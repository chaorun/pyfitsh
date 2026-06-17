#!/usr/bin/env python3
"""test_10: fistar input_candidates 两步路径对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import run_test, RED, GREEN, NC
from astropy.io import fits
from astropy.table import Table

REF_IMG = "L1_GRB260604C-Ic-180s-20260606_210318.fits"

def do_test():
    from pyfitsh import Fistar

    ref_data = fits.getdata(REF_IMG).astype(np.float64)

    fs1 = Fistar(flux_threshold=10000, skysigma=5, only_candidates=True,
                 fields='cx,cy,cs,cd,ck')
    r1 = fs1.do_fistar(ref_data)

    
    r1.output.sort('cx')
    outcat = Table.read("cli_test_input_candidates.dat", format='ascii.fast_no_header',names=('cx','cy','cs','cd','ck'))
    outcat.sort('cx')
    for name in ['cx','cy','cs','cd','ck']:
        print(name,np.median(outcat[name]-r1.output[name]))
    #print(outcat[0])
    #print(r1.output[0])
    #      1  842.482   40.226     0.00     0.00  0.148  0.000  0.000 6.121 0.000   0.0       0.00     0.00        - 
    #print(outcat)
    cands = np.column_stack([
        r1.output['cx'], r1.output['cy'],
        r1.output['cs'], r1.output['cd'], r1.output['ck']
    ])
    
    #cands = np.column_stack([outcat['cx'], outcat['cy'],outcat['cs'], outcat['cd'], outcat['ck']])
    #sys.exit(0)
    # col_xy=(1,2) 0-indexed => cx 是 col 0, cy 是 col 1 => col_xy=(0,1)? 
    # CLI: --col-xy 1,2 (1-indexed) => 对应 cands 的第 1,2 列 (0-indexed)
    # 但 cands 布局是: cx,cy,cs,cd,ck (5列)
    # col_xy=(0,1) 0-indexed 对应 cx,cy; col_shape=(2,3,4) 0-indexed 对应 cs,cd,ck
    fs2 = Fistar(flux_threshold=10000, skysigma=5, gain=2,
                 model='elliptic',
                 fields='id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,s/n',
                 input_candidates=cands)
    r2 = fs2.do_fistar(ref_data)
    # 这里退出没事
    print(r2.output) #--- 加了这一行，或者说调用了r2.output，程序就会退出，原因不明
    # lldb？
    sys.exit(0)
    cli_names = 'id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,sn'.split(',')
    cli_tab = Table.read("cli_test_from_candidates.cat", format='ascii.fast_no_header', names=cli_names)

    cy = r2.output
    diffs = []
    for col in cli_names:
        if col not in cy.colnames:
            diffs.append(f"column {col} not in cy output")
            continue
        cy_arr = np.asarray(cy[col], dtype=np.float64)
        cli_arr = np.asarray(cli_tab[col], dtype=np.float64)
        max_diff = np.abs(cy_arr - cli_arr).max()
        if max_diff > 1e-3:
            diffs.append(f"  {col}: max_diff={max_diff:.6e}")
        if len(diffs) >= 30:
            break
    if diffs:
        return False, f"{len(diffs)} mismatch(es):\n" + "\n".join(diffs[:30])

    return True, f"table ok nstar={r2.nstar}"

if __name__ == "__main__":
    sys.exit(run_test(do_test))
