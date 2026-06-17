#!/usr/bin/env python3
"""test_07: fistar gauss full + PSF 对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import compare_fits, run_test, RED, GREEN, NC
from astropy.io import fits
from astropy.table import Table

REF_IMG = "L1_GRB260604C-Ic-180s-20260606_210318.fits"

def do_test():
    from pyfitsh import Fistar

    ref_data = fits.getdata(REF_IMG).astype(np.float64)
    fs = Fistar(flux_threshold=10000, skysigma=5, sort='x',
                model='gauss', model_order=2, it_sym=4, it_gen=2, gain=2,
                psf='native,order=2', psf_output=True, verbose=False,
                fields='id,ix,iy,cx,cy,cbg,camp,cmax,npix,cs,cd,ck,x,y,bg,amp,s,d,k,l,sigma,delta,kappa,fwhm,ellip,pa,flux,noise,s/n,magnitude,px,py,pbg,pamp,ps,pd,pk,pl')
    r = fs.do_fistar(ref_data)

    # 对比星表：CLI 文件无表头，按 CLI --format 顺序命名
    cli_names = 'id,ix,iy,cx,cy,cbg,camp,cmax,npix,cs,cd,ck,x,y,bg,amp,s,d,k,l,sigma,delta,kappa,fwhm,ellip,pa,flux,noise,sn,magnitude,px,py,pbg,pamp,ps,pd,pk,pl'.split(',')
    cli_tab = Table.read("cli_test_fistar_gauss_full.cat", format='ascii.fast_no_header', names=cli_names)

    cy = r.output
    cy_cols_map = list(cy.colnames)
    diffs = []
    for col in cli_names:
        if col not in cy_cols_map:
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

    # 对比 PSF 3D cube
    ok_psf, msg_psf = compare_fits(r.psf, "cli_test_native_o2.psf.fits")
    if not ok_psf:
        return False, f"psf: {msg_psf}"

    return True, f"table ok nstar={r.nstar}, psf ok"

if __name__ == "__main__":
    sys.exit(run_test(do_test))
