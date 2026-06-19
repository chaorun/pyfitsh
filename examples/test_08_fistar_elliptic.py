#!/usr/bin/env python3
"""test_08: fistar elliptic 对比"""
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
    fs = Fistar(flux_threshold=10000, skysigma=5, gain=2,
                model='elliptic',
                fields='id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,s/n,magnitude')
    r = fs.do_fistar(ref_data)

    cli_names = 'id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,sn,magnitude'.split(',')
    cli_tab = Table.read("cli_test_fistar_elliptic.cat", format='ascii.fast_no_header', names=cli_names)

    cy = r.output
    cy.sort('x')
    cli_tab.sort('x')
    matched = 0
    diffs = []
    for col in cli_names:
        if col not in cy.colnames:
            diffs.append(f"column {col} not in cy output")
            continue
        cy_arr = np.asarray(cy[col], dtype=np.float64)
        cli_arr = np.asarray(cli_tab[col], dtype=np.float64)
        max_diff = np.abs(cy_arr - cli_arr).max()
        tol = 0.05 if col in ('s','d','k','ellip','pa','sn','noise','magnitude') else 0.01
        if max_diff > tol:
            diffs.append(f"  {col}: max_diff={max_diff:.6e}")
        if len(diffs) >= 30:
            break
    if diffs:
        return False, f"{len(diffs)} mismatch(es):\n" + "\n".join(diffs[:30])

    return True, f"table ok nstar={r.nstar}"

if __name__ == "__main__":
    sys.exit(run_test(do_test))
