#!/usr/bin/env python3
"""test_07b: fistar deviated moment 列对比"""
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
                model='deviated', model_order=2,
                fields='id,x,y,bg,amp,s,d,k,mom,l,sigma,delta,kappa,fwhm,ellip,pa,flux,noise,s/n,magnitude')
    r = fs.do_fistar(ref_data)
    cy = r.output

    cli_names = 'id,x,y,bg,amp,s,d,k,mom,l,sigma,delta,kappa,fwhm,ellip,pa,flux,noise,sn,magnitude'.split(',')
    cli_tab = Table.read("cli_test_fistar_deviated_mom.cat", format='ascii.fast_no_header', names=cli_names)
    cy.sort('x')
    cli_tab.sort('x')
    diffs = []
    for col in cli_names[:8]:  # id,x,y,bg,amp,s,d,k
        cy_arr = np.asarray(cy[col], dtype=np.float64)
        cli_arr = np.asarray(cli_tab[col], dtype=np.float64)
        max_diff = np.abs(cy_arr - cli_arr).max()
        if max_diff > 1e-3:
            diffs.append(f"  {col}: max_diff={max_diff:.6e}")

    # mom: CY has single value, CLI has comma-separated; compare CY value against CLI first mom
    cy_mom = np.asarray(cy['mom'], dtype=np.float64)
    cli_mom_strs = np.asarray(cli_tab['mom'], dtype=str)
    cli_mom0 = np.array([float(s.split(',')[0].replace('+','')) for s in cli_mom_strs])
    if not np.allclose(cy_mom, cli_mom0, rtol=1e-3):
        diffs.append(f"  mom[0]: max_diff={np.abs(cy_mom - cli_mom0).max():.6e}")

    if diffs:
        return False, f"{len(diffs)} mismatch(es):\n" + "\n".join(diffs[:30])

    return True, f"table ok nstar={r.nstar}"


if __name__ == "__main__":
    sys.exit(run_test(do_test))
