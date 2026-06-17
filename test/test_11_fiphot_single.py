#!/usr/bin/env python3
"""test_11: fiphot single aperture 对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import run_test, RED, GREEN, NC
from astropy.io import fits
from astropy.table import Table

REF_IMG = "L1_GRB260604C-Ic-180s-20260606_210318.fits"

def do_test():
    from pyfitsh import Fiphot

    ref_data = fits.getdata(REF_IMG).astype(np.float64)
    stars = np.loadtxt("cli_test_fiphot_positions.dat", ndmin=2)
    # stars columns: id(0), x(1), y(2)
    fp = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0),
                sky_fit='median,iterations=2,sigma=3')
    r = fp.photometry(ref_data, stars=stars, col_id=0, col_xy=(1, 2))

    cy = r.output
    # CY column order: id,x,y,flux,flux_err,bg,bg_err,mag,mag_err
    # CLI column order: id,x,y,bg,bg_err,flux,flux_err,mag,mag_err
    # Map CY col -> CLI col index
    cy_colnames = cy.colnames  # ['id','x','y','flux','flux_err','bg','bg_err','mag','mag_err']
    cy_arr = np.column_stack([np.asarray(cy[c], dtype=np.float64) for c in cy_colnames])
    # build column index map
    cy_idx = {c: i for i, c in enumerate(cy_colnames)}
    cli_order = ['id','x','y','bg','bg_err','flux','flux_err','mag','mag_err']
    cy_reordered = np.column_stack([cy_arr[:, cy_idx[c]] for c in cli_order])

    cli_data = np.loadtxt("cli_test_fiphot_single.cat", dtype=np.float64)

    if cy_reordered.shape != cli_data.shape:
        return False, f"shape: CY={cy_reordered.shape} CLI={cli_data.shape}"

    for i, col_name in enumerate(cli_order):
        cy_col = cy_reordered[:, i]
        cli_col = cli_data[:, i]
        if not np.allclose(cy_col, cli_col, rtol=1e-3, atol=1e-2, equal_nan=True):
            max_diff = np.abs(cy_col - cli_col).max()
            return False, f"col {col_name}: max_diff={max_diff:.6e}"

    return True, f"table ok {len(cy)} rows"

if __name__ == "__main__":
    sys.exit(run_test(do_test))
