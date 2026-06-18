#!/usr/bin/env python3
"""test_12: fiphot multi aperture 对比"""
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
    stars_tab = Table.read("cli_test_fiphot_positions.dat",
                           format='ascii.no_header',
                           names=['id', 'x', 'y'])
    stars = np.column_stack([stars_tab['id'], stars_tab['x'], stars_tab['y']])

    fp = Fiphot(apertures='8:18:12,10:18:12,12:20:15,15:25:15',
                gain='2', mag_flux=(10.0, 10000.0))
    r = fp.photometry(ref_data, stars=stars, col_id=0, col_xy=(1, 2))

    nap = 4
    nstar = len(stars)

    # CY per-aperture order: id,x,y,flux,flux_err,bg,bg_err,mag,mag_err
    # CLI per-aperture order: id,x,y,bg,bg_err,flux,flux_err,mag,mag_err
    cy_cols = ['id','x','y','bg','bg_err','flux','flux_err','mag','mag_err']
    cy_raw_cols = ['id','x','y','flux','flux_err','bg','bg_err','mag','mag_err']
    cy_reorder = [0, 1, 2, 5, 6, 3, 4, 7, 8]  # cy_raw_cols index -> cy_cols order

    # CLI: first aperture 9 cols, subsequent 6 cols each → 9+6*3=27
    cli_names = ['id','x','y','bg','bg_err','flux','flux_err','mag','mag_err']
    for a in range(2, nap+1):
        cli_names += [f'{c}_a{a}' for c in ['bg','bg_err','flux','flux_err','mag','mag_err']]
    cli_tab = Table.read("cli_test_fiphot_multiap.cat", format='ascii.no_header',
                         names=cli_names)
    cli_data = np.column_stack([np.asarray(cli_tab[c], dtype=np.float64) for c in cli_names])
    nstar = cli_data.shape[0]

    cy_nap = np.zeros(cli_data.shape)
    for j in range(nap):
        for k, ci in enumerate(cy_reorder):
            cname = cy_raw_cols[ci]
            cy_flat = np.asarray(r.output[cname], dtype=np.float64).flatten()
            cy_val = cy_flat[j::nap]
            if j == 0:
                cy_nap[:, k] = cy_val
            else:
                if k >= 3:
                    cy_nap[:, 9 + (j-1)*6 + (k-3)] = cy_val

    if not np.allclose(cy_nap, cli_data, rtol=1e-3, atol=1e-2, equal_nan=True):
        max_diff = np.abs(cy_nap - cli_data).max()
        return False, f"max_diff={max_diff:.6e}"

    return True, f"table ok {nstar} stars × {nap} apertures"

if __name__ == "__main__":
    sys.exit(run_test(do_test))
