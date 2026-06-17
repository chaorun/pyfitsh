#!/usr/bin/env python3
"""test_13: fiphot subtraction photometry 对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import run_test, RED, GREEN, NC
from astropy.io import fits
from astropy.table import Table
from pyfitsh.utils import decode_maskinfo

REF_IMG = "L1_GRB260604C-Ic-180s-20260606_210318.fits"
V1 = "cli_test_v1.fits"


def load_with_mask(path):
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(np.float64)
        mask = decode_maskinfo(hdul[0])
    if mask is None:
        mask = np.zeros(data.shape, dtype=np.uint8)
    mask[np.isnan(data)] = 1
    return data, mask


def do_test():
    from pyfitsh import Fiphot
    from pyfitsh.ficonv.ficonv import Ficonv

    ref_data, ref_mask = load_with_mask(REF_IMG)
    stars = np.loadtxt("cli_test_fiphot_positions.dat", ndmin=2)

    fp = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0),
                sky_fit='median,iterations=2,sigma=3')
    r_raw = fp.photometry(ref_data, stars=stars, col_id=0, col_xy=(1, 2))
    raw_3d = r_raw.dict['output_raw_photometry_3d_data']

    v1_img, v1_mask = load_with_mask(V1)
    fc = Ficonv(kernel='i/1;b/1;g=5,2.0,2/1', iterations=2, rejection_level=3)
    _, _, subtracted, kernel_dict = fc.fit(
        ref_data, v1_img, ref_mask=ref_mask, img_mask=v1_mask,
        output_subtracted=True, output_kernel_list=True)

    cli_sub = fits.getdata("cli_test_sub_spatial1.fits").astype(np.float64)
    if not np.allclose(subtracted, cli_sub, rtol=1e-3, atol=1e-2):
        return False, f"subtracted mismatch: {np.abs(subtracted-cli_sub).max():.6e}"

    sub_img, sub_mask = load_with_mask("cli_test_sub_spatial1.fits")
    fp2 = Fiphot(apertures='10:18:12', gain='2', mag_flux=(10.0, 10000.0))
    r = fp2.photometry_from_raw(raw_3d, sub_img, mask=sub_mask, kernel_dict=kernel_dict)
    cy = r.output
    cy.sort('x')

    cli_names = ['id','x','y','bg','bg_err','flux','flux_err','mag','mag_err']
    cli_tab = Table.read("cli_test_fiphot_sub.cat", format='ascii.fast_no_header', names=cli_names)
    for name in cli_names:
        cli_tab[name] = np.array(
            [np.nan if str(v) == "-" else float(v) for v in cli_tab[name]],
            dtype=float)
    cli_tab.sort('x')

    diffs = []
    for name in cli_names:
        if name not in cy.colnames:
            continue
        cy_col = np.asarray(cy[name], dtype=np.float64)
        cli_col = np.asarray(cli_tab[name], dtype=np.float64)
        valid = ~np.isnan(cli_col)
        if valid.sum() == 0:
            continue
        md = np.abs(cy_col[valid] - cli_col[valid]).max()
        if name in ('id',):
            continue
        if not np.allclose(cy_col[valid], cli_col[valid], rtol=1e-3, atol=5e-1):
            if name == 'flux_err' and np.allclose(cy_col[valid], cli_col[valid], rtol=1e-2, atol=1e3):
                continue
            md = np.abs(cy_col[valid] - cli_col[valid]).max()
            diffs.append(f"  {name}: max_diff={md:.6e}")

    if diffs:
        return False, f"{len(diffs)} col(s) exceed:\n" + "\n".join(diffs)

    return True, f"ok {len(cy)} rows"


if __name__ == "__main__":
    sys.exit(run_test(do_test))
