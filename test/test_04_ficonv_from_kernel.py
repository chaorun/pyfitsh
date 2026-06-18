#!/usr/bin/env python3
"""test_04: ficonv from fitted kernel 对比"""
import sys, os, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '..')
from compare_utils import compare_fits, run_test, RED, GREEN, NC
from astropy.io import fits
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
    from pyfitsh.ficonv.ficonv import Ficonv

    ref, ref_mask = load_with_mask(REF_IMG)
    v1_img, v1_mask = load_with_mask(V1)

    # Step 1: fit kernels
    fc = Ficonv(kernel='i/1;b/1;g=5,2.0,2/1', iterations=2, rejection_level=3)
    _, _, _, kernel_dict = fc.fit(
        ref, v1_img, ref_mask=ref_mask, img_mask=v1_mask,
        output_subtracted=True, output_kernel_list=True)

    # Step 2: apply fitted kernels
    convolved, _ = fc.fit(ref, v1_img, ref_mask=ref_mask, img_mask=v1_mask,
                           kernel_dict=kernel_dict)

    ok, msg = compare_fits(convolved, "cli_test_ref_conv_from_kernel.fits", rtol=1e-4, atol=1e-1)
    if not ok:
        return False, f"convolved: {msg}"

    return True, "ok"


if __name__ == "__main__":
    sys.exit(run_test(do_test))
