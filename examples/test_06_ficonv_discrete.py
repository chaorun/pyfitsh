#!/usr/bin/env python3
"""test_06: ficonv discrete kernel 对比"""
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

    fc = Ficonv(kernel='i/0;b/0;d=2/0', iterations=2, rejection_level=3)
    convolved, mask, subtracted, kernel_dict = fc.fit(
        ref, v1_img, ref_mask=ref_mask, img_mask=v1_mask,
        output_subtracted=True, output_kernel_list=True)

    ok1, msg1 = compare_fits(convolved, "cli_test_ref_conv_discrete.fits")
    if not ok1:
        return False, f"convolved: {msg1}"

    ok2, msg2 = compare_fits(subtracted, "cli_test_sub_discrete.fits")
    if not ok2:
        return False, f"subtracted: {msg2}"

    if len(kernel_dict.get('kernels', [])) < 1:
        return False, "kernel_dict has no kernels"

    return True, "ok"


if __name__ == "__main__":
    sys.exit(run_test(do_test))
