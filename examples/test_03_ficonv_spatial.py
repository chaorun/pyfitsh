#!/usr/bin/env python3
"""test_03: ficonv spatial kernel 卷积/减影对比"""
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

    fc = Ficonv(kernel='i/1;b/1;g=5,2.0,2/1', iterations=2, rejection_level=3)
    convolved, mask, subtracted, kernel_dict = fc.fit(
        ref, v1_img, ref_mask=ref_mask, img_mask=v1_mask,
        output_subtracted=True, output_kernel_list=True)

    ok1, msg1 = compare_fits(convolved, "cli_test_ref_conv_spatial1.fits")
    if not ok1:
        return False, f"convolved: {msg1}"

    ok2, msg2 = compare_fits(subtracted, "cli_test_sub_spatial1.fits")
    if not ok2:
        return False, f"subtracted: {msg2}"

    n_kernels = len(kernel_dict.get('kernels', []))
    if n_kernels < 1:
        return False, "kernel_dict has no kernels"

    return True, f"convolved/subtracted ok, kernels={n_kernels}"


if __name__ == "__main__":
    sys.exit(run_test(do_test))
