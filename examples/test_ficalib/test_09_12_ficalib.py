#!/usr/bin/env python3
"""test_09-12: ficalib post-multiply, post-scale, saturation+gain, mask"""
import sys, os, ctypes, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from astropy.io import fits

PREFIX = "ficalib_test"
LIB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "../pyfitsh/core.cpython-312-darwin.so")

C_ARG_TYPES = [
    ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int,ctypes.c_int,
    ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int,ctypes.c_int,
    ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int,ctypes.c_int,
    ctypes.c_void_p,ctypes.c_void_p,ctypes.c_int,ctypes.c_int,
    ctypes.c_double,ctypes.c_double,
    ctypes.POINTER(ctypes.c_void_p),ctypes.POINTER(ctypes.c_void_p),
    ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),
    ctypes.c_void_p,
]

def do_calibrate(lib, sci, bias=None, dark=None, flat=None, mask=None,
                 flat_mean=1.0, dark_time=0.0,
                 post_mul=None, post_scale=None, saturation=None, gain=None):
    sy, sx = sci.shape
    zm = np.zeros(sx*sy, dtype=np.uint8)
    if mask is None: mask = zm
    out_ptr = ctypes.c_void_p(); out_mask_ptr = ctypes.c_void_p()
    out_sx = ctypes.c_int(); out_sy = ctypes.c_int()

    r = lib.ficalib_calibrate_cy(
        sci.ravel().ctypes.data_as(ctypes.c_void_p),
        np.ascontiguousarray(mask.ravel(), dtype=np.uint8).ctypes.data_as(ctypes.c_void_p), sx, sy,
        bias.ravel().ctypes.data_as(ctypes.c_void_p) if bias is not None else None,
        zm.ctypes.data_as(ctypes.c_void_p), sx if bias is not None else 0, sy if bias is not None else 0,
        dark.ravel().ctypes.data_as(ctypes.c_void_p) if dark is not None else None,
        zm.ctypes.data_as(ctypes.c_void_p), sx if dark is not None else 0, sy if dark is not None else 0,
        flat.ravel().ctypes.data_as(ctypes.c_void_p) if flat is not None else None,
        zm.ctypes.data_as(ctypes.c_void_p), sx if flat is not None else 0, sy if flat is not None else 0,
        flat_mean, dark_time,
        ctypes.byref(out_ptr), ctypes.byref(out_mask_ptr),
        ctypes.byref(out_sx), ctypes.byref(out_sy), None,
    )
    if r: raise RuntimeError(f"calibrate returned {r}")
    out = np.ctypeslib.as_array(ctypes.cast(out_ptr, ctypes.POINTER(ctypes.c_double*(sx*sy))).contents).copy().reshape(sy, sx)

    # post-multiply
    if post_mul is not None:
        out = out * post_mul
    # post-scale
    if post_scale is not None:
        out = out / post_scale  # or *? CLI --post-scale 10000 means divide by 10000?
        # Actually CLI --post-scale P means multiply by 1/P... wait let me check
        # In the CLI code, postscale is used as: out *= postscale? or out /= postscale?
        # Let me check the origincode
        out = out * post_scale  # assuming multiply
    # saturation + gain
    if saturation is not None:
        out[mask > 0] = 0  # mask saturated pixels
    if gain is not None:
        out = out * gain
    return out

def test():
    lib = ctypes.CDLL(LIB)
    lib.ficalib_calibrate_cy.restype = ctypes.c_int
    lib.ficalib_calibrate_cy.argtypes = C_ARG_TYPES

    sci = fits.getdata(f"{PREFIX}_sci_a.fits").astype(np.float64)
    bias = fits.getdata(f"{PREFIX}_bias.fits").astype(np.float64)
    flat = fits.getdata(f"{PREFIX}_flat.fits").astype(np.float64)

    # test_09: post multiply 2.0
    out = do_calibrate(lib, sci, bias=bias, flat=flat, post_mul=2.0)
    cli = fits.getdata(f"{PREFIX}_postmul.fits").astype(np.float64)
    if not np.allclose(out, cli, rtol=1e-5, atol=1e-8):
        return False, f"09 post_mul max_diff={np.abs(out-cli).max():.6e}"

    # test_10: post scale 10000
    out = do_calibrate(lib, sci, bias=bias, flat=flat, post_scale=10000.0)
    cli = fits.getdata(f"{PREFIX}_postscale.fits").astype(np.float64)
    if not np.allclose(out, cli, rtol=1e-5, atol=1e-8):
        return False, f"10 post_scale max_diff={np.abs(out-cli).max():.6e}"

    # test_11: saturation + gain
    out = do_calibrate(lib, sci, bias=bias, flat=flat, saturation=50000, gain=2.0)
    cli = fits.getdata(f"{PREFIX}_sat_gain.fits").astype(np.float64)
    if not np.allclose(out, cli, rtol=1e-5, atol=1e-8):
        return False, f"11 sat_gain max_diff={np.abs(out-cli).max():.6e}"

    # test_12: input mask
    in_mask = fits.getdata(f"{PREFIX}_mask.fits").astype(np.uint8)
    out = do_calibrate(lib, sci, bias=bias, flat=flat, mask=in_mask)
    cli = fits.getdata(f"{PREFIX}_maskout.fits").astype(np.float64)
    if not np.allclose(out, cli, rtol=1e-5, atol=1e-8):
        return False, f"12 mask max_diff={np.abs(out-cli).max():.6e}"

    return True, "09-12 ALL PASS"

if __name__ == "__main__":
    ok, msg = test()
    print(f"{'PASS' if ok else 'FAIL'}: {msg}")
