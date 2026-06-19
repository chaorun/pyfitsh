#!/usr/bin/env python3
"""run_fiign_tests.py — compare Python Fiign API against CLI reference outputs."""
import sys, os, numpy as np
TESTDIR = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.dirname(os.path.dirname(TESTDIR))
sys.path.insert(0, ROOT)
os.chdir(TESTDIR)
from pyfitsh import Fiign
from astropy.io import fits

P = "fiign_test"

results = []

def chk(name, py_data, cli_path, tol=0):
    cli = fits.getdata(cli_path).astype(np.float64)
    d = np.abs(py_data - cli)
    md = d.max()
    nd = (d > tol).sum()
    ok = nd == 0
    s = f"{name:40s} max={md:.3e} n_diff={nd:>8}/{py_data.size} tol={tol:.0e}"
    results.append((ok, s))
    print(s, "  OK" if ok else "  FAIL")

def chk_mask(name, py_mask, cli_path, tol=0):
    cli = fits.getdata(cli_path).astype(np.int32)
    d = np.abs(py_mask.astype(np.int32) - cli)
    md = d.max()
    nd = (d > tol).sum()
    ok = nd == 0
    s = f"{name:40s} max={md:.3e} n_diff={nd:>8}/{py_mask.size} tol={tol:.0e}"
    results.append((ok, s))
    print(s, "  OK" if ok else "  FAIL")

# ====== load input images ======
img_base = fits.getdata(f"{P}_base.fits").astype(np.float64)
img_signed = fits.getdata(f"{P}_signed.fits").astype(np.float64)
img_satlevel = fits.getdata(f"{P}_satlevel.fits").astype(np.float64)

print("=" * 80)

# [01] ignore negative + apply mask
r = Fiign().apply(img_signed, ignore_neg=True, apply_mask=True, mask_value=-999)
chk("[01] neg apply data", r['data'], f"{P}_ignore_negative_apply.fits")

# [02] ignore zero + apply mask
r = Fiign().apply(img_signed, ignore_zero=True, apply_mask=True, mask_value=-888)
chk("[02] zero apply data", r['data'], f"{P}_ignore_zero_apply.fits")

# [03] ignore nonpositive + apply mask
r = Fiign().apply(img_signed, ignore_nonpos=True, apply_mask=True, mask_value=-777)
chk("[03] nonpos apply data", r['data'], f"{P}_ignore_nonpositive_apply.fits")

# [04] saturation scalar + apply mask
r = Fiign().apply(img_base, saturation=1200, apply_mask=True, mask_value=-666)
chk("[04] sat scalar data", r['data'], f"{P}_saturation_scalar_apply.fits")

# [05] saturation image + apply mask
r = Fiign().apply(img_base, saturation_img=img_satlevel, apply_mask=True, mask_value=-555)
chk("[05] sat image data", r['data'], f"{P}_saturation_image_apply.fits")

# [06] block rectangle (hot:100,100:300,300)
r = Fiign().apply(img_base, mask_block_list=["block:hot:100,100:300,300"], apply_mask=True, mask_value=-444)
chk("[06] block rect data", r['data'], f"{P}_block_rect_apply.fits")

# [07] block circle (cosmic:1200,1200:100)
r = Fiign().apply(img_base, mask_block_list=["circle:cosmic:1200,1200:100"], apply_mask=True, mask_value=-333)
chk("[07] block circle data", r['data'], f"{P}_block_circle_apply.fits")

# [08] block pixel (fault:500,500)
r = Fiign().apply(img_base, mask_block_list=["pixel:fault:500,500"], apply_mask=True, mask_value=-222)
chk("[08] block pixel data", r['data'], f"{P}_block_pixel_apply.fits")

# [09] block line (cosmic:100,200:800,900)
r = Fiign().apply(img_base, mask_block_list=["line:cosmic:100,200:800,900"], apply_mask=True, mask_value=-111)
chk("[09] block line data", r['data'], f"{P}_block_line_apply.fits")

# [10] expand mask (pixel hot:1000,1000 expand 5)
r = Fiign().apply(img_base, mask_block_list=["pixel:hot:1000,1000"], expand_hsize=5, apply_mask=True, mask_value=-123)
chk("[10] expand apply data", r['data'], f"{P}_expand_apply.fits")

# [11] external mask creation reference — not a comparison test
# The CLI creates external_mask.fits by applying -q block:fault:400,400:600,600 to base
# We load it for test 12
print("[11] external mask file: prepared by CLI (no comparison)")

# [12] input mask co-add + apply mask
# CLI: -M external_mask.fits reads mask from FITS header (MASKINFO)
# The external mask was created by: -q block:fault:400,400:600,600
# Python: manually create the mask array (no FITS header I/O)
ext_mask = np.zeros((2560, 2560), dtype=np.uint8)
ext_mask[400:601, 400:601] = 1  # MASK_FAULT in block 400,400:600,600
r = Fiign().apply(img_base, mask=ext_mask, apply_mask=True, mask_value=-321)
chk("[12] input mask data", r['data'], f"{P}_input_mask_apply.fits")

# [13] convert hot → cosmic
r = Fiign().apply(img_base, mask_block_list=["pixel:hot:1300,1300"],
                  convert_list=["hot:hot:none:cosmic"],
                  apply_mask=True, mask_value=-231)
chk("[13] convert data", r['data'], f"{P}_convert_hot_cosmic_apply.fits")

# [14] ignore input mask
# CLI: --ignore-mask means don't use the embedded mask from the FITS file
# Python: pass mask=None (don't load from file)
r = Fiign().apply(img_base, mask=None, apply_mask=True, mask_value=-135)
chk("[14] ignore mask data", r['data'], f"{P}_ignore_mask_apply.fits")

# ====== summary ======
print("=" * 80)
passed = sum(1 for ok, _ in results if ok)
failed = sum(1 for ok, _ in results if not ok)
print(f"{passed}/{len(results)} PASS, {failed}/{len(results)} FAIL")

if failed > 0:
    print("\nFAILED tests:")
    for ok, s in results:
        if not ok:
            print(f"  {s}")
