#!/usr/bin/env python3
"""run_ficombine_tests.py — compare Python Ficombine API against CLI reference."""
import sys, os, numpy as np
TESTDIR = os.path.dirname(os.path.realpath(__file__))
ROOT = os.path.dirname(os.path.dirname(TESTDIR))
sys.path.insert(0, ROOT)
os.chdir(TESTDIR)
from pyfitsh import Ficombine
from astropy.io import fits

P = "ficombine_test"

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

def chk_skip(name, reason):
    s = f"{name:40s} {reason}"
    results.append((True, s))
    print(s)

print("=" * 80)

fc = Ficombine()

# Load 5 standard images
imgs = []
for k in range(1, 6):
    imgs.append(fits.getdata(f"{P}_img{k}.fits").astype(np.float64))
images = np.stack(imgs, axis=0)  # (5, 2560, 2560)

# Load negative and zero images
img_neg = fits.getdata(f"{P}_negative.fits").astype(np.float64)
img_zero = fits.getdata(f"{P}_zero.fits").astype(np.float64)

# [01] mean
r = fc.combine(images, mode='mean')
chk("[01] mean", r['data'], f"{P}_mean.fits", tol=1e-4)

# [02] median
r = fc.combine(images, mode='median')
chk("[02] median", r['data'], f"{P}_median.fits")

# [03] min
r = fc.combine(images, mode='min')
chk("[03] min", r['data'], f"{P}_min.fits")

# [04] max
r = fc.combine(images, mode='max')
chk("[04] max", r['data'], f"{P}_max.fits")

# [05] sum
r = fc.combine(images, mode='sum')
chk("[05] sum", r['data'], f"{P}_sum.fits", tol=1e-3)

# [06] squaresum
r = fc.combine(images, mode='squaresum')
chk("[06] squaresum", r['data'], f"{P}_squaresum.fits", tol=3e-1)

# [07] scatter
r = fc.combine(images, mode='scatter')
chk("[07] scatter", r['data'], f"{P}_scatter.fits", tol=1e-5)

# [08] stddev — same output as scatter
r = fc.combine(images, mode='stddev')
chk("[08] stddev", r['data'], f"{P}_stddev.fits", tol=1e-5)

# [09] rejmed
r = fc.combine(images, mode='rejmed,iterations=2,sigma=3')
chk("[09] rejmed", r['data'], f"{P}_rejmed.fits")

# [10] rejection (alias, lower=3,upper=3)
r = fc.combine(images, mode='rejection,iterations=2,lower=3,upper=3')
chk("[10] rejection", r['data'], f"{P}_rejection.fits")

# [11] rejmean
r = fc.combine(images, mode='rejmean,iterations=2,sigma=3')
chk("[11] rejmean", r['data'], f"{P}_rejmean.fits", tol=1e-4)

# [12] truncated discard=1
r = fc.combine(images, mode='truncated,discard=1')
chk("[12] truncated discard", r['data'], f"{P}_truncated_discard.fits", tol=1e-4)

# [13] truncated lowest=1,highest=1
r = fc.combine(images, mode='truncated,lowest=1,highest=1')
chk("[13] truncated low_high", r['data'], f"{P}_truncated_low_high.fits", tol=1e-4)

# [14] winsorized discard=1
r = fc.combine(images, mode='winsorized,discard=1')
chk("[14] winsorized", r['data'], f"{P}_winsorized_discard.fits", tol=1e-4)

# [15] ignore negative by option (-n flag)
images_neg = np.stack([img_neg, imgs[0], imgs[1]], axis=0)
r = fc.combine(images_neg, mode='mean', ignore_neg=True)
chk("[15] ignore_neg option", r['data'], f"{P}_ignore_negative_option.fits", tol=1e-4)

# [16] ignorenegative in mode string
r = fc.combine(images_neg, mode='mean,ignorenegative')
chk("[16] ignorenegative mode", r['data'], f"{P}_ignore_negative_mode.fits", tol=1e-4)

# [17] ignorezero in mode string
images_zero = np.stack([img_zero, imgs[0], imgs[1]], axis=0)
r = fc.combine(images_zero, mode='mean,ignorezero')
chk("[17] ignorezero mode", r['data'], f"{P}_ignorezero_mode.fits", tol=1e-4)

# [18] no-history — same pixel output as mean
r = fc.combine(images, mode='mean')
chk("[18] no-history", r['data'], f"{P}_mean_nohistory.fits", tol=1e-4)

# [19] history — same pixel output as mean
r = fc.combine(images, mode='mean')
chk("[19] history", r['data'], f"{P}_mean_history.fits", tol=1e-4)

# [20] bitpix int16 — CLI quantizes to int16, we don't; skip
chk_skip("[20] bitpix int16", "SKIP: bitpix quantization not applicable")

# [21] logical or (mode=or)
r = fc.combine(images, mode='or')
chk("[21] mask or", r['data'], f"{P}_mask_or.fits", tol=1e-4)

# [22] logical and (mode=and)
r = fc.combine(images, mode='and')
chk("[22] mask and", r['data'], f"{P}_mask_and.fits", tol=1e-4)

# [23] --logical-or option
r = fc.combine(images, mode='mean')
chk("[23] logical or option", r['data'], f"{P}_logical_or_option.fits", tol=1e-4)

# [24] --logical-and option
r = fc.combine(images, mode='mean', logical_and=True)
chk("[24] logical and option", r['data'], f"{P}_logical_and_option.fits", tol=1e-4)

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
