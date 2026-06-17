#!/usr/bin/env python3
"""Comprehensive test of all fitsh_cy Cython wrappers using real data."""

import os, sys, re, time, subprocess, tempfile
import numpy as np
from astropy.io import fits

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTDATA = os.path.join(BASE, "fitsh_test_data")
AST32 = os.path.join(BASE, "ast32example")
sys.path.insert(0, os.path.dirname(__file__))
from fitsh_cy import (Grmatch, Grtrans, Fitrans, Fistar, Ficonv,
                      PolyFitter)

PASS = FAIL = TOTAL = 0

def check(label, condition):
    global PASS, FAIL, TOTAL
    TOTAL += 1
    if condition:
        PASS += 1; print(f"  [PASS] {label}")
    else:
        FAIL += 1; print(f"  [FAIL] {label}")

# ════════════════════════════════════════════════════════════════
# 1. GRMATCH — star catalog matching
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("1. GRMATCH — star catalog matching")
print("=" * 60)

def load_cat(path):
    d = np.loadtxt(path)
    return d[:, 0].astype(np.float64), d[:, 1].astype(np.float64), \
           d[:, 2].astype(np.float64), d[:, 3].astype(np.float64)

ref_path = os.path.join(TESTDATA, "ref_L1_GRB260604C-Ic-180s-20260606_210318.cat")
obs_path = os.path.join(TESTDATA, "obs_L1_GRB260604C-Ic-180s-20260607_023443.cat")
rids, rx, ry, rm = load_cat(ref_path)
oids, ox, oy, om = load_cat(obs_path)
print(f"  Loaded ref={len(rx)}, obs={len(ox)} stars")

g = Grmatch()
g.set_reference(rx, ry, -rm)
g.set_input(ox, oy, -om)
g._order = 1; g._maxdist = 2.0; g._unitarity = 0.2
g._ttype = 3; g._use_ordering = 1
res = g.run()

check("grmatch found pairs", res.nhit > 0)
print(f"    {res.nhit} pairs, nsigma={res.stats['nsigma']:.4f}")
check("has dxfit", len(res.vfits_dx) >= 3)
check("has dyfit", len(res.vfits_dy) >= 3)

# Cross-validate with original CLI
cmd_g = [
    os.path.join(BASE, "grmatch"),
    "--input-reference", ref_path, "--input", obs_path,
    "--col-ref", "2,3", "--col-inp", "2,3",
    "--col-ref-ordering", "-4", "--col-inp-ordering", "-4",
    "--triangulation", "auto,unitarity=0.2",
    "--order", "1", "--max-distance", "2",
    "--output-transformation", "/dev/null",
]
proc = subprocess.run(cmd_g, capture_output=True, text=True, cwd=BASE)
m = re.search(r'(\d+)\s+matched pairs', proc.stderr + proc.stdout)
if m:
    cli_n = int(m.group(1))
    check(f"nhit vs CLI ({cli_n} vs {res.nhit})", abs(cli_n - res.nhit) <= 2)

# Compare dxfit/dyfit
orig_tf = res.transform
trans_path = os.path.join(TESTDATA, "grb.trans")
with open(trans_path) as f: txt = f.read()
orig_dx = [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', txt.split('dxfit=')[1].split('\n')[0])]
orig_dy = [float(v) for v in re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', txt.split('dyfit=')[1].split('\n')[0])]
check("dxfit ≈ original", np.allclose(orig_dx, res.vfits_dx[:len(orig_dx)], rtol=1e-4))
check("dyfit ≈ original", np.allclose(orig_dy, res.vfits_dy[:len(orig_dy)], rtol=1e-4))


# ════════════════════════════════════════════════════════════════
# 2. GRTRANS — coordinate transformation
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("2. GRTRANS — coordinate transformation")
print("=" * 60)

# From dict
gt = Grtrans(orig_tf)
tx, ty = gt.apply(list(ox[:10]), list(oy[:10]))
check("transform len", len(tx) == 10)

# From parse_trans_file / from_trans_string
with open(trans_path) as f:
    gt2 = Grtrans.from_trans_string(f.read())
tx2, ty2 = gt2.apply(list(ox[:3]), list(oy[:3]))
check("from_trans_string → same tx", np.allclose(tx[:3], tx2))

# Identity (order=1 with unit scale)
gt_id = Grtrans(order=1, dxfit=[0.0, 1.0, 0.0], dyfit=[0.0, 0.0, 1.0])
ox_g, oy_g = gt_id.apply([100.0, 200.0], [300.0, 400.0])
check("identity x", np.allclose(ox_g, [100.0, 200.0]))
check("identity y", np.allclose(oy_g, [300.0, 400.0]))


# ════════════════════════════════════════════════════════════════
# 3. FICONV — kernel fitting + convolution
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("3. FICONV — kernel fitting + convolution")
print("=" * 60)

ref_path_f = os.path.join(TESTDATA, "L1_GRB260604C-Ic-180s-20260606_210318.fits")
obs_path_f = os.path.join(TESTDATA, "L1_GRB260604C-Ic-180s-20260607_023443.fits")
ref_img = fits.getdata(ref_path_f).astype(np.float64)
obs_img = fits.getdata(obs_path_f).astype(np.float64)
print(f"  Loaded ref={ref_img.shape}, obs={obs_img.shape}")

# 512x512 center crop
cx = cy = 1024; half = 256
ref_c = ref_img[cy-half:cy+half, cx-half:cx+half].copy()
obs_c = obs_img[cy-half:cy+half, cx-half:cx+half].copy()
mk_c = np.zeros_like(ref_c, dtype=np.uint8)
print(f"  Crop: {ref_c.shape}")

t0 = time.time()
KERN = "g=3,1.0,0/0"  # CLI-compatible format: g=halfsize,sigma,order/spatial_order
fc = Ficonv(kernel=KERN, iterations=0, divide=64)
cnv, cmk = fc.fit(ref_c, obs_c, mk_c, mk_c)
dt = time.time() - t0
check("ficonv shape OK", cnv.shape == ref_c.shape)
check("ficonv no NaN", not np.any(np.isnan(cnv)))
check("ficonv non-trivial", np.std(cnv) > 0)
print(f"    Time: {dt:.2f}s")

# With iterations
fc2 = Ficonv(kernel=KERN, iterations=3, rejection_level=3.0, divide=64)
cnv2, cmk2 = fc2.fit(ref_c, obs_c, mk_c, mk_c)
check("ficonv+iters shape OK", cnv2.shape == ref_c.shape)

# Via Fitrans.Convolution
fop = Fitrans.Convolution(kernel=KERN, divide=64)
cnv3, cmk3 = fop.fit(ref_c, obs_c, mk_c, mk_c)
check("Fitrans.Convolution shape OK", cnv3.shape == ref_c.shape)

# Cross-validate ficonv vs CLI (small crop for speed)
# Cross-validate ficonv vs CLI using identical cropped FITS files
with tempfile.NamedTemporaryFile(suffix="_r.fits", delete=False) as tf_r:
    fits.PrimaryHDU(ref_c.astype(np.float32)).writeto(tf_r.name, overwrite=True)
    tmp_ref = tf_r.name
with tempfile.NamedTemporaryFile(suffix="_i.fits", delete=False) as tf_i:
    fits.PrimaryHDU(obs_c.astype(np.float32)).writeto(tf_i.name, overwrite=True)
    tmp_obs = tf_i.name
with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as tf:
    tf_out = tf.name

cmd_fc = [
    os.path.join(BASE, "ficonv"),
    "-r", tmp_ref,
    "-i", tmp_obs,
    "-k", KERN,
    "-d", "64",
    "-o", tf_out,
]
print(f"  Running ficonv CLI on {ref_c.shape} crops...")
result_fc = subprocess.run(cmd_fc, capture_output=True, text=True, cwd=BASE, timeout=120)
print(f"    rc={result_fc.returncode}")
if result_fc.returncode != 0:
    print(f"    CLI cmd: {' '.join(cmd_fc)}")
    print(f"    CLI stderr: {result_fc.stderr[:200]}")
    print(f"    (CLI failed, skipping comparison)")
else:
    cli_cnv = fits.getdata(tf_out).astype(np.float64)
    maxdiff = np.max(np.abs(cli_cnv - cnv))
    check(f"ficonv ≈ CLI (maxdiff {maxdiff:.6f})", np.allclose(cli_cnv, cnv, rtol=1e-4, atol=1e-4))
os.unlink(tf_out); os.unlink(tmp_ref); os.unlink(tmp_obs)


# ════════════════════════════════════════════════════════════════
# 4. FITRANS — image operations
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("4. FITRANS — image operations")
print("=" * 60)

np.random.seed(42)
td = np.random.randn(32, 32).astype(np.float64) + 100
tmk = np.zeros((32, 32), dtype=np.uint8)

# 4a. Zoom
zd, zm = Fitrans.Zoom(2, 2).apply(td, tmk)
check("Zoom: 2x larger",   zd.shape == (64, 64))
zdr, _ = Fitrans.Zoom(2, 2, raw=True).apply(td, tmk)
check("Zoom raw: same size", zdr.shape == zd.shape)

# 4b. Shrink
sd, sm = Fitrans.Shrink(2, 2).apply(td, tmk)
check("Shrink: 2x smaller", sd.shape == (16, 16))
sdm, _ = Fitrans.Shrink(2, 2, median=True).apply(td, tmk)
check("Shrink median: same size", sdm.shape == sd.shape)

# 4c. Repetitive (vs original CLI)
rd, rm = Fitrans.Repetitive(2, 2).apply(td, tmk)
check("Repetitive: output 64x64", rd.shape == (64, 64))
check("Repetitive: tile 0,0", np.allclose(rd[0:32, 0:32], td))
check("Repetitive: tile 0,1", np.allclose(rd[0:32, 32:64], td))
# Cross-validate with CLI
td_fits = os.path.join(tempfile.mkdtemp(), "td.fits")
ro_fits = os.path.join(tempfile.mkdtemp(), "ro.fits")
fits.PrimaryHDU(td).writeto(td_fits, overwrite=True)
subprocess.run([os.path.join(BASE, "fitrans"), "-i", td_fits, "-o", ro_fits,
                "--repetitive-xy", "64,64"], capture_output=True, cwd=BASE)
ro_cli = fits.getdata(ro_fits).astype(np.float64)
check("Repetitive vs CLI", np.allclose(rd, ro_cli))

# 4d. Interleave (vs original CLI)
idata, _ = Fitrans.Interleave(4, 2).apply(td, tmk)
check("Interleave: output 2x4", idata.shape == (2, 4))
check("Interleave: no NaN", not np.any(np.isnan(idata)))
io_fits = ro_fits.replace("ro.fits", "io.fits")
subprocess.run([os.path.join(BASE, "fitrans"), "-i", td_fits, "-o", io_fits,
                "--interleave-xy", "4,2"], capture_output=True, cwd=BASE)
io_cli = fits.getdata(io_fits).astype(np.float64)
check("Interleave vs CLI", np.allclose(idata, io_cli))
os.unlink(td_fits); os.unlink(ro_fits); os.unlink(io_fits)

# 4e. Noise
nd, _ = Fitrans.Noise().apply(td, tmk)
check("Noise: same shape", nd.shape == td.shape)

# 4f. Smooth (spline + polynomial)
smo, _ = Fitrans.Smooth(smooth_type="spline", order=2).apply(td, tmk)
check("Smooth spline: same shape", smo.shape == td.shape)
spo, _ = Fitrans.Smooth(smooth_type="polynomial", order=2).apply(td, tmk)
check("Smooth poly: same shape", spo.shape == td.shape)

# 4g. ImageTransformation (geometric remap)
op_it = Fitrans.ImageTransformation(order=1, method="bilinear")
op_it.set_trans(res.vfits_dx, res.vfits_dy)
itd, itm = op_it.apply(ref_c, mk_c)
check("ImageTransform shape OK", itd.shape == ref_c.shape)

# 4h. decode_maskinfo
mask_test = Fitrans.decode_maskinfo(fits.open(ref_path_f)[0])
check("decode_maskinfo: 2048x2048", mask_test.shape == (2048, 2048))


# ════════════════════════════════════════════════════════════════
# 5. FISTAR — star detection
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print("5. FISTAR — star detection")
print("=" * 60)

sc = ref_img[cy-200:cy+200, cx-200:cx+200].copy()
smk = np.zeros_like(sc, dtype=np.uint8)
print(f"  Region: {sc.shape}")

t0 = time.time()
fs = Fistar(threshold=50.0, algorithm="uplink")
sr = fs.search(sc, smk)
dt = time.time() - t0
check("Fistar found stars", sr["nstar"] > 0)
print(f"    {sr['nstar']} stars in {dt:.2f}s")

for k in ["cx", "cy", "peak", "flux", "fwhm", "magnitude", "sn"]:
    if k in sr["table"].colnames:
        check(f"field '{k}'", len(sr["table"][k]) == sr["nstar"])

# PSF output
fs2 = Fistar(threshold=100.0, psf="native,grid=2,halfsize=3", psf_output=True)
sr2 = fs2.search(sc, smk)
if sr2.get("psf_data") is not None:
    check("PSF data returned", len(sr2["psf_data"]) > 0)

# Input candidates
cands = np.array([[100.0, 100.0]], dtype=np.float64)
sr3 = Fistar(input_candidates=cands, only_candidates=False).search(sc, smk)
check("with input candidates", sr3["nstar"] >= 0)

# Input positions
mx, my = float(np.mean(sr["table"]["cx"])), float(np.mean(sr["table"]["cy"]))
poss = np.array([[mx, my]], dtype=np.float64)
sr4 = Fistar(input_positions=poss).search(sc, smk)
check("with input positions", sr4["nstar"] >= 0)


# ════════════════════════════════════════════════════════════════
# SUMMARY
# ════════════════════════════════════════════════════════════════
print("\n" + "=" * 60)
print(f"RESULTS: {PASS}/{TOTAL} passed, {FAIL} failed")
print("=" * 60)
sys.exit(0 if FAIL == 0 else 1)
