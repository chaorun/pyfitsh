#!/usr/bin/env python3
import sys, os, re
import numpy as np

TDIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(TDIR, '..', '..'))
from pyfitsh.gropt import Optics

def load_opt(name):
    with open(os.path.join(TDIR, name)) as f:
        return f.read()

def load_transfer(name):
    d = {}
    with open(os.path.join(TDIR, name)) as f:
        for line in f:
            m = re.match(r'(\w+):\s+([\d.eE+-]+)', line)
            if m:
                d[m.group(1)] = float(m.group(2))
    return d

def load_spot(name):
    rows = []
    center = None
    with open(os.path.join(TDIR, name)) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                parts = line[1:].split()
                if len(parts) >= 2:
                    center = (float(parts[0]), float(parts[1]))
            elif line:
                parts = line.split()
                if len(parts) >= 2:
                    rows.append([float(parts[0]), float(parts[1])])
    return np.array(rows), center

def load_psf(name):
    rows = []
    with open(os.path.join(TDIR, name)) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                rows.append([int(parts[0]), int(parts[1]), float(parts[2])])
    data = np.array(rows)
    px_max = int(data[:, 0].max()) + 1
    py_max = int(data[:, 1].max()) + 1
    psf = np.zeros((py_max, px_max))
    for r in data:
        psf[int(r[1]), int(r[0])] = r[2]
    return psf

passed = 0
failed = 0
total = 0

def report(num, name, ok, detail=""):
    global passed, failed, total
    total += 1
    tag = "PASS" if ok else "FAIL"
    if ok:
        passed += 1
    else:
        failed += 1
    msg = f"  [{num:02d}] {tag}  {name}"
    if detail:
        msg += f"  {detail}"
    print(msg)

min_text = load_opt('gropt_test_min.opt')
two_text = load_opt('gropt_test_twolens.opt')

TOL = 5e-3

opt = Optics.from_string(min_text)
report(0, "parse min.opt", opt.nglass == 1 and opt.nlens == 1,
       f"nglass={opt.nglass} nlens={opt.nlens}")

tm = opt.transfer(wavelength=0.6)
ref_tm = load_transfer('gropt_test_transfer.out')
d1 = abs(tm['focal_plane'] - ref_tm['focal_plane'])
d2 = abs(tm['effective_focus'] - ref_tm['effective_focus'])
report(1, "transfer", d1 < TOL and d2 < TOL,
       f"focal_diff={d1:.2e} eff_diff={d2:.2e}")

def test_spot(num, name, ref_file, **kw):
    o = Optics.from_string(min_text)
    if 'focus' in kw:
        o.z_focal = kw.pop('focus')
    r = o.spot(**kw)
    ref, ref_center = load_spot(ref_file)
    spots = r['spots']
    if spots.shape != ref.shape:
        report(num, name, False, f"shape {spots.shape} vs {ref.shape}")
        return
    maxd = np.abs(spots - ref).max()
    ok = maxd < TOL
    report(num, name, ok, f"max_diff={maxd:.2e}")

test_spot(2, "basic spot", 'gropt_test_spot.cat',
          wavelength=0.55, aperture_radius=10, nrings=3, pixel_scale=0.01)

test_spot(3, "spot offset", 'gropt_test_spot_offset.cat',
          wavelength=0.55, aperture_radius=10, nrings=3, zstart=2, pixel_scale=0.01)

test_spot(4, "spot angle scalar", 'gropt_test_spot_angle_scalar.cat',
          wavelength=0.55, aperture_radius=10, nrings=3, angle=0.01, pixel_scale=0.01)

test_spot(5, "spot angle vector", 'gropt_test_spot_angle_vector.cat',
          wavelength=0.55, aperture_radius=10, nrings=3, angle_xy=(0.0, 0.01), pixel_scale=0.01)

test_spot(6, "focus override", 'gropt_test_spot_focus120.cat',
          wavelength=0.55, aperture_radius=10, nrings=3, pixel_scale=0.01, focus=120)

test_spot(7, "wavelength 0.80", 'gropt_test_spot_lambda080.cat',
          wavelength=0.80, aperture_radius=10, nrings=3, pixel_scale=0.01)

opt2 = Optics.from_string(min_text)
eps = opt2.to_eps(aperture_radius=10, nrings=3, wavelength=0.55, pixel_scale=0.01)
report(8, "EPS export", len(eps) > 100, f"len={len(eps)}")

scad = opt2.to_openscad()
report(9, "SCAD export", len(scad) > 100, f"len={len(scad)}")

def test_psf(num, name, ref_file, half_size, **kw):
    o = Optics.from_string(min_text)
    psf = o.psf(half_size=half_size, **kw)
    ref = load_psf(ref_file)
    if psf.shape != ref.shape:
        report(num, name, False, f"shape {psf.shape} vs {ref.shape}")
        return
    maxd = np.abs(psf - ref).max()
    ok = maxd < TOL
    report(num, name, ok, f"shape={psf.shape} max_diff={maxd:.2e}")

test_psf(10, "PSF 7x7", 'gropt_test_psf_h3.cat', half_size=3,
         wavelength=0.55, aperture_radius=10, pixel_scale=0.01)

test_psf(11, "PSF 11x11", 'gropt_test_psf_h5.cat', half_size=5,
         wavelength=0.55, aperture_radius=10, pixel_scale=0.01)

opt3 = Optics.from_string(two_text)
report(13, "parse twolens.opt", opt3.nglass == 2 and opt3.nlens == 2,
       f"nglass={opt3.nglass} nlens={opt3.nlens}")

tm3 = opt3.transfer(wavelength=0.6)
ref_tm3 = load_transfer('gropt_test_twolens_transfer.out')
d1 = abs(tm3['focal_plane'] - ref_tm3['focal_plane'])
d2 = abs(tm3['effective_focus'] - ref_tm3['effective_focus'])
report(14, "twolens transfer", d1 < TOL and d2 < TOL,
       f"focal_diff={d1:.2e} eff_diff={d2:.2e}")

o4 = Optics.from_string(two_text)
r4 = o4.spot(wavelength=0.55, aperture_radius=10, nrings=3, pixel_scale=0.01)
ref4, _ = load_spot('gropt_test_twolens_spot.cat')
if r4['spots'].shape == ref4.shape:
    maxd4 = np.abs(r4['spots'] - ref4).max()
    report(15, "twolens spot", maxd4 < TOL, f"max_diff={maxd4:.2e}")
else:
    report(15, "twolens spot", False, f"shape {r4['spots'].shape} vs {ref4.shape}")

print()
print(f"Result: {passed} passed, {failed} failed / {total} total")
sys.exit(1 if failed else 0)
