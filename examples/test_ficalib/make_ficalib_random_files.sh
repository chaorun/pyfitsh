#!/usr/bin/env bash
set -euo pipefail

P="ficalib_test"

check_file() {
    local f="$1"
    test -s "$f" || { echo "ERROR: missing or empty file: $f"; exit 1; }
}

check_fits() {
    local f="$1"
    local sx="${2:-}"
    local sy="${3:-}"
    check_file "$f"
    python - "$f" "$sx" "$sy" <<'PY'
import sys
from astropy.io import fits

fn = sys.argv[1]
sx = sys.argv[2]
sy = sys.argv[3]

data = fits.getdata(fn)
if data is None:
    raise SystemExit(f"ERROR: no FITS data: {fn}")

print(f"CHECK FITS OK: {fn} shape={data.shape} dtype={data.dtype}")

if sx and sy:
    exp = (int(sy), int(sx))
    if data.shape != exp:
        raise SystemExit(f"ERROR: shape mismatch for {fn}: got {data.shape}, expected {exp}")
PY
}

echo "[01] make science A"
firandom \
  --seed 10101 \
  --size 2560,2560 \
  --sky 1000 \
  --sky-noise 3 \
  --list "40[x=r(0.05,0.95),y=r(0.05,0.95),f=5,e=0.03,p=r(0,180),i=r(1000,3000)]" \
  --fep \
  --output-list "${P}_sci_a.lst" \
  -o "${P}_sci_a.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_sci_a.fits" 2560 2560

echo "[02] make science B"
firandom \
  --seed 20202 \
  --size 2560,2560 \
  --sky 1100 \
  --sky-noise 4 \
  --list "35[x=r(0.05,0.95),y=r(0.05,0.95),f=6,e=0.05,p=r(0,180),i=r(800,2500)]" \
  --fep \
  --output-list "${P}_sci_b.lst" \
  -o "${P}_sci_b.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_sci_b.fits" 2560 2560

echo "[03] make bias"
fiarith \
  -s 2560,2560 \
  "100.0" \
  -o "${P}_bias.fits"
check_fits "${P}_bias.fits" 2560 2560

echo "[04] make dark"
fiarith \
  -s 2560,2560 \
  "10.0" \
  -o "${P}_dark.fits"
check_fits "${P}_dark.fits" 2560 2560

echo "[05] make flat"
fiarith \
  -s 2560,2560 \
  "1.0" \
  -o "${P}_flat.fits"
check_fits "${P}_flat.fits" 2560 2560

echo "[06] make mask"
fiarith \
  -s 2560,2560 \
  "0.0" \
  -o "${P}_mask.fits"
check_fits "${P}_mask.fits" 2560 2560

echo "ALL RANDOM FICALIB INPUT FILES CREATED"