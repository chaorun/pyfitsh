#!/usr/bin/env bash
set -euo pipefail

P="fiign_test"

check_file() {
    local f="$1"
    test -s "$f" || { echo "ERROR: missing or empty file: $f"; exit 1; }
}

check_fits() {
    local f="$1"
    check_file "$f"
    python - "$f" <<'PY'
import sys
from astropy.io import fits
fn = sys.argv[1]
data = fits.getdata(fn)
if data is None:
    raise SystemExit(f"ERROR: no FITS data: {fn}")
print(f"CHECK FITS OK: {fn} shape={data.shape} dtype={data.dtype}")
PY
}

echo "[01] make base positive image"
firandom \
  --seed 31001 \
  --size 2560,2560 \
  --sky 1000 \
  --sky-noise 3 \
  --list "30[x=r(0.05,0.95),y=r(0.05,0.95),f=5,e=0.03,p=r(0,180),i=r(1000,3000)]" \
  --fep \
  --output-list "${P}_base.lst" \
  -o "${P}_base.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_base.fits"

echo "[02] make image with negative/zero values"
fiarith \
  -s 2560,2560 \
  "['${P}_base.fits'](a-1000)" \
  -o "${P}_signed.fits"
check_fits "${P}_signed.fits"

echo "[03] make saturation-level image"
fiarith \
  -s 2560,2560 \
  "1200.0" \
  -o "${P}_satlevel.fits"
check_fits "${P}_satlevel.fits"

echo "ALL FIIGN INPUT FILES CREATED"