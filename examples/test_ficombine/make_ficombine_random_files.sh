#!/usr/bin/env bash
set -euo pipefail

P="ficombine_test"

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

echo "[01] make image 1"
firandom \
  --seed 41001 \
  --size 2560,2560 \
  --sky 1000 \
  --sky-noise 3 \
  --list "30[x=r(0.05,0.95),y=r(0.05,0.95),f=5,e=0.03,p=r(0,180),i=r(1000,3000)]" \
  --fep \
  --output-list "${P}_img1.lst" \
  -o "${P}_img1.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_img1.fits" 2560 2560

echo "[02] make image 2"
firandom \
  --seed 41002 \
  --size 2560,2560 \
  --sky 1010 \
  --sky-noise 3 \
  --list "30[x=r(0.05,0.95),y=r(0.05,0.95),f=5,e=0.04,p=r(0,180),i=r(900,2800)]" \
  --fep \
  --output-list "${P}_img2.lst" \
  -o "${P}_img2.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_img2.fits" 2560 2560

echo "[03] make image 3"
firandom \
  --seed 41003 \
  --size 2560,2560 \
  --sky 990 \
  --sky-noise 4 \
  --list "30[x=r(0.05,0.95),y=r(0.05,0.95),f=6,e=0.02,p=r(0,180),i=r(1200,3500)]" \
  --fep \
  --output-list "${P}_img3.lst" \
  -o "${P}_img3.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_img3.fits" 2560 2560

echo "[04] make image 4"
firandom \
  --seed 41004 \
  --size 2560,2560 \
  --sky 1020 \
  --sky-noise 5 \
  --list "30[x=r(0.05,0.95),y=r(0.05,0.95),f=4,e=0.05,p=r(0,180),i=r(800,2600)]" \
  --fep \
  --output-list "${P}_img4.lst" \
  -o "${P}_img4.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_img4.fits" 2560 2560

echo "[05] make image 5"
firandom \
  --seed 41005 \
  --size 2560,2560 \
  --sky 980 \
  --sky-noise 4 \
  --list "30[x=r(0.05,0.95),y=r(0.05,0.95),f=7,e=0.01,p=r(0,180),i=r(1000,3200)]" \
  --fep \
  --output-list "${P}_img5.lst" \
  -o "${P}_img5.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_fits "${P}_img5.fits" 2560 2560

echo "[06] make negative image"
fiarith \
  -s 2560,2560 \
  "['${P}_img1.fits'](a-1200)" \
  -o "${P}_negative.fits"
check_fits "${P}_negative.fits" 2560 2560

echo "[07] make zero image"
fiarith \
  -s 2560,2560 \
  "0.0" \
  -o "${P}_zero.fits"
check_fits "${P}_zero.fits" 2560 2560

cat > "${P}_inputs.lst" <<EOF
${P}_img1.fits
${P}_img2.fits
${P}_img3.fits
${P}_img4.fits
${P}_img5.fits
EOF

echo "ALL FICOMBINE INPUT FILES CREATED"