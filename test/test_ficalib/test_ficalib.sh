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

for f in \
  "${P}_sci_a.fits" \
  "${P}_sci_b.fits" \
  "${P}_bias.fits" \
  "${P}_dark.fits" \
  "${P}_flat.fits" \
  "${P}_mask.fits"
do
    check_fits "$f" 2560 2560
done

echo "[01] ficalib copy only"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_copy.fits" \
  -b -32
check_fits "${P}_copy.fits" 2560 2560

echo "[02] ficalib bias"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_biascorr.fits" \
  -B "${P}_bias.fits" \
  -b -32
check_fits "${P}_biascorr.fits" 2560 2560

echo "[03] ficalib bias flat"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_bias_flat.fits" \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  -b -32
check_fits "${P}_bias_flat.fits" 2560 2560

echo "[04] ficalib bias dark flat no exptime correction"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_bdf_noexp.fits" \
  -B "${P}_bias.fits" \
  -D "${P}_dark.fits" \
  -F "${P}_flat.fits" \
  --no-exptime-correction \
  -b -32
check_fits "${P}_bdf_noexp.fits" 2560 2560

echo "[05] ficalib bias dark flat exptime correction"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_bdf_exp.fits" \
  -B "${P}_bias.fits" \
  -D "${P}_dark.fits" \
  -F "${P}_flat.fits" \
  --exptime-correction \
  -b -32
check_fits "${P}_bdf_exp.fits" 2560 2560

echo "[06] ficalib multi input output"
ficalib \
  -i "${P}_sci_a.fits" "${P}_sci_b.fits" \
  -o "${P}_multi_a.fits" "${P}_multi_b.fits" \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  -b -32
check_fits "${P}_multi_a.fits" 2560 2560
check_fits "${P}_multi_b.fits" 2560 2560

echo "[07] ficalib rewrite output name"
rm -f "${P}_sci_a.calib.fits" "${P}_sci_b.calib.fits"
ficalib \
  -i "${P}_sci_a.fits" "${P}_sci_b.fits" \
  -r '*.fits|*.calib.fits' \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  -b -32
check_fits "${P}_sci_a.calib.fits" 2560 2560
check_fits "${P}_sci_b.calib.fits" 2560 2560

echo "[08] ficalib trim"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_trim.fits" \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  --image 100:100:1900:1900 \
  --trim \
  -b -32
check_fits "${P}_trim.fits" 1801 1801

echo "[09] ficalib post multiply"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_postmul.fits" \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  --post-multiply 2.0 \
  -b -32
check_fits "${P}_postmul.fits" 2560 2560

echo "[10] ficalib post scale"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_postscale.fits" \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  --post-scale 10000 \
  -b -32
check_fits "${P}_postscale.fits" 2560 2560

echo "[11] ficalib saturation gain"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_sat_gain.fits" \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  --saturation 50000 \
  -g 2 \
  -b -32
check_fits "${P}_sat_gain.fits" 2560 2560

echo "[12] ficalib input mask"
ficalib \
  -i "${P}_sci_a.fits" \
  -o "${P}_maskout.fits" \
  -M "${P}_mask.fits" \
  -B "${P}_bias.fits" \
  -F "${P}_flat.fits" \
  -b -32
check_fits "${P}_maskout.fits" 2560 2560

echo "[13] summary"
python - <<'PY'
from pathlib import Path
files = sorted(Path(".").glob("ficalib_test*"))
print("Generated files:", len(files))
for f in files:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL FICALIB TESTS PASSED"