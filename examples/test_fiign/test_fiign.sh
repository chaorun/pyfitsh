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

for f in "${P}_base.fits" "${P}_signed.fits" "${P}_satlevel.fits"; do
    check_fits "$f"
done

echo "[01] ignore negative + apply mask"
fiign \
  -i "${P}_signed.fits" \
  -o "${P}_ignore_negative_apply.fits" \
  --ignore-negative \
  --apply-mask \
  --mask-value -999
check_fits "${P}_ignore_negative_apply.fits"

echo "[02] ignore zero + apply mask"
fiign \
  -i "${P}_signed.fits" \
  -o "${P}_ignore_zero_apply.fits" \
  --ignore-zero \
  --apply-mask \
  --mask-value -888
check_fits "${P}_ignore_zero_apply.fits"

echo "[03] ignore nonpositive + apply mask"
fiign \
  -i "${P}_signed.fits" \
  -o "${P}_ignore_nonpositive_apply.fits" \
  --ignore-nonpositive \
  --apply-mask \
  --mask-value -777
check_fits "${P}_ignore_nonpositive_apply.fits"

echo "[04] saturation scalar + apply mask"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_saturation_scalar_apply.fits" \
  -s 1200 \
  --apply-mask \
  --mask-value -666
check_fits "${P}_saturation_scalar_apply.fits"

echo "[05] saturation image + apply mask"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_saturation_image_apply.fits" \
  -S "${P}_satlevel.fits" \
  --apply-mask \
  --mask-value -555
check_fits "${P}_saturation_image_apply.fits"

echo "[06] mask block rectangle"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_block_rect_apply.fits" \
  -q block:hot:100,100:300,300 \
  --apply-mask \
  --mask-value -444
check_fits "${P}_block_rect_apply.fits"

echo "[07] mask block circle"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_block_circle_apply.fits" \
  -q circle:cosmic:1200,1200:100 \
  --apply-mask \
  --mask-value -333
check_fits "${P}_block_circle_apply.fits"

echo "[08] mask block pixel"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_block_pixel_apply.fits" \
  -q pixel:fault:500,500 \
  --apply-mask \
  --mask-value -222
check_fits "${P}_block_pixel_apply.fits"

echo "[09] mask block line"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_block_line_apply.fits" \
  -q line:cosmic:100,200:800,900 \
  --apply-mask \
  --mask-value -111
check_fits "${P}_block_line_apply.fits"

echo "[10] expand mask"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_expand_apply.fits" \
  -q pixel:hot:1000,1000 \
  --expand 5 \
  --apply-mask \
  --mask-value -123
check_fits "${P}_expand_apply.fits"

echo "[11] create external mask carrier"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_external_mask.fits" \
  -q block:fault:400,400:600,600
check_fits "${P}_external_mask.fits"

echo "[12] input mask co-add + apply mask"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_input_mask_apply.fits" \
  -M "${P}_external_mask.fits" \
  --apply-mask \
  --mask-value -321
check_fits "${P}_input_mask_apply.fits"

echo "[13] convert hot to cosmic"
fiign \
  -i "${P}_base.fits" \
  -o "${P}_convert_hot_cosmic_apply.fits" \
  -q pixel:hot:1300,1300 \
  --convert hot:hot:none:cosmic \
  --apply-mask \
  --mask-value -231
check_fits "${P}_convert_hot_cosmic_apply.fits"

echo "[14] ignore input mask"
fiign \
  -i "${P}_external_mask.fits" \
  -o "${P}_ignore_mask_apply.fits" \
  --ignore-mask \
  --apply-mask \
  --mask-value -135
check_fits "${P}_ignore_mask_apply.fits"

echo "[15] summary"
python - <<'PY'
from pathlib import Path
files = sorted(Path(".").glob("fiign_test*"))
print("Generated files:", len(files))
for f in files:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL FIIGN TESTS PASSED"