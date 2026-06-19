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

for f in \
  "${P}_img1.fits" \
  "${P}_img2.fits" \
  "${P}_img3.fits" \
  "${P}_img4.fits" \
  "${P}_img5.fits" \
  "${P}_negative.fits" \
  "${P}_zero.fits"
do
    check_fits "$f" 2560 2560
done

IMGS="${P}_img1.fits ${P}_img2.fits ${P}_img3.fits ${P}_img4.fits ${P}_img5.fits"

echo "[01] mean"
ficombine $IMGS \
  -m mean \
  -o "${P}_mean.fits" \
  -b -32
check_fits "${P}_mean.fits" 2560 2560

echo "[02] median"
ficombine $IMGS \
  -m median \
  -o "${P}_median.fits" \
  -b -32
check_fits "${P}_median.fits" 2560 2560

echo "[03] min"
ficombine $IMGS \
  -m min \
  -o "${P}_min.fits" \
  -b -32
check_fits "${P}_min.fits" 2560 2560

echo "[04] max"
ficombine $IMGS \
  -m max \
  -o "${P}_max.fits" \
  -b -32
check_fits "${P}_max.fits" 2560 2560

echo "[05] sum"
ficombine $IMGS \
  -m sum \
  -o "${P}_sum.fits" \
  -b -32
check_fits "${P}_sum.fits" 2560 2560

echo "[06] squaresum"
ficombine $IMGS \
  -m squaresum \
  -o "${P}_squaresum.fits" \
  -b -32
check_fits "${P}_squaresum.fits" 2560 2560

echo "[07] scatter"
ficombine $IMGS \
  -m scatter \
  -o "${P}_scatter.fits" \
  -b -32
check_fits "${P}_scatter.fits" 2560 2560

echo "[08] stddev alias"
ficombine $IMGS \
  -m stddev \
  -o "${P}_stddev.fits" \
  -b -32
check_fits "${P}_stddev.fits" 2560 2560

echo "[09] rejmed"
ficombine $IMGS \
  -m rejmed,iterations=2,sigma=3 \
  -o "${P}_rejmed.fits" \
  -b -32
check_fits "${P}_rejmed.fits" 2560 2560

echo "[10] rejection alias"
ficombine $IMGS \
  -m rejection,iterations=2,lower=3,upper=3 \
  -o "${P}_rejection.fits" \
  -b -32
check_fits "${P}_rejection.fits" 2560 2560

echo "[11] rejmean"
ficombine $IMGS \
  -m rejmean,iterations=2,sigma=3 \
  -o "${P}_rejmean.fits" \
  -b -32
check_fits "${P}_rejmean.fits" 2560 2560

echo "[12] truncated discard"
ficombine $IMGS \
  -m truncated,discard=1 \
  -o "${P}_truncated_discard.fits" \
  -b -32
check_fits "${P}_truncated_discard.fits" 2560 2560

echo "[13] truncated lowest highest"
ficombine $IMGS \
  -m truncated,lowest=1,highest=1 \
  -o "${P}_truncated_low_high.fits" \
  -b -32
check_fits "${P}_truncated_low_high.fits" 2560 2560

echo "[14] winsorized discard"
ficombine $IMGS \
  -m winsorized,discard=1 \
  -o "${P}_winsorized_discard.fits" \
  -b -32
check_fits "${P}_winsorized_discard.fits" 2560 2560

echo "[15] ignore negative by option"
ficombine \
  "${P}_negative.fits" \
  "${P}_img1.fits" \
  "${P}_img2.fits" \
  -n \
  -m mean \
  -o "${P}_ignore_negative_option.fits" \
  -b -32
check_fits "${P}_ignore_negative_option.fits" 2560 2560

echo "[16] ignorenegative in mode"
ficombine \
  "${P}_negative.fits" \
  "${P}_img1.fits" \
  "${P}_img2.fits" \
  -m mean,ignorenegative \
  -o "${P}_ignore_negative_mode.fits" \
  -b -32
check_fits "${P}_ignore_negative_mode.fits" 2560 2560

echo "[17] ignorezero in mode"
ficombine \
  "${P}_zero.fits" \
  "${P}_img1.fits" \
  "${P}_img2.fits" \
  -m mean,ignorezero \
  -o "${P}_ignorezero_mode.fits" \
  -b -32
check_fits "${P}_ignorezero_mode.fits" 2560 2560

echo "[18] no-history"
ficombine $IMGS \
  -m mean \
  --no-history \
  -o "${P}_mean_nohistory.fits" \
  -b -32
check_fits "${P}_mean_nohistory.fits" 2560 2560

echo "[19] history"
ficombine $IMGS \
  -m mean \
  --history \
  -o "${P}_mean_history.fits" \
  -b -32
check_fits "${P}_mean_history.fits" 2560 2560

echo "[20] bitpix int16"
ficombine $IMGS \
  -m mean \
  -o "${P}_mean_i16.fits" \
  -b 16
check_fits "${P}_mean_i16.fits" 2560 2560

echo "[21] logical or"
ficombine $IMGS \
  -m or \
  -o "${P}_mask_or.fits" \
  -b -32
check_fits "${P}_mask_or.fits" 2560 2560

echo "[22] logical and"
ficombine $IMGS \
  -m and \
  -o "${P}_mask_and.fits" \
  -b -32
check_fits "${P}_mask_and.fits" 2560 2560

echo "[23] logical-or option"
ficombine $IMGS \
  --logical-or \
  -m mean \
  -o "${P}_logical_or_option.fits" \
  -b -32
check_fits "${P}_logical_or_option.fits" 2560 2560

echo "[24] logical-and option"
ficombine $IMGS \
  --logical-and \
  -m mean \
  -o "${P}_logical_and_option.fits" \
  -b -32
check_fits "${P}_logical_and_option.fits" 2560 2560

echo "[25] summary"
python - <<'PY'
from pathlib import Path
files = sorted(Path(".").glob("ficombine_test*"))
print("Generated files:", len(files))
for f in files:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL FICOMBINE TESTS PASSED"