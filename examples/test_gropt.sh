#!/usr/bin/env bash
set -euo pipefail

P="gropt_test"

check_file() {
    local f="$1"
    test -s "$f" || { echo "ERROR: missing or empty file: $f"; exit 1; }
}

check_table() {
    local f="$1"
    check_file "$f"
    local n
    n=$(awk 'NF>0 && $1 !~ /^#/ {n++} END{print n+0}' "$f")
    if [ "$n" -le 0 ]; then
        echo "ERROR: no data rows in table: $f"
        exit 1
    fi
    echo "CHECK TABLE OK: $f rows=$n"
}

echo "[00] make minimal optical descriptor"
cat > "${P}_min.opt" <<'EOF'
# minimal gropt optical descriptor
# units: mm
glass BK7 1.5168
lens 0 5 25 +1/100 -1/100 BK7 - -
focal 100
EOF
check_file "${P}_min.opt"

echo "[01] transfer matrix analysis"
gropt \
  -i "${P}_min.opt" \
  -t \
  > "${P}_transfer.out"
check_file "${P}_transfer.out"
cat "${P}_transfer.out"

echo "[02] basic spot diagram"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -l 0.55 \
  -f 100 \
  -x 0.01 \
  -o "${P}_spot.cat"
check_table "${P}_spot.cat"

echo "[03] spot diagram with aperture offset"
gropt \
  -i "${P}_min.opt" \
  -s 10,3,2 \
  -l 0.55 \
  -f 100 \
  -x 0.01 \
  -o "${P}_spot_offset.cat"
check_table "${P}_spot_offset.cat"

echo "[04] spot diagram with scalar angle"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -a 0.01 \
  -l 0.55 \
  -f 100 \
  -x 0.01 \
  -o "${P}_spot_angle_scalar.cat"
check_table "${P}_spot_angle_scalar.cat"

echo "[05] spot diagram with normal vector angle"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -a 0.0,0.01 \
  -l 0.55 \
  -f 100 \
  -x 0.01 \
  -o "${P}_spot_angle_vector.cat"
check_table "${P}_spot_angle_vector.cat"

echo "[06] focus override"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -l 0.55 \
  -f 120 \
  -x 0.01 \
  -o "${P}_spot_focus120.cat"
check_table "${P}_spot_focus120.cat"

echo "[07] wavelength override"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -l 0.80 \
  -f 100 \
  -x 0.01 \
  -o "${P}_spot_lambda080.cat"
check_table "${P}_spot_lambda080.cat"

echo "[08] EPS geometry export"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -l 0.55 \
  -f 100 \
  -x 0.01 \
  -e "${P}_geom.eps"
check_file "${P}_geom.eps"

echo "[09] SCAD geometry export"
gropt \
  -i "${P}_min.opt" \
  -d "${P}_geom.scad"
check_file "${P}_geom.scad"

echo "[10] PSF 7x7"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -l 0.55 \
  -f 100 \
  -x 0.01 \
  -z 3 \
  -p "${P}_psf_h3.cat"
check_table "${P}_psf_h3.cat"

echo "[11] PSF 11x11"
gropt \
  -i "${P}_min.opt" \
  -s 10,3 \
  -l 0.55 \
  -f 100 \
  -x 0.01 \
  -z 5 \
  -p "${P}_psf_h5.cat"
check_table "${P}_psf_h5.cat"

echo "[12] stdin input transfer"
cat "${P}_min.opt" | gropt \
  -i - \
  -t \
  > "${P}_stdin_transfer.out"
check_file "${P}_stdin_transfer.out"

echo "[13] make two-lens descriptor"
cat > "${P}_twolens.opt" <<'EOF'
# two-lens gropt optical descriptor
# units: mm
glass BK7 1.5168
glass F2 1.6200
lens 0 5 25 +1/100 -1/100 BK7 - -
lens 20 4 22 +1/80 -1/120 F2 - -
focal 150
EOF
check_file "${P}_twolens.opt"

echo "[14] two-lens transfer"
gropt \
  -i "${P}_twolens.opt" \
  -t \
  > "${P}_twolens_transfer.out"
check_file "${P}_twolens_transfer.out"
cat "${P}_twolens_transfer.out"

echo "[15] two-lens spot"
gropt \
  -i "${P}_twolens.opt" \
  -s 10,3 \
  -l 0.55 \
  -f 150 \
  -x 0.01 \
  -o "${P}_twolens_spot.cat"
check_table "${P}_twolens_spot.cat"

echo "[16] summary"
python - <<'PY'
from pathlib import Path
files = sorted(Path(".").glob("gropt_test_*"))
print("Generated gropt files:", len(files))
for f in files:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL GROPT TESTS PASSED"