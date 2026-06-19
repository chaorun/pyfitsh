#!/usr/bin/env bash
set -euo pipefail

P="lfit_test"

MCMC_N="${MCMC_N:-20000}"
XMMC_N="${XMMC_N:-10000}"
EMCE_N="${EMCE_N:-2000}"
FIMA_N="${FIMA_N:-20000}"

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

for f in \
  "${P}_poly2.dat" \
  "${P}_weighted.dat" \
  "${P}_outlier.dat" \
  "${P}_pairs.dat"
do
    check_table "$f"
done

echo "[01] long MCMC accepted"
time lfit \
  -M \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P accepted \
  -s 12345 \
  -i "$MCMC_N" \
  "${P}_poly2.dat" \
  -o "${P}_long_mcmc_accepted.out"
check_table "${P}_long_mcmc_accepted.out"

echo "[02] long MCMC gibbs weighted"
time lfit \
  -M \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2,e:3 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -e "e" \
  -P gibbs \
  -s 12345 \
  -i "$MCMC_N" \
  "${P}_weighted.dat" \
  -o "${P}_long_mcmc_gibbs_weighted.out"
check_table "${P}_long_mcmc_gibbs_weighted.out"

echo "[03] long XMMC adaptive"
time lfit \
  -X \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P adaptive,window=20 \
  -s 12345 \
  -i "$XMMC_N" \
  "${P}_poly2.dat" \
  -o "${P}_long_xmmc_adaptive.out"
check_table "${P}_long_xmmc_adaptive.out"

echo "[04] long XMMC skip weighted"
time lfit \
  -X \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2,e:3 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -e "e" \
  -P skip \
  -s 12345 \
  -i "$XMMC_N" \
  "${P}_weighted.dat" \
  -o "${P}_long_xmmc_skip_weighted.out"
check_table "${P}_long_xmmc_skip_weighted.out"

echo "[05] long EMCE CLLS"
time lfit \
  -E \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P clls \
  --perturbations 0.01 \
  -s 12345 \
  -i "$EMCE_N" \
  "${P}_poly2.dat" \
  -o "${P}_long_emce_clls.out"
check_table "${P}_long_emce_clls.out"

echo "[06] long EMCE DHSX"
time lfit \
  -E \
  -v a=1:0.5,b=1:0.5,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P dhsx \
  --perturbations 0.01 \
  -s 12345 \
  -i "$EMCE_N" \
  "${P}_poly2.dat" \
  -o "${P}_long_emce_dhsx.out"
check_table "${P}_long_emce_dhsx.out"

echo "[07] long FIMA mock"
time lfit \
  -A \
  -v a=1,b=2,c=0.5 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -g "sum=a+b+c" \
  -s 12345 \
  -i "$FIMA_N" \
  "${P}_poly2.dat" \
  -o "${P}_long_fima_mock.out"
check_table "${P}_long_fima_mock.out"

echo "[08] long catalog-pair MCMC dx"
time lfit \
  -M \
  -v a=0:1,b=0:0.01,c=0:0.01 \
  -c xr:1,yr:2,xo:3,yo:4,mr:5,mo:6 \
  -f "a+b*xr+c*yr" \
  -y "xo-xr" \
  -P accepted \
  -s 12345 \
  -i "$MCMC_N" \
  "${P}_pairs.dat" \
  -o "${P}_long_pairs_dx_mcmc.out"
check_table "${P}_long_pairs_dx_mcmc.out"

echo "[09] long catalog-pair MCMC dy"
time lfit \
  -M \
  -v a=0:1,b=0:0.01,c=0:0.01 \
  -c xr:1,yr:2,xo:3,yo:4,mr:5,mo:6 \
  -f "a+b*xr+c*yr" \
  -y "yo-yr" \
  -P accepted \
  -s 12345 \
  -i "$MCMC_N" \
  "${P}_pairs.dat" \
  -o "${P}_long_pairs_dy_mcmc.out"
check_table "${P}_long_pairs_dy_mcmc.out"

echo "[10] summary"
python - <<'PY'
from pathlib import Path

files = sorted(Path(".").glob("lfit_test_long_*.out"))
print("Generated long lfit files:", len(files))
for f in files:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL LONG LFIT TESTS PASSED"