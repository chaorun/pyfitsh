#!/usr/bin/env bash
set -euo pipefail

P="lfit_test"

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
  "${P}_outlier.dat"
do
    check_table "$f"
done

echo "[01] MCMC basic"
lfit \
  -M \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_mcmc_basic.out"
check_table "${P}_mcmc_basic.out"

echo "[02] MCMC accepted"
lfit \
  -M \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P accepted \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_mcmc_accepted.out"
check_table "${P}_mcmc_accepted.out"

echo "[03] MCMC nonaccepted"
lfit \
  -M \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P nonaccepted \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_mcmc_nonaccepted.out"
check_table "${P}_mcmc_nonaccepted.out"

echo "[04] MCMC gibbs"
lfit \
  -M \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P gibbs \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_mcmc_gibbs.out"
check_table "${P}_mcmc_gibbs.out"

echo "[05] XMMC basic"
lfit \
  -X \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_xmmc_basic.out"
check_table "${P}_xmmc_basic.out"

echo "[06] XMMC skip"
lfit \
  -X \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P skip \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_xmmc_skip.out"
check_table "${P}_xmmc_skip.out"

echo "[07] XMMC adaptive"
lfit \
  -X \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P adaptive \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_xmmc_adaptive.out"
check_table "${P}_xmmc_adaptive.out"

echo "[08] XMMC window"
lfit \
  -X \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P window=10 \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_xmmc_window.out"
check_table "${P}_xmmc_window.out"

echo "[09] EMCE CLLS primary"
lfit \
  -E \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P clls \
  -s 12345 \
  -i 20 \
  "${P}_poly2.dat" \
  -o "${P}_emce_clls.out"
check_table "${P}_emce_clls.out"

echo "[10] EMCE NLLM primary"
lfit \
  -E \
  -v a=1,b=1,c=0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P nllm \
  -s 12345 \
  -i 20 \
  "${P}_poly2.dat" \
  -o "${P}_emce_nllm.out"
check_table "${P}_emce_nllm.out"

echo "[11] EMCE LMND primary"
lfit \
  -E \
  -v a=1,b=1,c=0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -q a=0.001,b=0.001,c=0.001 \
  -P lmnd \
  -s 12345 \
  -i 20 \
  "${P}_poly2.dat" \
  -o "${P}_emce_lmnd.out"
check_table "${P}_emce_lmnd.out"

echo "[12] EMCE DHSX primary"
lfit \
  -E \
  -v a=1:0.5,b=1:0.5,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P dhsx \
  -s 12345 \
  -i 20 \
  "${P}_poly2.dat" \
  -o "${P}_emce_dhsx.out"
check_table "${P}_emce_dhsx.out"

echo "[13] EMCE perturbations"
lfit \
  -E \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P clls \
  --perturbations 0.01 \
  -s 12345 \
  -i 20 \
  "${P}_poly2.dat" \
  -o "${P}_emce_perturb.out"
check_table "${P}_emce_perturb.out"

echo "[14] FIMA basic"
lfit \
  -A \
  -v a=1,b=2,c=0.5 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  "${P}_poly2.dat" \
  -o "${P}_fima_basic.out"
check_file "${P}_fima_basic.out"

echo "[15] FIMA with mock Gaussian output"
lfit \
  -A \
  -v a=1,b=2,c=0.5 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -s 12345 \
  -i 50 \
  "${P}_poly2.dat" \
  -o "${P}_fima_mock.out"
check_table "${P}_fima_mock.out"

echo "[16] FIMA derived variable"
lfit \
  -A \
  -v a=1,b=2,c=0.5 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -g "sum=a+b+c" \
  "${P}_poly2.dat" \
  -o "${P}_fima_derived.out"
check_file "${P}_fima_derived.out"

echo "[17] MCMC weighted"
lfit \
  -M \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2,e:3 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -e "e" \
  -P accepted \
  -s 12345 \
  -i 50 \
  "${P}_weighted.dat" \
  -o "${P}_mcmc_weighted.out"
check_table "${P}_mcmc_weighted.out"

echo "[18] XMMC weighted"
lfit \
  -X \
  -v a=1:0.2,b=1:0.2,c=0.1:0.1 \
  -c x:1,y:2,e:3 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -e "e" \
  -P skip \
  -s 12345 \
  -i 50 \
  "${P}_weighted.dat" \
  -o "${P}_xmmc_weighted.out"
check_table "${P}_xmmc_weighted.out"

echo "[19] EMCE with outlier data"
lfit \
  -E \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P clls \
  -s 12345 \
  -i 20 \
  "${P}_outlier.dat" \
  -o "${P}_emce_outlier.out"
check_table "${P}_emce_outlier.out"

echo "[20] summary"
python - <<'PY'
from pathlib import Path
files = sorted(Path(".").glob("lfit_test_*mc*.out")) + \
        sorted(Path(".").glob("lfit_test_*mcmc*.out")) + \
        sorted(Path(".").glob("lfit_test_*xmmc*.out")) + \
        sorted(Path(".").glob("lfit_test_*emce*.out")) + \
        sorted(Path(".").glob("lfit_test_*fima*.out"))
seen = set()
uniq = []
for f in files:
    if f not in seen:
        uniq.append(f)
        seen.add(f)
print("Generated Monte-Carlo/Fisher files:", len(uniq))
for f in uniq:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL LFIT MONTECARLO TESTS PASSED"