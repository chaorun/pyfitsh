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
  "${P}_linear.dat" \
  "${P}_poly2.dat" \
  "${P}_weighted.dat" \
  "${P}_outlier.dat" \
  "${P}_pairs.dat"
do
    check_table "$f"
done

echo "[01] CLLS linear fit"
lfit \
  -L \
  -v a,b \
  -c x:1,y:2 \
  -f "a+b*x" \
  -y "y" \
  "${P}_linear.dat" \
  -o "${P}_clls_linear.out"
check_file "${P}_clls_linear.out"

echo "[02] CLLS linear fit with variables output"
lfit \
  -L \
  -v a,b \
  -c x:1,y:2 \
  -f "a+b*x" \
  -y "y" \
  "${P}_linear.dat" \
  -o "${P}_clls_linear_main.out" \
  -l "${P}_clls_linear.vars" \
  -p "${P}_clls_linear.expr"
check_file "${P}_clls_linear_main.out"
check_file "${P}_clls_linear.vars"
check_file "${P}_clls_linear.expr"

echo "[03] CLLS poly2 fit"
lfit \
  -L \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  "${P}_poly2.dat" \
  -o "${P}_clls_poly2.out"
check_file "${P}_clls_poly2.out"

echo "[04] CLLS poly2 with error columns"
lfit \
  -L \
  -v a,b,c \
  -c x:1,y:2,e:3 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -e "e" \
  "${P}_weighted.dat" \
  -o "${P}_clls_weighted_error.out" \
  --errors
check_file "${P}_clls_weighted_error.out"

echo "[05] CLLS poly2 with weight expression"
lfit \
  -L \
  -v a,b,c \
  -c x:1,y:2,e:3 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -w "1/e" \
  "${P}_weighted.dat" \
  -o "${P}_clls_weighted_weight.out"
check_file "${P}_clls_weighted_weight.out"

echo "[06] CLLS outlier rejection"
lfit \
  -L \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -r 3 \
  -n 3 \
  "${P}_outlier.dat" \
  -o "${P}_clls_reject.out" \
  -u "${P}_clls_reject.fitted" \
  -j "${P}_clls_reject.rejected" \
  -a "${P}_clls_reject.all"
check_file "${P}_clls_reject.out"
check_file "${P}_clls_reject.fitted"
check_file "${P}_clls_reject.all"

echo "[07] CLLS outlier rejection with delta"
lfit \
  -L \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -r 3 \
  -n 3 \
  --delta \
  "${P}_outlier.dat" \
  -o "${P}_clls_reject_delta.out" \
  -u "${P}_clls_reject_delta.fitted" \
  -j "${P}_clls_reject_delta.rejected" \
  -a "${P}_clls_reject_delta.all"
check_file "${P}_clls_reject_delta.out"
check_file "${P}_clls_reject_delta.fitted"
check_file "${P}_clls_reject_delta.all"

echo "[08] CLLS output format"
lfit \
  -L \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -F a=%16.10g,b=%16.10g,c=%16.10g \
  "${P}_poly2.dat" \
  -o "${P}_clls_poly2_format.out"
check_file "${P}_clls_poly2_format.out"

echo "[09] derived variable"
lfit \
  -L \
  -v a,b,c \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -g "sum=a+b+c" \
  "${P}_poly2.dat" \
  -o "${P}_clls_poly2_derived.out"
check_file "${P}_clls_poly2_derived.out"

echo "[10] constraint via := fixed c"
lfit \
  -L \
  -v a,b,c:=0.5 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  "${P}_poly2.dat" \
  -o "${P}_clls_poly2_fixed_c.out"
check_file "${P}_clls_poly2_fixed_c.out"

echo "[11] nonlinear LM analytic"
lfit \
  -N \
  -v a=1,b=1,c=0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -P iterations=50 \
  "${P}_poly2.dat" \
  -o "${P}_nllm_poly2.out"
check_file "${P}_nllm_poly2.out"

echo "[12] nonlinear LM numerical derivatives"
lfit \
  -U \
  -v a=1,b=1,c=0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  -q a=0.001,b=0.001,c=0.001 \
  -P iterations=50 \
  "${P}_poly2.dat" \
  -o "${P}_lmnd_poly2.out"
check_file "${P}_lmnd_poly2.out"

echo "[13] downhill simplex"
lfit \
  -D \
  -v a=1:0.5,b=1:0.5,c=0.1:0.1 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  "${P}_poly2.dat" \
  -o "${P}_dhsx_poly2.out"
check_file "${P}_dhsx_poly2.out"

echo "[14] chi2 grid"
lfit \
  -K \
  -v a=[0:1:3],b=[1:1:4],c=0.5 \
  -c x:1,y:2 \
  -f "a+b*x+c*x*x" \
  -y "y" \
  "${P}_poly2.dat" \
  -o "${P}_chi2_grid.out"
check_table "${P}_chi2_grid.out"

echo "[15] function evaluation default output"
lfit \
  -c x:1,y:2 \
  -f "2+3*x" \
  "${P}_linear.dat" \
  -o "${P}_eval_linear.out"
check_table "${P}_eval_linear.out"

echo "[16] function evaluation replace columns"
lfit \
  -c x:1,y:2 \
  -f "2+3*x,y-(2+3*x)" \
  -z 1,2 \
  "${P}_linear.dat" \
  -o "${P}_eval_replace_cols.out"
check_table "${P}_eval_replace_cols.out"

echo "[17] macro definition"
lfit \
  -L \
  -v a,b \
  -c x:1,y:2 \
  -x "line(u,v,t)=u+v*t" \
  -f "line(a,b,x)" \
  -y "y" \
  "${P}_linear.dat" \
  -o "${P}_macro_linear.out"
check_file "${P}_macro_linear.out"

echo "[18] paired catalog dx fit"
lfit \
  -L \
  -v a,b,c \
  -c xr:1,yr:2,xo:3,yo:4,mr:5,mo:6 \
  -f "a+b*xr+c*yr" \
  -y "xo-xr" \
  "${P}_pairs.dat" \
  -o "${P}_pairs_dx_linear.out"
check_file "${P}_pairs_dx_linear.out"

echo "[19] paired catalog dy fit"
lfit \
  -L \
  -v a,b,c \
  -c xr:1,yr:2,xo:3,yo:4,mr:5,mo:6 \
  -f "a+b*xr+c*yr" \
  -y "yo-yr" \
  "${P}_pairs.dat" \
  -o "${P}_pairs_dy_linear.out"
check_file "${P}_pairs_dy_linear.out"

echo "[20] summary"
python - <<'PY'
from pathlib import Path
files = sorted(Path(".").glob("lfit_test*"))
print("Generated files:", len(files))
for f in files:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL LFIT TESTS PASSED"