#!/usr/bin/env bash
set -euo pipefail

REF_CAT="ref_L1_GRB260604C-Ic-180s-20260606_210318.cat"
OBS_CAT="obs_L1_GRB260604C-Ic-180s-20260607_023443.cat"
P="lfit_test"

test -s "$REF_CAT"
test -s "$OBS_CAT"

echo "[01] make linear data"
awk 'NF>=4 {x=$2; y=2.0+3.0*x; print x,y}' "$REF_CAT" > "${P}_linear.dat"

echo "[02] make polynomial data"
awk 'NF>=4 {x=$2/1000.0; y=1.0+2.0*x+0.5*x*x; print x,y}' "$REF_CAT" > "${P}_poly2.dat"

echo "[03] make weighted data"
awk 'NF>=4 {x=$2/1000.0; y=1.0+2.0*x+0.5*x*x; e=0.01+0.001*NR; print x,y,e}' "$REF_CAT" > "${P}_weighted.dat"

echo "[04] make outlier data"
awk 'NF>=4 {x=$2/1000.0; y=1.0+2.0*x+0.5*x*x; if (NR==5 || NR==10) y+=10.0; print x,y}' "$REF_CAT" > "${P}_outlier.dat"

echo "[05] make paired catalog data"
paste "$REF_CAT" "$OBS_CAT" | awk 'NF>=8 {print $2,$3,$6,$7,$4,$8}' > "${P}_pairs.dat"

for f in \
  "${P}_linear.dat" \
  "${P}_poly2.dat" \
  "${P}_weighted.dat" \
  "${P}_outlier.dat" \
  "${P}_pairs.dat"
do
    test -s "$f"
    echo "$f rows=$(awk 'NF>0{n++} END{print n+0}' "$f")"
    head -3 "$f"
done

echo "ALL LFIT TEST DATA CREATED"