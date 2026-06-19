#!/usr/bin/env bash
set -euo pipefail

REF_IMG="L1_GRB260604C-Ic-180s-20260606_210318.fits"
OBS_IMG="L1_GRB260604C-Ic-180s-20260607_023443.fits"
REF_CAT="ref_L1_GRB260604C-Ic-180s-20260606_210318.cat"
OBS_CAT="obs_L1_GRB260604C-Ic-180s-20260607_023443.cat"

PREFIX="cli_test"

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

echo "[01] grmatch"
grmatch \
  --input "$OBS_CAT" \
  --input-reference "$REF_CAT" \
  --col-ref 2,3 \
  --col-inp 2,3 \
  --match-points \
  --triangulation auto,unitarity=0.2 \
  --col-ref-ordering -4 \
  --col-inp-ordering -4 \
  --max-distance 2 \
  --output "${PREFIX}_match.cat" \
  --output-transformation "${PREFIX}_grb.trans"
check_table "${PREFIX}_match.cat"
check_file "${PREFIX}_grb.trans"

echo "[02] fitrans inverse"
fitrans \
  --input "$OBS_IMG" \
  -T "${PREFIX}_grb.trans" \
  --output "${PREFIX}_v1.fits" \
  --inverse
check_fits "${PREFIX}_v1.fits"

echo "[03] ficonv spatial kernel"
ficonv \
  -r "$REF_IMG" \
  -i "${PREFIX}_v1.fits" \
  -k "i/1;b/1;g=5,2.0,2/1" \
  --output-kernel-list "${PREFIX}_kernel_spatial1.txt" \
  -o "${PREFIX}_ref_conv_spatial1.fits" \
  --output-subtracted "${PREFIX}_sub_spatial1.fits" \
  -n 2 \
  -s 3
check_file "${PREFIX}_kernel_spatial1.txt"
check_fits "${PREFIX}_ref_conv_spatial1.fits"
check_fits "${PREFIX}_sub_spatial1.fits"

echo "[04] ficonv read fitted kernel"
ficonv \
  -r "$REF_IMG" \
  --input-kernel-list "${PREFIX}_kernel_spatial1.txt" \
  -o "${PREFIX}_ref_conv_from_kernel.fits"
check_fits "${PREFIX}_ref_conv_from_kernel.fits"

echo "[05] ficonv gaussian multi-kernel"
ficonv \
  -r "$REF_IMG" \
  -i "${PREFIX}_v1.fits" \
  -k "i/0;b/1;g=3,1.0,2/0;g=5,2.0,2/0;g=7,4.0,2/0" \
  --output-kernel-list "${PREFIX}_kernel_gauss.txt" \
  -o "${PREFIX}_ref_conv_gauss.fits" \
  --output-subtracted "${PREFIX}_sub_gauss.fits" \
  -n 3 \
  -s 3
check_file "${PREFIX}_kernel_gauss.txt"
check_fits "${PREFIX}_ref_conv_gauss.fits"
check_fits "${PREFIX}_sub_gauss.fits"

echo "[06] ficonv discrete kernel"
ficonv \
  -r "$REF_IMG" \
  -i "${PREFIX}_v1.fits" \
  -k "i/0;b/0;d=2/0" \
  --output-kernel-list "${PREFIX}_kernel_discrete.txt" \
  -o "${PREFIX}_ref_conv_discrete.fits" \
  --output-subtracted "${PREFIX}_sub_discrete.fits" \
  -n 2 \
  -s 3
check_file "${PREFIX}_kernel_discrete.txt"
check_fits "${PREFIX}_ref_conv_discrete.fits"
check_fits "${PREFIX}_sub_discrete.fits"

echo "[07] fistar gauss full + PSF"
fistar \
  --input "$REF_IMG" \
  -d 5 \
  -f 10000 \
  -s x \
  --format id,ix,iy,cx,cy,cbg,camp,cmax,npix,cs,cd,ck,x,y,bg,amp,s,d,k,l,sigma,delta,kappa,fwhm,ellip,pa,flux,noise,s/n,magnitude,px,py,pbg,pamp,ps,pd,pk,pl \
  -g 2 \
  --model gauss \
  --output "${PREFIX}_fistar_gauss_full.cat" \
  --psf native,order=2 \
  --output-psf "${PREFIX}_native_o2.psf.fits" \
  -V
check_table "${PREFIX}_fistar_gauss_full.cat"
check_fits "${PREFIX}_native_o2.psf.fits"

echo "[07b] fistar deviated moment column"
fistar \
  --input "$REF_IMG" \
  -d 5 \
  -f 10000 \
  -g 2 \
  --model deviated,order=2 \
  --format id,x,y,bg,amp,s,d,k,mom,l,sigma,delta,kappa,fwhm,ellip,pa,flux,noise,s/n,magnitude \
  --output "${PREFIX}_fistar_deviated_mom.cat"
check_table "${PREFIX}_fistar_deviated_mom.cat"

echo "[08] fistar elliptic"
fistar \
  --input "$REF_IMG" \
  -d 5 \
  -f 10000 \
  -g 2 \
  --model elliptic \
  --format id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,s/n,magnitude \
  --output "${PREFIX}_fistar_elliptic.cat"
check_table "${PREFIX}_fistar_elliptic.cat"

echo "[09] fistar parabolapeak only-candidates"
fistar \
  --input "$REF_IMG" \
  -d 5 \
  -f 10000 \
  --algorithm parabolapeak \
  --only-candidates \
  --format id,ix,iy,cx,cy,cbg,camp,cmax,npix,cs,cd,ck \
  --output "${PREFIX}_parabolapeak_candidates.cat"
check_table "${PREFIX}_parabolapeak_candidates.cat"

echo "[10] fistar input-candidates path"
fistar \
  --input "$REF_IMG" \
  -d 5 \
  -f 10000 \
  --only-candidates \
  --format cx,cy,cs,cd,ck \
  --output "${PREFIX}_input_candidates.dat"
check_table "${PREFIX}_input_candidates.dat"

fistar \
  --input "$REF_IMG" \
  --input-candidates "${PREFIX}_input_candidates.dat" \
  --col-xy 1,2 \
  --col-shape 3,4,5 \
  -g 2 \
  --model elliptic \
  --format id,x,y,bg,amp,s,d,k,fwhm,ellip,pa,flux,noise,s/n \
  --output "${PREFIX}_from_candidates.cat"
check_table "${PREFIX}_from_candidates.cat"

echo "[11] fiphot single aperture + raw"
awk '{print $1, $13, $14}' "${PREFIX}_fistar_gauss_full.cat" > "${PREFIX}_fiphot_positions.dat"
check_table "${PREFIX}_fiphot_positions.dat"

fiphot \
  --input "$REF_IMG" \
  --input-list "${PREFIX}_fiphot_positions.dat" \
  --col-id 1 \
  --col-xy 2,3 \
  --apertures 10:18:12 \
  --gain 2 \
  --sky-fit median,iterations=2,sigma=3 \
  --format IXY,BbFfMm \
  --output "${PREFIX}_fiphot_single.cat" \
  --output-raw-photometry "${PREFIX}_fiphot_single.raw"
check_table "${PREFIX}_fiphot_single.cat"
check_file "${PREFIX}_fiphot_single.raw"

echo "[12] fiphot multi aperture"
fiphot \
  --input "$REF_IMG" \
  --input-list "${PREFIX}_fiphot_positions.dat" \
  --col-id 1 \
  --col-xy 2,3 \
  --apertures 8:18:12,10:18:12,12:20:15,15:25:15 \
  --gain 2 \
  --format IXY,BbFfMm \
  --output "${PREFIX}_fiphot_multiap.cat"
check_table "${PREFIX}_fiphot_multiap.cat"

echo "[13] fiphot subtraction photometry"
fiphot \
  -s "${PREFIX}_sub_spatial1.fits" \
  --input-raw-photometry "${PREFIX}_fiphot_single.raw" \
  --input-kernel "${PREFIX}_kernel_spatial1.txt" \
  --format IXY,BbFfMm \
  --output "${PREFIX}_fiphot_sub.cat"
check_file "${PREFIX}_fiphot_sub.cat"

echo "[14] firandom 2049 A"
firandom \
  --seed 11111 \
  --size 2049,2049 \
  --sky 1000 \
  --sky-noise 3 \
  --list "10[x=r(0.2,0.8),y=r(0.2,0.8),f=5,e=0.0,p=0,i=r(1000,3000)]" \
  --fep \
  --output-list "${PREFIX}_fr_a.lst" \
  -o "${PREFIX}_fr_a.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_table "${PREFIX}_fr_a.lst"
check_fits "${PREFIX}_fr_a.fits"

echo "[15] firandom 2049 B"
firandom \
  --seed 22222 \
  --size 2049,2049 \
  --sky 900 \
  --sky-noise 2 \
  --list "8[x=r(0.2,0.8),y=r(0.2,0.8),f=4,e=0.0,p=0,i=r(500,2000)]" \
  --fep \
  --output-list "${PREFIX}_fr_b.lst" \
  -o "${PREFIX}_fr_b.fits" \
  --bitpix -32 \
  --gain 2 \
  --photon-noise \
  --integral
check_table "${PREFIX}_fr_b.lst"
check_fits "${PREFIX}_fr_b.fits"

echo "[16] fiarith image subtraction"
fiarith \
  "'${PREFIX}_fr_a.fits'-'${PREFIX}_fr_b.fits'" \
  -o "${PREFIX}_fiarith_a_minus_b.fits"
check_fits "${PREFIX}_fiarith_a_minus_b.fits"

echo "[17] fiarith laplace"
fiarith \
  "laplace('${PREFIX}_fr_a.fits')" \
  -o "${PREFIX}_fiarith_laplace.fits"
check_fits "${PREFIX}_fiarith_laplace.fits"

echo "[18] fiarith sign"
fiarith \
  "sign('${PREFIX}_fr_a.fits')" \
  -o "${PREFIX}_fiarith_sign.fits"
check_fits "${PREFIX}_fiarith_sign.fits"

echo "[19] fiarith theta"
fiarith \
  "theta('${PREFIX}_fr_a.fits')" \
  -o "${PREFIX}_fiarith_theta.fits"
check_fits "${PREFIX}_fiarith_theta.fits"

echo "[20] summary"
python - <<'PY'
from pathlib import Path
files = sorted(Path(".").glob("cli_test_*"))
print(f"Generated files: {len(files)}")
for f in files:
    print(f"{f.name}\t{f.stat().st_size}")
PY

echo "ALL CLI TESTS PASSED"