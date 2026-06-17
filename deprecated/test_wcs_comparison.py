#!/usr/bin/env python3
"""test_wcs_comparison.py — Compare fitsh_cython WCS fitting vs FITS header WCS

Pipeline:
  1. Read ref_L1_...cat → pixel (x_ref, y_ref) of reference stars
  2. Use FITS header WCS (astropy) → (ra_ref, dec_ref) for each star
  3. Split into train/test sets
  4. Fit WCS via fitsh_cython Grtrans.WCSFitter
  5. Compare:
     a) Fitted WCS residual (xy vs radec)
     b) Fitted WCS params vs original FITS header WCS params
     c) Compare with astropy's own WCS fitting (if available)
"""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS as AstroWCS
from fitsh_cy import Grmatch, Grtrans

DATA = '../fitsh_test_data'
FITS_FILE = f'{DATA}/L1_GRB260604C-Ic-180s-20260606_210318.fits'
REF_CAT = f'{DATA}/ref_L1_GRB260604C-Ic-180s-20260606_210318.cat'

# ---- 1. Load reference stars (id, x, y, mag) ----
print('=== 1. Load reference stars ===')
ref_data = np.loadtxt(REF_CAT)
ref_id  = ref_data[:, 0].astype(int)
ref_x   = ref_data[:, 1]
ref_y   = ref_data[:, 2]
ref_mag = ref_data[:, 3]
print(f'  Stars: {len(ref_x)}')

# ---- 2. Use FITS header WCS to compute (ra, dec) ----
print('=== 2. Compute RA/Dec from FITS header WCS ===')
hdul = fits.open(FITS_FILE)
orig_hdr = hdul[0].header
orig_wcs = AstroWCS(orig_hdr)
ra_orig, dec_orig = orig_wcs.all_pix2world(ref_x, ref_y, 1)
print(f'  Original WCS: CTYPE={orig_hdr["CTYPE1"]}/{orig_hdr["CTYPE2"]}')
print(f'  CRPIX=({orig_hdr["CRPIX1"]:.2f},{orig_hdr["CRPIX2"]:.2f})')
print(f'  CRVAL=({orig_hdr["CRVAL1"]:.4f},{orig_hdr["CRVAL2"]:.4f})')
print(f'  SIP order: A={orig_hdr.get("A_ORDER","N/A")}')
print(f'  RA range: [{ra_orig.min():.4f}, {ra_orig.max():.4f}]')
print(f'  DEC range: [{dec_orig.min():.4f}, {dec_orig.max():.4f}]')

# ---- 3. Fit WCS via fitsh_cython Grtrans.WCSFitter ----
print('=== 3. Fit WCS via fitsh_cython ===')
ra0 = orig_hdr['CRVAL1']
dec0 = orig_hdr['CRVAL2']
proj_type = orig_hdr['CTYPE1'].split('---')[1].replace('-SIP', '').lower()

for order in [2, 3, 4]:
    wcs_fitter = Grtrans.WCSFitter(ra0=ra0, dec0=dec0, projection=proj_type, order=order)
    wcs_fitter.fit(ra_orig, dec_orig, ref_x, ref_y)

    # Residual: how well does the fitted WCS reproduce the input?
    x_fit, y_fit = wcs_fitter.apply_to_radec(ra_orig, dec_orig)
    x_resid = ref_x - x_fit
    y_resid = ref_y - y_fit
    rms = np.sqrt(np.mean(x_resid**2 + y_resid**2))

    # Also test on ra/dec → pixel → ra/dec roundtrip
    ra_back, dec_back = wcs_fitter.pixel_to_radec(ref_x, ref_y)
    ra_diff_arcsec = np.abs(ra_back - ra_orig) * 3600
    dec_diff_arcsec = np.abs(dec_back - dec_orig) * 3600

    print(f'\n  Order {order}:')
    print(f'    XY residual  RMS: {rms:.4f} px')
    print(f'    XY residual  max: {np.sqrt(x_resid**2+y_resid**2).max():.4f} px')
    print(f'    RA roundtrip RMS: {np.mean(ra_diff_arcsec):.4f} arcsec')
    print(f'    DEC roundtrip RMS: {np.mean(dec_diff_arcsec):.4f} arcsec')
    print(f'    CRPIX: ({wcs_fitter.crpix1:.2f}, {wcs_fitter.crpix2:.2f})')
    print(f'    CD: [[{wcs_fitter.cd11:.6e}, {wcs_fitter.cd12:.6e}],')
    print(f'         [{wcs_fitter.cd21:.6e}, {wcs_fitter.cd22:.6e}]]')

    # Generate FITS header and compare with astropy
    if order >= 2:
        fitted_hdr = wcs_fitter.to_fits_header()
        fitted_wcs = AstroWCS(fitted_hdr)
        ra_sip, dec_sip = fitted_wcs.all_pix2world(ref_x, ref_y, 1)
        ra_sip_diff = np.abs(ra_sip - ra_orig) * 3600
        dec_sip_diff = np.abs(dec_sip - dec_orig) * 3600
        print(f'    SIP RA roundtrip RMS: {np.mean(ra_sip_diff):.4f} arcsec')
        print(f'    SIP DEC roundtrip RMS: {np.mean(dec_sip_diff):.4f} arcsec')

# ---- 4. Compare fitted WCS params vs original header ----
print('\n=== 4. WCS parameter comparison (order=3) ===')
wcs_fitter = Grtrans.WCSFitter(ra0=ra0, dec0=dec0, projection=proj_type, order=3)
wcs_fitter.fit(ra_orig, dec_orig, ref_x, ref_y)

print(f'  {"Param":10s} {"Fitted":>16s} {"Original":>16s} {"Diff":>16s}')
for name, fit_val, orig_val in [
    ('CRPIX1', wcs_fitter.crpix1, orig_hdr['CRPIX1']),
    ('CRPIX2', wcs_fitter.crpix2, orig_hdr['CRPIX2']),
    ('CD1_1', wcs_fitter.cd11, orig_hdr.get('CD1_1', 0)),
    ('CD1_2', wcs_fitter.cd12, orig_hdr.get('CD1_2', 0)),
    ('CD2_1', wcs_fitter.cd21, orig_hdr.get('CD2_1', 0)),
    ('CD2_2', wcs_fitter.cd22, orig_hdr.get('CD2_2', 0)),
]:
    diff = fit_val - orig_val
    print(f'  {name:10s} {fit_val:16.6e} {orig_val:16.6e} {diff:16.6e}')

# ---- 5. Compare with astropy.utils.fit_wcs (from astropy v5.3+) ----
print('\n=== 5. Compare with astropy fit_wcs ===')
try:
    from astropy.wcs.utils import fit_wcs_from_points
    from astropy.coordinates import SkyCoord
    astro_wcs = fit_wcs_from_points(
        (ref_x, ref_y),
        SkyCoord(ra=ra_orig, dec=dec_orig, unit='deg'),
        proj_point=SkyCoord(ra=ra0, dec=dec0, unit='deg'),
        projection='TAN',
        sip_degree=3,
    )
    astro_ra, astro_dec = astro_wcs.all_pix2world(ref_x, ref_y, 1)
    astro_ra_diff = np.abs(astro_ra - ra_orig) * 3600
    astro_dec_diff = np.abs(astro_dec - dec_orig) * 3600
    print(f'  astropy.fit_wcs RA RMS: {np.mean(astro_ra_diff):.4f} arcsec')
    print(f'  astropy.fit_wcs DEC RMS: {np.mean(astro_dec_diff):.4f} arcsec')

    # Compare fitsh_cython vs astropy
    cy_ra, cy_dec = wcs_fitter.pixel_to_radec(ref_x, ref_y)
    cy_ra_diff = np.abs(cy_ra - ra_orig) * 3600
    cy_dec_diff = np.abs(cy_dec - dec_orig) * 3600
    print(f'  fitsh_cython RA RMS: {np.mean(cy_ra_diff):.4f} arcsec')
    print(f'  fitsh_cython DEC RMS: {np.mean(cy_dec_diff):.4f} arcsec')

    # Compare WCS parameters
    print(f'\n  {"Param":10s} {"fitsh_cython":>16s} {"astropy":>16s}')
    for name, cy_val, astro_val in [
        ('CRPIX1', wcs_fitter.crpix1, astro_wcs.wcs.crpix[0]),
        ('CRPIX2', wcs_fitter.crpix2, astro_wcs.wcs.crpix[1]),
        ('CD1_1', wcs_fitter.cd11, astro_wcs.wcs.cd[0,0]),
        ('CD1_2', wcs_fitter.cd12, astro_wcs.wcs.cd[0,1]),
        ('CD2_1', wcs_fitter.cd21, astro_wcs.wcs.cd[1,0]),
        ('CD2_2', wcs_fitter.cd22, astro_wcs.wcs.cd[1,1]),
    ]:
        print(f'  {name:10s} {cy_val:16.6e} {astro_val:16.6e}')

except ImportError:
    print('  astropy.wcs.utils.fit_wcs not available (requires astropy >= 5.3)')

hdul.close()
print('\n=== DONE ===')
