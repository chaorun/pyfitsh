#ifndef _WCS_CORE_H
#define _WCS_CORE_H

typedef struct {
    double ra0, de0, roll0;
    int type, order;
    int n_k, n_p;
    double k[4], p[4];
} wcs_fit_input;

typedef struct {
    double crpix1, crpix2;
    double cd11, cd12, cd21, cd22;
    double res_prj, res_pix;
    int order, nvar;
    double *pixpoly_dx;   /* forward: pix -> proj, nvar coeffs */
    double *pixpoly_dy;
    double *prjpoly_dx;   /* inverse: proj -> pix, nvar coeffs */
    double *prjpoly_dy;
} wcs_fit_output;

int wcs_fit_cy(
    double *ra, double *dec, double *px, double *py, int npts,
    const wcs_fit_input *cfg,
    wcs_fit_output *out);

void wcs_fit_output_free(wcs_fit_output *out);

#endif
