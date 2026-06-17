/*****************************************************************************/
/* _wcs_core.c — WCS fitting: extract grtrans_wcs_perform_fit               */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "math/poly.h"
#include "math/polyfit.h"
#include "projection.h"
#include "wcs_core.h"

/*****************************************************************************/
/* Structs from grtrans.c                                                     */
/*****************************************************************************/

typedef struct {
    double x, y;
    double weight, delta;
    double *vals;
} fitpoint;

typedef struct {
    double  ra0, de0, roll0;
    double  qr0, qi0, qj0, qk0;
    double  max_distance;
    int     type, order;
    matrix  mproj;
    double  zfactor;
    projection_distort distort;
} projection_data;

typedef struct {
    projection_data init;
    double *pixpoly[2];
    double *prjpoly[2];
    double crpix1, crpix2;
    double cd11, cd12, cd21, cd22;
    double res_prj, res_pix;
} wcs_data;

/*****************************************************************************/
/* grtrans_wcs_perform_fit (from grtrans.c:106-306, unchanged)               */
/*****************************************************************************/

static int wcs_perform_fit(fitpoint *fps, int nfp, projection_data *wi, wcs_data *wcs)
{
    int i, k, m, nvar, order;
    point *p_pix, *p_prj, *ps;
    double *vfits[4], cd11, cd12, cd21, cd22, crpix1, crpix2, det, sig2;
    double ra, de, x, y, z;

    if (wi == NULL || wcs == NULL) return (-1);
    memmove(&wcs->init, wi, sizeof(projection_data));

    order = wcs->init.order;
    if (order <= 0) return (-1);

    nvar = (order + 1) * (order + 2) / 2;

    for (k = 0; k < 4; k++)
        vfits[k] = (double *)malloc(sizeof(double) * nvar);

    p_pix = (point *)malloc(sizeof(point) * nfp);
    p_prj = (point *)malloc(sizeof(point) * nfp);
    ps = (point *)malloc(sizeof(point) * nfp);

    /* Stage 1: prepare projection */
    projection_get_matrix_rdr(wcs->init.ra0, wcs->init.de0, 0.0, wcs->init.mproj);
    for (i = 0; i < nfp; i++) {
        ra = fps[i].x; de = fps[i].y;
        projection_do_coord(ra, de, wcs->init.ra0, wcs->init.de0, 0.0, &x, &y, &z);
        projection_do_distortion(wcs->init.type, NULL, &x, &y, &z);
        p_prj[i].x = x * M_R2D;
        p_prj[i].y = y * M_R2D;
        p_pix[i].x = fps[i].vals[0];
        p_pix[i].y = fps[i].vals[1];
    }

    /* Stage 2: initial linear fit for CD and CRPIX */
    for (i = 0; i < nfp; i++) { ps[i].x = p_pix[i].x; ps[i].y = p_pix[i].y; ps[i].weight = 1.0; }
    for (i = 0; i < nfp; i++) ps[i].value = p_prj[i].x;
    fit_2d_poly(ps, nfp, 1, vfits[0], 0, 0, 1);
    for (i = 0; i < nfp; i++) ps[i].value = p_prj[i].y;
    fit_2d_poly(ps, nfp, 1, vfits[1], 0, 0, 1);
    cd11 = vfits[0][1]; cd12 = vfits[0][2];
    cd21 = vfits[1][1]; cd22 = vfits[1][2];
    det = (cd11 * cd22 - cd12 * cd21);
    crpix1 = -(+cd22 * vfits[0][0] - cd12 * vfits[1][0]) / det;
    crpix2 = -(-cd21 * vfits[0][0] + cd11 * vfits[1][0]) / det;

    /* Stage 3: higher-order distortion refinement */
    if (2 <= order) {
        for (k = 0; k < 3; k++) {
            for (i = 0; i < nfp; i++) {
                ps[i].x = p_pix[i].x - crpix1;
                ps[i].y = p_pix[i].y - crpix2;
                ps[i].weight = 1.0;
            }
            for (i = 0; i < nfp; i++) ps[i].value = p_prj[i].x;
            fit_2d_poly(ps, nfp, order, vfits[0], 0, 0, 1);
            for (i = 0; i < nfp; i++) ps[i].value = p_prj[i].y;
            fit_2d_poly(ps, nfp, order, vfits[1], 0, 0, 1);
            cd11 = vfits[0][1]; cd12 = vfits[0][2];
            cd21 = vfits[1][1]; cd22 = vfits[1][2];
            det = (cd11 * cd22 - cd12 * cd21);
            crpix1 += -(+cd22 * vfits[0][0] - cd12 * vfits[1][0]) / det;
            crpix2 += -(-cd21 * vfits[0][0] + cd11 * vfits[1][0]) / det;
        }

        for (k = 0; k < 1; k++) {
            double cx11, cx12, cx21, cx22;
            det = (cd11 * cd22 - cd12 * cd21);
            for (i = 0; i < nfp; i++) {
                ps[i].x = p_pix[i].x - crpix1;
                ps[i].y = p_pix[i].y - crpix2; ps[i].weight = 1.0;
            }
            for (i = 0; i < nfp; i++) ps[i].value = (+cd22 * p_prj[i].x - cd12 * p_prj[i].y) / det;
            fit_2d_poly(ps, nfp, order, vfits[0], 0, 0, 1);
            for (i = 0; i < nfp; i++) ps[i].value = (-cd21 * p_prj[i].x + cd11 * p_prj[i].y) / det;
            fit_2d_poly(ps, nfp, order, vfits[1], 0, 0, 1);

            cx11 = cd11 * vfits[0][1] + cd12 * vfits[1][1];
            cx12 = cd11 * vfits[0][2] + cd12 * vfits[1][2];
            cx21 = cd21 * vfits[0][1] + cd22 * vfits[1][1];
            cx22 = cd21 * vfits[0][2] + cd22 * vfits[1][2];
            cd11 = cx11; cd12 = cx12; cd21 = cx21; cd22 = cx22;
        }
    }

    /* Force standard convention: linear terms in forward polynomial */
    vfits[0][0] = 0.0; vfits[0][1] = 1.0; vfits[0][2] = 0.0;
    vfits[1][0] = 0.0; vfits[1][1] = 0.0; vfits[1][2] = 1.0;

    /* Residual of forward fit */
    sig2 = 0.0;
    for (i = 0; i < nfp; i++) {
        double px1, py1, fpx, fpy;
        px1 = p_pix[i].x - crpix1; py1 = p_pix[i].y - crpix2;
        x = eval_2d_poly(px1, py1, order, vfits[0], 0, 0, 1);
        y = eval_2d_poly(px1, py1, order, vfits[1], 0, 0, 1);
        fpx = cd11 * x + cd12 * y;
        fpy = cd21 * x + cd22 * y;
        fpx -= p_prj[i].x; fpy -= p_prj[i].y;
        sig2 += fpx * fpx + fpy * fpy;
    }
    if (0 < nfp) sig2 = sqrt(sig2 / (double)nfp);
    wcs->res_prj = sig2;

    /* Inverse fit: projected -> pixel */
    det = (cd11 * cd22 - cd12 * cd21);
    for (i = 0; i < nfp; i++) {
        ps[i].x = (+cd22 * p_prj[i].x - cd12 * p_prj[i].y) / det;
        ps[i].y = (-cd21 * p_prj[i].x + cd11 * p_prj[i].y) / det;
        ps[i].weight = 1.0;
    }
    for (i = 0; i < nfp; i++) ps[i].value = p_pix[i].x - crpix1;
    fit_2d_poly(ps, nfp, order, vfits[2], 0, 0, 1);
    for (i = 0; i < nfp; i++) ps[i].value = p_pix[i].y - crpix2;
    fit_2d_poly(ps, nfp, order, vfits[3], 0, 0, 1);
    vfits[2][0] = 0.0; vfits[3][0] = 0.0;

    /* Residual of inverse fit */
    det = (cd11 * cd22 - cd12 * cd21);
    sig2 = 0.0;
    for (i = 0; i < nfp; i++) {
        double px1, py1, fpx, fpy;
        px1 = (+cd22 * p_prj[i].x - cd12 * p_prj[i].y) / det;
        py1 = (-cd21 * p_prj[i].x + cd11 * p_prj[i].y) / det;
        x = eval_2d_poly(px1, py1, order, vfits[2], 0, 0, 1);
        y = eval_2d_poly(px1, py1, order, vfits[3], 0, 0, 1);
        fpx = x - (p_pix[i].x - crpix1);
        fpy = y - (p_pix[i].y - crpix2);
        sig2 += fpx * fpx + fpy * fpy;
    }
    if (0 < nfp) sig2 = sqrt(sig2 / (double)nfp);
    wcs->res_pix = sig2;

    /* Save */
    wcs->crpix1 = crpix1; wcs->crpix2 = crpix2;
    wcs->cd11 = cd11; wcs->cd12 = cd12;
    wcs->cd21 = cd21; wcs->cd22 = cd22;

    /* Rescale: divide A_ij by (i!j!) */
    k = 0; m = 1;
    for (i = 0; i <= order; i++) {
        int j, w;
        w = m;
        for (j = 0; j <= i; j++) {
            vfits[0][k] /= (double)w;
            vfits[1][k] /= (double)w;
            vfits[2][k] /= (double)w;
            vfits[3][k] /= (double)w;
            if (j < i) w = (w / (i - j)) * (j + 1);
            k++;
        }
        m = m * (i + 1);
    }
    /* Subtract unity from linear part of inverse (standard SIP convention) */
    vfits[2][1] -= 1.0;
    vfits[3][2] -= 1.0;

    wcs->pixpoly[0] = vfits[0];
    wcs->pixpoly[1] = vfits[1];
    wcs->prjpoly[0] = vfits[2];
    wcs->prjpoly[1] = vfits[3];

    free(ps);
    free(p_prj);
    free(p_pix);

    return 0;
}

/*****************************************************************************/
/* Flat wrapper                                                              */
/*****************************************************************************/

int wcs_fit_cy(
    double *ra, double *dec, double *px, double *py, int npts,
    const wcs_fit_input *cfg,
    wcs_fit_output *out)
{
    int i, nvar;
    fitpoint *fps;
    projection_data wi;
    wcs_data wcs_data_static, *wcs = &wcs_data_static;
    memset(wcs, 0, sizeof(wcs_data));

    memset(&wi, 0, sizeof(wi));
    wi.ra0 = cfg->ra0;
    wi.de0 = cfg->de0;
    wi.roll0 = cfg->roll0;
    wi.type = cfg->type;
    wi.order = cfg->order;
    wi.max_distance = 0.0;
    wi.zfactor = M_R2D;
    wi.distort.brown_conrady_radial_ncoeff = cfg->n_k;
    for (i = 0; i < cfg->n_k && i < 4; i++)
        wi.distort.brown_conrady_radial_coeffs[i] = cfg->k[i];
    wi.distort.brown_conrady_tangential_ncoeff = cfg->n_p;
    for (i = 0; i < cfg->n_p && i < 4; i++)
        wi.distort.brown_conrady_tangential_coeffs[i] = cfg->p[i];

    fps = (fitpoint *)malloc(sizeof(fitpoint) * npts);
    for (i = 0; i < npts; i++) {
        fps[i].x = ra[i];
        fps[i].y = dec[i];
        fps[i].weight = 1.0;
        fps[i].delta = 0.0;
        fps[i].vals = (double *)malloc(sizeof(double) * 2);
        fps[i].vals[0] = px[i];
        fps[i].vals[1] = py[i];
    }

    int r = wcs_perform_fit(fps, npts, &wi, wcs);
    if (r < 0) {
        for (i = 0; i < npts; i++) free(fps[i].vals);
        free(fps);
        return r;
    }

    nvar = (cfg->order + 1) * (cfg->order + 2) / 2;
    out->crpix1 = wcs->crpix1 + 0.5;
    out->crpix2 = wcs->crpix2 + 0.5;
    out->cd11 = wcs->cd11; out->cd12 = wcs->cd12;
    out->cd21 = wcs->cd21; out->cd22 = wcs->cd22;
    out->res_prj = wcs->res_prj;
    out->res_pix = wcs->res_pix;
    out->order = cfg->order;
    out->nvar = nvar;

    out->pixpoly_dx = (double *)malloc(sizeof(double) * nvar);
    out->pixpoly_dy = (double *)malloc(sizeof(double) * nvar);
    out->prjpoly_dx = (double *)malloc(sizeof(double) * nvar);
    out->prjpoly_dy = (double *)malloc(sizeof(double) * nvar);
    memcpy(out->pixpoly_dx, wcs->pixpoly[0], sizeof(double) * nvar);
    memcpy(out->pixpoly_dy, wcs->pixpoly[1], sizeof(double) * nvar);
    memcpy(out->prjpoly_dx, wcs->prjpoly[0], sizeof(double) * nvar);
    memcpy(out->prjpoly_dy, wcs->prjpoly[1], sizeof(double) * nvar);

    for (i = 0; i < 4; i++) {
        if (wcs->pixpoly[i / 2] != NULL) { free(wcs->pixpoly[i / 2]); wcs->pixpoly[i / 2] = NULL; }
    }
    for (i = 0; i < npts; i++)
        free(fps[i].vals);
    free(fps);

    return 0;
}

void wcs_fit_output_free(wcs_fit_output *out)
{
    free(out->pixpoly_dx); free(out->pixpoly_dy);
    free(out->prjpoly_dx); free(out->prjpoly_dy);
}
