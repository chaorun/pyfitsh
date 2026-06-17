/*****************************************************************************/
/* _polyfit_core.c — flat wrappers around fit_2d_poly and compose_affine    */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "math/poly.h"
#include "math/polyfit.h"
#include "polyfit_core.h"

int poly_fit_cy(
    double *x, double *y, double *vals, double *weights, int npts,
    int order, double ox, double oy, double scale,
    double *coeff_out)
{
    int i, r, nvar;
    point *data;

    nvar = (order + 1) * (order + 2) / 2;
    for (i = 0; i < nvar; i++) coeff_out[i] = 0.0;

    data = (point *)malloc((size_t)npts * sizeof(point));
    if (data == NULL) return -1;

    for (i = 0; i < npts; i++) {
        data[i].x = x[i];
        data[i].y = y[i];
        data[i].value = vals[i];
        data[i].weight = (weights != NULL) ? weights[i] : 1.0;
    }

    r = fit_2d_poly(data, npts, order, coeff_out, ox, oy, scale);
    free(data);
    return r;
}

int poly_compose_affine_cy(
    double *coeff, int order,
    double xlin0, double xlin1, double xlin2,
    double ylin0, double ylin1, double ylin2,
    double *rc)
{
    double xlin[3] = {xlin0, xlin1, xlin2};
    double ylin[3] = {ylin0, ylin1, ylin2};
    return compose_2d_poly_with_affine(coeff, order, xlin, ylin, rc);
}
