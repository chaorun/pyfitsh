#ifndef _POLYFIT_CORE_H
#define _POLYFIT_CORE_H

int poly_fit_cy(
    double *x, double *y, double *vals, double *weights, int npts,
    int order, double ox, double oy, double scale,
    double *coeff_out);

int poly_compose_affine_cy(
    double *coeff, int order,
    double xlin0, double xlin1, double xlin2,
    double ylin0, double ylin1, double ylin2,
    double *rc);

#endif
