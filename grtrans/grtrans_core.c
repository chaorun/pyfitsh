/*****************************************************************************/
/* _trans_core.c — thin wrapper: apply 2D polynomial transformation        */
/*                  to flat arrays of (x,y) points.                         */
/* Extracted from transform.c, depends only on math/poly.c (eval_2d_poly).  */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "math/poly.h"
#include "grtrans_core.h"

#define EPS 1e-10

/*****************************************************************************/
/* Jacobian of the polynomial transformation (from transform.c:395)          */
/*****************************************************************************/

void trans_get_jacobi(int order, double scale,
    double *dxfit, double *dyfit,
    double **rjxx, double **rjxy, double **rjyx, double **rjyy)
{
    int jnvar, i, j, k;
    double *jxx, *jxy, *jyx, *jyy;

    jnvar = (order + 0) * (order + 1) / 2;

    jxx = (double *)malloc(sizeof(double) * jnvar);
    jxy = (double *)malloc(sizeof(double) * jnvar);
    jyx = (double *)malloc(sizeof(double) * jnvar);
    jyy = (double *)malloc(sizeof(double) * jnvar);

    for (i = 0; i < jnvar; i++)
        jxx[i] = jxy[i] = jyx[i] = jyy[i] = 0.0;

    for (i = 0, k = 0; i <= order - 1; i++)
    {
        for (j = 0; j <= i; j++, k++)
        {
            jxx[k] += dxfit[k + i + 1] / scale;
            jxy[k] += dxfit[k + i + 2] / scale;
            jyx[k] += dyfit[k + i + 1] / scale;
            jyy[k] += dyfit[k + i + 2] / scale;
        }
    }

    *rjxx = jxx;
    *rjxy = jxy;
    *rjyx = jyx;
    *rjyy = jyy;
}

/*****************************************************************************/
/* Forward transformation (from transform.c:428)                             */
/*****************************************************************************/

static void trans_eval_normal(double x, double y,
    int order, double *dxfit, double *dyfit,
    double ox, double oy, double scale,
    double *rx, double *ry)
{
    *rx = eval_2d_poly(x, y, order, dxfit, ox, oy, scale);
    *ry = eval_2d_poly(x, y, order, dyfit, ox, oy, scale);
}

/*****************************************************************************/
/* Inverse transformation via Newton iteration (from transform.c:437)        */
/*****************************************************************************/

static void trans_eval_invert(double x, double y,
    int order, double *dxfit, double *dyfit,
    double ox, double oy, double scale,
    double *jxx, double *jxy, double *jyx, double *jyy,
    double *rx, double *ry)
{
    double wx, wy, mxx, mxy, myx, myy, det, imxx, imxy, imyx, imyy;
    double x0, y0, px0, py0, dx, dy;
    int i, n;

    if (order < 1) { *rx = x; *ry = y; return; }

    wx = x - dxfit[0];
    wy = y - dyfit[0];
    mxx = dxfit[1]; mxy = dxfit[2];
    myx = dyfit[1]; myy = dyfit[2];
    det = 1.0 / (mxx * myy - mxy * myx);
    imxx = +myy * det;
    imxy = -mxy * det;
    imyx = -myx * det;
    imyy = +mxx * det;
    x0 = imxx * wx + imxy * wy;
    y0 = imyx * wx + imyy * wy;

    n = 100;
    for (i = 0; i < n && order >= 2; i++)
    {
        mxx = eval_2d_poly(x0, y0, order - 1, jxx, ox, oy, scale);
        mxy = eval_2d_poly(x0, y0, order - 1, jxy, ox, oy, scale);
        myx = eval_2d_poly(x0, y0, order - 1, jyx, ox, oy, scale);
        myy = eval_2d_poly(x0, y0, order - 1, jyy, ox, oy, scale);

        det = 1.0 / (mxx * myy - mxy * myx);
        imxx = +myy * det; imxy = -mxy * det;
        imyx = -myx * det; imyy = +mxx * det;
        px0 = eval_2d_poly(x0, y0, order, dxfit, ox, oy, scale);
        py0 = eval_2d_poly(x0, y0, order, dyfit, ox, oy, scale);
        dx = x - px0;
        dy = y - py0;
        x0 += imxx * dx + imxy * dy;
        y0 += imyx * dx + imyy * dy;
        if ( fabs(x0-px0) < (fabs(x0)+fabs(px0))*EPS &&
             fabs(y0-py0) < (fabs(y0)+fabs(py0))*EPS ) break;
    }

    *rx = x0;
    *ry = y0;
}

/*****************************************************************************/
/* Main entry point: apply transformation to flat arrays                     */
/*****************************************************************************/

int grtrans_apply_cy(
    double *in_x, double *in_y, int npts,
    double *out_x, double *out_y,
    int order, int is_invert,
    double ox, double oy, double scale,
    double *dxfit, double *dyfit)
{
    int i;
    double *jxx = NULL, *jxy = NULL, *jyx = NULL, *jyy = NULL;

    if (is_invert && order >= 1) {
        int all_zero = 1;
        int nvar = (order + 1) * (order + 2) / 2;
        for (i = 0; i < nvar; i++) {
            if (dxfit[i] != 0.0 || dyfit[i] != 0.0) { all_zero = 0; break; }
        }
        if (all_zero) return 1;
        trans_get_jacobi(order, scale, dxfit, dyfit, &jxx, &jxy, &jyx, &jyy);
    }

    for (i = 0; i < npts; i++)
    {
        if (is_invert && order >= 1)
            trans_eval_invert(in_x[i], in_y[i], order, dxfit, dyfit,
                ox, oy, scale, jxx, jxy, jyx, jyy, &out_x[i], &out_y[i]);
        else
            trans_eval_normal(in_x[i], in_y[i], order, dxfit, dyfit,
                ox, oy, scale, &out_x[i], &out_y[i]);
    }

    free(jxx); free(jxy); free(jyx); free(jyy);
    return 0;
}
