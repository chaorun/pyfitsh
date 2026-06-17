#ifndef _TRANS_CORE_H
#define _TRANS_CORE_H

int trans_apply_cy(
    double *in_x, double *in_y, int npts,
    double *out_x, double *out_y,
    int order, int is_invert,
    double ox, double oy, double scale,
    double *dxfit, double *dyfit);

#endif
