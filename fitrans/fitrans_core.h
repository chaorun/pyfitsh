#ifndef _FITRANS_CORE_H
#define _FITRANS_CORE_H

int fitrans_apply_cy(
    double *in_data, char *in_mask, int sx, int sy,
    double *out_data, char *out_mask, int nsx, int nsy,
    int ofx, int ofy,
    int order, int int_method, int is_invert,
    double ox, double oy, double scale,
    double *dxfit, double *dyfit);

#endif
