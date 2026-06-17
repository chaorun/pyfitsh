# grtrans C API declarations
cdef extern from "grtrans_core.h":
    int grtrans_apply_cy(
        double *in_x, double *in_y, int npts,
        double *out_x, double *out_y,
        int order, int is_invert,
        double ox, double oy, double scale,
        double *dxfit, double *dyfit) nogil

    void trans_get_jacobi(int order, double scale,
        double *dxfit, double *dyfit,
        double **rjxx, double **rjxy, double **rjyx, double **rjyy) nogil

    int compose_2d_poly_with_affine(double *coeff, int order,
        double *xlin, double *ylin, double *rc) nogil

