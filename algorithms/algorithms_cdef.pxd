# algorithms C API declarations
cdef extern from "trans_core.h":
    int trans_apply_cy(
        double *in_x, double *in_y, int npts,
        double *out_x, double *out_y,
        int order, int is_invert,
        double ox, double oy, double scale,
        double *dxfit, double *dyfit) nogil



cdef extern from "wcs_core.h":
    ctypedef struct wcs_fit_input:
        double ra0, de0, roll0
        int type, order
        int n_k, n_p
        double k[4], p[4]

    ctypedef struct wcs_fit_output:
        double crpix1, crpix2
        double cd11, cd12, cd21, cd22
        double res_prj, res_pix
        int order, nvar
        double *pixpoly_dx
        double *pixpoly_dy
        double *prjpoly_dx
        double *prjpoly_dy

    int wcs_fit_cy(
        double *ra, double *dec, double *px, double *py, int npts,
        const wcs_fit_input *cfg,
        wcs_fit_output *out) nogil

    void wcs_fit_output_free(wcs_fit_output *out) nogil



cdef extern from "polyfit_core.h":
    int poly_fit_cy(
        double *x, double *y, double *vals, double *weights, int npts,
        int order, double ox, double oy, double scale,
        double *coeff_out) nogil

    int poly_compose_affine_cy(
        double *coeff, int order,
        double xlin0, double xlin1, double xlin2,
        double ylin0, double ylin1, double ylin2,
        double *rc) nogil



