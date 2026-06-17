# fitrans C API declarations
cdef extern from "fitrans_core.h":
    int fitrans_apply_cy(
        double *in_data, char *in_mask, int sx, int sy,
        double *out_data, char *out_mask, int nsx, int nsy,
        int ofx, int ofy,
        int order, int int_method, int is_invert,
        double ox, double oy, double scale,
        double *dxfit, double *dyfit) nogil



cdef extern from "fitrans_ops.h":
    int fitrans_noise_cy(
        double *in_data, char *in_mask, int sx, int sy,
        double *out_data, char *out_mask) nogil

    int fitrans_zoom_cy(
        double *in_data, char *in_mask, int sx, int sy,
        double *out_data, char *out_mask, int nsx, int nsy,
        int ofx, int ofy, int scalex, int scaley, int is_raw) nogil

    int fitrans_shrink_cy(
        double *in_data, char *in_mask, int sx, int sy,
        double *out_data, char *out_mask, int nsx, int nsy,
        int ofx, int ofy, int scalex, int scaley,
        int mode_median, int mode_truncated_mean, int mode_average_mask) nogil

    int fitrans_smooth_cy(
        double *in_data, char *in_mask, int sx, int sy,
        double *out_data, char *out_mask,
        int smooth_type, int xorder, int yorder,
        int prefilter, int fxhsize, int fyhsize,
        double frejratio, int niter, double lower, double upper,
        int is_mean_unity, int is_detrend) nogil


    int fitrans_repetitive_cy(
        double *in_data, char *in_mask, int sx, int sy,
        double *out_data, char *out_mask, int nsx, int nsy,
        int ofx, int ofy) nogil

    int fitrans_interleave_cy(
        double *in_data, char *in_mask, int sx, int sy,
        double *out_data, char *out_mask, int nsx, int nsy,
        int ofx, int ofy,
        int mode_median, int mode_truncated_mean, int mode_average_mask) nogil



