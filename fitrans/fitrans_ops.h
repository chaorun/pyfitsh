#ifndef _FITRANS_OPS_H
#define _FITRANS_OPS_H

int fitrans_noise_cy(
    double *in_data, char *in_mask, int sx, int sy,
    double *out_data, char *out_mask);

int fitrans_zoom_cy(
    double *in_data, char *in_mask, int sx, int sy,
    double *out_data, char *out_mask, int nsx, int nsy,
    int ofx, int ofy, int scalex, int scaley, int is_raw);

int fitrans_shrink_cy(
    double *in_data, char *in_mask, int sx, int sy,
    double *out_data, char *out_mask, int nsx, int nsy,
    int ofx, int ofy, int scalex, int scaley,
    int mode_median, int mode_truncated_mean, int mode_average_mask);

int fitrans_smooth_cy(
    double *in_data, char *in_mask, int sx, int sy,
    double *out_data, char *out_mask,
    int smooth_type, int xorder, int yorder,
    int prefilter, int fxhsize, int fyhsize,
    double frejratio, int niter, double lower, double upper,
    int is_mean_unity, int is_detrend);

int fitrans_repetitive_cy(
    double *in_data, char *in_mask, int sx, int sy,
    double *out_data, char *out_mask, int nsx, int nsy,
    int ofx, int ofy);

int fitrans_interleave_cy(
    double *in_data, char *in_mask, int sx, int sy,
    double *out_data, char *out_mask, int nsx, int nsy,
    int ofx, int ofy,
    int mode_median, int mode_truncated_mean, int mode_average_mask);

#endif
