#ifndef FISTAR_PIPELINE_H
#define FISTAR_PIPELINE_H

/* output flags */
#define FISTAR_OUT_MARK  1
#define FISTAR_OUT_AREA  2
#define FISTAR_OUT_PSF   4

/* result struct */
typedef struct {
    double *ix, *iy, *cx, *cy, *cbg, *camp, *cmax;
    int    *npix;
    double *cs, *cd, *ck;
    double *x, *y, *bg, *amp;
    double *s, *d, *k, *l;
    double *sigma, *delta, *kappa;
    double *fwhm, *ellip, *pa;
    double *flux, *noise, *sn, *magnitude;
    double *px, *py, *pbg, *pamp;
    double *ps, *pd, *pk, *pl;
} fistar_result;

/* PSF output */
typedef struct {
    int nvar, nside;
    double *data;     /* flat (nvar*nside*nside) */
    double ox, oy, scale;
    int hsize, grid, order;
} fistar_psf_out;

int fistar_search_cy(
    double *img_data, char *mask_data, int sx, int sy,
    /* detection */
    double threshold, double flux_threshold, double critical_prominence,
    double skysigma,
    int algorithm,
    /* fitting */
    int model, int model_order,
    int iter_symmetric, int iter_general,
    int is_fit_model,
    int collfit_niter, int collfit_refinelevel, double collfit_bhsize,
    /* PSF */
    int is_determine_psf, int psf_hsize, int psf_grid, int psf_order,
    int psf_type, int psf_symmetrize, int psf_use_biquad,
    double psf_integral_kappa, double psf_circle_width, int psf_circle_order,
    /* gain */
    double gain,
    double mag_intensity, double mag_magnitude,
    /* input candidates */
    double *cand_xy, int in_ncand, int cand_ncol,
    double cand_radius,
    /* input positions */
    double *pos_xy, int npos,
    /* source range */
    int src_xmin, int src_xmax, int src_ymin, int src_ymax,
    /* mark drawing */
    int mark_symbol, int mark_size,
    /* output control */
    int out_flags,
    int sort,
    int is_verbose,
    /* results */
    int *nstar_out,
    fistar_result *result,
    /* mark image */
    double *mark_data,
    /* area image */
    double *area_data,
    /* PSF output */
    fistar_psf_out *psf_out,
    /* position flux output */
    double *pos_flux, double *pos_ferr);
#endif
