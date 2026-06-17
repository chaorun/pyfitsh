# fistar C API declarations
cdef extern from "fistar_core.h":
    ctypedef struct fistar_result:
        double *ix, *iy, *cx, *cy, *cbg, *camp, *cmax
        int    *npix
        double *cs, *cd, *ck
        double *x, *y, *bg, *amp
        double *s, *d, *k, *l
        double *sigma, *delta, *kappa
        double *fwhm, *ellip, *pa
        double *flux, *noise, *sn, *magnitude
        double *px, *py, *pbg, *pamp
        double *ps, *pd, *pk, *pl

    ctypedef struct fistar_psf_out:
        int nvar, nside
        double *data
        double ox, oy, scale
        int hsize, grid, order

    int fistar_search_cy(
        double *img, char *mask, int sx, int sy,
        double th, double fth, double cp, double skysigma,
        int alg, int model, int morder,
        int it_sym, int it_gen, int only_cand,
        int cf_niter, int cf_refine, double cf_bhsize,
        int psf, int psf_hsize, int psf_grid, int psf_order,
        int psf_type, int psf_symm, int psf_biquad,
        double psf_ikappa, double psf_cwidth, int psf_corder,
        double gain,
        double mag_intensity, double mag_magnitude,
        double *cand_xy, int ncand, int cand_ncol, double cand_rad,
        double *pos_xy, int npos,
        int src_xmin, int src_xmax, int src_ymin, int src_ymax,
        int mark_symbol, int mark_size,
        int out_flags, int sort,
        int is_verbose,
        int *nstar, fistar_result *res,
        double *mark, double *area,
        fistar_psf_out *psf_out,
        double *pos_flux, double *pos_ferr) nogil


