# fiphot C API declarations
cdef extern from "weight.h":
    ctypedef struct weightlist:
        pass

cdef extern from "fiphot.h":
    ctypedef struct photstar:
        pass


cdef extern from "fiphot_core.h":
    int fiphot_photometry_cy(
        double *img_data, char *img_mask, int sx, int sy,
        double *star_x, double *star_y, int nstar,
        double *star_r0, double *star_ra, double *star_da, int nap,
        char *aperture_spec, int zoom,
        int bg_type, int bg_scatter, int bg_rejniter,
        double bg_rejlower, double bg_rejupper,
        int mask_ignore, int use_biquad, int use_sky, double sky_level,
        int is_disjoint_rings, int is_disjoint_apertures, double disjoint_radius,
        double **subpixel_data, int subg,
        int sg_order, double *sg_coeff, double sg_vmin,
        double correlation_length, double sigma,
        int is_calc_opt_apert,
        weightlist *wl, int weightusage,
        double *out_flux, double *out_fluxerr,
        double *out_bgarea, double *out_bgflux, double *out_bgmedian, double *out_bgsigma,
        double *out_cntr_x, double *out_cntr_y, double *out_cntr_width,
        double *out_cntr_w_d, double *out_cntr_w_k,
        double *out_cntr_x_err, double *out_cntr_y_err, double *out_cntr_w_err,
        int *out_flag, int *out_rtot, int *out_rbad, int *out_rign,
        int *out_atot, int *out_abad,
        double *out_optimal_r0, double *out_optimal_ra, double *out_optimal_da,
        double *out_raw) nogil

    int read_raw_photometry_cy(
        photstar **rps, int *rnp,
        double *raw_data, int nstar, int nap,
        double *ref_mag_arr, double *ref_col_arr, double *ref_err_arr) nogil

    int fiphot_photometry_from_raw_cy(
        double *raw_data, int nstar, int nap,
        double *ref_mag_arr, double *ref_col_arr, double *ref_err_arr,
        double *img_data, char *img_mask, int sx, int sy,
        int bg_type, int bg_scatter, int bg_rejniter,
        double bg_rejlower, double bg_rejupper,
        int mask_ignore, int use_biquad, int use_sky, double sky_level,
        int is_disjoint_rings, int is_disjoint_apertures, double disjoint_radius,
        double **subpixel_data, int subg,
        int sg_order, double *sg_coeff, double sg_vmin,
        double correlation_length, double sigma,
        weightlist *wl, int weightusage,
        int nkernel,
        double ox, double oy, double scale, int ktype,
        int *k_types, int *k_orders, int *k_ncoeffs,
        int *k_hsizes, double *k_sigmas,
        int *k_bx, int *k_by,
        double **k_coeffs,
        char *kernel_spec, int normalize_kernel,
        double *out_flux, double *out_fluxerr,
        double *out_bgarea, double *out_bgflux, double *out_bgmedian, double *out_bgsigma,
        double *out_cntr_x, double *out_cntr_y, double *out_cntr_width,
        double *out_cntr_w_d, double *out_cntr_w_k,
        double *out_cntr_x_err, double *out_cntr_y_err, double *out_cntr_w_err,
        int *out_flag, int *out_rtot, int *out_rbad, int *out_rign,
        int *out_atot, int *out_abad) nogil

    int fiphot_subtracted_photometry_cy(
        double *img_data, char *img_mask, int sx, int sy,
        double *star_x, double *star_y, int nstar,
        double *star_r0, double *star_ra, double *star_da, int nap,
        char *aperture_spec, int zoom,
        int bg_type, int bg_scatter, int bg_rejniter,
        double bg_rejlower, double bg_rejupper,
        int mask_ignore, int use_biquad, int use_sky, double sky_level,
        int is_disjoint_rings, int is_disjoint_apertures, double disjoint_radius,
        double **subpixel_data, int subg,
        int sg_order, double *sg_coeff, double sg_vmin,
        double correlation_length,
        char *kernel_spec, int normalize_kernel,
        weightlist *wl, int weightusage,
        double *out_flux, double *out_fluxerr,
        double *out_bgarea, double *out_bgflux, double *out_bgmedian, double *out_bgsigma,
        double *out_cntr_x, double *out_cntr_y, double *out_cntr_width,
        double *out_cntr_w_d, double *out_cntr_w_k,
        double *out_cntr_x_err, double *out_cntr_y_err, double *out_cntr_w_err,
        int *out_flag, int *out_rtot, int *out_rbad, int *out_rign,
        int *out_atot, int *out_abad) nogil

    int fiphot_magnitude_fit_cy(
        double *flux, int *flag, int nstar, int nap,
        double *ref_flux, double *ref_col_arr, double *ref_mag_arr, double *ref_err_arr,
        int *orders, int norder, int niter, double sigma,
        int sx, int sy,
        double *out_mag,
        int *out_ninit, int *out_nrejs, int *out_nstar, int *out_naperture) nogil


