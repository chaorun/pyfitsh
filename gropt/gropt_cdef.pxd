from libc.stddef cimport size_t

cdef extern from "math/matrixvector.h":
    ctypedef double vector[3]
    ctypedef double matrix[3][3]

cdef extern from "optcalc.h":
    struct glass:
        char g_name[32]
        double g_n
        int g_nncoeff
        double g_ncoeffs[8]

    struct surface:
        double s_curvature
        double s_conic
        int s_nalpha
        double s_alphas[8]

    struct lens:
        double l_offset
        double l_thickness
        double l_radius1
        double l_radius2
        surface l_s1
        surface l_s2
        int l_index_glass
        double l_n_lambda

    struct optics:
        glass *opt_glasses
        int opt_nglass
        lens *opt_lenses
        int opt_nlens
        double opt_z_focal
        surface opt_s_focal

    struct raytrace:
        int rt_npoint
        vector *rt_points

    ctypedef double transfer_matrix[2][2]

    int optcalc_reset(optics *opt) nogil
    int optcalc_read_optics_from_string(const char *text,optics *opt) nogil
    int optcalc_ray_trace(optics *opt,double lam,vector v0,vector n0,raytrace *rt) nogil
    int optcalc_glass_refraction_precompute(optics *opt,double lam) nogil
    int optcalc_raytrace_reset(raytrace *rt) nogil
    int optcalc_raytrace_free(raytrace *rt) nogil

cdef extern from "gropt_interface.h":
    int gropt_transfer(optics *opt,double wavelength,
        double *out_focal_plane,double *out_eff_focus) nogil

    int gropt_spot_diagram(optics *opt,
        double wavelength,double aperture_radius,int nrings,
        double angle_nx,double angle_ny,double zstart,double pixel_scale,
        double *out_xy,int out_max,int *out_n,
        double *out_center_x,double *out_center_y) nogil

    int gropt_single_raytrace(optics *opt,double wavelength,
        double x0,double y0,double z0,
        double nx,double ny,double nz,
        double *out_points,int out_max,int *out_n) nogil

    int gropt_compute_psf(optics *opt,
        double wavelength,double aperture_radius,
        double angle_nx,double angle_ny,double zstart,double pixel_scale,
        int psf_hsize,double *out_psf) nogil

    int gropt_export_scad(optics *opt,
        char **out_str,size_t *out_len) nogil

    int gropt_export_eps(optics *opt,
        double aperture_radius,int nrings,
        double angle_nx,double angle_ny,double zstart,
        double wavelength,double pixel_scale,
        char **out_str,size_t *out_len) nogil
