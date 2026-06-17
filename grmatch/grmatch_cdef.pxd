# grmatch C API declarations
cdef extern from "grmatch_lib.h":
    ctypedef struct matchpointstat:
        double wsigma
        double nsigma
        double unitarity
        double time_total
        double time_trimatch
        double time_symmatch
        int nmiter
        int tri_level
        double hull_coverage



cdef extern from "grmatch_core.h":
    int do_pointmatch_cy(
        double *ref_x, double *ref_y, double *ref_ord, int nref,
        double *inp_x, double *inp_y, double *inp_ord, int ninp,
        int ttype, double maxdist, double unitarity, int parity,
        int use_ordering, int maxnum_ref, int maxnum_inp,
        int nmiter, int friter, double rejlevel,
        int wcat, int w_magnitude, double wpower,
        int is_centering, double refcx, double refcy,
        double inpcx, double inpcy, double maxcenterdist,
        int order,
        double ox, double oy, double scale,
        double *hint_dxfit, double *hint_dyfit, int hint_order,
        int **hits_idx0, int **hits_idx1, int *rnhit,
        double **vfits_dx, double **vfits_dy,
        matchpointstat *mps) nogil

    int coord_match_cy(
        double *ref_x, double *ref_y, int nref,
        double *inp_x, double *inp_y, int ninp,
        double maxdist,
        int **hits_idx0, int **hits_idx1, int *rnhit) nogil

    int id_match_cy(
        char **ref_ids, int nref,
        char **inp_ids, int ninp,
        int ambig,
        int **hits_idx0, int **hits_idx1, int *rnhit) nogil



