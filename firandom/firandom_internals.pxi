# firandom_internals.pxi — C declarations + Python helpers for firandom

cdef extern from "image.h":
    ctypedef struct image:
        int sx, sy
        int bit
        double **data
        int dim
        int naxis[17]
        void *vdata
        void *allocdata
        double curr_bscale, curr_bzero
        double read_bscale, read_bzero

    void image_free(image *img)
    int image_duplicate(image *ret, image *src, int flag)

cdef extern from "magnitude.h":
    ctypedef struct magflux:
        double magnitude
        double intensity

cdef extern from "psf.h":
    ctypedef struct psf:
        int hsize, grid, order
        double ox, oy, scale
        double ***coeff

cdef extern from "firandom.h":
    ctypedef image fitsimage

    int LISTTYPE_FEP
    int LISTTYPE_SIG
    int LISTTYPE_SDK
    int LISTTYPE_SMM
    int LISTTYPE_PSF
    int LISTTYPE_PDS

    int VAR_X, VAR_Y, VAR_AX, VAR_AY, VAR_I, VAR_M
    int VAR_SH_F, VAR_SH_E, VAR_SH_P
    int VAR_SN_S, VAR_SN_D, VAR_SN_K
    int VAR_SI_S, VAR_SI_D, VAR_SI_K

    ctypedef struct stargenparam:
        int subg
        double **subpixeldata
        double gain
        int is_photnoise, method
        int is_intinelect, dontquantize
        double nsuppress
        psf *tpd

    ctypedef struct starlistparam:
        magflux mf0
        double ox, oy, scale
        int sx, sy, basetype

    ctypedef struct inlistparam:
        int colx, coly
        int colflux, colmag
        int colsh1, colsh2, colsh3, colsh4
        int listtype
        magflux mf0

    double get_gaussian(double mean, double stddev)
    int get_gaussian_2d(double x0, double y0, double sigma, double delta, double kappa, double *rx, double *ry)

    int fep_to_sdk(double f, double e, double p, double *s, double *d, double *k)
    int sdk_to_fep(double s, double d, double k, double *f, double *e, double *p)
    int sdk_to_isdk(double s, double d, double k, double *_is, double *_id, double *_ik)
    int isdk_to_sdk(double _is, double _id, double _ik, double *s, double *d, double *k)

    int create_background(fitsimage *img, char *bgarg, double stddev, double ox, double oy, double scale, int zoom)
    int replace_limiters(char *buff)
    int create_input_list(char *buff, starlistparam *lp, star **rstars, int *rnstar, int sseed)

    int draw_starlist(fitsimage *img, stargenparam *sgp, star *stars, int nstar, int zoom)
    int quantize_image(fitsimage *img)
    int divide_image(fitsimage *img, double d)

cdef extern from "random.h":
    int random_seed(int j)
    double random_double()

cdef extern from "stars.h":
    int SHAPE_GAUSS
    int SHAPE_ELLIPTIC
    int SHAPE_DEVIATED
    int SHAPE_PSF
    int MAX_DEVIATION_COEFF

    ctypedef struct starlocation:
        double gamp, gbg
        double gcx, gcy

    ctypedef struct starshape:
        int model, order
        double gs, gd, gk, gl
        double mom[15]
        double factor

    ctypedef struct candidate:
        pass

    ctypedef struct star:
        starlocation location
        starshape shape
        starlocation psf
        double gsig, gdel, gkap
        double gfwhm, gellip, gpa
        double flux
        int marked
        candidate *cand

cdef extern from "algorithms/tensor.h":
    void *tensor_alloc_2d(int typesize, int a, int b)
    int tensor_free(void *tensor)

cdef double _SIG_FWHM = 2.3548200450309493
