/*****************************************************************************/
/* image.h — simple 2D image struct, no fits library dependency              */
/*****************************************************************************/

#ifndef	__IMAGE_H_INCLUDED
#define	__IMAGE_H_INCLUDED	1

#include <stdio.h>

/*****************************************************************************/

#define	FITS_TAPE_CARDIMAGESIZE		80

#define	FITS_VSTR	1
#define	FITS_VINT	3
#define	FITS_VBOOLEAN	2
#define	FITS_VDOUBLE	4
#define	FITS_VCOMMENT	5

#define	FITS_SH_FIRST		0
#define	FITS_SH_LAST		1
#define	FITS_SH_ADD		2
#define	FITS_SH_INSERT		3

#define	FITS_EXT_IMAGE		1

/*****************************************************************************/

typedef struct {
    int		sx, sy;
    int		bit;
    double	**data;
    int		dim;
    int		naxis[17];
    void	*vdata;
    void	*allocdata;
    double	curr_bscale, curr_bzero;
    double	read_bscale, read_bzero;
} image;

/*****************************************************************************/

typedef struct {
    char	name[FITS_TAPE_CARDIMAGESIZE];
    char	comment[FITS_TAPE_CARDIMAGESIZE];
    char	vstr[FITS_TAPE_CARDIMAGESIZE];
    int		vtype;
    int		vint;
    double	vdouble;
} fitsheader;

typedef struct {
    fitsheader	*hdrs;
    int		nhdr;
    int		ahdr;
} fitsheaderset;

/*****************************************************************************/

typedef union {
    image		i;
    struct { int _pad[32]; } t;
    struct { int _pad[32]; } b;
} fitsextensiondata;

typedef struct {
    int			type;
    fitsheaderset	header;
    fitsextensiondata	x;
} fitsextension;

typedef struct {
    fitsheaderset	header;
    image		i;
    fitsextension	*xtns;
    int			nxtn;
    int			length;
    char		*rawdata;
} fits;

/*****************************************************************************/

int	image_duplicate(image *ret, image *src, int flag);
void	image_free(image *img);
int	image_draw_line(image *img, int x, int y, int dx, int dy, double col, int style);

/*****************************************************************************/
/* stubs for fits_headerset_* and fits I/O functions (provided by fits_stubs.c) */

int	fits_headerset_get_count(fitsheaderset *header, char *hdr);
fitsheader *fits_headerset_get_header(fitsheaderset *header, char *hdr, int cnt);
int	fits_headerset_delete_all(fitsheaderset *header, char *hdr);
int	fits_headerset_set_string(fitsheaderset *header, char *hdr, int rule, char *str, char *comment);

fits		*fits_create(void);
void		 fits_free(fits *img);
void		 fits_set_standard(fits *img, char *comment);
int		 fits_set_image_params(fits *img);
int		 fits_write(FILE *fw, fits *img);

int fits_header_add_int(void *a, char *b, int c, char *d);
int fits_header_add_float(void *a, char *b, double c, int d, char *e);
int fits_bintable_alloc(void *a, int b);
int fits_bintable_check_fields(void *a);
int fits_table_add_column(void *a, void *b);

/*****************************************************************************/

#endif
