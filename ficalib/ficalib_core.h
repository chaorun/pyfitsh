/*****************************************************************************/
/* ficalib_core.h — pyfitsh calibration interface                           */
/*****************************************************************************/

#ifndef __FICALIB_CORE_H_INCLUDED
#define __FICALIB_CORE_H_INCLUDED   1

#include "image.h"

typedef struct {
    double *coeff;
    int order;
    double resd;
    double vmin;
} ficalib_flatpoly;

typedef struct {
    double *data;
    int nscan;
    /* overscan params returned */
} ficalib_overscan;

/* Main calibration entry point.
 * All input images are numpy float64 2D arrays passed as flat pointers.
 * mask must always be provided (use all-zero if no masking needed).
 * Returns 0 on success, stores output in *out_data / *out_mask / *out_flat. */
int ficalib_calibrate_cy(
    double *img_data, unsigned char *mask_data, int sx, int sy,
    double *bias_data, unsigned char *bias_mask,
    int bias_sx, int bias_sy,
    double *dark_data, unsigned char *dark_mask,
    int dark_sx, int dark_sy,
    double *flat_data, unsigned char *flat_mask,
    int flat_sx, int flat_sy,
    double flat_mean, double dark_time,
    /* output */
    double **out_data, unsigned char **out_mask, int *out_sx, int *out_sy,
    ficalib_flatpoly *out_flat);

#endif
