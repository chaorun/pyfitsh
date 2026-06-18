/*****************************************************************************/
/* fiign_core.h — memory-based pixel mask operations (no FITS I/O)           */
/*****************************************************************************/

#ifndef __FIIGN_CORE_H_INCLUDED
#define __FIIGN_CORE_H_INCLUDED  1

/*****************************************************************************/

#include "image.h"
#include "mask.h"

/*****************************************************************************/

#define FIIGN_VERSION  "0.9e"

/*****************************************************************************/

typedef struct {
    double **data;
    char   **mask;
    int    sx, sy;
    char   *history;
} fiign_result;

/*****************************************************************************/

fiign_result *fiign_apply(
    double **img, char **mask, int sx, int sy,
    double saturation, double *sat_img,
    int leak_method,
    int ignore_nonpos, int ignore_neg, int ignore_zero,
    int ignore_cosmics, int replace_cosmics,
    double th_low, double th_high, double sky_sigma,
    int expand_hsize,
    int apply_mask, double mask_value,
    int bitpix,
    char **convert_list,
    char **mask_block_list
);

void fiign_result_free(fiign_result *r);

int parse_mask_flags_simple(char *list);

/*****************************************************************************/

#endif
