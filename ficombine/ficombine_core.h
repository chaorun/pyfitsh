/*****************************************************************************/
/* ficombine_core.h — memory-based image combination (no FITS I/O)           */
/*****************************************************************************/

#ifndef __FICOMBINE_CORE_H_INCLUDED
#define __FICOMBINE_CORE_H_INCLUDED  1

/*****************************************************************************/

#include "image.h"
#include "mask.h"
#include "combine.h"

/*****************************************************************************/

#define FICOMBINE_VERSION  "0.9"

/*****************************************************************************/

typedef struct {
    double **data;
    char   **mask;
    int    sx, sy;
    char   *history;
} ficombine_result;

/*****************************************************************************/

int combine_parse_mode_simple(char *modstr, combine_parameters *cp);

ficombine_result *ficombine_apply(
    double **images, char **masks,
    int nimg, int sx, int sy,
    combine_parameters *cp,
    int apply_mask
);

void ficombine_result_free(ficombine_result *r);

/*****************************************************************************/

#endif
