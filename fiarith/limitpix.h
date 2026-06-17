/*****************************************************************************/
/* limitpix.h                                                                */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Marks or corrects pixels that exceed integer representation limits.       */
/* Pure in-memory operation, no file I/O.                                     */
/*****************************************************************************/

#ifndef __LIMITPIX_H_INCLUDED
#define __LIMITPIX_H_INCLUDED   1

#include "image.h"

int mark_integerlimited_pixels(image *img, char **mask, int bitpix,
                               int is_corr, int mvlo, int mvhi);

#endif
