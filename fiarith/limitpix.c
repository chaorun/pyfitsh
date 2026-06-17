/*****************************************************************************/
/* limitpix.c                                                                */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Marks or corrects pixels that are outside the valid range of the given    */
/* integer bitpix representation.  Pure in-memory, no file I/O.              */
/* Extracted from fitsh common.c:mark_integerlimited_pixels().               */
/*****************************************************************************/

#include <math.h>
#include "image.h"
#include "mask.h"
#include "limitpix.h"

int mark_integerlimited_pixels(image *img, char **mask, int bitpix,
                               int is_corr, int mvlo, int mvhi)
{
 double  llo, lhi;
 int     k, l;

 if ( img == NULL || img->sx <= 0 || img->sy <= 0 || img->data == NULL )
        return -1;

 if ( bitpix < 0 )
        return 0;       /* floating-point data, nothing to limit */

 if ( bitpix == 8 )
   {    llo = -128.0;
        lhi = +127.0;
   }
 else if ( bitpix == 16 )
   {    llo = -32768.0;
        lhi = +32767.0;
   }
 else if ( bitpix == 32 )
   {    llo = -2147483648.0;
        lhi = +2147483647.0;
   }
 else
        return -1;      /* invalid bitpix */

 for ( k = 0 ; k < img->sy ; k++ )
   {    for ( l = 0 ; l < img->sx ; l++ )
         {      if ( img->data[k][l] < llo )
                 {      if ( mask != NULL )   mask[k][l] |= mvlo;
                        if ( is_corr )        img->data[k][l] = llo;
                 }
                else if ( img->data[k][l] > lhi )
                 {      if ( mask != NULL )   mask[k][l] |= mvhi;
                        if ( is_corr )        img->data[k][l] = lhi;
                 }
                else if ( is_corr )
                        img->data[k][l] = floor(img->data[k][l]);
         }
   }

 return 0;
}
