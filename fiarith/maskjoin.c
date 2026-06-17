/*****************************************************************************/
/* maskjoin.c                                                                */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Merge one mask into another: dst[i][j] |= src[i][j] for all pixels.       */
/* Pure in-memory operation.  Replaces the file-I/O-based                     */
/* join_masks_from_files() from the original CLI.                            */
/*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "maskjoin.h"

int mask_merge(char **dst, char **src, int sx, int sy)
{
 int i;

 if ( dst == NULL || src == NULL || sx <= 0 || sy <= 0 )
        return -1;

 for ( i = 0 ; i < sy ; i++ )
   {    int j;
        for ( j = 0 ; j < sx ; j++ )
                dst[i][j] |= src[i][j];
   }

 return 0;
}
