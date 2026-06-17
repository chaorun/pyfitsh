/*****************************************************************************/
/* maskjoin.h                                                                */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Pure in-memory mask merging (replaces CLI join_masks_from_files).         */
/*****************************************************************************/

#ifndef __MASKJOIN_H_INCLUDED
#define __MASKJOIN_H_INCLUDED   1

int mask_merge(char **dst, char **src, int sx, int sy);

#endif
