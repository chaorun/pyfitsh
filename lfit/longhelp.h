/*****************************************************************************/
/* longhelp.h — stripped for pyfitsh (longhelp_entry typedef is in fitsh.h) */
/*****************************************************************************/

#ifndef	__LONGHELP_H_INCLUDED
#define	__LONGHELP_H_INCLUDED	1

#include "fitsh.h"

/*****************************************************************************/

/* these are specific terms which are supported by the ``help2man'' utility. */
#define		LONGHELP_OPTIONS	{ "Options:", NULL }
#define		LONGHELP_EXAMPLES	{ "Examples:", NULL }

/*****************************************************************************/

/* longhelp_entry is defined in fitsh.h */

/*****************************************************************************/

int	longhelp_fprint(FILE *fw,longhelp_entry *entries,int flags,int width);
int	longhelp_fprint_mediawiki(FILE *fw,longhelp_entry *entry);

/*****************************************************************************/

#endif
