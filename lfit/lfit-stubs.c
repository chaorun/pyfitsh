/*****************************************************************************/
/* lfit-stubs.c — stubs for CLI/I/O functions used by lfit (no-IO port)      */
/*****************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "fitsh.h"
#include "longhelp.h"

/*****************************************************************************/

/* --- from iof.h (file I/O) --- */

char *freadline(FILE *fw)
{
    (void)fw;
    return NULL;
}

FILE *fopenread(const char *name)
{
    (void)name;
    return NULL;
}

void fclosewrite(FILE *fp)
{
    if (fp && fp != stdout && fp != stderr) fclose(fp);
}

int fcloseread(FILE *fp)
{
    if (fp && fp != stdin) fclose(fp);
    return 0;
}

FILE *fopenwrite(const char *name)
{
    (void)name;
    return NULL;
}

int is_fits(const char *name)
{
    (void)name;
    return 0;
}

/* --- from scanarg.h (CLI parsing) --- */

#define SCANARG_ALLOW_FLAGS  0x01

int scanarg(int argc, char *argv[], int flags, ...)
{
    (void)argc; (void)argv; (void)flags;
    return -1;
}

#define SCANPAR_DEFAULT  0x00

int scanpar(char *par, int flags, ...)
{
    (void)par; (void)flags;
    return 0;
}

/* --- from ui.c (help/version) --- */

int fprint_generic_long_help(FILE *fw, int is_wiki, longhelp_entry *help,
    char *synopsis, char *description)
{
    (void)fw; (void)is_wiki; (void)help;
    (void)synopsis; (void)description;
    return 0;
}

int fprint_generic_version(FILE *fw, char *arg0, char *name, char *pv, int type)
{
    (void)fw; (void)arg0; (void)name; (void)pv; (void)type;
    return 0;
}

/* --- from scanarg.h (scanflag) --- */

int scanflag(char *par, int flags, ...)
{
    (void)par; (void)flags;
    return 0;
}

/* --- from fitsmask.h (mask functions) --- */

char **fits_mask_create_empty(int sx, int sy)
{
    (void)sx; (void)sy;
    return NULL;
}

void fits_mask_free(char **mask)
{
    (void)mask;
}

int fits_mask_expand_logic(char **mask, int sx, int sy, int dis,
    int imv, int smv, int expand_border)
{
    (void)mask; (void)sx; (void)sy; (void)dis;
    (void)imv; (void)smv; (void)expand_border;
    return 0;
}

int fits_mask_export_as_header(void *header, int cl, char **mask,
    int sx, int sy, char *hdr)
{
    (void)header; (void)cl; (void)mask;
    (void)sx; (void)sy; (void)hdr;
    return 0;
}

int fits_mask_mark_nans(void *img, char **mask, int setmask)
{
    (void)img; (void)mask; (void)setmask;
    return 0;
}

int join_masks_from_files(char **mask, int sx, int sy, char **inmasklist)
{
    (void)mask; (void)sx; (void)sy; (void)inmasklist;
    return 0;
}

int parse_mask_flags(char *list)
{
    (void)list;
    return 0;
}

/* --- from history.h --- */

void fits_history_export_command_line(void *img, char *prg, char *vrs,
    int argc, char *argv[])
{
    (void)img; (void)prg; (void)vrs;
    (void)argc; (void)argv;
}

/*****************************************************************************/
