/*****************************************************************************/
/* fits_stubs.c — stubs for fits functions called by code that must compile  */
/* but whose callers are unused in the flat-pipeline context.               */
/*****************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "fitsh.h"
#include "image.h"

/* --- logmsg (used by ficonv_core.c, kernel-base.c) --- */

int logmsg(int flag, char *msg, ...)
{
    if (flag) {
        va_list ap;
        va_start(ap, msg);
        vfprintf(stderr, msg, ap);
        va_end(ap);
        fflush(stderr);
        return 0;
    }
    return 1;
}

/* --- header I/O stubs (called from fitsmask.c dead functions) --- */

int fits_headerset_get_count(fitsheaderset *header, char *hdr)
{ (void)header; (void)hdr; return 0; }

fitsheader *fits_headerset_get_header(fitsheaderset *header, char *hdr, int cnt)
{ (void)header; (void)hdr; (void)cnt; return NULL; }

int fits_headerset_delete_all(fitsheaderset *header, char *hdr)
{ (void)header; (void)hdr; return 0; }

int fits_headerset_set_string(fitsheaderset *header, char *hdr, int rule, char *str, char *comment)
{ (void)header; (void)hdr; (void)rule; (void)str; (void)comment; return 0; }

/* --- kernel I/O stubs (called from kernel-io.c dump_kernel, unused) --- */

fits *fits_create(void)
{ return NULL; }

void fits_free(fits *img)
{ if (img) free(img); }

void fits_set_standard(fits *img, char *comment)
{ (void)img; (void)comment; }

int fits_set_image_params(fits *img)
{ (void)img; return 0; }

int fits_write(FILE *fw, fits *img)
{ (void)fw; (void)img; return 0; }

/* --- table/header stubs (from previous fits_stubs.c) --- */

int fits_header_add_int(void *a, char *b, int c, char *d)
{ (void)a; (void)b; (void)c; (void)d; return 0; }

int fits_header_add_float(void *a, char *b, double c, int d, char *e)
{ (void)a; (void)b; (void)c; (void)d; (void)e; return 0; }

int fits_bintable_alloc(void *a, int b)
{ (void)a; (void)b; return 0; }

int fits_bintable_check_fields(void *a)
{ (void)a; return 0; }

int fits_table_add_column(void *a, void *b)
{ (void)a; (void)b; return 0; }

/* --- fprint_error / fprint_warning (common to lfit, fiphot, fistar, ficalib, grmatch) --- */

static char *progbasename = "pyfitsh";

int fprint_error(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: error: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

int fprint_warning(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: warning: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}
