/* minimal iof.h stub for pyfitsh lfit port */
#ifndef __IOF_H_STUB
#define __IOF_H_STUB
#include <stdio.h>
char *freadline(FILE *fw);
FILE *fopenread(const char *name);
void fclosewrite(FILE *fp);
int  fcloseread(FILE *fp);
FILE *fopenwrite(const char *name);
int is_fits(const char *name);
#endif
