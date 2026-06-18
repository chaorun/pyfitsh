/* minimal scanarg.h stub for pyfitsh lfit port */
#ifndef __SCANARG_H_STUB
#define __SCANARG_H_STUB
#define SCANARG_ALLOW_FLAGS  0x01
#define SCANPAR_DEFAULT      0x00
#define SNf(i)               "%SN" #i "f"
#define SNf_str(i)           "%SN" #i "f"
int scanarg(int argc, char *argv[], int flags, ...);
int scanpar(char *par, int flags, ...);
int scanflag(char *par, int flags, ...);
#endif
