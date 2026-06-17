/*****************************************************************************/
/* _ficonv_core.c — thin flat wrapper, uses kernel-base.c + parser.        */
/*****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fits/fits.h>
#include "kernel.h"
#include "_ficonv_core.h"

extern int create_kernels_from_kernelarg(char *, kernellist *);
extern int fit_kernel_poly_coefficients_block(fitsimage *, fitsimage *,
    char **, double **, int, int, kernellist *, kernellist *);
extern int convolve_with_kernel_set(fitsimage *, char **, kernellist *, fitsimage *);

extern void fitrans_profile_dump(void);

int ficonv_fit_flat(
    double *ref_data, char *ref_mask,
    double *img_data, char *img_mask,
    int sx, int sy,
    char *kernel_spec,
    int niter, double rejlevel, double gain,
    double *cnv_data, char *cnv_mask)
{
    int i;
    fitsimage refimg, img, outimg;
    kernellist kl, xl;
    char **mr, **mi, **mo;

    memset(&kl, 0, sizeof(kl)); memset(&xl, 0, sizeof(xl));
    if (create_kernels_from_kernelarg(kernel_spec, &kl)) return 1;
    /* set offset/scale for polynomial stability (same as original ficonv) */
    kl.ox = 0.5 * (double)sx;
    kl.oy = 0.5 * (double)sy;
    kl.scale = 0.5 * (double)sx;
    xl.type = 1; xl.kernels = kl.kernels; xl.nkernel = kl.nkernel;

    memset(&refimg, 0, sizeof(refimg)); refimg.sx = sx; refimg.sy = sy;
    refimg.data = (double **)malloc(sy * sizeof(double *));
    for (i = 0; i < sy; i++) refimg.data[i] = ref_data + i * sx;
    memset(&img, 0, sizeof(img)); img.sx = sx; img.sy = sy;
    img.data = (double **)malloc(sy * sizeof(double *));
    for (i = 0; i < sy; i++) img.data[i] = img_data + i * sx;
    memset(&outimg, 0, sizeof(outimg)); outimg.sx = sx; outimg.sy = sy;
    outimg.data = (double **)malloc(sy * sizeof(double *));
    for (i = 0; i < sy; i++) outimg.data[i] = cnv_data + i * sx;

    mr = (char **)malloc(sy * sizeof(char *));
    for (i = 0; i < sy; i++) mr[i] = ref_mask + i * sx;
    mi = (char **)malloc(sy * sizeof(char *));
    for (i = 0; i < sy; i++) mi[i] = img_mask + i * sx;
    mo = (char **)malloc(sy * sizeof(char *));
    for (i = 0; i < sy; i++) mo[i] = cnv_mask + i * sx;

    { kernel *k; for (i = 0, k = kl.kernels; i < kl.nkernel; i++, k++) k->flag = 1;
      for (i = 0, k = xl.kernels; i < xl.nkernel; i++, k++) k->flag = 1; }
    fit_kernel_poly_coefficients_block(&refimg, &img, mr, NULL, 1, 1, &kl, &xl);
    kl.type = 1;
    convolve_with_kernel_set(&refimg, mo, &kl, &outimg);

    free(refimg.data); free(img.data); free(outimg.data);
    free(mr); free(mi); free(mo);
    return 0;
}
