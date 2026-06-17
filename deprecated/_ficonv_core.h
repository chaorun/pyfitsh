#ifndef _FICONV_CORE_H
#define _FICONV_CORE_H

int ficonv_fit_flat(
    double *ref_data, char *ref_mask,
    double *img_data, char *img_mask,
    int sx, int sy,
    char *kernel_spec,
    int niter, double rejlevel, double gain,
    double *cnv_data, char *cnv_mask);

#endif
