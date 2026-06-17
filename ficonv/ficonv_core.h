#ifndef FICONV_PIPELINE_H
#define FICONV_PIPELINE_H

#include "kernel.h"

int ficonv_fit_cy(
    double *ref_data, char *ref_mask,
    double *img_data, char *img_mask,
    char *inmask_data,
    int sx, int sy,
    char *kernel_spec,
    int method,
    int bdc,
    int niter, double rejlevel, double gain,
    int is_verbose,
    double *cnv_data, char *cnv_mask,
    double *sub_data,
    double *add_data,
    int unity_kernels,
    char **kernel_list_out,
    char *stamp_arg,
    int psx, int psy, char *spline_stamp_arg,
    /* pre-fitted kernel dict (skip fitting if nkernels > 0) */
    int prefit_nkernels,
    int *prefit_types,
    int *prefit_orders,
    int *prefit_ncoeffs,
    double *prefit_coeffs,
    double prefit_ox, double prefit_oy, double prefit_scale);

int fitsh_ficonv_fit_cy(
    double *ref_data, char *ref_mask,
    double *img_data, char *img_mask,
    char *inmask_data,
    int sx, int sy,
    char *kernel_spec,
    int method, int bdc,
    int niter, double rejlevel, double gain,
    int is_verbose,
    double *cnv_data, char *cnv_mask,
    double *sub_data,
    double *add_data,
    int unity_kernels,
    int max_kernels,
    int *out_nkernel,
    double *out_ox, double *out_oy, double *out_scale,
    int *out_type, int *out_ktotal,
    double *out_kcoeffs,
    int *out_k_types, int *out_k_orders, int *out_k_ncoeffs,
    int *out_k_hsizes, double *out_k_sigmas,
    int *out_k_bx, int *out_k_by,
    char *stamp_arg,
    int psx, int psy, char *spline_stamp_arg,
    int prefit_nkernels,
    int *prefit_types,
    int *prefit_orders,
    int *prefit_ncoeffs,
    double *prefit_coeffs,
    int *prefit_hsizes,
    double *prefit_sigmas,
    int *prefit_bx,
    int *prefit_by,
    int prefit_ktype,
    double prefit_ox, double prefit_oy, double prefit_scale);

int kernel_info_read_dicts_in_ficonv(
    kernellist *kl,
    int nkernel,
    double ox, double oy, double scale, int ktype,
    int *k_types, int *k_orders, int *k_ncoeffs,
    int *k_hsizes, double *k_sigmas,
    int *k_bx, int *k_by,
    double **k_coeffs);

#endif
