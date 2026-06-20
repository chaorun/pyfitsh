/*****************************************************************************/
/* gropt_interface.h                                                         */
/* C interface for gropt geometrical optics (no FILE I/O)                    */
/*****************************************************************************/

#ifndef GROPT_INTERFACE_H
#define GROPT_INTERFACE_H

#include "optcalc.h"
#include <stddef.h>

int gropt_transfer(struct optics *opt,double wavelength,
	double *out_focal_plane,double *out_eff_focus);

int gropt_spot_diagram(struct optics *opt,
	double wavelength,double aperture_radius,int nrings,
	double angle_nx,double angle_ny,double zstart,double pixel_scale,
	double *out_xy,int out_max,int *out_n,
	double *out_center_x,double *out_center_y);

int gropt_single_raytrace(struct optics *opt,double wavelength,
	double x0,double y0,double z0,
	double nx,double ny,double nz,
	double *out_points,int out_max,int *out_n);

int gropt_compute_psf(struct optics *opt,
	double wavelength,double aperture_radius,
	double angle_nx,double angle_ny,double zstart,double pixel_scale,
	int psf_hsize,double *out_psf);

int gropt_export_scad(struct optics *opt,
	char **out_str,size_t *out_len);

int gropt_export_eps(struct optics *opt,
	double aperture_radius,int nrings,
	double angle_nx,double angle_ny,double zstart,
	double wavelength,double pixel_scale,
	char **out_str,size_t *out_len);

#endif
