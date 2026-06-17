/*****************************************************************************/
/* stars.h								     */
/*****************************************************************************/

#ifndef	__STARS_H_INCLUDED
#define	__STARS_H_INCLUDED	1

#include "image.h"
#include "background.h"

typedef struct
 {	int	xmin,xmax;
	int	ymin,ymax;
 } range;

#ifndef	__IPOINT_STRUCT_DEFINED
#define	__IPOINT_STRUCT_DEFINED
typedef struct
 {	int	x,y;
 } ipoint;
#endif

typedef struct
 {	ipoint	*ipoints;
	int	nipoint;
 } ipointlist;

typedef struct
 {	ipoint	*ipoints;
	double	*yvals;
	int	nipoint;
 } mpointlist;

typedef struct
 {	int	ix,iy;			/* Integer coordinates of a star can-*/
					/* didate, the image has a local max-*/
					/* imum in the point (ix,iy), or     */
					/* something like that...	     */

	double	cx,cy;			/* A better approximation for the    */
					/* center of a star candidate, deri- */
					/* ved by a second-order polynomial  */
					/* fit (or something else, see also  */
					/* refine_candidate_params()...	     */

	double	peak,amp,bg;			
	double	sxx,syy,sxy;		/* Other parameters derived from this*/
					/* second-order polynomial fit.	     */
					/* All of these parameters are set   */
					/* by search_star_candidates() and   */
					/* some of them are used by 	     */
					/* markout_stars() and 		     */
					/* fit_gaussians().		     */

	ipoint	*ipoints;		/* ipoints[] contains nipoint	     */
	int	nipoint;		/* elements, (x,y) pairs, which	     */
					/* possible belong to the star.      */
					/* Determined by markout_*() and     */
					/* used by fit_gaussians().	     */

	double	area,			/* area (it is always an integer)    */
		noise,			/* bg noise			     */
		flux;			/* flux				     */

	int	flags,marked;

 } candidate;

#include "psf.h"

#define		SHAPE_GAUSS		1
#define		SHAPE_ELLIPTIC		2
#define		SHAPE_DEVIATED		3
#define		SHAPE_PSF		4

#define		MAX_DEVIATION_ORDER	4
#define		MAX_DEVIATION_COEFF	15	/* (MDO+1)*(MDO+2)/2	*/

/* #define	STAR_MULTIMODEL		1 */	/* obsoleted		*/

typedef struct
 {	double	gamp,gbg;
	double	gcx,gcy;
 } starlocation;

typedef struct
 {	int	model,order;	/* `model` can be: SHAPE_GAUSS, ELLIPTIC,    */
				/* DEVIATED. `order` is only defined	     */
				/* for SHAPE_DEVIATED (act. between 2 and 4) */
				/* PSF fit results are stored in another     */
				/* form (see (starlocation)star->psf).	     */

	double	gs,gd,gk,gl;	/* coeff's for SHAPE_ELLIPTIC (gs, gd, gk),  */
				/* DEVIATED (gs) and PSF (gs, gd, gk and gl).*/

	double	mom[MAX_DEVIATION_COEFF];	/* coeff's for SHAPE_DEVIATED*/

	double	factor;		/* factor of multi-model fitting (should not */
				/* really expand out from the interval [0,1])*/

 } starshape;

typedef struct
 {	starlocation	location;	/* An approximation for the center,  */
					/* background and amplitude of the   */
					/* star	(`gamp` is model-dependent). */

	starshape	shape;		/* Shape parameters of the star.     */

	starlocation	psf;		/* PSF fitting yields only info's    */
					/* like this: centroid coordinates,  */
					/* background and amplitude (later is*/
					/* equivalent to the total flux...). */

	double		gsig,gdel,gkap;	/* Derived parameters from the Gau-  */
	double		gfwhm,		/* sian fit (FWHM, ellpticity and    */
			gellip,		/* position angle), also set by      */
			gpa;		/* fit_gaussians().		     */

	double		flux;		/* The total flux of the star.	     */
					/* Set by fit functions, and used    */
					/* by firandom also.		     */

	int		marked;		/* A flag for cleanup_starlist().    */

	candidate	*cand;		/* candidate pointer (if available). */
 } star;

/* temporary removed from typedef struct { ... } star; */
#ifdef	STAR_MULTIMODEL
	starshape	*mshapes;	/* Shape parameters of the star if   */
	int		nmshape;	/* multi-model fitting was used.     */
#endif

typedef struct
 {	int	ix,iy;
	double	value;
 } imgpoint;

typedef struct
 {	int	model;	/* can be SHAPE_{GAUSS,ELLIPTIC,DEVIATED}	     */
	int	order;	/* up to MAX_DEVIATION_ORDER, only for SHAPE_DEVIATED*/
