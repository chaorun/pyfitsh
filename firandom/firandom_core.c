/*****************************************************************************/
/* firandom.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line utility to create artificial images.			     */
/*****************************************************************************/
#define	FITSH_FIRANDOM_VERSION	"1.4.0"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <stdarg.h>
#include <sys/time.h>

#include "random.h"
#include "fitsh.h"

#include "tensor.h"
#include "common.h"
#include "stars.h"
#include "psf.h"
#include "magnitude.h"

#include "firandom.h"

#define		FIRANDOM_DEFAULT_NOISE_SUPPRESSION	10000.0

#ifdef  HAVE_NO_CC_EXTENSION 
#define __extension__ 
#endif 

/*****************************************************************************/

// CLI globals removed

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// fprint_error/warning removed



/*****************************************************************************/

/* gain: fotons/adu, usually > 1... */
int draw_star_montecarlo_gauss(fitsimage *img,stargenparam *sgp,double x0,double y0,double dflux,double is,double id,double ik)
{
 int	sx,sy,ix,iy,subx,suby,flux;
 double	x,y,fadu;

 if ( img==NULL )	return(-1);
 if ( img->data==NULL )	return(-1);
 sx=img->sx,sy=img->sy;

 if ( sgp->is_intinelect )	flux=(int)dflux;
 else				flux=(int)(dflux*sgp->gain+0.5);

 if ( flux<=0 )		return(1);
 fadu=1.0/sgp->gain;

 if ( sgp->subpixeldata==NULL || sgp->subg==0 )
  {	for ( ; flux>0 ; flux-- )
	 {	get_gaussian_2d(x0,y0,is,id,ik,&x,&y);
		ix=(int)floor(x),
		iy=(int)floor(y);
		if ( ix<0 || iy<0 || ix>=sx || iy>=sy )	continue;
		img->data[iy][ix]+=fadu;
	 }
  }
 else
  {	for ( ; flux>0 ; flux-- )
	 {	get_gaussian_2d(x0,y0,is,id,ik,&x,&y);
		ix=(int)floor(x),
		iy=(int)floor(y);
		if ( ix<0 || iy<0 || ix>=sx || iy>=sy )	continue;
		subx=(int)((x-(double)ix+2)*(double)sgp->subg)%sgp->subg;
		suby=(int)((y-(double)iy+2)*(double)sgp->subg)%sgp->subg;
		img->data[iy][ix]+=sgp->subpixeldata[suby][subx]*fadu;
 	 }
  }
 return(0);
}

int draw_star_from_array(fitsimage *img,int ix0,int iy0,double **iarr,int hsize,int sg,double dflux,double **subpix,double gain)
{
 int	i,j,k,l,fsize,is_subgrid,sx,sy,si,sj,is_photnoise;
 double	sum,w,sw;

 fsize=(2*hsize+1)*sg;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);

 sum=0.0;
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
 	 {	sum+=iarr[i][j];		}
  }
 w=dflux/sum;
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
	 {	iarr[i][j]*=w;		}
  }

 if ( sg>1 && subpix!=NULL )	is_subgrid=1;
 else				is_subgrid=0;

 if ( gain > 0.0 )	is_photnoise=1;
 else			is_photnoise=0;

 for ( i=iy0-hsize,si=0 ; i<=iy0+hsize ; i++,si+=sg )
  {	if ( i<0 || i>=sy )	continue;
	for ( j=ix0-hsize,sj=0 ; j<=ix0+hsize ; j++,sj+=sg )
	 {	if ( j<0 || j>=sx )	continue;
		if ( ! is_subgrid )	w=iarr[si][sj];
		else
		 {	w=0.0;
			for ( k=0 ; k<sg ; k++ )
			 {	for ( l=0 ; l<sg ; l++ )
				 {	sw=subpix[k][l];
					w+=sw*iarr[si+k][sj+l];
				 }
			 }
		 }
		if ( is_photnoise && w>0.0 )
			w=get_gaussian(w,sqrt(w/gain));
		if ( w<=0.0 )	continue;

		img->data[i][j]+=w;
	 }
  }

 return(0);
}

int draw_star_integral_gauss(fitsimage *img,stargenparam *sgp,double x0,double y0,double dflux,double is,double id,double ik)
{
 int		hsize,grid,ix0,iy0,fsize;
 double		vgain,dgrid;
 static double	**iarr=NULL;
 static int	afsize=0;

 if ( img==NULL || img->data==NULL )	return(-1);

 if ( sgp->is_intinelect )	dflux/=sgp->gain;
 if ( sgp->is_photnoise )	vgain=sgp->gain;
 else				vgain=0.0;
 if ( dflux*sgp->gain<=0.0 )	return(0);

 hsize=(int)(is*sqrt(2.0*log(dflux*sgp->gain*sgp->nsuppress/(is*is)))+1.0);
 if ( hsize>1023 )	hsize=1023;	/* extraordinary... */

 if ( sgp->subpixeldata==NULL || sgp->subg<=0 )	grid=1;
 else						grid=sgp->subg;

 dgrid=(double)grid;
 fsize=(2*hsize+1)*grid;

 if ( fsize>afsize )
  {	afsize=fsize;
	if ( iarr != NULL )	tensor_free(iarr);
	iarr=(double **)tensor_alloc_2d(double,afsize,afsize);
  }

 ix0=(int)floor(x0),x0-=(double)ix0,
 iy0=(int)floor(y0),y0-=(double)iy0;

 star_draw_gauss(iarr,fsize,fsize,(x0+hsize)*dgrid,(y0+hsize)*dgrid,
		 is*dgrid,id*dgrid,ik*dgrid);

 draw_star_from_array(img,ix0,iy0,iarr,hsize,grid,dflux,sgp->subpixeldata,vgain);

 return(0);
}

int draw_star_integral_deviated(fitsimage *img,stargenparam *sgp,double x0,double y0,double dflux,double gs,int order,double *mom)
{
 int		hsize,grid,ix0,iy0,fsize;
 double		vgain,dgrid;
 static double	**iarr=NULL;
 static int	afsize=0;
 double		cgs,gmom[MAX_DEVIATION_COEFF],*cmom;

 if ( img==NULL || img->data==NULL )	return(-1);

 if ( sgp->is_intinelect )	dflux/=sgp->gain;
 if ( sgp->is_photnoise )	vgain=sgp->gain;
 else				vgain=0.0;
 if ( dflux*sgp->gain<=0.0 )	return(0);
 if ( gs<=0.0 )			return(0);

 ix0=(int)floor(x0),x0-=(double)ix0,
 iy0=(int)floor(y0),y0-=(double)iy0;

 hsize=(int)(sqrt(2.0*log(gs*dflux*sgp->gain*sgp->nsuppress)/gs)+1.0);

 if ( sgp->subpixeldata==NULL || sgp->subg<=0 )	grid=1;
 else						grid=sgp->subg;

 fsize=(2*hsize+1)*grid;
 dgrid=(double)grid;

 if ( fsize>afsize )
  {	afsize=fsize;
	if ( iarr  != NULL )	tensor_free(iarr);
	iarr =(double **)tensor_alloc_2d(double,afsize,afsize);
  }

 if ( grid>1 )
  {	double	igrid,w;
	int	o,l,j;

	igrid=1.0/dgrid;
	w=igrid*igrid;
	cgs=w*gs;
	for ( o=2,l=0 ; o<=order ; o++ )
	 {	for ( j=0 ; j<=o ; j++,l++ )
		 {	gmom[l]=mom[l]*w;
			w=w*igrid;
		 }
	 }
	cmom=gmom;
  }
 else
  {	cgs =gs;
	cmom=mom;
  }

 star_draw_deviated(iarr,fsize,fsize,(x0+hsize)*dgrid,(y0+hsize)*dgrid,cgs,order,cmom);

 draw_star_from_array(img,ix0,iy0,iarr,hsize,grid,dflux,sgp->subpixeldata,vgain);
 
 return(0);
}

int draw_star_integral_psf(fitsimage *img,stargenparam *sgp,double x0,double y0,
	double px,double py,
	double is,double id,double ik,double il,double dflux)
{
 static double	**iarr=NULL;
 static int	afsize=0;
 int		hsize,fsize;
 int		ix0,iy0,grid;
 double		vgain,dgrid;
 double		det,nrm;
 psf		*p;

 if ( img==NULL || img->data==NULL )	return(-1);

 if ( sgp->is_intinelect )	dflux/=sgp->gain;
 if ( sgp->is_photnoise )	vgain=sgp->gain;
 else				vgain=0.0;
 if ( dflux*sgp->gain<=0.0 )	return(0);

 det=is*is+il*il-id*id-ik*ik;
 if ( is<=0.0 && det<=0.0 )	nrm=1.0;
 else				nrm=sqrt(is*is+il*il+id*id+ik*ik);

 ix0=(int)floor(x0),x0-=(double)ix0,
 iy0=(int)floor(y0),y0-=(double)iy0;

 if ( sgp->subpixeldata==NULL || sgp->subg<=0 )	grid=1;
 else						grid=sgp->subg;

 p=sgp->tpd;
 if ( p==NULL )	return(1);

 hsize=(int)(nrm*(double)(p->hsize+2));
 fsize=grid*(2*hsize+1);
 dgrid=(double)grid;

 if ( fsize>afsize )
  {	if ( iarr  != NULL )	tensor_free(iarr);
	iarr =(double **)tensor_alloc_2d(double,fsize,fsize);
	afsize=fsize;
  }

 star_draw_psf(iarr,fsize,fsize,(x0+hsize)*dgrid,(y0+hsize)*dgrid,p,px,py,
	is*dgrid,id*dgrid,ik*dgrid,il*dgrid);

 draw_star_from_array(img,ix0,iy0,iarr,hsize,grid,dflux,sgp->subpixeldata,vgain);
 
 return(0); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int draw_starlist(fitsimage *img,stargenparam *sgp,star *stars,int nstar,int zoom)
{
 int	o,j,l,n,order;
 double	s,d,k,dflux,zf,iz,gs,mom[MAX_DEVIATION_COEFF],*wmom,wz;
 star	*ws;

 if ( zoom<1 )	zf=1.0,zoom=1;
 else		zf=(double)zoom;
 iz=1.0/(double)zoom;

 for ( n=0 ; n<nstar ; n++ )
  {	ws=&stars[n];

	dflux=ws->flux;
	if ( dflux<=0.0 )	continue;

	if ( ws->shape.model==SHAPE_GAUSS || ws->shape.model==SHAPE_ELLIPTIC )
	 {	s=ws->gsig;
		if ( ws->shape.model==SHAPE_ELLIPTIC )
			d=ws->gdel,k=ws->gkap;
		else	
			d=k=0;

		if ( ! sgp->method )
		 {	if ( sgp->is_photnoise )
				dflux=get_gaussian(dflux,sqrt(dflux/sgp->gain));
			draw_star_montecarlo_gauss(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,dflux,s*zf,d*zf,k*zf);
		 }
		else
			draw_star_integral_gauss(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,dflux,s*zf,d*zf,k*zf);
	 }
	else if ( ws->shape.model==SHAPE_DEVIATED )
	 {	if ( zoom>1 )
	 	 {	wz=iz*iz;
			gs=ws->shape.gs*wz;
			order=ws->shape.order;
			for ( o=2,l=0 ; o<=order ; o++ )
			 {	for ( j=0 ; j<=o ; j++,l++ )
				 {	mom[l]=ws->shape.mom[l]*wz;	}
				wz=wz*iz;
			 }
			wmom=mom;
		 }
		else
		 {	gs=ws->shape.gs;
			order=ws->shape.order;
			wmom=ws->shape.mom;
		 }
		draw_star_integral_deviated(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,dflux,gs,order,wmom);
	 }
	else if ( ws->shape.model==SHAPE_PSF )
	 {	draw_star_integral_psf(img,sgp,ws->location.gcx*zf,ws->location.gcy*zf,
		ws->location.gcx,ws->location.gcy,
		ws->shape.gs*zf,ws->shape.gd*zf,ws->shape.gk*zf,ws->shape.gl*zf,dflux);
	 }
  }
 return(0);		
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int quantize_image(fitsimage *img)
{
 int	i,j,sx,sy;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);
 
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[i][j]=floor(img->data[i][j]);		}
  }

 return(0);
}

int divide_image(fitsimage *img,double d)
{
 int	i,j,sx,sy;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);
 if ( d<=0.0 )				return(1);
 d=1.0/d;
 
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[i][j]=d*img->data[i][j];		}
  }

 return(0);
}

/*****************************************************************************/

// read_subpixel_file removed (CLI file I/O)
// normalize_subpixeldata removed (CLI-only)

// read_input_list removed (CLI-only)

// write_output_list removed (CLI-only)

/*****************************************************************************/

#define		MAX_COL	32

// read_input_list removed (CLI-only)

// write_output_list removed (CLI-only)

// fprint_firandom_usage removed (CLI-only)

// get_seed_arg removed (CLI-only)

// get_random_seed removed (CLI-only)

// set_seed removed (CLI-only)

/*****************************************************************************/

