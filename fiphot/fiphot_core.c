/*****************************************************************************/
/* fiphot_core.c — algorithm functions from fiphot.c + flat wrappers.   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool for performing photometry on FITS images.		     */
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>

#include "image.h"


#include "fitsh.h"

#include "mask.h"


#include "tokenize.h"
#include "math/spline/biquad.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "math/point.h"
#include "math/polygon.h"
#include "statistics.h"
#include "magnitude.h"

#include "tensor.h"
#include "common.h"
#include "apphot.h"
#include "kernel.h"
#include "weight.h"

#include "fiphot.h"
#include "fiphot_core.h"
#include "dump_helper.h"

#ifdef	HAVE_NO_CC_EXTENSION
#define	__extension__
#endif

/*****************************************************************************/

int	is_verbose,is_comment;
char	*progbasename;

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

/*****************************************************************************/

int get_kernel_block_coords(int id,int *rbx,int *rby)
{
 int	bx,by,fs,bl,hs;

 if ( id<=0 )	bx=by=0;
 else
  {	hs=1,fs=3;bl=8;id--;
	while ( 1 )
	 {	if ( id<bl )
		 {	if ( id<fs )		bx=-hs+id,by=-hs;
			else if ( id>=bl-fs )	bx=-hs+id-bl+fs,by=hs;
			else
			 {	id-=fs;
				by=-hs+1+id/2;
				if ( id%2==0 )	bx=-hs;
				else		bx=+hs;
			 }
			break;
		 }
		else
		 {	id-=bl,hs++,
			fs+=2,bl+=8;
		 }
	 };
  }
 *rbx=bx,
 *rby=by;
 return(0);
}
int get_kernel_block_id(int bx,int by)
{
 int	hs,b;
 if ( bx==0 && by==0 )	return(0);
 hs=0;
 if ( +bx>hs )	hs=+bx;
 if ( -bx>hs )	hs=-bx;
 if ( +by>hs )	hs=+by;
 if ( -by>hs )	hs=-by;
 b=(2*hs-1)*(2*hs-1);
      if ( by==-hs )	return(b+hs+bx);
 else if ( by==+hs )	return(b+7*hs-1+bx);
 else if ( bx==-hs )	return(b+2*hs+1+2*(hs+by-1)+0);
 else if ( bx==+hs )	return(b+2*hs+1+2*(hs+by-1)+1);
 else			return(-1);
}

/*****************************************************************************/

double optimal_aperture_xi_approx(double m)
{
 double	x,xi,b,d,c,q,w;
 x=log(m);
 b=2.316,d=-1.914,q=3.0;
 w=fabs(x-d);
 c=2.512862417;
 xi=sqrt(pow(b*b+pow(w,q),1.0/q)-(x-d)+c);
 return(xi);
}
double optimal_aperture_xi(double m,int n)
{
 double	xi,k,t;
 int	i;
 xi=optimal_aperture_xi_approx(m);
 for ( i=0 ; i<n ; i++ )
  {	t=2.0*m*(1.0+xi*xi)+1.0;
	k=0.5*(t-sqrt(t*t-8.0*m));
	xi=sqrt(-2.0*log(k));
  }
 return(xi);
}

double optimal_aperture(double g,double s,double bg,double bgar,double flux)
{
 double	m,xi,r;
 m=g*(s*s*bg*bg)*(1.0+1.0/bgar)*M_PI/flux;
 xi=optimal_aperture_xi(m,1);
 r=s*xi;
 return(r);
}

/*****************************************************************************/

int read_subpixel_file(FILE *fr,double ***rsubpixeldata,int *rsubg)
{
 double **subpixeldata,x,y,v;
 int	subg,n,i,j,ix,iy;
 char	buff[256],*cmd[4];
 point	*points;
 int	npoint;
 points=NULL;npoint=0;
 while ( ! feof(fr) )
  {	if ( fgets(buff,255,fr)==NULL )	break;
	remove_newlines_and_comments(buff);
	n=tokenize_spaces(buff,cmd,3);
	if ( n<3 )	continue;
	ix=iy=-1;
	sscanf(cmd[0],"%lg",&x);ix=(int)floor(x);
	sscanf(cmd[1],"%lg",&y);iy=(int)floor(y);
	if ( ! isfinite(x) || ! isfinite(y) )	continue;
	if ( ix<0 || iy<0 )			continue;
	sscanf(cmd[2],"%lg",&v);
	if ( ! isfinite(v) )			continue;
	points=(point *)realloc(points,sizeof(point)*(npoint+1));
	points[npoint].x=(double)ix,
	points[npoint].y=(double)iy;
	points[npoint].value=v;
	npoint++;
  };
 subg=0;
 for ( i=0 ; i<npoint ; i++ )
  {	ix=(int)floor(points[i].x),
	iy=(int)floor(points[i].y);
	if ( ix>subg )	subg=ix;
	if ( iy>subg )	subg=iy;
  }
 subg++;
 subpixeldata=tensor_alloc_2d(double,subg,subg);
 for ( i=0 ; i<subg ; i++ )
  { for ( j=0 ; j<subg ; j++ )
     {	subpixeldata[i][j]=1.0;		}
  }
 for ( i=0 ; i<npoint ; i++ )
  {	ix=(int)floor(points[i].x),
	iy=(int)floor(points[i].y);
	subpixeldata[iy][ix]=points[i].value;
  }
 free(points);

 *rsubpixeldata=subpixeldata;
 *rsubg=subg;
	
 return(0);
}
int normalize_subpixeldata(double **d,int g)
{
 int	i,j;
 double	s,n;
 for ( i=0,s=0.0,n=0.0 ; i<g ; i++ )
  {	for ( j=0 ; j<g ; j++,n=n+1.0 )
	 {	s+=d[i][j];		}
  }
 if ( s<=0.0 )	return(1);
 for ( i=0 ; i<g ; i++ )
  {	for ( j=0 ; j<g ; j++ )
	 {	d[i][j]=d[i][j]*n/s;	}
  }

 return(0);
}


/*****************************************************************************/

int calculate_magnitudes(photstar *ps,int np,magflux *mf0)
{
 int		i,j;
 photflux	*pf;
 magflux	*mf,mf1;

 for ( i=0 ; i<np ; i++ )
  {	for ( j=0 ; j<ps[i].n ; j++ )
	 {	pf=&ps[i].fluxes[j];
		if ( ps[i].use_ref && ps[i].rfflux != NULL )
		 {	mf1.magnitude=ps[i].ref_mag;
			mf1.intensity=ps[i].rfflux[j].flux;
			mf=&mf1;
		 }
		else
			mf=mf0;
		flux_to_mag_magerr(pf->flux,pf->fluxerr,mf,&pf->mag,&pf->magerr);
	 }
  }

 return(0); 
}

/*****************************************************************************/

int add_to_data_weight(double **data,char **mask,int sx,int sy,
	int hsize,int grid,weight *wg,double mul)
{
 int	i,j,fsize,ix,iy;

 if ( data==NULL || wg==NULL )	return(1);

 fsize=grid*(2*hsize+1);

 for ( i=0 ; i<fsize ; i++ )
  {	iy=(i/grid)-hsize+wg->iy;
	if ( iy<0 || iy>=sy )	continue;
	for ( j=0 ; j<fsize ; j++ )
	 {	ix=(j/grid)-hsize+wg->ix;
		if ( ix<0 || ix>=sx )			continue;
		if ( mask != NULL && mask[iy][ix] )	continue;
		data[iy][ix] += mul * wg->iarr[i][j];
	 }
  }

 return(0);
}
int add_to_image_weights(image *img,char **mask,weightlist *wl,double mul)
{
 int	i;

 if ( img==NULL || img->data==NULL )	return(1);
 if ( img->sx<=0 || img->sy<=0 )	return(1);
 if ( wl==NULL || wl->weights==NULL )	return(0);

 for ( i=0 ; i<wl->nweight ; i++ )
  {	add_to_data_weight(img->data,mask,img->sx,img->sy,
		wl->hsize,wl->grid,&wl->weights[i],mul);
  }

 return(0);
}

/*****************************************************************************/

int ringmask_subtract(char **ringmask,int sx,int sy,double x,double y,double r)
{
 int	i,j,i0,i1,j0,j1;
 double	r2,dx,dy;

 r2=r*r;
 i0=(int)(y-r-1.0); if ( i0<0  )	i0=0;
 i1=(int)(y+r+1.0); if ( i1>sy )	i1=sy;
 j0=(int)(x-r-1.0); if ( j0<0  )	j0=0;
 j1=(int)(x+r+1.0); if ( j1>sx )	j1=sx;
 for ( i=i0 ; i<i1 ; i++ )
  {	dy=(double)(i+0.5)-y;
	for ( j=j0 ; j<j1 ; j++ )
	 {	dx=(double)(j+0.5)-x;
		if ( dx*dx+dy*dy<=r2 && ringmask[i][j] )
			ringmask[i][j]--;
	 }
  }
 return(0);
}

int ringmask_coadd(char **ringmask,int sx,int sy,double x,double y,double r)
{
 int	i,j,i0,i1,j0,j1;
 double	r2,dx,dy;

 r2=r*r;
 i0=(int)(y-r-1.0); if ( i0<0  )	i0=0;
 i1=(int)(y+r+1.0); if ( i1>sy )	i1=sy;
 j0=(int)(x-r-1.0); if ( j0<0  )	j0=0;
 j1=(int)(x+r+1.0); if ( j1>sx )	j1=sx;
 for ( i=i0 ; i<i1 ; i++ )
  {	dy=(double)(i+0.5)-y;
	for ( j=j0 ; j<j1 ; j++ )
	 {	dx=(double)(j+0.5)-x;
		if ( dx*dx+dy*dy<=r2 && ringmask[i][j]<2 )
			ringmask[i][j]++;
	 }
  }
 return(0);
}

/* 
 dradius:
	- positive: use this value as disjoint area radius
	- zero: use the aperture radius as disjoint area radius
	- negative: use the inner radius of the annulus as disjoint area radius
*/

char **	ringmask_create(int sx,int sy,photstar *ps,int np,apgeom *inaps,int a,double dradius)
{
 char **ringmask;
 int	i,n;
 apgeom	*aps;
 double	r;

 ringmask=(char **)tensor_alloc_2d(char,sx,sy);
 for ( i=0 ; i<sy ; i++ )
  {	memset(ringmask[i],0,sx);		}

 for ( n=0 ; n<np ; n++ )
  {	
	if ( dradius>0.0 )
		r=dradius;
	else 
	 {	if ( inaps != NULL )	aps=inaps;
		else			aps=ps[n].inaps;
		if ( aps==NULL )	continue; /* however, unexpected. */
		if ( dradius<0.0 )	r=aps[a].ra;
		else			r=aps[a].r0;
	 }

 	ringmask_coadd(ringmask,sx,sy,ps[n].x,ps[n].y,r);
  }

 return(ringmask);
}

int do_photometry(image *img,char **mask,
	photstar *ps,int np,weightlist *wl,int weightusage,int is_calc_opt_apert,
	apgeom *inaps,int numap,
	spatialgain *sg,double correlation_length,double sigma,xphotpar *xpp)
{
 int		i,j,k,jp,sx,sy,r;
 apphotpar	ap;
 apgeom		*aps;
 photflux	*pf,*wpf;

 double		bgarea,bgflux,bgmedian,bgsigma,gain;
 double		area,flux,fluxerr,cx,cy,dx,dy,cfd2;
 double		**bqc;
 weight		*ww;

 double		**aphdata,**aphwsum,grid2,igrid2;
 char		**aphmask;
 char		***ringmasks;
 int		hsize,fsize,grid,wfsize;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);

 if ( xpp->use_biquad )
  {	bqc=tensor_alloc_2d(double,sx*2+1,sy*2+1);
	if ( bqc==NULL )	return(-1);
	biquad_coeff(img->data,sx,sy,bqc,NULL);
	weightusage &= ~(USE_WEIGHT_SUBTRACTED|USE_WEIGHT_WEIGHTED);
  }
 else
	bqc=NULL;

 if ( wl==NULL || ! weightusage )	wl=NULL,weightusage=0;

 memcpy(&ap.bgm,&xpp->bgm,sizeof(bgmode));

 if ( wl != NULL )
  {	grid=wl->grid;
	wfsize=grid*(2*wl->hsize+1);
  }
 else
  {	grid=0;
	wfsize=0;
  }

 if ( weightusage & USE_WEIGHT_SUBTRACTED )
  { 	add_to_image_weights(img,mask,wl,-1.0);	
	cfd2=xpp->wconfdist*xpp->wconfdist;
  }
 else	
	cfd2=0.0;

 if ( weightusage & USE_WEIGHT_WEIGHTED )
  {	hsize=0;
	if ( inaps != NULL )
	 {	aps=inaps;
		for ( j=0 ; j<numap ; j++ )
		 {	k=(int)(aps[j].r0+1.0);
			if ( k>hsize )	hsize=k;
		 }
	 }
	else
	 {	for ( i=0 ; i<np ; i++ )
		 {	aps=ps[i].inaps;
			for ( j=0 ; j<numap && aps != NULL ; j++ )
			 {	k=(int)(aps[j].r0+1.0);
				if ( k>hsize )	hsize=k;
			 }
		 }
	 }

	fsize=grid*(2*hsize+1);
	aphdata=(double **)tensor_alloc_2d(double,fsize,fsize);
	aphmask=(char **)tensor_alloc_2d(char,fsize,fsize);
	aphwsum=(double **)tensor_alloc_2d(double,2*wl->hsize+1,2*wl->hsize+1);
  }
 else
  {	hsize=fsize=0;
	grid=1;
	aphdata=NULL;
	aphwsum=NULL;
	aphmask=NULL;
  }

 if ( xpp->is_disjoint_rings )
  {	ringmasks=(char ***)malloc(sizeof(char **)*numap);
	for ( j=0 ; j<numap ; j++ )
	 {	ringmasks[j]=ringmask_create(sx,sy,ps,np,inaps,j,-1.0);		}
  }
 else if ( xpp->is_disjoint_apertures )
  {	ringmasks=(char ***)malloc(sizeof(char **)*numap);
	for ( j=0 ; j<numap ; j++ )
	 {	ringmasks[j]=ringmask_create(sx,sy,ps,np,inaps,j,0.0);		}
  }
 else if ( xpp->disjoint_radius > 0.0 )
  {	ringmasks=(char ***)malloc(sizeof(char **)*numap);
	for ( j=0 ; j<numap ; j++ )
	 {	ringmasks[j]=ringmask_create(sx,sy,ps,np,inaps,j,xpp->disjoint_radius);		}
  }
 else
	ringmasks=NULL;

 grid2=(double)(grid*grid);
 igrid2=1.0/grid2;

 for ( i=0 ; i<np ; i++ )
  {	
	if ( inaps != NULL )	aps=inaps;
	else			aps=ps[i].inaps;

	if ( aps==NULL )	continue;	/* however, unexpected. */

	pf=(photflux *)malloc(sizeof(photflux)*numap);
	ps[i].fluxes=pf;
	ps[i].n=numap;

	cx=ps[i].x,
	cy=ps[i].y;

	gain=eval_2d_poly(cx,cy,sg->order,sg->coeff,0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);
	if ( 0<sg->vmin && gain<sg->vmin )	gain=sg->vmin;
	ap.gain=gain;

	if ( weightusage & USE_WEIGHT_SUBTRACTED )
	 {	ww=weight_get_closest(wl,cx,cy);
		dx=cx-ww->x,dy=cy-ww->y;
		if ( cfd2<=0.0 || dx*dx+dy*dy<cfd2 )
			add_to_data_weight(img->data,mask,sx,sy,wl->hsize,wl->grid,ww,+1.0);
		else	
			ww=NULL;
	 }
	else
		ww=NULL;

/*
	if ( weightusage & USE_WEIGHT_WEIGHTED )
	 {	int	ix0,iy0,i,j;
		int	ix,iy,wx,wy;
		weight	*ww;

		ww=weight_get_closest(wl,cx,cy);

		ix0=((int)floor(cx))-hsize;
		iy0=((int)floor(cy))-hsize;

		for ( i=0 ; i<wfsize/grid ; i++ )
		 {	for ( j=0 ; j<wfsize/grid ; j++ )
			 {	aphwsum[i][j]=0.0;		}
		 }
		for ( i=0 ; i<wfsize ; i++ )
		 {	for ( j=0 ; j<wfsize ; j++ )
			 {	aphwsum[i/grid][j/grid]+=ww->iarr[i][j];	}
		 }
		
		for ( i=0 ; i<fsize ; i++ )
		 {	iy=iy0+i/grid;
			for ( j=0 ; j<fsize ; j++ )
			 {	ix=ix0+j/grid;
				if ( ix<0 || iy<0 || ix>=sx || iy>=sy )
				 {	aphmask[i][j] = MASK_OUTER;
					aphdata[i][j] = 0.0;
					continue;
				 }
				else if ( mask != NULL )
				 {	aphmask[i][j] = mask[iy][ix];	}
				else
				 {	aphmask[i][j] = 0;		}

				aphdata[i][j]=img->data[iy][ix]*igrid2;

				wx=j+grid*(-hsize+wl->hsize);
				wy=i+grid*(-hsize+wl->hsize);
				if ( wx>=0 && wy>=0 && wx<wfsize && wy<wfsize )
				 {	aphdata[i][j]+=ww->iarr[wy][wx];
					aphdata[i][j]-=aphwsum[wy/grid][wx/grid]*igrid2;
				 }
			
			 }
		 }
	 }
*/

	for ( j=0 ; j<numap ; j++ )
	 {	int		rtot,rbad,rign,atot,abad;
		int		apgeom_type;
		apphot_out	out;
		double		background;

		apgeom_type=aps[j].apgeom_type;

		ap.r0=aps[j].r0,
		ap.ra=aps[j].ra,
		ap.da=aps[j].da;

		ap.r0_poly=aps[j].r0_poly,ap.nr0=aps[j].nr0;
		ap.ra_poly=aps[j].ra_poly,ap.nra=aps[j].nra;
		ap.da_poly=aps[j].da_poly,ap.nda=aps[j].nda;

		wpf=&pf[j];		
		wpf->ag.r0=ap.r0;
		wpf->ag.ra=ap.ra;
		wpf->ag.da=ap.da;

		for ( k=0,jp=-1 ; k<j && jp<0 ; k++ )
		 {	if ( pf[k].ag.ra==ap.ra && pf[k].ag.da==ap.da )
			 {	jp=k;break;		}
		 }

		if ( jp>=0 )
		 {	bgarea=pf[jp].bgarea,
			bgflux=pf[jp].bgflux,
			bgmedian=pf[jp].bgmedian,
			bgsigma=pf[jp].bgsigma;
			atot=pf[jp].atot,
			abad=pf[jp].abad;
			r=0;
		 }
		else if ( (! xpp->use_sky) && ringmasks != NULL )
		 {	r=aperture_photometry_back_ring(img->data,mask,sx,sy,cx,cy,
				&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,
				&atot,&abad,ringmasks[j]);
		 }
		else if ( (! xpp->use_sky) )
		 {	if ( apgeom_type==APGEOM_TYPE_CIRCULAR )
			 {	r=aperture_photometry_back_ring(img->data,mask,sx,sy,cx,cy,
					&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,
					&atot,&abad,NULL);
			 }
			else if ( apgeom_type==APGEOM_TYPE_POLYGON )
			 {	r=aperture_photometry_back_polygons(img->data,mask,sx,sy,cx,cy,
					&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,
					&atot,&abad,NULL);
			 }
			else
				r=1;
		 }
		else
		 {	if ( apgeom_type==APGEOM_TYPE_CIRCULAR )
				bgarea=M_PI*(2.0*ap.ra*ap.da+ap.da*ap.da);
			else
				bgarea=polygon_area(ap.da_poly,ap.nda)-polygon_area(ap.ra_poly,ap.nra);

			bgflux=xpp->sky*bgarea;
			bgmedian=xpp->sky;
			bgsigma=0.0;
			atot=(int)bgarea;
			abad=0;
			r=0;
		 }

		background=bgmedian;

		if ( r )
		 {	wpf=&pf[j];
			wpf->flag=MASK_NOBACKGROUND;
			wpf->rtot=0;
			wpf->rbad=0;
			wpf->rign=0;
			wpf->atot=0;
			wpf->abad=0;
			continue;
		 }

		wpf=&pf[j];
		wpf->bgmedian=bgmedian,
		wpf->bgsigma=bgsigma;
		wpf->bgflux=bgflux;
		wpf->bgarea=bgarea;
		wpf->atot=atot;
		wpf->abad=abad;

		/*
		if ( xpp->is_disjoint_apertures )
			ringmask_subtract(ringmasks[j],sx,sy,cx,cy,ap.ra);
		*/

		wpf=&pf[j];
		wpf->flag=0;

		if ( bqc != NULL )
		 {	/*
			if ( xpp->is_disjoint_apertures || xpp->disjoint_radius>0.0 )
				r=aperture_photometry_flux_biquad(bqc,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&rtot,&rbad,&rign,xpp->maskignore,ringmasks[j]);
			else
			*/
			r=aperture_photometry_flux_biquad(bqc,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&rtot,&rbad,&rign,xpp->maskignore,NULL);
		 }
		else if ( ! (weightusage & USE_WEIGHT_WEIGHTED) )
		 {	/*
			if ( xpp->is_disjoint_apertures || xpp->disjoint_radius>0.0 )
				r=aperture_photometry_flux(img->data,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&rtot,&rbad,&rign,xpp->maskignore,xpp->subpixeldata,xpp->subg,ringmasks[j]);
			else
			*/
			if ( apgeom_type==APGEOM_TYPE_CIRCULAR )
				r=aperture_photometry_flux_circle(img->data,mask,sx,sy,cx,cy,ap.r0,&area,&flux,&out,background,&rtot,&rbad,&rign,xpp->maskignore,xpp->subpixeldata,xpp->subg,NULL);
			else if ( apgeom_type==APGEOM_TYPE_POLYGON )
				r=aperture_photometry_flux_polygon(img->data,mask,sx,sy,cx,cy,ap.r0_poly,ap.nr0,&area,&flux,&out,background,&rtot,&rbad,&rign,xpp->maskignore,NULL);
			else
			 {	r=-1;
				rtot=rbad=rign=0;
			 }
		 }
		else
		 {	double	r0;
			double	ccx,ccy,cx0,cy0,bg;
			int	ix0,iy0,i,j;
			int	ix,iy,wx,wy;
			weight	*ww;

			cx0=floor(cx);
			cy0=floor(cy);
			ccx=(double)grid*((cx-cx0)+(double)hsize);
			ccy=(double)grid*((cy-cy0)+(double)hsize);
			r0 =(double)grid*ap.r0;

			ww=weight_get_closest(wl,cx,cy);

			ix0=(int)cx0-hsize;
			iy0=(int)cy0-hsize;

			for ( i=0 ; i<wfsize/grid ; i++ )
			 {	for ( j=0 ; j<wfsize/grid ; j++ )
				 {	aphwsum[i][j]=0.0;		}
			 }
			for ( i=0 ; i<wfsize ; i++ )
			 {	for ( j=0 ; j<wfsize ; j++ )
				 {	aphwsum[i/grid][j/grid]+=ww->iarr[i][j];	}
			 }

			bg=bgmedian*igrid2;
		
			for ( i=0 ; i<fsize ; i++ )
			 {	iy=iy0+i/grid;
				for ( j=0 ; j<fsize ; j++ )
				 {	ix=ix0+j/grid;
					if ( ix<0 || iy<0 || ix>=sx || iy>=sy )
					 {	aphmask[i][j] = MASK_OUTER;
						aphdata[i][j] = 0.0;
						continue;
					 }
					else if ( mask != NULL )
					 {	aphmask[i][j] = mask[iy][ix];	}
					else
					 {	aphmask[i][j] = 0;		}

					aphdata[i][j]=img->data[iy][ix]*igrid2-bg;

					wx=j+grid*(-hsize+wl->hsize);
					wy=i+grid*(-hsize+wl->hsize);
					if ( wx>=0 && wy>=0 && wx<wfsize && wy<wfsize )
					 {	aphdata[i][j] *= grid2*ww->iarr[wy][wx]/aphwsum[wy/grid][wx/grid];	}

					aphdata[i][j]+=bg;
			
				 }
			 }

			r=aperture_photometry_flux_circle(aphdata,aphmask,fsize,fsize,ccx,ccy,r0,&area,&flux,&out,background,&rtot,&rbad,&rign,xpp->maskignore,xpp->subpixeldata,xpp->subg,NULL);
		 }
		wpf->rtot=rtot;
		wpf->rbad=rbad;
		wpf->rign=rign;
		wpf->flag |= r;	

		flux -= bgmedian*area*igrid2;
		wpf->flux=flux;
		if ( flux<=0.0 )
		 { 	wpf->fluxerr=0.0;
		 }
		else
		 {	fluxerr=sqrt((0.0<gain?flux/gain:0.0)+area*(bgsigma*bgsigma*correlation_length*correlation_length)*(1.0+1.0/bgarea)),
			wpf->fluxerr=fluxerr;
		 }

		if ( 0.0 < out.fw && 0 < flux )
		 {	double	mx,my,mxx,mxy,myy,sxx,syy,sxy,ac0,w;
			mx=out.fwx/out.fw;
			my=out.fwy/out.fw;
			wpf->cntr_x=mx+cx;
			wpf->cntr_y=my+cy;
			mxx=out.fwxx/out.fw;
			mxy=out.fwxy/out.fw;
			myy=out.fwyy/out.fw;
			sxx=mxx-mx*mx;
			sxy=mxy-mx*my;
			syy=myy-my*my;
			ac0=area*bgsigma/flux;
			wpf->cntr_x_err=sqrt((0.0<gain?(sxx)/(flux*gain):0.0)+(ac0*ac0)/(4.0*M_PI));
			wpf->cntr_y_err=sqrt((0.0<gain?(syy)/(flux*gain):0.0)+(ac0*ac0)/(4.0*M_PI));
			w=sxx*syy-sxy*sxy;
			if ( w <= 0.0 )	
			 {	wpf->cntr_width=0.0;
				wpf->cntr_w_err=-1.0;
			 }
			else	
			 {	double	sig2;
				double	a,b,c,ws,wd,wk;
				a=(sxx+syy)/2;
				b=(sxx-syy)/2;
				c=sxy;
				ws=sqrt(0.5*(a+sqrt(a*a-b*b-c*c)));
				wd=b/(2*ws);
				wk=c/(2*ws);
				wpf->cntr_width=ws;
				wpf->cntr_w_d=wd;
				wpf->cntr_w_k=wk;
				sig2=sqrt(w);
				wpf->cntr_w_err=sqrt((0.0<gain?(2*sig2*sig2/(flux*gain)):0.0)+ac0*ac0*area/(8.0*M_PI));
			 }
		 }
		else
		 {	wpf->cntr_x=0.0;
			wpf->cntr_y=0.0;
			wpf->cntr_x_err=-1.0;
			wpf->cntr_y_err=-1.0;
			wpf->cntr_width=-1.0;
			wpf->cntr_w_err=-1.0;
		 }


		/*
		if ( xpp->is_disjoint_apertures )
			ringmask_coadd(ringmasks[j],sx,sy,cx,cy,ap.ra);
		*/

	 } /* for ( j=0 ; j<numap ; j++ ) */

	if ( ww != NULL )
		add_to_data_weight(img->data,mask,sx,sy,wl->hsize,wl->grid,ww,-1.0);

	if ( is_calc_opt_apert && numap==1 && ( ! pf[0].flag ) )
	 {	double	r0opt;

		if ( 0.0<gain )
			r0opt=optimal_aperture(gain,sigma,pf[0].bgsigma,pf[0].bgarea,pf[0].flux);
		else
			r0opt=0.0;

		ps[i].optimal.r0=r0opt;
		ps[i].optimal.ra=pf[0].ag.ra;
		ps[i].optimal.da=pf[0].ag.da;
	 }
	else
	 {	ps[i].optimal.r0=0.0;
		ps[i].optimal.ra=0.0;
		ps[i].optimal.da=0.0;
	 }
  }

 if ( ringmasks != NULL )
  {	for ( j=numap-1 ; j>=0 ; j-- )
	 {	if ( ringmasks[j] != NULL )
			tensor_free(ringmasks[j]);	
	 }
	free(ringmasks);
	ringmasks=NULL;
  }

 if ( aphmask != NULL )	tensor_free(aphmask);
 if ( aphdata != NULL )	tensor_free(aphdata);
 if ( bqc != NULL )	tensor_free(bqc);

 if ( weightusage & USE_WEIGHT_SUBTRACTED )
 	add_to_image_weights(img,mask,wl,+1.0);

 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int do_subtracted_photometry(image *img,char **mask,photstar *ps,int np,
	apgeom *inaps,int numap,
	spatialgain *sg,double correlation_length,xphotpar *xpp,kernellist *klist,int normalize_kernel)
{
 int		i,j,k,jp,bx,by,sx,sy,r,atot,abad,rtot,rbad,rign,flag;
 apphotpar	ap;
 apgeom		*aps;
 int		n;
 photflux	*pf,*wpf;
 double		bgarea,bgflux,bgmedian,bgsigma,gain;
 double		area,flux,fluxerr,cx,cy,kx,ky,flux0,sflux,kflux;
 double		**bqc;
 photflux	**fluxarr;
 double		**kernarr;
 int		hsize,fsize;

 if ( img==NULL || img->data==NULL )	return(1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(1);

 if ( xpp->use_biquad )
  {	bqc=tensor_alloc_2d(double,sx*2+1,sy*2+1);
	if ( bqc==NULL )	return(-1);
	biquad_coeff(img->data,sx,sy,bqc,NULL);
  }
 else
	bqc=NULL;

 memcpy(&ap.bgm,&xpp->bgm,sizeof(bgmode));

 if ( klist != NULL )
  {	hsize=0;
	for ( k=0 ; k<klist->nkernel ; k++ )
	 {	if ( klist->kernels[k].hsize>hsize )
			hsize=klist->kernels[k].hsize;
	 }
  }
 else
	hsize=0;

 fsize=2*hsize+1;
 fluxarr=(photflux **)tensor_alloc_2d(photflux,fsize,fsize);
 kernarr=(double   **)tensor_alloc_2d(double  ,fsize,fsize);

 for ( i=0 ; i<np ; i++ )
  {	if ( inaps != NULL )	aps=inaps,n=numap;
	else			aps=ps[i].inaps,n=ps[i].n;

	ps[i].n=n;
	pf=(photflux *)malloc(sizeof(photflux)*n);
	ps[i].fluxes=pf;

	cx=ps[i].x ,cy=ps[i].y;
	kx=ps[i].x ,ky=ps[i].y;

	gain=eval_2d_poly(cx,cy,sg->order,sg->coeff,0.5*(double)sx,0.5*(double)sy,0.5*(double)sx);
	if ( sg->vmin>0 && gain<sg->vmin )	gain=sg->vmin;
	ap.gain=gain;

	for ( j=0 ; j<n ; j++ )
	 {	
		int	apgeom_type;

		wpf=&pf[j];

		wpf->flag=0;
		wpf->atot=0;
		wpf->abad=0;
		wpf->rtot=0;	
		wpf->rbad=0;
		wpf->rign=0;

		apgeom_type=aps[j].apgeom_type;

                ap.r0=aps[j].r0,
                ap.ra=aps[j].ra,
                ap.da=aps[j].da;

                ap.r0_poly=aps[j].r0_poly,ap.nr0=aps[j].nr0;
                ap.ra_poly=aps[j].ra_poly,ap.nra=aps[j].nra;
                ap.da_poly=aps[j].da_poly,ap.nda=aps[j].nda;

		for ( k=0,jp=-1 ; k<j && jp<0 ; k++ )
		 {	if ( pf[k].ag.ra==ap.ra && pf[k].ag.da==ap.da )
			 {	jp=k;break;		}
		 }
		if ( jp>=0 )
		 {	bgarea  =pf[jp].bgarea,
			bgflux  =pf[jp].bgflux,
			bgmedian=pf[jp].bgmedian,
			bgsigma =pf[jp].bgsigma;
			atot    =pf[jp].atot;
			abad    =pf[jp].abad;
			r=0;
		 }
		else if ( apgeom_type==APGEOM_TYPE_CIRCULAR )
			r=aperture_photometry_back_ring(img->data,mask,sx,sy,cx,cy,&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,&atot,&abad,NULL);
		else if ( apgeom_type==APGEOM_TYPE_POLYGON )
			r=aperture_photometry_back_polygons(img->data,mask,sx,sy,cx,cy,&ap,&bgarea,&bgflux,&bgmedian,&bgsigma,&atot,&abad,NULL);
		else
			r=1;

		if ( r )
		 {	wpf->flag=MASK_NOBACKGROUND;
			continue;
		 }

		flag=0;

		rtot=0;	
		rbad=0;
		rign=0;

		for ( by=-hsize ; by<=hsize ; by++ )
		 {	for ( bx=-hsize ; bx<=hsize ; bx++ )
			 {	double	a,f;
				int	rt,rb,ri;
				if ( xpp->use_biquad )
					r=aperture_photometry_flux_biquad(bqc,mask,sx,sy,cx+bx,cy+by,ap.r0,&a,&f,&rt,&rb,&ri,xpp->maskignore,NULL);
				else if ( apgeom_type==APGEOM_TYPE_CIRCULAR )
					r=aperture_photometry_flux_circle(img->data,mask,sx,sy,cx+bx,cy+by,ap.r0,&a,&f,NULL,0.0,&rt,&rb,&ri,xpp->maskignore,NULL,0,NULL);
				else if ( apgeom_type==APGEOM_TYPE_POLYGON )
					r=aperture_photometry_flux_polygon(img->data,mask,sx,sy,cx+bx,cy+by,ap.r0_poly,ap.nr0,&a,&f,NULL,0.0,&rt,&rb,&ri,xpp->maskignore,NULL);
				else
				 {	r=-1;
					a=f=0.0;
					rt=rb=ri=0;
				 }
					
				if ( rt>rtot )	rtot=rt;
				if ( rb>rbad )	rbad=rb;
				if ( ri>rign )	rign=ri;

				fluxarr[hsize+by][hsize+bx].flux=f-bgmedian*a;
				fluxarr[hsize+by][hsize+bx].flag=r;
				fluxarr[hsize+by][hsize+bx].bgarea=a;

				flag |= r;
			 }
		 }

		if ( flag )
		 {	wpf->flag |= flag;
			wpf->flux=0.0;
			wpf->fluxerr=0.0;
			continue;
		 }

		if ( klist == NULL )
		 {	flux=fluxarr[0][0].flux;
			kernarr[0][0]=1.0;
			kflux=1.0;
		 }

		else
		 {	int	c;
			kernel	*k;
			double	kcoeff;

			for ( by=0 ; by<fsize ; by++ )
			 {	for ( bx=0 ; bx<fsize ; bx++ )
				 {	kernarr[by][bx]=0.0;		}
			 }

			for ( c=0,k=klist->kernels ; c<klist->nkernel ; c++,k++ )
			 {	if ( k->type==KERNEL_BACKGROUND )
					continue;
				kcoeff=eval_2d_poly(kx,ky,k->order,k->coeff,klist->ox,klist->oy,klist->scale);

				if ( k->type==KERNEL_IDENTITY )
					kernarr[hsize][hsize]+=kcoeff;
				else
				 {	for ( by=-k->hsize ; by<=+k->hsize ; by++ )
					 {  for ( bx=-k->hsize ; bx<=+k->hsize ; bx++ )
					     {	kernarr[hsize+by][hsize+bx]+=kcoeff*k->image[k->hsize+by][k->hsize+bx];	}
					 }
				 }				
			 }

			kflux=0.0;
			for ( by=0 ; by<fsize ; by++ )
			 {	for ( bx=0 ; bx<fsize ; bx++ )
				 {	kflux += kernarr[by][bx];	}
			 }
			if ( normalize_kernel && 0<kflux )
			 {	for ( by=0 ; by<fsize ; by++ )
				 {	for ( bx=0 ; bx<fsize ; bx++ )
					 {	kernarr[by][bx] /= kflux;	}
				 }
				kflux=1.0;
			 }
		 }

		sflux=0.0;
		area=0.0;
		for ( by=0 ; by<fsize ; by++ )
		 {	for ( bx=0 ; bx<fsize ; bx++ )
			 {	sflux+=kernarr[by][bx]*fluxarr[by][bx].flux;
				area +=kernarr[by][bx]*fluxarr[by][bx].bgarea;
			 }
		 }
		sflux/=(kflux*kflux);
		area /=kflux;

		wpf->atot=atot;
		wpf->abad=abad;
		wpf->rtot=rtot;	
		wpf->rbad=rbad;
		wpf->rign=rign;

		wpf->bgmedian=bgmedian,
		wpf->bgsigma=bgsigma;
		wpf->bgflux=bgflux;
		wpf->bgarea=bgarea;

		wpf->ag.r0=ap.r0;
		wpf->ag.ra=ap.ra;
		wpf->ag.da=ap.da;

		if ( ps[i].rfflux[j].flag )
		 {	wpf->flag |= ps[i].rfflux[j].flag;
			wpf->flux=0.0;
			wpf->fluxerr=0.0;
			continue;
		 }

		flux0=ps[i].rfflux[j].flux;

		flux = flux0 + sflux;

		if ( flux<=0.0 || flux0<=0.0 )
		 {	wpf->flux=0.0,
			wpf->fluxerr=0.0;
		 }
		else
		 {	fluxerr=sqrt((0.0<gain?flux/gain:0.0)+area*(bgsigma*bgsigma*correlation_length*correlation_length)*(1.0+1.0/bgarea));
			wpf->flux=flux,
			wpf->fluxerr=fluxerr;
		 }

		wpf->flag=0;
	 }

	ps[i].n=n;
  }

 if ( bqc != NULL )	tensor_free(bqc);

 tensor_free(kernarr);
 tensor_free(fluxarr);

 return(0);
}

/*****************************************************************************/

int do_magnitude_fit(photstar *stars,int nstar,
magfitparams *mfp,int sx,int sy,magfitstat *mfs)
{
 int		a,nvar,i,j,n,maxorder,k,l,v,iiter;
 double		**amatrix,*bvector,*fvars,*monoms,yvar;
 double		ox,oy,scale,dmag,cmul,w;
 photstar	*s;
 photflux	*fi,*fr;
 double		*dmags,*wghts,*smags;
 
 nvar=0;
 maxorder=0;
 for ( j=0 ; j<=mfp->norder ; j++ )
  {	nvar+=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
	if ( maxorder<mfp->orders[j] )
		maxorder=mfp->orders[j];
  }

 /* n: total (maximal) number of apertures */
 n=0;
 for ( i=0 ; i<nstar ; i++ )
  {	if ( n < stars[i].n )	n=stars[i].n;		}

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);
 monoms =vector_alloc((maxorder+1)*(maxorder+2)/2);

 ox=(double)sx/2.0;
 oy=(double)sy/2.0;
 scale=ox;

 dmags=(double *)malloc(sizeof(double)*nstar);
 smags=(double *)malloc(sizeof(double)*nstar);
 wghts=(double *)malloc(sizeof(double)*nstar);

 if ( mfs != NULL )
  {	mfs->nstar=nstar;
	mfs->naperture=n;
	mfs->ninit=(int *)malloc(sizeof(int)*n);
	mfs->nrejs=(int *)malloc(sizeof(int)*n);
	for ( a=0 ; a<n ; a++ )
	 {	mfs->ninit[a]=0;
		mfs->nrejs[a]=0;
	 }
  }

 for ( a=0 ; a<n ; a++ )
  {
	for ( i=0 ; i<nstar ; i++ )
	 {	wghts[i]=-1.0;
		dmags[i]=0.0;

		s=&stars[i];
		if ( s->n < a || s->fluxes==NULL || s->rfflux==NULL )
			continue;
		fi=&s->fluxes[a];
		fr=&s->rfflux[a];

		if ( fi->flag || fr->flag || fi->flux<=0.0 || fr->flux<=0.0 )
			continue;

		dmag=log(fi->flux/fr->flux);
		dmags[i]=dmag;
		w=sqrt(fr->flux);
		wghts[i]=w;

		mfs->ninit[a]++;
	 }

	for ( iiter=0 ; iiter<=mfp->niter ; iiter++ )

	 {	for ( k=0 ; k<nvar ; k++ )
		 {	for ( l=0 ; l<nvar ; l++ )
			 {	amatrix[k][l]=0.0;		}
			bvector[k]=0.0;
		 }

		for ( i=0 ; i<nstar ; i++ )
		 {	
			if ( wghts[i] <= 0.0 )
				continue;

			s=&stars[i];
			fi=&s->fluxes[a];
			fr=&s->rfflux[a];

			dmag=dmags[i];

			v=0;
			cmul=1.0;

			w=wghts[i];	/* weight */
	
			for ( j=0 ; j<=mfp->norder ; j++ )
			 {	eval_2d_monoms(s->x,s->y,mfp->orders[j],fvars+v,ox,oy,scale);
				l=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
				for ( k=0 ; k<l ; k++ )
				 {	fvars[v] *= cmul;
					v++;
				 }
				cmul=cmul*s->ref_col;
			 }
			yvar=dmag;
	
			for ( k=0 ; k<nvar ; k++ )
			 {	for ( l=0 ; l<nvar ; l++ )
				 {	amatrix[k][l] += w*fvars[k]*fvars[l];	}
				bvector[k] += w*fvars[k]*yvar;
			 }
		 }

		solve_gauss(amatrix,bvector,nvar);

		if ( iiter < mfp->niter )
		 {	double	s0,s1,s2,d;

			for ( i=0 ; i<nstar ; i++ )
			 {	s=&stars[i];
				if ( wghts[i] <= 0.0 )
					continue;
				fi=&s->fluxes[a];
				fr=&s->rfflux[a];
	
				v=0;
				cmul=1.0;
				for ( j=0 ; j<=mfp->norder ; j++ )
				 {	eval_2d_monoms(s->x,s->y,mfp->orders[j],fvars+v,ox,oy,scale);
					l=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
					for ( k=0 ; k<l ; k++ )
					 {	fvars[v] *= cmul;
						v++;
					 }
					cmul=cmul*s->ref_col;
				 }
	
				dmag=0.0;
				for ( k=0 ; k<nvar ; k++ )
				 {	dmag += fvars[k]*bvector[k];		}
				smags[i]=dmag;
			 }

			s0=s1=s2=0.0;			
			for ( i=0 ; i<nstar ; i++ )
			 {	if ( wghts[i] <= 0.0 )
					continue;
				d=dmags[i]-smags[i];		
				w=1.0;
				s0+=w;
				s1+=w*d;
				s2+=w*d*d;
			 }
			s1/=s0;
			s2/=s0;
			s2=s2-s1*s1;
			if ( s2<=0.0 )	s2=0.0;
			else		s2=sqrt(s2);

			for ( i=0 ; i<nstar ; i++ )
		 	 {	if ( wghts[i] <= 0.0 )
					continue;
				d=dmags[i]-smags[i];
				if ( fabs(d) >= s2*mfp->sigma )
				 {	wghts[i]=-1.0;
					mfs->nrejs[a]++;
				 }
			 }
		 }
	 }

	for ( i=0 ; i<nstar ; i++ )
	 {	s=&stars[i];
		if ( s->n < a || s->fluxes==NULL || s->rfflux==NULL )
			continue;
		fi=&s->fluxes[a];
		fr=&s->rfflux[a];

		if ( fi->flux<=0.0 || fr->flux<=0.0 )
			continue;

		v=0;
		cmul=1.0;
		for ( j=0 ; j<=mfp->norder ; j++ )
		 {	eval_2d_monoms(s->x,s->y,mfp->orders[j],fvars+v,ox,oy,scale);
			l=(mfp->orders[j]+1)*(mfp->orders[j]+2)/2;
			for ( k=0 ; k<l ; k++ )
			 {	fvars[v] *= cmul;
				v++;
			 }
			cmul=cmul*s->ref_col;
		 }

		dmag=dmags[i];
		for ( k=0 ; k<nvar ; k++ )
		 {	dmag -= fvars[k]*bvector[k];		}

		fi->flux=fr->flux*exp(dmag);
	
	 }

  }

 free(wghts);
 free(smags);
 free(dmags);

 vector_free(monoms);
 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 return(0);
}

/*****************************************************************************/

/*****************************************************************************/
/* Flat-array wrappers                                                       */
/*****************************************************************************/

static void fiphot_fill_photometry_output(photstar *ps, int nstar, int numap,
    double *out_flux, double *out_fluxerr,
    double *out_bgarea, double *out_bgflux, double *out_bgmedian, double *out_bgsigma,
    double *out_cntr_x, double *out_cntr_y, double *out_cntr_width,
    double *out_cntr_w_d, double *out_cntr_w_k,
    double *out_cntr_x_err, double *out_cntr_y_err, double *out_cntr_w_err,
    int *out_flag, int *out_rtot, int *out_rbad, int *out_rign,
    int *out_atot, int *out_abad,
    double *out_optimal_r0, double *out_optimal_ra, double *out_optimal_da)
{
    int i, j, idx;
    if (!out_flux) return;
    for (i = 0; i < nstar; i++) {
        for (j = 0; j < numap; j++) {
            idx = i * numap + j;
            photflux *pf = &ps[i].fluxes[j];
            if (out_flux) out_flux[idx] = pf->flux;
            if (out_fluxerr) out_fluxerr[idx] = pf->fluxerr;
            if (out_bgarea) out_bgarea[idx] = pf->bgarea;
            if (out_bgflux) out_bgflux[idx] = pf->bgflux;
            if (out_bgmedian) out_bgmedian[idx] = pf->bgmedian;
            if (out_bgsigma) out_bgsigma[idx] = pf->bgsigma;
            if (out_cntr_x) out_cntr_x[idx] = pf->cntr_x;
            if (out_cntr_y) out_cntr_y[idx] = pf->cntr_y;
            if (out_cntr_width) out_cntr_width[idx] = pf->cntr_width;
            if (out_cntr_w_d) out_cntr_w_d[idx] = pf->cntr_w_d;
            if (out_cntr_w_k) out_cntr_w_k[idx] = pf->cntr_w_k;
            if (out_cntr_x_err) out_cntr_x_err[idx] = pf->cntr_x_err;
            if (out_cntr_y_err) out_cntr_y_err[idx] = pf->cntr_y_err;
            if (out_cntr_w_err) out_cntr_w_err[idx] = pf->cntr_w_err;
            if (out_flag) out_flag[idx] = pf->flag;
            if (out_rtot) out_rtot[idx] = pf->rtot;
            if (out_rbad) out_rbad[idx] = pf->rbad;
            if (out_rign) out_rign[idx] = pf->rign;
            if (out_atot) out_atot[idx] = pf->atot;
            if (out_abad) out_abad[idx] = pf->abad;
        }
        if (out_optimal_r0) out_optimal_r0[i] = ps[i].optimal.r0;
        if (out_optimal_ra) out_optimal_ra[i] = ps[i].optimal.ra;
        if (out_optimal_da) out_optimal_da[i] = ps[i].optimal.da;
    }
}

static photstar *fiphot_build_stars(double *x, double *y, int nstar,
    double *r0, double *ra, double *da, int nap, int has_per_star)
{
    int i, j;
    photstar *ps = (photstar *)calloc(nstar, sizeof(photstar));
    for (i = 0; i < nstar; i++) {
        ps[i].x = x[i]; ps[i].y = y[i];
        if (has_per_star && r0) {
            ps[i].inaps = (apgeom *)malloc(nap * sizeof(apgeom));
            for (j = 0; j < nap; j++) {
                memset(&ps[i].inaps[j], 0, sizeof(apgeom));
                ps[i].inaps[j].apgeom_type = APGEOM_TYPE_CIRCULAR;
                ps[i].inaps[j].r0 = r0[i * nap + j];
                ps[i].inaps[j].ra = ra[i * nap + j];
                ps[i].inaps[j].da = da[i * nap + j];
            }
            ps[i].ninap = nap;
        }
    }
    return ps;
}

static int fiphot_build_inaps(apgeom **pinaps, int *pninap, int *pnumap,
    char *aperture_spec, double *star_r0, double *star_ra, double *star_da, int nap, int zoom)
{
    int j;
    if (aperture_spec && aperture_spec[0]) {
        return create_input_ap_data(aperture_spec, pinaps, pninap, zoom) ? 1 : (*pnumap = *pninap, 0);
    } else if (star_r0) {
        *pinaps = (apgeom *)malloc(nap * sizeof(apgeom));
        for (j = 0; j < nap; j++) {
            memset(&(*pinaps)[j], 0, sizeof(apgeom));
            (*pinaps)[j].apgeom_type = APGEOM_TYPE_CIRCULAR;
            (*pinaps)[j].r0 = star_r0[j];
            (*pinaps)[j].ra = star_ra[j];
            (*pinaps)[j].da = star_da[j];
        }
        *pnumap = nap; *pninap = nap;
        return 0;
    }
    return 1;
}

static void fiphot_free_photstar(photstar *ps, int nstar)
{
    int i;
    for (i = 0; i < nstar; i++) {
        if (ps[i].fluxes) free(ps[i].fluxes);
        if (ps[i].rfflux && ps[i].rfflux != ps[i].fluxes) free(ps[i].rfflux);
        if (ps[i].inaps) free(ps[i].inaps);
    }
    free(ps);
}

static void fiphot_free_inaps(apgeom *inaps, int ninap)
{
    if (inaps) free(inaps);
}

/*****************************************************************************/

int fiphot_photometry_cy(
    double *img_data, char *img_mask, int sx, int sy,
    double *star_x, double *star_y, int nstar,
    double *star_r0, double *star_ra, double *star_da, int nap,
    char *aperture_spec, int zoom,
    int bg_type, int bg_scatter, int bg_rejniter,
    double bg_rejlower, double bg_rejupper,
    int mask_ignore, int use_biquad, int use_sky, double sky_level,
    int is_disjoint_rings, int is_disjoint_apertures, double disjoint_radius,
    double **subpixel_data, int subg,
    int sg_order, double *sg_coeff, double sg_vmin,
    double correlation_length, double sigma,
    int is_calc_opt_apert,
    weightlist *wl, int weightusage,
    double *out_flux, double *out_fluxerr,
    double *out_bgarea, double *out_bgflux, double *out_bgmedian, double *out_bgsigma,
    double *out_cntr_x, double *out_cntr_y, double *out_cntr_width,
    double *out_cntr_w_d, double *out_cntr_w_k,
    double *out_cntr_x_err, double *out_cntr_y_err, double *out_cntr_w_err,
    int *out_flag, int *out_rtot, int *out_rbad, int *out_rign,
    int *out_atot, int *out_abad,
    double *out_optimal_r0, double *out_optimal_ra, double *out_optimal_da,
    double *out_raw)
{
    int i, r;
    image img;
    char **mk = NULL;
    photstar *ps;
    xphotpar xpp;
    spatialgain sg;
    apgeom *inaps = NULL;
    int ninap = 0, numap = 0;

    if (!img_data || !star_x || !star_y || nstar <= 0 || sx <= 0 || sy <= 0) return 1;

    img.sx = sx; img.sy = sy;
    img.data = (double **)malloc(sy * sizeof(double *));
    for (i = 0; i < sy; i++) img.data[i] = img_data + i * sx;
    if (img_mask) {
        mk = (char **)malloc(sy * sizeof(char *));
        for (i = 0; i < sy; i++) mk[i] = img_mask + i * sx;
    }

    progbasename = "fitsh_cy.fiphot";
    is_verbose = 0; is_comment = 0;

    memset(&xpp, 0, sizeof(xphotpar));
    xpp.bgm.type = bg_type; xpp.bgm.scatter = bg_scatter;
    xpp.bgm.rejniter = bg_rejniter;
    xpp.bgm.rejlower = bg_rejlower; xpp.bgm.rejupper = bg_rejupper;
    xpp.maskignore = mask_ignore; xpp.use_biquad = use_biquad;
    xpp.wconfdist = 0.0; xpp.sky = sky_level; xpp.use_sky = use_sky;
    xpp.is_disjoint_rings = is_disjoint_rings;
    xpp.is_disjoint_apertures = is_disjoint_apertures;
    xpp.disjoint_radius = disjoint_radius;
    xpp.subpixeldata = subpixel_data; xpp.subg = subg;
    /* normalize and invert subpixel data (matches CLI main()) */
    if (xpp.subpixeldata && xpp.subg > 0) {
        normalize_subpixeldata(xpp.subpixeldata, xpp.subg);
        { int si, sj;
          for (si = 0; si < xpp.subg; si++)
          for (sj = 0; sj < xpp.subg; sj++)
              xpp.subpixeldata[si][sj] = 1.0 / xpp.subpixeldata[si][sj];
        }
    }

    sg.order = sg_order; sg.vmin = sg_vmin;
    if (sg_coeff && sg_order >= 0) {
        int ncoeff = (sg_order + 1) * (sg_order + 2) / 2;
        sg.coeff = (double *)malloc(ncoeff * sizeof(double));
        memcpy(sg.coeff, sg_coeff, ncoeff * sizeof(double));
    } else { sg.order = 0; sg.coeff = (double *)malloc(sizeof(double)); sg.coeff[0] = 1.0; }

    r = fiphot_build_inaps(&inaps, &ninap, &numap, aperture_spec, star_r0, star_ra, star_da, nap, zoom);
    if (r) { r = 1; goto cleanup_sg; }

    ps = fiphot_build_stars(star_x, star_y, nstar, star_r0, star_ra, star_da, nap, star_r0 != NULL);
    r = do_photometry(&img, mk, ps, nstar, wl, weightusage, is_calc_opt_apert, inaps, numap, &sg, correlation_length, sigma, &xpp);

    if (r == 0) {
        fiphot_fill_photometry_output(ps, nstar, numap, out_flux, out_fluxerr, out_bgarea, out_bgflux, out_bgmedian, out_bgsigma, out_cntr_x, out_cntr_y, out_cntr_width, out_cntr_w_d, out_cntr_w_k, out_cntr_x_err, out_cntr_y_err, out_cntr_w_err, out_flag, out_rtot, out_rbad, out_rign, out_atot, out_abad, out_optimal_r0, out_optimal_ra, out_optimal_da);
        if (out_raw) {
            int j;
            for (i = 0; i < nstar; i++) {
                double ar0 = 0.0, ara = 0.0, ada = 0.0;
                for (j = 0; j < numap; j++) {
                    int idx = (i * numap + j) * 12;
                    /* aperture params: per-star or global */
                    if (star_r0 && star_ra && star_da && nap > 0) {
                        ar0 = star_r0[i * nap + j];
                        ara = star_ra[i * nap + j];
                        ada = star_da[i * nap + j];
                    } else if (ninap > 0) {
                        /* from inaps structure */
                        ar0 = inaps[j].r0; ara = inaps[j].ra; ada = inaps[j].da;
                    }
                    out_raw[idx + 0] = (double)i;
                    out_raw[idx + 1] = ps[i].x;
                    out_raw[idx + 2] = ps[i].y;
                    out_raw[idx + 3] = (double)j;
                    out_raw[idx + 4] = (double)ps[i].fluxes[j].flag;
                    out_raw[idx + 5] = ps[i].fluxes[j].bgmedian;
                    out_raw[idx + 6] = ps[i].fluxes[j].bgsigma;
                    out_raw[idx + 7] = ar0;
                    out_raw[idx + 8] = ara;
                    out_raw[idx + 9] = ada;
                    out_raw[idx +10] = ps[i].fluxes[j].flux;
                    out_raw[idx +11] = ps[i].fluxes[j].fluxerr;
                }
            }
        }
    }

    fiphot_free_photstar(ps, nstar);
    fiphot_free_inaps(inaps, ninap);

cleanup_sg:
    if (sg.coeff) free(sg.coeff);
    free(img.data);
    if (mk) free(mk);
    return r;
}

int fiphot_subtracted_photometry_cy(
    double *img_data, char *img_mask, int sx, int sy,
    double *star_x, double *star_y, int nstar,
    double *star_r0, double *star_ra, double *star_da, int nap,
    char *aperture_spec, int zoom,
    int bg_type, int bg_scatter, int bg_rejniter,
    double bg_rejlower, double bg_rejupper,
    int mask_ignore, int use_biquad, int use_sky, double sky_level,
    int is_disjoint_rings, int is_disjoint_apertures, double disjoint_radius,
    double **subpixel_data, int subg,
    int sg_order, double *sg_coeff, double sg_vmin,
    double correlation_length,
    char *kernel_spec, int normalize_kernel,
    weightlist *wl, int weightusage,
    double *out_flux, double *out_fluxerr,
    double *out_bgarea, double *out_bgflux, double *out_bgmedian, double *out_bgsigma,
    double *out_cntr_x, double *out_cntr_y, double *out_cntr_width,
    double *out_cntr_w_d, double *out_cntr_w_k,
    double *out_cntr_x_err, double *out_cntr_y_err, double *out_cntr_w_err,
    int *out_flag, int *out_rtot, int *out_rbad, int *out_rign,
    int *out_atot, int *out_abad)
{
    int i, r;
    image img;
    char **mk = NULL;
    photstar *ps;
    xphotpar xpp;
    spatialgain sg;
    apgeom *inaps = NULL;
    int ninap = 0, numap = 0;
    kernellist kl_data, *kl = NULL;

    if (!img_data || !star_x || !star_y || nstar <= 0 || sx <= 0 || sy <= 0) return 1;

    img.sx = sx; img.sy = sy;
    img.data = (double **)malloc(sy * sizeof(double *));
    for (i = 0; i < sy; i++) img.data[i] = img_data + i * sx;
    if (img_mask) {
        mk = (char **)malloc(sy * sizeof(char *));
        for (i = 0; i < sy; i++) mk[i] = img_mask + i * sx;
    }

    progbasename = "fitsh_cy.fiphot";
    is_verbose = 0; is_comment = 0;

    memset(&xpp, 0, sizeof(xphotpar));
    xpp.bgm.type = bg_type; xpp.bgm.scatter = bg_scatter;
    xpp.bgm.rejniter = bg_rejniter;
    xpp.bgm.rejlower = bg_rejlower; xpp.bgm.rejupper = bg_rejupper;
    xpp.maskignore = mask_ignore; xpp.use_biquad = use_biquad;
    xpp.wconfdist = 0.0; xpp.sky = sky_level; xpp.use_sky = use_sky;
    xpp.is_disjoint_rings = is_disjoint_rings;
    xpp.is_disjoint_apertures = is_disjoint_apertures;
    xpp.disjoint_radius = disjoint_radius;
    xpp.subpixeldata = subpixel_data; xpp.subg = subg;
    /* normalize and invert subpixel data (matches CLI main()) */
    if (xpp.subpixeldata && xpp.subg > 0) {
        normalize_subpixeldata(xpp.subpixeldata, xpp.subg);
        { int si, sj;
          for (si = 0; si < xpp.subg; si++)
          for (sj = 0; sj < xpp.subg; sj++)
              xpp.subpixeldata[si][sj] = 1.0 / xpp.subpixeldata[si][sj];
        }
    }

    sg.order = sg_order; sg.vmin = sg_vmin;
    if (sg_coeff && sg_order >= 0) {
        int ncoeff = (sg_order + 1) * (sg_order + 2) / 2;
        sg.coeff = (double *)malloc(ncoeff * sizeof(double));
        memcpy(sg.coeff, sg_coeff, ncoeff * sizeof(double));
    } else { sg.order = 0; sg.coeff = (double *)malloc(sizeof(double)); sg.coeff[0] = 1.0; }

    if (kernel_spec && kernel_spec[0]) {
        memset(&kl_data, 0, sizeof(kernellist));
        kl = &kl_data;
        if (create_kernels_from_kernelarg(kernel_spec, kl)) { r = 1; goto cleanup_sg2; }
        kernel_init_images(kl);
    }

    r = fiphot_build_inaps(&inaps, &ninap, &numap, aperture_spec, star_r0, star_ra, star_da, nap, zoom);
    if (r) { r = 1; goto cleanup_kernels; }

    ps = fiphot_build_stars(star_x, star_y, nstar, star_r0, star_ra, star_da, nap, star_r0 != NULL);

    r = do_photometry(&img, mk, ps, nstar, wl, weightusage, 0, inaps, numap, &sg, correlation_length, 0.0, &xpp);
    if (!r) {
        for (i = 0; i < nstar; i++) {
            ps[i].rfflux = ps[i].fluxes;
            ps[i].fluxes = NULL;
        }
    }

    r = do_subtracted_photometry(&img, mk, ps, nstar, inaps, numap, &sg, correlation_length, &xpp, kl, normalize_kernel);

    if (r == 0)
        fiphot_fill_photometry_output(ps, nstar, numap, out_flux, out_fluxerr, out_bgarea, out_bgflux, out_bgmedian, out_bgsigma, out_cntr_x, out_cntr_y, out_cntr_width, out_cntr_w_d, out_cntr_w_k, out_cntr_x_err, out_cntr_y_err, out_cntr_w_err, out_flag, out_rtot, out_rbad, out_rign, out_atot, out_abad, NULL, NULL, NULL);

    fiphot_free_photstar(ps, nstar);
    fiphot_free_inaps(inaps, ninap);

cleanup_kernels:
    if (kl && kl->nkernel > 0) {
        for (i = 0; i < kl->nkernel; i++) {
            if (kl->kernels[i].image) free(kl->kernels[i].image);
            if (kl->kernels[i].coeff) free(kl->kernels[i].coeff);
        }
        if (kl->kernels) free(kl->kernels);
    }
cleanup_sg2:
    if (sg.coeff) free(sg.coeff);
    free(img.data);
    if (mk) free(mk);
    return r;
}

int fiphot_magnitude_fit_cy(
    double *flux, int *flag, int nstar, int nap,
    double *ref_flux, double *ref_col_arr, double *ref_mag_arr, double *ref_err_arr,
    int *orders, int norder, int niter, double sigma,
    int sx, int sy,
    double *out_mag,
    int *out_ninit, int *out_nrejs, int *out_nstar, int *out_naperture)
{
    int i, j, r;
    photstar *ps;
    magfitparams mfp;
    magfitstat mfs;
    magflux mf0;

    if (!flux || !out_mag || nstar <= 0 || nap <= 0) return 1;
    if (norder < 0) return 0;

    for (i = 0; i < 4 && i <= norder; i++) mfp.orders[i] = orders ? orders[i] : 0;
    for (; i < 4; i++) mfp.orders[i] = 0;
    mfp.norder = norder; mfp.niter = niter; mfp.sigma = sigma;
    mf0.magnitude = 10.0; mf0.intensity = 10000.0;

    ps = (photstar *)calloc(nstar, sizeof(photstar));
    for (i = 0; i < nstar; i++) {
        ps[i].n = nap;
        ps[i].fluxes = (photflux *)calloc(nap, sizeof(photflux));
        if (ref_flux) {
            ps[i].rfflux = (photflux *)calloc(nap, sizeof(photflux));
            for (j = 0; j < nap; j++) {
                int idx = i * nap + j;
                ps[i].rfflux[j].flux = ref_flux[idx];
            }
        }
        ps[i].use_ref = (ref_mag_arr && ref_mag_arr[i] != 0.0) ? 1 : 0;
        ps[i].ref_mag = ref_mag_arr ? ref_mag_arr[i] : 0.0;
        ps[i].ref_col = ref_col_arr ? ref_col_arr[i] : 0.0;
        ps[i].ref_err = ref_err_arr ? ref_err_arr[i] : 0.0;
        for (j = 0; j < nap; j++) {
            int idx = i * nap + j;
            ps[i].fluxes[j].flux = flux[idx];
            ps[i].fluxes[j].flag = flag ? flag[idx] : 0;
        }
    }

    memset(&mfs, 0, sizeof(magfitstat));
    r = do_magnitude_fit(ps, nstar, &mfp, sx, sy, &mfs);

    if (r == 0) {
        calculate_magnitudes(ps, nstar, &mf0);
        for (i = 0; i < nstar; i++)
            for (j = 0; j < nap; j++)
                out_mag[i * nap + j] = ps[i].fluxes[j].mag;
    }
    if (out_ninit && mfs.ninit) for (j = 0; j < nap && j < 1; j++) out_ninit[j] = mfs.ninit[j];
    if (out_nrejs && mfs.nrejs) for (j = 0; j < nap && j < 1; j++) out_nrejs[j] = mfs.nrejs[j];
    if (out_nstar) *out_nstar = mfs.nstar;
    if (out_naperture) *out_naperture = mfs.naperture;

    for (i = 0; i < nstar; i++) {
        if (ps[i].fluxes) free(ps[i].fluxes);
        if (ps[i].rfflux) free(ps[i].rfflux);
    }
    free(ps);
    if (mfs.ninit) free(mfs.ninit);
    if (mfs.nrejs) free(mfs.nrejs);

    return r;
}
/* ---- commented out: dump for photstar comparison ----
static void _dump_photstar_bin_cmp(const char *path, photstar *ps, int np)
{
    ...
}
int fiphot_dump_raw_photstar_cython(...) { ... }
-------------------------------------------------------- */

int read_raw_photometry_cy(
    photstar **rps, int *rnp,
    double *raw_data, int nstar, int nap,
    double *ref_mag_arr, double *ref_col_arr, double *ref_err_arr)
{
    int i, j, idx;
    photstar *ps;
    photflux *rfflux;
    apgeom *inaps;

    if (!raw_data || !rps || !rnp || nstar <= 0 || nap <= 0) return 1;

    ps = (photstar *)calloc(nstar, sizeof(photstar));
    if (!ps) return 1;

    for (i = 0; i < nstar; i++) {
        photstar *p = &ps[i];

        p->x = raw_data[i * nap * 12 + 0 * 12 + 1];
        p->y = raw_data[i * nap * 12 + 0 * 12 + 2];
        p->n = nap;
        p->ninap = nap;

        if (ref_mag_arr && ref_mag_arr[i] != 0.0) {
            p->use_ref = 1;
            p->ref_mag = ref_mag_arr[i];
        } else {
            p->use_ref = 0;
            p->ref_mag = 0.0;
        }
        p->ref_col = ref_col_arr ? ref_col_arr[i] : 0.0;
        p->ref_err = ref_err_arr ? ref_err_arr[i] : 0.0;

        inaps = (apgeom *)malloc(sizeof(apgeom) * nap);
        for (j = 0; j < nap; j++) {
            memset(&inaps[j], 0, sizeof(apgeom));
            inaps[j].apgeom_type = APGEOM_TYPE_CIRCULAR;
            inaps[j].r0 = raw_data[i * nap * 12 + j * 12 + 7];
            inaps[j].ra = raw_data[i * nap * 12 + j * 12 + 8];
            inaps[j].da = raw_data[i * nap * 12 + j * 12 + 9];
        }
        p->inaps = inaps;

        rfflux = (photflux *)malloc(sizeof(photflux) * nap);
        for (j = 0; j < nap; j++) {
            idx = i * nap + j;
            rfflux[j].flux = raw_data[idx * 12 + 10];
            rfflux[j].fluxerr = raw_data[idx * 12 + 11];
            rfflux[j].flag = (int)raw_data[idx * 12 + 4];
            rfflux[j].bgarea = 0.0;
            rfflux[j].bgflux = 0.0;
            rfflux[j].bgmedian = 0.0;
            rfflux[j].bgsigma = 0.0;
            rfflux[j].mag = 0.0;
            rfflux[j].magerr = 0.0;
            rfflux[j].rtot = 0;
            rfflux[j].rbad = 0;
            rfflux[j].rign = 0;
            rfflux[j].atot = 0;
            rfflux[j].abad = 0;
        }
        p->rfflux = rfflux;

        p->fluxes = NULL;
        p->id = NULL;
    }

    *rps = ps;
    *rnp = nstar;
    return 0;
}

int fiphot_photometry_from_raw_cy(
    double *raw_data, int nstar, int nap,
    double *ref_mag_arr, double *ref_col_arr, double *ref_err_arr,
    double *img_data, char *img_mask, int sx, int sy,
    int bg_type, int bg_scatter, int bg_rejniter,
    double bg_rejlower, double bg_rejupper,
    int mask_ignore, int use_biquad, int use_sky, double sky_level,
    int is_disjoint_rings, int is_disjoint_apertures, double disjoint_radius,
    double **subpixel_data, int subg,
    int sg_order, double *sg_coeff, double sg_vmin,
    double correlation_length, double sigma,
    weightlist *wl, int weightusage,
    int nkernel,
    double ox, double oy, double scale, int ktype,
    int *k_types, int *k_orders, int *k_ncoeffs,
    int *k_hsizes, double *k_sigmas,
    int *k_bx, int *k_by,
    double **k_coeffs,
    char *kernel_spec, int normalize_kernel,
    double *out_flux, double *out_fluxerr,
    double *out_bgarea, double *out_bgflux, double *out_bgmedian, double *out_bgsigma,
    double *out_cntr_x, double *out_cntr_y, double *out_cntr_width,
    double *out_cntr_w_d, double *out_cntr_w_k,
    double *out_cntr_x_err, double *out_cntr_y_err, double *out_cntr_w_err,
    int *out_flag, int *out_rtot, int *out_rbad, int *out_rign,
    int *out_atot, int *out_abad)
{
    int i, j, r = 0;
    image img;
    char **mk = NULL;
    photstar *ps;
    int np;
    xphotpar xpp;
    spatialgain sg;
    kernellist kl_data, *kl = NULL;

    if (!img_data || !raw_data || nstar <= 0 || nap <= 0 || sx <= 0 || sy <= 0) return 1;

    img.sx = sx; img.sy = sy;
    img.data = (double **)malloc(sy * sizeof(double *));
    for (i = 0; i < sy; i++) img.data[i] = img_data + i * sx;
    if (img_mask) {
        mk = (char **)malloc(sy * sizeof(char *));
        for (i = 0; i < sy; i++) mk[i] = img_mask + i * sx;
    }

    progbasename = "fitsh_cy.fiphot";
    is_verbose = 0; is_comment = 0;

    memset(&xpp, 0, sizeof(xphotpar));
    xpp.bgm.type = bg_type; xpp.bgm.scatter = bg_scatter;
    xpp.bgm.rejniter = bg_rejniter;
    xpp.bgm.rejlower = bg_rejlower; xpp.bgm.rejupper = bg_rejupper;
    xpp.maskignore = mask_ignore; xpp.use_biquad = use_biquad;
    xpp.wconfdist = 0.0; xpp.sky = sky_level; xpp.use_sky = use_sky;
    xpp.is_disjoint_rings = is_disjoint_rings;
    xpp.is_disjoint_apertures = is_disjoint_apertures;
    xpp.disjoint_radius = disjoint_radius;
    xpp.subpixeldata = subpixel_data; xpp.subg = subg;
    if (xpp.subpixeldata && xpp.subg > 0) {
        normalize_subpixeldata(xpp.subpixeldata, xpp.subg);
        for (i = 0; i < xpp.subg; i++)
            for (j = 0; j < xpp.subg; j++)
                xpp.subpixeldata[i][j] = 1.0 / xpp.subpixeldata[i][j];
    }

    sg.order = sg_order; sg.vmin = sg_vmin;
    if (sg_coeff && sg_order >= 0) {
        int ncoeff = (sg_order + 1) * (sg_order + 2) / 2;
        sg.coeff = (double *)malloc(ncoeff * sizeof(double));
        memcpy(sg.coeff, sg_coeff, ncoeff * sizeof(double));
    } else { sg.order = 0; sg.coeff = (double *)malloc(sizeof(double)); sg.coeff[0] = 1.0; }

    if (nkernel > 0 && k_coeffs != NULL) {
        if (kernel_info_read_dicts(&kl_data, nkernel, ox, oy, scale, ktype,
            k_types, k_orders, k_ncoeffs, k_hsizes, k_sigmas,
            k_bx, k_by, k_coeffs))
            { r = 1; goto cleanup_sg; }
        kl = &kl_data;
    } else if (kernel_spec && kernel_spec[0]) {
        memset(&kl_data, 0, sizeof(kernellist));
        if (create_kernels_from_kernelarg(kernel_spec, &kl_data)) { r = 1; goto cleanup_sg; }
        kl_data.type = 0;
        kernel_init_images(&kl_data);
        kl = &kl_data;
    }

    read_raw_photometry_cy(&ps, &np, raw_data, nstar, nap, ref_mag_arr, ref_col_arr, ref_err_arr);

    // { FILE *_df = fopen("/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp/cy_do_sub_dump.bin", "wb");
    //   if (_df) { _dump_do_sub_input(_df, sx, sy, img.data, mk, ps, np, NULL, nap, &sg, correlation_length, &xpp, kl, normalize_kernel); fclose(_df); }
    // }

    r = do_subtracted_photometry(&img, mk, ps, np, NULL, nap, &sg, correlation_length, &xpp, kl, normalize_kernel);

    // { FILE *_df = fopen("/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp/cy_do_sub_dump.bin", "ab");
    //   if (_df) { _dump_do_sub_output(_df, ps, np); fclose(_df); }
    // }

    if (r == 0)
        fiphot_fill_photometry_output(ps, np, nap, out_flux, out_fluxerr, out_bgarea, out_bgflux, out_bgmedian, out_bgsigma, out_cntr_x, out_cntr_y, out_cntr_width, out_cntr_w_d, out_cntr_w_k, out_cntr_x_err, out_cntr_y_err, out_cntr_w_err, out_flag, out_rtot, out_rbad, out_rign, out_atot, out_abad, NULL, NULL, NULL);

    fiphot_free_photstar(ps, np);

cleanup_sg:
    if (sg.coeff) free(sg.coeff);
    if (kl && kl->nkernel > 0) {
        for (i = 0; i < kl->nkernel; i++) {
            if (kl->kernels[i].image) free(kl->kernels[i].image);
            if (kl->kernels[i].coeff) free(kl->kernels[i].coeff);
        }
        if (kl->kernels) free(kl->kernels);
    }
    free(img.data);
    if (mk) free(mk);
    return r;
}

int kernel_info_read_dicts(
    kernellist *kl,
    int nkernel,
    double ox, double oy, double scale, int ktype,
    int *k_types, int *k_orders, int *k_ncoeffs,
    int *k_hsizes, double *k_sigmas,
    int *k_bx, int *k_by,
    double **k_coeffs)
{
    int i, j;

    memset(kl, 0, sizeof(kernellist));
    kl->nkernel = nkernel;
    kl->ox = ox;
    kl->oy = oy;
    kl->scale = scale;
    kl->type = ktype;
    kl->kernels = (kernel *)calloc(nkernel, sizeof(kernel));

    for (i = 0; i < nkernel; i++)
    {
        kernel *k = &kl->kernels[i];
        k->type = k_types[i];
        k->order = k_orders[i];
        k->bx = k_bx ? k_bx[i] : 0;
        k->by = k_by ? k_by[i] : 0;
        k->flag = 0;
        k->target = 0;

        if (k->type == KERNEL_GAUSSIAN)
        {
            k->hsize = k_hsizes ? k_hsizes[i] : 0;
            k->sigma = k_sigmas ? k_sigmas[i] : 0.0;
        }

        if (k_ncoeffs && k_coeffs && k_coeffs[i])
        {
            int nv = k_ncoeffs[i];
            k->coeff = (double *)malloc(nv * sizeof(double));
            for (j = 0; j < nv; j++)
                k->coeff[j] = k_coeffs[i][j];
        }
    }

    kernel_init_images(kl);
    return 0;
}
