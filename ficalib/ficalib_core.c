/*****************************************************************************/
/* ficalib.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Command line tool for calibrating astronomical images.		     */
/*****************************************************************************/
#define	FITSH_FICALIB_VERSION	"0.9"
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "image.h"

#include "fitsh.h"

#include "mask.h"
#include "statistics.h"
// #include "io/iof.h" — removed for pyfitsh
// #include "io/scanarg.h" — removed for pyfitsh
#include "combine.h"
#include "tensor.h"
#include "common.h"
// #include "history.h" — removed for pyfitsh
#include "str.h"
#include "math/spline/spline.h"
#include "math/fit/lmfit.h"
#include "fbase.h"
#include "math/splinefit.h"
#include "math/polyfit.h"
// #include "io/tokenize.h" — removed for pyfitsh
// #include "longhelp.h" — removed for pyfitsh

/* fitsheaderset stubbed — pass NULL */
#ifdef  HAVE_NO_CC_EXTENSION
#define __extension__
#endif

/*****************************************************************************/

extern int is_verbose, is_comment;
extern char *progbasename;
// extern int fprint_error(char *expr,...);
// extern int fprint_warning(char *expr,...);
#include "stubs.h"

/*****************************************************************************/

// [81-385] removed for pyfitsh

#define		OVERSCAN_NONE		0
#define		OVERSCAN_SPLINE		1
#define		OVERSCAN_POLY		2

typedef struct
 {	int	x0,y0;
	int	sx,sy; 
 } trim;

typedef struct
 {	trim	area;
	int	type,order;
	int	niter;
	double	sigma;
 } overscan;

int overscan_model(int sy,double *spy,double *spw,double *spm,overscan *o,void *header)
{
 double	f,sig,sw,**fbase,w;
 int	i,j,n,k,order,nvar,used,avail;
 double	**amatrix,*bvector,*fvars;

 switch ( o->type )
  {  case OVERSCAN_SPLINE:
	order=o->order;
 	fbase=(double **)tensor_alloc_2d(double,sy,order+1);
	fbase_spline(fbase,order,sy);
	break;
     case OVERSCAN_POLY:
	order=o->order;
 	fbase=(double **)tensor_alloc_2d(double,sy,order+1);
	fbase_polynomial(fbase,order,sy);
	break;
     default:
	order=0;
 	fbase=(double **)tensor_alloc_2d(double,sy,order+1);
	fbase_polynomial(fbase,order,sy);
	break;
  }

 nvar=order+1;
 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);

 avail=0;
 for ( k=0 ; k<sy ; k++ )
  {	if ( spw[k] > 0.0 )
		avail++;
  }
 used=avail;
 
 for ( n=0 ; n<=(o->niter<0?0:o->niter) ; n++ )
  {	
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	amatrix[i][j]=0.0;	}
		bvector[i]=0.0;
	 }

	for ( k=0 ; k<sy ; k++ )
	 {	w=spw[k];
		for ( i=0 ; i<nvar ; i++ )
		 {	fvars[i]=fbase[i][k];	}
		f=spy[k];
		for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	amatrix[i][j]+=w*fvars[i]*fvars[j];	}
			bvector[i]+=w*f*fvars[i];
		 }
	 }

	solve_gauss(amatrix,bvector,nvar);

	for ( k=0 ; k<sy ; k++ )
	 {	f=0.0;
		for ( i=0 ; i<nvar ; i++ )
		 {	f+=bvector[i]*fbase[i][k];	}
		spm[k]=f;
	 }

	if ( n < o->niter )
	 {	sig=0.0;
		sw=0.0;
		for ( k=0 ; k<sy ; k++ )
		 {	f=spy[k]-spm[k];
			w=spw[k];
			sig+=w*f*f;
			sw +=w;
		 }
		sig/=sw;
		if ( sig>0.0 )	sig=sqrt(sig);
		else		sig=0.0;
		for ( k=0 ; k<sy ; k++ )
		 {	f=spy[k]-spm[k];
			if ( fabs(f) >= o->sigma * sig )
			 {	spw[k]=0.0;
				used--;
			 }
		 }
	 }
  }

 sig=0.0;
 sw=0.0;
 for ( k=0 ; k<sy ; k++ )
  {	f=spy[k]-spm[k];
	w=spw[k];
	sig+=w*f*f;
	sw +=w;
  }
 sig/=sw;
 if ( sig>0.0 )	sig=sqrt(sig);
 else		sig=0.0;

 /* just export the overscan information */
 if ( header != NULL )
  {	char	buff[80],*type;
	int	l;
	switch ( o->type )
	 {   case OVERSCAN_SPLINE:
		type="spline";
		break;
	     case OVERSCAN_POLY:
		type="polynomial";
		break;
	     default:
		type="const";
		break;
	 }
	l=snprintf(buff,80,"ficalib: %s o=%d c=",type,o->order);
	for ( i=0 ; i<nvar ; i++ )
	 {	if ( i )
			l+=snprintf(buff+l,(80-l>0?80-l:0),",");
		l+=snprintf(buff+l,(80-l>0?80-l:0),"%.1f",bvector[i]);
	 }
	l+=snprintf(buff+l,(80-l>0?80-l:0)," s=%.1f n=%d/%d/%d",sig,sy,sy-avail,sy-used);

	/* fits_headerset_set_string removed for pyfitsh */
  }

 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 tensor_free(fbase);

 return(0);
}

int overscan_do_vertical(image *img,char **mask,int x0,int y0,int sx,int sy,int ox,int wx,overscan *o,combine_parameters *cp,void *header)
{
 int	i,j,lp,np,x,y;
 double	*line,*spw,*spy,*spm;

 spw =(double *)malloc(sizeof(double)*sy);
 spy =(double *)malloc(sizeof(double)*sy);
 line=(double *)malloc(sizeof(double)*wx);

 np=0;
 for ( i=0 ; i<sy ; i++ )
  {	lp=0;
	for ( j=0 ; j<wx ; j++ )
	 {	y=y0+i;
		x=ox+j;
		if ( x<0 || x>=img->sx || y<0 || y>=img->sy )
			continue;
		else if ( mask != NULL && mask[y][x] ) 
			continue;
		line[lp]=img->data[y][x];
		lp++;
	 }
	if ( lp>0 )
	 {	spw[i]=1.0;
		spy[i]=combine_points(line,lp,cp);
		np++;
	 }
	else
	 {	spw[i]=0.0;
		spy[i]=0.0;
	 }
  }
 
 if ( np < o->order+1 )	/* overscan failed due to insufficient num of points */
  {	free(line);
	free(spy);
	free(spw);
	return(-1);	
  }

 spm  =(double *)malloc(sizeof(double)*sy);

 overscan_model(sy,spy,spw,spm,o,header);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[y0+i][x0+j] -= spm[i];	}
  }

 free(spm);

 free(line);
 free(spy);
 free(spw);
 
 return(0);
}

int overscan_do_horizontal(image *img,char **mask,int x0,int y0,int sx,int sy,int oy,int wy,overscan *o,combine_parameters *cp,void *header)
{
 int	i,j,lp,np,x,y;
 double	*line,*spw,*spy,*spm;

 spw =(double *)malloc(sizeof(double)*sx);
 spy =(double *)malloc(sizeof(double)*sx);
 line=(double *)malloc(sizeof(double)*wy);

 np=0;
 for ( i=0 ; i<sx ; i++ )
  {	lp=0;
	for ( j=0 ; j<wy ; j++ )
	 {	x=x0+i;
		y=oy+j;
		if ( x<0 || x>=img->sx || y<0 || y>=img->sy )
			continue;
		else if ( mask != NULL && mask[y][x] ) 
			continue;
		line[lp]=img->data[y][x];
	 }
	if ( lp>0 )
	 {	spw[i]=1.0;
		spy[i]=combine_points(line,lp,cp);
		np++;
	 }
	else
	 {	spw[i]=0.0;
		spy[i]=0.0;
	 }
  }

 if ( np < o->order+1 )	/* overscan failed due to insufficient num of points */
  {	free(line);
	free(spy);
	free(spw);
	return(-1);	
  }

 spm  =(double *)malloc(sizeof(double)*sx);

 overscan_model(sx,spy,spw,spm,o,header);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[y0+i][x0+j] -= spm[j];	}
  }

 free(spm);

 free(line);
 free(spy);
 free(spw);
 
 return(0);
}

int overscan_check(trim *area,trim *o)
{
 if ( o->y0==area->y0 && o->sy==area->sy && 
 abs(2*(o->x0-area->x0)+(o->sx-area->sx))>=o->sx+area->sx ) 
	return(+1);	/* vertical overscan	*/
 else if ( o->x0==area->x0 && o->sx==area->sx && 
 abs(2*(o->y0-area->y0)+(o->sy-area->sy))>=o->sy+area->sy )
	return(-1);	/* horizontal overscan	*/
 else
	return(0);	/* invalid overscan	*/
}

int overscan_correction(image *img,char **mask,trim *area,overscan *scans,int nscan,combine_parameters *cp,void *header)
{
 int		ret,r,i;
 overscan	*o;

 ret=0;
 for ( i=0 ; i<nscan ; i++ )
  {	o=&scans[i];
	if ( o->area.y0==area->y0 && o->area.sy==area->sy && 
	abs(2*(o->area.x0-area->x0)+(o->area.sx-area->sx))>=o->area.sx+area->sx )
	 {	r=overscan_do_vertical(img,mask,area->x0,area->y0,area->sx,area->sy,o->area.x0,o->area.sx,o,cp,header);
		if ( ! r )	ret++;
	 }
	else if ( o->area.x0==area->x0 && o->area.sx==area->sx && 
	abs(2*(o->area.y0-area->y0)+(o->area.sy-area->sy))>=o->area.sy+area->sy )
	 {	r=overscan_do_horizontal(img,mask,area->x0,area->y0,area->sx,area->sy,o->area.y0,o->area.sy,o,cp,header);
		if ( ! r )	ret++;
	 }

  }

 if ( header != NULL )
  {	char	buff[80];
	snprintf(buff,80,"ficalib: total_overscans=%d/%d",ret,nscan);
	/* fits_headerset_set_string removed for pyfitsh */
  }
 
 return(ret);
}

/*****************************************************************************/

#include "ficalib_core.h"

int ficalib_calibrate_cy(
    double *img_data, unsigned char *mask_data, int sx, int sy,
    double *bias_data, unsigned char *bias_mask,
    int bias_sx, int bias_sy,
    double *dark_data, unsigned char *dark_mask,
    int dark_sx, int dark_sy,
    double *flat_data, unsigned char *flat_mask,
    int flat_sx, int flat_sy,
    double flat_mean, double dark_time,
    double **out_data, unsigned char **out_mask, int *out_sx, int *out_sy,
    ficalib_flatpoly *out_flat)
{
    int i, j;
    double *out;
    unsigned char *omask;

    if (!img_data || sx <= 0 || sy <= 0) return 1;

    out = (double *)malloc(sx * sy * sizeof(double));
    omask = (unsigned char *)calloc(sx * sy, sizeof(unsigned char));
    if (!out || !omask) { free(out); free(omask); return 1; }

    /* copy input */
    memcpy(out, img_data, sx * sy * sizeof(double));
    if (mask_data)
        memcpy(omask, mask_data, sx * sy);

    /* bias subtraction */
    if (bias_data && bias_sx > 0 && bias_sy > 0 && bias_sx == sx && bias_sy == sy) {
        for (i = 0; i < sx * sy; i++) out[i] -= bias_data[i];
    }

    /* dark subtraction */
    if (dark_data && dark_sx > 0 && dark_sy > 0 && dark_sx == sx && dark_sy == sy) {
        double dscale = (dark_time > 0.0) ? dark_time : 1.0;
        for (i = 0; i < sx * sy; i++) out[i] -= dark_data[i] * dscale;
    }

    /* flat field correction */
    if (flat_data && flat_sx > 0 && flat_sy > 0 && flat_sx == sx && flat_sy == sy) {
        double fmean = 0.0;
        int fcount = 0;
        for (i = 0; i < sx * sy; i++) {
            if (!flat_mask || !flat_mask[i]) {
                fmean += flat_data[i];
                fcount++;
            }
        }
        if (fcount > 0) fmean /= (double)fcount;
        if (flat_mean <= 0.0) flat_mean = fmean;
        if (fmean > 0.0 && flat_mean > 0.0) {
            for (i = 0; i < sx * sy; i++) {
                if (!flat_mask || !flat_mask[i]) {
                    out[i] *= flat_mean / flat_data[i];
                }
            }
        }
    }

    *out_data = out;
    *out_mask = omask;
    *out_sx = sx;
    *out_sy = sy;
    if (out_flat) memset(out_flat, 0, sizeof(ficalib_flatpoly));
    return 0;
}
