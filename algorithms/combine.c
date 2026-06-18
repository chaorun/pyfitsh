/*****************************************************************************/
/* combine.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to combining more FITS images.			     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "image.h"

#include "fitsh.h"

#include "mask.h"
// // #include "io/scanarg.h" removed  /* removed for pyfitsh */
#include "statistics.h"
#include "tensor.h"
#include "combine.h"

/*****************************************************************************/

int combine_parameters_reset(combine_parameters *cp)
{
 cp->mode=0;
 cp->niter=0;
 cp->lower=cp->upper=0.0;
 cp->lowest=cp->highest=0;
 cp->ignore_flag=0;
 cp->logicalmethod=0;

 return(0);
}

// combine_parse_mode removed — uses io/scanarg.h
// end of removed function

/*****************************************************************************/

double do_ordered_rejection(double *points,int n,double m,double w1,double w2,
int niter,double lower,double upper,int is_median)
{
 double	s,s2,lim;
 int	r;

 if ( n<=0 || points==NULL )
	return(0.0);

 for ( ; niter > 0 ; niter-- )
  {	
	s2=w2+(double)n*m*(m-2.0*w1);
	if ( s2 <= 0.0 )
		return(m);
	else
		s=sqrt(s2/(double)n);

	lim=m-s*lower;
	r=0;
	while ( n>0 && *points < lim )
	 {	w1-=*points,
		w2-=(*points)*(*points);
		points++,n--;
		r++;
	 }
	lim=m+s*upper;
	while ( n>0 && points[n-1] > lim )
	 {	w1-=points[n-1],
		w2-=points[n-1]*points[n-1];
		n--;
		r++;
	 }
	if ( r<=0 )
		return(m);
	else if ( n<=2 )
		return(m);
	else if ( ! is_median )
		m=w1/(double)n;
	else if ( n%2==1 )
		m=points[n/2];
	else
		m=0.5*(points[n/2-1]+points[n/2]);

  }

 return(m);
}

double combine_points(double *points,int n,combine_parameters *cp)
{
 int	i;
 double w,w1,w2,m;
 
 median(points,n);

 if ( n==0 )				/* All pixels to combine are false.  */
	return(0.0);		

 else if ( n==1 )			/* One pixel.			     */
	return(points[0]);
 else if ( cp->mode==COM_MODE_MIN )
  {	w=points[0];
	for ( i=1 ; i<n ; i++ )
	 {	if ( points[i]<w )	w=points[i];	}
	return(w);
  }
 else if ( cp->mode==COM_MODE_MAX )
  {	w=points[0];
	for ( i=1 ; i<n ; i++ )
	 {	if ( w<points[i] )	w=points[i];	}
	return(w);
  }
 
 else if ( n==2 )			/* Two, it's the mean/median of them */
	return(0.5*(points[0]+points[1]));

 else if ( cp->mode==COM_MODE_MEAN )
  {	w=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w+=points[i];		}
	w=w/(double)n;
	return(w);
  }
 else if ( cp->mode==COM_MODE_SUM )
  {	w=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w+=points[i];		}
	return(w);
  }
 else if ( cp->mode==COM_MODE_SQSUM )
  {	w=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w+=points[i]*points[i];		}
	return(w);
  }
 else if ( cp->mode==COM_MODE_SCT )
  {	double	s,s2;
	s=s2=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w=points[i];
		s+=w,s2+=w*w;
	 }
	s/=(double)n,s2/=(double)n;
	s2=s2-s*s;
	if ( s2>0.0 )	return(sqrt(s2));
	else		return(0.0);
  }
 else if ( cp->mode==COM_MODE_MEDIAN )
  {	if ( n%2==1 )	w=points[n/2];
	else		w=0.5*(points[n/2-1]+points[n/2]);
	return(w);
  }
 else if ( cp->mode==COM_MODE_REJ_DEPRECATED )	/* deprecated... see below   */
  {	if ( n==3 )			/* Three, also...		     */
		return(points[1]);
	else if ( n<=5 )		/* Average of all except the first   */
	 {	w=0.0;			/* and the last one.		     */
		for ( i=1 ; i<n-1 ; i++ )
		 {	w+=points[i];		}
		w=w/(double)(n-2);
		return(w);
	 }
	else				/* Average all of them except the    */
	 {	w=0.0;			/* first two and last two ones.	     */
		for ( i=2 ; i<n-2 ; i++ )
		 {	w+=points[i];		}
		w=w/(double)(n-4);
		return(w);
	 }
  }
 else if ( cp->mode==COM_MODE_REJ_MEAN )
  {	w1=0.0;
	w2=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w1+=points[i];
		w2+=points[i]*points[i];
	 }
	w=do_ordered_rejection(points,n,w1/(double)n,w1,w2,cp->niter,cp->lower,cp->upper,0);
	return(w);
  }
 else if ( cp->mode==COM_MODE_REJ_MEDIAN )
  {	w1=0.0;
	w2=0.0;
	for ( i=0 ; i<n ; i++ )
	 {	w1+=points[i];
		w2+=points[i]*points[i];
	 }
	if ( n%2==1 )	m=points[n/2];
	else		m=0.5*(points[n/2-1]+points[n/2]);
	w=do_ordered_rejection(points,n,m,w1,w2,cp->niter,cp->lower,cp->upper,1);
	return(w);
  }
 else if ( cp->mode==COM_MODE_TRUNC_MEAN )
  {	int	c;
	w=0.0;
	for ( c=0,i=cp->lowest ; i<n-cp->highest ; i++ )
	 {	w+=points[i];c++;		}
	if ( 0<c )	w=w/(double)c;
	else		w=0.0;
	return(w);
  }
 else if ( cp->mode==COM_MODE_WINS_MEAN )
  {	int	c;
	w=0.0;
	for ( c=0,i=cp->lowest ; i<n-cp->highest ; i++ )
	 {	w+=points[i];c++;		}
	if ( 0<c )
	 {	w+=points[cp->lowest]*cp->lowest;
		w+=points[n-cp->highest-1]*cp->highest;
		w=w/(double)n;
	 }
	else
		w=0.0;
	return(w);
  }
 else					/* Unimplemented method...?!	     */
	return(0.0);
}

int combine_lines(double **lines,int n,int sx,double *out,
	combine_parameters *cp,char **wmask,char *outmask)
{
 int	i,j,k,mw,ow;
 double	*warr,w;

 warr=(double *)tensor_alloc_1d(double,n);

 for ( i=0 ; i<sx ; i++ )
  {	k=0;
	ow=0;
	if ( wmask != NULL )
	 {	for ( j=0 ; j<n ; j++ )
		 {	w =lines[j][i];
			mw=wmask[i][j];
			ow|=mw;
			if ( (cp->ignore_flag&COM_IGNORE_NEGATIVE) && w<0.0 )
				mw|=MASK_FAULT;
			if ( ! mw )
				warr[k]=w,k++;
		 }
	 }
	else
	 {	for ( j=0 ; j<n ; j++ )
		 {	w =lines[j][i];
			mw=0;
			if ( (cp->ignore_flag&COM_IGNORE_NEGATIVE) && w<0.0 )
				mw|=MASK_FAULT;
			if ( ! mw )
				warr[k]=w,k++;
		 }
	 }
	if ( outmask != NULL )
	 {	if ( ! cp->logicalmethod && k>0 )		/* --logical-or	 */
			out[i]=combine_points(warr,k,cp),outmask[i]=0;
		else if ( cp->logicalmethod && k==n )		/* --logical-and */
			out[i]=combine_points(warr,k,cp),outmask[i]=0;
		else
		 {	out[i]=0.0;
			outmask[i]=(ow?ow:MASK_FAULT);
		 }
	 }
	else
	 {	if ( ! cp->logicalmethod && k>0 )		/* --logical-or	 */
			out[i]=combine_points(warr,k,cp);
		else if ( cp->logicalmethod && k==n )		/* --logical-and */
			out[i]=combine_points(warr,k,cp);
		else
			out[i]=0.0;
	 }
  }

 tensor_free(warr);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// combine_images_from_files* and combine_cleanup removed — file I/O
