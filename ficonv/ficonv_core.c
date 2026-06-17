/*****************************************************************************/
/* ficonv_core.c -- complete kernel fitting pipeline from ficonv.c       */
/*                                                                           */
/* Changes from original ficonv.c:                                           */
/*   - Removed CLI main() / scanarg / usage / long-help                      */
/*   - Added flat-array entry point ficonv_fit_cy_v2                       */
/*   - Removed file I/O includes (io/iof.h, io/scanarg.h)                    */
/*   - Removed longhelp.h, history.h, common.h (not needed)                  */
/*   - All algorithm functions verbatim from ficonv.c                        */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h> /* for mkdir, used by dump routines below */
#include <stdarg.h>
#include <time.h>

#include "image.h"

#include "fitsh.h"

#include "math/spline/spline.h"
#include "math/spline/biquad.h"
#include "math/spline/bicubic.h"
#include "statistics.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "mask.h"
#include "tokenize.h"

#include "tensor.h"
#include "kernel.h"

#ifdef  HAVE_NO_CC_EXTENSION
#define __extension__
#endif

int kernel_info_read_dicts_in_ficonv(
    kernellist *kl,
    int nkernel,
    double ox, double oy, double scale, int ktype,
    int *k_types, int *k_orders, int *k_ncoeffs,
    int *k_hsizes, double *k_sigmas,
    int *k_bx, int *k_by,
    double **k_coeffs);

/*****************************************************************************/

extern int	is_comment, is_verbose;
extern char	*progbasename;



/*****************************************************************************/

#define		FM_NORMAL		0x00
#define		FM_MASKED		0x01
#define		FM_WEIGHTED		0x02
#define		FM_BGITERATIVE		0x04

/*****************************************************************************/

int affine_img_transformation_spline_fit(image *ref,image *img,char **mask,int bx,int by,double *coeff)
{
 int	i,j,k,l,v,sx,sy,nvar;
 double	**xbase,**ybase;
 double	**xnspc,**ynspc;
 double	*tmp,x,y;
 double	**amatrix,*bvector,*fvars,f,w;

 if ( ref==NULL || img==NULL || coeff==NULL )
	return(-1);
 else if ( (sx=ref->sx) != img->sx || (sy=ref->sy) != img->sy )
	return(-1);
 else if ( bx<=0 || by<=0 )
	return(-1);

 nvar=2*(bx+1)*(by+1);

 xbase=tensor_alloc_2d(double,sx,bx+1);
 ybase=tensor_alloc_2d(double,sy,by+1);
 xnspc=tensor_alloc_2d(double,bx+1,bx+1);
 ynspc=tensor_alloc_2d(double,by+1,by+1);
 tmp  =tensor_alloc_1d(double,(bx>by?bx:by)+1);

 for ( i=0 ; i<=by ; i++ )
  {	for ( k=0 ; k<=by ; k++ )
 	 {	tmp[k]=0.0;		}
 	tmp[i]=1.0;
 	natspline_coeff(tmp,by+1,ynspc[i]);
 	for ( k=0 ; k<sy ; k++ )
 	 {	y=(double)by*((double)k+0.5)/(double)sy;
 		ybase[i][k]=natspline_inter(tmp,ynspc[i],by+1,y);
 	 }
  }

 for ( j=0 ; j<=bx ; j++ )
  {	for ( l=0 ; l<=bx ; l++ )
 	 {	tmp[l]=0.0;		}
 	tmp[j]=1.0;
 	natspline_coeff(tmp,bx+1,xnspc[j]);
 	for ( l=0 ; l<sx ; l++ )
 	 {	x=(double)bx*((double)l+0.5)/(double)sx;
 		xbase[j][l]=natspline_inter(tmp,xnspc[j],bx+1,x);
 	 }
  }

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 fvars  =vector_alloc(nvar);

 for ( i=0 ; i<nvar ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
 	 {	amatrix[i][j]=0.0;		}
 	bvector[i]=0.0;
  }

 for ( k=0 ; k<sy ; k++ )
  {	for ( l=0 ; l<sx ; l++ )
 	 {	if ( mask != NULL && mask[k][l] )
 			continue;
 		v=0;
 		w=ref->data[k][l];
 		for ( i=0 ; i<=by ; i++ )
 		 {	for ( j=0 ; j<=bx ; j++ )
 			 {	f=xbase[j][l]*ybase[i][k];
 				fvars[v+0]=f*w;
 				fvars[v+1]=f;
 				v+=2;
 			 }
 		 }
 		for ( i=0 ; i<nvar ; i++ )
 		 {	for ( j=0 ; j<nvar ; j++ )
 			 {	amatrix[i][j] += fvars[i]*fvars[j];	}
 			bvector[i] += fvars[i]*img->data[k][l];
 		 }
 	 }
  }

 solve_gauss(amatrix,bvector,nvar);

 for ( i=0 ; i<nvar ; i++ )
 	coeff[i]=bvector[i];

 vector_free(fvars);
 vector_free(bvector);
 matrix_free(amatrix);

 tensor_free(tmp);
 tensor_free(ynspc);
 tensor_free(xnspc);
 tensor_free(ybase);
 tensor_free(xbase);

 return(0);
}

int affine_img_transformation_spline_eval(image *ref,char **mask,int bx,int by,double *coeff)
{
 double	**bcscl,**bcoff,**tbscl,**tboff;
 int	i,j,k,l,v,sx,sy;
 double	scl,off,x,y;

 bcscl=(double **)tensor_alloc_2d(double,2*(bx+1),2*(by+1));
 bcoff=(double **)tensor_alloc_2d(double,2*(bx+1),2*(by+1));
 tbscl=(double **)tensor_alloc_2d(double,bx+1,by+1);
 tboff=(double **)tensor_alloc_2d(double,bx+1,by+1);

 v=0;
 for ( i=0 ; i<=by ; i++ )
  {	for ( j=0 ; j<=bx ; j++ )
 	 {	tbscl[i][j]=coeff[v+0];
 		tboff[i][j]=coeff[v+1];
 		v+=2;
 	 }
  }

 bicubic_coeff(tbscl,bx+1,by+1,bcscl,NULL);
 bicubic_coeff(tboff,bx+1,by+1,bcoff,NULL);

 tensor_free(tboff);
 tensor_free(tbscl);

 sx=ref->sx;
 sy=ref->sy;

 for ( k=0 ; k<sy ; k++ )
  {	y=(double)by*((double)k+0.5)/(double)sy;

 	for ( l=0 ; l<sx ; l++ )
 	 {	if ( mask != NULL && mask[k][l] )
 			continue;

 		x=(double)bx*((double)l+0.5)/(double)sx;

 		scl=bicubic_inter(bcscl,x,y);
 		off=bicubic_inter(bcoff,x,y);
 	
 		ref->data[k][l]=ref->data[k][l]*scl+off;
 	 }
  }
 
 tensor_free(bcoff);
 tensor_free(bcscl);

 return(0);
}

/*****************************************************************************/

char ** create_level_mask(image *img,char **inmask,double mlev1,double klev1,double mlev2,double klev2)
{
 char	**mask;
 double	**bqc,*ps,sc,th1,th2;
 int	i,j,sx,sy,np;

 sx=img->sx,
 sy=img->sy;
 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 biquad_coeff(img->data,sx,sy,bqc,inmask);
 ps=tensor_alloc_1d(double,sx*sy);
 mask=mask_duplicate(inmask,sx,sy);
 np=0;
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
 	 {	if ( mask[i][j] )	continue;
 		sc=biquad_scatter(bqc,j,i);
 		bqc[2*i+1][2*j+1]=ps[np]=sc;
 		np++;
 	 }
  }
 median(ps,np);

 if ( klev1>0.0 ) i=(int)((double)(np-1)*mlev1),th1=ps[i]*klev1;
 else		  th1=0.0;
 if ( klev2>0.0 ) i=(int)((double)(np-1)*mlev2),th2=ps[i]*klev2;
 else		  th2=0.0;

 tensor_free(ps);

 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
 	 {	if ( mask[i][j] )	continue;
 		sc=bqc[2*i+1][2*j+1];
 		if ( sc<th1 && th1>0.0 )	mask[i][j]=-1;
 		else if ( sc>th2 && th2>0.0 )	mask[i][j]=+1;
 		else				mask[i][j]= 0;
 	 }
  }	
 tensor_free(bqc);
 return(mask);
}

char ** create_foreground_mask(image *img,char **inmask,double mlev,double klev)
{
 char	**ret;
 ret=create_level_mask(img,inmask,0.0,0.0,mlev,klev);
 return(ret); 
}
char ** create_background_mask(image *img,char **inmask,double mlev,double klev)
{
 char	**ret;
 int	i,j;
 ret=create_level_mask(img,inmask,mlev,klev,0.0,0.0);
 for ( i=0 ; i<img->sy ; i++ )
  {	for ( j=0 ; j<img->sx ; j++ )
 	 {	if ( ret[i][j]<0 )	ret[i][j]=1;	}
  }
 return(ret); 
}

double **create_image_weight(image *img,char **inmask)
{
 double	**bqc,**wght;
 int	i,j,sx,sy;
 sx=img->sx,sy=img->sy;
 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 biquad_coeff(img->data,sx,sy,bqc,inmask);
 wght=(double **)tensor_alloc_2d(double,sx,sy);
 if ( inmask==NULL )
  {	for ( i=0 ; i<sy ; i++ )
 	 {	for ( j=0 ; j<sx ; j++ )
 		 {	wght[i][j]=biquad_scatter(bqc,j,i);	 }
 	 }
  }
 else
  {	for ( i=0 ; i<sy ; i++ )
 	 {	for ( j=0 ; j<sx ; j++ )
 		 {	if ( inmask[i][j] )
 				wght[i][j]=0.0;
 			else
 				wght[i][j]=biquad_scatter(bqc,j,i);
 		 }	
 	 }
  }
 tensor_free(bqc);
 return(wght);
}

/*****************************************************************************/

int mark_stamps(image *img,rectan *rcts,int nrct)
{
 int	i;
 rectan	*r; 
 double	color;

 color=10000.0;
 for ( i=0 ; i<nrct ; i++ )
  {	r=&rcts[i];
 	image_draw_line(img,r->x1  ,r->y1  ,r->x2-1,r->y1  ,color,0x3333);
 	image_draw_line(img,r->x2-1,r->y1  ,r->x2-1,r->y2-1,color,0x3333);
 	image_draw_line(img,r->x2-1,r->y2-1,r->x1  ,r->y2-1,color,0x3333);
 	image_draw_line(img,r->x1  ,r->y2-1,r->x1  ,r->y1  ,color,0x3333);
  }
 return(0);
}

/*****************************************************************************/

int make_subtracted_image(image *cnv,char **mask,image *img,char **mask_img,image *out,char **mask_out,kernellist *xlist)
{
 image	xkcimg;
 int		i,j,sx,sy;
 char		**imask,**aimask;

 if ( cnv==NULL || cnv->data==NULL )	return(1);
 if ( img==NULL || img->data==NULL )	return(1);
 if ( out==NULL || out->data==NULL )	return(1);
 sx=cnv->sx,sy=cnv->sy;

 /* this xlist != NULL part is needed to be obsoleted sooner-or-later: */ 
 if ( xlist != NULL && xlist->kernels != 0 && xlist->nkernel>0 )
  {	kernel	*k;
 	int	ihsize;

 	ihsize=0;
 	for ( i=0,k=xlist->kernels ; i<xlist->nkernel ; i++,k++ )
 	 {	if ( k->hsize > ihsize ) ihsize=k->hsize;	}
 	imask=aimask=mask_expand_false(mask_img,sx,sy,ihsize,-1,-1,1);
 	image_duplicate(&xkcimg,img,1);
 	convolve_with_kernel_set(img,imask,xlist,&xkcimg);
  }
 else	
  {	xkcimg.data=NULL;
 	aimask=NULL;
 	imask=mask_img;
  }
		
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
 	 {	if ( ! mask[i][j] && ! imask[i][j] )
 		 {	out->data[i][j]=img->data[i][j]-cnv->data[i][j];
 			if ( xkcimg.data != NULL )
 				out->data[i][j]+=xkcimg.data[i][j];
 			mask_out[i][j]=0;
 		 }
 		else
 		 {	out->data[i][j]=0.0;
 			mask_out[i][j]=-1;
 		 }
 	 }
  }
 if ( xkcimg.data != NULL )	image_free(&xkcimg);
 if ( aimask != NULL )		mask_free(aimask);

 return(0);
}

/*****************************************************************************/

int create_weights(image *img,char **mask,char **levmask,double **wcc,int sign)
{
 int	i,j,sx,sy;
 sx=img->sx,sy=img->sy;
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
 	 {	if ( mask[i][j] )			wcc[i][j]=0.0;
 		else if ( sign>0 && levmask[i][j]>0 )	wcc[i][j]=1.0;
 		else if ( sign<0 && levmask[i][j]<0 )	wcc[i][j]=1.0;
 		else if ( sign==0 && levmask[i][j] )	wcc[i][j]=1.0;
 		else					wcc[i][j]=0.0;
 	 }
  }
 return(0);
}

int fit_kernels_native(image *ref,image *img,char **mask_conv,
	kernellist *klist,kernellist *xlist,int method,int bdc)
{
 int	i,j,l,sx,sy;
 kernel	*bgkernel,*k;

 sx=ref->sx,sy=ref->sy;

 bgkernel=NULL;
 for ( i=0,k=klist->kernels ; i<klist->nkernel ; i++,k++ )
  {	if ( k->type==KERNEL_BACKGROUND )	bgkernel=k;
 	k->flag=1;
  }
 if ( xlist != NULL )
  {	for ( i=0,k=xlist->kernels ; i<xlist->nkernel ; i++,k++ )
 	 {	k->flag=1;	}
  }

 if ( method & FM_WEIGHTED )	method|=FM_MASKED;	/* -w implies -m ! */

 if ( bgkernel != NULL && ( method & FM_MASKED ) )
  {	double		**wcc,w;
 	kernel		*k;
 	char		**fgmask;
 	
 	wcc=tensor_alloc_2d(double,sx,sy);

 	fgmask=create_foreground_mask(ref,mask_conv,0.5,3.0);

 	for ( i=0 ; i<klist->nkernel ; i++ )
 	 {	k=&klist->kernels[i];k->flag=1;		}

 	l=0;
 	for ( i=0 ; i<sy ; i++ )
 	 {	for ( j=0 ; j<sx ; j++ )
 		 {	if ( mask_conv[i][j] )
 				wcc[i][j]=0.0;
 			else if ( fgmask[i][j] )
 		 	 {	if ( method & FM_WEIGHTED )
 					w=ref->data[i][j];
 				else	
 					w=1.0;					
 				wcc[i][j]=w;
 				l++;
 			 }
 			else	wcc[i][j]=0.0;
 		 }
 	  }		
 	logmsg(is_verbose>=1,"Foreground pixels: %d/%d\n",l,sx*sy);
 	fit_kernel_poly_coefficients_block(ref,img,mask_conv,wcc,bdc,bdc,klist,xlist);
 	mask_free(fgmask);
 	tensor_free(wcc);
  }

 else
  	fit_kernel_poly_coefficients_block(ref,img,mask_conv,NULL,bdc,bdc,klist,xlist);

 return(0);
}


int fit_kernels(image *ref,char **mask_ref,image *img,char **mask_img,
	char **inmask,kernellist *klist,kernellist *xlist,int method,int bdc,
	int niter,double rlevel,double gain)
{
 int	hsize,i,sx,sy;
 char	**mask,**mask_conv;
 kernel	*k;

 if ( ref==NULL || img==NULL )		return(1);
 if ( ref->data==NULL || img->data==NULL )	return(1);
 sx=ref->sx,sy=ref->sy;
 if ( sx != img->sx || sy != img->sy )		return(1);

 if ( xlist==NULL || xlist->nkernel<=0 || xlist->kernels==NULL )
	xlist=NULL;

 logmsg(is_verbose>=1,"Fitting kernel coefficients ...\n");

 hsize=0;
 for ( i=0,k=klist->kernels ; i<klist->nkernel ; i++,k++ )
  {	if ( k->hsize > hsize )		hsize=k->hsize;		}
 if ( xlist != NULL )
  {	for ( i=0,k=xlist->kernels ; k != NULL && i<xlist->nkernel ; i++,k++ )
	 {	if ( k->hsize > hsize )	hsize=k->hsize;		}
  }

 mask=mask_create_empty(sx,sy);
 mask_and(mask,sx,sy,mask_img);
 mask_and(mask,sx,sy,mask_ref);
 if ( inmask != NULL )
 	mask_and(mask,sx,sy,inmask);
 mask_conv=mask_expand_false(mask,sx,sy,hsize,-1,-1,1);

 if ( niter<=0 )	niter=1;
 else			niter++;

 while ( niter>0 )
  {	fit_kernels_native(ref,img,mask_conv,klist,xlist,method,bdc);
	niter--;
	if ( niter>0 )
	 {	image	sub;
		int		i,j,rej,tot;
		char		**ms;
		double		s,s2,w,n,sig,nz;

		ms=mask_conv;

		/* image_duplicate(&sub,ref,1) → contiguous alloc + memcpy */
    sub.sx = ref->sx;
    sub.sy = ref->sy;
    {	size_t _rb = sub.sx * sizeof(double);
        double *_fb = malloc(sub.sx * sub.sy * sizeof(double));
        sub.data = malloc(sub.sy * sizeof(double *));
        int _ii;
        for (_ii=0 ; _ii<sub.sy ; _ii++) {
            sub.data[_ii] = _fb + _ii * sub.sx;
            memcpy(sub.data[_ii], ref->data[_ii], _rb);
        }
    }
		convolve_to_subtracted(ref,img,ms,klist,xlist,&sub);
		n=s=s2=0.0;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	if ( ! ms[i][j] )	continue;
				w=sub.data[i][j];
				s+=w,s2+=w*w,n+=1.0;
			 }
		 }
		s/=n,s2/=n;
		sig=s2-s*s;
				rej=tot=0;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ;j<sx ;  j++ )
			 {	if ( ! ms[i][j] )	continue;
				tot++;
				nz=sqrt(sig+ref->data[i][j]/gain);
				if ( fabs(sub.data[i][j])>rlevel*nz )
					ms[i][j]=0,rej++;
			 }
		 }
		logmsg(is_verbose>=1,"Rejected: %d from %d\n",rej,tot);
		    /* fits_image_free → free flat + row pointers */
    free(sub.data[0]);
    free(sub.data);
	 }
  };

 tensor_free(mask_conv);
 mask_free(mask);

 return(0);
}

/*****************************************************************************/

int stamp_parse_argument(char *stamparg,char **fitmask,int sx,int sy)
{
 char	*tmpstamparg,**cmd;
 int	k,ix,iy,is,i,j;
 double	x,y,s;

 tmpstamparg=strdup(stamparg);

 cmd=tokenize_char_dyn(tmpstamparg,':');
 for ( k=0 ; cmd != NULL && cmd[k] != NULL ; k++ )
  {	if ( sscanf(cmd[k],"%lg,%lg,%lg",&x,&y,&s)<3 )
 		continue;
 	if ( s<=0.0 )
 		continue;
 	ix=(int)x,
 	iy=(int)y;
 	is=(int)s;
 	for ( i=iy-is ; i<=iy+is ; i++ )
 	 {	if ( i<0 || i>=sy )	continue;
 		for ( j=ix-is ; j<=ix+is ; j++ )
 		 {	if ( j<0 || j>=sx )	continue;
 			fitmask[i][j] |= 0x80;
 		 }
 	 }
  }
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
 	 {	if ( fitmask[i][j] & 0x80 )
 			fitmask[i][j] &= 0x7f;
 		else
 			fitmask[i][j] = (fitmask[i][j] & 0x7f) | MASK_OUTER;
 	 }
  }

 if ( cmd != NULL )	free(cmd);
 free(tmpstamparg);

 return(0);
}

/*****************************************************************************/

/* stamp_read_file removed — file I/O not supported in Python pipeline */

/*****************************************************************************/

/*****************************************************************************/
/* dump functions for checkpointing (commented out)                         */
/*****************************************************************************/
/*
#define FICONV_DUMP_DIR "ficonv_dump"

static FILE *ficonv_dump_open(char *stage, char *ext)
{
    char fname[256];
    sprintf(fname, FICONV_DUMP_DIR "/ficonv_%s.%s", stage, ext);
    return fopen(fname, ext[0] == 'b' ? "wb" : "w");
}

static void ficonv_dump_image(char *stage, double **data, int sx, int sy)
{
    FILE *fb = ficonv_dump_open(stage, "bin");
    int i, j;
    if (fb != NULL) {
        fwrite(&sx, sizeof(sx), 1, fb); fwrite(&sy, sizeof(sy), 1, fb);
        for (i = 0; i < sy; i++)
            for (j = 0; j < sx; j++)
                fwrite(&data[i][j], sizeof(double), 1, fb);
        fclose(fb);
    }
}

static void ficonv_dump_mask(char *stage, char **mask, int sx, int sy)
{
    FILE *fb = ficonv_dump_open(stage, "mask.bin");
    int i, j;
    if (fb != NULL) {
        fwrite(&sx, sizeof(sx), 1, fb); fwrite(&sy, sizeof(sy), 1, fb);
        for (i = 0; i < sy; i++)
            for (j = 0; j < sx; j++)
                fwrite(&mask[i][j], sizeof(char), 1, fb);
        fclose(fb);
    }
}
*/

 /*****************************************************************************/
 /*        ficonv_fit_cy 已删除。请使用 fitsh_ficonv_fit_cy                    */
 /*****************************************************************************/
 
int fitsh_ficonv_fit_cy(
    double *ref_data, char *ref_mask,
    double *img_data, char *img_mask,
    char *inmask_data,
    int sx, int sy,
    char *kernel_spec,
    int method,
    int bdc,
    int niter, double rejlevel, double gain,
    int is_verbose_arg,
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
    double prefit_ox, double prefit_oy, double prefit_scale)
{
 int	i,r = 0;
 image	refimg, img, outimg;
 char		**mr, **mi, **mo, **inmask = NULL;
 kernellist	kl_data, xl_data, *kl = &kl_data, *xl = &xl_data;

 progbasename = "fitsh_cy.ficonv";
 is_verbose = is_verbose_arg;
 is_comment = 0;

 refimg.sx = sx; refimg.sy = sy;
 refimg.data = malloc(sy * sizeof(double *));
 for (i=0; i<sy; i++) refimg.data[i] = ref_data + i * sx;

 img.sx = sx; img.sy = sy;
 img.data = malloc(sy * sizeof(double *));
 for (i=0; i<sy; i++) img.data[i] = img_data + i * sx;

 outimg.sx = sx; outimg.sy = sy;
 outimg.data = malloc(sy * sizeof(double *));
 for (i=0; i<sy; i++) outimg.data[i] = cnv_data + i * sx;

 mr = malloc(sy * sizeof(char *));
 for (i=0; i<sy; i++) mr[i] = ref_mask + i * sx;
 mi = malloc(sy * sizeof(char *));
 for (i=0; i<sy; i++) mi[i] = img_mask + i * sx;
 mo = malloc(sy * sizeof(char *));
 for (i=0; i<sy; i++) mo[i] = cnv_mask + i * sx;
  if (inmask_data) {
      inmask = malloc(sy * sizeof(char *));
      for (i=0; i<sy; i++) inmask[i] = inmask_data + i * sx;
  }
  if (stamp_arg && stamp_arg[0]) {
      if (!inmask) {
          inmask = mask_create_empty(sx, sy);
      }
      stamp_parse_argument(stamp_arg, inmask, sx, sy);
  }

  /* pre-spline affine correction */
  if (psx > 0 && psy > 0) {
      char **spfmask;
      double *coeff;
      if (spline_stamp_arg && spline_stamp_arg[0]) {
          spfmask = mask_create_empty(sx, sy);
          stamp_parse_argument(spline_stamp_arg, spfmask, sx, sy);
      } else {
          spfmask = mask_create_empty(sx, sy);
      }
      if (inmask != NULL) mask_and(spfmask, sx, sy, inmask);
      mask_and(spfmask, sx, sy, mr);
      mask_and(spfmask, sx, sy, mi);
      coeff = malloc(sizeof(double) * 2 * (psx+1) * (psy+1));
      affine_img_transformation_spline_fit(&img, &refimg, spfmask, psx, psy, coeff);
      affine_img_transformation_spline_eval(&img, mi, psx, psy, coeff);
      free(coeff);
      mask_free(spfmask);
  }

  memset(kl, 0, sizeof(kernellist));
 memset(xl, 0, sizeof(kernellist));
  /* 两条路径互斥：预拟合数据（来自 kernel_dict）与内核规格字符串（-k 参数） */
  if (prefit_nkernels <= 0      /* 无预拟合数据                            */
      && kernel_spec             /* 内核规格字符串已提供                    */
      && create_kernels_from_kernelarg(kernel_spec, kl))  /* 解析规格构建 klist */
   {	r = 1; goto cleanup; }

 kl->type = 0;

 for (i=0 ; i<kl->nkernel ; i++)
 	kl->kernels[i].target = 0;

 xl->nkernel = 0;
 xl->kernels = NULL;

  logmsg(is_verbose>=1, "Number of kernels: %d.\n", kl->nkernel);

  /* DUMP: kernel_init_images input (commented out)
   { char p[512]; FILE *df; int ki, oo, nv;
     snprintf(p,sizeof(p),"%s/pipe_kernelinit_in.bin",
       "/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp");
     df=fopen(p,"wb"); if(df){
       fwrite(&kl->nkernel,4,1,df);
       fwrite(&kl->ox,8,1,df); fwrite(&kl->oy,8,1,df); fwrite(&kl->scale,8,1,df);
       fwrite(&kl->type,4,1,df);
       for(ki=0;ki<kl->nkernel;ki++){
         kernel *kk=&kl->kernels[ki];
         fwrite(&kk->type,4,1,df); fwrite(&kk->order,4,1,df);
         fwrite(&kk->hsize,4,1,df); fwrite(&kk->sigma,8,1,df);
         fwrite(&kk->bx,4,1,df); fwrite(&kk->by,4,1,df);
         fwrite(&kk->target,4,1,df); fwrite(&kk->flag,4,1,df);
         nv=(kk->order+1)*(kk->order+2)/2;
         fwrite(&nv,4,1,df);
         if(kk->coeff) fwrite(kk->coeff,8,nv,df);
       }
       fclose(df);
     }
   }
   */

  kernel_init_images(kl);
  kernel_init_images(xl);

  /* DUMP: kernel_init_images output (commented out)
   { char p[512]; FILE *df; int ki, oo, nv, hsz, _i, _j;
     snprintf(p,sizeof(p),"%s/pipe_kernelinit_out.bin",
       "/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp");
     df=fopen(p,"wb"); if(df){
       fwrite(&kl->nkernel,4,1,df);
       fwrite(&kl->ox,8,1,df); fwrite(&kl->oy,8,1,df); fwrite(&kl->scale,8,1,df);
       fwrite(&kl->type,4,1,df);
       for(ki=0;ki<kl->nkernel;ki++){
         kernel *kk=&kl->kernels[ki];
         fwrite(&kk->type,4,1,df); fwrite(&kk->order,4,1,df);
         fwrite(&kk->hsize,4,1,df); fwrite(&kk->sigma,8,1,df);
         fwrite(&kk->bx,4,1,df); fwrite(&kk->by,4,1,df);
         fwrite(&kk->target,4,1,df); fwrite(&kk->flag,4,1,df);
         nv=(kk->order+1)*(kk->order+2)/2;
         fwrite(&nv,4,1,df);
         if(kk->coeff) fwrite(kk->coeff,8,nv,df);
         hsz = 2*kk->hsize + 1;
         if(kk->image) { for(_i=0;_i<hsz;_i++) fwrite(kk->image[_i],8,hsz,df); }
       }
       fclose(df);
     }
   }
    return 0;  (commented out — remove to dump and exit early)
    */

  if (prefit_nkernels > 0 && prefit_coeffs != NULL)
    {	/* 从预拟合数据构建 klist：构造系数指针数组，调 kernel_info_read_dicts_in_ficonv */
        double **coeff_ptrs;
        int ki, ci;
        coeff_ptrs = (double **)malloc(prefit_nkernels * sizeof(double *));
        for (ki = 0, ci = 0; ki < prefit_nkernels; ki++) {
            coeff_ptrs[ki] = &prefit_coeffs[ci];
            ci += prefit_ncoeffs[ki];
        }
        kernel_info_read_dicts_in_ficonv(kl,
            prefit_nkernels, prefit_ox, prefit_oy, prefit_scale, prefit_ktype,
            prefit_types, prefit_orders, prefit_ncoeffs,
            prefit_hsizes, prefit_sigmas,
            prefit_bx, prefit_by,
            coeff_ptrs);
        free(coeff_ptrs);
    }
  else

  r = fit_kernels(&refimg, mr, &img, mi, inmask,
  	kl, xl, method, bdc, niter, rejlevel, gain);
  if ( r )
  {	r = 2; goto cleanup; }

   /* output kernel data via arrays */
   if ( out_nkernel != NULL )
    {	int ki, oo, nv, nk = kl->nkernel, ci = 0;
        if (nk > max_kernels) nk = max_kernels;
        *out_nkernel = nk;
        *out_ox = kl->ox;
        *out_oy = kl->oy;
        *out_scale = kl->scale;
        *out_type = kl->type;
        for (ki=0; ki<nk; ki++)
         {	kernel *kk = &kl->kernels[ki];
            out_k_types[ki] = kk->type;
            out_k_orders[ki] = kk->order;
            nv = (kk->order+1)*(kk->order+2)/2;
            out_k_ncoeffs[ki] = nv;
            if ( kk->type == KERNEL_GAUSSIAN ) {
                out_k_hsizes[ki] = kk->hsize;
                out_k_sigmas[ki] = kk->sigma;
            } else {
                out_k_hsizes[ki] = 0;
                out_k_sigmas[ki] = 0.0;
            }
            out_k_bx[ki] = kk->bx;
            out_k_by[ki] = kk->by;
            if ( kk->coeff && kl->type==1 ) {
                for (oo=0; oo<nv; oo++)
                    out_kcoeffs[ci++] = kk->coeff[oo];
            } else {
                for (oo=0; oo<nv; oo++)
                    out_kcoeffs[ci++] = 0.0;
            }
         }
        *out_ktotal = ci;
    }

 /* build proper convolution mask (same as ficonv.c main()) */
  {	int	 hsize = 0;
  	kernel	*k;
  	for (i=0,k=kl->kernels; i<kl->nkernel; i++,k++)
  	 {	if ( k->type==KERNEL_BACKGROUND || k->type==KERNEL_IDENTITY )
  			continue;
  		if (k->hsize > hsize) hsize = k->hsize;
  	 }
 	{	char **mask_base, **cmask;
 		mask_base = mask_create_empty(sx,sy);
 		mask_and(mask_base,sx,sy,mr);
 		if (inmask != NULL) mask_and(mask_base,sx,sy,inmask);
  	cmask = mask_expand_false(mask_base,sx,sy,hsize,-1,-1,1);
  	mask_free(mask_base);
  	r = convolve_with_kernel_set(&refimg, cmask, kl, &outimg);
  		for (i=0; i<sy; i++)
 		memcpy(mo[i], cmask[i], sx);
 		if (add_data != NULL)
 		 {	int k;
 		 	for (i=0; i<sy; i++)
 		 	 {	for (k=0; k<sx; k++)
 		 	 	 {	if ( ! cmask[i][k] )
 		 	 	 		outimg.data[i][k] += add_data[i*sx + k],
 		 	 	 		mo[i][k] = 0;
 		 	 	 	else
 		 	 	 		mo[i][k] = (char)(-1);
 		 	 	 }
 		 	 }
 		 }
 		if (sub_data != NULL)
 		 {	image sub_img;
 		 	char **submask_out;
 		 	sub_img.sx = sx; sub_img.sy = sy;
 		 	sub_img.data = malloc(sy * sizeof(double *));
 		 	{	double *fb = malloc(sx * sy * sizeof(double));
 		 		for (i=0; i<sy; i++) sub_img.data[i] = fb + i * sx; }
 		 	submask_out = mask_create_empty(sx,sy);
 		 	make_subtracted_image(&outimg, cmask, &img, mi, &sub_img, submask_out, xl);
 		 	for (i=0; i<sy; i++)
 		 		memcpy(sub_data + i * sx, sub_img.data[i], sx * sizeof(double));
 		 	free(sub_img.data[0]); free(sub_img.data);
 		 	mask_free(submask_out);
 		 }
 		mask_free(cmask);
 	}
 }
 if ( r )
  {	r = 3; goto cleanup; }

 /* unity kernels — after convolution, matching CLI */
 if ( unity_kernels )
  {	int ki; double norm = 1.0;
  	for (ki=0; ki<kl->nkernel; ki++)
  	 {	if (kl->kernels[ki].type == KERNEL_IDENTITY)
  	 	 {	int oo, nv;
  	 	 	norm = kl->kernels[ki].coeff[0];
  	 	 	nv = (kl->kernels[ki].order+1)*(kl->kernels[ki].order+2)/2;
  	 	 	kl->kernels[ki].coeff[0] = 1.0;
  	 	 	for (oo=1; oo<nv; oo++) kl->kernels[ki].coeff[oo] = 0.0;
  	 	 }
  	 }
  	for (ki=0; ki<kl->nkernel; ki++)
   	 {	if (kl->kernels[ki].type == KERNEL_DDELTA && 0 < norm)
  	 	 {	int oo, nv;
  	 	 	nv = (kl->kernels[ki].order+1)*(kl->kernels[ki].order+2)/2;
  	 	 	for (oo=0; oo<nv; oo++) kl->kernels[ki].coeff[oo] /= norm;
  	 	 }
  	 }
  }

  /* dump: final output (commented out)
   { char p[512]; FILE *df; int ii;
     snprintf(p,sizeof(p),"%s/ficonv_dump_out.bin",
       "/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp");
     df=fopen(p,"wb"); if(df){
       fwrite(&sx,4,1,df); fwrite(&sy,4,1,df);
       for(ii=0;ii<sy;ii++) fwrite(cnv_data+ii*sx,8,sx,df);
       fclose(df);
     }
   }
   */
  r = 0;

cleanup:
 free(refimg.data);
 free(img.data);
 free(outimg.data);
 free(mr); free(mi); free(mo); if (inmask) free(inmask);
 return r;
}

int kernel_info_read_dicts_in_ficonv(
    kernellist *kl,
    int nkernel,                         /* 对应传入 prefit_nkernels */
    double ox, double oy, double scale,
    int ktype,                           /* 对应传入 prefit_ktype */
    int *k_types,                        /* 对应传入 prefit_types */
    int *k_orders,                       /* 对应传入 prefit_orders */
    int *k_ncoeffs,                      /* 对应传入 prefit_ncoeffs */
    int *k_hsizes,                       /* 对应传入 prefit_hsizes */
    double *k_sigmas,                    /* 对应传入 prefit_sigmas */
    int *k_bx,                           /* 对应传入 prefit_bx */
    int *k_by,                           /* 对应传入 prefit_by */
    double **k_coeffs)                   /* 由 prefit_coeffs 按 ncoeffs 偏移构造的指针数组 */
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
