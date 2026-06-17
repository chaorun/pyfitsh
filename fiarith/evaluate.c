/*****************************************************************************/
/* evaluate.c — fiarith 表达式求值器                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* 从 origincode fiarith.c 提取并适配：fits* → image*，去文件 I/O。        */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <psn/psn.h>
#include "image.h"
#include "tensor.h"
#include "mask.h"

#include "math/spline/biquad.h"
#include "math/dft/pbfft.h"
#include "psn/psn-general.h"

#include "evaluate.h"

/*****************************************************************************/

#define		IO_ADD		1
#define		IO_SUB		2
#define		IO_MUL		3
#define		IO_DIV		4
#define		IO_CHS		6

#define		IF_MIN		32
#define		IF_MAX		33
#define		IF_MEA		34

#define		IF_ABS		36
#define		IF_SQRT		37
#define		IF_SQ		38
#define		IF_NORM		39
#define		IF_SIGN		40
#define		IF_THETA	41

#define		IF_LAP		56
#define		IF_SCT		57
#define		IF_SMT		58
#define		IF_WGH		59
#define		IF_CORR		60

#define		UDF_OFFSET	64

/*****************************************************************************/

psnsym		psn_subimg_var[] = {	{ T_VAR, 0, "x", 0 },
					{ T_VAR, 1, "y", 0 },
					{ T_VAR, 2, "X", 0 },
					{ T_VAR, 3, "Y", 0 },
					{ T_VAR, 4, "a", 0 },
					{ T_VAR, 5, "b", 0 },
					{ T_VAR, 6, "c", 0 },
					{ T_VAR, 7, "d", 0 },
					{ T_VAR, 8, "e", 0 },
					{ T_VAR, 9, "f", 0 },
					{ T_VAR,10, "g", 0 },
					{ T_VAR,11, "h", 0 },
					{ 0,0,NULL,0 }
				   };

psnsym		psn_img_op[]	= {	{ T_OP, IO_ADD, "+", TO_INFIX  },
					{ T_OP, 0     , "+", TO_PREFIX },
					{ T_OP, IO_SUB, "-", TO_INFIX  },
					{ T_OP, IO_CHS, "-", TO_PREFIX },
					{ T_OP, IO_MUL, "*", TO_INFIX  },
					{ T_OP, IO_DIV, "/", TO_INFIX  },

					{0,0,NULL,0}
				  };

psnsym		psn_img_fn[]	= {	
					{ T_FN, IF_MIN, "min"    , -1  },
					{ T_FN, IF_MAX, "max"    , -1  },
					{ T_FN, IF_MEA, "mean"   , -1  },
					{ T_FN, IF_ABS, "abs"	 ,  1  },
					{ T_FN, IF_SIGN,"sign"	 ,  1  },
					{ T_FN, IF_THETA,"theta" ,  1  },
					{ T_FN, IF_NORM, "norm"	 ,  1  },
					{ T_FN, IF_SQRT,"sqrt"	 ,  1  },
{ T_FN, IF_SQ,  "sq"    ,  1  },
					{ T_FN, IF_LAP, "laplace",  1  },
					{ T_FN, IF_SCT, "scatter",  1  },
					{ T_FN, IF_SMT, "smooth" ,  3  },
					{ T_FN, IF_WGH, "weight" ,  3  },
					{ T_FN, IF_CORR,"corr"	 ,  2  },
					{ 0,0,NULL,0 } 
				 };


psnprop		psn_img_prop[] = {	{ IO_ADD,  2, 20, ASSOC_LEFT },
					{ IO_SUB,  2, 20, ASSOC_LEFT },
					{ IO_MUL,  2, 21, ASSOC_LEFT },
					{ IO_DIV,  2, 21, ASSOC_LEFT },
					{ IO_CHS,  1, 22 },
					{ IF_MIN, -1, 0 },
					{ IF_MAX, -1, 0 },
					{ IF_MEA, -1, 0 },	
					{ IF_ABS,  1, 0 },
					{ IF_SIGN, 1, 0 },
					{ IF_THETA,1, 0 },
					{ IF_NORM, 1, 0 },
					{ IF_SQRT, 1, 0 },
					{ IF_SQ,   1, 0 },
					{ IF_LAP,  1, 0 },
					{ IF_SCT,  1, 0 },
					{ IF_SMT,  3, 0 },
					{ IF_WGH,  3, 0 },
					{ IF_CORR, 2, 0 },
					{ 0, 0, 0,0  } };

/*****************************************************************************/
/* 图像内存分配（单次 malloc，row 指针 + 像素数据合一）                    */

image * image_new_empty(int sx, int sy)
{
 image  *img;
 int    i;

 if ( sx <= 0 || sy <= 0 )       return NULL;

 img = (image *)malloc(sizeof(image));
 if ( img == NULL )              return NULL;
 memset(img, 0, sizeof(image));

 img->sx = sx;
 img->sy = sy;
 img->data = (double **)malloc(sizeof(double *) * sy + sizeof(double) * sx * sy);
 if ( img->data == NULL )
   {    free(img);
        return NULL;
   }

 double *row0 = (double *)(img->data + sy);
 for ( i = 0 ; i < sy ; i++ )
        img->data[i] = row0 + i * sx;

 return img;
}

image * image_new_constant(int sx, int sy, double d)
{
 image  *img;
 int    i, j;

 img = image_new_empty(sx, sy);
 if ( img == NULL )      return NULL;

 for ( i = 0 ; i < sy ; i++ )
   {    for ( j = 0 ; j < sx ; j++ )
                img->data[i][j] = d;
   }

 return img;
}

/* 与 image_new_empty 的单次分配配套的释放函数，不逐行 free */
static void img_free(image *img)
{
 if ( img != NULL )
  {    if ( img->data != NULL )   free(img->data);
       free(img);
  }
}

/*****************************************************************************/
/* scatter: 双二次噪声水平估计                                              */

int evaluate_scatter(image *img)
{ 
 double	**bqc;
 int	i,j,sx,sy;

 sx=img->sx,
 sy=img->sy;

 bqc=(double **)tensor_alloc_2d(double,2*sx+1,2*sy+1);
 biquad_coeff(img->data,sx,sy,bqc,NULL);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	img->data[i][j]=biquad_scatter(bqc,j,i);		}
  }
 tensor_free(bqc);
 return(0);
}

/*****************************************************************************/
/* cyclic_laplace: 循环拉普拉斯变换（从 imgtrans.c 移入）                  */

static int cyclic_laplace_of_image_ign_flag(image *img, int flag)
{
 int     i, j, sx, sy;
 double  *d0, *d1, *di, *dd;

 if ( img == NULL || img->data == NULL ) return 1;

 sx = img->sx; sy = img->sy;

 dd = (double *)malloc(sizeof(double) * sx * 3);
 if ( dd == NULL ) return -1;

 di = dd; d0 = dd + sx; d1 = d0 + sx;
 for ( j = 0 ; j < sx ; j++ )
   {    di[j] = d1[j] = img->data[0][j];
        d0[j] = img->data[sy - 1][j];
   }

 for ( i = 0 ; i < sy ; i++ )
   {    double  *dprev, *dcurr, *dnext;
        dprev = d0; dcurr = d1;
        if ( i < sy - 1 )       dnext = img->data[i + 1];
        else                    dnext = di;

        img->data[i][0] = 4 * dcurr[0] - (dprev[0] + dnext[0] + dcurr[sx - 1] + dcurr[1]);
        for ( j = 1 ; j < sx - 1 ; j++ )
                img->data[i][j] = 4 * dcurr[j] - (dprev[j] + dnext[j] + dcurr[j - 1] + dcurr[j + 1]);
        img->data[i][sx - 1] = 4 * dcurr[sx - 1] - (dprev[sx - 1] + dnext[sx - 1] + dcurr[sx - 2] + dcurr[0]);

        if ( i < sy - 1 )
         {    for ( j = 0 ; j < sx ; j++ )
               {    d0[j] = d1[j]; d1[j] = img->data[i + 1][j]; }
         }
   }

 free(dd);
 return 0;
}

/*****************************************************************************/

int evaluate_min(imgstack *stack,int narg,int sx,int sy)
{
 int	i,j,k;
 double	w,c;
 image *img;
 if ( narg <= 0 )
  {	stack[0].img=NULL;stack[0].val=0.0;return(0);	  }
 for ( k=0,img=NULL ; k<narg && img==NULL ; k++ )
  {	img=stack[k].img;	}
 if ( img==NULL )
  {	w=stack[0].val;
	for ( k=1 ; k<narg ; k++ )
	 {	if ( stack[k].val<w )	w=stack[k].val;		}
	stack[0].val=w;
	return(0);
  }
 for ( i=0,w=0.0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	for ( k=0,w=0.0 ; k<narg ; k++ )
		 {	if ( stack[k].img != NULL )
				c=stack[k].img->data[i][j];
			else
				c=stack[k].val;

			if ( k==0 )	w=c;
			else if ( c<w )	w=c;
		 }
		img->data[i][j]=w;
	 } 
  }
 for ( k=0 ; k<narg ; k++ )
  {	if ( stack[k].img != NULL && stack[k].img != img )
		img_free(stack[k].img);
  }
 stack[0].img=img;stack[0].val=0.0;
 return(0);
}

int evaluate_max(imgstack *stack,int narg,int sx,int sy)
{
 int	i,j,k;
 double	w,c;
 image *img;
 if ( narg <= 0 )
  {	stack[0].img=NULL;stack[0].val=0.0;return(0);	  }
 for ( k=0,img=NULL ; k<narg && img==NULL ; k++ )
  {	img=stack[k].img;	}
 if ( img==NULL )
  {	w=stack[0].val;
	for ( k=1 ; k<narg ; k++ )
	 {	if ( stack[k].val>w )	w=stack[k].val;		}
	stack[0].val=w;
	return(0);
  }
 for ( i=0,w=0.0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	for ( k=0,w=0.0 ; k<narg ; k++ )
		 {	if ( stack[k].img != NULL )
				c=stack[k].img->data[i][j];
			else
				c=stack[k].val;
			if ( k==0 )	w=c;
			else if ( c>w )	w=c;
		 }
		img->data[i][j]=w;
	 } 
  }
 for ( k=0 ; k<narg ; k++ )
 {	if ( stack[k].img != NULL && stack[k].img != img )
		img_free(stack[k].img);
  }
 stack[0].img=img;stack[0].val=0.0;
 return(0);
}

int evaluate_mean(imgstack *stack,int narg,int sx,int sy)
{
 int	i,j,k;
 double	w,c,div;
 image *img;
 if ( narg <= 0 )
  {	stack[0].img=NULL;stack[0].val=0.0;return(0);	  }
 div=1.0/(double)narg;
 for ( k=0,img=NULL ; k<narg && img==NULL ; k++ )
  {	img=stack[k].img;	}
 if ( img==NULL )
  {	w=0.0;
	for ( k=0 ; k<narg ; k++ )
	 {	w+=stack[k].val;		}
	stack[0].val=w*div;
	return(0);
  }
 for ( i=0,w=0.0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	for ( k=0,w=0.0 ; k<narg ; k++ )
		 {	if ( stack[k].img != NULL )
				c=stack[k].img->data[i][j];
			else
				c=stack[k].val;
			w+=c;
		 }
		img->data[i][j]=w*div;
	 } 
  }
 for ( k=0 ; k<narg ; k++ )
  {	if ( stack[k].img != NULL && stack[k].img != img )
		img_free(stack[k].img);
  }
 stack[0].img=img;stack[0].val=0.0;
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int evaluate_smooth(evaldata *ed,imgstack *stack,double sigma,double dhsize)
{
 int	hsize,fsize,i,j,sx,sy,k,l;
 double	**kernel,sum,w,km,x0,x1,y0,y1;
 char	**nmask;
 image *img,*nimg;
 if ( stack->img==NULL )	return(0);	/* nothing to be done */
 hsize=(int)dhsize;
 if ( hsize<=0 )		return(0);

 fsize=2*hsize+1;
 kernel=tensor_alloc_2d(double,fsize,fsize);

 km=1.0/(sqrt(2.0)*sigma);
 for ( i=-hsize ; i<=hsize ; i++ )
  {	y0=km*((double)i-0.5);
	y1=km*((double)i+0.5);
	for ( j=-hsize ; j<=hsize ; j++ )
	 {	x0=km*((double)j-0.5);
		x1=km*((double)j+0.5);
		kernel[i+hsize][j+hsize]=(erf(x1)-erf(x0))*(erf(y1)-erf(y0));
	 }
  }

 sum=0.0;
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
	 {	sum+=kernel[i][j];		}
  }
 for ( i=0 ; i<fsize ; i++ )
  {	for ( j=0 ; j<fsize ; j++ )
	 {	kernel[i][j]/=sum;		}
  }

 img=stack->img;
 sx=ed->sx,
 sy=ed->sy;
 nmask=mask_expand_false(ed->mask,sx,sy,hsize,-1,-1,1);
 nimg=image_new_constant(ed->sx,ed->sy,0.0);
 for ( i=0 ; i<sy ; i++ )
  {	for ( j=0 ; j<sx ; j++ )
	 {	if ( nmask[i][j] )
		 {	if ( ! ed->mask[i][j] )	nimg->data[i][j]=0.0;
			continue;
		 }
		w=0.0;
		for ( k=-hsize ; k<=hsize ; k++ )
		 {	for ( l=-hsize ; l<=hsize ; l++ )
			 {	w+=img->data[i-k][j-l]*
					kernel[hsize+k][hsize+l];
			 }
		 }
		nimg->data[i][j]=w;
	 }
  }
 mask_free(ed->mask);
 ed->mask=nmask;
 img_free(img);
 stack->img=nimg;

 tensor_free(kernel);
 return(0);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int evaluate_weight(evaldata *ed,imgstack *simg,imgstack *swgh,double dbsize)
{
 int	bsize,sx,sy,i,j,bx,by,k,l;
 image *img,*wgh;
 double	sm,ii,iw,sw;

 /* nothing to be done if any of the images (target image and/or weight)   */
 /* are constant, so (imgstack *)->img == NULL (no associated FITS image): */
 if ( simg->img==NULL || swgh->img==NULL )	return(0);

 sx=ed->sx,
 sy=ed->sy;
 img=simg->img;
 wgh=swgh->img;
 
 bsize=(int)dbsize;

 /* if block size is less than or equal to 1, do nothing. */
 if ( bsize<=1 )	return(0);

 for ( i=0 ; i<sy ; i+=bsize )
  {	by=bsize;
	if ( i+by>sy )	by=sy-i;
	for ( j=0 ; j<sx ; j+=bsize )
	 {	bx=bsize;
		if ( j+bx>sx )	bx=sx-j;
		sm=sw=0.0;
		for ( k=0 ; k<by ; k++ )
		 {	for ( l=0 ; l<bx ; l++ )
			 {	ii=img->data[i+k][j+l];
				iw=wgh->data[i+k][j+l];
				sm+=(ii-iw);
				sw+=iw;
			 }
		 }
		if ( sw<=0.0 )	continue;
		sm=sm/(double)(bx*by);
		for ( k=0 ; k<by ; k++ )
		 {	for ( l=0 ; l<bx ; l++ )
			 {	iw=wgh->data[i+k][j+l];
				img->data[i+k][j+l] = iw+sm;
			 }
		 }
	 }
  }
	
 return(0); 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int evaluate_correlation(evaldata *ed,imgstack *s1,imgstack *s2)
{
 int		sx,sy,i,j;
 double		val;

 if ( s1->img==NULL && s2->img==NULL )
  {	s1->val = s1->val * s2->val;
	return(0);
  }
 else if ( s1->img==NULL || s2->img==NULL )
  {	if ( s1->img != NULL )
	 {	s2->img=s1->img;
		s1->val=s2->val;
		s1->img=NULL;
	 }
	sx=ed->sx,
	sy=ed->sy;
	val=0.0;
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	val+=s2->img->data[i][j];		}
	 }
	s1->val = s1->val * val;
	return(0);
  }
 else
  {	complex		**c1,**c2,*c,w;
	image *i1,*i2;
	double		norm;

	i1= s1->img;
	i2= s2->img;

	sx=ed->sx;
	sy=ed->sy;

	c1=tensor_alloc_2d(complex,sx,sy);
	c2=tensor_alloc_2d(complex,sx,sy);
	c =tensor_alloc_1d(complex,sy);

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
	 	 {	c1[i][j].re=i1->data[i][j];
			c1[i][j].im=0.0;
			c2[i][j].re=i2->data[i][j];
			c2[i][j].im=0.0;
		 }
		pbfft_conv(c1[i],sx,0);
		pbfft_conv(c2[i],sx,0);
	 }

	for ( j=0 ; j<sx ; j++ )
	 {	for ( i=0 ; i<sy ; i++ )
		 {	c[i].re=c1[i][j].re;
			c[i].im=c1[i][j].im;
	  	 }
		pbfft_conv(c,sy,0);
		for ( i=0 ; i<sy ; i++ )
		 {	c1[i][j].re=c[i].re;
			c1[i][j].im=c[i].im;
	  	 }

		for ( i=0 ; i<sy ; i++ )
		 {	c[i].re=c2[i][j].re;
			c[i].im=c2[i][j].im;
	  	 }
		pbfft_conv(c,sy,0);
		for ( i=0 ; i<sy ; i++ )
		 {	c2[i][j].re=c[i].re;
			c2[i][j].im=c[i].im;
	  	 }
	 }
	
	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	w.re=+c1[i][j].re*c2[i][j].re+c1[i][j].im*c2[i][j].im;
			w.im=+c1[i][j].re*c2[i][j].im-c1[i][j].im*c2[i][j].re;
			c1[i][j].re=w.re;
			c1[i][j].im=w.im;
		 }
	 }
	/*
	c1[0][0].re=0.0;
	c1[0][0].im=0.0;
	*/

	for ( j=0 ; j<sx ; j++ )
	 {	for ( i=0 ; i<sy ; i++ )
		 {	c[i].re=c1[i][j].re;
			c[i].im=c1[i][j].im;
	  	 }
		pbfft_conv(c,sy,1);
		for ( i=0 ; i<sy ; i++ )
		 {	c1[i][j].re=c[i].re;
			c1[i][j].im=c[i].im;
	  	 }
	 }

	for ( i=0 ; i<sy ; i++ )
	 {	pbfft_conv(c1[i],sx,1);		}

	norm=1.0/((double)sx*(double)sy);

	for ( i=0 ; i<sy ; i++ )
	 {	for ( j=0 ; j<sx ; j++ )
		 {	w.re=c1[i][j].re;
			w.im=c1[i][j].im;
			i1->data[i][j]=norm*w.re;
		 }
	 }
	
	tensor_free(c);
	tensor_free(c2);
	tensor_free(c1);

	return(0);
  }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

image *evaluate(evaldata *ed,psn *pseq)
{
 imgstack stack[32]; 
 image *img;
 psnterm  *seq;
 int	  i,j,k,sp,sx,sy,udfid,narg,iac;
 psn	  *udf;
 double	  uvars[16],res,crc,w;

 if ( pseq==NULL || pseq->terms==NULL ) return NULL;

 /* 确保 mask 不为 NULL（smooth/weight 需要），与原始 CLI 保持一致 */
 if ( ed->mask == NULL )
        ed->mask = mask_create_empty(ed->sx,ed->sy);

 sx=ed->sx,
 sy=ed->sy;

 sp=0;
 for ( seq=pseq->terms ; seq->type ; seq++ )
  { narg=seq->minor;
    switch( seq->type )
    {  
     case T_CONST:
	stack[sp].img=NULL;
	stack[sp].val=pseq->cons[seq->major];
	sp++;
  	break;
     case T_SCONST:
	stack[sp].img=NULL;
	stack[sp].val=(double)(seq->major);
	sp++;
	break;
     case T_VAR:
	k=seq->major;
	img=ed->ops[k].img;
	stack[sp].img=img;
	stack[sp].val=0.0;
	sp++;
	break;
     case T_OP: case T_FN:
	switch ( seq->major )
	 {   case IO_ADD:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->data[i][j]+=
					stack[sp-1].img->data[i][j];
				 }
			 }
			img_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-1].img->data[i][j]+=stack[sp-2].val;	}
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->data[i][j]+=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val += stack[sp-1].val;	}
		sp--;
		break;
	     case IO_SUB:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->data[i][j]-=
					stack[sp-1].img->data[i][j];
				 }
			 }
			img_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-1].img->data[i][j]=stack[sp-2].val-stack[sp-1].img->data[i][j];	}
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->data[i][j]-=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val -= stack[sp-1].val;	}
		sp--;
		break;
	     case IO_MUL:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->data[i][j]*=
					stack[sp-1].img->data[i][j];
				 }
			 }
			img_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-1].img->data[i][j]*=stack[sp-2].val;	}
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->data[i][j]*=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val *= stack[sp-1].val;	}
		sp--;
		break;
	     case IO_DIV:
		if ( stack[sp-2].img!=NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	w=stack[sp-1].img->data[i][j];
					if ( w != 0.0 )
						stack[sp-2].img->data[i][j]/=w;
					else
						stack[sp-2].img->data[i][j]=0.0;
				 }
			 }
			img_free(stack[sp-1].img);
		 }
		else if ( stack[sp-2].img==NULL && stack[sp-1].img!=NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	w=stack[sp-1].img->data[i][j];
					if ( w != 0.0 )
						stack[sp-1].img->data[i][j]=stack[sp-2].val/w;
					else
						stack[sp-1].img->data[i][j]=0.0;
				 }
			 }
			stack[sp-2].img=stack[sp-1].img;
		 }
		else if ( stack[sp-2].img!=NULL && stack[sp-1].img==NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	stack[sp-2].img->data[i][j]/=stack[sp-1].val;	}
			 }
		 }
		else
		 {	stack[sp-2].val /= stack[sp-1].val;	}
		sp--;
		break;
	     case IO_CHS:
		if ( stack[sp-1].img != NULL )
	 	 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
					stack[sp-1].img->data[i][j]=-stack[sp-1].img->data[i][j];
			 }
		 }
		else
		 {	stack[sp-1].val=-stack[sp-1].val;	}
		break;
	     case IF_MIN:
		evaluate_min(&stack[sp-narg],narg,sx,sy);
		sp+=1-narg;
		break;
	     case IF_MAX:
		evaluate_max(&stack[sp-narg],narg,sx,sy);
		sp+=1-narg;
		break;
	     case IF_MEA:
		evaluate_mean(&stack[sp-narg],narg,sx,sy);
		break;
	     case IF_ABS:
		if ( stack[sp-1].img != NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
					stack[sp-1].img->data[i][j]=fabs(stack[sp-1].img->data[i][j]);
			 }
		 }
		else
			stack[sp-1].val=fabs(stack[sp-1].val);
		break;
	     case IF_SIGN:
		if ( stack[sp-1].img != NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	double	d;
					d=stack[sp-1].img->data[i][j];
					stack[sp-1].img->data[i][j]=(d<0?-1.0:d>0?+1.0:0.0);
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			stack[sp-1].val=(d<0?-1.0:d>0?+1.0:0.0);
		 }
		break;
	     case IF_THETA:
		if ( stack[sp-1].img != NULL )
		 {	for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	double	d;
					d=stack[sp-1].img->data[i][j];
					stack[sp-1].img->data[i][j]=(d<0?-1.0:+1.0);
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			stack[sp-1].val=(d<0?-1.0:+1.0);
		 }
		break;
	     case IF_NORM:
		if ( stack[sp-1].img != NULL )
		 {	double	norm,n;
			norm=0.0,n=0.0;
			for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	norm += fabs(stack[sp-1].img->data[i][j]);
					n    += 1.0;
				 }
			 }
			img_free(stack[sp-1].img);
			stack[sp-1].img=NULL;
			stack[sp-1].val=norm/n;
		 }
		else
			stack[sp-1].val=fabs(stack[sp-1].val);
		break;
	     case IF_SQRT:
		if ( stack[sp-1].img != NULL )
		 {	double	d;
			for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	d=stack[sp-1].img->data[i][j];
					if ( d>0.0 )	d=sqrt(d);
					else		d=0.0;
					stack[sp-1].img->data[i][j]=d;
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			if ( d>0.0 )	d=sqrt(d);
			else		d=0.0;
			stack[sp-1].val=d;
		 }
		break;
	     case IF_SQ:
		if ( stack[sp-1].img != NULL )
		 {	double	d;
			for ( i=0 ; i<sy ; i++ )
			 {	for ( j=0 ; j<sx ; j++ )
				 {	d=stack[sp-1].img->data[i][j];
					stack[sp-1].img->data[i][j]=d*d;
				 }
			 }
		 }
		else
		 {	double	d;
			d=stack[sp-1].val;
			stack[sp-1].val=d*d;
		 }
		break;
	     case IF_LAP:
		if ( stack[sp-1].img != NULL )
			cyclic_laplace_of_image_ign_flag(stack[sp-1].img,0);
		else
			stack[sp-1].val=0.0;
		break;
	     case IF_SCT:
		if ( stack[sp-1].img != NULL )
			evaluate_scatter(stack[sp-1].img);
		else
			stack[sp-1].val=0.0;
		break;
	     case IF_SMT:
		evaluate_smooth(ed,&stack[sp-3],stack[sp-2].val,stack[sp-1].val);
		if ( stack[sp-2].img != NULL )	img_free(stack[sp-2].img);
		if ( stack[sp-1].img != NULL )	img_free(stack[sp-1].img);
		sp+=1-3;
		break;
	     case IF_WGH:
		evaluate_weight(ed,&stack[sp-3],&stack[sp-2],stack[sp-1].val);
		if ( stack[sp-2].img != NULL )	img_free(stack[sp-2].img);
		if ( stack[sp-1].img != NULL )	img_free(stack[sp-1].img);
		sp+=1-3;
		break;
	     case IF_CORR:
		evaluate_correlation(ed,&stack[sp-2],&stack[sp-1]);
		if ( stack[sp-1].img != NULL )	img_free(stack[sp-1].img);
		sp+=1-2;
		break;
	     default :
		if ( seq->major < UDF_OFFSET )	break;
		udfid=seq->major-UDF_OFFSET;			
		udf=ed->udfs[udfid];
		narg=seq->minor;
		img=image_new_constant(ed->sx,ed->sy,0.0);
		iac=-1;crc=0.0;
		for ( i=0 ; i<sy ; i++ )
		 {	for ( j=0 ; j<sx ; j++ )
			 {	uvars[0]=(double)(2*j-sx)/(double)sx;
				uvars[1]=(double)(2*i-sy)/(double)sx;
				uvars[2]=(double)j;
				uvars[3]=(double)i;
				for ( k=0 ; k<narg ; k++ )
				 {	if ( stack[sp-narg+k].img != NULL )
						uvars[4+k]=stack[sp-narg+k].img->data[i][j];
					else
						uvars[4+k]=stack[sp-narg+k].val;
				 }
				k=psn_double_calc(udf,psn_general_funct,&res,uvars);
				if ( k )	res=0.0;
				if ( iac<0 )	crc=res,iac=0;
				else if ( ! iac && crc != res )	iac=1;
				img->data[i][j]=res;
			 }
		 }
		for ( k=0 ; k<narg ; k++ )
		 {	if ( stack[sp-narg+k].img != NULL )
				img_free(stack[sp-narg+k].img);
		 }
		sp-=narg;
		if ( ! iac )
		 {	img_free(img);
			stack[sp].img=NULL;
			stack[sp].val=crc;
		 }
		else
		 {	stack[sp].img=img;
			stack[sp].val=0.0;
		 }
		sp++;
	 }
	break;
     }
  }
 if ( sp != 1 )				return(NULL);
 else if ( stack[0].img != NULL )	return(stack[0].img);
 else					return(image_new_constant(ed->sx,ed->sy,stack[0].val));
}
