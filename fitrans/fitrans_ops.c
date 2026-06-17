/*****************************************************************************/
/* _fitrans_ops.c — fitrans image operations: noise, zoom, shrink, smooth.  */
/* Extracted from fitrans.c.                                                 */
/*****************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "math/spline/spline.h"
#include "fbase.h"
#include "statistics.h"
#include "tensor.h"
#include "fitrans_ops.h"

#define FOP_MASK_OUTER 0x01
#define FOP_MASK_FAULT 0x04

typedef struct { int sx,sy; double **data; } fop_img;

/* ---- helpers from fitrans.c ---- */
static int fop_min4(int a,int b,int c,int d){
    if(b<a)a=b;if(c<a)a=c;if(d<a)a=d;return a;}
static int fop_max4(int a,int b,int c,int d){
    if(b>a)a=b;if(c>a)a=c;if(d>a)a=d;return a;}

/* ---- noise (scatter_image) ---- */
int fitrans_noise_cy(
    double *in_data,char *in_mask,int sx,int sy,double *out_data,char *out_mask)
{
    int i,j;
    fop_img img,out; char **mk,**om;
    double **bqc;
    img.sx=out.sx=sx;img.sy=out.sy=sy;
    img.data=(double**)malloc(sy*sizeof(double*));
    out.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++){img.data[i]=in_data+i*sx;out.data[i]=out_data+i*sx;}
    mk=(char**)malloc(sy*sizeof(char*));om=(char**)malloc(sy*sizeof(char*));
    for(i=0;i<sy;i++){mk[i]=in_mask+i*sx;om[i]=out_mask+i*sx;}
    bqc=(double**)tensor_alloc_2d(double,2*sx+1,2*sy+1);
    biquad_coeff(img.data,sx,sy,bqc,mk);
    for(i=0;i<sy;i++)for(j=0;j<sx;j++){
        if(!mk[i][j]){om[i][j]=0;out.data[i][j]=biquad_scatter(bqc,j,i);}
        else{om[i][j]=mk[i][j];out.data[i][j]=0;}}
    tensor_free(bqc);free(img.data);free(out.data);free(mk);free(om);
    return 0;
}

/* ---- zoom (zoom_image / zoom_raw_image) ---- */
int fitrans_zoom_cy(
    double *in_data,char *in_mask,int sx,int sy,
    double *out_data,char *out_mask,int nsx,int nsy,
    int ofx,int ofy,int scalex,int scaley,int is_raw)
{
    int i,j,k,l,ii,jj,m;
    double **bqc,**wret,sfactor;
    fop_img img,out; char **mk,**om;
    img.sx=sx;img.sy=sy;img.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++)img.data[i]=in_data+i*sx;
    out.sx=nsx;out.sy=nsy;out.data=(double**)malloc(nsy*sizeof(double*));
    for(i=0;i<nsy;i++)out.data[i]=out_data+i*nsx;
    mk=(char**)malloc(sy*sizeof(char*));om=(char**)malloc(nsy*sizeof(char*));
    for(i=0;i<sy;i++)mk[i]=in_mask+i*sx;
    for(i=0;i<nsy;i++)om[i]=out_mask+i*nsx;
    if(!scalex||!scaley||!nsx||!nsy||!sx||!sy){free(img.data);free(out.data);free(mk);free(om);return 0;}
    if(is_raw){sfactor=1.0/(double)(scalex*scaley);
        for(i=0;i<nsy/scaley;i++){ii=i+ofy;
            for(j=0;j<nsx/scalex;j++){jj=j+ofx;m=0;
                if(ii<0||jj<0||ii>sy-1||jj>sx-1){m|=FOP_MASK_OUTER;
                    for(k=0;k<scaley;k++)for(l=0;l<scalex;l++){out.data[i*scaley+k][j*scalex+l]=0;om[i*scaley+k][j*scalex+l]=m;}}
                else{m=mk[ii][jj];
                    for(k=0;k<scaley;k++)for(l=0;l<scalex;l++){out.data[i*scaley+k][j*scalex+l]=m?0:img.data[ii][jj]*sfactor;om[i*scaley+k][j*scalex+l]=m;}}}}}
    else{bqc=(double**)tensor_alloc_2d(double,2*sx+1,2*sy+1);biquad_coeff(img.data,sx,sy,bqc,mk);
        wret=(double**)tensor_alloc_2d(double,scalex,scaley);
        for(i=0;i<nsy/scaley;i++){ii=i+ofy;
            for(j=0;j<nsx/scalex;j++){jj=j+ofx;m=0;
                if(ii<0||jj<0||ii>sy-1||jj>sx-1){m|=FOP_MASK_OUTER;
                    for(k=0;k<scaley;k++)for(l=0;l<scalex;l++)wret[k][l]=0;}
                else if(mk[ii][jj]){m|=mk[ii][jj];
                    for(k=0;k<scaley;k++)for(l=0;l<scalex;l++)wret[k][l]=0;}
                else{biquad_isc_int_block_subpixels(bqc,jj,ii,scalex,scaley,wret);m=0;}
                for(k=0;k<scaley;k++)for(l=0;l<scalex;l++){out.data[i*scaley+k][j*scalex+l]=wret[k][l];om[i*scaley+k][j*scalex+l]=m;}}}
        tensor_free(bqc);tensor_free(wret);}
    free(img.data);free(out.data);free(mk);free(om);
    return 0;
}

/* ---- shrink (shrink_image) ---- */
typedef struct {int mode_median,mode_truncated_mean,mode_average_mask;} fop_shrink;
int fitrans_shrink_cy(
    double *in_data,char *in_mask,int sx,int sy,
    double *out_data,char *out_mask,int nsx,int nsy,
    int ofx,int ofy,int scalex,int scaley,
    int mode_median,int mode_truncated_mean,int mode_average_mask)
{
    int i,j,k,l,ii,jj,m,t;
    double w,*medarr;
    fop_shrink sm;sm.mode_median=mode_median;sm.mode_truncated_mean=mode_truncated_mean;sm.mode_average_mask=mode_average_mask;
    fop_img img,out; char **mk,**om;
    img.sx=sx;img.sy=sy;img.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++)img.data[i]=in_data+i*sx;
    out.sx=nsx;out.sy=nsy;out.data=(double**)malloc(nsy*sizeof(double*));
    for(i=0;i<nsy;i++)out.data[i]=out_data+i*nsx;
    mk=(char**)malloc(sy*sizeof(char*));om=(char**)malloc(nsy*sizeof(char*));
    for(i=0;i<sy;i++)mk[i]=in_mask+i*sx;
    for(i=0;i<nsy;i++)om[i]=out_mask+i*nsx;
    if(!sx||!sy||!nsx||!nsy){free(img.data);free(out.data);free(mk);free(om);return 0;}
    medarr=(sm.mode_median||sm.mode_truncated_mean>0)?(double*)malloc(scalex*scaley*8):NULL;
    for(i=0;i<nsy;i++){ii=i*scaley+ofy;
        for(j=0;j<nsx;j++){jj=j*scalex+ofx;
            if(ii<0||jj<0||ii>sy-scaley||jj>sx-scalex){out.data[i][j]=0;om[i][j]=FOP_MASK_OUTER;}
            else if(medarr){
                m=0;
                if(!sm.mode_average_mask){
                    for(k=0,t=0;k<scaley;k++)for(l=0;l<scalex;l++){m|=mk[ii+k][jj+l];if(!mk[ii+k][jj+l])medarr[t++]=img.data[ii+k][jj+l];}
                    om[i][j]=m;
                    if(sm.mode_median)out.data[i][j]=median(medarr,t)*(double)(scalex*scaley);
                    else if(sm.mode_truncated_mean>0)out.data[i][j]=truncated_mean(medarr,t,sm.mode_truncated_mean)*(double)(scalex*scaley);}
                else{for(k=0,t=0;k<scaley;k++)for(l=0;l<scalex;l++){m|=~mk[ii+k][jj+l];if(!mk[ii+k][jj+l])medarr[t++]=img.data[ii+k][jj+l];}
                    om[i][j]=(~m)&0xFF;
                    if(sm.mode_median)out.data[i][j]=median(medarr,t)*(double)(scalex*scaley);
                    else if(sm.mode_truncated_mean>0)out.data[i][j]=truncated_mean(medarr,t,sm.mode_truncated_mean)*(double)(scalex*scaley);}}
            else{w=0;m=0;
                for(k=0;k<scaley;k++)for(l=0;l<scalex;l++){m|=mk[ii+k][jj+l];w+=img.data[ii+k][jj+l];}
                om[i][j]=m;out.data[i][j]=w;}}}
    if(medarr)free(medarr);free(img.data);free(out.data);free(mk);free(om);
    return 0;
}

/* ---- smooth (combined_smooth_image + helpers) ---- */
#define FOP_SMOOTH_NONE 0
#define FOP_SMOOTH_SPLINE 1
#define FOP_SMOOTH_POLY 2
#define FOP_PREF_NONE 0
#define FOP_PREF_MEAN 1
#define FOP_PREF_MEDIAN 2

typedef struct {int x,y;double val;} fop_ppixel;
static int fop_ppix_cmp(const void *a,const void *b){
    double da=((fop_ppixel*)a)->val,db=((fop_ppixel*)b)->val;
    return da<db?-1:1;}
static int fop_ppix_sort(fop_ppixel *p,int n){if(p&&n>0)qsort(p,n,sizeof(fop_ppixel),fop_ppix_cmp);return 0;}

static int fop_prefilter(double **data,int sx,int sy,double **fltd,
    int filter,int fxh,int fyh,double frejratio,char **mask,char **fmsk)
{
    int i,j,k,l,k0,k1,hx=fxh,hy=fyh;
    fop_ppixel *pps,*ppt,*ppn,**apr,*ppw;int npp,tpp,nnp,*ran;double sumpp;
    if(!filter||(hx<=0&&hy<=0)){for(i=0;i<sy;i++)for(j=0;j<sx;j++){fltd[i][j]=data[i][j];if(fmsk&&mask)fmsk[i][j]=mask[i][j];else if(fmsk)fmsk[i][j]=0;}return 0;}
    pps=(fop_ppixel*)malloc((2*hx+1)*(2*hy+1)*sizeof(fop_ppixel));
    ppt=(fop_ppixel*)malloc((2*hx+1)*(2*hy+1)*sizeof(fop_ppixel));
    ppn=(fop_ppixel*)malloc((2*hy+1)*sizeof(fop_ppixel));
    apr=(fop_ppixel**)tensor_alloc_2d(fop_ppixel,2*hy+1,sx);
    ran=(int*)malloc(sx*sizeof(int));
    for(j=0;j<sx;j++){l=0;for(k=0;k<=hy&&k<sy;k++){if(mask&&mask[k][j])continue;apr[j][l].x=j;apr[j][l].y=k;apr[j][l].val=data[k][j];l++;}fop_ppix_sort(apr[j],l);ran[j]=l;}
    for(i=0;i<sy;i++){npp=0;sumpp=0;k0=-hy;if(i+k0<0)k0=-i;k1=hy;if(i+k1>=sy)k1=sy-1-i;
        for(k=k0;k<=k1;k++)for(l=0;l<=hx&&l<sx;l++){if(mask&&mask[i+k][l])continue;pps[npp].x=l;pps[npp].y=i+k;sumpp+=(pps[npp].val=data[i+k][l]);npp++;}
        fop_ppix_sort(pps,npp);
        for(j=0;j<sx;j++){
            if(npp>0){if(filter==FOP_PREF_MEDIAN)fltd[i][j]=0.5*(pps[(npp-1)/2].val+pps[npp/2].val);
            else if(filter==FOP_PREF_MEAN&&frejratio<=0)fltd[i][j]=sumpp/npp;
            else if(filter==FOP_PREF_MEAN){int c0=(int)(npp*(frejratio/2.0)),c1=npp-c0;fltd[i][j]=0;for(k=c0;k<c1;k++)fltd[i][j]+=pps[k].val;if(c1>c0)fltd[i][j]/=(c1-c0);else{fltd[i][j]=0;fmsk[i][j]=FOP_MASK_FAULT;}}
            else{fltd[i][j]=0;fmsk[i][j]=FOP_MASK_FAULT;}}else{fltd[i][j]=0;fmsk[i][j]=FOP_MASK_FAULT;}
            if(!(j<sx-1))continue;
            if(j+hx+1<sx){ppw=apr[j+hx+1];nnp=ran[j+hx+1];}else{ppw=NULL;nnp=0;}
            tpp=0;sumpp=0;for(k=0,l=0;k<npp;k++){if(pps[k].x>j-hx){while(l<nnp&&ppw[l].val<pps[k].val){ppt[tpp]=ppw[l];sumpp+=ppw[l].val;tpp++;l++;}ppt[tpp]=pps[k];sumpp+=pps[k].val;tpp++;}}
            while(l<nnp){ppt[tpp]=ppw[l];sumpp+=ppw[l].val;tpp++;l++;}
            ppw=pps;pps=ppt;ppt=ppw;npp=tpp;}
        if(!(i<sy-1))continue;
        for(j=0;j<sx;j++){l=0;
            if(i+hy+1<sy&&(!mask||!mask[i+hy+1][j])){ppn[0].x=j;ppn[0].y=i+hy+1;ppn[0].val=data[i+hy+1][j];nnp=1;}else nnp=0;
            tpp=0;for(k=0,l=0;k<ran[j];k++){if(apr[j][k].y>i-hy){while(l<nnp&&ppn[l].val<apr[j][k].val){ppt[tpp]=ppn[l];tpp++;l++;}ppt[tpp]=apr[j][k];tpp++;}}
            while(l<nnp){ppt[tpp]=ppn[l];tpp++;l++;}memcpy(apr[j],ppt,tpp*sizeof(fop_ppixel));ran[j]=tpp;}}
    free(ran);tensor_free(apr);free(ppn);free(ppt);free(pps);return 0;}

#include "math/fit/lmfit.h"
static int fop_smooth_img(double **data,int sx,int sy,double **fltd,int st,int xo,int yo,
    int filter,int fxh,int fyh,double fr,int nit,double lo,double up,char **mask,char **omsk)
{
    int xorder=xo,yorder=yo,nxvar=xo+1,nyvar=yo+1,nvar=nxvar*nyvar,i,j,iiter,k,l;
    double **fxbase,**fybase,*fvars,*bvector,**amatrix,w,f;
    char **tmsk;
    if(!st){xo=0;yo=0;nxvar=1;nyvar=1;nvar=1;}
    fxbase=(double**)tensor_alloc_2d(double,sx,xo+1);fybase=(double**)tensor_alloc_2d(double,sy,yo+1);
    if(st==FOP_SMOOTH_SPLINE){fbase_spline(fxbase,xo,sx);fbase_spline(fybase,yo,sy);}
    else{fbase_polynomial(fxbase,xo,sx);fbase_polynomial(fybase,yo,sy);}
    tmsk=(char**)tensor_alloc_2d(char,sx,sy);
    for(i=0;i<sy;i++){if(mask)memcpy(tmsk[i],mask[i],sx);else memset(tmsk[i],0,sx);}
    fvars=vector_alloc(nvar);bvector=vector_alloc(nvar);amatrix=matrix_alloc(nvar);
    for(iiter=0;iiter<=(nit>0?nit:0);iiter++){
        for(k=0;k<nvar;k++){for(l=0;l<nvar;l++)amatrix[k][l]=0;bvector[k]=0;}
        for(i=0;i<sy;i++)for(j=0;j<sx;j++){if(tmsk[i][j])continue;
            for(k=0;k<nyvar;k++)for(l=0;l<nxvar;l++)fvars[k*nxvar+l]=fybase[k][i]*fxbase[l][j];
            w=1;f=data[i][j];
            for(k=0;k<nvar;k++){for(l=0;l<nvar;l++)amatrix[k][l]+=w*fvars[k]*fvars[l];bvector[k]+=w*f*fvars[k];}}
        solve_gauss(amatrix,bvector,nvar);
        for(i=0;i<sy;i++)for(j=0;j<sx;j++){for(k=0;k<nyvar;k++)for(l=0;l<nxvar;l++)fvars[k*nxvar+l]=fybase[k][i]*fxbase[l][j];
            f=0;for(k=0;k<nvar;k++)f+=bvector[k]*fvars[k];fltd[i][j]=f;}
        if(iiter<nit){double s2=0,sig;for(i=0;i<sy;i++)for(j=0;j<sx;j++){f=data[i][j]-fltd[i][j];s2+=f*f;}s2/=(double)(sx*sy);
            if(s2>0)sig=sqrt(s2);else break;
            for(i=0;i<sy;i++)for(j=0;j<sx;j++){if(fabs(data[i][j]-fltd[i][j])>sig*((data[i][j]<fltd[i][j])?lo:up))tmsk[i][j]=1;}}}
    for(i=0;i<sy;i++)memset(omsk[i],0,sx);
    matrix_free(amatrix);vector_free(bvector);vector_free(fvars);tensor_free(tmsk);free(fybase);free(fxbase);
    return 0;
}

int fitrans_smooth_cy(
    double *in_data,char *in_mask,int sx,int sy,
    double *out_data,char *out_mask,
    int smooth_type,int xorder,int yorder,
    int prefilter,int fxhsize,int fyhsize,
    double frejratio,int niter,double lower,double upper,
    int is_mean_unity,int is_detrend)
{
    int i,j;double **data,**fltd,**fmsk2;char **mk,**om,**fmsk;double s,n;
    fop_img in,out;
    in.sx=sx;in.sy=sy;in.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++)in.data[i]=in_data+i*sx;
    out.sx=sx;out.sy=sy;out.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++)out.data[i]=out_data+i*sx;
    mk=(char**)malloc(sy*sizeof(char*));om=(char**)malloc(sy*sizeof(char*));
    for(i=0;i<sy;i++){mk[i]=in_mask+i*sx;om[i]=out_mask+i*sx;}
    if(prefilter){data=(double**)tensor_alloc_2d(double,sx,sy);fmsk=(char**)tensor_alloc_2d(char,sx,sy);
        fop_prefilter(in.data,sx,sy,data,prefilter,fxhsize,fyhsize,frejratio,mk,fmsk);}
    else{data=in.data;fmsk=mk;}
    if(smooth_type)fop_smooth_img(data,sx,sy,out.data,smooth_type,xorder,yorder,prefilter,fxhsize,fyhsize,frejratio,niter,lower,upper,fmsk,om);
    else for(i=0;i<sy;i++)for(j=0;j<sx;j++){out.data[i][j]=data[i][j];om[i][j]=fmsk[i][j];}
    if(is_mean_unity){s=n=0;for(i=0;i<sy;i++)for(j=0;j<sx;j++){if(om[i][j])continue;n+=1;s+=out.data[i][j];}
        if(s>0){s=n/s;for(i=0;i<sy;i++)for(j=0;j<sx;j++)out.data[i][j]*=s;}}
    else if(is_detrend){s=n=0;for(i=0;i<sy;i++)for(j=0;j<sx;j++){if(om[i][j])continue;n+=1;s+=out.data[i][j];}s/=n;
        for(i=0;i<sy;i++)for(j=0;j<sx;j++){double d=out.data[i][j];out.data[i][j]=d>0?in.data[i][j]*s/d:0;}}
    if(prefilter){tensor_free(fmsk);tensor_free(data);}
    free(in.data);free(out.data);free(mk);free(om);
    return 0;
}

/* ---- repetitive (repetitive_image from fitrans.c) ---- */
int fitrans_repetitive_cy(
    double *in_data,char *in_mask,int sx,int sy,
    double *out_data,char *out_mask,int nsx,int nsy,
    int ofx,int ofy)
{
    int i,j,ii,jj;
    fop_img img,out; char **mk,**om;
    if(!sx||!sy||!nsx||!nsy) return 0;
    img.sx=sx;img.sy=sy;img.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++)img.data[i]=in_data+i*sx;
    out.sx=nsx;out.sy=nsy;out.data=(double**)malloc(nsy*sizeof(double*));
    for(i=0;i<nsy;i++)out.data[i]=out_data+i*nsx;
    mk=(char**)malloc(sy*sizeof(char*));om=(char**)malloc(nsy*sizeof(char*));
    for(i=0;i<sy;i++)mk[i]=in_mask+i*sx;
    for(i=0;i<nsy;i++)om[i]=out_mask+i*nsx;
    if(ofx<0) ofx+=((-ofx)/nsx+2)*nsx;
    if(ofy<0) ofy+=((-ofy)/nsy+2)*nsy;
    ofx%=nsx; ofy%=nsy;
    for(i=0;i<nsy;i++){ii=(i+ofy)%sy;
        for(j=0;j<nsx;j++){jj=(j+ofx)%sx;
            out.data[i][j]=img.data[ii][jj];
            om[i][j]=mk[ii][jj];}}
    free(img.data);free(out.data);free(mk);free(om);
    return 0;
}

/* ---- interleave (interleave_image from fitrans.c) ---- */
int fitrans_interleave_cy(
    double *in_data,char *in_mask,int sx,int sy,
    double *out_data,char *out_mask,int nsx,int nsy,
    int ofx,int ofy,
    int mode_median,int mode_truncated_mean,int mode_average_mask)
{
    int i,j,k,l,ii,jj,m,t,bwx,bwy;
    double *medarr;
    fop_img img,out; char **mk,**om;
    if(!sx||!sy||!nsx||!nsy) return 0;
    bwx=(sx+nsx-1)/nsx; bwy=(sy+nsy-1)/nsy;
    img.sx=sx;img.sy=sy;img.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++)img.data[i]=in_data+i*sx;
    out.sx=nsx;out.sy=nsy;out.data=(double**)malloc(nsy*sizeof(double*));
    for(i=0;i<nsy;i++)out.data[i]=out_data+i*nsx;
    mk=(char**)malloc(sy*sizeof(char*));om=(char**)malloc(nsy*sizeof(char*));
    for(i=0;i<sy;i++)mk[i]=in_mask+i*sx;
    for(i=0;i<nsy;i++)om[i]=out_mask+i*nsx;
    if(ofx<0) ofx+=((-ofx)/nsx+2)*nsx;
    if(ofy<0) ofy+=((-ofy)/nsy+2)*nsy;
    ofx%=nsx; ofy%=nsy;
    medarr=(double*)malloc(sizeof(double)*(bwx*bwy));
    for(i=0;i<nsy;i++){int i0=(i+ofy)%nsy;
        for(j=0;j<nsx;j++){int j0=(j+ofx)%nsx;
            t=0;m=0;
            for(k=0;k<bwy;k++){ii=i0+k*nsy;if(sy<=ii)continue;
                for(l=0;l<bwx;l++){jj=j0+l*nsx;if(sx<=jj)continue;
                    medarr[t]=img.data[ii][jj];
                    if(!mode_average_mask) m|=mk[ii][jj];
                    else m|=~mk[ii][jj];
                    if(mode_average_mask&&mk[ii][jj]) continue;
                    t++;}}
            if(mode_median) out.data[i][j]=median(medarr,t);
            else out.data[i][j]=mean(medarr,t);
            if(!mode_average_mask) om[i][j]=m;
            else om[i][j]=(~m)&0xFF;}}
    if(medarr) free(medarr);
    free(img.data);free(out.data);free(mk);free(om);
    return 0;
}
