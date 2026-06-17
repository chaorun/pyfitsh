/*****************************************************************************/
/* _fitrans_core.c — fitrans image transformation: 6 interpolation modes.   */
/*****************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "math/poly.h"
#include "math/spline/bicubic.h"
#include "math/spline/biquad.h"
#include "math/spline/biquad-isc.h"
#include "tensor.h"
#include "fitrans_core.h"

typedef struct { int sx,sy; double **data; } image_lite;
typedef struct { int type,order,nval; double ox,oy,scale,**vfits; } fitrans_tf;

#define FT_IP 0
#define FT_IG 1
#define FT_SP 2
#define FT_L3 4
#define FT_L4 5
#define FT_OUTER 0x01
#define FT_EPS 1e-10

static int min4(int a,int b,int c,int d){
    if(b<a)a=b;if(c<a)a=c;if(d<a)a=d;return a;
}
static int max4(int a,int b,int c,int d){
    if(b>a)a=b;if(c>a)a=c;if(d>a)a=d;return a;
}
/*****************************************************************************/
static void jacobi(fitrans_tf *tf,double **xx,double **xy,double **yx,double **yy){
    int n=(tf->order)*(tf->order+1)/2,i,j,k;
    double s=tf->scale,*jx,*jy,*ix,*iy;
    jx=(double*)malloc(n*8);jy=(double*)malloc(n*8);
    ix=(double*)malloc(n*8);iy=(double*)malloc(n*8);
    for(i=0;i<n;i++)jx[i]=jy[i]=ix[i]=iy[i]=0;
    for(i=0,k=0;i<=tf->order-1;i++)
        for(j=0;j<=i;j++,k++){
            jx[k]+=tf->vfits[0][k+i+1]/s;jy[k]+=tf->vfits[0][k+i+2]/s;
            ix[k]+=tf->vfits[1][k+i+1]/s;iy[k]+=tf->vfits[1][k+i+2]/s;
        }
    *xx=jx;*xy=jy;*yx=ix;*yy=iy;
}
static void fwd(double x,double y,fitrans_tf *tf,double *rx,double *ry){
    *rx=eval_2d_poly(x,y,tf->order,tf->vfits[0],tf->ox,tf->oy,tf->scale);
    *ry=eval_2d_poly(x,y,tf->order,tf->vfits[1],tf->ox,tf->oy,tf->scale);
}
static void inv_tf(double x,double y,fitrans_tf *tf,double *jx,double *jy,double *ix,double *iy,double *rx,double *ry){
    int o=tf->order,i,n;double ox=tf->ox,oy=tf->oy,s=tf->scale;
    double *dx=tf->vfits[0],*dy=tf->vfits[1];
    if(o<1){*rx=x;*ry=y;return;}
    double wx=x-dx[0],wy=y-dy[0];
    double a=dx[1],b=dx[2],c=dy[1],d=dy[2];
    double dt=1.0/(a*d-b*c);
    double ia=+d*dt,ib=-b*dt,ic=-c*dt,id=+a*dt;
    double x0=ia*wx+ib*wy,y0=ic*wx+id*wy,px,py,px0,py0;
    n=100;
    for(i=0;i<n&&o>=2;i++){
        a=eval_2d_poly(x0,y0,o-1,jx,ox,oy,s);b=eval_2d_poly(x0,y0,o-1,jy,ox,oy,s);
        c=eval_2d_poly(x0,y0,o-1,ix,ox,oy,s);d=eval_2d_poly(x0,y0,o-1,iy,ox,oy,s);
        dt=1.0/(a*d-b*c);ia=+d*dt;ib=-b*dt;ic=-c*dt;id=+a*dt;
        px0=eval_2d_poly(x0,y0,o,dx,ox,oy,s);py0=eval_2d_poly(x0,y0,o,dy,ox,oy,s);
        px=x-px0;py=y-py0;x0+=ia*px+ib*py;y0+=ic*px+id*py;
        if(sqrt(px*px+py*py)<FT_EPS)break;
    }
    *rx=x0;*ry=y0;
}
/*****************************************************************************/
static int lin(image_lite *img,char **mask,double x,double y,double *ret){
    int ix=(int)floor(x),iy=(int)floor(y),flag=0;
    double u=x-ix,v=y-iy;
    if(x<0||y<0||ix>=img->sx-2||iy>=img->sy-2)return FT_OUTER;
    if(mask)flag=mask[iy][ix]|mask[iy+1][ix+1]|mask[iy][ix+1]|mask[iy+1][ix];
    *ret=img->data[iy][ix]*(1-u)*(1-v)+img->data[iy][ix+1]*u*(1-v)
        +img->data[iy+1][ix]*(1-u)*v+img->data[iy+1][ix+1]*u*v;
    return flag;
}
static int bic(double **c,int sx,int sy,char **mask,double fx,double fy,double *ret){
    int ix=(int)floor(fx),iy=(int)floor(fy),flag=0;
    if(fx<0||fy<0||ix>=sx-2||iy>=sy-2)return FT_OUTER;
    if(mask)flag=mask[iy][ix]|mask[iy+1][ix+1]|mask[iy][ix+1]|mask[iy+1][ix];
    *ret=bicubic_inter(c,fx,fy);return flag;
}
static double sinc(double x){return fabs(x)<1e-15?1.0:sin(M_PI*x)/(M_PI*x);}
static double lk(double x,int a){
    if(fabs(x)>=(double)a)return 0;return sinc(x)*sinc(x/(double)a);
}
static int lcz(image_lite *img,char **mask,double fx,double fy,int a,double *ret){
    int ix0=(int)floor(fx)-a+1,iy0=(int)floor(fy)-a+1,i,j,ix,iy,flag=0;
    double s=0,w=0,wx,wy,ww;
    if(ix0<0||iy0<0||ix0+2*a>img->sx||iy0+2*a>img->sy)return FT_OUTER;
    for(j=0;j<2*a;j++){iy=iy0+j;wy=lk(fy-iy,a);
        for(i=0;i<2*a;i++){ix=ix0+i;wx=lk(fx-ix,a);ww=wx*wy;
            s+=img->data[iy][ix]*ww;w+=ww;
            if(mask)flag|=mask[iy][ix];
    }}
    *ret=w>0?s/w:0;return flag;
}
/*****************************************************************************/
static int interpolate(
    image_lite *img,char **mask,image_lite *out,char **omask,
    int ofx,int ofy,fitrans_tf *tf,int m,int inv)
{
    int i,j,sx=img->sx,sy=img->sy,nsx=out->sx,nsy=out->sy,o=tf->order,flag;
    double x,y,nx,ny,w,ox=tf->ox,oy=tf->oy,s=tf->scale,**sc=NULL;
    double *xx,*xy,*yx,*yy,cxx,cxy,cyx,cyy,cj;
    int sp=(m==2||m==3),l3=(m==FT_L3),l4=(m==FT_L4),la=l3?3:(l4?4:0);
    if(sp){sc=(double**)tensor_alloc_2d(double,2*sx,2*sy);
        bicubic_coeff(img->data,sx,sy,sc,mask);}
    jacobi(tf,&xx,&xy,&yx,&yy);
    for(i=0;i<nsy;i++)for(j=0;j<nsx;j++){
        x=j-ofx;y=i-ofy;
        if(!inv)fwd(x,y,tf,&nx,&ny);else inv_tf(x,y,tf,xx,xy,yx,yy,&nx,&ny);
        cxx=eval_2d_poly(x,y,o-1,xx,ox,oy,s);cxy=eval_2d_poly(x,y,o-1,xy,ox,oy,s);
        cyx=eval_2d_poly(x,y,o-1,yx,ox,oy,s);cyy=eval_2d_poly(x,y,o-1,yy,ox,oy,s);
        cj=inv?1.0/(cxx*cyy-cxy*cyx):(cxx*cyy-cxy*cyx);
        if(sp)flag=bic(sc,sx,sy,mask,nx,ny,&w);
        else if(l3||l4)flag=lcz(img,mask,nx,ny,la,&w);
        else flag=lin(img,mask,nx,ny,&w);
        if(!flag){omask[i][j]=0;out->data[i][j]=w*cj;}
        else{omask[i][j]=(char)flag;out->data[i][j]=0;}
    }
    if(sc)tensor_free(sc);free(xx);free(xy);free(yx);free(yy);
    return 0;
}
/*****************************************************************************/
static int integrate(
    image_lite *img,char **mask,image_lite *out,char **omask,
    int ofx,int ofy,fitrans_tf *tf,int m,int inv)
{
    int i,j,sx=img->sx,sy=img->sy,nsx=out->sx,nsy=out->sy,flag;
    double x,y,nx,ny,dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4,w;
    double **bqc=NULL,**dl,*wd,*jx,*jy,*ix,*iy;
    int sp=(m&FT_SP),i1,i2,i3,i4,j1,j2,j3,j4,il,ih,jl,jh,iu,ju;
    if(sp){bqc=(double**)tensor_alloc_2d(double,2*sx+1,2*sy+1);
        biquad_coeff(img->data,sx,sy,bqc,mask);}
    if(inv)jacobi(tf,&jx,&jy,&ix,&iy);else jx=jy=ix=iy=NULL;
    dl=(double**)tensor_alloc_2d(double,nsx+1,4);
    y=0-ofy;for(j=0;j<=nsx;j++){x=j-ofx;
        if(!inv)fwd(x,y,tf,&nx,&ny);else inv_tf(x,y,tf,jx,jy,ix,iy,&nx,&ny);
        dl[0][j]=nx;dl[1][j]=ny;}
    for(i=0;i<nsy;i++){
        y=i+1-ofy;for(j=0;j<=nsx;j++){x=j-ofx;
            if(!inv)fwd(x,y,tf,&nx,&ny);else inv_tf(x,y,tf,jx,jy,ix,iy,&nx,&ny);
            dl[2][j]=nx;dl[3][j]=ny;}
        for(j=0;j<nsx;j++){
            dx1=dl[0][j];dx2=dl[0][j+1];dy1=dl[1][j];dy2=dl[1][j+1];
            dx3=dl[2][j];dx4=dl[2][j+1];dy3=dl[3][j];dy4=dl[3][j+1];
            i1=(int)dx1;i2=(int)dx2;i3=(int)dx3;i4=(int)dx4;
            j1=(int)dy1;j2=(int)dy2;j3=(int)dy3;j4=(int)dy4;
            flag=0;
            if(dx1<0||dy1<0||dx2<0||dy2<0)flag|=FT_OUTER;
            else if(dx3<0||dy3<0||dx4<0||dy4<0)flag|=FT_OUTER;
            else if(dx1>=sx||dy1>=sy||dx2>=sx||dy2>=sy)flag|=FT_OUTER;
            else if(dx3>=sx||dy3>=sy||dx4>=sx||dy4>=sy)flag|=FT_OUTER;
            else if(mask){il=min4(i1,i2,i3,i4);ih=max4(i1,i2,i3,i4);
                jl=min4(j1,j2,j3,j4);jh=max4(j1,j2,j3,j4);
                if(il<0){flag|=FT_OUTER;il=0;}if(ih>=sx){flag|=FT_OUTER;ih=sx-1;}
                if(jl<0){flag|=FT_OUTER;jl=0;}if(jh>=sy){flag|=FT_OUTER;jh=sy-1;}
                for(iu=jl;iu<=jh;iu++)for(ju=il;ju<=ih;ju++)flag|=mask[iu][ju];}
            if(!flag){omask[i][j]=0;
                if(sp&&bqc)w=biquad_isc_int_triangle(bqc,1,dx1,dy1,dx2,dy2,dx3,dy3,sx,sy)
                    +biquad_isc_int_triangle(bqc,1,dx2,dy2,dx4,dy4,dx3,dy3,sx,sy);
                else w=biquad_isc_int_triangle(img->data,0,dx1,dy1,dx2,dy2,dx3,dy3,sx,sy)
                    +biquad_isc_int_triangle(img->data,0,dx2,dy2,dx4,dy4,dx3,dy3,sx,sy);
                out->data[i][j]=w;}
            else{omask[i][j]=(char)flag;out->data[i][j]=0;}
        }
        wd=dl[0];dl[0]=dl[2];dl[2]=wd;wd=dl[1];dl[1]=dl[3];dl[3]=wd;
    }
    if(bqc)tensor_free(bqc);tensor_free(dl);free(jx);free(jy);free(ix);free(iy);
    return 0;
}
/*****************************************************************************/
int fitrans_apply_cy(
    double *in_data,char *in_mask,int sx,int sy,
    double *out_data,char *out_mask,int nsx,int nsy,
    int ofx,int ofy,int order,int method,int invert,
    double ox,double oy,double scale,double *dxfit,double *dyfit)
{
    int nv=(order+1)*(order+2)/2,i;
    image_lite img,out; char **mk,**om; fitrans_tf tf;
    memset(&tf,0,sizeof(tf));tf.type=1;tf.order=order;
    tf.ox=ox;tf.oy=oy;tf.scale=scale;tf.nval=2;
    tf.vfits=(double**)malloc(2*sizeof(double*));
    tf.vfits[0]=(double*)malloc(nv*8);tf.vfits[1]=(double*)malloc(nv*8);
    memcpy(tf.vfits[0],dxfit,nv*8);memcpy(tf.vfits[1],dyfit,nv*8);
    img.sx=sx;img.sy=sy;img.data=(double**)malloc(sy*sizeof(double*));
    for(i=0;i<sy;i++)img.data[i]=in_data+i*sx;
    out.sx=nsx;out.sy=nsy;out.data=(double**)malloc(nsy*sizeof(double*));
    for(i=0;i<nsy;i++)out.data[i]=out_data+i*nsx;
    mk=(char**)malloc(sy*sizeof(char*));
    for(i=0;i<sy;i++)mk[i]=in_mask+i*sx;
    om=(char**)malloc(nsy*sizeof(char*));
    for(i=0;i<nsy;i++)om[i]=out_mask+i*nsx;
    if(method&FT_IG)integrate(&img,mk,&out,om,ofx,ofy,&tf,method,invert);
    else interpolate(&img,mk,&out,om,ofx,ofy,&tf,method,invert);
    free(img.data);free(out.data);free(mk);free(om);
    free(tf.vfits[0]);free(tf.vfits[1]);free(tf.vfits);
    return 0;
}
