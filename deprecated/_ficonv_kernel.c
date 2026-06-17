/*****************************************************************************/
/* _ficonv_kernel.c — parser only (create_kernels_from_kernelarg).         */
/*****************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "kernel.h"

static void _rm_spaces(char *s) {
    int i,j; for(i=j=0;s[i];i++) if(s[i]!=' '&&s[i]!='\t'&&s[i]!='\n'&&s[i]!='\r')s[j++]=s[i];s[j]=0;
}
static int _tok(char *s,char **c,char sep,int n) {
    int m=0;c[0]=s;m++;for(;*s;s++)if(*s==sep){*s=0;if(m<n)c[m++]=s+1;}for(;m<n;m++)c[m]=NULL;return m;
}

/* ---- add_kernel_* functions (from kernel-io.c) ---- */
static int _add_bg(kernellist *kl,int o) {
    kernel *k;kl->kernels=(kernel*)realloc(kl->kernels,(kl->nkernel+1)*sizeof(kernel));
    k=&kl->kernels[kl->nkernel];kl->nkernel++;
    memset(k,0,sizeof(kernel));k->type=KERNEL_BACKGROUND;k->order=o;k->hsize=0;k->flag=1;return 0;
}
static int _add_id(kernellist *kl,int o) {
    kernel *k;kl->kernels=(kernel*)realloc(kl->kernels,(kl->nkernel+1)*sizeof(kernel));
    k=&kl->kernels[kl->nkernel];kl->nkernel++;
    memset(k,0,sizeof(kernel));k->type=KERNEL_IDENTITY;k->order=o;k->hsize=0;k->flag=1;return 0;
}
static int _add_shift(kernellist *kl,int o){return 0;}
static int _add_gauss_set(kernellist *kl,double sigma,int order,int hsize,int sp) {
    int i,j,nv=(order+1)*(order+2)/2,di,dj;
    kernel *k;kl->kernels=(kernel*)realloc(kl->kernels,(kl->nkernel+nv)*sizeof(kernel));
    for(di=0,k=&kl->kernels[kl->nkernel];di<=order;di++)for(dj=0;dj<=di;dj++,k++) {
        memset(k,0,sizeof(kernel));k->type=KERNEL_GAUSSIAN;k->order=0;
        k->hsize=hsize;k->sigma=sigma;k->bx=di;k->by=dj;k->flag=1;k->target=0;}
    /* calculate kernel images via numerical integration over subpixels */
    int subg=10,hs=hsize;for(di=0,k=&kl->kernels[kl->nkernel];di<nv;di++,k++){
        if(!k->type){continue;}
        k->image=(double**)malloc((2*hs+1)*sizeof(double*));
        for(i=0;i<=2*hs;i++){k->image[i]=(double*)malloc((2*hs+1)*sizeof(double));}
        for(i=0;i<=2*hs;i++)for(j=0;j<=2*hs;j++){
            double w=0;int m,n;
            for(m=0;m<subg;m++){double dy=(1+2*m-subg)/(2.0*subg);
                for(n=0;n<subg;n++){double dx=(1+2*n-subg)/(2.0*subg);
                    double x=j-hs+dx,y=i-hs+dy;
                    double g=exp(-(x*x+y*y)/(2*sigma*sigma));
                    double bx=pow(hs?x:1,di),by=pow(hs?y:1,dj);
                    w+=g*bx*by;}}
            k->image[i][j]=w/(subg*subg);}
        if(di==0&&dj==0){/* normalize and subtract identity at center */
            double s=0;for(i=0;i<=2*hs;i++)for(j=0;j<=2*hs;j++)s+=k->image[i][j];
            s=1.0/s;for(i=0;i<=2*hs;i++)for(j=0;j<=2*hs;j++)k->image[i][j]*=s;
            k->image[hs][hs]-=1.0;}}
    kl->nkernel+=nv;return 0;
}
static int _add_linear_set(kernellist *kl,int hsize,int sp) {
    int kr;for(kr=1;kr<=hsize;kr++){
        kernel *k;kl->kernels=(kernel*)realloc(kl->kernels,(kl->nkernel+2)*sizeof(kernel));
        k=&kl->kernels[kl->nkernel];kl->nkernel++;
        memset(k,0,sizeof(kernel));k->type=KERNEL_DDELTA;k->hsize=kr;k->bx=kr;k->by=0;k->order=0;k->flag=1;
        k=&kl->kernels[kl->nkernel];kl->nkernel++;
        memset(k,0,sizeof(kernel));k->type=KERNEL_DDELTA;k->hsize=kr;k->bx=0;k->by=kr;k->order=0;k->flag=1;}
    return 0;
}

int create_kernels_from_kernelarg(char *kernelarg, kernellist *kl) {
    char *ka,*cmd[32],*icmd[4],*jcmd[4];
    int bg=-1,id=-1,g=-1,d=-1,s=-1,i,j,ls=0;
    kl->nkernel=0;kl->kernels=NULL;kl->ox=kl->oy=0;kl->scale=1;kl->type=0;
    ka=(char*)malloc(strlen(kernelarg)+1);strcpy(ka,kernelarg);_rm_spaces(ka);
    _tok(ka,cmd,';',31);
    for(i=0;cmd[i];i++){if(cmd[i][0]=='i')id=i;else if(cmd[i][0]=='b')bg=i;
        else if(cmd[i][0]=='g')g=i;else if(cmd[i][0]=='d')d=i;
        else if(cmd[i][0]=='s')s=i;else{free(ka);return 1;}}
    if(bg>=0){int o=-1;_tok(cmd[bg],jcmd,'/',2);o=jcmd[1]?atoi(jcmd[1]):0;if(o<0){free(ka);return 1;}_add_bg(kl,o);}
    if(id>=0){int o=-1;_tok(cmd[id],jcmd,'/',2);o=jcmd[1]?atoi(jcmd[1]):0;if(o<0){free(ka);return 1;}_add_id(kl,o);}
    if(s>=0){int o=-1;_tok(cmd[s],jcmd,'/',2);o=jcmd[1]?atoi(jcmd[1]):0;if(o<0){free(ka);return 1;}_add_shift(kl,o);}
    for(i=0;cmd[i];i++){double si;int hs,or,sp;if(cmd[i][0]!='g')continue;
        _tok(cmd[i],icmd,'=',2);if(!icmd[1]){free(ka);return 1;}sp=0;
        if(sscanf(icmd[1],"%d,%lg,%d/%d",&hs,&si,&or,&sp)<3){free(ka);return 1;}
        if(sp<0){free(ka);return 1;}_add_gauss_set(kl,si,or,hs,sp);}
    for(i=0;cmd[i];i++){int hs,sp;if(cmd[i][0]!='d')continue;
        _tok(cmd[i],icmd,'=',2);if(!icmd[1]){free(ka);return 1;}sp=0;
        if(sscanf(icmd[1],"%d/%d",&hs,&sp)<1){free(ka);return 1;}if(sp<0){free(ka);return 1;}
        for(j=ls+1;j<=hs;j++)_add_linear_set(kl,j,sp);ls=hs;}
    free(ka);return 0;
}
