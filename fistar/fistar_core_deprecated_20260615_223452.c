/*****************************************************************************/
/* fistar_core.c — extracted from fistar.c, star search + PSF pipeline   */
/* v0.01                                                                     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/stat.h>

#include "image.h"

#include "fitsh.h"
#include "mask.h"
#include "tensor.h"
#include "math/fit/lmfit.h"
#include "math/poly.h"
#include "statistics.h"
#include "index/sort.h"

#include "stars.h"
#include "psf.h"
#include "psf-base.h"
#include "psf-determine.h"
#include "psf-io.h"
#include "magnitude.h"
#include "link/linkpoint.h"
#include "link/floodfill.h"
#include "link/linkblock.h"

#include "fistar.h"
#include "fistar_core.h"

/*****************************************************************************/

extern int	is_verbose, is_comment;
extern char	*progbasename;
extern int fprint_error(char *expr,...);
extern int fprint_warning(char *expr,...);

int fprint_info(char *expr,...)
{
 va_list	ap;
 fprintf(stderr,"%s: ",progbasename);
 va_start(ap,expr);
 vfprintf(stderr,expr,ap);
 va_end(ap);
 fprintf(stderr,"\n");
 return(0);
}

/*****************************************************************************/
/* types extracted from fistar.c                                             */
/*****************************************************************************/

#define		ALG_PP		1
#define		ALG_TRB		2
#define		ALG_LNK		3
#define		ALG_BIQ		4

typedef struct
 {	int	niter;
	int	refinelevel;
	double	bhsize;
 } collectivefit;

typedef struct
 {	int	is_fit_model;
	int	is_determine_psf;
	starfit		sfp;
	starmodelfit	smfs[16];
	int		nsmf;
	range		srcrange;
	spatial		imgbg;
	int		algorithm;
	candidate	*cands;
	int		ncand;
	initcand	*icands;
	int		nicand;
	star		*stars;
	int		nstar;
	collectivefit	cf;
 } starsearch;

/*****************************************************************************/
/* directly copied from fistar.c                                             */
/*****************************************************************************/

int make_star_candidates(initcand *icands,int nicand,double rad,
	image *img,char **mask,candidate **rcands,int *rncand)
{
 candidate	*cands,*wc;
 int		ncand,i,ix1,ix2,iy1,iy2,ix,iy,irad,k,l,t0,t1;
 int		sx,sy;
 double		dx,dy,rad2,dy2,f0;
 initcand	*ic;
 ipoint		*ipoints;
 int		nipoint;
 starshape	shp;

 if ( img==NULL || img->data==NULL )	return(-1);
 sx=img->sx,sy=img->sy;
 if ( sx<=0 || sy<=0 )			return(-1);

 cands=NULL;
 ncand=0;
 
 irad=2+(int)rad;
 rad2=rad*rad;

 for ( i=0 ; i<nicand && icands != NULL ; i++ )
  {	ic=&icands[i];
	ix=(int)floor(ic->x);
	iy=(int)floor(ic->y);
	ix1=ix-irad,ix2=ix+irad;
	iy1=iy-irad,iy2=iy+irad;
	ipoints=(ipoint *)malloc(sizeof(ipoint)*((2*irad+1)*(2*irad+1)));
	nipoint=0,
	t0=t1=0;
	for ( k=iy1 ; k<=iy2 ; k++ )
	 {	
		dy=(ic->y)-((double)k+0.5);
		dy2=dy*dy;
		for ( l=ix1 ; l<=ix2 ; l++ )
		 {	dx=(ic->x)-((double)l+0.5);
			if ( dx*dx+dy2>rad2 )	continue;

			t0++;
			if ( k<0 || sy<=k )	continue;
			if ( l<0 || sx<=l )	continue;
			if ( mask != NULL && mask[k][l] )	continue;
			t1++;

			ipoints[nipoint].x=l,
			ipoints[nipoint].y=k;
			nipoint++;
		 }
	 }
	if ( nipoint==0 || t1<t0 )
	 {	free(ipoints);
		continue;
	 }
	cands=(candidate *)realloc(cands,sizeof(candidate)*(ncand+1));
	wc=&cands[ncand];

	wc->ix=ix,wc->cx=ic->x,
	wc->iy=iy,wc->cy=ic->y;
	
	wc->bg=ic->bg;

	wc->sxx=ic->s+ic->d;
	wc->syy=ic->s-ic->d;
	wc->sxy=ic->k;

	shp.model=SHAPE_ELLIPTIC;
	shp.gs=ic->s,
	shp.gd=ic->d,
	shp.gk=ic->k;
	f0=star_get_unity_flux(&shp);
	if ( f0>0.0 )	wc->peak=wc->amp=ic->flux/f0;
	else		wc->peak=wc->amp=0.0;
	wc->flux=ic->flux;

	wc->ipoints=ipoints;
	wc->nipoint=nipoint;

	wc->area=0.0;

	wc->flags=0;
	wc->marked=0;
	ncand++;
  }

 if ( rcands != NULL )	*rcands=cands;
 if ( rncand != NULL )	*rncand=ncand;

 return(0);
}

/*****************************************************************************/

int fistar_determine_psf(image *img,char **mask,starsearch *ss,psfdetermine *pdparam,psf *tpd)
{
 psfcandidate	*pcands,*wc;
 int		npcand,n,i,nipoint,ii,niter;
 star		*ws;

 if ( is_verbose )
  {	i=pdparam->hsize*2+1;
	fprintf(stderr,"Determination of PSF [%dx%d]x(%dx%d), "
		"spatial order: %d ... ",
		i,i,pdparam->grid,pdparam->grid,pdparam->order);
	fflush(stderr);
  }

 npcand=ss->nstar;
 pcands=(psfcandidate *)malloc(sizeof(psfcandidate)*npcand);

 for ( n=0 ; n<npcand ; n++ )
  {	wc=&pcands[n];
	ws=&ss->stars[n];

	wc->nipoint=nipoint=ws->cand->nipoint;
	wc->ipoints=ws->cand->ipoints;
	wc->yvals=(double *)malloc(sizeof(double)*nipoint);
	memset(wc->yvals,0,sizeof(double)*nipoint);
	drawback_model(wc->ipoints,nipoint,wc->yvals,
		&ws->location,&ws->shape,+1.0);
	wc->bg =ws->location.gbg;
	wc->amp=ws->location.gamp;
	wc->cx =ws->location.gcx;
	wc->cy =ws->location.gcy;
  }

 niter=0;

 for ( ii=0 ; ii<=niter ; ii++ )
  {	psf_determine(img,mask,pcands,npcand,1,pdparam,tpd);
	psf_bgamp_fit(img,mask,pcands,npcand,1,tpd);
	if ( ii<niter )
	 {	for ( n=0 ; n<npcand ; n++ )
		 {	wc=&pcands[n];
			memset(wc->yvals,0,sizeof(double)*wc->nipoint);
			drawback_psf(wc->ipoints,wc->nipoint,wc->yvals,
				wc->cx,wc->cy,wc->amp,tpd,+1.0);
		 }
	 }
  }

 if ( is_verbose )
  {	fprintf(stderr,"done.\n");		}

 for ( n=0 ; n<npcand ; n++ )
  {	wc=&pcands[n];
	ws=&ss->stars[n];
	ws->psf.gcx=wc->cx,
	ws->psf.gcy=wc->cy;
	ws->psf.gamp=wc->amp;
	ws->psf.gbg =wc->bg ;
  }

 for ( n=npcand-1 ; n>=0 ; n-- )
  {	wc=&pcands[n];
	if ( wc->yvals != NULL )	free(wc->yvals);
  }

 free(pcands);

 return(0);
}

/*****************************************************************************/
/* ======================================================================== */
/* flat-array entry point                                                    */
/* ======================================================================== */

/* comparison functions for sort */
int compare_x(int i,int j,void *p)
{ if (((star *)p)[i].location.gcx<((star *)p)[j].location.gcx) return(1);
  if (((star *)p)[i].location.gcx>((star *)p)[j].location.gcx) return(-1); return(0); }
int compare_y(int i,int j,void *p)
{ if (((star *)p)[i].location.gcy<((star *)p)[j].location.gcy) return(1);
  if (((star *)p)[i].location.gcy>((star *)p)[j].location.gcy) return(-1); return(0); }
int compare_amp(int i,int j,void *p)
{ if (((star *)p)[i].location.gamp<((star *)p)[j].location.gamp) return(1);
  if (((star *)p)[i].location.gamp>((star *)p)[j].location.gamp) return(-1); return(0); }
int compare_flux(int i,int j,void *p)
{ if (((star *)p)[i].flux<((star *)p)[j].flux) return(1);
  if (((star *)p)[i].flux>((star *)p)[j].flux) return(-1); return(0); }
int compare_peak(int i,int j,void *p)
{ if (((star *)p)[i].cand->peak<((star *)p)[j].cand->peak) return(1);
  if (((star *)p)[i].cand->peak>((star *)p)[j].cand->peak) return(-1); return(0); }
int compare_fwhm(int i,int j,void *p)
{ if (((star *)p)[i].gfwhm<((star *)p)[j].gfwhm) return(1);
  if (((star *)p)[i].gfwhm>((star *)p)[j].gfwhm) return(-1); return(0); }
int compare_noise(int i,int j,void *p)
{ if (((star *)p)[i].cand==NULL||((star *)p)[j].cand==NULL) return(0);
  if (((star *)p)[i].cand->noise<((star *)p)[j].cand->noise) return(1);
  if (((star *)p)[i].cand->noise>((star *)p)[j].cand->noise) return(-1); return(0); }
int compare_sn(int i,int j,void *p)
{ double sn1,sn2;
  if (((star *)p)[i].cand==NULL||((star *)p)[j].cand==NULL) return(0);
  else if (((star *)p)[i].cand->noise<=0.0||((star *)p)[j].cand->noise<=0.0) return(0);
  else { sn1=((star *)p)[i].flux/((star *)p)[i].cand->noise;
         sn2=((star *)p)[j].flux/((star *)p)[j].cand->noise;
         if (sn1<sn2) return(1); else if (sn1>sn2) return(-1); else return(0); } }

extern int search_star_candidates(image *img, char **mask,
    candidate **rcands, int *rncand, range *srcrange,
    double threshold, spatial *bg, double skysigma);
extern int determine_background(image *img, spatial *bg, int m, int n, int k);
extern int markout_candidates(image *img, char **mask, candidate *cands, int ncand);
extern int cleanup_candlist(candidate **rcands, int *rncand);
extern int convert_candidates(candidate *cands, int ncand, star **rstars, int *rnstar);
extern int cleanup_starlist(star **rstars, int *rnstar);
extern int collective_fit_star_single_model_iterative(image *img, char **mask,
    star *stars, int nstar, ipointlist *ipl,
    starfit *sfp, int level, int niter);
extern int collective_fit_star_single_model_blocked(image *img, char **mask,
    star *stars, int nstar, ipointlist *ipl, double bhsize);
extern int search_star_candidates_link(image *img, char **mask,
    candidate **rcands, int *rncand, range *sr,
    double threshold, double flux_threshold, double critical_prominence);
extern int refine_candidate_params(image *img, candidate *cands, int ncand);
extern int fit_star_single_model(image *img, char **mask,
    candidate *cands, int ncand, star **rstars, int *rnstar,
    starfit *sfp, int model, int order);
extern int free_stars(star *stars, int nstar);

int fistar_search_flat(
    double *img_data, char *mask_data, int sx, int sy,
    double threshold, double flux_threshold, double critical_prominence,
    double skysigma,
    int algorithm, int model, int model_order,
    int iter_symmetric, int iter_general,
    int is_fit_model,
    int collfit_niter, int collfit_refinelevel, double collfit_bhsize,
    int is_determine_psf, int psf_hsize, int psf_grid, int psf_order,
    int psf_type, int psf_symmetrize, int psf_use_biquad,
    double psf_integral_kappa, double psf_circle_width, int psf_circle_order,
    double gain,
    double mag_intensity, double mag_magnitude,
    double *cand_xy, int in_ncand, double cand_radius,
    double *pos_xy, int npos,
    int src_xmin, int src_xmax, int src_ymin, int src_ymax,
    int mark_symbol, int mark_size,
    int out_flags, int sort,
    int is_verbose_arg,
    int *nstar_out, fistar_result *result,
    double *mark_data, double *area_data,
    fistar_psf_out *psf_out,
    double *pos_flux, double *pos_ferr)
{
    int i, r;
    image img;
    char **mask;
    candidate *cands = NULL;
    int ncand = 0, nstar = 0;
    star *stars = NULL;
    starsearch ss;
    starfit sfp;
    psfdetermine pdparam;
    psf tpd;

    progbasename = "fitsh_cy.fistar";
    is_verbose = is_verbose_arg;
    is_comment = 0;

    /* build image row pointers from flat data */
    memset(&img, 0, sizeof(img));
    img.sx = sx; img.sy = sy;
    img.data = malloc(sy * sizeof(double *));
    for (i = 0; i < sy; i++) img.data[i] = img_data + i * sx;

    /* build mask row pointers */
    mask = malloc(sy * sizeof(char *));
    for (i = 0; i < sy; i++) mask[i] = mask_data + i * sx;

    /* mark NaN pixels in mask (same as original ficonv.c main()) */
    mask_mark_nans(&img, mask, MASK_NAN);

    /* ---- initialize star search parameters ---- */
    /* original does not zero ss — fields set individually below
    memset(&ss, 0, sizeof(ss));
    */
    ss.is_fit_model = 1;
    ss.is_determine_psf = is_determine_psf;
    ss.algorithm = algorithm;
    ss.nsmf = 1;
    ss.smfs[0].model = model;
    ss.smfs[0].order = model_order;

    /* original only sets iter_sym and iter_gen on ss.sfp directly
    memset(&sfp, 0, sizeof(sfp));
    sfp.iter_symmetric = 4;
    sfp.iter_general = 2;
    sfp.fit_flags = FIT_XY | FIT_AB | FIT_WIDTH | FIT_DEVIATION;
    ss.sfp = sfp;
    */
    ss.sfp.iter_symmetric = 4;
    ss.sfp.iter_general = 2;

    /* ---- star search ---- */
    /* dump: all input params + image + mask */

    if (cand_xy != NULL && in_ncand > 0) {
        /* input candidates mode: build initcand from flat arrays */
        int ci;
        initcand *icands = (initcand *)malloc(sizeof(initcand) * in_ncand);
        for (ci = 0; ci < in_ncand; ci++) {
            icands[ci].x = cand_xy[ci * 2];
            icands[ci].y = cand_xy[ci * 2 + 1];
            icands[ci].s = 2.0;   /* default shape */
            icands[ci].d = 0.0;
            icands[ci].k = 0.0;
            icands[ci].bg = 0.0;
            icands[ci].flux = 1.0;
        }
        make_star_candidates(icands, in_ncand, cand_radius,
            &img, mask, &cands, &ncand);
        free(icands);
    } else switch (algorithm) {
        /* original does not override sfp or is_fit_model here
        ss.sfp.iter_symmetric = iter_symmetric;
        ss.sfp.iter_general = iter_general;
        ss.is_fit_model = is_fit_model;
        */

        /* build source range */
        range srcrange, *sr = NULL;
        if (src_xmax > src_xmin && src_ymax > src_ymin) {
            srcrange.xmin = src_xmin; srcrange.xmax = src_xmax;
            srcrange.ymin = src_ymin; srcrange.ymax = src_ymax;
            sr = &srcrange;
        }

        case ALG_PP: {
            spatial imgbg;
            /* dump: determine_background input */

            determine_background(&img, &imgbg, 3, 3, 2);

            /* dump: determine_background output */

            /* dump: search_star_candidates_pp input */

            r = search_star_candidates(&img, mask, &cands, &ncand, sr,
                threshold, &imgbg, skysigma);

            /* dump: search_star_candidates_pp output */

            if (r) { r = 1; goto cleanup; }
            markout_candidates(&img, mask, cands, ncand);
            cleanup_candlist(&cands, &ncand);

            /* dump: after cleanup */
            break;
        }
        case ALG_LNK:
            if (flux_threshold > 0.0) threshold = 0.0;

            /* dump: search_star_candidates_link input (same format as original)
            {   char path[512]; FILE *df; int ii;
                mkdir("/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp",0755);
                snprintf(path,sizeof(path),"%s/pipe_search_lnk_in.bin",
                    "/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp");
                df=fopen(path,"wb");
                if(df){
                    fwrite(&sx,4,1,df); fwrite(&sy,4,1,df);
                    fwrite(&threshold,8,1,df);
                    fwrite(&flux_threshold,8,1,df);
                    fwrite(&critical_prominence,8,1,df);
                    for(ii=0;ii<sy;ii++) fwrite(img_data+ii*sx,8,sx,df);
                    for(ii=0;ii<sy;ii++) fwrite(mask_data+ii*sx,1,sx,df);
                    fclose(df);
                }
            }
            */

            r = search_star_candidates_link(&img, mask, &cands, &ncand, sr,
                threshold, flux_threshold, critical_prominence);

            /* dump output: ncand + candidates
            {   char path[512]; FILE *df; int ci;
                snprintf(path,sizeof(path),"%s/pipe_search_lnk_out.bin",
                    "/Users/chaorun/Code/Githubs/fitsh-0.9.4/testgrmatch/tmp");
                df=fopen(path,"wb");
                if(df){
                  fwrite(&ncand,4,1,df);
                  for(ci=0;ci<ncand;ci++){
                    candidate *c=&cands[ci];
                    double a[16]={c->ix,c->iy,c->cx,c->cy,c->peak,c->amp,c->bg,
                      c->sxx,c->syy,c->sxy,c->flux,c->noise,(double)c->nipoint,
                      (double)c->area,(double)c->marked,(double)c->flags};
                    fwrite(a,8,16,df);
                  }
                  fclose(df);
                }
            }
            */

            if (r) { r = 1; goto cleanup; }
            refine_candidate_params(&img, cands, ncand);

            /* dump: after refine */
            break;
        default:
            fprint_error("unsupported algorithm %d", algorithm);
            r = 5; goto cleanup;
    }

    /* ---- fit star models (or collective fit) ---- */

    ss.cf.niter = collfit_niter;
    ss.cf.refinelevel = collfit_refinelevel;
    ss.cf.bhsize = collfit_bhsize;

    if (ss.cf.niter < 0) {
        if (is_fit_model) {
            ss.sfp.fit_flags = FIT_XY | FIT_AB | FIT_WIDTH | FIT_DEVIATION;

            /* DUMP: fit_star_single_model input — commented out
            { ... dump code ... }
            */

            r = fit_star_single_model(&img, mask, cands, ncand, &stars, &nstar,
                &ss.sfp, ss.smfs[0].model, ss.smfs[0].order);

            /* DUMP: fit_star_single_model output — commented out
            { ... dump code ... }
            return 0;
            */

            if (r) { r = 2; goto cleanup; }
        }
    } else {
        ipointlist *ipls;
        int nn;
        convert_candidates(cands, ncand, &stars, &nstar);
        ipls = (ipointlist *)malloc(sizeof(ipointlist) * nstar);
        for (nn = 0; nn < nstar; nn++) {
            ipls[nn].ipoints = cands[nn].ipoints;
            ipls[nn].nipoint = cands[nn].nipoint;
        }
        if (ss.cf.niter > 0) {
            collective_fit_star_single_model_iterative(&img, mask,
                stars, nstar, ipls, &ss.sfp, ss.cf.refinelevel, ss.cf.niter - 1);
        } else {
            collective_fit_star_single_model_blocked(&img, mask,
                stars, nstar, ipls, ss.cf.bhsize);
        }
        cleanup_starlist(&stars, &nstar);
        free(ipls);
    }

    ss.stars = stars;
    ss.nstar = nstar;
    ss.cands = cands;
    ss.ncand = ncand;

    /* dump: fit_star_model output */

    /* ---- PSF determination (optional) ---- */
    if (is_determine_psf) {
        memset(&pdparam, 0, sizeof(pdparam));
        pdparam.hsize = psf_hsize;
        pdparam.grid = psf_grid;
        pdparam.order = psf_order;
        pdparam.type = psf_type;
        pdparam.is_symmetrize = psf_symmetrize;
        pdparam.param.native.use_biquad = psf_use_biquad;
        pdparam.param.integral.kappa = psf_integral_kappa;
        pdparam.param.circle.width = psf_circle_width;
        pdparam.param.circle.order = psf_circle_order;
        memset(&tpd, 0, sizeof(tpd));

        fistar_determine_psf(&img, mask, &ss, &pdparam, &tpd);
    }

    /* ---- extract PSF output ---- */
    if (is_determine_psf && psf_out != NULL) {
        int nvar = (pdparam.order + 1) * (pdparam.order + 2) / 2;
        int nside = (2 * pdparam.hsize + 1) * pdparam.grid;
        psf_out->nvar = nvar;
        psf_out->nside = nside;
        psf_out->hsize = pdparam.hsize;
        psf_out->grid = pdparam.grid;
        psf_out->order = pdparam.order;
        psf_out->ox = 0.5 * (double)sx;
        psf_out->oy = 0.5 * (double)sy;
        psf_out->scale = 0.5 * (double)sx;
        psf_out->data = (double *)malloc(nvar * nside * nside * sizeof(double));
        if (tpd.coeff != NULL) {
            int vi, yi, xi;
            for (vi = 0; vi < nvar; vi++)
                for (yi = 0; yi < nside; yi++)
                    for (xi = 0; xi < nside; xi++)
                        psf_out->data[vi * nside * nside + yi * nside + xi] =
                            tpd.coeff[vi][yi][xi];
        } else {
            memset(psf_out->data, 0, nvar * nside * nside * sizeof(double));
        }
    }

    /* ---- extract results to flat arrays ---- */
    *nstar_out = nstar;

    /* ---- position matching ---- */
    if (pos_xy != NULL && npos > 0 && pos_flux != NULL) {
        int pi, ci, jj;
        int **refarr = (int **)tensor_alloc_2d(int, sx, sy);
        for (i = 0; i < sy; i++)
            for (jj = 0; jj < sx; jj++)
                refarr[i][jj] = -1;
        for (ci = 0; ci < ncand; ci++)
            for (jj = 0; jj < cands[ci].nipoint; jj++)
                refarr[cands[ci].ipoints[jj].y][cands[ci].ipoints[jj].x] = ci;
        for (pi = 0; pi < npos; pi++) {
            int ix = (int)floor(pos_xy[pi*2]), iy = (int)floor(pos_xy[pi*2+1]);
            pos_flux[pi] = 0.0;
            if (pos_ferr) pos_ferr[pi] = 0.0;
            if (ix >= 0 && ix < sx && iy >= 0 && iy < sy && refarr[iy][ix] >= 0) {
                candidate *wc = &cands[refarr[iy][ix]];
                if (wc->flux > 0.0) {
                    pos_flux[pi] = wc->flux;
                    if (pos_ferr)
                        pos_ferr[pi] = sqrt(wc->flux / gain + wc->noise * wc->nipoint);
                }
            }
        }
        tensor_free(refarr);
    }

    /* ---- sort if requested ---- */
    if (sort >= 0 && nstar > 0) {
        int *indx = (int *)malloc(sizeof(int) * nstar);
        for (i = 0; i < nstar; i++) indx[i] = i;
        typedef int (*cmp_fn)(int, int, void *);
        cmp_fn cmps[] = { (cmp_fn)compare_x, (cmp_fn)compare_y, (cmp_fn)compare_peak, (cmp_fn)compare_fwhm, (cmp_fn)compare_amp, (cmp_fn)compare_flux, (cmp_fn)compare_noise, (cmp_fn)compare_sn };
        if (sort < 8 && cmps[sort])
            index_qsort(indx, nstar, cmps[sort], (void *)stars);
        /* reorder stars array */
        star *sorted_stars = (star *)malloc(sizeof(star) * nstar);
        for (i = 0; i < nstar; i++) sorted_stars[i] = stars[indx[i]];
        free(stars); stars = sorted_stars;
        free(indx);
    }

    #define ALLOC(N) calloc(nstar, sizeof(*(result->N)))
    result->ix    = (double *)ALLOC(ix);
    result->iy    = (double *)ALLOC(iy);
    result->cx    = ALLOC(cx);
    result->cy    = ALLOC(cy);
    result->cbg   = ALLOC(cbg);
    result->camp  = ALLOC(camp);
    result->cmax  = ALLOC(cmax);
    result->npix  = (int *)calloc(nstar, sizeof(int));
    result->cs    = ALLOC(cs);
    result->cd    = ALLOC(cd);
    result->ck    = ALLOC(ck);
    result->x     = ALLOC(x);
    result->y     = ALLOC(y);
    result->bg    = ALLOC(bg);
    result->amp   = ALLOC(amp);
    result->s     = ALLOC(s);
    result->d     = ALLOC(d);
    result->k     = ALLOC(k);
    result->l     = ALLOC(l);
    result->sigma = ALLOC(sigma);
    result->delta = ALLOC(delta);
    result->kappa = ALLOC(kappa);
    result->fwhm  = ALLOC(fwhm);
    result->ellip = ALLOC(ellip);
    result->pa    = ALLOC(pa);
    result->flux  = ALLOC(flux);
    result->noise = ALLOC(noise);
    result->sn    = ALLOC(sn);
    result->magnitude = ALLOC(magnitude);
    result->px    = ALLOC(px);
    result->py    = ALLOC(py);
    result->pbg   = ALLOC(pbg);
    result->pamp  = ALLOC(pamp);
    result->ps    = ALLOC(ps);
    result->pd    = ALLOC(pd);
    result->pk    = ALLOC(pk);
    result->pl    = ALLOC(pl);
    #undef ALLOC

    for (i = 0; i < nstar; i++) {
        star     *ws  = &stars[i];
        candidate *wc = ws->cand;

        /* candidate fields */
        result->ix[i]   = wc ? wc->ix + 1 : 0;
        result->iy[i]   = wc ? wc->iy + 1 : 0;
        result->cx[i]   = wc ? wc->cx : 0.0;
        result->cy[i]   = wc ? wc->cy : 0.0;
        result->cbg[i]  = wc ? wc->bg : 0.0;
        result->camp[i] = wc ? wc->amp : 0.0;
        result->cmax[i] = wc ? wc->peak : 0.0;
        result->npix[i] = wc ? wc->nipoint : 0;
        result->cs[i]   = wc ? 0.5 * (wc->sxx + wc->syy) : 0.0;
        result->cd[i]   = wc ? 0.5 * (wc->sxx - wc->syy) : 0.0;
        result->ck[i]   = wc ? wc->sxy : 0.0;

        /* fitted location */
        result->x[i]   = ws->location.gcx;
        result->y[i]   = ws->location.gcy;
        result->bg[i]  = ws->location.gbg;
        result->amp[i] = ws->location.gamp;

        /* fitted shape */
        result->s[i]     = ws->shape.gs;
        result->d[i]     = ws->shape.gd;
        result->k[i]     = ws->shape.gk;
        result->l[i]     = ws->shape.gl;
        result->sigma[i] = ws->gsig;
        result->delta[i] = ws->gdel;
        result->kappa[i] = ws->gkap;
        result->fwhm[i]  = ws->gfwhm;
        result->ellip[i] = ws->gellip;
        result->pa[i]    = ws->gpa;

        /* flux / noise */
        result->flux[i]  = ws->flux;
        result->noise[i] = wc ? wc->noise : 0.0;
        result->sn[i]    = (wc && wc->noise > 0.0) ? ws->flux / wc->noise : 0.0;
        /* magnitude (match fistar-io.c fprint_star_mag) */
        result->magnitude[i] = (ws->flux > 0.0)
            ? mag_magnitude - 2.5 * log10(ws->flux / mag_intensity) : 0.0;

        /* PSF */
        result->px[i]   = ws->psf.gcx;
        result->py[i]   = ws->psf.gcy;
        result->pbg[i]  = ws->psf.gbg;
        result->pamp[i] = ws->psf.gamp;
        result->ps[i]   = 1.0;
        result->pd[i]   = 0.0;
        result->pk[i]   = 0.0;
        result->pl[i]   = 0.0;
    }

    /* ---- extract PSF data ---- */
    if (is_determine_psf && psf_out != NULL) {
        psf_out->hsize = pdparam.hsize;
        psf_out->grid = pdparam.grid;
        psf_out->order = pdparam.order;
        psf_out->ox = 0.5 * (double)sx;
        psf_out->oy = 0.5 * (double)sy;
        psf_out->scale = 0.5 * (double)sx;
    }

    /* ---- free internal structures ---- */
    if (stars) free_stars(stars, nstar);
    /* cands will be freed by free_stars since stars point to cands */
    /* PSF data owned by tpd, freed by psf_free */

    /* ---- mark / area output ---- */
    if (mark_data != NULL && (out_flags & 1)) {
        memcpy(mark_data, img_data, sx * sy * sizeof(double));
    }
    if (area_data != NULL && (out_flags & 2)) {
        memcpy(area_data, img_data, sx * sy * sizeof(double));
    }

    /* ---- fallback: if stars is NULL after fit, convert from cands ---- */
    if (stars == NULL && is_fit_model)
        convert_candidates(cands, ncand, &stars, &nstar);

    /* ---- draw marks ---- */
    if (mark_data != NULL && (out_flags & 1) && nstar > 0) {
        int mi, ix, iy, d, dj;
        double color = 0.0;
        for (mi = 0; mi < nstar; mi++) {
            if (stars[mi].cand == NULL) continue;
            ix = stars[mi].cand->ix; iy = stars[mi].cand->iy;
            /* mark drawing: 0=dot, 1=square, 2=circle */
            if (mark_symbol == 0) {
                mark_data[iy * sx + ix] = color;
            } else if (mark_symbol == 1) {
                for (d = -mark_size; d <= mark_size; d++) {
                    if (iy+d >= 0 && iy+d < sy) {
                        for (dj = -mark_size; dj <= mark_size; dj++) {
                            if (ix+dj >= 0 && ix+dj < sx)
                                mark_data[(iy+d) * sx + ix+dj] = color;
                        }
                    }
                }
            } else {
                /* default: cross marker (mark_symbol == 2 or unknown) */
                mark_data[iy * sx + ix] = color;
                for (d = 1; d <= mark_size; d++) {
                    if (ix+d < sx) mark_data[iy * sx + ix+d] = color;
                    if (ix-d >= 0) mark_data[iy * sx + ix-d] = color;
                    if (iy+d < sy) mark_data[(iy+d) * sx + ix] = color;
                    if (iy-d >= 0) mark_data[(iy-d) * sx + ix] = color;
                }
            }
        }
    }
    if (area_data != NULL && (out_flags & 2) && nstar > 0) {
        int ii, jj, ix, iy;
        for (ii = 0; ii < nstar; ii++) {
            if (stars[ii].cand == NULL) continue;
            for (jj = 0; jj < stars[ii].cand->nipoint; jj++) {
                ix = stars[ii].cand->ipoints[jj].x;
                iy = stars[ii].cand->ipoints[jj].y;
                area_data[iy * sx + ix] = 0.0;
            }
        }
    }

    r = 0;

cleanup:
    free(img.data);
    free(mask);
    return r;
}

