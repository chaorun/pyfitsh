/*****************************************************************************/
/* lfit_interface.c — Python-callable wrapper for lfit (no FITS I/O)         */
/*                                                                           */
/* Strategy: #include "lfit.c" to bring in all structs + algorithms.         */
/* Replace main() + FILE* outputs with lfit_python_apply() + fmemopen().      */
/*****************************************************************************/

/* Pull in the entire lfit core (structs, fit_* functions, PSN helpers) */
#include "lfit.c"

/*****************************************************************************/
/* New types for Python interface                                            */
/*****************************************************************************/

typedef struct
 {	double	**array_ptrs;
	int	*array_nrows;
	int	*array_ncols;
 } lfit_input_arrays;

typedef struct
 {	double	*params;
	double	*errors;
	double	chi2;
	double	*chain;
	int	chain_count;
	int	chain_max;
	int	nvar_chain;
	int	*used_mask;
	int	nrow;
	int	nused;
	double	residual_sigma;
	double	acceptance;
	double	**cov_matrix;
	double	**corr_matrix;
	double	*eval_data;
	int	eval_nrow;
	int	eval_ncol;
	int	error_code;
	char	error_msg[256];
 } lfit_result;

typedef void (*lfit_step_callback)(double *bvector,int nvar,double chi2,void *user_data);

/*****************************************************************************/
/* load_fit_data_from_array: single datablock, replaces read_fit_data        */
/*****************************************************************************/

int load_fit_data_from_array(
	double *array_data,int nrow_in,int ncol_in,
	int dbidx,int nvar,int maxncol,
	int *rnrow,
	double **rfitdata,double **rfitdep,double **rfitwght,
	fitinputrow **rfitinputrows,
	psn *pdep,psn *perr,int errtype)
{
 int		i,j,nrow;
 double		*fitdata,*fitdep,*fitwght;
 fitinputrow	*fitinputrows;
 double		f,weight;
 double		*wvars;

 wvars=(double *)malloc(sizeof(double)*(nvar+maxncol));

 nrow=*rnrow;
 fitdata =*rfitdata;
 fitdep  =*rfitdep;
 fitwght =*rfitwght;
 fitinputrows=*rfitinputrows;

 fitdata =(double *)realloc(fitdata ,(nrow+nrow_in)*sizeof(double)*maxncol);
 fitdep  =(double *)realloc(fitdep  ,(nrow+nrow_in)*sizeof(double));
 fitwght =(double *)realloc(fitwght ,(nrow+nrow_in)*sizeof(double));
 fitinputrows=(fitinputrow *)realloc(fitinputrows,(nrow+nrow_in)*sizeof(fitinputrow));

 for ( i=0 ; i<nrow_in ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	wvars[j]=0;		}
	for ( j=0 ; j<ncol_in && j<maxncol ; j++ )
	 {	wvars[j+nvar]=array_data[i*ncol_in+j];	}

	lfit_psn_double_calc(pdep,NULL,lpg.pl_funct,&f,wvars);

	if ( perr != NULL )
	 {	lfit_psn_double_calc(perr,NULL,lpg.pl_funct,&weight,wvars);
		if ( errtype )		weight=weight;
		else if ( weight>0.0 )	weight=1.0/weight;
		else			weight=0.0;
	 }
	else	weight=1.0;

	j=nrow*maxncol;
	for ( j=0 ; j<ncol_in && j<maxncol ; j++ )
	 {	fitdata[nrow*maxncol+j]=array_data[i*ncol_in+j];	}
	for ( ; j<maxncol ; j++ )
	 {	fitdata[nrow*maxncol+j]=0.0;	}

	fitdep [nrow]=f;
	fitwght[nrow]=weight;
	fitinputrows[nrow].line =NULL;
	fitinputrows[nrow].dbidx=dbidx;
	fitinputrows[nrow].x    =NULL;

	nrow++;
  }

 free(wvars);

 *rnrow=nrow;
 *rfitdata =fitdata;
 *rfitdep  =fitdep;
 *rfitwght =fitwght;
 *rfitinputrows=fitinputrows;

 return(0);
}

/*****************************************************************************/
/* load_all_fit_data_from_arrays: all datablocks, replaces read_all_fit_data */
/*****************************************************************************/

int load_all_fit_data_from_arrays(lfitdata *lf,int nvar,
	lfit_input_arrays *input,
	double **rfitdata,double **rfitdep,double **rfitwght,
	fitinputrow **rfitinputrows,int *rnrow)
{
 int		d;
 datablock	*db;

 for ( d=0 ; d<lf->ndatablock ; d++ )
  {	db=&lf->datablocks[d];
	db->offset=*rnrow;

	load_fit_data_from_array(
		input->array_ptrs[d],input->array_nrows[d],input->array_ncols[d],
		d,nvar,lf->maxncol,
		rnrow,rfitdata,rfitdep,rfitwght,rfitinputrows,
		db->dep,db->err,db->errtype);

	db->size=(*rnrow)-db->offset;
  }

 return(0);
}

/*****************************************************************************/
/* Functions moved from lfit.c (contain FILE* I/O, to be modified later)    */
/*****************************************************************************/

int fit_linear_or_nonlinear(lfitdata *lf,lfit_result *result,
	variable *vars,int nvar,int errtype,
	int errdump,int resdump,int is_dump_delta,int is_weighted_sigma,int is_numeric_derivs,
	lfit_input_arrays *input)
{
 int		i,j,k,n,isig,nreject,nc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fiteval,*fitwght,lam,lam_mpy;
 double		*wvars;
 fitinputrow	*fitinputrows;
 int		nrow;
 double		**amatrix,*bvector,*evector,*xvector;
 fitfunctdata	ffd;
 constraint	ccstat,*cc;
 fitparam	*fp;

 fp=&lf->parameters;

 if ( lf->fconstraint.nc>0 )
  {	nc= ccstat.nc= lf->fconstraint.nc;
	ccstat.cmatrix=lf->fconstraint.cmatrix;
	ccstat.cvector=lf->fconstraint.cvector;
 	cc=&ccstat;
  }
 else
  {	cc=NULL;
	nc=0;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

  // fprintf(stderr,"DBG: fit_lon: calling load_all nvar=%d ndatablock=%d maxncol=%d\n",nvar,lf->ndatablock,lf->maxncol); fflush(stderr);
 load_all_fit_data_from_arrays(lf,nvar,input,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);
  // fprintf(stderr,"DBG: fit_lon: load_all done nrow=%d\n",nrow); fflush(stderr);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	snprintf(result->error_msg,sizeof(result->error_msg),"too few lines for fitting");
	result->error_code=-1;
	return(-1);
  }

 for ( i=0 ; i<nrow ; i++ )
  {	fitwght[i] = fitwght[i]*fitwght[i];		}
  // fprintf(stderr,"DBG: fit_lon: fitwght squared, allocating matrices\n"); fflush(stderr);

 amatrix=matrix_alloc(nvar);
 bvector=vector_alloc(nvar);
 evector=vector_alloc(nvar);
 xvector=vector_alloc(nvar);

  // fprintf(stderr,"DBG: fit_lon: matrices allocated, setting initial values\n"); fflush(stderr);

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=0.0;
	xvector[i]=vars[i].diff;
	for ( j=0 ; j<nvar ; j++ )	amatrix[i][j]=0.0;
  } 
  // fprintf(stderr,"DBG: fit_lon: initial values set\n"); fflush(stderr);

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);
  // fprintf(stderr,"DBG: fit_lon: entering isig loop niter=%d is_linear=%d\n",fp->niter,lf->is_linear); fflush(stderr);

 for ( isig=0 ; isig<=fp->niter ; isig++ )
  {	
  // fprintf(stderr,"DBG: fit_lon: isig=%d\n",isig); fflush(stderr);
	if ( lf->is_linear )
	 {	if ( errdump<=0 )
 {	// fprintf(stderr,"DBG: fit_lon: calling lin_fit_con nvar=%d nrow=%d cc=%p\n",nvar,nrow,(void*)cc); fflush(stderr);
			i=lin_fit_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,NULL);
  // fprintf(stderr,"DBG: fit_lon: lin_fit_con returned i=%d\n",i); fflush(stderr);
 }
		else
			i=lin_fit_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,evector);
		if ( i )		
		 {	fprint_error("singular matrix");
			snprintf(result->error_msg,sizeof(result->error_msg),"singular matrix");
			result->error_code=-1;
			return(-1);
		 }
	 }
	else 
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	bvector[i]=vars[i].init;		}
		lam    =fp->tune.nllm.lambda;
		lam_mpy=fp->tune.nllm.lambda_mpy;
		for ( n=0 ; n<fp->tune.nllm.max_iter ; n++ )
		 {	if ( is_verbose>=2 )
			 {	int vi;
				for ( vi=0 ; vi<nvar ; vi++ )
				 {	fprintf(stderr,vars[vi].format,bvector[vi]);
					fprintf(stderr," ");
				 }
				fprintf(stderr,"(%11g)\n",lam);
			 }
			if ( ! is_numeric_derivs)
				lam=nlm_fit_base_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,lam,lam_mpy);
			else
				lam=nlm_fit_nmdf_con(fitpnt,fitdep,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,lam,lam_mpy,xvector);
		 }
	 }

	if ( isig<fp->niter )
	 {	double	d,s,sdd,sig,maxdev,w;
		
		s=sdd=0.0;
		for ( i=0 ; i<nrow ; i++ )
		 {	fit_function(fitpnt[i],bvector,&fiteval[i],NULL,&ffd);
			d=fiteval[i]-fitdep[i];
			if ( is_weighted_sigma )
			 {	w=fitwght[i];
				s  +=w;
				sdd+=d*d*w;
			 }
			else
			 {	s+=1.0;
				sdd+=d*d;
			 }
		 }
		sdd/=s;
		sig =sqrt(sdd);

		maxdev=sig*fp->sigma;
		for ( i=0,j=0,k=0 ; i<nrow ; i++ )
		 {	if ( fabs(fiteval[i]-fitdep[i]) > maxdev && fitwght[i]>0.0 )
			 {	fitwght[i]=0.0;
				j++;
			 }
			if ( fitwght[i]>0.0 )
				k++;
		 }
		nreject=j;
		if ( is_verbose>0 )
			fprintf(stderr,"Iteration#%3d: rejected/be used/all: %4d/%4d/%4d\n",isig+1,j,k,nrow);
	 }
	else
		nreject=0;

	if ( is_verbose>0 && fp->niter>0 && ( isig==fp->niter || nreject==0 ) )
	 {	for ( i=0,k=0 ; i<nrow ; i++ )
		 {	if ( fitwght[i]>0.0 )	k++;	}
		fprintf(stderr,"In the last  : used/all: %4d/%4d\n",k,nrow);
	 }
	
	if ( ! ( nreject>0 ) )	break;
  }

 if ( ! lf->is_linear )	errdump=0;

 for ( i=0 ; i<nvar ; i++ )
  {	result->params[i]=bvector[i];	}
 if ( errdump>0 )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	result->errors[i]=evector[i];	}
  }
 {	double	d,sdd;
	sdd=0.0;
	for ( i=0 ; i<nrow ; i++ )
	 {	fit_function(fitpnt[i],bvector,&fiteval[i],NULL,&ffd);
		d=fiteval[i]-fitdep[i];
		sdd+=d*d;
	 }
	if ( nrow>nvar )	sdd/=(double)(nrow-nvar);
	else			sdd=0.0;
	result->residual_sigma=sqrt(sdd);
	result->chi2=sdd;
 }
 result->nrow=nrow;
 for ( i=0,k=0 ; i<nrow ; i++ )
  {	if ( fitwght[i]>0.0 )	k++;	}
 result->nused=k;
 if ( result->used_mask != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	result->used_mask[i]=(fitwght[i]>0.0);	}
  }
 result->error_code=0;

 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);

 vector_free(xvector);
 vector_free(evector);
 vector_free(bvector);
 matrix_free(amatrix);

 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 return(0);
}

int fit_markov_chain_monte_carlo(lfitdata *lf,lfit_result *result,variable *vars,int nvar,
	lfit_input_arrays *input,lfit_step_callback callback,void *callback_data)
{
 int		i,n,nc,nrc,ccnt,is_not_accepted;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fiteval,*fitwght;
 double		*wvars,chi2,nchi,mchi,vx,u,alpha,dalpha;
 fitinputrow	*fitinputrows;
 int		nrow,cgibbs;
 double		*bvector,*dvector,*evector,*mvector;
 fitfunctdata	ffd;
 fitparam	*fp;
 
 fp=&lf->parameters;

 nc=lf->fconstraint.nc;
 nrc=nc;
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);
 dvector=vector_alloc(nvar);
 evector=vector_alloc(nvar);
 mvector=vector_alloc(nvar);

 load_all_fit_data_from_arrays(lf,nvar,input,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	snprintf(result->error_msg,sizeof(result->error_msg),"too few lines for fitting");
	result->error_code=-1;
	return(-1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;			}

 ccnt=1;
 result->chain_count=0;

 if ( nc > 0 )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].flags & VAR_IS_CONSTANT )
	 	 {	bvector[i]=vars[i].init;		}
	 }
	if ( nrc>0 )
	 {	constraint_force(bvector,nvar,nc,
		lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
	 }
  }

 chi2=data_get_chisquare(fitpnt,fitdep,fitwght,fiteval,nrow,bvector,&ffd,0);
 mchi=chi2;
 for ( i=0 ; i<nvar ; i++ )
  {	mvector[i]=bvector[i];		}

 if ( result->chain != NULL && result->chain_count < result->chain_max )
  {	int stride=nvar+1;
	memcpy(&result->chain[result->chain_count*stride],bvector,nvar*sizeof(double));
	result->chain[result->chain_count*stride+nvar]=chi2;
	result->chain_count++;
  }
 if ( callback != NULL )
	callback(bvector,nvar,chi2,callback_data);

 cgibbs=0;

 for ( n=0 ; (lf->parameters.tune.mcmc.count_accepted?ccnt:n)<fp->mc_iterations ; n++ )
  {	
	is_not_accepted=0;
	if ( ! fp->tune.mcmc.use_gibbs )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( ! (vars[i].flags & VAR_IS_CONSTANT) )
			 {	dvector[i]=lfit_get_gaussian(bvector[i],vars[i].ierr);
				if ( (vars[i].flags & VAR_MIN) && dvector[i] < vars[i].imin )
				 {	is_not_accepted=1;
					dvector[i]=vars[i].imin;
				 }
				if ( (vars[i].flags & VAR_MAX) && dvector[i] > vars[i].imax )
				 {	is_not_accepted=1;
					dvector[i]=vars[i].imax;
				 }
			 }
			else
				dvector[i]=bvector[i];
		 }
	 }
	else 
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	dvector[i]=bvector[i];		}

		while ( vars[cgibbs].flags & VAR_IS_CONSTANT )
			cgibbs=(cgibbs+1)%nvar;

		i=cgibbs;
		dvector[i]=lfit_get_gaussian(bvector[i],vars[i].ierr);
		if ( (vars[i].flags & VAR_MIN) && dvector[i] < vars[i].imin )
		 {	is_not_accepted=1;
			dvector[i]=vars[i].imin;
		 }
		if ( (vars[i].flags & VAR_MAX) && dvector[i] > vars[i].imax )
		 {	is_not_accepted=1;
			dvector[i]=vars[i].imax;
		 }

		cgibbs=(cgibbs+1)%nvar;
	 }

	if ( is_not_accepted )
		continue;

	if ( nc > 0 )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( vars[i].flags & VAR_IS_CONSTANT )
		 	 {	dvector[i]=vars[i].init;		}
		 }
		if ( nrc>0 )
		 {	constraint_force(dvector,nvar,nc,
			lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
		 }
	 }

	nchi=data_get_chisquare(fitpnt,fitdep,fitwght,fiteval,nrow,dvector,&ffd,0);

	if ( nchi<mchi )
	 {	mchi=nchi;
		for ( i=0 ; i<nvar ; i++ )
		 {	mvector[i]=dvector[i];		}
	 }

	vx=exp(-0.5*(nchi-chi2));
	alpha=(vx<1.0?vx:1.0);

	if ( alpha>0.0 && (u=random_double()) <= alpha )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	bvector[i]=dvector[i];		}
		chi2=nchi;
		ccnt++;
		if ( result->chain != NULL && result->chain_count < result->chain_max )
		 {	int stride=nvar+1;
			memcpy(&result->chain[result->chain_count*stride],bvector,nvar*sizeof(double));
			result->chain[result->chain_count*stride+nvar]=chi2;
			result->chain_count++;
		 }
		if ( callback != NULL )
			callback(bvector,nvar,chi2,callback_data);
	 }
	else
		is_not_accepted=1;	
  }

 if ( n>0 )
  {	alpha=(double)ccnt/(double)n;
	dalpha=sqrt((double)ccnt)/(double)n;
  }
 else
  {	alpha=0.0;
	dalpha=0.0;
  }

 mchi=data_get_chisquare(fitpnt,fitdep,fitwght,fiteval,nrow,mvector,&ffd,1);
 result->chi2=mchi;
 result->acceptance=alpha;
 for ( i=0 ; i<nvar ; i++ )
  {	result->params[i]=mvector[i];	}
 result->nrow=nrow;
 result->error_code=0;

 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(mvector);
 vector_free(evector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}

int fit_map_chi2(lfitdata *lf,lfit_result *result,variable *vars,int nvar,
	lfit_input_arrays *input,lfit_step_callback callback,void *callback_data)
{
 int		i,j,n;
 void		**fitpnt;
 double		*fitdata,*fitdep,*fitwght;
 double		*wvars,chi2;
 fitinputrow	*fitinputrows;
 int		nrow,ndim,ngrid,igrid;
 mapaxis	*axes,*cm;
 gridpoint	*griddata;
 
 double		*bvector;
 fitfunctdata	ffd;
 
 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 ndim=0;
 ngrid=0;
 axes=NULL;
 for ( i=0 ; i<nvar ; i++ )
  {	if ( ! (vars[i].flags & VAR_STEP) )
		continue;
	n=1+(int)((vars[i].imax+0.5*vars[i].istp-vars[i].imin)/vars[i].istp);
	if ( n <= 0 )
		continue;
	axes=(mapaxis *)realloc(axes,sizeof(mapaxis)*(ndim+1));
	cm=&axes[ndim];
	cm->n=n;
	cm->curr=0;
	cm->imin=vars[i].imin;
	cm->istp=vars[i].istp;
	cm->varindex=i;
	ndim++;
	if ( ngrid )	ngrid=ngrid*cm->n;
	else		ngrid=cm->n;
  }
 if ( ngrid>0 && ndim>0 )
	griddata=(gridpoint *)malloc(sizeof(gridpoint)*ngrid);
 else
	griddata=NULL;

 bvector=vector_alloc(nvar);

 load_all_fit_data_from_arrays(lf,nvar,input,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 for ( i=0 ; i<nvar ; i++ )
  {	if ( vars[i].flags & VAR_STEP )
		bvector[i]=vars[i].imin;
	else
		bvector[i]=vars[i].init;
  }

 igrid=0;
 result->chain_count=0;

 while ( 1 )
  {	
	chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,bvector,&ffd,1);

	if ( result->chain != NULL && result->chain_count < result->chain_max )
	 {	int stride=nvar+1;
		memcpy(&result->chain[result->chain_count*stride],bvector,nvar*sizeof(double));
		result->chain[result->chain_count*stride+nvar]=chi2;
		result->chain_count++;
	 }
	if ( callback != NULL )
		callback(bvector,nvar,chi2,callback_data);

	if ( griddata != NULL )
		griddata[igrid].chi2=chi2;	

	for ( i=0 ; i<ndim ; i++ )
	 {	axes[i].curr++;
		j=axes[i].varindex;
		bvector[j]=axes[i].imin+axes[i].istp*axes[i].curr;
		if ( axes[i].curr < axes[i].n )
			break;
		else
		 {	bvector[j]=axes[i].imin;
			axes[i].curr=0;
		 }
	 }
	if ( i >= ndim )
		break;

	igrid++;
  };

 if ( griddata != NULL )
  {	int	k,s,v;
	int	*mins,nmin;
	int	bdist,b;
	double	**mmatrix,*mvector,mscalar;

	map_chi2_downlink(axes,ndim,griddata);
	mins=NULL;
	map_chi2_partitions(griddata,ngrid,&mins,&nmin);

	mmatrix=matrix_alloc(ndim);
	mvector=vector_alloc(ndim);

	for ( n=0 ; n<nmin ; n++ )
	 {	k=mins[n];
		bdist=-1;
		for ( i=0 ; i<ndim ; i++ )
		 {	s=k%axes[i].n;
			b=axes[i].n-1-s;
			if ( s<b )	b=s;
			if ( bdist<0 || b<bdist )	bdist=b;
			v=axes[i].varindex;
			bvector[v]=axes[i].imin+s*axes[i].istp;
			k/=axes[i].n;
		 }
		if ( ! ( bdist>0 ) )
			continue;

		map_chi2_fit_minimum(axes,ndim,griddata,mins[n],3,
			mmatrix,mvector,&mscalar);
		for ( i=0 ; i<ndim ; i++ )
		 {	v=axes[i].varindex;
			bvector[v]=mvector[i];
		 }
		for ( i=0 ; i<nvar ; i++ )
		 {	result->params[i]=bvector[i];	}
		result->chi2=griddata[mins[n]].chi2;
		break;
	 }

	vector_free(mvector);
	matrix_free(mmatrix);

	if ( mins != NULL )	free(mins);
  }

 result->nrow=nrow;
 result->error_code=0;

 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
 	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(bvector);

 if ( griddata != NULL )	free(griddata);
 if ( axes != NULL )		free(axes);

 return(0);
 
}

int fit_error_monte_carlo_estimation(lfitdata *lf,lfit_result *result,variable *vars,int nvar,
	lfit_input_arrays *input,lfit_step_callback callback,void *callback_data)
{
 int		i,n,nc,nrc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fitvary,*fiteval,*fitwght;
 double		*wvars,mchi,nchi;
 fitinputrow	*fitinputrows;
 int		nrow,downlink_max_iter,max_decrease_rejection;
 double		decrease_error_ratio;
 double		*bvector,*dvector,*evector,*mvector,**ematrix,*xvector;
 fitfunctdata	ffd;
 fitparam	*fp;
 int		use_dft_rednoise=0;
 complex	*resdft;
 constraint	ccstat,*cc;

 fp=&lf->parameters;

 nc=lf->fconstraint.nc;
 nrc=nc;
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 if ( lf->fconstraint.nc>0 )
  {	ccstat.nc=     lf->fconstraint.nc;
	ccstat.cmatrix=lf->fconstraint.cmatrix;
	ccstat.cvector=lf->fconstraint.cvector;
	cc=&ccstat;
  }
 else
	cc=NULL;

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);
 dvector=vector_alloc(nvar);
 evector=vector_alloc(nvar);
 mvector=vector_alloc(nvar);
 xvector=vector_alloc(nvar);

 load_all_fit_data_from_arrays(lf,nvar,input,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	snprintf(result->error_msg,sizeof(result->error_msg),"too few lines for fitting");
	result->error_code=-1;
	return(-1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitvary=(double *)malloc(nrow*sizeof(double));
 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( i=0 ; i<nvar ; i++ )
  {	mvector[i]=vars[i].init;
  }

 downlink_max_iter=1000;
 max_decrease_rejection=10; 
 decrease_error_ratio=0.5;

 for ( i=0 ; i<nvar ; i++ )
  {	evector[i]=vars[i].ierr;
	xvector[i]=vars[i].diff;
  }

 if ( fp->tune.emce.sub_method==FIT_METHOD_DHSX )
  {	ematrix=matrix_alloc(nvar);
	get_error_matrix(ematrix,nvar,evector,NULL,lf->fconstraint.nc,lf->fconstraint.cproject);
  }
 else
	ematrix=NULL;

 if ( fp->tune.emce.skip_initial_fit )
	mchi=-1;
 else if ( fp->tune.emce.sub_method==FIT_METHOD_CLLS )
  {	lin_fit_con(fitpnt,fitdep,mvector,fitwght,
	fit_function,nvar,nrow,&ffd,cc,NULL);
	mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,&ffd,1);
  }
 else if ( fp->tune.emce.sub_method==FIT_METHOD_NLLM )
  {	double	lam;
	lam=fp->tune.nllm.lambda;
	for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
	 {	lam=nlm_fit_base_con(fitpnt,fitdep,mvector,fitwght,
		fit_function,nvar,nrow,&ffd,cc,
		lam,fp->tune.nllm.lambda_mpy);
	 }
	mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,&ffd,1);
  }
 else if ( fp->tune.emce.sub_method==FIT_METHOD_LMND )
  {	double	lam;
	lam=fp->tune.nllm.lambda;
	for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
	 {	lam=nlm_fit_nmdf_con(fitpnt,fitdep,mvector,fitwght,
		fit_function,nvar,nrow,&ffd,cc,
		lam,fp->tune.nllm.lambda_mpy,xvector);
	 }
	mchi=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,mvector,&ffd,1);
  }
 else if ( fp->tune.emce.sub_method==FIT_METHOD_MCMC )
  {	mchi=find_monte_carlo_minimum(fitpnt,fitdep,fitwght,nrow,mvector,&ffd,
		evector,downlink_max_iter,max_decrease_rejection,decrease_error_ratio,
		lf,vars,nvar);
  }
 else /* FIT_METHOD_DHSX: default */
  {	mchi=find_downhill_simplex_minimum(fitpnt,fitdep,fitwght,nrow,
		mvector,&ffd,ematrix,nvar,nvar-nc,NULL,0);
  }

 data_get_scatters(fitpnt,fitdep,fitwght,fiteval,
	nrow,mvector,&ffd,0);

 if ( use_dft_rednoise )
	resdft=data_get_noise_dft(fitpnt,fitdep,fitwght,fiteval,
		nrow,mvector,&ffd);
 else
	resdft=NULL;
 
 for ( i=0 ; i<nvar ; i++ )
  {	result->params[i]=mvector[i];	}
 result->chi2=mchi;
 result->chain_count=0;

 if ( result->chain != NULL && result->chain_count < result->chain_max )
  {	int stride=nvar+1;
	memcpy(&result->chain[result->chain_count*stride],mvector,nvar*sizeof(double));
	result->chain[result->chain_count*stride+nvar]=(mchi>=0.0?mchi:0.0);
	result->chain_count++;
  }
 if ( callback != NULL )
	callback(mvector,nvar,(mchi>=0.0?mchi:0.0),callback_data);

 for ( n=0 ; n<fp->mc_iterations ; n++ )
  {	
	if ( resdft != NULL )
		data_perturb_dft(fitpnt,fiteval,fitvary,nrow,&ffd,resdft);
	else
		data_perturb_wtn(fitpnt,fiteval,fitvary,nrow,&ffd);

	for ( i=0 ; i<nvar ; i++ )
	 {	evector[i]=vars[i].ierr;
		bvector[i]=mvector[i];
	 }

	if ( fp->tune.emce.sub_method==FIT_METHOD_CLLS )
	 {	lin_fit_con(fitpnt,fitvary,bvector,fitwght,fit_function,nvar,nrow,&ffd,cc,NULL);
		nchi=data_get_chisquare(fitpnt,fitvary,fitwght,NULL,nrow,bvector,&ffd,1);
	 }
	else if ( fp->tune.emce.sub_method==FIT_METHOD_NLLM )
	 {	double	lam;
		lam=fp->tune.nllm.lambda;
		for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
		 {	lam=nlm_fit_base_con(fitpnt,fitvary,bvector,fitwght,
			fit_function,nvar,nrow,&ffd,cc,lam,fp->tune.nllm.lambda_mpy);
		 }
		nchi=data_get_chisquare(fitpnt,fitvary,fitwght,NULL,nrow,bvector,&ffd,1);
	 }
	else if ( fp->tune.emce.sub_method==FIT_METHOD_LMND )
	 {	double	lam;
		lam=fp->tune.nllm.lambda;
		for ( i=0 ; i<fp->tune.nllm.max_iter ; i++ )
		 {	lam=nlm_fit_nmdf_con(fitpnt,fitvary,bvector,fitwght,
			fit_function,nvar,nrow,&ffd,cc,
			lam,fp->tune.nllm.lambda_mpy,xvector);
		 }
		nchi=data_get_chisquare(fitpnt,fitvary,fitwght,NULL,nrow,bvector,&ffd,1);
	 }
	else if ( fp->tune.emce.sub_method==FIT_METHOD_MCMC )
	 {	nchi=find_monte_carlo_minimum(fitpnt,fitvary,fitwght,nrow,bvector,&ffd,
			evector,downlink_max_iter,max_decrease_rejection,decrease_error_ratio,
			lf,vars,nvar);
	 }	
	else /* FIT_METHOD_DHSX: default */
	 {	nchi=find_downhill_simplex_minimum(fitpnt,fitvary,fitwght,nrow,
			bvector,&ffd,ematrix,nvar,nvar-nc,NULL,0);
	 }

	if ( result->chain != NULL && result->chain_count < result->chain_max )
	 {	int stride=nvar+1;
		memcpy(&result->chain[result->chain_count*stride],bvector,nvar*sizeof(double));
		result->chain[result->chain_count*stride+nvar]=nchi;
		result->chain_count++;
	 }
	if ( callback != NULL )
		callback(bvector,nvar,nchi,callback_data);
  }

 result->nrow=nrow;
 result->error_code=0;

 if ( resdft != NULL )		free(resdft);
	
 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitvary  != NULL )	free(fitvary);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(xvector);
 vector_free(mvector);
 if ( ematrix != NULL )		free(ematrix);
 vector_free(evector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}

int fit_extended_markov_chain_mc(lfitdata *lf,lfit_result *result,variable *vars,int nvar,
	lfit_input_arrays *input,lfit_step_callback callback,void *callback_data)
{
 int		i,j,n,nc,nrc,ccnt,is_not_accepted;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fiteval,*fitwght;
 double		*wvars,chi2,nchi,mchi,vx,u,alpha,dalpha;
 fitinputrow	*fitinputrows;
 int		nrow,window,corr_c,iiter,niter;
 double		*bvector,*dvector,*mvector,*gvector,*tvector;
 double		**rcmatrix,**cmatrix;
 double		**dcmatrix,*dcvector,*ddvector,*dsvector;
 double		**corrlendd,**corr_ftmp,**corr_prev,**corr_carr,*corr_sum1,*corr_sum2;
 fitfunctdata	ffd;
 fitparam	*fp;
 int		nseparated,*idx;
 fitconstraint	ifc;

 fp=&lf->parameters;

 nc=lf->fconstraint.nc;
 nrc=nc;
 nseparated=0;
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
	if ( vars[i].is_separated )
		nseparated++;
  }
 
 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);
 dvector=vector_alloc(nvar);
 mvector=vector_alloc(nvar);
 gvector=vector_alloc(nvar);
 tvector=vector_alloc(nvar);

 dcmatrix=matrix_alloc(nvar);
 dcvector=vector_alloc(nvar);
 ddvector=vector_alloc(nvar);
 dsvector=vector_alloc(nvar);

 window=lf->parameters.tune.xmmc.window;
 corrlendd=matrix_alloc_gen(window,3*nvar);
 corr_ftmp=corrlendd+0*nvar;
 corr_prev=corrlendd+1*nvar;
 corr_carr=corrlendd+2*nvar;
 corr_sum1=vector_alloc(nvar);
 corr_sum2=vector_alloc(nvar);

 load_all_fit_data_from_arrays(lf,nvar,input,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	snprintf(result->error_msg,sizeof(result->error_msg),"too few lines for fitting");
	result->error_code=-1;
	return(-1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 fiteval=(double *)malloc(sizeof(double)*nrow);

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;			}

 nseparated=get_separate_index(vars,nvar,&idx);

 result->chain_count=0;

 ifc.nc=nc+nseparated;
 ifc.cmatrix=matrix_alloc(nvar);
 for ( i=0 ; i<nc ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	ifc.cmatrix[i][j]=lf->fconstraint.cmatrix[i][j];		}
  }
 for ( i=0 ; i<nseparated ; i++ )
  {	for ( j=0 ; j<nvar ; j++ )
 	 {	ifc.cmatrix[nc+i][j]=0.0;		}
	ifc.cmatrix[nc+i][idx[i]]=1.0;
  }
 ifc.cproject=matrix_alloc(nvar);
 constraint_initialize_proj_matrix(ifc.nc,nvar,ifc.cmatrix,ifc.cproject);

 if ( ! lf->parameters.tune.xmmc.skip_initial_fit )
  {	double	**ematrix;
	double	**icmatrix;

	ematrix=matrix_alloc(nvar);
	icmatrix=data_get_inversecovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd);
	get_error_matrix(ematrix,nvar,NULL,icmatrix,ifc.nc,ifc.cproject);
	matrix_free(icmatrix);
	find_downhill_simplex_minimum(fitpnt,fitdep,fitwght,nrow,bvector,&ffd,ematrix,nvar,nvar-ifc.nc,idx,nseparated);
	matrix_free(ematrix);
  }

 if ( ! lf->parameters.tune.xmmc.is_adaptive )
	rcmatrix=data_get_rootcovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd,ifc.nc,ifc.cproject);
 else
	rcmatrix=NULL;

 if ( nc > 0 )
  {	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].flags & VAR_IS_CONSTANT )
	 	 {	bvector[i]=vars[i].init;		}
	 }
	if ( nrc>0 )
	 {	constraint_force(bvector,nvar,nc,
		lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
	 }
  }

 chi2=data_get_chisquare_separated_linear(fitpnt,fitdep,fitwght,fiteval,nrow,bvector,&ffd,0,nvar,idx,nseparated,bvector);
 mchi=chi2;
 for ( i=0 ; i<nvar ; i++ )
  {	 mvector[i]=bvector[i];
	ddvector[i]=bvector[i];
  }

 if ( ! lf->parameters.tune.xmmc.skip_initial_fit )
  {	if ( result->chain != NULL && result->chain_count < result->chain_max )
	 {	int stride=nvar+1;
		memcpy(&result->chain[result->chain_count*stride],bvector,nvar*sizeof(double));
		result->chain[result->chain_count*stride+nvar]=chi2;
		result->chain_count++;
	 }
	if ( callback != NULL )
		callback(bvector,nvar,chi2,callback_data);
  }

 niter=(lf->parameters.tune.xmmc.niter>0?lf->parameters.tune.xmmc.niter:0);

 for ( iiter=0 ; iiter <= niter ; iiter++ )
  {

	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	dcmatrix[i][j]=0.0;		}
		dsvector[i]=0.0;
	 }

	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<window ; j++ )
		 {	corr_ftmp[i][j]=0.0;
			corr_prev[i][j]=0.0;
			corr_carr[i][j]=0.0;
		 }
		corr_sum1[i]=0.0;
		corr_sum2[i]=0.0;
	 }
	corr_c=0;

	for ( n=0,ccnt=0 ; (lf->parameters.tune.xmmc.count_accepted?ccnt:n)<fp->mc_iterations ; n++ )
	 {	
		is_not_accepted=0;
		for ( j=0 ; j<nvar ; j++ )
		 {	gvector[j]=lfit_get_gaussian(0.0,1.0);		}

		if ( lf->parameters.tune.xmmc.is_adaptive )
			rcmatrix=data_get_rootcovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd,ifc.nc,ifc.cproject);

		for ( i=0 ; i<nvar ; i++ )
		 {	dvector[i]=bvector[i];
			for ( j=0 ; j<nvar ; j++ )
			 {	dvector[i]+=rcmatrix[i][j]*gvector[j];		}
		 }

		if ( lf->parameters.tune.xmmc.is_adaptive )
		 {	matrix_free(rcmatrix);
			rcmatrix=NULL;
		 }
	
		for ( i=0 ; i<nvar ; i++ )
		 {
			if ( (vars[i].flags & VAR_MIN) && dvector[i] < vars[i].imin )
			 {	is_not_accepted=1;
				dvector[i]=vars[i].imin;
			 }
			if ( (vars[i].flags & VAR_MAX) && dvector[i] > vars[i].imax )
			 {	is_not_accepted=1;
				dvector[i]=vars[i].imax;
			 }
		 }

		for ( i=0 ; i<lf->dconstraint.ndcfunct && (!is_not_accepted) ; i++ )
		 {	double	w;
			lfit_psn_double_calc(lf->dconstraint.dcfuncts[i].funct,NULL,lpg.pl_funct,&w,dvector);
			if ( w<=0.0 )
				is_not_accepted=1;
		 }
	
		if ( is_not_accepted )
		 {	continue;
		 }

		if ( nc > 0 )
		 {	for ( i=0 ; i<nvar ; i++ )
			 {	if ( vars[i].flags & VAR_IS_CONSTANT )
			 	 {	dvector[i]=vars[i].init;		}
			 }
			if ( nrc>0 )
			 {	constraint_force(dvector,nvar,nc,
				lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
			 }
		 }

		nchi=data_get_chisquare_separated_linear(fitpnt,fitdep,fitwght,fiteval,nrow,dvector,&ffd,0,nvar,idx,nseparated,dvector);

		if ( nchi<mchi )
		 {	mchi=nchi;
			for ( i=0 ; i<nvar ; i++ )
			 {	mvector[i]=dvector[i];		}
		 }
	
		vx=exp(-0.5*(nchi-chi2));
		alpha=(vx<1.0?vx:1.0);

		if ( alpha>0.0 && (u=random_double()) <= alpha )
		 {	for ( i=0 ; i<nvar ; i++ )
			 {	bvector[i]=dvector[i];
				dcvector[i]=bvector[i]-ddvector[i];
			 }
			chi2=nchi;

			for ( i=0 ; i<nvar ; i++ )
			 {	double	d;
				d=dcvector[i];
				if ( ccnt<window )
				 {	corr_ftmp[i][ccnt]=d;		}
				corr_prev[i][corr_c]=d;
				for ( j=0 ; j<window && j<=ccnt ; j++ )
				 {	corr_carr[i][j]+=d*corr_prev[i][(corr_c+window-j)%window];	}
				corr_sum1[i]+=d;
				corr_sum2[i]+=d*d;
			 }
			corr_c=(corr_c+1)%window;

			ccnt++;

			if ( result->chain != NULL && result->chain_count < result->chain_max )
			 {	int stride=nvar+1;
				memcpy(&result->chain[result->chain_count*stride],bvector,nvar*sizeof(double));
				result->chain[result->chain_count*stride+nvar]=chi2;
				result->chain_count++;
			 }
			if ( callback != NULL )
				callback(bvector,nvar,chi2,callback_data);

			for ( i=0 ; i<nvar ; i++ )
			 {	for ( j=0 ; j<nvar ; j++ )
				 {	dcmatrix[i][j]+=dcvector[i]*dcvector[j];	}
				dsvector[i]+=dcvector[i];
			 }
		 }
		else
			is_not_accepted=1;	
	 }

	for ( i=0 ; i<nvar ; i++ )
	 {	double	s1,s2,sig2;
		int	j,n;
	
	 	if ( ccnt>0 )
		 {	corr_sum1[i]/=(double)ccnt;
			corr_sum2[i]/=(double)ccnt;
		 }
		s1=corr_sum1[i];
		s2=corr_sum2[i];

		sig2=s2-s1*s1;
		if ( sig2<0.0 )	sig2=0.0;
		
		for ( n=1 ; n<window ; n++ )
		 {	for ( j=0 ; j<n ; j++ )
			 {	corr_carr[i][n]+=corr_ftmp[i][j]*corr_prev[i][(corr_c+window-n)%window];	}
		 }
		for ( n=0 ; n<window ; n++ )
		 {	if ( ccnt>0 && sig2>0.0 )
				corr_carr[i][n]=(corr_carr[i][n]-s1*s1*ccnt)/(ccnt*sig2);
			else
				corr_carr[i][n]=0.0;
		 }
	 }

	for ( i=0 ; i<nvar && ccnt>0 ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	dcmatrix[i][j]/=(double)ccnt;		}
		dsvector[i]/=(double)ccnt;
	 }	
	for ( i=0 ; i<nvar ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	dcmatrix[i][j] -= dsvector[i]*dsvector[j];	}
	 }

	if ( n>0 )
	 {	alpha=(double)ccnt/(double)n;
		dalpha=sqrt((double)ccnt)/(double)n;
	 }
	else
	 {	alpha=0.0;
		dalpha=0.0;
	 }

	mchi=data_get_chisquare_separated_linear(fitpnt,fitdep,fitwght,fiteval,nrow,mvector,&ffd,1,nvar,idx,nseparated,mvector);

	cmatrix=data_get_covariance(fitpnt,fitdep,fitwght,nrow,nvar,mvector,&ffd);
	if ( cmatrix != NULL )
	 {	if ( result->cov_matrix != NULL )
		 {	for ( i=0 ; i<nvar ; i++ )
			 {	for ( j=0 ; j<nvar ; j++ )
				 {	result->cov_matrix[i][j]=cmatrix[i][j];		}
			 }
		 }
		matrix_free(cmatrix);
	 }

	if ( iiter<niter && rcmatrix != NULL )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	for ( j=0 ; j<nvar ; j++ )
			 {	rcmatrix[i][j]=dcmatrix[i][j];		}
		 }
		matrix_root(rcmatrix,nvar,NULL);
	 }

  }

 result->chi2=mchi;
 result->acceptance=alpha;
 for ( i=0 ; i<nvar ; i++ )
  {	result->params[i]=mvector[i];	}
 result->nrow=nrow;
 result->error_code=0;

 if ( ! lf->parameters.tune.xmmc.is_adaptive )
	matrix_free(rcmatrix);

 if ( fiteval != NULL )		free(fiteval);
 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(corr_sum2);
 vector_free(corr_sum1);
 matrix_free(corrlendd);

 if ( ifc.cmatrix  != NULL )	matrix_free(ifc.cmatrix);
 if ( ifc.cproject != NULL )	matrix_free(ifc.cproject);

 vector_free(dsvector);
 vector_free(ddvector);
 vector_free(dcvector);
 matrix_free(dcmatrix);

 vector_free(tvector);
 vector_free(gvector);
 vector_free(mvector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}

int fit_fisher_matrix_analysis(lfitdata *lf,lfit_result *result,
	variable *vars,int nvar,dvariable *dvars,int ndvar,
	lfit_input_arrays *input,lfit_step_callback callback,void *callback_data)
{
 int		i,j,n,nc,nrc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fitwght;
 double		*wvars,chi2;
 fitinputrow	*fitinputrows;
 int		nrow;
 double		*bvector,*dvector,*gvector;
 double		**rcmatrix;
 fitfunctdata	ffd;
 fitparam	*fp;

 if ( ndvar<=0 || dvars==NULL )
	ndvar=0,dvars=NULL;

 fp=&lf->parameters;

 nc=lf->fconstraint.nc;
 nrc=nc;
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar+ndvar);
 dvector=vector_alloc(nvar);
 gvector=vector_alloc(nvar);
 
 load_all_fit_data_from_arrays(lf,nvar,input,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;			}
 for ( i=0 ; i<ndvar ; i++ )
  {	double	d;
	lfit_psn_double_calc(dvars[i].funct,NULL,lpg.pl_funct,&d,bvector);
	bvector[i+nvar]=d;
  }

 data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,bvector,&ffd,1);
 data_get_scatters(fitpnt,fitdep,fitwght,NULL,nrow,bvector,&ffd,0);

 for ( i=0 ; i<nvar ; i++ )
  {	result->params[i]=bvector[i];	}
 result->chain_count=0;

 if ( result->chain != NULL && result->chain_count < result->chain_max )
  {	int stride=nvar+1;
	memcpy(&result->chain[result->chain_count*stride],bvector,nvar*sizeof(double));
	result->chain[result->chain_count*stride+nvar]=0.0;
	result->chain_count++;
  }
 if ( callback != NULL )
	callback(bvector,nvar,0.0,callback_data);

 rcmatrix=data_get_rootcovariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd,lf->fconstraint.nc,lf->fconstraint.cproject);

 for ( n=0 ; n<fp->mc_iterations && fp->tune.fima.do_montecarlo ; n++ )
  {	for ( j=0 ; j<nvar ; j++ )
	 {	gvector[j]=lfit_get_gaussian(0.0,1.0);		}
	for ( i=0 ; i<nvar ; i++ )
	 {	dvector[i]=bvector[i];
		for ( j=0 ; j<nvar ; j++ )
		 {	dvector[i]+=rcmatrix[i][j]*gvector[j];		}
	 }
	if ( nc > 0 )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	if ( vars[i].flags & VAR_IS_CONSTANT )
		 	 {	dvector[i]=vars[i].init;		}
		 }
		if ( nrc>0 )
		 {	constraint_force(dvector,nvar,nc,
			lf->fconstraint.invlmatrix,lf->fconstraint.cvector);
		 }
	 }
	chi2=data_get_chisquare(fitpnt,fitdep,fitwght,NULL,nrow,dvector,&ffd,1);
	if ( result->chain != NULL && result->chain_count < result->chain_max )
	 {	int stride=nvar+1;
		memcpy(&result->chain[result->chain_count*stride],dvector,nvar*sizeof(double));
		result->chain[result->chain_count*stride+nvar]=chi2;
		result->chain_count++;
	 }
	if ( callback != NULL )
		callback(dvector,nvar,chi2,callback_data);
  }

 matrix_free(rcmatrix);

 {	double	**cmatrix,**ocmatrix,*tvector;

	tvector=vector_alloc(nvar+ndvar);
	cmatrix=matrix_alloc(nvar+ndvar);
	ocmatrix=data_get_covariance(fitpnt,fitdep,fitwght,nrow,nvar,bvector,&ffd);

	if ( ocmatrix != NULL )
	 {	int	k,l,m,n;
		double	c,p1,p2,**diff;

		diff=matrix_alloc(nvar+ndvar);
		for ( k=0 ; k<ndvar ; k++ )
		 {	for ( m=0 ; m<nvar ; m++ )
			 {	lfit_psn_double_calc(dvars[k].diff[m],NULL,lpg.pl_funct,&p1,bvector);
				diff[k][m]=p1;
			 }
		 }

		for ( k=0 ; k<nvar+ndvar ; k++ )
		 {  for ( l=0 ; l<=k ; l++ )
		     {	c=0.0;
			for ( m=0 ; m<nvar ; m++ )
			 {  for ( n=0 ; n<nvar ; n++ )
			     {	if ( k<nvar && k==m )	p1=1.0;
				else if ( k<nvar )	p1=0.0;
				else			p1=diff[k-nvar][m];
				if ( l<nvar && l==n )	p2=1.0;
				else if ( l<nvar )	p2=0.0;
				else			p2=diff[l-nvar][n];
				c+=p1*ocmatrix[m][n]*p2;
			     }
			 }
			cmatrix[k][l]=c;
		     }
		 }
		for ( k=0 ; k<nvar+ndvar ; k++ )
		 {  for ( l=k+1 ; l<nvar+ndvar ; l++ )
		     {	cmatrix[k][l]=cmatrix[l][k];	}
		 }

		matrix_free(diff);
		matrix_free(ocmatrix);
	 }
	else
		cmatrix=NULL;

	if ( cmatrix != NULL )
	 {	for ( i=0 ; i<nvar ; i++ )
		 {	tvector[i]=sqrt(cmatrix[i][i]);
			result->errors[i]=tvector[i];
		 }
		if ( result->cov_matrix != NULL )
		 {	for ( i=0 ; i<nvar ; i++ )
			 {	for ( j=0 ; j<nvar ; j++ )
				 {	result->cov_matrix[i][j]=cmatrix[i][j];	}
			 }
		 }
		matrix_free(cmatrix);
	 }
	vector_free(tvector);
 }

 result->nrow=nrow;
 result->error_code=0;

 if ( fitpnt  != NULL )		free(fitpnt);
 if ( fitinputrows != NULL )
  {	for ( i=0 ; i<nrow ; i++ )
	 {	free(fitinputrows[i].line);	}
	free(fitinputrows);
  }
 if ( fitwght  != NULL )	free(fitwght);
 if ( fitdep   != NULL )	free(fitdep);
 if ( fitdata  != NULL )	free(fitdata);

 vector_free(gvector);
 vector_free(dvector);
 vector_free(bvector);

 return(0);
}

int fit_downhill_simplex(lfitdata *lf,lfit_result *result,variable *vars,int nvar,
	lfit_input_arrays *input)
{
 int		i,nc,nrc;

 void		**fitpnt;
 double		*fitdata,*fitdep,*fitwght;
 double		*wvars,chi2;
 fitinputrow	*fitinputrows;
 int		nrow;
 double		*bvector,*evector;
 double		**ematrix;
 fitfunctdata	ffd;
 int		*idx,nseparated;

 nc=lf->fconstraint.nc;
 nrc=nc;
 for ( i=0 ; i<nvar ; i++ )	
  {	if ( vars[i].flags & VAR_IS_CONSTANT )
		nrc--;
  }

 fitdata =NULL;
 fitdep  =NULL;
 fitwght =NULL;
 fitinputrows=NULL;
 nrow=0;

 bvector=vector_alloc(nvar);
 evector=vector_alloc(nvar);
 ematrix=matrix_alloc(nvar);

 load_all_fit_data_from_arrays(lf,nvar,input,&fitdata,&fitdep,&fitwght,&fitinputrows,&nrow);

 if ( nrow < nvar-nc )	
  {	fprint_error("too few lines for fitting");
	snprintf(result->error_msg,sizeof(result->error_msg),"too few lines for fitting");
	result->error_code=-1;
	return(-1);
  }

 wvars=(double *)malloc(sizeof(double)*(nvar+lf->maxncol));

 fitpnt=(void **)malloc(nrow*sizeof(void *));
 for ( i=0 ; i<nrow ; i++ )
  {	fitinputrows[i].x=&fitdata[i*lf->maxncol];
	fitpnt[i]=(void *)(&fitinputrows[i]);
  }

 ffd.nvar=nvar;
 ffd.wvars=wvars;
 ffd.functs=lpg.pl_funct;
 ffd.lf=lf;

 for ( i=0 ; i<nvar ; i++ )
  {	bvector[i]=vars[i].init;
	evector[i]=vars[i].ierr;
  }

 ematrix=matrix_alloc(nvar);

 nseparated=get_separate_index(vars,nvar,&idx);

 if ( idx != NULL && 0<nseparated )
  {	fitconstraint	ifc;
	int		i,j;

	ifc.nc=nc+nseparated;
	ifc.cmatrix=matrix_alloc(nvar);
	for ( i=0 ; i<nc ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	ifc.cmatrix[i][j]=lf->fconstraint.cmatrix[i][j];		}
	 }
	for ( i=0 ; i<nseparated ; i++ )
	 {	for ( j=0 ; j<nvar ; j++ )
		 {	ifc.cmatrix[nc+i][j]=0.0;		}
		ifc.cmatrix[nc+i][idx[i]]=1.0;
	 }
	ifc.cproject=matrix_alloc(nvar);
	constraint_initialize_proj_matrix(ifc.nc,nvar,ifc.cmatrix,ifc.cproject);
 	get_error_matrix(ematrix,nvar,evector,NULL,ifc.nc,ifc.cproject);
	matrix_free(ifc.cproject);
	matrix_free(ifc.cmatrix);
  }
 else
 	get_error_matrix(ematrix,nvar,evector,NULL,lf->fconstraint.nc,lf->fconstraint.cproject);

 chi2=find_downhill_simplex_minimum(fitpnt,fitdep,fitwght,nrow,bvector,&ffd,ematrix,nvar,nvar-nc-nseparated,NULL,0);

 result->chi2=chi2;
 for ( i=0 ; i<nvar ; i++ )
  {	result->params[i]=bvector[i];	}
 result->nrow=nrow;
 result->error_code=0;

 if ( idx != NULL )	free(idx);

 matrix_free(ematrix);
 vector_free(evector);
 vector_free(bvector);

 return(0);
}

/*****************************************************************************/
/* lfit_python_apply: Python-callable entry point                            */
/*****************************************************************************/

static int lfit_builtins_registered = 0;

static void lfit_free_psn(void)
{
 if ( lpg.pl_sym )      { free(lpg.pl_sym);      }
 if ( lpg.pl_prop )     { free(lpg.pl_prop);     }
 if ( lpg.pl_funct )    { free(lpg.pl_funct);    }
 if ( lpg.pl_diff )     { free(lpg.pl_diff);     }
 if ( lpg.pl_symeval )  { free(lpg.pl_symeval);  }
 if ( lpg.pl_macro )    { free(lpg.pl_macro);    }
 if ( lpg.symbols )     { free(lpg.symbols);     }
 if ( lpg.lffregs )     { free(lpg.lffregs);     }
 memset(&lpg,0,sizeof(lfitpsnglobal));
 lpg.pl_simp=psn_lfit_simp;
 lfit_builtins_registered=0;
}

static int lfit_ensure_builtins(void)
{
 int	i,r;
 if ( lfit_builtins_registered )
	return(0);
 for ( i=0 ; psnlfit_list_builtin_normal_operators[i].name != NULL ; i++ )
  {	r=lfit_register_internal(&psnlfit_list_builtin_normal_operators[i]);
	if ( r ) return(-1);
  }
 for ( i=0 ; psnlfit_list_builtin_elementary_functions[i].name != NULL ; i++ )
  {	r=lfit_register_internal(&psnlfit_list_builtin_elementary_functions[i]);
	if ( r ) return(-1);
  }
 for ( i=0 ; psnlfit_list_builtin_interpolators[i].name != NULL ; i++ )
  {	r=lfit_register_internal(&psnlfit_list_builtin_interpolators[i]);
	if ( r ) return(-1);
  }
 for ( i=0 ; psnlfit_list_builtin_aa_functions[i].name != NULL ; i++ )
  {	r=lfit_register_internal(&psnlfit_list_builtin_aa_functions[i]);
	if ( r ) return(-1);
  }
 lfit_builtins_registered=1;
 return(0);
}

int lfit_python_apply(
	double *array_data,int nrow_in,int ncol_in,
	char *variables,
	char *columns,
	char *function,
	char *dependent,
	char *error,
	char *weight,
	int fit_method,
	char *parameters,
	char *differences,
	char *separate,
	char *perturbations,
	int seed,
	int mc_iterations,
	double rejection_level,
	int rejection_niter,
	int weighted_sigma,
	char **macros,
	char *constraints,
	int errdump,
	char *format,
	char *correlation_format,
	char *derived_variables,
	int is_dump_delta,
	int resdump,
	int force_nonlinear,
	char *columns_output,
	lfit_result *result,
	double *chain,int chain_max,
	lfit_step_callback callback,void *callback_data)
{
 int		i,j,r,calc_derivatives;
 variable	*vars;
 int		nvar;
 dvariable	*dvars;
 int		ndvar;
 column		*cols;
 int		ncol;
 int		errtype;
 char		*errstr;
 psnsym		*varsym,*colsym,*mysyms[8];
 lfitdata	lf_static,*lf=&lf_static;
 datablock	db_static,*db=&db_static;
 lfit_input_arrays	input;
 double		*array_ptrs[1];
 int		array_nrows[1],array_ncols[1];
 int		is_fit;

 memset(lf,0,sizeof(lfitdata));
 memset(db,0,sizeof(datablock));
 memset(result,0,sizeof(lfit_result));

 lfit_free_psn();
 lfit_ensure_builtins();

 if ( error != NULL && weight != NULL )
  {	snprintf(result->error_msg,sizeof(result->error_msg),"both error and weight specified, use only one");
	result->error_code=-1;
	return(-1);
  }
 if ( error != NULL )
	errtype=0, errstr=error;
 else if ( weight != NULL )
	errtype=1, errstr=weight;
 else
	errtype=0, errstr=NULL;

 char *s_col  = columns    ? strdup(columns)    : NULL;
 char *s_var  = variables  ? strdup(variables)  : NULL;
 char *s_func = function   ? strdup(function)   : NULL;
 char *s_dep  = dependent  ? strdup(dependent)  : NULL;
 char *s_err  = errstr     ? strdup(errstr)     : NULL;
 char *s_diff = differences? strdup(differences): NULL;
 char *s_sep  = separate   ? strdup(separate)   : NULL;
 char *s_pert = perturbations ? strdup(perturbations) : NULL;
 char *s_cnt  = constraints? strdup(constraints): NULL;
 char *s_fmt  = format     ? strdup(format)     : NULL;
 char *s_cfmt = correlation_format ? strdup(correlation_format) : NULL;
 char *s_dvr  = derived_variables ? strdup(derived_variables) : NULL;
 char *s_par  = parameters ? strdup(parameters) : NULL;
 char *s_co   = columns_output ? strdup(columns_output) : NULL;

 #define FREE_ALL_STRDUP do { \
	free(s_col); free(s_var); free(s_func); free(s_dep); free(s_err); \
	free(s_diff); free(s_sep); free(s_pert); free(s_cnt); free(s_fmt); \
	free(s_cfmt); free(s_dvr); free(s_par); free(s_co); } while(0)

 #define RETURN_ERROR(msg_fmt,...) do { \
	snprintf(result->error_msg,sizeof(result->error_msg),msg_fmt,##__VA_ARGS__); \
	result->error_code=-1; FREE_ALL_STRDUP; return(-1); } while(0)

 if ( lfit_ensure_builtins() )
	RETURN_ERROR("failed to register builtins");

 if ( macros != NULL )
  {	char	**xdef;
	int	omajor=0;
	for ( xdef=macros ; *xdef ; xdef++ )
	 {	lfit_define_macro(lf,*xdef,&omajor);	}
  }

 is_fit = (s_dep != NULL) ? 1 : 0;

 if ( s_var != NULL )
  {	remove_spaces_and_comments(s_var);
	extract_variables(s_var,&vars,&nvar);
  }
 else
  {	vars=NULL; nvar=0;	}

 if ( vars==NULL || nvar<=0 )
  {	if ( is_fit )
		RETURN_ERROR("no variables defined");
	nvar=0;
  }

 for ( i=0 ; i<nvar ; i++ )
  {	lfit_register_symbol_add(&lpg,vars[i].name);	}

 if ( s_dvr != NULL )
  {	remove_spaces_and_comments(s_dvr);
	extract_derived_variables(s_dvr,&dvars,&ndvar);
  }
 else
  {	dvars=NULL; ndvar=0;	}
 for ( i=0 ; dvars != NULL && i<ndvar ; i++ )
  {	lfit_register_symbol_add(&lpg,dvars[i].name);	}

 if ( s_col == NULL )
	RETURN_ERROR("column definitions missing");
 if ( extract_columns(s_col,&cols,&ncol) )
	RETURN_ERROR("invalid column specification");

 db->key="default";
 db->cols=cols;
 db->ncol=ncol;
 db->colarg=s_col;
 db->fncarg=s_func;
 db->deparg=s_dep;
 db->errarg=s_err;
 db->errtype=errtype;
 db->inparg=NULL;
 db->xnoise=0.0;

 lf->ndatablock=1;
 lf->datablocks=db;
 lf->maxncol=ncol;

 lf->parameters.fit_method=fit_method;
 lf->parameters.niter=rejection_niter;
 lf->parameters.sigma=rejection_level;
 lf->parameters.mc_iterations=mc_iterations;
 lf->parameters.tune.nllm.lambda=0.001;
 lf->parameters.tune.nllm.lambda_mpy=10.0;
 lf->parameters.tune.nllm.max_iter=10;
 lf->parameters.tune.nllm.numeric_derivs=0;
 lf->parameters.tune.dhsx.use_fisher_sx=0;
 lf->parameters.tune.mcmc.use_gibbs=0;
 lf->parameters.tune.mcmc.count_accepted=1;
 lf->parameters.tune.xmmc.skip_initial_fit=0;
 lf->parameters.tune.xmmc.is_adaptive=0;
 lf->parameters.tune.xmmc.count_accepted=1;
 lf->parameters.tune.xmmc.window=16;
 lf->parameters.tune.xmmc.niter=0;
 lf->parameters.tune.emce.skip_initial_fit=0;
 lf->parameters.tune.emce.sub_method=FIT_METHOD_DHSX;
 lf->parameters.tune.fima.do_montecarlo=0;
 lf->parameters.tune.fima.write_original=1;
 lf->parameters.tune.fima.write_uncert=1;
 lf->parameters.tune.fima.write_correl=1;

 if ( s_par != NULL )
  {	char *ptmp=strdup(s_par),*tok,*eq;
	for ( tok=strtok(ptmp,",") ; tok != NULL ; tok=strtok(NULL,",") )
	 {	eq=strchr(tok,'=');
		if ( eq ) *eq=0;
		if ( strcmp(tok,"accepted")==0 )
		 {	lf->parameters.tune.mcmc.count_accepted=1;
			lf->parameters.tune.xmmc.count_accepted=1;
		 }
		else if ( strcmp(tok,"nonaccepted")==0 )
		 {	lf->parameters.tune.mcmc.count_accepted=0;
			lf->parameters.tune.xmmc.count_accepted=0;
		 }
		else if ( strcmp(tok,"gibbs")==0 )
			lf->parameters.tune.mcmc.use_gibbs=1;
		else if ( strcmp(tok,"skip")==0 )
		 {	lf->parameters.tune.xmmc.skip_initial_fit=1;
			lf->parameters.tune.emce.skip_initial_fit=1;
		 }
		else if ( strcmp(tok,"adaptive")==0 )
			lf->parameters.tune.xmmc.is_adaptive=1;
		else if ( strcmp(tok,"fisher")==0 )
			lf->parameters.tune.dhsx.use_fisher_sx=1;
		else if ( strcmp(tok,"clls")==0 || strcmp(tok,"linear")==0 )
			lf->parameters.tune.emce.sub_method=FIT_METHOD_CLLS;
		else if ( strcmp(tok,"nllm")==0 || strcmp(tok,"nonlinear")==0 )
			lf->parameters.tune.emce.sub_method=FIT_METHOD_NLLM;
		else if ( strcmp(tok,"lmnd")==0 )
			lf->parameters.tune.emce.sub_method=FIT_METHOD_LMND;
		else if ( strcmp(tok,"dhsx")==0 || strcmp(tok,"downhill")==0 )
			lf->parameters.tune.emce.sub_method=FIT_METHOD_DHSX;
		else if ( strcmp(tok,"mc")==0 || strcmp(tok,"montecarlo")==0 )
		 {	if ( fit_method==FIT_METHOD_FIMA )
				lf->parameters.tune.fima.do_montecarlo=1;
			else
				lf->parameters.tune.emce.sub_method=FIT_METHOD_MCMC;
		 }
		else if ( eq != NULL )
		 {	double dval; int ival;
			if ( strcmp(tok,"lambda")==0 && sscanf(eq+1,"%lg",&dval)==1 )
				lf->parameters.tune.nllm.lambda=dval;
			else if ( strcmp(tok,"multiply")==0 && sscanf(eq+1,"%lg",&dval)==1 )
				lf->parameters.tune.nllm.lambda_mpy=dval;
			else if ( strcmp(tok,"iterations")==0 && sscanf(eq+1,"%d",&ival)==1 )
			 {	lf->parameters.tune.nllm.max_iter=ival;
				lf->parameters.tune.xmmc.niter=ival;
			 }
			else if ( strcmp(tok,"window")==0 && sscanf(eq+1,"%d",&ival)==1 )
				lf->parameters.tune.xmmc.window=ival;
		 }
	 }
	free(ptmp);
  }

 if ( fit_method==FIT_METHOD_LMND )
	lf->parameters.tune.nllm.numeric_derivs=1;

 if ( ! force_nonlinear )
  {	if ( fit_method==FIT_METHOD_CLLS )
		force_nonlinear=0;
	else if ( fit_method==FIT_METHOD_EMCE && lf->parameters.tune.emce.sub_method==FIT_METHOD_CLLS )
		force_nonlinear=0;
	else if ( fit_method==FIT_METHOD_NONE )
		force_nonlinear=0;
	else
		force_nonlinear=1;
  }

 if ( s_cfmt != NULL )
  {	char *cf=s_cfmt;
	if ( cf[0]=='%' ) cf++;
	lf->corrfm[0]='%';
	strncpy(lf->corrfm+1,cf,15);
  }
 else
	strcpy(lf->corrfm,LFIT_DEFAULT_CORR_FORMAT);

 if ( s_fmt != NULL && is_fit )
  {	i=extract_variable_formats(s_fmt,vars,nvar,dvars,ndvar);
	if ( i )	RETURN_ERROR("invalid variable format string");
  }

 if ( s_diff != NULL )
  {	i=extract_variable_differences(s_diff,vars,nvar);
	if ( i )	RETURN_ERROR("invalid variable differences");
	for ( i=0 ; i<nvar ; i++ )
	 {	if ( vars[i].diff <= 0.0 )
			RETURN_ERROR("negative or zero difference for variable");
	 }
  }

 if ( seed != 0 )
  {	if ( seed < 0 )
	 {	struct timeval tv;
		gettimeofday(&tv,NULL);
		seed=tv.tv_usec | ((tv.tv_sec & 0xFFF) << 20);
	 }
	random_seed(seed);
  }
 else
	random_seed(0);

 if ( s_pert != NULL )
  {	double xnoise=0.0;
	if ( sscanf(s_pert,"%lg",&xnoise)==1 && xnoise>=0.0 )
		db->xnoise=xnoise;
  }

 varsym=(psnsym *)malloc(sizeof(psnsym)*(nvar+1));
 for ( i=0 ; i<nvar ; i++ )
  {	varsym[i].type=T_VAR;
	varsym[i].major=i;
	varsym[i].name=vars[i].name;
  }
 varsym[nvar].type=0;
 varsym[nvar].major=0;
 varsym[nvar].name=NULL;
 lf->varsym=varsym;

 colsym=(psnsym *)malloc(sizeof(psnsym)*(ncol+1));
 for ( j=0 ; j<ncol ; j++ )
  {	colsym[j].type=T_VAR;
	colsym[j].major=j+nvar;
	colsym[j].name=cols[j].name;
  }
 colsym[ncol].type=0;
 colsym[ncol].major=0;
 colsym[ncol].name=NULL;
 db->colsym=colsym;

 mysyms[0]=varsym;
 mysyms[1]=db->colsym;
 mysyms[2]=lpg.pl_sym;
 mysyms[3]=NULL;

 for ( i=0 ; dvars != NULL && i<ndvar ; i++ )
  {	r=build_psn_derived_variable_expression(&dvars[i],mysyms,nvar,1);
	if ( r )	RETURN_ERROR("error in derived variable expression");
  }

 if ( is_fit && s_cnt != NULL )
  {	mysyms[1]=lpg.pl_sym;
	mysyms[2]=NULL;
	i=build_psn_constraint_sequences(lf,s_cnt,mysyms,vars,nvar);
	mysyms[1]=db->colsym;
	mysyms[2]=lpg.pl_sym;
	mysyms[3]=NULL;
	if ( i )	RETURN_ERROR("constraint error");
  }
 else
  {	lf->fconstraint.nc=0;
	lf->fconstraint.cmatrix=NULL;
	lf->fconstraint.cproject=NULL;
  }

 calc_derivatives=lfit_is_method_requries_derivatives(fit_method);
 if ( fit_method==FIT_METHOD_EMCE )
	calc_derivatives|=lfit_is_method_requries_derivatives(lf->parameters.tune.emce.sub_method);

 lf->is_linear=1;
 if ( is_fit )
  {	r=build_psn_fit_sequences(db,db->fncarg,db->deparg,
		NULL,0,
		db->errarg,
		mysyms,vars,nvar,force_nonlinear,calc_derivatives);
	if ( r )
	 {	free(varsym); free(colsym);
		RETURN_ERROR("PSN compile error %d",r);
	 }
	if ( ! db->is_linear )	lf->is_linear=0;
  }
 else
  {	j=build_psn_base_sequences(db,db->fncarg,NULL,0,mysyms,nvar);
	if ( j )
	 {	free(varsym); free(colsym);
		RETURN_ERROR("PSN base compile error %d",j);
	 }
  }

 if ( s_sep != NULL && is_fit )
	lfit_set_separated_linears(&lf->fconstraint,vars,nvar,s_sep);

 array_ptrs[0]=array_data;
 array_nrows[0]=nrow_in;
 array_ncols[0]=ncol_in;
 input.array_ptrs=array_ptrs;
 input.array_nrows=array_nrows;
 input.array_ncols=array_ncols;

 result->chain=chain;
 result->chain_max=chain_max;
 result->chain_count=0;
 result->nvar_chain=nvar+1;

 if ( fit_method==FIT_METHOD_FIMA || fit_method==FIT_METHOD_XMMC )
  {	int ci;
	result->cov_matrix=(double **)calloc(nvar>0?nvar:1,sizeof(double *));
	for ( ci=0 ; ci<nvar ; ci++ )
		result->cov_matrix[ci]=(double *)calloc(nvar>0?nvar:1,sizeof(double));
  }

 result->params=(double *)calloc(nvar>0?nvar:1,sizeof(double));
 result->errors=(double *)calloc(nvar>0?nvar:1,sizeof(double));
 result->used_mask=(int *)calloc(nrow_in>0?nrow_in:1,sizeof(int));

 if ( is_fit )
  {
	switch ( fit_method )
	 {
		case FIT_METHOD_CLLS:
		case FIT_METHOD_NLLM:
		case FIT_METHOD_LMND:
			r=fit_linear_or_nonlinear(lf,result,vars,nvar,
				errtype,
				errdump,resdump,is_dump_delta,weighted_sigma,
				(fit_method==FIT_METHOD_LMND?1:0),
				&input);
			break;
		case FIT_METHOD_MCMC:
			r=fit_markov_chain_monte_carlo(lf,result,vars,nvar,
				&input,callback,callback_data);
			break;
		case FIT_METHOD_MCHI:
			r=fit_map_chi2(lf,result,vars,nvar,
				&input,callback,callback_data);
			break;
		case FIT_METHOD_EMCE:
			r=fit_error_monte_carlo_estimation(lf,result,vars,nvar,
				&input,callback,callback_data);
			break;
		case FIT_METHOD_DHSX:
			r=fit_downhill_simplex(lf,result,vars,nvar,&input);
			break;
		case FIT_METHOD_XMMC:
			r=fit_extended_markov_chain_mc(lf,result,vars,nvar,
				&input,callback,callback_data);
			break;
		case FIT_METHOD_FIMA:
			r=fit_fisher_matrix_analysis(lf,result,vars,nvar,
				dvars,ndvar,&input,callback,callback_data);
			break;
		default:
			snprintf(result->error_msg,sizeof(result->error_msg),"unknown fit method %d",fit_method);
			result->error_code=-1;
			r=-1;
			break;
	 }
  }
 else
  {	double		*eval_wvars;
	fitinputrow	eval_fir;
	fitfunctdata	eval_ffd;

	eval_wvars=(double *)calloc(nvar+lf->maxncol,sizeof(double));

	eval_ffd.nvar=nvar;
	eval_ffd.wvars=eval_wvars;
	eval_ffd.functs=lpg.pl_funct;
	eval_ffd.lf=lf;

	result->eval_ncol=1;
	result->eval_nrow=nrow_in;
	result->eval_data=(double *)calloc(nrow_in>0?nrow_in:1,sizeof(double));
	result->nrow=nrow_in;

	for ( i=0 ; i<nrow_in ; i++ )
	 {	for ( j=0 ; j<ncol_in && j<lf->maxncol ; j++ )
			eval_wvars[nvar+j]=array_data[i*ncol_in+j];

		memset(&eval_fir,0,sizeof(fitinputrow));
		eval_fir.x=&eval_wvars[nvar];
		eval_fir.dbidx=0;

		fit_function((void *)&eval_fir,NULL,
			&result->eval_data[i],NULL,&eval_ffd);
	 }

	free(eval_wvars);
	r=0;
	result->error_code=0;
  }

 free(varsym);
 free(colsym);
 FREE_ALL_STRDUP;

 return(r);
}

/*****************************************************************************/

int lfit_python_ready(void)
{
    return 1;
}

/*****************************************************************************/
