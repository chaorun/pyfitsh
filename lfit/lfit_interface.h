#ifndef LFIT_INTERFACE_H
#define LFIT_INTERFACE_H

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
	lfit_step_callback callback,void *callback_data);

int lfit_python_ready(void);

#ifndef FIT_METHOD_NONE
#define FIT_METHOD_NONE 0
#define FIT_METHOD_CLLS 1
#define FIT_METHOD_NLLM 2
#define FIT_METHOD_MCMC 3
#define FIT_METHOD_MCHI 4
#define FIT_METHOD_EMCE 5
#define FIT_METHOD_DHSX 6
#define FIT_METHOD_XMMC 7
#define FIT_METHOD_LMND 8
#define FIT_METHOD_FIMA 9
#endif

#endif
