LFIT_METHODS = {
    'clls': FIT_METHOD_CLLS,
    'nllm': FIT_METHOD_NLLM,
    'mcmc': FIT_METHOD_MCMC,
    'mchi': FIT_METHOD_MCHI,
    'emce': FIT_METHOD_EMCE,
    'dhsx': FIT_METHOD_DHSX,
    'xmmc': FIT_METHOD_XMMC,
    'lmnd': FIT_METHOD_LMND,
    'fima': FIT_METHOD_FIMA,
}

def lfit_ready():
    return lfit_python_ready()

def lfit_fit(data,
             str variables not None,
             str columns not None,
             str function not None,
             str dependent=None,
             str error=None,
             str weight=None,
             str method='clls',
             str parameters=None,
             str differences=None,
             str separate=None,
             str perturbations=None,
             int seed=0,
             int mc_iterations=1000,
             double rejection_level=3.0,
             int rejection_niter=0,
             int weighted_sigma=0,
             list macros=None,
             str constraints=None,
             int errdump=0,
             str format=None,
             str correlation_format=None,
             str derived_variables=None,
             int is_dump_delta=0,
             int resdump=0,
             int force_nonlinear=0,
             str columns_output=None):
    import numpy as np

    cdef int nrow, ncol, r, i, nm
    cdef int fit_method
    cdef lfit_result result
    cdef double *data_ptr
    cdef int nvar_fit, stride, chain_capacity
    cdef double *chain_buf

    if method in LFIT_METHODS:
        fit_method = LFIT_METHODS[method]
    else:
        raise ValueError(f"unknown fit method '{method}', use one of: {list(LFIT_METHODS.keys())}")

    arr = np.ascontiguousarray(data, dtype=np.float64)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    elif arr.ndim != 2:
        raise ValueError("data must be 1D or 2D array")

    cdef double[:, ::1] mv = arr
    nrow = mv.shape[0]
    ncol = mv.shape[1]
    data_ptr = &mv[0, 0]

    cdef bytes b_var = variables.encode('utf-8')
    cdef bytes b_col = columns.encode('utf-8')
    cdef bytes b_func = function.encode('utf-8')
    cdef bytes b_dep
    cdef bytes b_err, b_wgt, b_par, b_diff, b_sep, b_pert
    cdef bytes b_cnt, b_fmt, b_cfmt, b_dvr, b_co

    cdef char *c_var  = b_var
    cdef char *c_col  = b_col
    cdef char *c_func = b_func
    cdef char *c_dep  = NULL
    cdef char *c_err  = NULL
    cdef char *c_wgt  = NULL
    cdef char *c_par  = NULL
    cdef char *c_diff = NULL
    cdef char *c_sep  = NULL
    cdef char *c_pert = NULL
    cdef char *c_cnt  = NULL
    cdef char *c_fmt  = NULL
    cdef char *c_cfmt = NULL
    cdef char *c_dvr  = NULL
    cdef char *c_co   = NULL

    if dependent is not None:
        b_dep = dependent.encode('utf-8'); c_dep = b_dep
    if error is not None:
        b_err = error.encode('utf-8'); c_err = b_err
    if weight is not None:
        b_wgt = weight.encode('utf-8'); c_wgt = b_wgt
    if parameters is not None:
        b_par = parameters.encode('utf-8'); c_par = b_par
    if differences is not None:
        b_diff = differences.encode('utf-8'); c_diff = b_diff
    if separate is not None:
        b_sep = separate.encode('utf-8'); c_sep = b_sep
    if perturbations is not None:
        b_pert = perturbations.encode('utf-8'); c_pert = b_pert
    if constraints is not None:
        b_cnt = constraints.encode('utf-8'); c_cnt = b_cnt
    if format is not None:
        b_fmt = format.encode('utf-8'); c_fmt = b_fmt
    if correlation_format is not None:
        b_cfmt = correlation_format.encode('utf-8'); c_cfmt = b_cfmt
    if derived_variables is not None:
        b_dvr = derived_variables.encode('utf-8'); c_dvr = b_dvr
    if columns_output is not None:
        b_co = columns_output.encode('utf-8'); c_co = b_co

    cdef char **c_macros = NULL
    cdef list macro_bytes = []
    if macros is not None and len(macros) > 0:
        nm = len(macros)
        c_macros = <char **>malloc((nm + 1) * sizeof(char *))
        for i in range(nm):
            mb = macros[i].encode('utf-8') if isinstance(macros[i], str) else macros[i]
            macro_bytes.append(mb)
            c_macros[i] = mb
        c_macros[nm] = NULL

    nvar_fit = len(variables.split(','))
    stride = nvar_fit + 1
    chain_capacity = 0
    chain_buf = NULL
    if fit_method in (FIT_METHOD_MCMC, FIT_METHOD_XMMC, FIT_METHOD_EMCE,
                      FIT_METHOD_FIMA, FIT_METHOD_MCHI):
        chain_capacity = mc_iterations + 2
        chain_buf = <double *>malloc(chain_capacity * stride * sizeof(double))

    memset(&result, 0, sizeof(lfit_result))

    with nogil:
        r = lfit_python_apply(
            data_ptr, nrow, ncol,
            c_var, c_col, c_func, c_dep,
            c_err, c_wgt,
            fit_method,
            c_par, c_diff, c_sep, c_pert,
            seed, mc_iterations,
            rejection_level, rejection_niter, weighted_sigma,
            c_macros,
            c_cnt,
            errdump, c_fmt, c_cfmt, c_dvr,
            is_dump_delta, resdump, force_nonlinear, c_co,
            &result,
            chain_buf, chain_capacity,
            NULL, NULL)

    if c_macros != NULL:
        free(c_macros)

    out = {}
    out['return_code'] = r
    out['error_code'] = result.error_code
    out['error_msg'] = result.error_msg.decode('utf-8', errors='replace').rstrip('\x00')
    out['chi2'] = result.chi2
    out['nrow'] = result.nrow
    out['nused'] = result.nused
    out['residual_sigma'] = result.residual_sigma
    out['acceptance'] = result.acceptance

    cdef int nvar_out = 0
    if b_var:
        nvar_out = len(variables.split(','))

    if result.params != NULL and nvar_out > 0:
        out['params'] = [result.params[i] for i in range(nvar_out)]
    else:
        out['params'] = []

    if result.errors != NULL and nvar_out > 0:
        out['errors'] = [result.errors[i] for i in range(nvar_out)]
    else:
        out['errors'] = []

    if result.used_mask != NULL and result.nrow > 0:
        out['used_mask'] = [bool(result.used_mask[i]) for i in range(result.nrow)]
    else:
        out['used_mask'] = []

    if chain_buf != NULL and result.chain_count > 0:
        chain_list = [chain_buf[i] for i in range(result.chain_count * stride)]
        out['chain'] = np.array(chain_list).reshape(result.chain_count, stride)
    else:
        out['chain'] = None

    if chain_buf != NULL:
        free(chain_buf)

    if result.cov_matrix != NULL and nvar_fit > 0:
        cov_rows = []
        for i in range(nvar_fit):
            if result.cov_matrix[i] != NULL:
                cov_rows.append([result.cov_matrix[i][j] for j in range(nvar_fit)])
                free(result.cov_matrix[i])
        free(result.cov_matrix)
        if len(cov_rows) == nvar_fit:
            out['cov_matrix'] = np.array(cov_rows)
        else:
            out['cov_matrix'] = None
    else:
        out['cov_matrix'] = None

    if result.eval_data != NULL and result.eval_nrow > 0 and result.eval_ncol > 0:
        eval_list = [result.eval_data[i] for i in range(result.eval_nrow * result.eval_ncol)]
        out['eval_data'] = np.array(eval_list).reshape(result.eval_nrow, result.eval_ncol)
        free(result.eval_data)
    else:
        out['eval_data'] = None

    return out
