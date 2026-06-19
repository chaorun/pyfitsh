"""Lfit — general-purpose curve fitting (Python wrapper around C lfit engine)."""

import numpy as np


METHODS = {
    'clls': 'Classical Linear Least Squares',
    'nllm': 'Nonlinear Levenberg-Marquardt',
    'mcmc': 'Markov Chain Monte-Carlo',
    'mchi': 'Chi-squared grid mapping',
    'emce': 'Monte-Carlo error estimation',
    'dhsx': 'Downhill simplex',
    'xmmc': 'Extended Markov Chain Monte-Carlo',
    'lmnd': 'Levenberg-Marquardt (numerical derivatives)',
    'fima': 'Fisher information matrix analysis',
}


class LfitResult:
    """Container for lfit fitting results."""

    def __init__(self, raw, meta=None):
        self.params = np.array(raw.get('params', []))
        self.errors = np.array(raw.get('errors', []))
        self.chi2 = raw.get('chi2', 0.0)
        self.nrow = raw.get('nrow', 0)
        self.nused = raw.get('nused', 0)
        self.residual_sigma = raw.get('residual_sigma', 0.0)
        self.acceptance = raw.get('acceptance', 0.0)
        self.used_mask = raw.get('used_mask', [])
        self.chain = raw.get('chain', None)
        self.cov_matrix = raw.get('cov_matrix', None)
        self.eval_data = raw.get('eval_data', None)
        self.return_code = raw.get('return_code', -1)
        self.error_code = raw.get('error_code', -1)
        self.error_msg = raw.get('error_msg', '')
        self.raw = raw
        self.meta = meta or {}

    @property
    def ok(self):
        return self.return_code == 0 and self.error_code == 0

    @property
    def nrejected(self):
        return self.nrow - self.nused

    @property
    def variables(self):
        return self.meta.get('variables', [])

    @property
    def param_dict(self):
        vnames = self.variables
        if len(vnames) == len(self.params):
            return dict(zip(vnames, self.params))
        return {}

    @property
    def corr_matrix(self):
        if self.cov_matrix is None:
            return None
        cov = self.cov_matrix
        n = cov.shape[0]
        diag = np.sqrt(np.diag(cov))
        diag[diag == 0] = 1.0
        return cov / np.outer(diag, diag)

    def to_dict(self):
        out = dict(self.raw)
        out['meta'] = self.meta
        return out

    KEYS = (
        'params', 'errors', 'chi2', 'nrow', 'nused',
        'residual_sigma', 'acceptance', 'used_mask', 'chain',
        'cov_matrix', 'corr_matrix', 'eval_data',
        'return_code', 'error_code', 'error_msg',
        'ok', 'nrejected', 'variables', 'param_dict', 'meta',
    )

    def __getitem__(self, key):
        if key in self.KEYS:
            return getattr(self, key)
        raise KeyError(key)

    def __contains__(self, key):
        return key in self.KEYS

    def keys(self):
        return list(self.KEYS)

    def values(self):
        return [getattr(self, k) for k in self.KEYS]

    def items(self):
        return [(k, getattr(self, k)) for k in self.KEYS]

    def __repr__(self):
        if not self.ok:
            return f"LfitResult(error={self.error_msg!r})"
        parts = []
        vnames = self.variables
        for i, p in enumerate(self.params):
            name = vnames[i] if i < len(vnames) else f"p{i}"
            if i < len(self.errors) and self.errors[i] != 0:
                parts.append(f"{name}={p:.6g} +/- {self.errors[i]:.6g}")
            else:
                parts.append(f"{name}={p:.6g}")
        pstr = ', '.join(parts)
        return (f"LfitResult({pstr}, chi2={self.chi2:.6g}, "
                f"nrow={self.nrow}, nused={self.nused})")

    def __str__(self):
        return self.__repr__()


class Lfit:
    """General-purpose curve fitting engine.

    Usage:
        lf = Lfit(method='clls', rejection_niter=3, rejection_level=3.0)
        result = lf.fit(data, variables='a,b', columns='x:1,y:2',
                        function='a+b*x', dependent='y')
        print(result.params)
    """

    def __init__(self, method='clls',
                 rejection_level=3.0, rejection_niter=0,
                 weighted_sigma=0,
                 parameters=None, differences=None, separate=None,
                 perturbations=None, seed=0, mc_iterations=1000,
                 macros=None, constraints=None,
                 errdump=0, format=None, correlation_format=None,
                 derived_variables=None,
                 is_dump_delta=0, resdump=0, force_nonlinear=0,
                 columns_output=None):
        if method not in METHODS:
            raise ValueError(
                f"unknown method '{method}', use one of: {list(METHODS.keys())}")
        self.method = method
        self.rejection_level = rejection_level
        self.rejection_niter = rejection_niter
        self.weighted_sigma = weighted_sigma
        self.parameters = parameters
        self.differences = differences
        self.separate = separate
        self.perturbations = perturbations
        self.seed = seed
        self.mc_iterations = mc_iterations
        self.macros = macros
        self.constraints = constraints
        self.errdump = errdump
        self.format = format
        self.correlation_format = correlation_format
        self.derived_variables = derived_variables
        self.is_dump_delta = is_dump_delta
        self.resdump = resdump
        self.force_nonlinear = force_nonlinear
        self.columns_output = columns_output

    def fit(self, data, variables, columns, function, dependent=None,
            error=None, weight=None):
        from pyfitsh.core import lfit_fit as c_lfit_fit

        raw = c_lfit_fit(
            data,
            variables=variables,
            columns=columns,
            function=function,
            dependent=dependent,
            error=error,
            weight=weight,
            method=self.method,
            parameters=self.parameters,
            differences=self.differences,
            separate=self.separate,
            perturbations=self.perturbations,
            seed=self.seed,
            mc_iterations=self.mc_iterations,
            rejection_level=self.rejection_level,
            rejection_niter=self.rejection_niter,
            weighted_sigma=self.weighted_sigma,
            macros=self.macros,
            constraints=self.constraints,
            errdump=self.errdump,
            format=self.format,
            correlation_format=self.correlation_format,
            derived_variables=self.derived_variables,
            is_dump_delta=self.is_dump_delta,
            resdump=self.resdump,
            force_nonlinear=self.force_nonlinear,
            columns_output=self.columns_output,
        )

        varlist = [v.split('=')[0].split(':')[0].strip()
                   for v in variables.split(',')]

        meta = {
            'method': self.method,
            'method_name': METHODS.get(self.method, ''),
            'function': function,
            'dependent': dependent,
            'columns': columns,
            'variables': varlist,
            'error': error,
            'weight': weight,
        }

        return LfitResult(raw, meta)


def lfit_fit(data, variables, columns, function, dependent=None,
             error=None, weight=None, method='clls',
             parameters=None, differences=None, separate=None,
             perturbations=None, seed=0, mc_iterations=1000,
             rejection_level=3.0, rejection_niter=0, weighted_sigma=0,
             macros=None, constraints=None,
             errdump=0, format=None, correlation_format=None,
             derived_variables=None,
             is_dump_delta=0, resdump=0, force_nonlinear=0,
             columns_output=None):
    """Convenience function for one-shot curve fitting.

    Returns LfitResult.
    """
    lf = Lfit(method=method,
              rejection_level=rejection_level, rejection_niter=rejection_niter,
              weighted_sigma=weighted_sigma,
              parameters=parameters, differences=differences,
              separate=separate, perturbations=perturbations,
              seed=seed, mc_iterations=mc_iterations,
              macros=macros, constraints=constraints,
              errdump=errdump, format=format,
              correlation_format=correlation_format,
              derived_variables=derived_variables,
              is_dump_delta=is_dump_delta, resdump=resdump,
              force_nonlinear=force_nonlinear,
              columns_output=columns_output)
    return lf.fit(data, variables, columns, function, dependent,
                  error=error, weight=weight)
