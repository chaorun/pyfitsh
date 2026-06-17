# poly_helpers.pxi — polynomial evaluation helper
def eval_2d_poly_py(x, y, coeff, order):
    """Equivalent to C eval_2d_poly(x, y, order, coeff, 0, 0, 1)."""
    idx = 0; result = 0.0
    for k in range(order + 1):
        for j in range(k + 1):
            term = 1.0
            for p in range(1, k - j + 1): term *= x / p
            for p in range(1, j + 1): term *= y / p
            result += coeff[idx] * term
            idx += 1
    return result

# polyfit.pxi — polynomial fitting
cdef class PolyFitter:
    """
    Fit 2D polynomial coefficients from (x, y) → target value data,
    with optional iterative sigma rejection.

    Usage:
        fitter = PolyFitter(order=2, niter=3, rejlevel=3.0)
        fitter.fit(ref_x, ref_y, target_x, weights=err_list)
        print(fitter.coeff, fitter.residual)

        # weights from stellar magnitudes (mag → flux^power)
        fitter = PolyFitter(order=2, magnitude=True, power=2.0)
        fitter.fit(x, y, vals, weights=mag_list)
    """

    cdef int _order, _nvar, _niter
    cdef double _rejlevel, _ox, _oy, _scale
    cdef double *_coeff
    cdef double _residual
    cdef int _npts, _nrej
    cdef bint _w_magnitude
    cdef double _w_power

    def __cinit__(self):
        self._coeff = NULL

    def __init__(self, order=2, niter=0, rejlevel=3.0,
                 ox=0.0, oy=0.0, scale=1.0,
                 magnitude=False, power=1.0):
        self._order = order
        self._nvar = (order + 1) * (order + 2) // 2
        self._niter = niter
        self._rejlevel = rejlevel
        self._ox = ox; self._oy = oy; self._scale = scale
        self._w_magnitude = magnitude
        self._w_power = power
        self._coeff = <double*>malloc(self._nvar * sizeof(double))
        self._residual = 0.0
        self._nrej = 0

    def fit(self, x, y, vals, weights=None):
        """Fit polynomial. Returns self for chaining."""
        cdef int n = len(x)
        assert len(y) == n == len(vals)
        self._npts = n

        cdef double *cx = <double*>malloc(n * sizeof(double))
        cdef double *cy = <double*>malloc(n * sizeof(double))
        cdef double *cv = <double*>malloc(n * sizeof(double))
        cdef double *cw = <double*>malloc(n * sizeof(double))
        cdef int i, irej
        cdef double s, sdd, w, fd, sig, ss, ssdd
        for i in range(n):
            cx[i] = x[i]; cy[i] = y[i]; cv[i] = vals[i]
            if self._w_magnitude and weights is not None:
                flux = 10.0 ** (-0.4 * weights[i])
                cw[i] = flux ** self._w_power
            elif weights is not None:
                cw[i] = weights[i]
            else:
                cw[i] = 1.0

        cdef int niter = self._niter
        if niter < 0: niter = 0

        for irej in range(niter + 1):
            poly_fit_cy(cx, cy, cv, cw, n, self._order,
                self._ox, self._oy, self._scale, self._coeff)

            if irej >= niter:
                break

            # Compute residual
            s = 0.0; sdd = 0.0
            for i in range(n):
                if cw[i] <= 0: continue
                fd = eval_2d_poly_py(cx[i], cy[i], self.coeff_list, self._order) - cv[i]
                w = cw[i]
                s += w; sdd += fd * fd * w
            if s <= 0: break
            sig = (sdd / s) ** 0.5

            # Reject outliers
            self._nrej = 0
            for i in range(n):
                if cw[i] <= 0: continue
                fd = abs(eval_2d_poly_py(cx[i], cy[i], self.coeff_list, self._order) - cv[i])
                if fd > self._rejlevel * sig:
                    cw[i] = 0.0
                    self._nrej += 1

        # Final residual
        ss = 0.0; ssdd = 0.0
        for i in range(n):
            w = cw[i]
            if w <= 0: continue
            fd = eval_2d_poly_py(cx[i], cy[i], self.coeff_list, self._order) - cv[i]
            ss += w; ssdd += fd * fd * w
        self._residual = (ssdd / ss) ** 0.5 if ss > 0 else 0.0

        free(cx); free(cy); free(cv); free(cw)
        return self

    @property
    def coeff(self):
        if self._coeff == NULL: return []
        return [self._coeff[i] for i in range(self._nvar)]

    @property
    def coeff_list(self):
        return self.coeff

    @property
    def residual(self): return self._residual

    @property
    def nrej(self): return self._nrej

    def __dealloc__(self):
        if self._coeff: free(self._coeff)


