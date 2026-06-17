# fitrans_internals.pxi — Fitrans internal operation classes
cdef class ImageTransform:
    """Internal: 2D polynomial image transformation engine."""
    cdef int _order, _method, _invert
    cdef double _ox, _oy, _scale
    cdef double *_dxfit
    cdef double *_dyfit
    cdef int _nvar
    cdef bint _trans_set

    _METHOD_MAP = {'bilinear':0,'integrate':1,'bicubic':2,
                   'spline_integrate':3,'lanczos3':4,'lanczos4':5}

    def __cinit__(self):
        self._dxfit = NULL; self._dyfit = NULL; self._trans_set = False

    def __init__(self, order=1, method='bilinear', ox=0.0, oy=0.0, scale=1.0,
                 inverse=False, T=None, trans=None, transformation=None):
        t = T or trans or transformation
        if t is not None:
            order = t.get('order', order)
            ox = t.get('offset', [ox, oy])[0]
            oy = t.get('offset', [ox, oy])[1]
            scale = t.get('scale', scale)
        if order < 1:
            raise ValueError(f"order must be >= 1, got {order}")
        self._order = order
        if isinstance(method, str):
            method = self._METHOD_MAP.get(method, 0)
        self._method = method
        self._ox = ox; self._oy = oy; self._scale = scale
        self._invert = 0 if inverse else 1
        self._nvar = (order + 1) * (order + 2) // 2
        self._dxfit = <double*>malloc(self._nvar * sizeof(double))
        self._dyfit = <double*>malloc(self._nvar * sizeof(double))
        if t is not None:
            self.set_trans(t['dxfit'], t['dyfit'])

    def set_trans(self, dxfit, dyfit):
        """Set transformation coefficients."""
        cdef int i
        for i in range(self._nvar):
            self._dxfit[i] = dxfit[i]; self._dyfit[i] = dyfit[i]
        self._trans_set = True

    def compose_shift(self, dx, dy):
        """Compose a pixel shift into the transformation (modifies coefficients in-place).
        Equivalent to --shift dx,dy in fitrans CLI."""
        if not self._trans_set:
            self._dxfit[0] = 0.0; self._dyfit[0] = 0.0
            self._dxfit[1] = 1.0; self._dyfit[2] = 1.0
            self._trans_set = True
        self._dxfit[0] += dx
        self._dyfit[0] += dy

    def compose_flip(self, x=False, y=False):
        """Compose X/Y flipping into the transformation (modifies coefficients in-place).
        Equivalent to --flip-x / --flip-y in fitrans CLI."""
        if not self._trans_set:
            self._dxfit[0] = 0.0; self._dyfit[0] = 0.0
            self._dxfit[1] = 1.0; self._dyfit[2] = 1.0
            self._trans_set = True
        if x:
            for i in range(self._nvar): self._dxfit[i] = -self._dxfit[i]
        if y:
            for i in range(self._nvar): self._dyfit[i] = -self._dyfit[i]

    def apply(self, img_data, img_mask=None, out_shape=None, ofx=0, ofy=0,
              weight=None, verbose=False):
        """Apply transformation to numpy 2D arrays.

        Args:
            img_data: numpy 2D float64 array, shape (sy, sx)
            img_mask: numpy 2D uint8 array or list of arrays (optional).
                      Multiple masks are OR'd together.
            out_shape: (nsy, nsx) output shape, defaults to input shape
            ofx, ofy: output offset in input coordinates
            weight: numpy 2D float64 array, same shape (optional).
                    Multiplied into image data before transformation.
            verbose: print progress info
        Returns:
            (out_data, out_mask) as numpy float64 / uint8 2D arrays
        """
        import numpy as np
        if not self._trans_set:
            raise RuntimeError("Call set_trans() first")
        sy, sx = img_data.shape
        if out_shape is None:
            nsy, nsx = sy, sx
        else:
            nsy, nsx = out_shape

        if weight is not None:
            img_data = np.ascontiguousarray(img_data, dtype=np.float64)
            img_data = img_data * weight

        if verbose:
            print(f'fitrans: {sx}x{sy} → {nsx}x{nsy}, order={self._order}, inverse={not bool(self._invert)}')
        # ensure contiguous C-order and correct dtypes
        in_arr = np.ascontiguousarray(img_data, dtype=np.float64)
        if isinstance(img_mask, list):
            mk_arr = np.ascontiguousarray(img_mask[0], dtype=np.uint8).copy()
            for m in img_mask[1:]:
                np.bitwise_or(mk_arr, np.ascontiguousarray(m, dtype=np.uint8), out=mk_arr)
        elif img_mask is not None:
            mk_arr = np.ascontiguousarray(img_mask, dtype=np.uint8)
        else:
            mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        # output buffers
        out_data = np.empty((nsy, nsx), dtype=np.float64)
        out_mask = np.empty((nsy, nsx), dtype=np.uint8)
        cdef double [:,::1] in_view = in_arr
        cdef unsigned char [:,::1] mk_view = mk_arr
        cdef double [:,::1] out_view = out_data
        cdef unsigned char [:,::1] om_view = out_mask
        cdef double *c_out = &out_view[0,0]
        cdef char   *c_mout = <char*>&om_view[0,0]
        fitrans_apply_cy(&in_view[0,0], <char*>&mk_view[0,0], sx, sy,
            c_out, c_mout, nsx, nsy, ofx, ofy,
            self._order, self._method, self._invert,
            self._ox, self._oy, self._scale, self._dxfit, self._dyfit)
        return out_data, out_mask

    def to_dict(self):
        """Export transformation params as dict."""
        method_str = {v:k for k,v in self._METHOD_MAP.items()}.get(self._method, 'bilinear')
        return {
            'order': self._order, 'method': method_str,
            'ox': self._ox, 'oy': self._oy, 'scale': self._scale,
            'inverse': not bool(self._invert),
            'dxfit': [self._dxfit[i] for i in range(self._nvar)],
            'dyfit': [self._dyfit[i] for i in range(self._nvar)],
        }

    def __dealloc__(self):
        if self._dxfit: free(self._dxfit)
        if self._dyfit: free(self._dyfit)

cdef class ShiftOp:
    """Simple pixel shift (matching --shift dx,dy in fitrans CLI)."""
    cdef double _dx, _dy
    def __init__(self, dx, dy):
        self._dx = dx; self._dy = dy
    def apply(self, img_data, img_mask=None, weight=None):
        import numpy as np
        result = img_data.copy()
        if weight is not None:
            np.multiply(result, weight, out=result)
        result = np.roll(np.roll(result, int(round(self._dx)), axis=1), int(round(self._dy)), axis=0)
        if isinstance(img_mask, list):
            mask_out = img_mask[0].copy()
            for m in img_mask[1:]: np.bitwise_or(mask_out, m, out=mask_out)
        elif img_mask is not None:
            mask_out = img_mask.copy()
        else:
            mask_out = np.zeros_like(img_data, dtype=np.uint8)
        mask_out = np.roll(np.roll(mask_out, int(round(self._dx)), axis=1), int(round(self._dy)), axis=0)
        return result, mask_out



cdef class FlipOp:
    """Flip image in X/Y (matching --flip-x, --flip-y in fitrans CLI)."""
    cdef bint _x, _y
    def __init__(self, x, y):
        self._x = x; self._y = y
    def apply(self, img_data, img_mask=None, weight=None):
        import numpy as np
        result = img_data.copy()
        if weight is not None:
            np.multiply(result, weight, out=result)
        if isinstance(img_mask, list):
            mask_out = img_mask[0].copy()
            for m in img_mask[1:]: np.bitwise_or(mask_out, m, out=mask_out)
        elif img_mask is not None:
            mask_out = img_mask.copy()
        else:
            mask_out = np.zeros_like(result, dtype=np.uint8)
        if self._x:
            result = np.fliplr(result)
            mask_out = np.fliplr(mask_out)
        if self._y:
            result = np.flipud(result)
            mask_out = np.flipud(mask_out)
        return result, mask_out


cdef class RepetitiveOp:
    """Repetitive image expansion (absolute output size, matching CLI --repetitive-xy)."""
    cdef int _nsx, _nsy, _ox, _oy
    def __init__(self, nsx, nsy, ox=0, oy=0):
        self._nsx = nsx; self._nsy = nsy; self._ox = ox; self._oy = oy
    def apply(self, img_data, img_mask=None, weight=None):
        import numpy as np
        sy, sx = img_data.shape
        nsx, nsy = self._nsx, self._nsy
        if nsx < sx: nsx = sx
        if nsy < sy: nsy = sy
        in_arr = np.ascontiguousarray(img_data, dtype=np.float64)
        if weight is not None:
            in_arr = in_arr * np.ascontiguousarray(weight, dtype=np.float64)
        if img_mask is not None: mk_arr = np.ascontiguousarray(img_mask, dtype=np.uint8)
        else: mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        out_data = np.empty((nsy, nsx), dtype=np.float64)
        out_mask = np.empty((nsy, nsx), dtype=np.uint8)
        cdef double [:,::1] in_v = in_arr
        cdef double [:,::1] out_v = out_data
        cdef unsigned char [:,::1] mk_v = mk_arr, om_v = out_mask
        fitrans_repetitive_cy(&in_v[0,0], <char*>&mk_v[0,0], sx, sy,
            &out_v[0,0], <char*>&om_v[0,0], nsx, nsy, self._ox, self._oy)
        return out_data, out_mask


cdef class InterleaveOp:
    """Interleave image expansion."""
    cdef int _fx, _fy, _ox, _oy, _median, _optimistic
    def __init__(self, fx, fy, ox=0, oy=0, median=False, optimistic_mask=False):
        self._fx = fx; self._fy = fy; self._ox = ox; self._oy = oy
        self._median = 1 if median else 0
        self._optimistic = 1 if optimistic_mask else 0
    def apply(self, img_data, img_mask=None):
        import numpy as np
        sy, sx = img_data.shape
        nsx, nsy = self._fx, self._fy
        in_arr = np.ascontiguousarray(img_data, dtype=np.float64)
        if img_mask is not None: mk_arr = np.ascontiguousarray(img_mask, dtype=np.uint8)
        else: mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        out_data = np.empty((nsy, nsx), dtype=np.float64)
        out_mask = np.empty((nsy, nsx), dtype=np.uint8)
        cdef double [:,::1] in_v = in_arr
        cdef double [:,::1] out_v = out_data
        cdef unsigned char [:,::1] mk_v = mk_arr, om_v = out_mask
        fitrans_interleave_cy(&in_v[0,0], <char*>&mk_v[0,0], sx, sy,
            &out_v[0,0], <char*>&om_v[0,0], nsx, nsy,
            self._ox, self._oy, self._median, 0, self._optimistic)
        return out_data, out_mask


cdef class NoiseOp:
    """Noise estimation operator."""
    def apply(self, img_data, img_mask=None):
        import numpy as np
        sy, sx = img_data.shape
        in_arr = np.ascontiguousarray(img_data, dtype=np.float64)
        if img_mask is not None: mk_arr = np.ascontiguousarray(img_mask, dtype=np.uint8)
        else: mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        out_data = np.empty((sy, sx), dtype=np.float64)
        out_mask = np.empty((sy, sx), dtype=np.uint8)
        cdef double [:,::1] in_view = in_arr, out_view = out_data
        cdef unsigned char [:,::1] mk_view = mk_arr, om_view = out_mask
        fitrans_noise_cy(&in_view[0,0], <char*>&mk_view[0,0], sx, sy,
            &out_view[0,0], <char*>&om_view[0,0])
        return out_data, out_mask


cdef class ZoomOp:
    """Zoom/magnify operator."""
    cdef int _sx, _sy, _ox, _oy, _raw
    def __init__(self, sx, sy, ox=0, oy=0, raw=False):
        self._sx = sx; self._sy = sy; self._ox = ox; self._oy = oy; self._raw = 1 if raw else 0
    def apply(self, img_data, img_mask=None):
        import numpy as np
        sy, sx = img_data.shape
        nsx, nsy = sx * self._sx, sy * self._sy
        in_arr = np.ascontiguousarray(img_data, dtype=np.float64)
        if img_mask is not None: mk_arr = np.ascontiguousarray(img_mask, dtype=np.uint8)
        else: mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        out_data = np.empty((nsy, nsx), dtype=np.float64)
        out_mask = np.empty((nsy, nsx), dtype=np.uint8)
        cdef double [:,::1] in_v = in_arr, out_v = out_data
        cdef unsigned char [:,::1] mk_v = mk_arr, om_v = out_mask
        fitrans_zoom_cy(&in_v[0,0], <char*>&mk_v[0,0], sx, sy,
            &out_v[0,0], <char*>&om_v[0,0], nsx, nsy,
            self._ox, self._oy, self._sx, self._sy, self._raw)
        return out_data, out_mask


cdef class ShrinkOp:
    """Shrink operator."""
    cdef int _sx, _sy, _ox, _oy, _med, _tm, _am
    def __init__(self, sx, sy, ox=0, oy=0, median=False, truncated_mean=0, optimistic_mask=False):
        self._sx = sx; self._sy = sy; self._ox = ox; self._oy = oy
        self._med = 1 if median else 0; self._tm = truncated_mean; self._am = 1 if optimistic_mask else 0
    def apply(self, img_data, img_mask=None):
        import numpy as np
        sy, sx = img_data.shape
        nsx = (sx + self._sx - 1) // self._sx
        nsy = (sy + self._sy - 1) // self._sy
        in_arr = np.ascontiguousarray(img_data, dtype=np.float64)
        if img_mask is not None: mk_arr = np.ascontiguousarray(img_mask, dtype=np.uint8)
        else: mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        out_data = np.empty((nsy, nsx), dtype=np.float64)
        out_mask = np.empty((nsy, nsx), dtype=np.uint8)
        cdef double [:,::1] in_v = in_arr, out_v = out_data
        cdef unsigned char [:,::1] mk_v = mk_arr, om_v = out_mask
        fitrans_shrink_cy(&in_v[0,0], <char*>&mk_v[0,0], sx, sy,
            &out_v[0,0], <char*>&om_v[0,0], nsx, nsy,
            self._ox, self._oy, self._sx, self._sy,
            self._med, self._tm, self._am)
        return out_data, out_mask


cdef class SmoothOp:
    """Smoothing operator."""
    cdef int _st, _xo, _yo, _pf, _fhx, _fhy, _nit, _mu, _de
    cdef double _fr, _lo, _up
    def __init__(self, st, xo, yo, pf, fhx, fhy, fr, nit, lo, up, mu, de):
        self._st = st; self._xo = xo; self._yo = yo
        self._pf = pf; self._fhx = fhx; self._fhy = fhy
        self._fr = fr; self._nit = nit; self._lo = lo; self._up = up
        self._mu = 1 if mu else 0; self._de = 1 if de else 0
    def apply(self, img_data, img_mask=None):
        import numpy as np
        sy, sx = img_data.shape
        in_arr = np.ascontiguousarray(img_data, dtype=np.float64)
        if img_mask is not None: mk_arr = np.ascontiguousarray(img_mask, dtype=np.uint8)
        else: mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        out_data = np.empty((sy, sx), dtype=np.float64)
        out_mask = np.empty((sy, sx), dtype=np.uint8)
        cdef double [:,::1] in_v = in_arr, out_v = out_data
        cdef unsigned char [:,::1] mk_v = mk_arr, om_v = out_mask
        fitrans_smooth_cy(&in_v[0,0], <char*>&mk_v[0,0], sx, sy,
            &out_v[0,0], <char*>&om_v[0,0],
            self._st, self._xo, self._yo,
            self._pf, self._fhx, self._fhy,
            self._fr, self._nit, self._lo, self._up,
            self._mu, self._de)
        return out_data, out_mask

