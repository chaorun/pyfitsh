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


def _parse_kernel_list(raw):
    """Parse kernel list text into structured dict matching CLI kernel.txt values."""
    result = {'global': {}, 'kernels': []}
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('type='):
            result['global']['type'] = int(line.split('=')[1])
        elif line.startswith('offset='):
            parts = line.split('=')[1].split(',')
            result['global']['offset'] = [float(parts[0]), float(parts[1])]
        elif line.startswith('scale='):
            result['global']['scale'] = float(line.split('=')[1])
        elif line.startswith('kernel='):
            kd = {}
            parts = line.split()
            in_coeff = False
            for p in parts:
                if '=' in p:
                    k, _, v = p.partition('=')
                    if k == 'kernel':
                        kd['index'] = int(v)
                    elif k == 'type':
                        kd['type'] = int(v)
                    elif k == 'order':
                        kd['order'] = int(v)
                    elif k == 'hsize':
                        kd['hsize'] = int(v)
                    elif k == 'sigma':
                        kd['sigma'] = float(v)
                    elif k == 'coeff':
                        kd['coeff'] = [float(v)]
                        in_coeff = True
                elif in_coeff:
                    # bare number tokens after coeff= are additional coefficients
                    try:
                        kd['coeff'].append(float(p))
                    except ValueError:
                        in_coeff = False
            if 'coeff' not in kd:
                kd['coeff'] = []
            result['kernels'].append(kd)
    return result

cdef class ConvolutionOp:
    """Kernel fitting + convolution (ficonv)."""
    cdef bytes _kernel_spec
    cdef int _iterations, _method, _divide
    cdef double _rejection_level, _gain
    cdef bint _verbose, _unity_kernels

    # ficonv fitting modes (bit flags from ficonv_core.c)
    FM_NORMAL       = 0
    FM_MASKED       = 1
    FM_WEIGHTED     = 2
    FM_BGITERATIVE  = 4

    def __init__(self, kernel, iterations=0, rejection_level=3.0,
                 masked=False, weighted=False, background_iterative=False,
                 divide=32, gain=1.0, verbose=False,
                 unity_kernels=False):
        self._kernel_spec = (kernel or "").encode('utf-8')
        self._iterations = iterations
        self._rejection_level = rejection_level
        self._divide = divide; self._gain = gain
        self._verbose = verbose
        self._unity_kernels = unity_kernels
        # method: combine flags, WEIGHTED implies MASKED
        self._method = 0
        if masked:                self._method |= ConvolutionOp.FM_MASKED
        if weighted:              self._method |= ConvolutionOp.FM_WEIGHTED
        if background_iterative:  self._method |= ConvolutionOp.FM_BGITERATIVE

    def fit(self, ref_data, img_data, ref_mask=None, img_mask=None,
            stamps=None, add_to=None, maskinfo=None,
            output_subtracted=False, output_kernel_list=False,
            kernel_dict=None):
        """Fit kernel and convolve.

        Parameters
        ----------
        ref_data, img_data : float64 ndarray (2D)   — -r / -i
        ref_mask, img_mask : uint8 ndarray (optional) — -M
        stamps : [(x,y,sx,sy),...]                    — -t
        add_to : float64 ndarray                      — -a
        output_subtracted : bool                      — --output-subtracted
        output_kernel_list : bool                     — --output-kernel-list
        kernel_dict : dict                            — pre-fitted kernel list, skip fitting

        Returns
        -------
        (convolved, mask[, subtracted][, add_to_result][, kernel_list])
        """
        import numpy as np
        sy, sx = ref_data.shape
        in_ref = np.ascontiguousarray(ref_data, dtype=np.float64)
        in_img = np.ascontiguousarray(img_data, dtype=np.float64)
        if ref_mask is not None:
            mr = np.ascontiguousarray(ref_mask, dtype=np.uint8)
        else:
            mr = np.zeros((sy,sx), dtype=np.uint8)
        if maskinfo is not None:
            mi_mask = parse_maskinfo(maskinfo, sx, sy)
            mr = np.bitwise_or(mr, mi_mask)

        if stamps is not None:
            # Original CLI: stamps define INCLUSION zones for fitting.
            # fitmask: stamp pixels = 0 (allowed), non-stamp pixels = MASK_OUTER (excluded).
            inmask = np.full((sy, sx), 0x08, dtype=np.uint8)  # 0x08 = MASK_OUTER
            for sx0, sy0, ssx, ssy in stamps:
                x0, y0 = int(sx0), int(sy0)
                inmask[y0:y0+ssy, x0:x0+ssx] = 0
            if img_mask is not None:
                mi = np.ascontiguousarray(img_mask, dtype=np.uint8)
            else:
                mi = np.zeros((sy, sx), dtype=np.uint8)
        elif img_mask is not None:
            mi = np.ascontiguousarray(img_mask, dtype=np.uint8)
            inmask = None
        else:
            mi = np.zeros((sy, sx), dtype=np.uint8)
            inmask = None

        if self._verbose:
            print(f'ficonv: fitting {sx}x{sy}, method={self._method}, divide={self._divide}, gain={self._gain}')

        out_data = np.empty((sy,sx), dtype=np.float64)
        out_mask = np.zeros((sy,sx), dtype=np.uint8)
        sub_data = np.empty((sy,sx), dtype=np.float64) if output_subtracted else None
        cdef double [:,::1] ir = in_ref, ii = in_img, od = out_data
        cdef unsigned char [:,::1] omr = mr, omi = mi, omsk = out_mask
        cdef unsigned char [:,::1] oinmask_view
        cdef char *inmask_ptr = NULL
        if inmask is not None:
            oinmask_view = np.ascontiguousarray(inmask, dtype=np.uint8)
            inmask_ptr = <char*>&oinmask_view[0,0]
        cdef char *klist_out = NULL
        cdef double *sub_ptr = NULL
        cdef double [:,::1] sd
        if sub_data is not None:
            sd = sub_data
            sub_ptr = &sd[0,0]
        cdef double *add_ptr = NULL
        cdef double [:,::1] ad
        if add_to is not None:
            ad = np.ascontiguousarray(add_to, dtype=np.float64)
            add_ptr = &ad[0,0]

        # prefit kernel dict → flat arrays
        cdef int prefit_nk = 0
        cdef int [:] prefit_types_v
        cdef int [:] prefit_orders_v
        cdef int [:] prefit_ncoeffs_v
        cdef double [:] prefit_coeffs_v
        cdef int *pf_types = NULL, *pf_orders = NULL, *pf_ncoeffs = NULL
        cdef double *pf_coeffs = NULL
        cdef double pf_ox = 0.0, pf_oy = 0.0, pf_scale = 1.0
        if kernel_dict is not None and 'kernels' in kernel_dict:
            kds = kernel_dict['kernels']
            prefit_nk = len(kds)
            prefit_types_v = np.zeros(prefit_nk, dtype=np.int32)
            prefit_orders_v = np.zeros(prefit_nk, dtype=np.int32)
            prefit_ncoeffs_v = np.zeros(prefit_nk, dtype=np.int32)
            total_coeffs = sum(len(kd['coeff']) for kd in kds)
            prefit_coeffs_v = np.zeros(total_coeffs, dtype=np.float64)
            pc_idx = 0
            for ki, kd in enumerate(kds):
                prefit_types_v[ki] = kd['type']
                prefit_orders_v[ki] = kd['order']
                prefit_ncoeffs_v[ki] = len(kd['coeff'])
                for cv in kd['coeff']:
                    prefit_coeffs_v[pc_idx] = cv
                    pc_idx += 1
            pf_types = &prefit_types_v[0]
            pf_orders = &prefit_orders_v[0]
            pf_ncoeffs = &prefit_ncoeffs_v[0]
            pf_coeffs = &prefit_coeffs_v[0]
            g = kernel_dict.get('global', {})
            pf_ox = g.get('offset', [0.0, 0.0])[0]
            pf_oy = g.get('offset', [0.0, 0.0])[1]
            pf_scale = g.get('scale', 1.0)

        ficonv_fit_cy(&ir[0,0], <char*>&omr[0,0], &ii[0,0], <char*>&omi[0,0],
            inmask_ptr,
            sx, sy, self._kernel_spec,
            self._method, self._divide,
            self._iterations, self._rejection_level, self._gain,
            1 if self._verbose else 0,
            &od[0,0], <char*>&omsk[0,0],
            sub_ptr,
            add_ptr,
            self._unity_kernels, &klist_out,
            <char*>NULL, 0, 0, <char*>NULL,
            prefit_nk, pf_types, pf_orders, pf_ncoeffs, pf_coeffs,
            pf_ox, pf_oy, pf_scale)

        results = [out_data, out_mask]

        if output_subtracted:
            results.append(sub_data)

        if add_to is not None:
            results.append(out_data)

        if output_kernel_list:
            if klist_out != NULL:
                raw = klist_out.decode('utf-8')
                results.append(_parse_kernel_list(raw))
                free(klist_out)
            else:
                results.append({})

        return tuple(results) if len(results) > 2 else (out_data, out_mask)

