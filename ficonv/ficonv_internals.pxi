# ficonv_internals.pxi — ficonv kernel list parsing + ConvolutionOp
def parse_kernel_list(raw):
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

cdef class Ficonv_ConvolutionOp:
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
        if masked:                self._method |= Ficonv_ConvolutionOp.FM_MASKED
        if weighted:              self._method |= Ficonv_ConvolutionOp.FM_WEIGHTED
        if background_iterative:  self._method |= Ficonv_ConvolutionOp.FM_BGITERATIVE

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

        # output kernel arrays
        cdef int max_k = 32, max_c = 512
        cdef int [:] out_nk_arr = np.zeros(1, dtype=np.int32)
        cdef double [:] out_ox_arr = np.zeros(1, dtype=np.float64)
        cdef double [:] out_oy_arr = np.zeros(1, dtype=np.float64)
        cdef double [:] out_sc_arr = np.zeros(1, dtype=np.float64)
        cdef int [:] out_tp_arr = np.zeros(1, dtype=np.int32)
        cdef int [:] out_ct_arr = np.zeros(1, dtype=np.int32)
        cdef double [:] out_kc_arr = np.zeros(max_c, dtype=np.float64)
        cdef int [:] out_kt_arr = np.zeros(max_k, dtype=np.int32)
        cdef int [:] out_ko_arr = np.zeros(max_k, dtype=np.int32)
        cdef int [:] out_kn_arr = np.zeros(max_k, dtype=np.int32)
        cdef int [:] out_kh_arr = np.zeros(max_k, dtype=np.int32)
        cdef double [:] out_ks_arr = np.zeros(max_k, dtype=np.float64)
        cdef int [:] out_kx_arr = np.zeros(max_k, dtype=np.int32)
        cdef int [:] out_ky_arr = np.zeros(max_k, dtype=np.int32)

        fitsh_ficonv_fit_cy(&ir[0,0], <char*>&omr[0,0], &ii[0,0], <char*>&omi[0,0],
            inmask_ptr,
            sx, sy, self._kernel_spec,
            self._method, self._divide,
            self._iterations, self._rejection_level, self._gain,
            1 if self._verbose else 0,
            &od[0,0], <char*>&omsk[0,0],
            sub_ptr,
            add_ptr,
            self._unity_kernels,
            max_k,
            &out_nk_arr[0],
            &out_ox_arr[0], &out_oy_arr[0], &out_sc_arr[0],
            &out_tp_arr[0], &out_ct_arr[0],
            &out_kc_arr[0],
            &out_kt_arr[0], &out_ko_arr[0], &out_kn_arr[0],
            &out_kh_arr[0], &out_ks_arr[0],
            &out_kx_arr[0], &out_ky_arr[0],
            <char*>NULL, 0, 0, <char*>NULL,
            prefit_nk, pf_types, pf_orders, pf_ncoeffs, pf_coeffs,
            pf_ox, pf_oy, pf_scale)

        results = [out_data, out_mask]

        if output_subtracted:
            results.append(sub_data)

        if add_to is not None:
            results.append(out_data)

        if output_kernel_list:
            nk = out_nk_arr[0]
            kt = out_ct_arr[0]
            kd = {'global': {
                'type': int(out_tp_arr[0]),
                'offset': [float(out_ox_arr[0]), float(out_oy_arr[0])],
                'scale': float(out_sc_arr[0]),
            }, 'kernels': []}
            ci = 0
            for ki in range(nk):
                kn = out_kn_arr[ki]
                kern = {
                    'index': ki,
                    'type': int(out_kt_arr[ki]),
                    'order': int(out_ko_arr[ki]),
                }
                if out_kt_arr[ki] == 3:
                    kern['hsize'] = int(out_kh_arr[ki])
                    kern['sigma'] = float(out_ks_arr[ki])
                kern['bx'] = int(out_kx_arr[ki])
                kern['by'] = int(out_ky_arr[ki])
                kern['coeff'] = [float(out_kc_arr[ci + oo]) for oo in range(kn)]
                ci += kn
                kd['kernels'].append(kern)
            results.append(kd)

        return tuple(results) if len(results) > 2 else (out_data, out_mask)
