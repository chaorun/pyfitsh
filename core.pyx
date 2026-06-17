# cython: language_level=3
"""
grmatch_cy.pyx — Python wrapper around do_pointmatch_cy

Provides the Matcher class which takes flat coordinate arrays and
returns matched pairs, transformation coefficients, and statistics.
"""

from libc.stdlib cimport malloc, free
from libc.string cimport memset, memcpy
# from libc.stdio cimport FILE, fopen, fwrite, fclose  # dump helper
include "grmatch/grmatch_result.pxi"
include "wcs/wcs.pxi"
include "polyfit/polyfit.pxi"
include "ficonv/ficonv_internals.pxi"
include "fitrans/fitrans_internals.pxi"
include "fiphot/fiphot_internals.pxi"
include "firandom/firandom_internals.pxi"
include "algorithms/utils.pxi"


cdef class Grmatch:
    """
    Wrapper around the C do_pointmatch function.

    Usage:
        m = Grmatch(order=2, maxdist=2.0, unitarity=0.01)
        result = m.matchpoints(ref_x, ref_y, inp_x, inp_y, -rmag, -omag)
        print(result.hits_ref, result.vfits_dx, result.stats)
    """

    cdef double *ref_x
    cdef double *ref_y
    cdef double *ref_ord
    cdef double *inp_x
    cdef double *inp_y
    cdef double *inp_ord
    cdef int nref
    cdef int ninp

    # matchpointtune fields (defaults matching grmatch.c:79-89)
    cdef public int _order
    cdef public double _maxdist, _unitarity
    cdef public int _ttype, _parity, _use_ordering
    cdef public int _maxnum_ref, _maxnum_inp
    cdef public int _nmiter, _friter
    cdef public double _rejlevel
    cdef public int _wcat, _w_magnitude
    cdef public double _wpower
    cdef public int _is_centering
    cdef public double _refcx, _refcy, _inpcx, _inpcy, _maxcenterdist
    # offset/scale for output transformation
    cdef public double _ox, _oy, _scale
    # hint transformation
    cdef public list _hint_dxfit, _hint_dyfit
    cdef public int _hint_order
    cdef object _logger
    cdef bint _verbose

    def __cinit__(self):
        self.ref_x = NULL
        self.ref_y = NULL
        self.ref_ord = NULL
        self.inp_x = NULL
        self.inp_y = NULL
        self.inp_ord = NULL
        self.nref = 0
        self.ninp = 0
        self._logger = None
        self._verbose = False

    def __init__(self,
                 order=2,
                 maxdist=2.0,
                 unitarity=0.01,
                 ttype=2,
                 parity=0,
                 use_ordering=1,
                 maxnum_ref=0,
                 maxnum_inp=0,
                 nmiter=0,
                 friter=0,
                 rejlevel=3.0,
                 wcat=-1,
                 w_magnitude=0,
                 wpower=0.0,
                 is_centering=0,
                 refcx=0.0, refcy=0.0,
                 inpcx=0.0, inpcy=0.0,
                 maxcenterdist=0.0,
                 ox=0.0, oy=0.0, scale=1.0,
                 hint_dxfit=None, hint_dyfit=None, hint_order=0,
                 verbose=False,
                 logger=None):

        # store params as Python attributes for the run() call
        self._order = order
        self._maxdist = maxdist
        self._unitarity = unitarity
        self._ttype = ttype
        self._parity = parity
        self._use_ordering = use_ordering
        self._maxnum_ref = maxnum_ref
        self._maxnum_inp = maxnum_inp
        self._nmiter = nmiter
        self._friter = friter
        self._rejlevel = rejlevel
        self._wcat = wcat
        self._w_magnitude = w_magnitude
        self._wpower = wpower
        self._is_centering = is_centering
        self._refcx = refcx; self._refcy = refcy
        self._inpcx = inpcx; self._inpcy = inpcy
        self._maxcenterdist = maxcenterdist
        self._ox = ox; self._oy = oy; self._scale = scale
        self._hint_dxfit = hint_dxfit if hint_dxfit is not None else []
        self._hint_dyfit = hint_dyfit if hint_dyfit is not None else []
        self._hint_order = hint_order
        self._verbose = verbose
        self._logger = logger

    def matchpoints(self, ref_x, ref_y, inp_x, inp_y,
                    ref_order=None, inp_order=None):
        """Run point matching end-to-end. Returns MatchResult.
        
        ref_order/inp_order controls brightness sorting: provide
        negative magnitude (e.g. -mag) to sort brightest first.
        None skips internal sort.
        """
        cdef int i

        # set reference
        self.nref = len(ref_x)
        assert len(ref_y) == self.nref
        assert ref_order is None or len(ref_order) == self.nref
        if self.ref_x: free(self.ref_x)
        if self.ref_y: free(self.ref_y)
        if self.ref_ord: free(self.ref_ord)
        self.ref_x = <double*>malloc(self.nref * sizeof(double))
        self.ref_y = <double*>malloc(self.nref * sizeof(double))
        self.ref_ord = <double*>malloc(self.nref * sizeof(double))
        for i in range(self.nref):
            self.ref_x[i] = ref_x[i]
            self.ref_y[i] = ref_y[i]
            self.ref_ord[i] = ref_order[i] if ref_order is not None else 0.0
        if ref_order is None:
            self._use_ordering = 0

        # set input
        self.ninp = len(inp_x)
        assert len(inp_y) == self.ninp
        assert inp_order is None or len(inp_order) == self.ninp
        if self.inp_x: free(self.inp_x)
        if self.inp_y: free(self.inp_y)
        if self.inp_ord: free(self.inp_ord)
        self.inp_x = <double*>malloc(self.ninp * sizeof(double))
        self.inp_y = <double*>malloc(self.ninp * sizeof(double))
        self.inp_ord = <double*>malloc(self.ninp * sizeof(double))
        for i in range(self.ninp):
            self.inp_x[i] = inp_x[i]
            self.inp_y[i] = inp_y[i]
            self.inp_ord[i] = inp_order[i] if inp_order is not None else 0.0
        if inp_order is None:
            self._use_ordering = 0

        return self.run()

    def _log(self, msg, level='info'):
        """Internal: send message to logger if verbose=True or logger configured."""
        if not self._verbose:
            return
        if self._logger is None:
            print(msg)
            return
        if hasattr(self._logger, level):
            getattr(self._logger, level)(msg)
        elif callable(self._logger):
            self._logger(msg)

    cpdef MatchResult run(self):
        """Execute do_pointmatch and return a MatchResult."""
        if self.nref == 0 or self.ninp == 0:
            self._log("grmatch: no data to match", 'error')
            raise RuntimeError("No data set. Call matchpoints() first")
        self._log(f"grmatch: matching {self.nref} ref × {self.ninp} inp points")

        cdef:
            int *hits_idx0 = NULL
            int *hits_idx1 = NULL
            int nhit = 0
            double *vfits_dx = NULL
            double *vfits_dy = NULL
            matchpointstat mps
            int i, nvar, order

        memset(&mps, 0, sizeof(matchpointstat))
        order = self._order
        nvar = (order + 1) * (order + 2) // 2

        # Build hint arrays
        cdef double *hdx = NULL
        cdef double *hdy = NULL
        cdef int hnvar = 0
        if self._hint_dxfit and self._hint_dyfit and self._hint_order > 0:
            hnvar = (self._hint_order + 1) * (self._hint_order + 2) // 2
            hdx = <double*>malloc(hnvar * sizeof(double))
            hdy = <double*>malloc(hnvar * sizeof(double))
            for i in range(hnvar):
                hdx[i] = self._hint_dxfit[i]
                hdy[i] = self._hint_dyfit[i]

        do_pointmatch_cy(
            self.ref_x, self.ref_y, self.ref_ord, self.nref,
            self.inp_x, self.inp_y, self.inp_ord, self.ninp,
            self._ttype, self._maxdist, self._unitarity, self._parity,
            self._use_ordering, self._maxnum_ref, self._maxnum_inp,
            self._nmiter, self._friter, self._rejlevel,
            self._wcat, self._w_magnitude, self._wpower,
            self._is_centering, self._refcx, self._refcy,
            self._inpcx, self._inpcy, self._maxcenterdist,
            order,
            self._ox, self._oy, self._scale,
            hdx, hdy, self._hint_order,
            &hits_idx0, &hits_idx1, &nhit,
            &vfits_dx, &vfits_dy,
            &mps)

        if hdx: free(hdx)
        if hdy: free(hdy)

        self._log(f"grmatch: {nhit} pairs found, nsigma={mps.nsigma:.4f}")

        # Build result
        cdef MatchResult result = MatchResult()
        result.nhit = nhit
        result.order = order
        result.nvar = nvar

        result.hits_ref = [hits_idx0[i] for i in range(nhit)] if hits_idx0 else []
        result.hits_inp = [hits_idx1[i] for i in range(nhit)] if hits_idx1 else []

        # compute excluded indices (not in hits)
        cdef set _hset_ref = set(result.hits_ref)
        cdef set _hset_inp = set(result.hits_inp)
        result.excluded_ref = sorted([i for i in range(self.nref) if i not in _hset_ref])
        result.excluded_inp = sorted([i for i in range(self.ninp) if i not in _hset_inp])
        self._log(f"grmatch: {nhit} matched, {len(result.excluded_ref)}/{len(result.excluded_inp)} excluded ref/inp")

        result.vfits_dx = [vfits_dx[i] for i in range(nvar)] if vfits_dx else []
        result.vfits_dy = [vfits_dy[i] for i in range(nvar)] if vfits_dy else []

        result.stats = {
            'wsigma': mps.wsigma,
            'nsigma': mps.nsigma,
            'unitarity': mps.unitarity,
            'time_total': mps.time_total,
            'time_trimatch': mps.time_trimatch,
            'time_symmatch': mps.time_symmatch,
            'nmiter': mps.nmiter,
            'tri_level': mps.tri_level,
            'hull_coverage': mps.hull_coverage,
        }

        # free C allocations from do_pointmatch_cy
        if hits_idx0: free(hits_idx0)
        if hits_idx1: free(hits_idx1)
        # vfits_dx and vfits_dy are already freed by do_pointmatch_cy? 
        # No, they're allocated by calloc in do_pointmatch_cy and returned.
        # We need to free them too.
        if vfits_dx: free(vfits_dx)
        if vfits_dy: free(vfits_dy)

        return result

    def match_by_coords(self, ref_x, ref_y, inp_x, inp_y, maxdist=-1):
        """Simple nearest-neighbor coordinate matching. Returns (idx0, idx1) lists."""
        cdef int nref = len(ref_x), ninp = len(inp_x)
        cdef double *rx = <double*>malloc(nref * sizeof(double))
        cdef double *ry = <double*>malloc(nref * sizeof(double))
        cdef double *ix = <double*>malloc(ninp * sizeof(double))
        cdef double *iy = <double*>malloc(ninp * sizeof(double))
        cdef int *h0 = NULL, *h1 = NULL, nhit = 0
        cdef int i
        for i in range(nref): rx[i] = ref_x[i]; ry[i] = ref_y[i]
        for i in range(ninp): ix[i] = inp_x[i]; iy[i] = inp_y[i]

        coord_match_cy(rx, ry, nref, ix, iy, ninp, maxdist, &h0, &h1, &nhit)

        cdef list r0 = [h0[i] for i in range(nhit)] if h0 else []
        cdef list r1 = [h1[i] for i in range(nhit)] if h1 else []

        free(rx); free(ry); free(ix); free(iy)
        if h0: free(h0)
        if h1: free(h1)
        return r0, r1

    def match_by_id(self, ref_ids, inp_ids, ambiguity='none'):
        """Match by star identifier strings. Returns (idx0, idx1) lists.
        ambiguity: 'none'=require 1:1, 'first'=first match, 'any'=any, 'full'=cartesian product."""
        cdef int nref = len(ref_ids), ninp = len(inp_ids)
        cdef int ambig = 0  # AMBIG_NONE
        if ambiguity == 'first': ambig = 1
        elif ambiguity == 'any': ambig = 2
        elif ambiguity == 'full': ambig = 3

        cdef char **r_ids = <char**>malloc(nref * sizeof(char*))
        cdef char **i_ids = <char**>malloc(ninp * sizeof(char*))
        cdef int *h0 = NULL, *h1 = NULL, nhit = 0
        cdef int j
        cdef list ref_bytes = [], inp_bytes = []
        for j in range(nref):
            ref_bytes.append(str(ref_ids[j]).encode('utf-8'))
            r_ids[j] = <bytes>ref_bytes[j]
        for j in range(ninp):
            inp_bytes.append(str(inp_ids[j]).encode('utf-8'))
            i_ids[j] = <bytes>inp_bytes[j]

        id_match_cy(r_ids, nref, i_ids, ninp, ambig, &h0, &h1, &nhit)

        cdef list r0 = [h0[j] for j in range(nhit)] if h0 else []
        cdef list r1 = [h1[j] for j in range(nhit)] if h1 else []

        free(r_ids); free(i_ids)
        if h0: free(h0)
        if h1: free(h1)
        return r0, r1

    def __dealloc__(self):
        if self.ref_x: free(self.ref_x)
        if self.ref_y: free(self.ref_y)
        if self.ref_ord: free(self.ref_ord)
        if self.inp_x: free(self.inp_x)
        if self.inp_y: free(self.inp_y)
        if self.inp_ord: free(self.inp_ord)

    def to_dict(self):
        """Export Matcher parameters as JSON-serializable dict."""
        d = {
            'order': self._order, 'maxdist': self._maxdist,
            'unitarity': self._unitarity, 'ttype': self._ttype,
            'parity': self._parity, 'use_ordering': self._use_ordering,
            'rejlevel': self._rejlevel,
        }
        if self._ox != 0 or self._oy != 0 or self._scale != 1:
            d['ox'] = self._ox; d['oy'] = self._oy; d['scale'] = self._scale
        if self._hint_order > 0:
            d['hint_order'] = self._hint_order
            d['hint_dxfit'] = self._hint_dxfit
            d['hint_dyfit'] = self._hint_dyfit
        return d

    @classmethod
    def from_dict(cls, d):
        """Create Matcher from a dict of parameters."""
        return cls(**{k: d[k] for k in d
            if k in ('order','maxdist','unitarity','ttype','parity',
                     'use_ordering','rejlevel','maxnum_ref','maxnum_inp',
                     'nmiter','friter','wcat','w_magnitude','wpower',
                     'is_centering','refcx','refcy','inpcx','inpcy','maxcenterdist')})

cdef class Grtrans:
    """
    Apply a 2D polynomial transformation to arrays of (x,y) points.

    The transformation parameters typically come from a grmatch result
    (MatchResult.vfits_dx, .vfits_dy) or from a .trans file.

    Usage:
        t = Transformer(order=2, dxfit=vfits_dx, dyfit=vfits_dy, invert=True)
        new_x, new_y = t.apply(x_list, y_list)
    """

    cdef int _order
    cdef int _is_invert
    cdef double _ox, _oy, _scale
    cdef double *_dxfit
    cdef double *_dyfit
    cdef int _nvar

    property invert:
        """Get/set inversion flag."""
        def __get__(self): return bool(self._is_invert)
        def __set__(self, v): self._is_invert = 1 if v else 0

    def __cinit__(self):
        self._dxfit = NULL
        self._dyfit = NULL

    def __init__(self, order_or_dict=None, dxfit=None, dyfit=None,
                 ox=0.0, oy=0.0, scale=1.0, invert=False, order=None):
        # Accept either:
        #   Transformer(dict)   — from parse_trans_file / MatchResult.transformation / JSON
        #   Transformer(order=2, dxfit=..., dyfit=...) — explicit params
        #   Transformer(df, dxfit, dyfit) — positional (order, dxfit, dyfit)
        if dxfit is None and dyfit is None and isinstance(order_or_dict, dict):
            d = order_or_dict
            order = d['order']
            dxfit = d['dxfit']; dyfit = d['dyfit']
            offset = d.get('offset', [0.0, 0.0])
            if isinstance(offset, list):
                ox, oy = offset[0], offset[1]
            else:
                ox = d.get('ox', 0.0); oy = d.get('oy', 0.0)
            scale = d.get('scale', 1.0)
            invert = d.get('invert', invert)
        else:
            order = order_or_dict if order is None else order
            if order is None:
                raise TypeError('Transformer requires either a dict or order= parameter')

        self._order = order
        self._ox = ox; self._oy = oy; self._scale = scale
        self._is_invert = 1 if invert else 0
        self._nvar = (order + 1) * (order + 2) // 2
        assert len(dxfit) == self._nvar
        assert len(dyfit) == self._nvar

        self._dxfit = <double*>malloc(self._nvar * sizeof(double))
        self._dyfit = <double*>malloc(self._nvar * sizeof(double))
        cdef int i
        for i in range(self._nvar):
            self._dxfit[i] = dxfit[i]
            self._dyfit[i] = dyfit[i]

    def apply(self, x, y):
        """Apply transformation to (x, y) lists. Returns (new_x, new_y)."""
        cdef int n = len(x)
        assert len(y) == n
        cdef double *in_x  = <double*>malloc(n * sizeof(double))
        cdef double *in_y  = <double*>malloc(n * sizeof(double))
        cdef double *out_x = <double*>malloc(n * sizeof(double))
        cdef double *out_y = <double*>malloc(n * sizeof(double))
        cdef int i
        for i in range(n):
            in_x[i] = x[i]; in_y[i] = y[i]
        trans_apply_cy(in_x, in_y, n, out_x, out_y,
            self._order, self._is_invert, self._ox, self._oy, self._scale,
            self._dxfit, self._dyfit)
        rx = [out_x[i] for i in range(n)]
        ry = [out_y[i] for i in range(n)]
        free(in_x); free(in_y); free(out_x); free(out_y)
        return rx, ry

    def jacobian(self):
        """Return Jacobian matrices (jxx, jxy, jyx, jyy) as lists of coefficients.
        Useful for transforming ellipticity/shape parameters."""
        cdef int jnvar = (self._order + 0) * (self._order + 1) // 2
        cdef double *jxx = NULL
        cdef double *jxy = NULL
        cdef double *jyx = NULL
        cdef double *jyy = NULL
        trans_get_jacobi(self._order, self._scale, self._dxfit, self._dyfit,
                          &jxx, &jxy, &jyx, &jyy)
        result = (
            [jxx[i] for i in range(jnvar)] if jxx else [],
            [jxy[i] for i in range(jnvar)] if jxy else [],
            [jyx[i] for i in range(jnvar)] if jyx else [],
            [jyy[i] for i in range(jnvar)] if jyy else [],
        )
        free(jxx); free(jxy); free(jyx); free(jyy)
        return result

    def __dealloc__(self):
        if self._dxfit: free(self._dxfit)
        if self._dyfit: free(self._dyfit)

    def to_dict(self):
        """Export transformation as JSON-serializable dict."""
        return {
            'type': 'polynomial',
            'order': self._order,
            'offset': [self._ox, self._oy],
            'scale': self._scale,
            'dxfit': [self._dxfit[i] for i in range(self._nvar)],
            'dyfit': [self._dyfit[i] for i in range(self._nvar)],
        }

    @classmethod
    def from_dict(cls, d):
        """Create Transformer from a dict (e.g. parsed JSON or .trans)."""
        return cls(
            order=d['order'],
            dxfit=d['dxfit'], dyfit=d['dyfit'],
            ox=d.get('offset', [0.0, 0.0])[0],
            oy=d.get('offset', [0.0, 0.0])[1] if isinstance(d.get('offset', 0), list) else d.get('oy', 0.0),
            scale=d.get('scale', 1.0),
            invert=d.get('invert', False))

    @classmethod
    def WCS(cls, *args, **kwargs):
        return GrtransWCS(*args, **kwargs)

    @staticmethod
    def Proj(ra, dec, ra0, dec0, proj_type='tan'):
        """Project (RA,Dec) to tangent plane (xi,eta) in degrees."""
        return project_radec(ra, dec, ra0, dec0, {'sin':0,'arc':1,'tan':2}[proj_type])

    @staticmethod
    def Deproj(xi, eta, ra0, dec0, proj_type='tan'):
        """Deproject (xi,eta) back to (RA,Dec) in degrees."""
        return deproject_radec(xi, eta, ra0, dec0, {'sin':0,'arc':1,'tan':2}[proj_type])

    @classmethod
    def from_trans_string(cls, s):
        """Parse a .trans format string, return a Grtrans instance.
        The string follows the C grmatch --output-transformation format:
            type = <type>
            order = <order>
            offset = <ox>, <oy>
            scale = <scale>
            dxfit= <c0>, <c1>, ...
            dyfit= <c0>, <c1>, ...
        """
        params = {}
        for line in s.splitlines():
            line = line.split('#')[0].strip()
            if not line or '=' not in line: continue
            key, val = line.split('=', 1)
            key = key.strip(); val = val.strip()
            if key in ('dxfit', 'dyfit'):
                params[key] = [float(x) for x in val.split(',') if x.strip()]
            elif key == 'offset':
                params[key] = list(map(float, val.split(',')))
            else:
                try: params[key] = float(val) if '.' in val else int(val)
                except ValueError: params[key] = val
        return cls(
            order=params.get('order', 1),
            dxfit=params.get('dxfit', []), dyfit=params.get('dyfit', []),
            ox=params.get('offset', [0.0, 0.0])[0],
            oy=params.get('offset', [0.0, 0.0])[1] if isinstance(params.get('offset', 0), list) else params.get('oy', 0.0),
            scale=params.get('scale', 1.0),
            invert=False)

    def to_trans_string(self):
        """Serialize to .trans format string compatible with C grmatch."""
        d = self.to_dict()
        lines = [
            f"type = {d.get('type', 'polynomial')}",
            f"order = {d['order']}",
            f"offset = {d['offset'][0]}, {d['offset'][1]}",
            f"scale = {d.get('scale', 1.0)}",
            f"dxfit= {', '.join(str(x) for x in d['dxfit'])}",
            f"dyfit= {', '.join(str(x) for x in d['dyfit'])}",
        ]
        return '\n'.join(lines) + '\n'

    @staticmethod
    def compose_affine(dxfit, dyfit, offset_x, offset_y, scale=1.0):
        """Compose polynomial transformation with affine (bake offset/scale into coeffs)."""
        cdef int nvar = len(dxfit)
        cdef int order = (int((1 + 8 * nvar) ** 0.5) - 3) // 2
        cdef double *coeff_in  = <double*>malloc(nvar * sizeof(double))
        cdef double *rc_dx = <double*>malloc(nvar * sizeof(double))
        cdef double *rc_dy = <double*>malloc(nvar * sizeof(double))
        cdef int i
        for i in range(nvar): coeff_in[i] = dxfit[i]
        poly_compose_affine_cy(coeff_in, order,
            -offset_x / scale, 1.0 / scale, 0.0,
            -offset_y / scale, 0.0, 1.0 / scale, rc_dx)
        for i in range(nvar): coeff_in[i] = dyfit[i]
        poly_compose_affine_cy(coeff_in, order,
            -offset_x / scale, 1.0 / scale, 0.0,
            -offset_y / scale, 0.0, 1.0 / scale, rc_dy)
        result = ([rc_dx[i] for i in range(nvar)], [rc_dy[i] for i in range(nvar)])
        free(coeff_in); free(rc_dx); free(rc_dy)
        return result

class Fitrans:
    """Image transformation. All operations composed into a single pipeline.

    Usage:
        ff = Fitrans(transformation={...}, shift=(3,-2), flip=(True,False))
        result, mask = ff.apply(img)
    """
    def __init__(self, transformation=None,
                 shift=None, flip=None,
                 method='lanczos3', inverse=True,
                 offset=None, size=None,
                 zoom=None, shrink=None,
                 shrink_median=False, shrink_truncated_mean=0, shrink_optimistic=False,
                 interleave=None, interleave_median=False, interleave_optimistic=False,
                 repetitive=None,
                 smooth=None, noise=False):
        self._offset = offset
        self._size = size
        self._ops = []
        if noise:
            self._ops.append(('noise', Fitrans.Noise()))
        if smooth:
            self._ops.append(('smooth', Fitrans.Smooth(**smooth) if isinstance(smooth, dict) else smooth))
        if transformation or shift or flip:
            tf_kw = dict(method=method, inverse=inverse, transformation=transformation,
                         shift=shift, flip=flip)
            self._ops.append(('transform', Fitrans.ImageTransformation(**tf_kw)))
        if zoom:
            zx, zy, ox, oy = self._parse_xy(zoom)
            self._ops.append(('zoom', Fitrans.Zoom(zx, zy, ox, oy)))
        if shrink:
            sx, sy, ox, oy = self._parse_xy(shrink)
            self._ops.append(('shrink', Fitrans.Shrink(sx, sy, ox, oy,
                median=shrink_median, truncated_mean=shrink_truncated_mean,
                optimistic_mask=shrink_optimistic)))
        if interleave:
            ix, iy, ox, oy = self._parse_xy(interleave)
            self._ops.append(('interleave', Fitrans.Interleave(ix, iy, ox, oy,
                median=interleave_median, optimistic_mask=interleave_optimistic)))
        if repetitive:
            rx, ry = repetitive if isinstance(repetitive, tuple) else (repetitive, repetitive)
            self._ops.append(('repetitive', Fitrans.Repetitive(rx, ry)))

    @staticmethod
    def _parse_xy(v):
        if isinstance(v, (int, float)): return (int(v), int(v), 0, 0)
        if len(v) == 2: return (int(v[0]), int(v[1]), 0, 0)
        return (int(v[0]), int(v[1]), int(v[2]) if len(v)>2 else 0, int(v[3]) if len(v)>3 else 0)

    def apply(self, img_data, img_mask=None, weight=None, verbose=False):
        import numpy as np
        data = img_data.astype(np.float64, copy=False)
        mask = img_mask
        if isinstance(mask, list):
            mask = mask[0].copy()
            for m in mask[1:]: np.bitwise_or(mask, np.ascontiguousarray(m, dtype=np.uint8), out=mask)
        elif mask is not None:
            mask = np.ascontiguousarray(mask, dtype=np.uint8)
        else:
            mask = np.zeros_like(data, dtype=np.uint8)
        if weight is not None:
            data = data * np.ascontiguousarray(weight, dtype=np.float64)
        for name, op in self._ops:
            if name == 'transform' and (self._offset or self._size):
                ofx, ofy = (self._offset[0], self._offset[1]) if self._offset else (0, 0)
                data, mask = op.apply(data, mask, out_shape=self._size, ofx=ofx, ofy=ofy)
            else:
                data, mask = op.apply(data, mask)
        return data, mask

    # ---- Static helpers (kept for granular use) ----
    @staticmethod
    def ImageTransformation(order=1, method='lanczos3', inverse=True,
                             T=None, trans=None, transformation=None, shift=None, flip=None):
        """Factory for ImageTransform. inverse=True matches CLI --inverse."""
        t = T or trans or transformation
        tr = ImageTransform(order=order, method=method, inverse=inverse, T=t)
        if shift and (shift[0] or shift[1]):
            tr.compose_shift(shift[0], shift[1])
        if flip and (flip[0] or flip[1]):
            tr.compose_flip(flip[0], flip[1])
        return tr

    @staticmethod
    def Noise():
        """Noise level estimation via biquadratic scatter."""
        return NoiseOp()

    @staticmethod
    def Zoom(factor_x, factor_y=None, offset_x=0, offset_y=0, raw=False):
        """Integer zoom with biquadratic subpixel interpolation.
        Set raw=True for simple magnification without subpixel."""
        if factor_y is None: factor_y = factor_x
        return ZoomOp(factor_x, factor_y, offset_x, offset_y, raw)

    @staticmethod
    def Shrink(factor_x, factor_y=None, offset_x=0, offset_y=0,
               median=False, truncated_mean=0, optimistic_mask=False):
        """Integer shrink with flux conservation.
        median=True uses median averaging; truncated_mean=N rejects N lower/upper."""
        if factor_y is None: factor_y = factor_x
        return ShrinkOp(factor_x, factor_y, offset_x, offset_y,
                          median, truncated_mean, optimistic_mask)

    @staticmethod
    def Repetitive(nsx, nsy, offset_x=0, offset_y=0):
        """Repetitive expansion to absolute output size.
        Corresponds to --repetitive-xy <nsx>,<nsy> in fitrans CLI."""
        return RepetitiveOp(nsx, nsy, offset_x, offset_y)

    @staticmethod
    def Interleave(factor_x, factor_y, offset_x=0, offset_y=0,
                   median=False, optimistic_mask=False):
        """Interleave pixels: expand by interleaving original pixels.
        Corresponds to --interleave-xy <x>,<y> in fitrans CLI.
        median=True uses median; optimistic_mask=True for complement mask."""
        return InterleaveOp(factor_x, factor_y, offset_x, offset_y,
                              median, optimistic_mask)

    @staticmethod
    def Smooth(smooth_type='spline', order=2, xorder=None, yorder=None,
               prefilter=None, fxhsize=0, fyhsize=0, frejratio=0,
               niter=0, lower=3.0, upper=3.0,
               mean_unity=False, detrend=False):
        """Large-scale image smoothing via spline or polynomial fit.
        prefilter: 'mean' or 'median' for prefilter box filtering."""
        if xorder is None: xorder = order
        if yorder is None: yorder = order
        st = 1 if smooth_type == 'spline' else 2 if smooth_type == 'polynomial' else 0
        pf = 1 if prefilter == 'mean' else 2 if prefilter == 'median' else 0
        return SmoothOp(st, xorder, yorder, pf, fxhsize, fyhsize,
                          frejratio, niter, lower, upper, mean_unity, detrend)

    @staticmethod
    def Convolution(kernel=None, iterations=0, rejection_level=3.0,
                    masked=False, weighted=False, background_iterative=False,
                    divide=32, gain=1.0, verbose=False,
                    unity_kernels=False):
        """Kernel fitting + image convolution (ficonv).

        Parameters (matching ficonv --long-help keywords):
        kernel              : -k/--kernel       kernel specification string
        iterations          : -n/--iterations   number of rejection iterations
        rejection_level     : -s/--rejection-level  rejection sigma
        masked              : -m/--masked       masked fitting mode
        weighted            : -w/--weighted     weighted fitting (implies masked)
        background_iterative: -b/--bg-iterative background-iterative fitting
        divide              : -d/--divide       block divide factor
        gain                : -g/--gain         detector gain (e-/ADU)
        verbose             : --verbose         print progress
        unity_kernels       : -u/--unity-kernels normalize identity kernel to flux=1

        Returns a ConvolutionOp; call .fit(ref, img, ...) to execute.
        """
        return Ficonv_ConvolutionOp(kernel, iterations, rejection_level,
                              masked=masked, weighted=weighted,
                              background_iterative=background_iterative,
                              divide=divide, gain=gain,
                              verbose=verbose, unity_kernels=unity_kernels)

    @staticmethod
    def decode_maskinfo(hdu):
        """Decode MASKINFO from FITS HDU header into uint8 mask array."""
        from pyfitsh.utils import decode_maskinfo as _decode
        return _decode(hdu)

from pyfitsh.ficonv.ficonv import Ficonv  # noqa: F401
from pyfitsh.fiarith import Fiarith  # noqa: F401


cdef class Firandom:
    """Artificial image generator (firandom).

    Usage:
        fir = Firandom(sx=512, sy=512, gain=2.5, seed=42, background='100+5*x')
        result = fir.generate('100[x=256;y=256;i=5000;f=3.0;e=0.2;p=45]')
        img = result.image  # numpy 2D float64 array
    """
    cdef int _sx, _sy
    cdef double _gain
    cdef int _seed, _seed_noise, _seed_spatial, _seed_photon
    cdef bytes _bgexpr
    cdef int _zoom, _subg, _method
    cdef int _is_photnoise, _is_intinelect, _dontquantize
    cdef double _nsuppress
    # cdef starlistparam _dump_slp  # dump helper

    def __init__(self, sx=512, sy=512, gain=1.0, seed=0,
                 background=None,
                 zoom=1, subgrid=1, method=1,
                 is_photnoise=False, is_intinelect=False,
                 dontquantize=True, nsuppress=10000.0,
                 seed_noise=0, seed_spatial=0, seed_photon=0):
        self._sx = sx
        self._sy = sy
        self._gain = gain
        self._seed = seed
        self._seed_noise = seed_noise if seed_noise else seed
        self._seed_spatial = seed_spatial if seed_spatial else seed
        self._seed_photon = seed_photon if seed_photon else seed
        self._bgexpr = (background or "").encode('utf-8')
        self._zoom = zoom if zoom >= 1 else 1
        self._subg = subgrid if subgrid >= 1 else 1
        self._method = 0 if method == 0 else 1
        self._is_photnoise = 1 if is_photnoise else 0
        self._is_intinelect = 1 if is_intinelect else 0
        self._dontquantize = 1 if dontquantize else 0
        self._nsuppress = nsuppress

    cdef inline void _setup_sgp(self, stargenparam *sgp):
        memset(sgp, 0, sizeof(stargenparam))
        sgp.subg = self._subg
        sgp.subpixeldata = NULL
        sgp.gain = self._gain
        sgp.is_photnoise = self._is_photnoise
        sgp.method = self._method
        sgp.is_intinelect = self._is_intinelect
        sgp.dontquantize = self._dontquantize
        sgp.nsuppress = self._nsuppress
        sgp.tpd = NULL

    cdef inline int _alloc_image(self, image *img):
        cdef int i
        memset(img, 0, sizeof(image))
        img.sx = self._sx
        img.sy = self._sy
        img.data = <double**>malloc(self._sy * sizeof(double*))
        if img.data == NULL:
            return 1
        for i in range(self._sy):
            img.data[i] = <double*>malloc(self._sx * sizeof(double))
            if img.data[i] == NULL:
                return 1
            memset(img.data[i], 0, self._sx * sizeof(double))
        return 0

    cdef inline void _free_image(self, image *img):
        cdef int i
        if img.data != NULL:
            for i in range(img.sy):
                if img.data[i] != NULL:
                    free(img.data[i])
            free(img.data)
            img.data = NULL

    cdef inline void _img2np(self, image *img, double[:,::1] out):
        cdef int i, j
        for i in range(img.sy):
            for j in range(img.sx):
                out[i, j] = img.data[i][j]

    cdef inline int _make_stars_str(self, str spec, double ox, double oy, double scale,
                                      double mag, double flux,
                                      star **out_stars, int *out_n):
        cdef bytes buff = spec.encode('utf-8')
        cdef bytearray mb = bytearray(buff)
        cdef char *mb_ptr = <char*>mb
        cdef starlistparam lp
        cdef int ret
        replace_limiters(mb_ptr)
        memset(&lp, 0, sizeof(starlistparam))
        lp.mf0.magnitude = mag
        lp.mf0.intensity = flux
        lp.ox = ox; lp.oy = oy
        lp.scale = scale if scale != 0 else 1.0
        lp.sx = self._sx; lp.sy = self._sy
        lp.basetype = LISTTYPE_FEP
        # self._dump_slp = lp  # dump helper
        ret = create_input_list(mb_ptr, &lp, out_stars, out_n, self._seed_spatial)
        return ret

    cdef inline int _make_stars_list(self, stars_list, double ox, double oy, double scale,
                                       star **out_stars, int *out_n):
        cdef int n = len(stars_list)
        cdef int i, j
        cdef dict sd
        cdef double ff, fe, fp, fs, fd, fk
        cdef list mom_list
        cdef str model_str
        cdef star *stars_ptr

        if n == 0:
            out_stars[0] = NULL; out_n[0] = 0
            return 0

        stars_ptr = <star*>malloc(n * sizeof(star))
        if stars_ptr == NULL:
            return 1
        memset(stars_ptr, 0, n * sizeof(star))

        for i in range(n):
            sd = dict(stars_list[i])
            stars_ptr[i].location.gcx = float(sd.get('x', 0.0)) * scale + ox
            stars_ptr[i].location.gcy = float(sd.get('y', 0.0)) * scale + oy
            stars_ptr[i].flux = float(sd.get('flux', 0.0))
            stars_ptr[i].shape.model = SHAPE_ELLIPTIC

            if 's' in sd and 'd' in sd and 'k' in sd:
                stars_ptr[i].gsig = float(sd['s'])
                stars_ptr[i].gdel = float(sd['d'])
                stars_ptr[i].gkap = float(sd['k'])
            elif 'fwhm' in sd:
                ff = float(sd['fwhm'])
                fe = float(sd.get('ellip', 0.0))
                fp = float(sd.get('pa', 0.0))
                fep_to_sdk(ff, fe, fp, &fs, &fd, &fk)
                stars_ptr[i].gsig = fs
                stars_ptr[i].gdel = fd
                stars_ptr[i].gkap = fk
            else:
                stars_ptr[i].gsig = 1.0
                stars_ptr[i].gdel = 0.0
                stars_ptr[i].gkap = 0.0

            model_str = sd.get('model', 'elliptic').lower()
            if model_str == 'gauss':
                stars_ptr[i].shape.model = SHAPE_GAUSS
            elif model_str == 'deviated':
                stars_ptr[i].shape.model = SHAPE_DEVIATED
                stars_ptr[i].shape.order = int(sd.get('order', 2))
                stars_ptr[i].shape.gs = float(sd.get('gs', 1.0))
                mom_list = sd.get('mom', [])
                for j in range(min(len(mom_list), MAX_DEVIATION_COEFF)):
                    stars_ptr[i].shape.mom[j] = float(mom_list[j])

        out_stars[0] = stars_ptr; out_n[0] = n
        return 0

#     cdef void _dump_bytes(self, char *path, bytes data):
#         cdef int l = len(data)
#         cdef FILE *f = fopen(path, "wb")
#         fwrite(&l, 4, 1, f)
#         fwrite(<char*>data, 1, l, f)
#         fclose(f)
# 
#     cdef void _dump_struct(self, char *path, void *p, int sz):
#         cdef FILE *f = fopen(path, "wb")
#         fwrite(p, sz, 1, f)
#         fclose(f)
# 
#     cdef void _dump_stars(self, char *path, star *stars, int nstar, int seed):
#         cdef FILE *f = fopen(path, "wb")
#         cdef int i
#         fwrite(&nstar, 4, 1, f)
#         fwrite(&seed, 4, 1, f)
#         for i in range(nstar):
#             fwrite(&stars[i], sizeof(star), 1, f)
#         fclose(f)
# 
#     cdef void _dump_img(self, char *path, image *img):
#         cdef FILE *f = fopen(path, "wb")
#         cdef int i
#         fwrite(&img.sx, 4, 1, f)
#         fwrite(&img.sy, 4, 1, f)
#         for i in range(img.sy):
#             fwrite(img.data[i], 8, img.sx, f)
#         fclose(f)

    cdef inline int _prepare_stars(self, stars, double ox, double oy, double scale,
                                     double mag, double flux,
                                     star **out_stars, int *out_n):
        if stars is None:
            out_stars[0] = NULL; out_n[0] = 0
            return 0
        if isinstance(stars, str):
            return self._make_stars_str(stars, ox, oy, scale, mag, flux, out_stars, out_n)
        elif isinstance(stars, (list, tuple)):
            return self._make_stars_list(stars, ox, oy, scale, out_stars, out_n)
        else:
            raise TypeError(f"stars must be str or list of dicts, got {type(stars)}")

    def generate(self, stars=None, ox=0.0, oy=0.0, scale=1.0,
                 magnitude=10.0, intensity=10000.0,
                 background_stddev=0.0):
        cdef int ret, nstar
        cdef image img
        cdef star *cstars
        cdef stargenparam sgp
        cdef char *bg_ptr
        cdef double[:,::1] out_img

        import numpy as np

        if ox == 0.0 and oy == 0.0 and scale == 1.0:
            ox = 0.5 * self._sx
            oy = 0.5 * self._sy
            scale = 0.5 * self._sx

        # if isinstance(stars, str):
        #     self._dump_bytes("tmp/dump_cy/03_stararg.bin", stars.encode('utf-8'))

        nstar = 0
        cstars = NULL
        if stars is not None:
            ret = self._prepare_stars(stars, ox, oy, scale, magnitude, intensity, &cstars, &nstar)
            if ret != 0:
                raise RuntimeError("prepare_stars failed")

        # if cstars != NULL and nstar > 0:
        #     self._dump_struct("tmp/dump_cy/04_slp.bin", &self._dump_slp, sizeof(starlistparam))
        #     self._dump_stars("tmp/dump_cy/04_stars.bin", cstars, nstar, self._seed_spatial)

        random_seed(self._seed_noise)
        # np.int32(self._seed_noise).tofile("tmp/dump_cy/06_nseed.bin")

        ret = self._alloc_image(&img)
        if ret != 0:
            raise MemoryError("Failed to allocate image")

        bg_ptr = <char*>self._bgexpr if len(self._bgexpr) > 0 else NULL
        ret = create_background(&img, bg_ptr, background_stddev, ox, oy, scale, self._zoom)
        if ret != 0:
            self._free_image(&img)
            raise RuntimeError("create_background failed")

        # self._dump_img("tmp/dump_cy/07_bg.bin", &img)

        if self._is_intinelect:
            divide_image(&img, self._gain)

        if cstars != NULL and nstar > 0:
            random_seed(self._seed_photon)
            # np.int32(self._seed_photon).tofile("tmp/dump_cy/08_pseed.bin")
            self._setup_sgp(&sgp)
            draw_starlist(&img, &sgp, cstars, nstar, self._zoom)
            # self._dump_img("tmp/dump_cy/09_final.bin", &img)

        if not self._dontquantize:
            quantize_image(&img)

        out_img = np.zeros((self._sy, self._sx), dtype=np.float64)
        self._img2np(&img, out_img)
        self._free_image(&img)
        if cstars != NULL:
            free(cstars)

        import types
        result = types.SimpleNamespace()
        result.image = np.asarray(out_img)
        result.dict = {'sx': self._sx, 'sy': self._sy, 'nstar': nstar}
        return result

    def draw_stars(self, image_data, stars, ox=0.0, oy=0.0, scale=1.0,
                   magnitude=10.0, intensity=10000.0, zoom=None):
        cdef int sy, sx, i, nstar, z, ret
        cdef double[:,::1] img_v
        cdef double **row_ptrs
        cdef image img
        cdef star *cstars
        cdef stargenparam sgp

        import numpy as np

        z = self._zoom if zoom is None else zoom
        if z < 1:
            z = 1

        sy = image_data.shape[0]
        sx = image_data.shape[1]

        img_v = np.ascontiguousarray(image_data, dtype=np.float64)
        image_data = np.asarray(img_v)

        row_ptrs = <double**>malloc(sy * sizeof(double*))
        for i in range(sy):
            row_ptrs[i] = &img_v[i, 0]

        memset(&img, 0, sizeof(image))
        img.sx = sx; img.sy = sy; img.data = row_ptrs

        ret = self._prepare_stars(stars, ox, oy, scale, magnitude, intensity, &cstars, &nstar)
        if ret != 0:
            free(row_ptrs)
            raise RuntimeError("prepare_stars failed")
        if cstars != NULL and nstar > 0:
            self._setup_sgp(&sgp)
            draw_starlist(&img, &sgp, cstars, nstar, z)

        free(row_ptrs)
        if cstars != NULL:
            free(cstars)
        return image_data


cdef class Fistar:
    cdef double _th, _fth, _cp, _skysigma, _gain, _cand_rad
    cdef int _alg, _model, _morder, _it_sym, _it_gen
    cdef int _only_cand, _psf, _psf_hsize, _psf_grid, _psf_order
    cdef int _psf_type, _psf_symm, _psf_biquad
    cdef double _psf_ikappa, _psf_cwidth
    cdef int _psf_corder
    cdef int _out_mark, _out_area, _out_psf, _mark_sym, _mark_size
    cdef int _sort, _verbose
    cdef object _fields, collfit_str, _psf_str
    cdef object _cands, _poss, _magflux
    cdef int _src_xmin, _src_xmax, _src_ymin, _src_ymax

    def __init__(self, threshold=100.0, flux_threshold=0.0, critical_prominence=0.0,
                 skysigma=0.0, algorithm='uplink', model='elliptic', model_order=2,
                 it_sym=4, it_gen=2, only_candidates=False, collective_fit=None,
                 psf=None, gain=1.0,
                 input_candidates=None, candidate_radius=2.0,
                 input_positions=None, mag_flux=(10.0, 10000.0),
                 output_mark=True, output_area=True, psf_output=False,
                 mark_symbol='dot', mark_size=2,
                 sort='x', fields=None, section=None, verbose=False):
        self._th = threshold; self._fth = flux_threshold; self._cp = critical_prominence
        self._skysigma = skysigma; self._gain = gain
        self._alg = 1 if algorithm == 'parabolapeak' else 3
        self._model = {'gauss':1, 'elliptic':2, 'deviated':3}.get(model, 2)
        self._morder = model_order; self._it_sym = it_sym; self._it_gen = it_gen
        self._only_cand = 0 if only_candidates else 1
        self._cand_rad = candidate_radius
        self._psf = 1 if psf else 0
        self._psf_str = psf; self.collfit_str = collective_fit
        self._out_mark = 1 if output_mark else 0
        self._out_area = 1 if output_area else 0
        self._out_psf = 1 if psf_output else 0
        self._mark_sym = {'dot':0, 'square':1, 'circle':2}.get(mark_symbol, 0)
        self._mark_size = mark_size
        if section is not None:
            self._src_xmin, self._src_xmax, self._src_ymin, self._src_ymax = section
        else:
            self._src_xmin = self._src_xmax = self._src_ymin = self._src_ymax = 0
        s_map = {'x':0,'y':1,'peak':2,'fwhm':3,'amp':4,'flux':5,'noise':6,'sn':7}
        self._sort = s_map.get(sort, 0)
        self._fields = fields
        self._cands = input_candidates
        self._poss = input_positions
        self._magflux = mag_flux
        self._verbose = verbose

        # parse PSF string
        self._psf_hsize = 4; self._psf_grid = 4; self._psf_order = 0
        self._psf_type = 0; self._psf_symm = 0; self._psf_biquad = 0
        self._psf_ikappa = 0.0; self._psf_cwidth = 0.0; self._psf_corder = 0
        if psf and isinstance(psf, str):
            for kv in psf.split(','):
                k, _, v = kv.partition('=')
                k = k.strip()
                if k == 'native': self._psf_type = 0
                elif k == 'integral': self._psf_type = 1
                elif k == 'circle': self._psf_type = 2
                elif k == 'biquad' or k == 'spline': self._psf_type = 0; self._psf_biquad = 1
                elif k == 'order': self._psf_order = int(v)
                elif k == 'grid': self._psf_grid = int(v)
                elif k == 'halfsize': self._psf_hsize = int(v)
                elif k == 'symmetrize': self._psf_symm = 1
                elif k == 'kappa': self._psf_ikappa = float(v)
                elif k == 'circlewidth' or k == 'width': self._psf_cwidth = float(v)
                elif k == 'circleorder': self._psf_corder = int(v)

    def do_fistar(self, img_data, mask=None, maskinfo=None):
        """Run star detection. Returns dict with numpy arrays for all 38 format fields.

        Parameters:
            maskinfo: str or list of str — MASKINFO header strings from a FITS file.
                      Parsed via parse_maskinfo() and merged into the mask.
        """
        import numpy as np
        sy, sx = img_data.shape
        in_img = np.ascontiguousarray(img_data, dtype=np.float64)
        if mask is not None:
            mk_arr = np.ascontiguousarray(mask, dtype=np.uint8)
        else:
            mk_arr = np.zeros((sy, sx), dtype=np.uint8)
        if maskinfo is not None:
            mi_mask = parse_maskinfo(maskinfo, sx, sy)
            mk_arr = np.bitwise_or(mk_arr, mi_mask)
            mk_arr = np.zeros((sy, sx), dtype=np.uint8)

        cdef double [:,::1] img_view = in_img
        cdef unsigned char [:,::1] mk_view = mk_arr

        cdef int nstar = 0
        cdef fistar_result res
        memset(&res, 0, sizeof(fistar_result))
        cdef double *psf_data = NULL
        cdef int psf_nx = 0, psf_ny = 0, psf_nvar = 0

        cdef double th = self._th
        # threshold zeroing is handled in fistar_core.c (ALG_LNK branch)

        # allocate output buffers for mark/area
        cdef double [:,::1] mark_view, area_view
        cdef double *mark_ptr = NULL, *area_ptr = NULL
        if self._out_mark:
            out_mark = np.empty((sy, sx), dtype=np.float64)
            mark_view = out_mark; mark_ptr = &mark_view[0,0]
        if self._out_area:
            out_area = np.empty((sy, sx), dtype=np.float64)
            area_view = out_area; area_ptr = &area_view[0,0]

        # input_candidates: numpy (N,2) or (N,3) → flat double*
        cdef double *cand_ptr = NULL
        cdef int nc = 0, nc_cols = 2
        cdef double [:,::1] cand_view
        if self._cands is not None:
            cand_arr = np.ascontiguousarray(self._cands, dtype=np.float64)
            nc = cand_arr.shape[0]
            nc_cols = cand_arr.shape[1] if cand_arr.ndim > 1 else 1
            cand_view = cand_arr; cand_ptr = &cand_view[0,0]

        # input_positions: list of (x,y) → flat double*
        cdef double *pos_ptr = NULL
        cdef int npos = 0
        cdef double [:] pos_view
        cdef double [:] pos_flux_view
        cdef double [:] pos_ferr_view
        cdef double *pflux = NULL, *pferr = NULL
        if self._poss is not None:
            npos = len(self._poss)
            pos_xy = np.array(self._poss, dtype=np.float64).ravel()
            pos_view = pos_xy
            pos_ptr = &pos_view[0]
            pos_flux_out = np.zeros(npos, dtype=np.float64)
            pos_ferr_out = np.zeros(npos, dtype=np.float64)
            pos_flux_view = pos_flux_out; pos_ferr_view = pos_ferr_out
            pflux = &pos_flux_view[0]; pferr = &pos_ferr_view[0]

        # PSF output
        cdef fistar_psf_out pso
        memset(&pso, 0, sizeof(pso))
        cdef fistar_psf_out *psoptr = &pso if self._out_psf else NULL

        cdef int collfit_niter = -1
        cdef int collfit_refinelevel = 0
        cdef double collfit_bhsize = 10.0
        cdef int collfit_niter_raw
        if self.collfit_str:
            for tok in self.collfit_str.split(','):
                kv = tok.strip()
                if kv.startswith('iterations='):
                    collfit_niter_raw = int(kv.split('=')[1])
                    collfit_niter = collfit_niter_raw + 1
                elif kv == 'position':
                    collfit_refinelevel |= 1
                elif kv == 'shape':
                    collfit_refinelevel |= 2
                elif kv.startswith('bhsize=') or kv.startswith('blockhalfsize='):
                    collfit_bhsize = float(kv.split('=')[1])
                elif kv.isdigit():
                    collfit_niter_raw = int(kv)
                    collfit_niter = collfit_niter_raw + 1

        fistar_search_cy(&img_view[0,0], <char*>&mk_view[0,0], sx, sy,
            th, self._fth, self._cp, self._skysigma,
            self._alg, self._model, self._morder,
            self._it_sym, self._it_gen, self._only_cand,
            collfit_niter, collfit_refinelevel, collfit_bhsize,
            self._psf, self._psf_hsize, self._psf_grid, self._psf_order,
            self._psf_type, self._psf_symm, self._psf_biquad,
            self._psf_ikappa, self._psf_cwidth, self._psf_corder,
            self._gain,
            self._magflux[1], self._magflux[0],   # mag_intensity, mag_magnitude
            cand_ptr, nc, nc_cols, self._cand_rad,
            pos_ptr, npos,
            self._src_xmin, self._src_xmax, self._src_ymin, self._src_ymax,
            self._mark_sym, self._mark_size,
            (1 if self._out_mark else 0) | (2 if self._out_area else 0) | (4 if self._out_psf else 0),
            self._sort,
            1 if self._verbose else 0,
            &nstar, &res,
            mark_ptr, area_ptr,
            psoptr,
            pflux, pferr)

        if self._fields: wanted = set(self._fields.split(','))
        else: wanted = {'id','ix','iy','cx','cy','cbg','camp','cmax','npix','cs','cd','ck',
                        'x','y','bg','amp','s','d','k','l','mom',
                        'sigma','delta','kappa','fwhm','ellip','pa',
                        'flux','noise','sn','magnitude',
                        'px','py','pbg','pamp','ps','pd','pk','pl'}

        # map aliases
        if 'mag' in wanted: wanted.discard('mag'); wanted.add('magnitude')
        if 's/n' in wanted: wanted.discard('s/n'); wanted.add('sn')

        result = {}
        if nstar > 0:
            if 'id' in wanted: result['id'] = np.array([res.id[i] for i in range(nstar)], dtype=np.int32)
            if 'ix' in wanted: result['ix'] = np.array([res.ix[i] for i in range(nstar)])
            if 'iy' in wanted: result['iy'] = np.array([res.iy[i] for i in range(nstar)])
            if 'cx' in wanted: result['cx'] = np.array([res.cx[i] for i in range(nstar)])
            if 'cy' in wanted: result['cy'] = np.array([res.cy[i] for i in range(nstar)])
            if 'cbg' in wanted: result['cbg'] = np.array([res.cbg[i] for i in range(nstar)])
            if 'camp' in wanted: result['camp'] = np.array([res.camp[i] for i in range(nstar)])
            if 'cmax' in wanted: result['cmax'] = np.array([res.cmax[i] for i in range(nstar)])
            if 'npix' in wanted: result['npix'] = np.array([res.npix[i] for i in range(nstar)], dtype=np.int32)
            if 'cs' in wanted: result['cs'] = np.array([res.cs[i] for i in range(nstar)])
            if 'cd' in wanted: result['cd'] = np.array([res.cd[i] for i in range(nstar)])
            if 'ck' in wanted: result['ck'] = np.array([res.ck[i] for i in range(nstar)])
            if 'x' in wanted: result['x'] = np.array([res.x[i] for i in range(nstar)])
            if 'y' in wanted: result['y'] = np.array([res.y[i] for i in range(nstar)])
            if 'bg' in wanted: result['bg'] = np.array([res.bg[i] for i in range(nstar)])
            if 'amp' in wanted: result['amp'] = np.array([res.amp[i] for i in range(nstar)])
            if 's' in wanted: result['s'] = np.array([res.s[i] for i in range(nstar)])
            if 'd' in wanted: result['d'] = np.array([res.d[i] for i in range(nstar)])
            if 'k' in wanted: result['k'] = np.array([res.k[i] for i in range(nstar)])
            if 'l' in wanted: result['l'] = np.array([res.l[i] for i in range(nstar)])
            if 'mom' in wanted: result['mom'] = np.array([res.mom[i * 15] for i in range(nstar)])
            if 'sigma' in wanted: result['sigma'] = np.array([res.sigma[i] for i in range(nstar)])
            if 'delta' in wanted: result['delta'] = np.array([res.delta[i] for i in range(nstar)])
            if 'kappa' in wanted: result['kappa'] = np.array([res.kappa[i] for i in range(nstar)])
            if 'fwhm' in wanted: result['fwhm'] = np.array([res.fwhm[i] for i in range(nstar)])
            if 'ellip' in wanted: result['ellip'] = np.array([res.ellip[i] for i in range(nstar)])
            if 'pa' in wanted: result['pa'] = np.array([res.pa[i] for i in range(nstar)])
            if 'flux' in wanted: result['flux'] = np.array([res.flux[i] for i in range(nstar)])
            if 'noise' in wanted: result['noise'] = np.array([res.noise[i] for i in range(nstar)])
            if 'sn' in wanted: result['sn'] = np.array([res.sn[i] for i in range(nstar)])
            if 'magnitude' in wanted: result['magnitude'] = np.array([res.magnitude[i] for i in range(nstar)])
            if 'px' in wanted: result['px'] = np.array([res.px[i] for i in range(nstar)])
            if 'py' in wanted: result['py'] = np.array([res.py[i] for i in range(nstar)])
            if 'pbg' in wanted: result['pbg'] = np.array([res.pbg[i] for i in range(nstar)])
            if 'pamp' in wanted: result['pamp'] = np.array([res.pamp[i] for i in range(nstar)])
            if 'ps' in wanted: result['ps'] = np.array([res.ps[i] for i in range(nstar)])
            if 'pd' in wanted: result['pd'] = np.array([res.pd[i] for i in range(nstar)])
            if 'pk' in wanted: result['pk'] = np.array([res.pk[i] for i in range(nstar)])
            if 'pl' in wanted: result['pl'] = np.array([res.pl[i] for i in range(nstar)])

        from astropy.table import Table
        table = Table(result)
        out_dict = {'table': table, 'nstar': nstar}
        if self._out_mark:
            out_dict['output_mark'] = out_mark
        if self._out_area:
            out_dict['output_area'] = out_area
        if self._out_psf and pso.data != NULL:
            nv = pso.nvar; ns = pso.nside
            psf_cube = np.array([pso.data[i] for i in range(nv*ns*ns)], dtype=np.float32).reshape(nv, ns, ns)
            from astropy.io import fits
            hdr = fits.Header()
            hdr['PSFHSIZE'] = pso.hsize; hdr['PSFSGRID'] = pso.grid; hdr['PSFORDER'] = pso.order
            hdr['PSFOFFSX'] = pso.ox; hdr['PSFOFFSY'] = pso.oy; hdr['PSFSCALE'] = pso.scale
            hdr['ORIGIN'] = 'pyfitsh/fistar_pipeline'
            psf_hdu = fits.ImageHDU(data=psf_cube, header=hdr)
            out_dict['psf'] = psf_hdu
            if pso.data: free(pso.data)
        if self._poss is not None:
            m0 = self._magflux[0]
            f0 = self._magflux[1]
            pos_mag = np.where(pos_flux_out > 0, m0 - 2.5 * np.log10(pos_flux_out / f0), 0.0)
            pos_merr = np.where(pos_flux_out > 0,
                                np.abs(-2.5 * pos_ferr_out / (pos_flux_out * np.log(10))), 0.0)
            pos_x = [p[0] for p in self._poss]
            pos_y = [p[1] for p in self._poss]
            pos_table = Table({'input_id': range(len(self._poss)),
                               'input_x': pos_x, 'input_y': pos_y,
                               'mag': pos_mag, 'merr': pos_merr,
                               'flux': pos_flux_out, 'ferr': pos_ferr_out})
            out_dict['table'] = pos_table

        free(res.id); free(res.ix); free(res.iy); free(res.cx); free(res.cy)
        free(res.cbg); free(res.camp); free(res.cmax)
        if res.npix: free(res.npix)
        free(res.cs); free(res.cd); free(res.ck)
        free(res.x); free(res.y); free(res.bg); free(res.amp)
        free(res.s); free(res.d); free(res.k); free(res.l); free(res.mom)
        free(res.sigma); free(res.delta); free(res.kappa)
        free(res.fwhm); free(res.ellip); free(res.pa)
        free(res.flux); free(res.noise); free(res.sn); free(res.magnitude)
        free(res.px); free(res.py); free(res.pbg); free(res.pamp)
        free(res.ps); free(res.pd); free(res.pk); free(res.pl)

        import types
        ns = types.SimpleNamespace()
        ns.output = out_dict.pop('table')
        ns.nstar = out_dict.pop('nstar')
        ns.output_mark = out_dict.pop('output_mark', None)
        ns.output_area = out_dict.pop('output_area', None)
        ns.psf = out_dict.pop('psf', None)
        ns.dict = {'output': ns.output, 'nstar': ns.nstar}
        if ns.output_mark is not None: ns.dict['output_mark'] = ns.output_mark
        if ns.output_area is not None: ns.dict['output_area'] = ns.output_area
        if ns.psf is not None: ns.dict['psf'] = ns.psf
        return ns

# ========================================================================
# Pure-Python helpers for WCS projection (mirrors projection.c)
# ========================================================================

class Fiphot:
    """Aperture photometry (fiphot).

    Parameters (matching fiphot --long-help):
      apertures     : str or list — e.g. '4.0:6.0:2.0' or ['4.0:6.0:2.0','2.0:3.0:1.0']
      sky_fit       : str — 'mean', 'median', 'mode,iterations=3,lower=3.0,upper=3.0'
      gain          : str — '1.0' or '0.9,-0.01,0.02' (polynomial)
      gain_vmin     : float — minimal gain value
      fwhm          : float — FWHM for optimal aperture
      skynoise      : float — sky noise for optimal aperture
      mag_flux      : (mag, flux) — magnitude zero point
      correlation_length : float — background correlation length
      mask_ignore   : str — e.g. 'saturated,hot' or int mask flags
      spline        : bool — use biquadratic spline interpolation (-k)
      disjoint_rings : bool — disjoint background annuli (-j)
      disjoint_apertures : bool — disjoint apertures (-p)
      disjoint_radius : float — exclusion radius around other centroids (-x)
      zoom          : int — multiply coordinates and radii
      serial        : str — serial identifier
    """
    def __init__(self, apertures=None, sky_fit=None,
                 gain='1.0', gain_vmin=0.0, fwhm=5.0, skynoise=0.0,
                 mag_flux=(10.0, 10000.0), correlation_length=1.0,
                 mask_ignore=0, spline=False,
                 disjoint_rings=False, disjoint_apertures=False, disjoint_radius=-1.0,
                 zoom=1, serial=''):
        self._apertures = apertures
        self._sky_fit = sky_fit
        self._gain = gain
        self._gain_vmin = gain_vmin
        self._fwhm = fwhm
        self._skynoise = skynoise
        self._mag_flux = mag_flux
        self._correlation_length = correlation_length
        self._mask_ignore = mask_ignore
        self._spline = spline
        self._disjoint_rings = disjoint_rings
        self._disjoint_apertures = disjoint_apertures
        self._disjoint_radius = disjoint_radius
        self._zoom = zoom
        self._serial = serial

    def photometry(self, img_data, mask=None, stars=None,
                   col_xy=(0, 1), col_ap=None,
                   col_id=-1, col_mag=-1, col_col=-1, col_err=-1,
                   calc_optimal=False):
        """Perform aperture photometry.

        Args:
            img_data: float64 ndarray (2D)
            mask: uint8 ndarray (optional)
            stars: ndarray shape (N, M) — star positions + optional columns.
                   Columns indexed by col_xy, col_ap, col_id, col_mag, col_col, col_err.
                   Can also be a list of (x,y) or (x,y,r0,ra,da,...) tuples.
            col_xy: (col_x, col_y) — default (0, 1)
            col_ap: list of column indices containing 'r0:ra:da' aperture specs
            col_id: column for object identifier (-1 = none)
            col_mag: column for reference magnitude (-1 = none)
            col_col: column for photometric color (-1 = none)
            col_err: column for magnitude uncertainty (-1 = none)
            calc_optimal: compute optimal aperture radius

        Returns:
            dict with keys: flux, fluxerr, magnitude, magerr, bgarea, bgflux,
            bgmedian, bgsigma, cntr_x, cntr_y, cntr_width, cntr_w_d, cntr_w_k,
            cntr_x_err, cntr_y_err, cntr_w_err, flag, rtot, rbad, rign, atot, abad,
            optimal_r0, optimal_ra, optimal_da (if calc_optimal=True),
            id (if col_id >= 0).
            Each per-star/per-aperture value is a 2D ndarray shape (nstar, nap).
            id is a 1D array shape (nstar,).
        """
        import numpy as np
        cdef int i, j

        sy, sx = img_data.shape
        in_img = np.ascontiguousarray(img_data, dtype=np.float64)
        if mask is not None:
            mk_arr = np.ascontiguousarray(mask, dtype=np.uint8)
        else:
            mk_arr = np.zeros((sy, sx), dtype=np.uint8)

        cdef double [:,::1] img_view = in_img
        cdef unsigned char [:,::1] mk_view = mk_arr

        # Parse stars — support both ndarray and list-of-tuples
        cdef int nstar, nap, has_per_star_aps
        cdef double [:] sx_arr, sy_arr
        cdef double [:] sr0_arr, sra_arr, sda_arr
        cdef double [:] srefmag_arr, srefcol_arr, sreferr_arr
        cdef list ids_list

        if stars is None:
            raise ValueError("stars is required")

        # Convert to ndarray if list/tuple
        cdef double [:,::1] star_mat
        cdef bint is_ndarray = isinstance(stars, np.ndarray)
        if is_ndarray:
            # ndarray path — use column indices
            star_mat = np.ascontiguousarray(stars, dtype=np.float64)
            nstar = star_mat.shape[0]
            ids_list = None
            cx, cy = col_xy[0], col_xy[1]

            sx_arr = np.asarray(star_mat[:, cx]) * self._zoom
            sy_arr = np.asarray(star_mat[:, cy]) * self._zoom

            # Per-star apertures from columns
            if col_ap is not None and len(col_ap) > 0:
                has_per_star_aps = True
                nap = len(col_ap)
                sr0_arr = np.zeros(nstar * nap, dtype=np.float64)
                sra_arr = np.zeros(nstar * nap, dtype=np.float64)
                sda_arr = np.zeros(nstar * nap, dtype=np.float64)
                for j, acol in enumerate(col_ap):
                    for i in range(nstar):
                        apspec = str(stars[i, acol])
                        parts = apspec.replace(',', ':').split(':')
                        idx = i * nap + j
                        if len(parts) >= 3:
                            sr0_arr[idx] = float(parts[0]) * self._zoom
                            sra_arr[idx] = float(parts[1]) * self._zoom
                            sda_arr[idx] = float(parts[2]) * self._zoom
                        elif len(parts) == 1:
                            sr0_arr[idx] = float(parts[0]) * self._zoom
                            sra_arr[idx] = 2.0 * self._zoom
                            sda_arr[idx] = 1.0 * self._zoom
            else:
                has_per_star_aps = False
                nap = 0
                sr0_arr = np.zeros(0, dtype=np.float64)
                sra_arr = np.zeros(0, dtype=np.float64)
                sda_arr = np.zeros(0, dtype=np.float64)

            # Reference mag / color / err from columns
            srefmag_arr = np.zeros(nstar, dtype=np.float64)
            srefcol_arr = np.zeros(nstar, dtype=np.float64)
            sreferr_arr = np.zeros(nstar, dtype=np.float64)
            col_id_used = col_id
            if col_mag >= 0 and col_mag < star_mat.shape[1]:
                srefmag_arr = np.asarray(star_mat[:, col_mag])
            if col_col >= 0 and col_col < star_mat.shape[1]:
                srefcol_arr = np.asarray(star_mat[:, col_col])
            if col_err >= 0 and col_err < star_mat.shape[1]:
                sreferr_arr = np.asarray(star_mat[:, col_err])
            if col_id >= 0 and col_id < star_mat.shape[1]:
                ids_list = [stars[i, col_id] for i in range(nstar)]

        else:
            # list-of-tuples path (backward compat)
            nstar = len(stars)
            first = stars[0]
            if isinstance(first, (list, tuple)) and len(first) >= 5:
                has_per_star_aps = True
                nap = (len(first) - 2) // 3
            else:
                has_per_star_aps = False
                nap = 0

            sx_arr = np.zeros(nstar, dtype=np.float64)
            sy_arr = np.zeros(nstar, dtype=np.float64)
            srefmag_arr = np.zeros(nstar, dtype=np.float64)
            srefcol_arr = np.zeros(nstar, dtype=np.float64)
            sreferr_arr = np.zeros(nstar, dtype=np.float64)

            if has_per_star_aps:
                sr0_arr = np.zeros(nstar * nap, dtype=np.float64)
                sra_arr = np.zeros(nstar * nap, dtype=np.float64)
                sda_arr = np.zeros(nstar * nap, dtype=np.float64)
            else:
                sr0_arr = np.zeros(0, dtype=np.float64)
                sra_arr = np.zeros(0, dtype=np.float64)
                sda_arr = np.zeros(0, dtype=np.float64)

            for i in range(nstar):
                s = stars[i]
                sx_arr[i] = s[0] * self._zoom
                sy_arr[i] = s[1] * self._zoom
                if has_per_star_aps:
                    for j in range(nap):
                        idx = i * nap + j
                        sr0_arr[idx] = s[2 + 3*j] * self._zoom
                        sra_arr[idx] = s[3 + 3*j] * self._zoom
                        sda_arr[idx] = s[4 + 3*j] * self._zoom
            ids_list = None
            col_id_used = -1

        # Build aperture spec from global apertures
        cdef bytes ap_spec_b = b''
        if not has_per_star_aps and self._apertures is not None:
            if isinstance(self._apertures, (list, tuple)):
                ap_spec = ','.join(str(a) for a in self._apertures)
            else:
                ap_spec = str(self._apertures)
            ap_spec_b = ap_spec.encode('utf-8')
            nap = ap_spec.count(',') + 1
        elif not has_per_star_aps:
            raise ValueError("Either provide per-star apertures (col_ap or 5+ element tuples) or global apertures")

        cdef char *ap_spec_ptr = NULL
        if len(ap_spec_b) > 0:
            ap_spec_ptr = ap_spec_b

        # Parse sky_fit / background params (unchanged)
        cdef int bg_type = 2
        cdef int bg_scatter = 0
        cdef int bg_rejniter = 0
        cdef double bg_rejlower = 3.0
        cdef double bg_rejupper = 3.0
        cdef int use_sky = 0
        cdef double sky_level = 0.0
        if self._sky_fit is not None:
            for tok in self._sky_fit.replace(',', ' ').split():
                if tok == 'mean': bg_type = 1
                elif tok == 'median': bg_type = 2
                elif tok == 'mode': bg_type = 3
                elif tok == 'mad': bg_scatter = 1
                elif tok == 'stddev': bg_scatter = 0
                elif tok.startswith('iterations='): bg_rejniter = int(tok.split('=')[1])
                elif tok.startswith('lower='): bg_rejlower = float(tok.split('=')[1])
                elif tok.startswith('upper='): bg_rejupper = float(tok.split('=')[1])
                elif tok.startswith('sigma='):
                    s = float(tok.split('=')[1])
                    bg_rejlower = bg_rejupper = s
                elif tok.startswith('force='):
                    sky_level = float(tok.split('=')[1])
                    use_sky = 1

        cdef int mask_ignore = 0
        if isinstance(self._mask_ignore, str) and self._mask_ignore:
            mask_map = {'saturated': 0x02, 'hot': 0x04, 'outer': 0x01}
            for tok in self._mask_ignore.replace(',', ' ').split():
                mask_ignore |= mask_map.get(tok.lower(), 0)
        elif isinstance(self._mask_ignore, int):
            mask_ignore = self._mask_ignore

        cdef int sg_order = 0
        cdef double sg_vmin = self._gain_vmin
        cdef double [:] sg_coeff_view
        cdef double *sg_coeff_ptr = NULL
        parts = [float(x) for x in self._gain.replace(',', ' ').split()]
        sg_order = int((1 + 8 * len(parts)) ** 0.5 - 3) // 2
        ncoeff = (sg_order + 1) * (sg_order + 2) // 2
        sg_coeff_arr = np.array(parts[:ncoeff], dtype=np.float64)
        sg_coeff_view = sg_coeff_arr
        sg_coeff_ptr = &sg_coeff_view[0]

        # Output arrays
        cdef double [:,::1] out_flux_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_fluxerr_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgarea_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgflux_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgmedian_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgsigma_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_x_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_y_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_width_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_d_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_k_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_x_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_y_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef int [:,::1] out_flag_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rtot_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rbad_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rign_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_atot_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_abad_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef double [:] out_opt_r0_arr = np.zeros(nstar, dtype=np.float64)
        cdef double [:] out_opt_ra_arr = np.zeros(nstar, dtype=np.float64)
        cdef double [:] out_opt_da_arr = np.zeros(nstar, dtype=np.float64)

        cdef double [:] out_raw_arr = np.zeros(nstar * max(nap,1) * 12, dtype=np.float64)

        cdef double sigma = self._fwhm / 2.354820045

        cdef double *sr0_ptr = NULL
        cdef double *sra_ptr = NULL
        cdef double *sda_ptr = NULL
        if has_per_star_aps and nap > 0:
            sr0_ptr = &sr0_arr[0]
            sra_ptr = &sra_arr[0]
            sda_ptr = &sda_arr[0]

        fiphot_photometry_cy(
            &img_view[0,0], <char*>&mk_view[0,0], sx, sy,
            &sx_arr[0], &sy_arr[0], nstar,
            sr0_ptr, sra_ptr, sda_ptr,
            nap if has_per_star_aps else 0,
            ap_spec_ptr, self._zoom,
            bg_type, bg_scatter, bg_rejniter, bg_rejlower, bg_rejupper,
            mask_ignore, 1 if self._spline else 0, use_sky, sky_level,
            1 if self._disjoint_rings else 0,
            1 if self._disjoint_apertures else 0,
            self._disjoint_radius,
            NULL, 0,
            sg_order, sg_coeff_ptr, sg_vmin,
            self._correlation_length, sigma,
            1 if calc_optimal else 0,
            NULL, 0,
            &out_flux_arr[0,0], &out_fluxerr_arr[0,0],
            &out_bgarea_arr[0,0], &out_bgflux_arr[0,0],
            &out_bgmedian_arr[0,0], &out_bgsigma_arr[0,0],
            &out_cntr_x_arr[0,0], &out_cntr_y_arr[0,0],
            &out_cntr_width_arr[0,0], &out_cntr_w_d_arr[0,0],
            &out_cntr_w_k_arr[0,0], &out_cntr_x_err_arr[0,0],
            &out_cntr_y_err_arr[0,0], &out_cntr_w_err_arr[0,0],
            <int*>&out_flag_arr[0,0], <int*>&out_rtot_arr[0,0],
            <int*>&out_rbad_arr[0,0], <int*>&out_rign_arr[0,0],
            <int*>&out_atot_arr[0,0], <int*>&out_abad_arr[0,0],
            &out_opt_r0_arr[0], &out_opt_ra_arr[0], &out_opt_da_arr[0],
            &out_raw_arr[0])

        # Compute magnitudes
        cdef double [:,::1] flux_arr = out_flux_arr
        cdef double [:,::1] fluxerr_arr = out_fluxerr_arr
        m0, f0 = self._mag_flux
        mag_arr = np.zeros((nstar, nap), dtype=np.float64)
        magerr_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] mag_v = mag_arr
        cdef double [:,::1] magerr_v = magerr_arr
        cdef int mi, mj
        for mi in range(nstar):
            mag_m0 = srefmag_arr[mi] if srefmag_arr is not None and srefmag_arr[mi] != 0 else m0
            mag_f0 = f0
            for mj in range(nap):
                if flux_arr[mi, mj] > 0:
                    mag_v[mi, mj] = mag_m0 - 2.5 * np.log10(flux_arr[mi, mj] / mag_f0)
                    magerr_v[mi, mj] = np.abs(-2.5 * fluxerr_arr[mi, mj] / (flux_arr[mi, mj] * np.log(10)))

        result_dict = {
            'flux': np.asarray(out_flux_arr),
            'fluxerr': np.asarray(out_fluxerr_arr),
            'magnitude': mag_arr,
            'magerr': magerr_arr,
            'bgarea': np.asarray(out_bgarea_arr),
            'bgflux': np.asarray(out_bgflux_arr),
            'bgmedian': np.asarray(out_bgmedian_arr),
            'bgsigma': np.asarray(out_bgsigma_arr),
            'cntr_x': np.asarray(out_cntr_x_arr),
            'cntr_y': np.asarray(out_cntr_y_arr),
            'cntr_width': np.asarray(out_cntr_width_arr),
            'cntr_w_d': np.asarray(out_cntr_w_d_arr),
            'cntr_w_k': np.asarray(out_cntr_w_k_arr),
            'cntr_x_err': np.asarray(out_cntr_x_err_arr),
            'cntr_y_err': np.asarray(out_cntr_y_err_arr),
            'cntr_w_err': np.asarray(out_cntr_w_err_arr),
            'flag': np.asarray(out_flag_arr),
            'rtot': np.asarray(out_rtot_arr),
            'rbad': np.asarray(out_rbad_arr),
            'rign': np.asarray(out_rign_arr),
            'atot': np.asarray(out_atot_arr),
            'abad': np.asarray(out_abad_arr),
        }
        if calc_optimal:
            result_dict['optimal_r0'] = np.asarray(out_opt_r0_arr)
            result_dict['optimal_ra'] = np.asarray(out_opt_ra_arr)
            result_dict['optimal_da'] = np.asarray(out_opt_da_arr)
        if ids_list is not None:
            result_dict['id'] = ids_list
        if nap > 0:
            result_dict['output_raw_photometry_3d_data'] = np.asarray(out_raw_arr).reshape(nstar, nap, 12)

        import types
        output_ns = types.SimpleNamespace()
        output_ns.dict = result_dict

        # Build .output Table (CLI --output format IXY,BbFfMm)
        from astropy.table import Table as Table
        _nrows = nstar * nap
        _flat_ids = []
        _flat_x = []
        _flat_y = []
        for _mi in range(nstar):
            _sid = ids_list[_mi] if ids_list is not None else (_mi + 1)
            for _mj in range(nap):
                _flat_ids.append(_sid)
                _flat_x.append(sx_arr[_mi])
                _flat_y.append(sy_arr[_mi])
        output_ns.output = Table({
            'id': _flat_ids,
            'x': np.array(_flat_x, dtype=np.float64),
            'y': np.array(_flat_y, dtype=np.float64),
            'flux': result_dict['flux'].ravel(),
            'flux_err': result_dict['fluxerr'].ravel(),
            'bg': result_dict['bgmedian'].ravel(),
            'bg_err': result_dict['bgsigma'].ravel(),
            'mag': result_dict['magnitude'].ravel(),
            'mag_err': result_dict['magerr'].ravel(),
        })

        # Build .output_raw_photometry Table (CLI --output-raw-photometry)
        _raw_data_3d = result_dict['output_raw_photometry_3d_data']
        _raw_flat = _raw_data_3d.reshape(-1, 12)
        output_ns.output_raw_photometry = Table({
            'id': _flat_ids,
            'x': np.array(_flat_x, dtype=np.float64),
            'y': np.array(_flat_y, dtype=np.float64),
            'nap': [1] * _nrows,
            'ref_mag': [np.nan] * _nrows,
            'ref_col': [np.nan] * _nrows,
            'r0': _raw_flat[:, 7],
            'ra': _raw_flat[:, 8],
            'da': _raw_flat[:, 9],
            'flux': _raw_flat[:, 10],
            'flux_err': _raw_flat[:, 11],
            'flag': _raw_flat[:, 4].astype(np.int32),
        })

        return output_ns

    def photometry_from_raw(self, raw_data, img_data, mask=None,
                            ref_mag=None, ref_col=None, ref_err=None,
                            kernel_spec=None, normalize_kernel=False,
                            kernel_dict=None):
        """Execute subtracted photometry using pre-computed raw photometry data.

        Corresponds to CLI: --input-raw-photometry <file> -s <subtracted_image>

        Parameters
        ----------
        raw_data : ndarray shape (nstar, nap, 12)
            Raw photometry from a previous photometry() call (r['raw']).
        img_data : float64 ndarray (2D)
            The subtracted image.
        mask : uint8 ndarray (optional)
        ref_mag : ndarray shape (nstar,) (optional)
            Reference magnitudes for calibration.
        ref_col : ndarray shape (nstar,) (optional)
            Photometric color index.
        ref_err : ndarray shape (nstar,) (optional)
            Reference magnitude uncertainties.
        kernel_spec : str (optional)
        normalize_kernel : bool

        Returns
        -------
        dict with keys: flux, fluxerr, magnitude, magerr, bgarea, bgflux,
            bgmedian, bgsigma, cntr_x, cntr_y, cntr_width, cntr_w_d, cntr_w_k,
            cntr_x_err, cntr_y_err, cntr_w_err, flag, rtot, rbad, rign, atot, abad.
        """
        import numpy as np
        cdef int i, j

        # Accept Table or numpy 3D array as input_raw_photometry
        if hasattr(raw_data, 'columns'):
            _tbl = raw_data
            _nstar = len(_tbl)
            _nap = 1  # TODO: extract from Table if multi-aperture support added
            raw_arr = np.zeros((_nstar, _nap, 12), dtype=np.float64)
            for _i in range(_nstar):
                raw_arr[_i, 0, 0] = _i  # index, not used by read_raw_photometry_cy
                raw_arr[_i, 0, 1] = _tbl['x'][_i]
                raw_arr[_i, 0, 2] = _tbl['y'][_i]
                raw_arr[_i, 0, 3] = _tbl['nap'][_i]
                raw_arr[_i, 0, 4] = float(_tbl['flag'][_i])
                raw_arr[_i, 0, 5] = 0.0
                raw_arr[_i, 0, 6] = 0.0
                raw_arr[_i, 0, 7] = _tbl['r0'][_i]
                raw_arr[_i, 0, 8] = _tbl['ra'][_i]
                raw_arr[_i, 0, 9] = _tbl['da'][_i]
                raw_arr[_i, 0, 10] = _tbl['flux'][_i]
                raw_arr[_i, 0, 11] = _tbl['flux_err'][_i]
        else:
            raw_arr = np.ascontiguousarray(raw_data, dtype=np.float64)
        cdef double [:, :, ::1] raw_view = raw_arr
        cdef int nstar = raw_view.shape[0]
        cdef int nap = raw_view.shape[1]

        sy, sx = img_data.shape
        in_img = np.ascontiguousarray(img_data, dtype=np.float64)
        if mask is not None:
            mk_arr = np.ascontiguousarray(mask, dtype=np.uint8)
        else:
            mk_arr = np.zeros((sy, sx), dtype=np.uint8)

        cdef double [:,::1] img_view = in_img
        cdef unsigned char [:,::1] mk_view = mk_arr

        cdef int bg_type = 2, bg_scatter = 0, bg_rejniter = 0
        cdef double bg_rejlower = 3.0, bg_rejupper = 3.0
        cdef int use_sky = 0
        cdef double sky_level = 0.0
        if self._sky_fit is not None:
            for tok in self._sky_fit.replace(',', ' ').split():
                if tok == 'mean': bg_type = 1
                elif tok == 'median': bg_type = 2
                elif tok == 'mode': bg_type = 3
                elif tok == 'mad': bg_scatter = 1
                elif tok == 'stddev': bg_scatter = 0
                elif tok.startswith('iterations='): bg_rejniter = int(tok.split('=')[1])
                elif tok.startswith('lower='): bg_rejlower = float(tok.split('=')[1])
                elif tok.startswith('upper='): bg_rejupper = float(tok.split('=')[1])
                elif tok.startswith('sigma='):
                    s = float(tok.split('=')[1])
                    bg_rejlower = bg_rejupper = s
                elif tok.startswith('force='):
                    sky_level = float(tok.split('=')[1])
                    use_sky = 1

        cdef int mask_ignore = 0
        if isinstance(self._mask_ignore, str) and self._mask_ignore:
            mask_map = {'saturated': 0x02, 'hot': 0x04, 'outer': 0x01}
            for tok in self._mask_ignore.replace(',', ' ').split():
                mask_ignore |= mask_map.get(tok.lower(), 0)
        elif isinstance(self._mask_ignore, int):
            mask_ignore = self._mask_ignore

        cdef int sg_order = 0
        cdef double sg_vmin = self._gain_vmin
        cdef double [:] sg_coeff_view
        cdef double *sg_coeff_ptr = NULL
        parts = [float(x) for x in self._gain.replace(',', ' ').split()]
        sg_order = int((1 + 8 * len(parts)) ** 0.5 - 3) // 2
        ncoeff = (sg_order + 1) * (sg_order + 2) // 2
        sg_coeff_arr = np.array(parts[:ncoeff], dtype=np.float64)
        sg_coeff_view = sg_coeff_arr
        sg_coeff_ptr = &sg_coeff_view[0]

        cdef bytes kernel_spec_b
        cdef char *kernel_spec_ptr = NULL
        if kernel_spec:
            kernel_spec_b = kernel_spec.encode('utf-8')
            kernel_spec_ptr = kernel_spec_b
        cdef int normalize_kernel_int = 1 if normalize_kernel else 0

        cdef double [:] ref_mag_v
        cdef double *ref_mag_ptr = NULL
        if ref_mag is not None:
            ref_mag_v = np.ascontiguousarray(ref_mag, dtype=np.float64)
            ref_mag_ptr = &ref_mag_v[0]

        cdef double [:] ref_col_v
        cdef double *ref_col_ptr = NULL
        if ref_col is not None:
            ref_col_v = np.ascontiguousarray(ref_col, dtype=np.float64)
            ref_col_ptr = &ref_col_v[0]

        cdef double [:] ref_err_v
        cdef double *ref_err_ptr = NULL
        if ref_err is not None:
            ref_err_v = np.ascontiguousarray(ref_err, dtype=np.float64)
            ref_err_ptr = &ref_err_v[0]

        # Output arrays
        cdef double [:,::1] out_flux_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_fluxerr_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgarea_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgflux_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgmedian_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgsigma_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_x_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_y_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_width_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_d_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_k_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_x_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_y_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef int [:,::1] out_flag_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rtot_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rbad_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rign_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_atot_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_abad_arr = np.zeros((nstar, nap), dtype=np.int32)

        cdef double sigma = self._fwhm / 2.354820045

        cdef int kd_nkernel = 0
        cdef double kd_ox = 0.0, kd_oy = 0.0, kd_scale = 1.0
        cdef int kd_ktype = 0
        cdef int [:] kd_ktv, kd_kov, kd_knv, kd_khv
        cdef double [:] kd_ksv, kd_kxv_d, kd_kyv_d, kd_kcfv
        cdef int [:] kd_kxv, kd_kyv
        cdef double **kd_coeff_ptrs = NULL
        cdef int kd_ki
        if kernel_dict is not None:
            kd_nkernel, kd_ox, kd_oy, kd_scale, kd_ktype, \
            _kt, _ko, _kn, _kh, _ks, _kx, _ky, _kcf, kd_offsets = \
                fiphot_build_kernel_arrays(kernel_dict)
            kd_ktv = _kt; kd_kov = _ko; kd_knv = _kn
            kd_khv = _kh; kd_ksv = _ks
            kd_kxv = _kx; kd_kyv = _ky
            kd_kcfv = _kcf
            kd_coeff_ptrs = <double**>malloc(kd_nkernel * sizeof(double*))
            for kd_ki in range(kd_nkernel):
                kd_coeff_ptrs[kd_ki] = &kd_kcfv[kd_offsets[kd_ki]]

        fiphot_photometry_from_raw_cy(
            &raw_view[0,0,0], nstar, nap,
            ref_mag_ptr, ref_col_ptr, ref_err_ptr,
            &img_view[0,0], <char*>&mk_view[0,0], sx, sy,
            bg_type, bg_scatter, bg_rejniter, bg_rejlower, bg_rejupper,
            mask_ignore, 1 if self._spline else 0, use_sky, sky_level,
            1 if self._disjoint_rings else 0,
            1 if self._disjoint_apertures else 0,
            self._disjoint_radius,
            NULL, 0,
            sg_order, sg_coeff_ptr, sg_vmin,
            self._correlation_length, sigma,
            NULL, 0,
            kd_nkernel,
            kd_ox, kd_oy, kd_scale, kd_ktype,
            <int*>&kd_ktv[0] if kd_nkernel > 0 else NULL,
            <int*>&kd_kov[0] if kd_nkernel > 0 else NULL,
            <int*>&kd_knv[0] if kd_nkernel > 0 else NULL,
            <int*>&kd_khv[0] if kd_nkernel > 0 else NULL,
            <double*>&kd_ksv[0] if kd_nkernel > 0 else NULL,
            <int*>&kd_kxv[0] if kd_nkernel > 0 else NULL,
            <int*>&kd_kyv[0] if kd_nkernel > 0 else NULL,
            kd_coeff_ptrs,
            kernel_spec_ptr, normalize_kernel_int,
            &out_flux_arr[0,0], &out_fluxerr_arr[0,0],
            &out_bgarea_arr[0,0], &out_bgflux_arr[0,0],
            &out_bgmedian_arr[0,0], &out_bgsigma_arr[0,0],
            &out_cntr_x_arr[0,0], &out_cntr_y_arr[0,0],
            &out_cntr_width_arr[0,0], &out_cntr_w_d_arr[0,0],
            &out_cntr_w_k_arr[0,0], &out_cntr_x_err_arr[0,0],
            &out_cntr_y_err_arr[0,0], &out_cntr_w_err_arr[0,0],
            <int*>&out_flag_arr[0,0], <int*>&out_rtot_arr[0,0],
            <int*>&out_rbad_arr[0,0], <int*>&out_rign_arr[0,0],
            <int*>&out_atot_arr[0,0], <int*>&out_abad_arr[0,0])

        if kd_coeff_ptrs != NULL:
            free(kd_coeff_ptrs)

        # Compute magnitudes
        cdef double [:,::1] flux_arr = out_flux_arr
        cdef double [:,::1] fluxerr_arr = out_fluxerr_arr
        m0, f0 = self._mag_flux
        mag_arr = np.zeros((nstar, nap), dtype=np.float64)
        magerr_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] mag_v = mag_arr
        cdef double [:,::1] magerr_v = magerr_arr
        cdef int mi, mj
        for mi in range(nstar):
            for mj in range(nap):
                if flux_arr[mi, mj] > 0:
                    mag_v[mi, mj] = m0 - 2.5 * np.log10(flux_arr[mi, mj] / f0)
                    magerr_v[mi, mj] = np.abs(-2.5 * fluxerr_arr[mi, mj] / (flux_arr[mi, mj] * np.log(10)))

        result_dict = {
            'flux': np.asarray(out_flux_arr),
            'fluxerr': np.asarray(out_fluxerr_arr),
            'magnitude': mag_arr,
            'magerr': magerr_arr,
            'bgarea': np.asarray(out_bgarea_arr),
            'bgflux': np.asarray(out_bgflux_arr),
            'bgmedian': np.asarray(out_bgmedian_arr),
            'bgsigma': np.asarray(out_bgsigma_arr),
            'cntr_x': np.asarray(out_cntr_x_arr),
            'cntr_y': np.asarray(out_cntr_y_arr),
            'cntr_width': np.asarray(out_cntr_width_arr),
            'cntr_w_d': np.asarray(out_cntr_w_d_arr),
            'cntr_w_k': np.asarray(out_cntr_w_k_arr),
            'cntr_x_err': np.asarray(out_cntr_x_err_arr),
            'cntr_y_err': np.asarray(out_cntr_y_err_arr),
            'cntr_w_err': np.asarray(out_cntr_w_err_arr),
            'flag': np.asarray(out_flag_arr),
            'rtot': np.asarray(out_rtot_arr),
            'rbad': np.asarray(out_rbad_arr),
            'rign': np.asarray(out_rign_arr),
            'atot': np.asarray(out_atot_arr),
            'abad': np.asarray(out_abad_arr),
        }

        import types
        output_ns = types.SimpleNamespace()
        output_ns.dict = result_dict

        # Build .output Table (CLI --output format IXY,BbFfMm)
        from astropy.table import Table as Table
        _nrows = nstar * nap
        _flat_ids = []
        _flat_x = []
        _flat_y = []
        for _mi in range(nstar):
            for _mj in range(nap):
                _flat_ids.append(int(raw_arr[_mi, _mj, 0]))
                _flat_x.append(raw_arr[_mi, _mj, 1])
                _flat_y.append(raw_arr[_mi, _mj, 2])

        output_ns.output = Table({
            'id': _flat_ids,
            'x': np.array(_flat_x, dtype=np.float64),
            'y': np.array(_flat_y, dtype=np.float64),
            'flux': result_dict['flux'].ravel(),
            'flux_err': result_dict['fluxerr'].ravel(),
            'bg': result_dict['bgmedian'].ravel(),
            'bg_err': result_dict['bgsigma'].ravel(),
            'mag': result_dict['magnitude'].ravel(),
            'mag_err': result_dict['magerr'].ravel(),
        })
        return output_ns

    def subtracted_photometry(self, img_data, mask=None, stars=None,
                              col_xy=(0, 1), col_ap=None,
                              col_id=-1, col_mag=-1, col_col=-1, col_err=-1,
                              kernel_spec=None, normalize_kernel=False):
        """Perform aperture photometry on a subtracted image.

        Similar to photometry, but operates on a pre-subtracted image where
        neighbour stars have already been removed via kernel convolution.

        Additional parameters:
            kernel_spec: str — kernel specification string used for subtraction
            normalize_kernel: bool — normalize the identity kernel to flux=1
        """
        import numpy as np
        cdef int i, j

        sy, sx = img_data.shape
        in_img = np.ascontiguousarray(img_data, dtype=np.float64)
        if mask is not None:
            mk_arr = np.ascontiguousarray(mask, dtype=np.uint8)
        else:
            mk_arr = np.zeros((sy, sx), dtype=np.uint8)

        cdef double [:,::1] img_view = in_img
        cdef unsigned char [:,::1] mk_view = mk_arr

        cdef int nstar, nap, has_per_star_aps
        cdef double [:] sx_arr, sy_arr
        cdef double [:] sr0_arr, sra_arr, sda_arr
        cdef double [:] srefmag_arr, srefcol_arr, sreferr_arr
        cdef list ids_list

        if stars is None:
            raise ValueError("stars is required")

        cdef double [:,::1] star_mat
        cdef bint is_ndarray = isinstance(stars, np.ndarray)
        if is_ndarray:
            star_mat = np.ascontiguousarray(stars, dtype=np.float64)
            nstar = star_mat.shape[0]
            ids_list = None
            cx, cy = col_xy[0], col_xy[1]

            sx_arr = np.asarray(star_mat[:, cx]) * self._zoom
            sy_arr = np.asarray(star_mat[:, cy]) * self._zoom

            if col_ap is not None and len(col_ap) > 0:
                has_per_star_aps = True
                nap = len(col_ap)
                sr0_arr = np.zeros(nstar * nap, dtype=np.float64)
                sra_arr = np.zeros(nstar * nap, dtype=np.float64)
                sda_arr = np.zeros(nstar * nap, dtype=np.float64)
                for j, acol in enumerate(col_ap):
                    for i in range(nstar):
                        apspec = str(stars[i, acol])
                        parts = apspec.replace(',', ':').split(':')
                        idx = i * nap + j
                        if len(parts) >= 3:
                            sr0_arr[idx] = float(parts[0]) * self._zoom
                            sra_arr[idx] = float(parts[1]) * self._zoom
                            sda_arr[idx] = float(parts[2]) * self._zoom
                        elif len(parts) == 1:
                            sr0_arr[idx] = float(parts[0]) * self._zoom
                            sra_arr[idx] = 2.0 * self._zoom
                            sda_arr[idx] = 1.0 * self._zoom
            else:
                has_per_star_aps = False
                nap = 0
                sr0_arr = np.zeros(0, dtype=np.float64)
                sra_arr = np.zeros(0, dtype=np.float64)
                sda_arr = np.zeros(0, dtype=np.float64)

            srefmag_arr = np.zeros(nstar, dtype=np.float64)
            srefcol_arr = np.zeros(nstar, dtype=np.float64)
            sreferr_arr = np.zeros(nstar, dtype=np.float64)
            col_id_used = col_id
            if col_mag >= 0 and col_mag < star_mat.shape[1]:
                srefmag_arr = np.asarray(star_mat[:, col_mag])
            if col_col >= 0 and col_col < star_mat.shape[1]:
                srefcol_arr = np.asarray(star_mat[:, col_col])
            if col_err >= 0 and col_err < star_mat.shape[1]:
                sreferr_arr = np.asarray(star_mat[:, col_err])
            if col_id >= 0 and col_id < star_mat.shape[1]:
                ids_list = [stars[i, col_id] for i in range(nstar)]

        else:
            nstar = len(stars)
            first = stars[0]
            if isinstance(first, (list, tuple)) and len(first) >= 5:
                has_per_star_aps = True
                nap = (len(first) - 2) // 3
            else:
                has_per_star_aps = False
                nap = 0

            sx_arr = np.zeros(nstar, dtype=np.float64)
            sy_arr = np.zeros(nstar, dtype=np.float64)
            srefmag_arr = np.zeros(nstar, dtype=np.float64)
            srefcol_arr = np.zeros(nstar, dtype=np.float64)
            sreferr_arr = np.zeros(nstar, dtype=np.float64)

            if has_per_star_aps:
                sr0_arr = np.zeros(nstar * nap, dtype=np.float64)
                sra_arr = np.zeros(nstar * nap, dtype=np.float64)
                sda_arr = np.zeros(nstar * nap, dtype=np.float64)
            else:
                sr0_arr = np.zeros(0, dtype=np.float64)
                sra_arr = np.zeros(0, dtype=np.float64)
                sda_arr = np.zeros(0, dtype=np.float64)

            for i in range(nstar):
                s = stars[i]
                sx_arr[i] = s[0] * self._zoom
                sy_arr[i] = s[1] * self._zoom
                if has_per_star_aps:
                    for j in range(nap):
                        idx = i * nap + j
                        sr0_arr[idx] = s[2 + 3*j] * self._zoom
                        sra_arr[idx] = s[3 + 3*j] * self._zoom
                        sda_arr[idx] = s[4 + 3*j] * self._zoom
            ids_list = None
            col_id_used = -1

        cdef bytes ap_spec_b = b''
        if not has_per_star_aps and self._apertures is not None:
            if isinstance(self._apertures, (list, tuple)):
                ap_spec = ','.join(str(a) for a in self._apertures)
            else:
                ap_spec = str(self._apertures)
            ap_spec_b = ap_spec.encode('utf-8')
            nap = ap_spec.count(',') + 1
        elif not has_per_star_aps:
            raise ValueError("Either provide per-star apertures or global apertures")

        cdef char *ap_spec_ptr = NULL
        if len(ap_spec_b) > 0:
            ap_spec_ptr = ap_spec_b

        cdef int bg_type = 2
        cdef int bg_scatter = 0
        cdef int bg_rejniter = 0
        cdef double bg_rejlower = 3.0
        cdef double bg_rejupper = 3.0
        cdef int use_sky = 0
        cdef double sky_level = 0.0
        if self._sky_fit is not None:
            for tok in self._sky_fit.replace(',', ' ').split():
                if tok == 'mean': bg_type = 1
                elif tok == 'median': bg_type = 2
                elif tok == 'mode': bg_type = 3
                elif tok == 'mad': bg_scatter = 1
                elif tok == 'stddev': bg_scatter = 0
                elif tok.startswith('iterations='): bg_rejniter = int(tok.split('=')[1])
                elif tok.startswith('lower='): bg_rejlower = float(tok.split('=')[1])
                elif tok.startswith('upper='): bg_rejupper = float(tok.split('=')[1])
                elif tok.startswith('sigma='):
                    s = float(tok.split('=')[1])
                    bg_rejlower = bg_rejupper = s
                elif tok.startswith('force='):
                    sky_level = float(tok.split('=')[1])
                    use_sky = 1

        cdef int mask_ignore = 0
        if isinstance(self._mask_ignore, str) and self._mask_ignore:
            mask_map = {'saturated': 0x02, 'hot': 0x04, 'outer': 0x01}
            for tok in self._mask_ignore.replace(',', ' ').split():
                mask_ignore |= mask_map.get(tok.lower(), 0)
        elif isinstance(self._mask_ignore, int):
            mask_ignore = self._mask_ignore

        cdef int sg_order = 0
        cdef double sg_vmin = self._gain_vmin
        cdef double [:] sg_coeff_view
        cdef double *sg_coeff_ptr = NULL
        parts = [float(x) for x in self._gain.replace(',', ' ').split()]
        sg_order = int((1 + 8 * len(parts)) ** 0.5 - 3) // 2
        ncoeff = (sg_order + 1) * (sg_order + 2) // 2
        sg_coeff_arr = np.array(parts[:ncoeff], dtype=np.float64)
        sg_coeff_view = sg_coeff_arr
        sg_coeff_ptr = &sg_coeff_view[0]

        cdef bytes kernel_spec_b
        cdef char *kernel_spec_ptr = NULL
        if kernel_spec:
            kernel_spec_b = kernel_spec.encode('utf-8')
            kernel_spec_ptr = kernel_spec_b
        cdef int normalize_kernel_int = 1 if normalize_kernel else 0

        # Output arrays
        cdef double [:,::1] out_flux_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_fluxerr_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgarea_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgflux_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgmedian_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_bgsigma_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_x_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_y_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_width_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_d_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_k_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_x_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_y_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] out_cntr_w_err_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef int [:,::1] out_flag_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rtot_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rbad_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_rign_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_atot_arr = np.zeros((nstar, nap), dtype=np.int32)
        cdef int [:,::1] out_abad_arr = np.zeros((nstar, nap), dtype=np.int32)

        cdef double *sr0_ptr = NULL
        cdef double *sra_ptr = NULL
        cdef double *sda_ptr = NULL
        if has_per_star_aps and nap > 0:
            sr0_ptr = &sr0_arr[0]
            sra_ptr = &sra_arr[0]
            sda_ptr = &sda_arr[0]

        fiphot_subtracted_photometry_cy(
            &img_view[0,0], <char*>&mk_view[0,0], sx, sy,
            &sx_arr[0], &sy_arr[0], nstar,
            sr0_ptr, sra_ptr, sda_ptr,
            nap if has_per_star_aps else 0,
            ap_spec_ptr, self._zoom,
            bg_type, bg_scatter, bg_rejniter, bg_rejlower, bg_rejupper,
            mask_ignore, 1 if self._spline else 0, use_sky, sky_level,
            1 if self._disjoint_rings else 0,
            1 if self._disjoint_apertures else 0,
            self._disjoint_radius,
            NULL, 0,
            sg_order, sg_coeff_ptr, sg_vmin,
            self._correlation_length,
            kernel_spec_ptr, normalize_kernel_int,
            NULL, 0,
            &out_flux_arr[0,0], &out_fluxerr_arr[0,0],
            &out_bgarea_arr[0,0], &out_bgflux_arr[0,0],
            &out_bgmedian_arr[0,0], &out_bgsigma_arr[0,0],
            &out_cntr_x_arr[0,0], &out_cntr_y_arr[0,0],
            &out_cntr_width_arr[0,0], &out_cntr_w_d_arr[0,0],
            &out_cntr_w_k_arr[0,0], &out_cntr_x_err_arr[0,0],
            &out_cntr_y_err_arr[0,0], &out_cntr_w_err_arr[0,0],
            <int*>&out_flag_arr[0,0], <int*>&out_rtot_arr[0,0],
            <int*>&out_rbad_arr[0,0], <int*>&out_rign_arr[0,0],
            <int*>&out_atot_arr[0,0], <int*>&out_abad_arr[0,0])

        # Compute magnitudes
        cdef double [:,::1] flux_arr = out_flux_arr
        cdef double [:,::1] fluxerr_arr = out_fluxerr_arr
        m0, f0 = self._mag_flux
        mag_arr = np.zeros((nstar, nap), dtype=np.float64)
        magerr_arr = np.zeros((nstar, nap), dtype=np.float64)
        cdef double [:,::1] mag_v = mag_arr
        cdef double [:,::1] magerr_v = magerr_arr
        cdef int mi, mj
        for mi in range(nstar):
            mag_m0 = srefmag_arr[mi] if srefmag_arr is not None and srefmag_arr[mi] != 0 else m0
            mag_f0 = f0
            for mj in range(nap):
                if flux_arr[mi, mj] > 0:
                    mag_v[mi, mj] = mag_m0 - 2.5 * np.log10(flux_arr[mi, mj] / mag_f0)
                    magerr_v[mi, mj] = np.abs(-2.5 * fluxerr_arr[mi, mj] / (flux_arr[mi, mj] * np.log(10)))

        result = {
            'flux': np.asarray(out_flux_arr),
            'fluxerr': np.asarray(out_fluxerr_arr),
            'magnitude': mag_arr,
            'magerr': magerr_arr,
            'bgarea': np.asarray(out_bgarea_arr),
            'bgflux': np.asarray(out_bgflux_arr),
            'bgmedian': np.asarray(out_bgmedian_arr),
            'bgsigma': np.asarray(out_bgsigma_arr),
            'cntr_x': np.asarray(out_cntr_x_arr),
            'cntr_y': np.asarray(out_cntr_y_arr),
            'cntr_width': np.asarray(out_cntr_width_arr),
            'cntr_w_d': np.asarray(out_cntr_w_d_arr),
            'cntr_w_k': np.asarray(out_cntr_w_k_arr),
            'cntr_x_err': np.asarray(out_cntr_x_err_arr),
            'cntr_y_err': np.asarray(out_cntr_y_err_arr),
            'cntr_w_err': np.asarray(out_cntr_w_err_arr),
            'flag': np.asarray(out_flag_arr),
            'rtot': np.asarray(out_rtot_arr),
            'rbad': np.asarray(out_rbad_arr),
            'rign': np.asarray(out_rign_arr),
            'atot': np.asarray(out_atot_arr),
            'abad': np.asarray(out_abad_arr),
        }
        if ids_list is not None:
            result['id'] = ids_list

        import types
        output_ns = types.SimpleNamespace()
        output_ns.dict = result

        from astropy.table import Table as Table
        _nrows = nstar * nap
        _flat_ids = []
        _flat_x = []
        _flat_y = []
        for _mi in range(nstar):
            _sid = ids_list[_mi] if ids_list is not None else (_mi + 1)
            for _mj in range(nap):
                _flat_ids.append(_sid)
                _flat_x.append(sx_arr[_mi])
                _flat_y.append(sy_arr[_mi])
        output_ns.output = Table({
            'id': _flat_ids,
            'x': np.array(_flat_x, dtype=np.float64),
            'y': np.array(_flat_y, dtype=np.float64),
            'flux': result['flux'].ravel(),
            'flux_err': result['fluxerr'].ravel(),
            'bg': result['bgmedian'].ravel(),
            'bg_err': result['bgsigma'].ravel(),
            'mag': result['magnitude'].ravel(),
            'mag_err': result['magerr'].ravel(),
        })
        return output_ns

    def magfit(self, flux, flag=None, ref_flux=None, ref_col=None, ref_mag=None,
                      ref_err=None, orders=None, niter=3, sigma=3.0):
        """Fit instrumental magnitudes to reference magnitudes (magfit).

        Performs a polynomial fit between instrumental and reference magnitudes
        to derive calibrated magnitudes, with iterative sigma-clipping.

        Parameters
        ----------
        flux : ndarray shape (nstar, nap)
            Measured flux values per star per aperture.
        flag : ndarray shape (nstar, nap), optional
            Photometry flags (zero = valid). Stars with any non-zero flag
            in any aperture are excluded from the fit.
        ref_flux : ndarray shape (nstar, nap), optional
            Reference frame photometry fluxes. Required for proper calibration.
            Should come from photometry() on the reference image.
        ref_col : ndarray shape (nstar,), optional
            Photometric color index for each star.
        ref_mag : ndarray shape (nstar,), optional
            Reference magnitudes for each star.
        ref_err : ndarray shape (nstar,), optional
            Reference magnitude uncertainties.
        orders : list of int, optional
            Polynomial orders per fitting term. E.g. [2, 1, 0] for
            magnitude^2, color^1, constant. Default: [0] (constant offset).
        niter : int
            Number of sigma-clipping iterations. Default 3.
        sigma : float
            Rejection threshold in sigma. Default 3.0.

        Returns
        -------
        dict with keys:
            mag : ndarray shape (nstar, nap) — calibrated magnitudes
            ninit : int — initial number of stars used
            nrejs : int — number of stars rejected
            nstar : int — final number of stars used
            naperture : int — number of apertures used
        """
        import numpy as np
        cdef int i, j

        cdef double [:,::1] flux_v = np.ascontiguousarray(flux, dtype=np.float64)
        cdef int nstar = flux_v.shape[0]
        cdef int nap = flux_v.shape[1]

        cdef int [:,::1] flag_v
        if flag is not None:
            flag_v = np.ascontiguousarray(flag, dtype=np.int32)
        else:
            flag_v = np.zeros((nstar, nap), dtype=np.int32)

        cdef double [:] ref_col_v
        if ref_col is not None:
            ref_col_v = np.ascontiguousarray(ref_col, dtype=np.float64)
        else:
            ref_col_v = np.zeros(nstar, dtype=np.float64)

        cdef double [:] ref_mag_v
        if ref_mag is not None:
            ref_mag_v = np.ascontiguousarray(ref_mag, dtype=np.float64)
        else:
            ref_mag_v = np.zeros(nstar, dtype=np.float64)

        cdef double [:] ref_err_v
        if ref_err is not None:
            ref_err_v = np.ascontiguousarray(ref_err, dtype=np.float64)
        else:
            ref_err_v = np.zeros(nstar, dtype=np.float64)

        cdef double *ref_flux_ptr = NULL
        cdef double [:] ref_flux_v
        if ref_flux is not None:
            ref_flux_v = np.ascontiguousarray(ref_flux, dtype=np.float64).ravel()
            ref_flux_ptr = &ref_flux_v[0]

        if orders is None:
            orders = [0]
        cdef int [:] orders_v = np.ascontiguousarray(orders, dtype=np.int32)
        cdef int norder = len(orders)

        cdef double [:] out_mag_v = np.zeros(nstar * nap, dtype=np.float64)
        cdef int out_ninit = 0, out_nrejs = 0, out_nstar = 0, out_naperture = 0

        fiphot_magnitude_fit_cy(
            &flux_v[0,0], &flag_v[0,0], nstar, nap,
            ref_flux_ptr, &ref_col_v[0], &ref_mag_v[0], &ref_err_v[0],
            &orders_v[0], norder, niter, sigma,
            nstar, nap,
            &out_mag_v[0],
            &out_ninit, &out_nrejs, &out_nstar, &out_naperture)

        result = {
            'mag': np.asarray(out_mag_v).reshape(nstar, nap),
            'ninit': out_ninit,
            'nrejs': out_nrejs,
            'nstar': out_nstar,
            'naperture': out_naperture,
        }

        import types
        output_ns = types.SimpleNamespace()
        output_ns.dict = result

        from astropy.table import Table as Table
        output_ns.output = Table({'mag': np.asarray(out_mag_v).reshape(nstar, nap).ravel()})
        return output_ns

