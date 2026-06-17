# grmatch/grmatch.pxi — star catalog matching
cdef class MatchResult:
    """Result container for do_pointmatch output."""
    cdef public int nhit
    cdef public list hits_ref
    cdef public list hits_inp
    cdef public list excluded_ref
    cdef public list excluded_inp
    cdef public list vfits_dx
    cdef public list vfits_dy
    cdef public dict stats
    cdef public int order
    cdef public int nvar

    property transform:
        def __get__(self):
            return {
                'type': 'polynomial',
                'order': self.order,
                'offset': [0.0, 0.0],
                'scale': 1.0,
                'dxfit': self.vfits_dx,
                'dyfit': self.vfits_dy,
            }

    def to_dict(self):
        return {
            'nhit': self.nhit,
            'order': self.order,
            'nvar': self.nvar,
            'hits_ref': self.hits_ref,
            'hits_inp': self.hits_inp,
            'excluded_ref': self.excluded_ref,
            'excluded_inp': self.excluded_inp,
            'vfits_dx': self.vfits_dx,
            'vfits_dy': self.vfits_dy,
            'stats': self.stats,
        }



cdef class Grmatch:
    """
    Wrapper around the C do_pointmatch function.

    Usage:
        m = Matcher(order=2, maxdist=2.0, unitarity=0.01)
        m.set_reference(ref_x, ref_y, ref_ord)
        m.set_input(inp_x, inp_y, inp_ord)
        result = m.run()
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

    def set_reference(self, x, y, ord_=None):
        """Set reference point arrays. ord_=None means skip internal sort."""
        self.nref = len(x)
        assert len(y) == self.nref
        assert ord_ is None or len(ord_) == self.nref

        if self.ref_x: free(self.ref_x)
        if self.ref_y: free(self.ref_y)
        if self.ref_ord: free(self.ref_ord)

        self.ref_x = <double*>malloc(self.nref * sizeof(double))
        self.ref_y = <double*>malloc(self.nref * sizeof(double))
        self.ref_ord = <double*>malloc(self.nref * sizeof(double))

        cdef int i
        for i in range(self.nref):
            self.ref_x[i] = x[i]
            self.ref_y[i] = y[i]
            self.ref_ord[i] = ord_[i] if ord_ is not None else 0.0
        if ord_ is None:
            self._use_ordering = 0

    def set_input(self, x, y, ord_=None):
        """Set input point arrays. ord_=None means skip internal sort."""
        self.ninp = len(x)
        assert len(y) == self.ninp
        assert ord_ is None or len(ord_) == self.ninp

        if self.inp_x: free(self.inp_x)
        if self.inp_y: free(self.inp_y)
        if self.inp_ord: free(self.inp_ord)

        self.inp_x = <double*>malloc(self.ninp * sizeof(double))
        self.inp_y = <double*>malloc(self.ninp * sizeof(double))
        self.inp_ord = <double*>malloc(self.ninp * sizeof(double))

        cdef int i
        for i in range(self.ninp):
            self.inp_x[i] = x[i]
            self.inp_y[i] = y[i]
            self.inp_ord[i] = ord_[i] if ord_ is not None else 0.0
        if ord_ is None:
            self._use_ordering = 0

    def matchxy(self, ref_x, ref_y, inp_x, inp_y,
                ref_ord=None, inp_ord=None):
        """Set reference and input data. ref_ord/inp_ord=None skips
        internal brightness sort; provide (e.g. -mag) to enable."""
        self.set_reference(ref_x, ref_y, ref_ord)
        self.set_input(inp_x, inp_y, inp_ord)

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
            raise RuntimeError("Call matchxy() or set_reference()+set_input() first")
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


