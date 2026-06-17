# grmatch_result.pxi — MatchResult container
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

    property transformation:
        def __get__(self):
            """Export transformation matching .trans file format."""
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
            'transformation': {
                'type': 'polynomial',
                'order': self.order,
                'offset': [0.0, 0.0],
                'scale': 1.0,
                'dxfit': self.vfits_dx,
                'dyfit': self.vfits_dy,
            },
            'stats': self.stats,
        }


