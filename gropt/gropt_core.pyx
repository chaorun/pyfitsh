import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free
from libc.string cimport memset

include "gropt_cdef.pxd"

cdef class _OpticsHandle:
    cdef optics _opt
    cdef bint _valid

    def __cinit__(self):
        optcalc_reset(&self._opt)
        self._valid = False

    def __dealloc__(self):
        if self._opt.opt_glasses != NULL:
            free(self._opt.opt_glasses)
            self._opt.opt_glasses = NULL
        if self._opt.opt_lenses != NULL:
            free(self._opt.opt_lenses)
            self._opt.opt_lenses = NULL

    @staticmethod
    def from_string(str text):
        cdef _OpticsHandle h = _OpticsHandle()
        cdef bytes btext = text.encode('utf-8')
        cdef int ret
        ret = optcalc_read_optics_from_string(btext, &h._opt)
        if ret > 0:
            raise ValueError(f"Parse error at line {ret}")
        if ret < 0:
            raise ValueError("Invalid input")
        h._valid = True
        return h

    @property
    def nglass(self):
        return self._opt.opt_nglass

    @property
    def nlens(self):
        return self._opt.opt_nlens

    @property
    def z_focal(self):
        return self._opt.opt_z_focal

    @z_focal.setter
    def z_focal(self, double val):
        self._opt.opt_z_focal = val

    def transfer(self, double wavelength):
        cdef double focal_plane = 0.0, eff_focus = 0.0
        gropt_transfer(&self._opt, wavelength, &focal_plane, &eff_focus)
        return {'focal_plane': focal_plane, 'effective_focus': eff_focus}

    def spot(self, double wavelength, double aperture_radius, int nrings,
             double angle=0.0, tuple angle_xy=None,
             double zstart=0.0, double pixel_scale=1.0):
        cdef double angle_nx = 0.0, angle_ny = 0.0
        cdef double center_x = 0.0, center_y = 0.0
        cdef int npoints = 0
        cdef int maxpts

        if angle_xy is not None:
            angle_nx = angle_xy[0]
            angle_ny = angle_xy[1]
        elif angle != 0.0:
            angle_nx = 0.0
            angle_ny = np.sin(angle)

        maxpts = 3 * nrings * (nrings + 1)
        cdef np.ndarray[double, ndim=2] out = np.zeros((maxpts, 2), dtype=np.float64)

        gropt_spot_diagram(&self._opt,
            wavelength, aperture_radius, nrings,
            angle_nx, angle_ny, zstart, pixel_scale,
            &out[0, 0], maxpts, &npoints,
            &center_x, &center_y)

        return {
            'spots': out[:npoints],
            'center_x': center_x,
            'center_y': center_y,
        }

    def psf(self, double wavelength, double aperture_radius,
            int half_size=3, double pixel_scale=1.0,
            double angle=0.0, tuple angle_xy=None, double zstart=0.0):
        cdef double angle_nx = 0.0, angle_ny = 0.0
        cdef int psf_size

        if angle_xy is not None:
            angle_nx = angle_xy[0]
            angle_ny = angle_xy[1]
        elif angle != 0.0:
            angle_nx = 0.0
            angle_ny = np.sin(angle)

        psf_size = 2 * half_size + 1
        cdef np.ndarray[double, ndim=2] out = np.zeros((psf_size, psf_size), dtype=np.float64)

        gropt_compute_psf(&self._opt,
            wavelength, aperture_radius,
            angle_nx, angle_ny, zstart, pixel_scale,
            half_size, &out[0, 0])

        return out

    def raytrace(self, double wavelength,
                 double x0, double y0, double z0,
                 double nx, double ny, double nz):
        cdef int npoints = 0
        cdef int maxpts = 256
        cdef np.ndarray[double, ndim=2] out = np.zeros((maxpts, 3), dtype=np.float64)

        gropt_single_raytrace(&self._opt, wavelength,
            x0, y0, z0, nx, ny, nz,
            &out[0, 0], maxpts, &npoints)

        return out[:npoints]

    def to_openscad(self):
        cdef char *buf = NULL
        cdef size_t buflen = 0
        gropt_export_scad(&self._opt, &buf, &buflen)
        if buf == NULL:
            return ""
        result = buf[:buflen].decode('utf-8')
        free(buf)
        return result

    def to_eps(self, double aperture_radius=0.0, int nrings=0,
               double angle=0.0, tuple angle_xy=None,
               double zstart=0.0, double wavelength=0.6,
               double pixel_scale=1.0):
        cdef double angle_nx = 0.0, angle_ny = 0.0
        cdef char *buf = NULL
        cdef size_t buflen = 0

        if angle_xy is not None:
            angle_nx = angle_xy[0]
            angle_ny = angle_xy[1]
        elif angle != 0.0:
            angle_nx = 0.0
            angle_ny = np.sin(angle)

        gropt_export_eps(&self._opt,
            aperture_radius, nrings,
            angle_nx, angle_ny, zstart,
            wavelength, pixel_scale,
            &buf, &buflen)

        if buf == NULL:
            return ""
        result = buf[:buflen].decode('utf-8')
        free(buf)
        return result
