from pyfitsh.gropt_core import _OpticsHandle


class Optics:

    @classmethod
    def from_string(cls, text):
        obj = cls.__new__(cls)
        obj._handle = _OpticsHandle.from_string(text)
        return obj

    @property
    def nglass(self):
        return self._handle.nglass

    @property
    def nlens(self):
        return self._handle.nlens

    @property
    def z_focal(self):
        return self._handle.z_focal

    @z_focal.setter
    def z_focal(self, val):
        self._handle.z_focal = val

    def transfer(self, wavelength=0.6):
        return self._handle.transfer(wavelength)

    def spot(self, wavelength=0.6, aperture_radius=10.0, nrings=5, **kw):
        return self._handle.spot(wavelength, aperture_radius, nrings, **kw)

    def psf(self, wavelength=0.6, aperture_radius=10.0, **kw):
        return self._handle.psf(wavelength, aperture_radius, **kw)

    def raytrace(self, wavelength, x0, y0, z0, nx, ny, nz):
        return self._handle.raytrace(wavelength, x0, y0, z0, nx, ny, nz)

    def to_openscad(self):
        return self._handle.to_openscad()

    def to_eps(self, **kw):
        return self._handle.to_eps(**kw)
