# wcs_helpers.pxi — WCS projection helpers
cdef extern from "wcs_core.h":
    ctypedef struct wcs_fit_input:
        double ra0, de0, roll0
        int type, order
        int n_k, n_p
        double k[4], p[4]

    ctypedef struct wcs_fit_output:
        double crpix1, crpix2
        double cd11, cd12, cd21, cd22
        double res_prj, res_pix
        int order, nvar
        double *pixpoly_dx
        double *pixpoly_dy
        double *prjpoly_dx
        double *prjpoly_dy

    int wcs_fit_cy(
        double *ra, double *dec, double *px, double *py, int npts,
        const wcs_fit_input *cfg,
        wcs_fit_output *out) nogil

    void wcs_fit_output_free(wcs_fit_output *out) nogil

import math

M_R2D = 180.0 / math.pi
M_D2R = math.pi / 180.0

def project_radec(ra, dec, ra0, dec0, proj_type):
    """Mirror projection_do_coord + projection_do_distortion from projection.c."""
    import math
    M_D2R = math.pi / 180.0
    M_R2D = 180.0 / math.pi

    xi_out, eta_out = [], []
    ra0_rad = ra0 * M_D2R
    de0_rad = dec0 * M_D2R

    for r, d in zip(ra, dec):
        # projection_do_coord
        dra = M_D2R * (r - ra0)
        drad = M_D2R * d
        sinda = math.sin(dra); cosda = math.cos(dra)
        sind0 = math.sin(de0_rad); cosd0 = math.cos(de0_rad)
        sind = math.sin(drad); cosd = math.cos(drad)
        x = +cosd * sinda
        y = -sind0 * cosd * cosda + cosd0 * sind
        rz = +cosd0 * cosd * cosda + sind0 * sind

        # projection_do_distortion
        if proj_type == 2:  # TAN
            m = 1.0 / math.sqrt(1.0 - x * x - y * y)
            x *= m; y *= m
        elif proj_type == 1:  # ARC
            d_val = math.sqrt(x * x + y * y)
            if d_val > 0 and rz < 0: m = math.asin(d_val) / d_val
            elif d_val > 0: m = (math.pi - math.asin(d_val)) / d_val
            else: m = 1.0
            x *= m; y *= m
        # SIN: nothing

        xi_out.append(x * M_R2D)
        eta_out.append(y * M_R2D)

    return xi_out, eta_out


def deproject_radec(xi, eta, ra0, dec0, proj_type):
    """Inverse of project_radec: projected (xi,eta) in deg → (RA,Dec) in deg."""
    import math
    M_D2R = math.pi / 180.0
    M_R2D = 180.0 / math.pi
    ra0_rad = ra0 * M_D2R; de0_rad = dec0 * M_D2R
    sa0 = math.sin(ra0_rad); ca0 = math.cos(ra0_rad)
    sd0 = math.sin(de0_rad); cd0 = math.cos(de0_rad)
    ra_out, de_out = [], []
    for xi_i, eta_i in zip(xi, eta):
        x = xi_i * M_D2R; y = eta_i * M_D2R
        if proj_type == 2:  # TAN
            x /= math.sqrt(1.0 + x*x + y*y)
            y /= math.sqrt(1.0 + x*x + y*y)
        elif proj_type == 1:  # ARC: forward x *= asin(d)/d; inv: d = sin(d_out), then x_in = x_out * d / d_out
            d_val = math.sqrt(x*x + y*y)
            if d_val > 0 and d_val < math.pi:
                x *= math.sin(d_val) / d_val
                y *= math.sin(d_val) / d_val
        # SIN: nothing
        z2 = 1.0 - x*x - y*y
        if z2 < 0:
            ra_out.append(float('nan')); de_out.append(float('nan'))
            continue
        z = math.sqrt(z2)  # positive z (Python forward convention, opposite of C)
        # inverse rotation: R^T where R = [[-sa0, ca0, 0], [-sd0*ca0, -sd0*sa0, cd0], [cd0*ca0, cd0*sa0, sd0]]
        x_cel = -sa0*x - sd0*ca0*y + cd0*ca0*z
        y_cel =  ca0*x - sd0*sa0*y + cd0*sa0*z
        z_cel =            cd0*y + sd0*z
        de_out.append(math.asin(max(-1.0, min(1.0, z_cel))) * M_R2D)
        ra = math.atan2(y_cel, x_cel) * M_R2D
        ra_out.append(ra % 360.0)
    return ra_out, de_out


def pix_to_proj(px, py, crpix1, crpix2, dxpoly, dypoly, order, nvar, cd11, cd12, cd21, cd22):
    """Forward: pixel → projected coords (SIP convention)."""
    xi, eta = [], []
    for i in range(len(px)):
        u = px[i] - crpix1
        v = py[i] - crpix2
        fx = eval_2d_poly_py(u, v, dxpoly, order)
        fy = eval_2d_poly_py(u, v, dypoly, order)
        xi.append(cd11 * fx + cd12 * fy)
        eta.append(cd21 * fx + cd22 * fy)
    return xi, eta


def proj_to_pix(xi, eta, crpix1, crpix2, dxpoly, dypoly, order, nvar, cd11, cd12, cd21, cd22):
    """Inverse: projected coords → pixel (SIP convention)."""
    det = cd11 * cd22 - cd12 * cd21
    px_out, py_out = [], []
    for i in range(len(xi)):
        u = (+cd22 * xi[i] - cd12 * eta[i]) / det
        v = (-cd21 * xi[i] + cd11 * eta[i]) / det
        fx = eval_2d_poly_py(u, v, dxpoly, order)
        fy = eval_2d_poly_py(u, v, dypoly, order)
        # SIP convention: prjpoly includes the identity term subtracted
        px_out.append(fx + u + crpix1)
        py_out.append(fy + v + crpix2)
    return px_out, py_out

# wcs.pxi — WCS fitting
cdef class GrtransWCS:
    """
    Fit a FITS WCS (World Coordinate System) from (RA, Dec, pixel_x, pixel_y)
    star list data.

    Usage:
        wcs = WCSFitter(ra0=189.0, dec0=62.0, projection='tan', order=3)
        wcs.fit(ra_list, dec_list, x_list, y_list)
        print(wcs.crpix1, wcs.cd11, wcs.res_prj)
        xi, eta = wcs.project(ra_list, dec_list)
        x, y = wcs.pix_to_proj(px_list, py_list)
    """

    cdef wcs_fit_input _inp
    cdef wcs_fit_output _out
    cdef int _fitted

    PROJECTION_TYPES = {'sin': 0, 'arc': 1, 'tan': 2}
    PROJ_NAMES = {0: 'SIN', 1: 'ARC', 2: 'TAN'}

    def __cinit__(self):
        memset(&self._inp, 0, sizeof(wcs_fit_input))
        memset(&self._out, 0, sizeof(wcs_fit_output))
        self._fitted = 0

    def __init__(self, ra0, dec0, projection='tan', order=3,
                 roll0=0.0, k=None, p=None):
        if order < 1:
            raise ValueError(f"order must be >= 1, got {order}")
        self._inp.ra0 = ra0
        self._inp.de0 = dec0
        self._inp.roll0 = roll0
        self._inp.type = self.PROJECTION_TYPES.get(projection, 2)
        self._inp.order = order
        if k:
            self._inp.n_k = min(len(k), 4)
            for i in range(self._inp.n_k):
                self._inp.k[i] = k[i]
        if p:
            self._inp.n_p = min(len(p), 4)
            for i in range(self._inp.n_p):
                self._inp.p[i] = p[i]

    @property
    def crpix1(self): return self._out.crpix1
    @property
    def crpix2(self): return self._out.crpix2
    @property
    def cd11(self): return self._out.cd11
    @property
    def cd12(self): return self._out.cd12
    @property
    def cd21(self): return self._out.cd21
    @property
    def cd22(self): return self._out.cd22
    @property
    def res_prj(self): return self._out.res_prj
    @property
    def res_pix(self): return self._out.res_pix
    @property
    def fit_order(self): return self._out.order

    @property
    def nvar(self): return self._out.nvar

    def get_pixpoly_dx(self):
        cdef int i
        if self._out.pixpoly_dx == NULL: return []
        return [self._out.pixpoly_dx[i] for i in range(self._out.nvar)]

    def get_pixpoly_dy(self):
        cdef int i
        if self._out.pixpoly_dy == NULL: return []
        return [self._out.pixpoly_dy[i] for i in range(self._out.nvar)]

    def get_prjpoly_dx(self):
        cdef int i
        if self._out.prjpoly_dx == NULL: return []
        return [self._out.prjpoly_dx[i] for i in range(self._out.nvar)]

    def get_prjpoly_dy(self):
        cdef int i
        if self._out.prjpoly_dy == NULL: return []
        return [self._out.prjpoly_dy[i] for i in range(self._out.nvar)]

    def fit(self, ra, dec, px, py):
        """Fit WCS from matched (RA, Dec) → (X, Y) data points."""
        cdef int n = len(ra)
        assert len(dec) == n == len(px) == len(py)

        cdef double *ra_arr  = <double*>malloc(n * sizeof(double))
        cdef double *de_arr  = <double*>malloc(n * sizeof(double))
        cdef double *px_arr  = <double*>malloc(n * sizeof(double))
        cdef double *py_arr  = <double*>malloc(n * sizeof(double))
        cdef int i
        for i in range(n):
            ra_arr[i] = ra[i]; de_arr[i] = dec[i]
            px_arr[i] = px[i]; py_arr[i] = py[i]

        wcs_fit_cy(ra_arr, de_arr, px_arr, py_arr, n, &self._inp, &self._out)
        self._fitted = 1

        free(ra_arr); free(de_arr); free(px_arr); free(py_arr)

    def project(self, ra, dec):
        """Project (RA, Dec) → projected plane coords (xi, eta) in degrees.
        Equivalent to projection_do_coord + projection_do_distortion(TAN/SIN/ARC)."""
        return project_radec(ra, dec, self._inp.ra0, self._inp.de0, self._inp.type)

    def pix_to_proj(self, px, py):
        """Forward transform: pixel (X,Y) → projected coords (xi, eta)."""
        return pix_to_proj(px, py,
            self._out.crpix1, self._out.crpix2,
            self.get_pixpoly_dx(), self.get_pixpoly_dy(),
            self._out.order, self._out.nvar,
            self._out.cd11, self._out.cd12, self._out.cd21, self._out.cd22)

    def proj_to_pix(self, xi, eta):
        """Inverse transform: projected coords (xi, eta) → pixel (X,Y)."""
        return proj_to_pix(xi, eta,
            self._out.crpix1, self._out.crpix2,
            self.get_prjpoly_dx(), self.get_prjpoly_dy(),
            self._out.order, self._out.nvar,
            self._out.cd11, self._out.cd12, self._out.cd21, self._out.cd22)

    def apply_to_radec(self, ra, dec):
        """Full pipeline: (RA, Dec) → pixel (X,Y)."""
        xi, eta = self.project(ra, dec)
        return self.proj_to_pix(xi, eta)

    def apply_to_pixel(self, px, py):
        """Full pipeline: pixel (X,Y) → (RA, Dec)."""
        xi, eta = self.pix_to_proj(px, py)
        return deproject_radec(xi, eta, self._inp.ra0, self._inp.de0, self._inp.type)

    def deproject(self, xi, eta):
        """Inverse spherical projection: projected coords (xi,eta) → (RA,Dec)."""
        return deproject_radec(xi, eta, self._inp.ra0, self._inp.de0, self._inp.type)

    def pixel_to_radec(self, px, py):
        """Full pipeline: pixel (X,Y) → (RA, Dec).  Alias for apply_to_pixel."""
        return self.apply_to_pixel(px, py)

    def to_fits_header(self):
        """Return an astropy.io.fits.Header with WCS keywords incl. SIP distortion.
        
        Includes CRPIX, CD, CTYPE, CRVAL, and forward SIP (A_ORDER/A_p_q/B_p_q) 
        computed from pixpoly.  Original pixpoly coeffs are stored in comment cards.
        """
        from astropy.io import fits as _pfts
        import numpy as np
        hdr = _pfts.Header()
        hdr['CRPIX1'] = self.crpix1
        hdr['CRPIX2'] = self.crpix2
        hdr['CD1_1']  = self.cd11
        hdr['CD1_2']  = self.cd12
        hdr['CD2_1']  = self.cd21
        hdr['CD2_2']  = self.cd22
        pix_dx = self.get_pixpoly_dx()
        pix_dy = self.get_pixpoly_dy()
        order = self.fit_order
        pname = self.PROJ_NAMES.get(self._inp.type, 'TAN')
        has_sip = order >= 2 and pix_dx and pix_dy and abs(self.cd11*self.cd22 - self.cd12*self.cd21) > 1e-30
        sip_suffix = '-SIP' if has_sip else ''
        hdr['CTYPE1'] = f'RA---{pname}{sip_suffix}'
        hdr['CTYPE2'] = f'DEC--{pname}{sip_suffix}'
        hdr['CRVAL1'] = self._inp.ra0
        hdr['CRVAL2'] = self._inp.de0
        hdr['RADESYS'] = 'ICRS'
        hdr['HISTORY'] = 'WCS fit by pyfitsh'

        if has_sip:
            hdr['A_ORDER'] = order
            hdr['B_ORDER'] = order
            # pixpoly stores coeffs by total order k: (k=0,q=0), (k=1,q=0), (k=1,q=1), (k=2,q=0)...
            idx = 0
            for k in range(order + 1):
                for q in range(k + 1):
                    p = k - q
                    cx = pix_dx[idx] if idx < len(pix_dx) else 0.0
                    cy = pix_dy[idx] if idx < len(pix_dy) else 0.0
                    if p + q >= 2:
                        import math
                        fact = math.factorial(p) * math.factorial(q)
                        hdr[f'A_{p}_{q}'] = cx / fact
                        hdr[f'B_{p}_{q}'] = cy / fact
                    idx += 1
        # Store original pixpoly as comment cards
        if pix_dx:
            dx_str = ', '.join(['%.6e' % c for c in pix_dx])
            hdr['COMMENT'] = 'FITSH PIXPOLY DX: ' + dx_str
        if pix_dy:
            dy_str = ', '.join(['%.6e' % c for c in pix_dy])
            hdr['COMMENT'] = 'FITSH PIXPOLY DY: ' + dy_str
        return hdr

    def to_dict(self):
        """Export WCS parameters as JSON-serializable dict."""
        return {
            'ra0': self._inp.ra0, 'dec0': self._inp.de0,
            'projection': self.PROJ_NAMES.get(self._inp.type, 'TAN'),
            'order': self.fit_order,
            'crpix': [self.crpix1, self.crpix2],
            'cd': [[self.cd11, self.cd12],
                   [self.cd21, self.cd22]],
            'res_prj': self.res_prj,
            'res_pix': self.res_pix,
            'pixpoly_dx': self.get_pixpoly_dx(),
            'pixpoly_dy': self.get_pixpoly_dy(),
            'prjpoly_dx': self.get_prjpoly_dx(),
            'prjpoly_dy': self.get_prjpoly_dy(),
            'k': [self._inp.k[i] for i in range(self._inp.n_k)],
            'p': [self._inp.p[i] for i in range(self._inp.n_p)],
        }

    def __dealloc__(self):
        if self._fitted:
            wcs_fit_output_free(&self._out)


# Module-level polynomial evaluation (needed by PolyFitter)

