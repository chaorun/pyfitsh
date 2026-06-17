# wcs_helpers.pxi — WCS projection helpers
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
