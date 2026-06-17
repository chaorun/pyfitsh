# fiphot_internals.pxi — fiphot internal helper functions

cdef tuple fiphot_build_kernel_arrays(dict kernel_dict):
    """Flatten kernel_dict to C arrays for fiphot_photometry_from_raw_cy."""
    import numpy as np
    g = kernel_dict.get('global', {})
    kds = kernel_dict['kernels']
    nk = len(kds)
    ox = g.get('offset', [0.0, 0.0])[0]
    oy = g.get('offset', [0.0, 0.0])[1]
    sc = g.get('scale', 1.0)
    kt = g.get('type', 0)

    k_types_v = np.zeros(nk, dtype=np.int32)
    k_orders_v = np.zeros(nk, dtype=np.int32)
    k_ncoeffs_v = np.zeros(nk, dtype=np.int32)
    k_hsizes_v = np.zeros(nk, dtype=np.int32)
    k_sigmas_v = np.zeros(nk, dtype=np.float64)
    k_bx_v = np.zeros(nk, dtype=np.int32)
    k_by_v = np.zeros(nk, dtype=np.int32)

    total_coeffs = sum(len(kd.get('coeff', [])) for kd in kds)
    k_coeffs_flat_v = np.zeros(total_coeffs, dtype=np.float64)

    ci = 0
    cdef list k_coeff_offsets = []
    for i, kd in enumerate(kds):
        k_types_v[i] = kd.get('type', 0)
        k_orders_v[i] = kd.get('order', 0)
        clist = kd.get('coeff', [])
        nv = len(clist)
        k_ncoeffs_v[i] = nv
        k_hsizes_v[i] = kd.get('hsize', 0)
        k_sigmas_v[i] = kd.get('sigma', 0.0)
        k_bx_v[i] = kd.get('bx', 0)
        k_by_v[i] = kd.get('by', 0)
        for cv in clist:
            k_coeffs_flat_v[ci] = cv
            ci += 1
        k_coeff_offsets.append(ci - nv)

    return (nk, ox, oy, sc, kt,
            k_types_v, k_orders_v, k_ncoeffs_v,
            k_hsizes_v, k_sigmas_v,
            k_bx_v, k_by_v,
            k_coeffs_flat_v, k_coeff_offsets)
