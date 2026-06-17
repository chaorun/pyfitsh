# poly_helpers.pxi — polynomial evaluation helper
def eval_2d_poly_py(x, y, coeff, order):
    """Equivalent to C eval_2d_poly(x, y, order, coeff, 0, 0, 1)."""
    idx = 0; result = 0.0
    for k in range(order + 1):
        for j in range(k + 1):
            term = 1.0
            for p in range(1, k - j + 1): term *= x / p
            for p in range(1, j + 1): term *= y / p
            result += coeff[idx] * term
            idx += 1
    return result

