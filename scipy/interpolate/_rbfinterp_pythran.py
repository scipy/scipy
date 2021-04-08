import numpy as np


# pythran export _distance_vector(float[:], float[:, :])
def _distance_vector(x, y):
    """
    Returns the distance between the point `x` and each point in `y`
    """
    n = y.shape[0]
    out = np.empty((n,), dtype=float)
    for i in range(n):
        out[i] = np.linalg.norm(x - y[i])

    return out


# pythran export _distance_matrix(float[:, :])
def _distance_matrix(x):
    """
    Returns the distance between each pair of points in `x`
    """
    n = x.shape[0]
    out = np.empty((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1):
            out[i, j] = np.linalg.norm(x[i] - x[j])
            out[j, i] = out[i, j]

    return out


# pythran export _polynomial_vector(float[:], int[:, :])
def _polynomial_vector(x, powers):
    """
    Returns monomials with exponents from `powers` evaluated at the point `x`
    """
    n = powers.shape[0]
    out = np.empty((n,), dtype=float)
    for i in range(n):
        out[i] = np.prod(x**powers[i])

    return out


# pythran export _polynomial_matrix(float[:, :], int[:, :])
def _polynomial_matrix(x, powers):
    """
    Returns monomials with exponents from `powers` evaluated at each point in
    `x`
    """
    n = x.shape[0]
    m = powers.shape[0]
    out = np.empty((n, m), dtype=float)
    for i in range(n):
        for j in range(m):
            out[i, j] = np.prod(x[i]**powers[j])

    return out


# pythran export _linear(float[:] or float[:, :])
def _linear(r):
    """Linear / 1st order polyharmonic spline"""
    return -r


# pythran export _tps(float[:] or float[:, :])
def _tps(r):
    """Thin plate spline / 2nd order polyharmonic spline"""
    # this is equivalent to xlogy(r**2, r)
    out = r**2*np.log(r)
    out[r == 0.0] = 0.0
    return out


# pythran export _cubic(float[:] or float[:, :])
def _cubic(r):
    """Cubic / 3rd order polyharmonic spline"""
    return r**3


# pythran export _quintic(float[:] or float[:, :])
def _quintic(r):
    """Quintic / 5th order polyharmonic spline"""
    return -r**5


# pythran export _mq(float[:] or float[:, :])
def _mq(r):
    """Multiquadratic"""
    return -np.sqrt(r**2 + 1)


# pythran export _imq(float[:] or float[:, :])
def _imq(r):
    """Inverse multiquadratic"""
    return 1/np.sqrt(r**2 + 1)


# pythran export _iq(float[:] or float[:, :])
def _iq(r):
    """Inverse quadratic"""
    return 1/(r**2 + 1)


# pythran export _ga(float[:] or float[:, :])
def _ga(r):
    """Gaussian"""
    return np.exp(-r**2)


_NAME_TO_FUNC = {
    'linear': _linear,
    'tps': _tps,
    'cubic': _cubic,
    'quintic': _quintic,
    'mq': _mq,
    'imq': _imq,
    'iq': _iq,
    'ga': _ga
    }


# pythran export _build_system(float[:, :],
#                              float[:, :] or complex[:, :],
#                              float[:],
#                              str,
#                              float,
#                              int[:, :])
def _build_system(y, d, smoothing, kernel, epsilon, powers):
    """
    Create the left-hand-side and right-hand-side of the system used to solve
    for the RBF interpolant coefficients

    Parameters
    ----------
    y : (P, N) float ndarray
        Data points coordinates
    d : (P, S) float or complex ndarray
        Data values at `y`
    smoothing : (P,) float ndarray
        Smoothing parameter for each data point
    kernel : str
        Name of the RBF
    epsilon : float
        Shape parameter
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial

    Returns
    -------
    lhs : (P + R, P + R) float ndarray
    rhs : (P + R, S) float or complex ndarray
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix
    scale : float
        Domain scaling used to create the polynomial matrix

    """
    p = d.shape[0]
    s = d.shape[1]
    r = powers.shape[0]

    yeps = y*epsilon
    kernel_func = _NAME_TO_FUNC[kernel]
    kmat = kernel_func(_distance_matrix(yeps))

    # shift and scale the polynomial domain for numerical stability
    shift = np.mean(y, axis=0)
    if p > 1:
        scale = np.ptp(y, axis=0).max()
    else:
        scale = 1.0

    yhat = (y - shift)/scale
    pmat = _polynomial_matrix(yhat, powers)

    lhs = np.empty((p + r, p + r), dtype=float)
    lhs[:p, :p] = kmat
    lhs[:p, p:] = pmat
    lhs[p:, :p] = pmat.T
    lhs[p:, p:] = 0.0
    for i in range(p):
        lhs[i, i] += smoothing[i]

    rhs = np.empty((p + r, s), dtype=d.dtype)
    rhs[:p] = d
    rhs[p:] = 0.0

    return lhs, rhs, shift, scale


# pythran export _evaluate(float[:, :],
#                          float[:, :],
#                          str,
#                          float,
#                          int[:, :],
#                          float[:],
#                          float,
#                          float[:, :] or complex[:, :])
def _evaluate(x, y, kernel, epsilon, powers, shift, scale, coeffs):
    '''
    Evaluates the RBF interpolant at `x`

    Parameters
    ----------
    x : (Q, N) float ndarray
        Interpolation point coordinates
    y : (P, N) float ndarray
        Data point coordinates
    kernel : str
        Name of the RBF
    epsilon : float
        Shape parameter
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial
    shift : (N,) float ndarray
        Shift the polynomial domain for numerical stability
    scale : float
        Scale the polynomial domain for numerical stability
    coeffs : (P + R, S) float or complex ndarray
        RBF and polynomial coefficients

    Returns
    -------
    (Q, S) float or complex ndarray

    '''
    q = x.shape[0]
    p = y.shape[0]
    s = coeffs.shape[1]

    yeps = y*epsilon
    kernel_func = _NAME_TO_FUNC[kernel]
    out = np.empty((q, s), dtype=coeffs.dtype)
    for i in range(q):
        xeps = x[i]*epsilon
        kvec = kernel_func(_distance_vector(xeps, yeps))

        xhat = (x[i] - shift)/scale
        pvec = _polynomial_vector(xhat, powers)

        for j in range(s):
            out[i, j] = np.sum(coeffs[:p, j]*kvec)
            out[i, j] += np.sum(coeffs[p:, j]*pvec)

    return out
