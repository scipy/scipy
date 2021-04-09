import numpy as np


def _linear(r):
    """Linear / 1st order polyharmonic spline"""
    return -r


def _tps(r):
    """Thin plate spline / 2nd order polyharmonic spline"""
    if r == 0:
        return 0.0
    else:
        return r**2*np.log(r)


def _cubic(r):
    """Cubic / 3rd order polyharmonic spline"""
    return r**3


def _quintic(r):
    """Quintic / 5th order polyharmonic spline"""
    return -r**5


def _mq(r):
    """Multiquadratic"""
    return -np.sqrt(r**2 + 1)


def _imq(r):
    """Inverse multiquadratic"""
    return 1/np.sqrt(r**2 + 1)


def _iq(r):
    """Inverse quadratic"""
    return 1/(r**2 + 1)


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


def _kernel_func_vector(x, y, kernel_func):
    """
    Returns RBFs with centers at `y` evaluated at the single point `x`
    """
    n = y.shape[0]
    out = np.empty((n,), dtype=float)
    for i in range(n):
        out[i] = kernel_func(np.linalg.norm(x - y[i]))

    return out


def _polynomial_vector(x, powers):
    """
    Returns monomials with exponents from `powers` evaluated at the single
    point `x`
    """
    n = powers.shape[0]
    out = np.empty((n,), dtype=float)
    for i in range(n):
        out[i] = np.prod(x**powers[i])

    return out


def _kernel_func_matrix(x, kernel_func):
    """
    Returns RBFs with centers at `x` evaluated at `x`
    """
    n = x.shape[0]
    out = np.empty((n, n), dtype=float)
    for i in range(n):
        for j in range(i+1):
            out[i, j] = kernel_func(np.linalg.norm(x[i] - x[j]))
            out[j, i] = out[i, j]

    return out


# pythran export _kernel_matrix(float[:, :], str)
def _kernel_matrix(x, kernel):
    '''
    Returns RBFs with centers at `x` evaluated at `x`. The RBF is specified as
    a string
    '''
    kernel_func = _NAME_TO_FUNC[kernel]
    return _kernel_func_matrix(x, kernel_func)


# pythran export _polynomial_matrix(float[:, :], int[:, :])
def _polynomial_matrix(x, powers):
    """
    Returns monomials with exponents from `powers` evaluated at the points `x`
    """
    n = x.shape[0]
    m = powers.shape[0]
    out = np.empty((n, m), dtype=float)
    for i in range(n):
        for j in range(m):
            out[i, j] = np.prod(x[i]**powers[j])

    return out


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
    kmat = _kernel_func_matrix(yeps, kernel_func)

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
        kvec = _kernel_func_vector(xeps, yeps, kernel_func)

        xhat = (x[i] - shift)/scale
        pvec = _polynomial_vector(xhat, powers)

        for j in range(s):
            out[i, j] = np.sum(coeffs[:p, j]*kvec)
            out[i, j] += np.sum(coeffs[p:, j]*pvec)

    return out
