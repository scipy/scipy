import numpy as np


def linear(r):
    """Linear / 1st order polyharmonic spline"""
    return -r


def tps(r):
    """Thin plate spline / 2nd order polyharmonic spline"""
    if r == 0:
        return 0.0
    else:
        return r**2*np.log(r)


def cubic(r):
    """Cubic / 3rd order polyharmonic spline"""
    return r**3


def quintic(r):
    """Quintic / 5th order polyharmonic spline"""
    return -r**5


def mq(r):
    """Multiquadric"""
    return -np.sqrt(r**2 + 1)


def imq(r):
    """Inverse multiquadric"""
    return 1/np.sqrt(r**2 + 1)


def iq(r):
    """Inverse quadratic"""
    return 1/(r**2 + 1)


def ga(r):
    """Gaussian"""
    return np.exp(-r**2)


NAME_TO_FUNC = {
   'linear': linear,
   'tps': tps,
   'cubic': cubic,
   'quintic': quintic,
   'mq': mq,
   'imq': imq,
   'iq': iq,
   'ga': ga
   }


def kernel_vector(x, y, kernel_func, out):
    """
    Returns RBFs with centers at `y` evaluated at the single point `x`.
    """
    for i in range(y.shape[0]):
        out[i] = kernel_func(np.linalg.norm(x - y[i]))


def polynomial_vector(x, powers, out):
    """
    Returns monomials with exponents from `powers` evaluated at the single
    point `x`.
    """
    for i in range(powers.shape[0]):
        out[i] = np.prod(x**powers[i])


def kernel_matrix(x, kernel_func, out):
    """
    Returns RBFs with centers at `x` evaluated at `x`.
    """
    for i in range(x.shape[0]):
        for j in range(i+1):
            out[i, j] = kernel_func(np.linalg.norm(x[i] - x[j]))
            out[j, i] = out[i, j]


def polynomial_matrix(x, powers, out):
    """
    Returns monomials with exponents from `powers` evaluated at the points `x`.
    """
    for i in range(x.shape[0]):
        for j in range(powers.shape[0]):
            out[i, j] = np.prod(x[i]**powers[j])


# pythran export _kernel_matrix(float[:, :], str)
def _kernel_matrix(x, kernel):
    """
    Returns RBFs with centers at `x` evaluated at `x`.
    """
    out = np.empty((x.shape[0], x.shape[0]), dtype=float)
    kernel_func = NAME_TO_FUNC[kernel]
    kernel_matrix(x, kernel_func, out)
    return out


# pythran export _polynomial_matrix(float[:, :], int[:, :])
def _polynomial_matrix(x, powers):
    """
    Returns monomials with exponents from `powers` evaluated at the points `x`.
    """
    out = np.empty((x.shape[0], powers.shape[0]), dtype=float)
    polynomial_matrix(x, powers, out)
    return out


# pythran export _build_system(float[:, :],
#                              float[:, :],
#                              float[:],
#                              str,
#                              float,
#                              int[:, :])
def _build_system(y, d, smoothing, kernel, epsilon, powers):
    """
    Create the left-hand-side and right-hand-side of the system used to solve
    for the RBF interpolant coefficients.

    Parameters
    ----------
    y : (P, N) float ndarray
        Data point coordinates.
    d : (P, S) float ndarray
        Data values at `y`.
    smoothing : (P,) float ndarray
        Smoothing parameter for each data point.
    kernel : str
        Name of the RBF.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.

    Returns
    -------
    lhs : (P + R, P + R) float ndarray
    rhs : (P + R, S) float ndarray
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix.
    scale : (N,) float ndarray
        Domain scaling used to create the polynomial matrix.

    """
    p = d.shape[0]
    s = d.shape[1]
    r = powers.shape[0]
    kernel_func = NAME_TO_FUNC[kernel]

    # Shift and scale the polynomial domain to the unit hypercube.
    mins = np.min(y, axis=0)
    maxs = np.max(y, axis=0)
    shift = (maxs + mins)/2
    scale = (maxs - mins)/2
    # The scale may be zero if there is a single point or all the points have
    # the same value for some dimension. Avoid division by zero by replacing
    # zeros with ones.
    scale[scale == 0.0] = 1.0

    yeps = y*epsilon
    yhat = (y - shift)/scale

    # Transpose to make the array fortran contiguous. This is required for
    # dgesv to not make a copy of lhs.
    lhs = np.empty((p + r, p + r), dtype=float).T
    kernel_matrix(yeps, kernel_func, lhs[:p, :p])
    polynomial_matrix(yhat, powers, lhs[:p, p:])
    lhs[p:, :p] = lhs[:p, p:].T
    lhs[p:, p:] = 0.0
    for i in range(p):
        lhs[i, i] += smoothing[i]

    # Transpose to make the array fortran contiguous.
    rhs = np.empty((s, p + r), dtype=float).T
    rhs[:p] = d
    rhs[p:] = 0.0

    return lhs, rhs, shift, scale


# pythran export _evaluate(float[:, :],
#                          float[:, :],
#                          str,
#                          float,
#                          int[:, :],
#                          float[:],
#                          float[:],
#                          float[:, :])
def _evaluate(x, y, kernel, epsilon, powers, shift, scale, coeffs):
    """
    Evaluates the RBF interpolant at `x`.

    Parameters
    ----------
    x : (Q, N) float ndarray
        Interpolation point coordinates.
    y : (P, N) float ndarray
        Data point coordinates.
    kernel : str
        Name of the RBF.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.
    shift : (N,) float ndarray
        Shifts the polynomial domain for numerical stability.
    scale : (N,) float ndarray
        Scales the polynomial domain for numerical stability.
    coeffs : (P + R, S) float ndarray
        RBF and polynomial coefficients.

    Returns
    -------
    (Q, S) float ndarray

    """
    q = x.shape[0]
    p = y.shape[0]
    r = powers.shape[0]
    s = coeffs.shape[1]
    kernel_func = NAME_TO_FUNC[kernel]

    yeps = y*epsilon
    xeps = x*epsilon
    xhat = (x - shift)/scale

    out = np.empty((q, s), dtype=float)
    vec = np.empty((p + r,), dtype=float)
    for i in range(q):
        kernel_vector(xeps[i], yeps, kernel_func, vec[:p])
        polynomial_vector(xhat[i], powers, vec[p:])
        for j in range(s):
            out[i, j] = np.dot(coeffs[:, j], vec)

    return out

