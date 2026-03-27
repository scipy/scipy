import numpy as np

# pythran export capsule linear(float64)
def linear(r):
    return -r

# pythran export capsule thin_plate_spline(float64)
def thin_plate_spline(r):
    if r == 0:
        return 0.0
    else:
        return r**2*np.log(r)

# pythran export capsule cubic(float64)
def cubic(r):
    return r**3

# pythran export capsule quintic(float64)
def quintic(r):
    return -r**5

# pythran export capsule multiquadric(float64)
def multiquadric(r):
    return -np.sqrt(r**2 + 1)

# pythran export capsule inverse_multiquadric(float64)
def inverse_multiquadric(r):
    return 1/np.sqrt(r**2 + 1)

# pythran export capsule inverse_quadratic(float64)
def inverse_quadratic(r):
    return 1/(r**2 + 1)

# pythran export capsule gaussian(float64)
def gaussian(r):
    return np.exp(-r**2)

# pythran export capsule matern1_2(float64)
def matern1_2(r):
    return np.exp(-r)


# pythran export capsule matern3_2(float64)
def matern3_2(r):
    term = np.sqrt(3.0) * r
    return (1.0 + term) * np.exp(-term)


# pythran export capsule matern5_2(float64)
def matern5_2(r):
    term = np.sqrt(5.0) * r
    return (1.0 + term + 5.0 * r**2 / 3.0) * np.exp(-term)


NAME_TO_FUNC = {
   "linear": linear,
   "thin_plate_spline": thin_plate_spline,
   "cubic": cubic,
   "quintic": quintic,
   "multiquadric": multiquadric,
   "inverse_multiquadric": inverse_multiquadric,
   "inverse_quadratic": inverse_quadratic,
   "gaussian": gaussian,
   "matern1_2": matern1_2,
   "matern3_2": matern3_2,
   "matern5_2": matern5_2
}


def kernel_vector(x, y, kernel_func, out):
    """Evaluate RBFs, with centers at `y`, at the point `x`."""
    for i in range(y.shape[0]):
        out[i] = kernel_func(np.linalg.norm(x - y[i]))


def polynomial_vector(x, powers, out):
    """Evaluate monomials, with exponents from `powers`, at the point `x`."""
    for i in range(powers.shape[0]):
        out[i] = np.prod(x**powers[i])


def kernel_matrix(x, kernel_func, out):
    """Evaluate RBFs, with centers at `x`, at `x`."""
    for i in range(x.shape[0]):
        for j in range(i+1):
            out[i, j] = kernel_func(np.linalg.norm(x[i] - x[j]))
            out[j, i] = out[i, j]


def polynomial_matrix(x, powers, out):
    """Evaluate monomials, with exponents from `powers`, at `x`."""
    for i in range(x.shape[0]):
        for j in range(powers.shape[0]):
            out[i, j] = np.prod(x[i]**powers[j])


# pythran export _kernel_matrix(float[:, :], str)
def _kernel_matrix(x, kernel):
    """Return RBFs, with centers at `x`, evaluated at `x`."""
    out = np.empty((x.shape[0], x.shape[0]), dtype=float)
    kernel_func = NAME_TO_FUNC[kernel]
    kernel_matrix(x, kernel_func, out)

    return out


# pythran export _polynomial_matrix(float[:, :], int64[:, :])
def _polynomial_matrix(x, powers):
    """Return monomials, with exponents from `powers`, evaluated at `x`."""
    out = np.empty((x.shape[0], powers.shape[0]), dtype=float)
    polynomial_matrix(x, powers, out)
    return out


# pythran export _build_system_with_kernel(float[:, :],
#                                          float[:, :],
#                                          float[:],
#                                          float64(float64),
#                                          float,
#                                          int64[:, :])
def _build_system_with_kernel(y, d, smoothing, kernel_func, epsilon, powers):
    """Build the system used to solve for the RBF interpolant coefficients.

    Parameters
    ----------
    y : (P, N) float ndarray
        Data point coordinates.
    d : (P, S) float ndarray
        Data values at `y`.
    smoothing : (P,) float ndarray
        Smoothing parameter for each data point.
    kernel_func : float64(float64) capsule
        Compiled RBF kernel: maps scalar distance r to scalar value.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.

    Returns
    -------
    lhs : (P + R, P + R) float ndarray
        Left-hand side matrix.
    rhs : (P + R, S) float ndarray
        Right-hand side matrix.
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix.
    scale : (N,) float ndarray
        Domain scaling used to create the polynomial matrix.
    """
    p = d.shape[0]
    s = d.shape[1]
    r = powers.shape[0]

    mins = np.min(y, axis=0)
    maxs = np.max(y, axis=0)
    shift = (maxs + mins) / 2
    scale = (maxs - mins) / 2
    scale[scale == 0.0] = 1.0

    yeps = y * epsilon
    yhat = (y - shift) / scale

    lhs = np.empty((p + r, p + r), dtype=float).T
    kernel_matrix(yeps, kernel_func, lhs[:p, :p])
    polynomial_matrix(yhat, powers, lhs[:p, p:])
    lhs[p:, :p] = lhs[:p, p:].T
    lhs[p:, p:] = 0.0
    for i in range(p):
        lhs[i, i] += smoothing[i]

    rhs = np.empty((s, p + r), dtype=float).T
    rhs[:p] = d
    rhs[p:] = 0.0

    return lhs, rhs, shift, scale


# pythran export _build_evaluation_coefficients_with_kernel(float[:, :],
#                                                           float[:, :],
#                                                           float64(float64),
#                                                           float,
#                                                           int64[:, :],
#                                                           float[:],
#                                                           float[:])
def _build_evaluation_coefficients_with_kernel(x, y, kernel_func, epsilon,
                                               powers, shift, scale):
    """
    Construct the coefficients needed to evaluate the RBF.

    Parameters
    ----------
    x : (Q, N) float ndarray
        Evaluation point coordinates.
    y : (P, N) float ndarray
        Data point coordinates.
    kernel_func : float64(float64) capsule
        Compiled RBF kernel: maps scalar distance r to scalar value.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.
    shift : (N,) float ndarray
        Shifts the polynomial domain for numerical stability.
    scale : (N,) float ndarray
        Scales the polynomial domain for numerical stability.

    Returns
    -------
    (Q, P + R) float ndarray
    """
    q = x.shape[0]
    p = y.shape[0]
    r = powers.shape[0]

    yeps = y * epsilon
    xeps = x * epsilon
    xhat = (x - shift) / scale

    vec = np.empty((q, p + r), dtype=float)
    for i in range(q):
        kernel_vector(xeps[i], yeps, kernel_func, vec[i, :p])
        polynomial_vector(xhat[i], powers, vec[i, p:])

    return vec
