""" Internal primitives for RBF interpolation.

Routines here are of dual purpose:
- AOT compiled by pythran (hence # pythran export comments)
- used for the "generic", Array API agnostic backend with non-numpy array libraries,
  See `_rbfinterp_xp.py`.
"""

def linear(r, xp):
    return -r


def thin_plate_spline(r, xp):
    return xp.where(r == 0, 0, r**2 * xp.log(r))


def cubic(r, xp):
    return r**3


def quintic(r, xp):
    return -r**5


def multiquadric(r, xp):
    return -xp.sqrt(r**2 + 1)


def inverse_multiquadric(r, xp):
    return 1.0 / xp.sqrt(r**2 + 1.0)


def inverse_quadratic(r, xp):
    return 1.0 / (r**2 + 1.0)


def gaussian(r, xp):
    return xp.exp(-r**2)


NAME_TO_FUNC = {
   "linear": linear,
   "thin_plate_spline": thin_plate_spline,
   "cubic": cubic,
   "quintic": quintic,
   "multiquadric": multiquadric,
   "inverse_multiquadric": inverse_multiquadric,
   "inverse_quadratic": inverse_quadratic,
   "gaussian": gaussian
   }


# pythran export polynomial_matrix(float[:, :], int64[:, :], numpy pkg)
def polynomial_matrix(x, powers, xp):
    """Evaluate monomials, with exponents from `powers`, at `x`."""
    return xp.prod(x[:, None, :] ** powers, axis=-1)


# NB: exported only to use in the tests
# pythran export kernel_matrix_at_centers(float[:, :], str, numpy pkg)
def kernel_matrix_at_centers(x, kernel, xp):
    """Return RBFs, with centers at `x`, evaluated at `x`."""
    kernel_func = NAME_TO_FUNC[kernel]
    return kernel_matrix(x, kernel_func, xp)


def kernel_matrix(x, kernel_func, xp):
    """Evaluate RBFs, with centers at `x`, at `x`."""
    return _kernel_matrix_impl(x, x, kernel_func, xp)


def _kernel_matrix_impl(x, y, kernel_func, xp):
    """Evaluate RBFs, with centerx at `y`, at `x`."""
    return kernel_func(
        xp.linalg.vector_norm(x[None, :, :] - y[:, None, :], axis=-1), xp
    )


# pythran export _build_system(float[:, :],
#                              float[:, :],
#                              float[:],
#                              str,
#                              float,
#                              int64[:, :],
#                              numpy pkg)
def _build_system(y, d, smoothing, kernel, epsilon, powers, xp):
    """Build the system used to solve for the RBF interpolant coefficients.

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
        Left-hand side matrix.
    rhs : (P + R, S) float ndarray
        Right-hand side matrix.
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix.
    scale : (N,) float ndarray
        Domain scaling used to create the polynomial matrix.

    """
    s = d.shape[1]
    r = powers.shape[0]
    kernel_func = NAME_TO_FUNC[kernel]

    # Shift and scale the polynomial domain to be between -1 and 1
    mins = xp.min(y, axis=0)
    maxs = xp.max(y, axis=0)
    shift = (maxs + mins)/2
    scale = (maxs - mins)/2
    # The scale may be zero if there is a single point or all the points have
    # the same value for some dimension. Avoid division by zero by replacing
    # zeros with ones.
    scale = xp.where(scale == 0.0, 1.0, scale)

    yeps = y*epsilon
    yhat = (y - shift)/scale

    out_kernels  = kernel_matrix(yeps, kernel_func, xp)
    out_poly = polynomial_matrix(yhat, powers, xp)

    lhs = xp.concat(
        [
            xp.concat((out_kernels, out_poly), axis=1),
            xp.concat((out_poly.T, xp.zeros((r, r))), axis=1)
        ]
    , axis=0) + xp.diag(xp.concat([smoothing, xp.zeros(r)]))

    rhs = xp.concat([d, xp.zeros((r, s))], axis=0)

    return lhs, rhs, shift, scale


# pythran export _build_evaluation_coefficients(float[:, :],
#                          float[:, :],
#                          str,
#                          float,
#                          int64[:, :],
#                          float[:],
#                          float[:],
#                          numpy pkg)
def _build_evaluation_coefficients(x, y, kernel, epsilon, powers,
                                   shift, scale, xp):
    """Construct the coefficients needed to evaluate
    the RBF.

    Parameters
    ----------
    x : (Q, N) float ndarray
        Evaluation point coordinates.
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

    Returns
    -------
    (Q, P + R) float ndarray

    """
    kernel_func = NAME_TO_FUNC[kernel]

    yeps = y*epsilon
    xeps = x*epsilon
    xhat = (x - shift)/scale

    vec = xp.concat(
        [
            _kernel_matrix_impl(yeps, xeps, kernel_func, xp),
            xp.prod(xhat[:, None, :] ** powers, axis=-1)
        ], axis=-1
    )

    return vec

