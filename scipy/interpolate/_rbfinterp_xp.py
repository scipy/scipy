import numpy as np


from ._rbfinterp_common import _monomial_powers_impl

def _monomial_powers(ndim, degree, xp):
    out = _monomial_powers_impl(ndim, degree)
    out = xp.asarray(out)
    if out.shape[0] == 0:
        out = xp.reshape(out, (0, ndim))
    return out


def _build_and_solve_system(y, d, smoothing, kernel, epsilon, powers, xp):
    """Build and solve the RBF interpolation system of equations.

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
    coeffs : (P + R, S) float ndarray
        Coefficients for each RBF and monomial.
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix.
    scale : (N,) float ndarray
        Domain scaling used to create the polynomial matrix.

    """
    lhs, rhs, shift, scale = _build_system(
        y, d, smoothing, kernel, epsilon, powers, xp
        )

    coeffs = xp.linalg.solve(lhs, rhs)

    '''
    # XXX cannot give same diagnostic in a backend-agnostic way?
    _, _, coeffs, info = dgesv(lhs, rhs, overwrite_a=True, overwrite_b=True)
    if info < 0:
        raise ValueError(f"The {-info}-th argument had an illegal value.")
    elif info > 0:
        msg = "Singular matrix."
        nmonos = powers.shape[0]
        if nmonos > 0:
            pmat = _polynomial_matrix((y - shift)/scale, powers)
            rank = np.linalg.matrix_rank(pmat)
            if rank < nmonos:
                msg = (
                    "Singular matrix. The matrix of monomials evaluated at "
                    "the data point coordinates does not have full column "
                    f"rank ({rank}/{nmonos})."
                    )

        raise LinAlgError(msg)
    '''

    return shift, scale, coeffs



def linear(r, xp):
    return -r


# NB: changed w.r.t. pythran
def thin_plate_spline(r, xp):
 #   if r == 0:
 #       return 0.0
 #   else:
 #       return r**2 * xp.log(r)
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


def kernel_vector(x, y, kernel_func, out, xp):
    """Evaluate RBFs, with centers at `y`, at the point `x`."""
    for i in range(y.shape[0]):
        out[i] = kernel_func(xp.linalg.norm(x - y[i]), xp)


def polynomial_vector(x, powers, out, xp):
    """Evaluate monomials, with exponents from `powers`, at the point `x`."""
    for i in range(powers.shape[0]):
        out[i] = xp.prod(x**powers[i])


def kernel_matrix(x, kernel_func, out, xp):
    """Evaluate RBFs, with centers at `x`, at `x`."""
    for i in range(x.shape[0]):
        for j in range(i+1):
            out[i, j] = kernel_func(xp.linalg.norm(x[i] - x[j]), xp)
            out[j, i] = out[i, j]


def polynomial_matrix(x, powers, out, xp):
    """Evaluate monomials, with exponents from `powers`, at `x`."""
    for i in range(x.shape[0]):
        for j in range(powers.shape[0]):
            out[i, j] = xp.prod(x[i]**powers[j])


# pythran export _kernel_matrix(float[:, :], str)
def _kernel_matrix(x, kernel, xp):
    """Return RBFs, with centers at `x`, evaluated at `x`."""
    out = xp.empty((x.shape[0], x.shape[0]), dtype=xp.float64)
    kernel_func = NAME_TO_FUNC[kernel]
    kernel_matrix(x, kernel_func, out, xp)
    return out


# pythran export _polynomial_matrix(float[:, :], int[:, :])
def _polynomial_matrix(x, powers, xp):
    """Return monomials, with exponents from `powers`, evaluated at `x`."""
    out = xp.empty((x.shape[0], powers.shape[0]), dtype=xp.float64)
    polynomial_matrix(x, powers, out, xp)
    return out


# pythran export _build_system(float[:, :],
#                              float[:, :],
#                              float[:],
#                              str,
#                              float,
#                              int[:, :])
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
    p = d.shape[0]
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
    scale[scale == 0.0] = 1.0

    yeps = y*epsilon
    yhat = (y - shift)/scale

    # Transpose to make the array fortran contiguous. This is required for
    # dgesv to not make a copy of lhs.
    lhs = xp.empty((p + r, p + r), dtype=xp.float64).T
    kernel_matrix(yeps, kernel_func, lhs[:p, :p], xp)
    polynomial_matrix(yhat, powers, lhs[:p, p:], xp)
    lhs[p:, :p] = lhs[:p, p:].T
    lhs[p:, p:] = 0.0
    for i in range(p):
        lhs[i, i] += smoothing[i]

    # Transpose to make the array fortran contiguous.
    rhs = xp.empty((s, p + r), dtype=xp.float64).T
    rhs[:p] = d
    rhs[p:] = 0.0

    return lhs, rhs, shift, scale


# pythran export _build_evaluation_coefficients(float[:, :],
#                          float[:, :],
#                          str,
#                          float,
#                          int[:, :],
#                          float[:],
#                          float[:])
import torch

torch._dynamo.config.cache_size_limit = 160

@torch.compile(fullgraph=True, dynamic=True)
def _build_evaluation_coefficients(x, y, kernel, epsilon, powers, shift, scale, xp):
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
    q = x.shape[0]
    p = y.shape[0]
    r = powers.shape[0]
    kernel_func = NAME_TO_FUNC[kernel]

    yeps = y*epsilon
    xeps = x*epsilon
    xhat = (x - shift)/scale

    # NB: changed w.r.t. pythran
###    vec = _coeffs_inner(xeps, xhat, yeps, powers, kernel_func, xp)
    vec = xp.empty((q, p + r), dtype=xp.float64)
    vec[:, :p] = kernel_func(xp.linalg.vector_norm(xeps[:, None, :] - yeps[None, :, :], axis=-1), xp)
    vec[:, p:] = xp.prod(xhat[:, None, :] ** powers, axis=-1)

#    for i in range(q):
#        kernel_vector(xeps[i], yeps, kernel_func, vec[i, :p], xp)
#        polynomial_vector(xhat[i], powers, vec[i, p:], xp)

    '''
    from numpy.testing import assert_allclose
    assert_allclose(vec[:, p:], xp.prod(xhat[:, None, :] ** powers, axis=-1), atol=1e-15)
    '''

    '''
    vec2 = xp.empty((q, p), dtype=xp.float64)
    for i in range(x.shape[0]):
        for j in range(y.shape[0]):
            vec2[i, j] = kernel_func(xp.linalg.vector_norm(xeps[i] - yeps[j]), xp)


    vec3 = kernel_func(xp.linalg.vector_norm(xeps[:, None, :] - yeps[None, :, :], axis=-1), xp)

    from numpy.testing import assert_allclose
    assert_allclose(vec[:, :p], vec2)

    assert_allclose(vec2, vec3)
    '''

    return vec


import torch

torch._dynamo.config.cache_size_limit = 160

@torch.compile(fullgraph=True, dynamic=True)
def _coeffs_inner(xeps, xhat, yeps, powers, kernel_func, xp):
    q = xeps.shape[0]
    p = yeps.shape[0]
    r = powers.shape[0]

    vec = xp.empty((q, p + r), dtype=xp.float64)
    vec[:, :p] = kernel_func(xp.linalg.vector_norm(xeps[:, None, :] - yeps[None, :, :], axis=-1), xp)
    vec[:, p:] = xp.prod(xhat[:, None, :] ** powers, axis=-1)

    return vec
