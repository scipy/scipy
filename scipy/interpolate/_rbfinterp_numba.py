"""RBF backend with numba JIT compiled evaluations.
"""
import numba
import numpy as np


from ._rbfinterp_np import *

# bring the _undersored names, too
from ._rbfinterp_np import (
    _build_and_solve_system, _build_evaluation_coefficients, _build_system,
    _monomial_powers, _monomial_powers_impl, polynomial_matrix
)


### Copy-paste from _rbfinterp_pythran.py

@numba.njit
def linear(r):
    return -r


@numba.njit
def thin_plate_spline(r):
    if r == 0:
        return 0.0
    else:
        return r**2*np.log(r)


@numba.njit
def cubic(r):
    return r**3


@numba.njit
def quintic(r):
    return -r**5


@numba.njit
def multiquadric(r):
    return -np.sqrt(r**2 + 1)


@numba.njit
def inverse_multiquadric(r):
    return 1/np.sqrt(r**2 + 1)


@numba.njit
def inverse_quadratic(r):
    return 1/(r**2 + 1)


@numba.njit
def gaussian(r):
    return np.exp(-r**2)


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


# copy-paste from _rbfinterp_np.py / _rbfinterp_pythran.py

def _build_evaluation_coefficients(x, y, kernel, epsilon, powers, shift, scale):
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
    return _build_evaluation_coefficients_impl(x, y, kernel_func, epsilon, powers, shift, scale)


@numba.njit
def kernel_vector(x, y, kernel_func, out):
    """Evaluate RBFs, with centers at `y`, at the point `x`."""
    # NB: explicit prange
    for i in numba.prange(y.shape[0]):
        out[i] = kernel_func(np.linalg.norm(x - y[i]))


@numba.njit
def polynomial_vector(x, powers, out):
    """Evaluate monomials, with exponents from `powers`, at the point `x`."""
    # NB: explicit prange
    for i in numba.prange(powers.shape[0]):
        out[i] = np.prod(x**powers[i])


@numba.jit(nopython=True, parallel=True)
def _build_evaluation_coefficients_impl(x, y, kernel_func, epsilon, powers, shift, scale):
    q = x.shape[0]
    p = y.shape[0]
    r = powers.shape[0]

    yeps = y*epsilon
    xeps = x*epsilon
    xhat = (x - shift)/scale

    vec = np.empty((q, p + r), dtype=float)

    # NB: explicit prange
    for i in numba.prange(q):
        kernel_vector(xeps[i], yeps, kernel_func, vec[i, :p])
        polynomial_vector(xhat[i], powers, vec[i, p:])

    return vec


def compute_interpolation_impl(x, y, kernel_func, epsilon, powers, shift, scale, coeffs):
    vec = _build_evaluation_coefficients_impl(x, y, kernel_func, epsilon, powers, shift, scale)
    return vec @ coeffs


def compute_interpolation(x, y, kernel, epsilon, powers, shift, scale, coeffs, xp):
    kernel_func = NAME_TO_FUNC[kernel]
    return compute_interpolation_impl(x, y, kernel_func, epsilon, powers, shift, scale, coeffs)


# XXX: apparently not used. Was added to try working around numba failing to
#      jit compile np.linal.norm

@numba.njit
def _norm(x):
    # numba 0.61.2 fails to compile np.linalg.vector_norm or np.linalg.norm(x, axis=-1)
    # Thus inline from https://github.com/numpy/numpy/blob/v2.2.0/numpy/linalg/_linalg.py#L2779C1-L2781C69
    s = (x * x)
    return np.sqrt(np.sum(s, axis=-1))
