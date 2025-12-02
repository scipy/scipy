"""'Generic' Array API backend for RBF interpolation."""
from numpy.linalg import LinAlgError
from ._rbfinterp_common import _monomial_powers_impl

from ._rbfinterp_pythran_src import (
    _build_evaluation_coefficients,
    _build_system,
    _kernel_matrix_impl,
    cubic,
    gaussian,
    inverse_multiquadric,
    inverse_quadratic,
    kernel_matrix,
    kernel_matrix_at_centers,
    linear,
    multiquadric,
    polynomial_matrix,
    quintic,
    thin_plate_spline
)

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
    try:
        coeffs = xp.linalg.solve(lhs, rhs)
    except Exception:
        # Best-effort attempt to emit a helpful message.
        # `_rbfinterp_np` backend gives better diagnostics; it is hard to
        # match it in a backend-agnostic way: e.g. jax emits no error at all,
        # and instead returns an array of nans for a singular `lhs`.
        msg = "Singular matrix"
        nmonos = powers.shape[0]
        if nmonos > 0:
            pmat = polynomial_matrix((y - shift)/scale, powers, xp=xp)
            rank = xp.linalg.matrix_rank(pmat)
            if rank < nmonos:
                msg = (
                    "Singular matrix. The matrix of monomials evaluated at "
                    "the data point coordinates does not have full column "
                    f"rank ({rank}/{nmonos})."
                    )
        raise LinAlgError(msg)

    return shift, scale, coeffs


def compute_interpolation(x, y, kernel, epsilon, powers, shift, scale, coeffs, xp):
    vec = _build_evaluation_coefficients(
        x, y, kernel, epsilon, powers, shift, scale, xp
    )
    return vec @ coeffs
