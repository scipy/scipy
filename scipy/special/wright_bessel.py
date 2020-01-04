import numbers
import numpy as np
import warnings

from ._ufuncs import rgamma
from scipy._lib._util import _asarray_validated


def wright_bessel(a, b, z, max_iter=1000, tol=1e-16):
    r"""Compute Wright's generalized Bessel function.

    Wright's generalized Bessel function is an entire function and defined as

    .. math:: \Phi(a, b; z) = \sum_{k=0}^\infty \frac{z^k}{k! \Gamma(a k + b)}

    Parameters
    ----------
    a : float
        a > -1
    b : complex
    z : array_like of complex
    max_iter : int, optional
        Maximum number of terms included in the series expansion.
    tol : float, optional
        Evaluation tolerance, tol > 0.

    References
    ----------
    .. [1] https://dlmf.nist.gov/10.46.E1
    """
    if tol <= 0:
        raise ValueError("The tolerance tol must be larger than zero.")

    max_iter = int(max_iter)
    if max_iter <= 0:
        raise ValueError("The number if terms included, max_iter, must be "
                         "larger a equal to 1.")

    a = float(a)
    if a <= -1:
        raise ValueError("Argument a must be larger than -1.")

    if not isinstance(b, numbers.Number):
        raise ValueError("Argument b must be a number.")
    elif not np.isfinite(b):
        raise ValueError("Argument b must be a finite number.")

    b = 1. * b  # make it at least a float

    z = _asarray_validated(z, check_finite=False)

    # determine dtype of result
    if np.iscomplexobj(z) or isinstance(b, complex):
        dtype = max(z.dtype, complex)
        result = np.zeros_like(z, dtype=dtype)
    else:
        dtype = max(z.dtype, float)
        result = np.zeros_like(z, dtype=dtype)

    # deal with np.nan and np.inf
    # TODO: np.isfinite(z) otherwise return np.nan

    # first some special cases
    if np.all(z == 0):
        return 0
    elif a == 0:
        result += np.exp(z) * rgamma(b)
    # now the general case
    else:
        _wright_series(a, b, z, result, eps_rel=tol, eps_abs=tol)

    return result


def _wright_series(a, b, z, result, max_iter=1000,
                   eps_rel=1e-17, eps_abs=1e-17):
    """Compute Wright's generalized Bessel function by series."""

    # maximum of abs(z)
    imax = np.argmax(np.abs(z))
    # zmax = z.flatten()[imax]

    # term k=0
    result += rgamma(b)
    z_k = np.ones_like(z, dtype=max(float, z.dtype))
    converged = False
    for k in range(1, max_iter):
        z_k *= z/k
        rg = rgamma(a*k + b)
        term = z_k * rg
        result += term

        # check convergence, only if rgamma != 0 (a*k+b = 0, -1, -2, ..)
        if rg != 0:
            term_max = np.abs(term.flatten()[imax])
            res_max = np.abs(result.flatten()[imax])
            if (term_max < eps_abs) and (term_max < eps_rel * res_max):
                converged = True
                break

    if not converged:
        warnings.warn("Wright's generalized Bessel function did not converge "
                      "to desired accuracy", RuntimeWarning)
