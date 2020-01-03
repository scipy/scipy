import numbers
import numpy as np
import warnings

from ._ufuncs import iv, rgamma
from scipy._lib._util import _asarray_validated


def wright_bessel(a, b, z, max_iter=1000, tol=1e-15):
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

    b = 1. * b  # make it at least a float

    z = _asarray_validated(z, check_finite=True)

    # first some special cases
    if np.all(z == 0):
        return 0
    elif a == 0:
        return np.exp(z) * rgamma(b)
    elif a == 1:
        # Phi(1, v+1, z**2/4) = (z/2)**(-v) * I_v(z)
        x = 2 * np.sqrt(z)
        v = b - 1
        return np.power(0.5*x, -v) * iv(v, x)
    # now the general case
    else:
        if np.iscomplexobj(z) or isinstance(b, complex):
            dtype = max(z.dtype, complex)
            result = np.zeros_like(z, dtype=dtype)
        else:
            dtype = max(z.dtype, float)
            result = np.zeros_like(z, dtype=dtype)

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
        term = z_k * rgamma(a*k + b)
        result += term

        # stop early?
        term_max = np.abs(term.flatten()[imax])
        res_max = np.abs(result.flatten()[imax])
        if (term_max < eps_abs) and (term_max < eps_rel * res_max):
            converged = True
            break

    if not converged:
        warnings.warn("Wright's generalized Bessel function did not converge "
                      "to desired accuracy", RuntimeWarning)
