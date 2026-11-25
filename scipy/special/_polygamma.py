"""This was in ._basic.py, but was moved here to avoid a circular import.

_basic.py consumes _support_alternative_backends.py and has a number of
functions which are array-agnostic modulo delegation for internal dependencies.
"""
import numpy as np

import scipy.special._ufuncs as _ufuncs


def polygamma(n, x):
    r"""Polygamma functions.

    Defined as :math:`\psi^{(n)}(x)` where :math:`\psi` is the
    `digamma` function. See [dlmf]_ for details.

    Parameters
    ----------
    n : array_like
        The order of the derivative of the digamma function; must be
        integral
    x : array_like
        Real valued input

    Returns
    -------
    ndarray
        Function results

    See Also
    --------
    digamma

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/5.15

    Examples
    --------
    >>> from scipy import special
    >>> x = [2, 3, 25.5]
    >>> special.polygamma(1, x)
    array([ 0.64493407,  0.39493407,  0.03999467])
    >>> special.polygamma(0, x) == special.psi(x)
    array([ True,  True,  True], dtype=bool)

    """
    n, x = np.asarray(n), np.asarray(x)
    fac2 = (-1.0)**(n+1) * _ufuncs.gamma(n+1.0) * _ufuncs._zeta(n+1, x)
    return np.where(n == 0, _ufuncs.psi(x), fac2)
