"""
Implementations of fitting algorithms for common functions that can be
computed analytically.

Functions
---------
.. autosummary::
   :toctree: generated/

    exp_fit
    pow_fit
"""

from __future__ import division, print_function, absolute_import

__all__ = ['exp_fit', 'pow_fit']

import numpy as np
from numpy import (
    argsort, asfarray, cumsum, diff, empty, empty_like, exp, log,
    square,
)
from numpy.linalg import inv


def exp_fit(x, y, sorted=True):
    """
    Fit an exponential curve to raveled 1D data.

    This algorithm does not require any a-priori knowledge of the data,
    such as the intercept. The fitting parameters are comptued for:

    .. math::

       y = A + Be^{Cx}

    Parameters
    ----------
    x : array-like
        The x-values of the data points. The fit will be performed on a
        raveled version of this array.
    y : array-like
        The y-values of the data points corresponding to `x`. Must be
        the same size as `x`. The fit will be performed on a raveled
        version of this array.
    sorted : bool
        Set to True if `x` is already monotonically increasing or
        decreasing. If False, x will be sorted into increasing order,
        and y will be sorted along with it.

    Return
    ------
    a, b, c : array
        A 3-element array of optimized fitting parameters. The first
        element is the additive bias, the second the multiplicative, and
        the third the exponential.

    Notes
    -----
    The fit is computed non-iteratively in a single pass. It can be used
    to initialize other regression methods based on different
    optimization criteria. The algorithm and the theory behind it is
    presented in the paper below.

    References
    ----------
    Jacquelin, Jean. “REGRESSIONS Et EQUATIONS INTEGRALES” 14 Jan. 2009, pp. 15–18., https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    x = asfarray(x).ravel()
    y = asfarray(y).ravel()
    if x.size != y.size:
        raise ValueError('x and y must be the same size')
    if not sorted:
        # Is there a better way to do this in scipy?
        ind = argsort(x)
        x = x[ind]
        y = y[ind]

    s = empty_like(y)
    s[0] = 0
    s[1:] = cumsum(0.5 * (y[1:] + y[:-1]) * diff(x))
    # This might be better: Needs a benchmark
    #s[1:] = y[:-1]
    #s[1:] += y[1:]
    #s[1:] *= np.diff(x)
    #s *= 0.5
    #s = np.cumsum(s)

    xn = x - x[0]
    yn = y - y[0]

    sx2 = square(xn).sum()
    sxs = (xn * s).sum()
    sys = (yn * s).sum()
    ss2 = square(s).sum()
    sxy = (xn * yn).sum()

    out = empty(3, dtype=float)

    _, out[2] = inv([[sx2, sxs], [sxs, ss2]]).dot([[sxy], [sys]])

    ex = exp(out[2] * x)

    se1 = ex.sum()
    se2 = square(ex).sum()
    sy0 = y.sum()
    sye = (y * ex).sum()

    out[0], out[1] = inv([[x.size, se1], [se1, se2]]).dot([[sy0], [sye]])

    return out


def pow_fit(x, y, sorted=True):
    """
    Fit a power curve to raveled 1D data.

    The fitting parameters are comptued for:

    .. math::

       y = A + Bx^C

    Parameters
    ----------
    x : array-like
        The x-values of the data points. The fit will be performed on a
        raveled version of this array. All elements must be positive.
    y : array-like
        The y-values of the data points corresponding to `x`. Must be
        the same size as `x`. The fit will be performed on a raveled
        version of this array.
    sorted : bool
        Set to True if `x` is already monotonically increasing or
        decreasing. If False, x will be sorted into increasing order,
        and y will be sorted along with it.

    Return
    ------
    a, b, c : array
        A 3-element array of optimized fitting parameters. The first
        element is the additive bias, the second the multiplicative, and
        the third is the power.

    Notes
    -----
    The fit is computed non-iteratively in a single pass. It can be used
    to initialize other regression methods based on different
    optimization criteria. ``pow_fit(x, y, sorted)`` is equivalent to
    ``exp_fit(log(x), y, sorted)`` since

    .. math::

       A + Be^{Cx} = A + B(e^x)^C

    The algorithm and the theory behind it is presented in the paper
    below.

    References
    ----------
    Jacquelin, Jean. “REGRESSIONS Et EQUATIONS INTEGRALES” 14 Jan. 2009, pp. 15–18., https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
    """
    return exp_fit(log(x), y, sorted)
