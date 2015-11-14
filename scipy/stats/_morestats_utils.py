"""
Utility functions related to statistics.

These have been moved into a separate module to break import cycles.

"""
from __future__ import division, print_function, absolute_import

import numpy as np

__all__ = ['filliben']


def filliben(n):
    """
    Approximations of uniform order statistic medians.

    Parameters
    ----------
    n : int
        Sample size.

    Returns
    -------
    v : 1d float array
        Approximations of the order statistic medians.

    References
    ----------
    .. [1] James J. Filliben, "The Probability Plot Correlation Coefficient
           Test for Normality", Technometrics, Vol. 17, pp. 111-117, 1975.

    Examples
    --------
    Order statistics of the uniform distribution on the unit interval
    are marginally distributed according to beta distributions.
    The expectations of these order statistic are evenly spaced across
    the interval, but the distributions are skewed in a way that
    pushes the medians slightly towards the endpoints of the unit interval:

    >>> n = 4
    >>> k = np.arange(1, n+1)
    >>> from scipy.stats import beta
    >>> a = k
    >>> b = n-k+1
    >>> beta.mean(a, b)
    array([ 0.2,  0.4,  0.6,  0.8])
    >>> beta.median(a, b)
    array([ 0.15910358,  0.38572757,  0.61427243,  0.84089642])

    The Filliben approximation uses the exact medians of the smallest
    and greatest order statistics, and the remaining medians are approximated
    by points spread evenly across a sub-interval of the unit interval:

    >>> from scipy.stats import filliben
    >>> filliben(n)
    array([ 0.15910358,  0.38545246,  0.61454754,  0.84089642])

    """
    v = np.zeros(n, dtype=np.float64)
    v[-1] = 0.5**(1.0 / n)
    v[0] = 1 - v[-1]
    i = np.arange(2, n)
    v[1:-1] = (i - 0.3175) / (n + 0.365)
    return v
