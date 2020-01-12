from __future__ import division, print_function, absolute_import

from ._ufuncs import _wright_bessel


def wright_bessel(a, b, x):
    r"""Compute Wright's generalized Bessel function.

    Wright's generalized Bessel function is an entire function and defined as

    .. math:: \Phi(a, b; x) = \sum_{k=0}^\infty \frac{x^k}{k! \Gamma(a k + b)}

    See also [1].

    Parameters
    ----------
    a : float
        a >= 0
    b : float
        b >= 0
    x : array_like of float
        x >= 0

    References
    ----------
    .. [1] https://dlmf.nist.gov/10.46.E1
    """
    return _wright_bessel(a, b, x)
