"""
Evaluate orthogonal polynomial values using special functions or
recurrence relations.

References
----------

.. [AMS55] Abramowitz & Stegun, Section 22.5.

.. [MH] Mason & Handscombe, Chebyshev Polynomials, CRC Press (2003).

"""
#
# Copyright (C) 2009 Pauli Virtanen
# Distributed under the same license as Scipy.
#

#------------------------------------------------------------------------------
# Actual evaluation functions
#------------------------------------------------------------------------------

import numpy as np
from scipy.special._ufuncs import gamma, hyp2f1, hyp1f1, gammaln, _eval_chebyt
from numpy import exp, sqrt

def binom(n, k):
    """
    binom(n, k)

    Binomial coefficient
    """
    return exp(gammaln(1+n) - gammaln(1+k) - gammaln(1+n-k))

def eval_jacobi(n, alpha, beta, x, out=None):
    """
    eval_jacobi(n, alpha, beta, x, out=None)

    Evaluate Jacobi polynomial at a point.
    """
    d = binom(n+alpha, n)
    a = -n
    b = n + alpha + beta + 1
    c = alpha + 1
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_sh_jacobi(n, p, q, x, out=None):
    """
    eval_sh_jacobi(n, p, q, x, out=None)

    Evaluate shifted Jacobi polynomial at a point.
    """
    factor = exp(gammaln(1+n) + gammaln(n+p) - gammaln(2*n+p))
    return factor * eval_jacobi(n, p-q, q-1, 2*x-1)

def eval_gegenbauer(n, alpha, x, out=None):
    """
    eval_gegenbauer(n, alpha, x, out=None)

    Evaluate Gegenbauer polynomial at a point.
    """
    d = gamma(n+2*alpha)/gamma(1+n)/gamma(2*alpha)
    a = -n
    b = n + 2*alpha
    c = alpha + 0.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_chebyt(n, x, out=None):
    """
    eval_chebyt(n, x, out=None)

    Evaluate Chebyshev T polynomial at a point.

    This routine is numerically stable for `x` in ``[-1, 1]`` at least
    up to order ``10000``.
    """
    return _eval_chebyt(n, x, out)

def eval_chebyu(n, x, out=None):
    """
    eval_chebyu(n, x, out=None)

    Evaluate Chebyshev U polynomial at a point.
    """
    d = n+1
    a = -n
    b = n+2
    c = 1.5
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_chebys(n, x, out=None):
    """
    eval_chebys(n, x, out=None)

    Evaluate Chebyshev S polynomial at a point.
    """
    return eval_chebyu(n, x/2, out=out)

def eval_chebyc(n, x, out=None):
    """
    eval_chebyc(n, x, out=None)

    Evaluate Chebyshev C polynomial at a point.
    """
    return 2*eval_chebyt(n, x/2.0, out)

def eval_sh_chebyt(n, x, out=None):
    """
    eval_sh_chebyt(n, x, out=None)

    Evaluate shifted Chebyshev T polynomial at a point.
    """
    return eval_chebyt(n, 2*x-1, out=out)

def eval_sh_chebyu(n, x, out=None):
    """
    eval_sh_chebyu(n, x, out=None)

    Evaluate shifted Chebyshev U polynomial at a point.
    """
    return eval_chebyu(n, 2*x-1, out=out)

def eval_legendre(n, x, out=None):
    """
    eval_legendre(n, x, out=None)

    Evaluate Legendre polynomial at a point.
    """
    d = 1
    a = -n
    b = n+1
    c = 1
    g = (1-x)/2.0
    return hyp2f1(a, b, c, g) * d

def eval_sh_legendre(n, x, out=None):
    """
    eval_sh_legendre(n, x, out=None)

    Evaluate shifted Legendre polynomial at a point.
    """
    return eval_legendre(n, 2*x-1, out=out)

def eval_genlaguerre(n, alpha, x, out=None):
    """
    eval_genlaguerre(n, alpha, x, out=None)

    Evaluate generalized Laguerre polynomial at a point.
    """
    d = binom(n+alpha, n)
    a = -n
    b = alpha + 1
    g = x
    return hyp1f1(a, b, g) * d

def eval_laguerre(n, x, out=None):
    """
    eval_laguerre(n, x, out=None)

    Evaluate Laguerre polynomial at a point.
    """
    return eval_genlaguerre(n, 0., x, out=out)

def eval_hermite(n, x, out=None):
    """
    eval_hermite(n, x, out=None)

    Evaluate Hermite polynomial at a point.
    """
    n, x = np.broadcast_arrays(n, x)
    n, x = np.atleast_1d(n, x)

    if out is None:
        out = np.zeros_like(0*n + 0*x)
    if (n % 1 != 0).any():
        raise ValueError("Order must be integer")

    even = (n % 2 == 0)

    m = n[even]/2
    out[even] = ((-1)**m * 2**(2*m) * gamma(1+m)
                 * eval_genlaguerre(m, -0.5, x[even]**2))

    m = (n[~even]-1)/2
    out[~even] = ((-1)**m * 2**(2*m+1) * gamma(1+m)
                  * x[~even] * eval_genlaguerre(m, 0.5, x[~even]**2))

    return out

def eval_hermitenorm(n, x, out=None):
    """
    eval_hermitenorm(n, x, out=None)

    Evaluate normalized Hermite polynomial at a point.
    """
    return eval_hermite(n, x/sqrt(2)) * 2**(-n/2.0)
