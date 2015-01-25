"""
A collection of functions to find the weights and abscissas for
Gaussian Quadrature.

These calculations are done by finding the eigenvalues of a
tridiagonal matrix whose entries are dependent on the coefficients
in the recursion formula for the orthogonal polynomials with the
corresponding weighting function over the interval.

Many recursion relations for orthogonal polynomials are given:

.. math::

    a1n f_{n+1} (x) = (a2n + a3n x ) f_n (x) - a4n f_{n-1} (x)

The recursion relation of interest is

.. math::

    P_{n+1} (x) = (x - A_n) P_n (x) - B_n P_{n-1} (x)

where :math:`P` has a different normalization than :math:`f`.

The coefficients can be found as:

.. math::

    A_n = -a2n / a3n
    \\qquad
    B_n = ( a4n / a3n \\sqrt{h_n-1 / h_n})^2

where

.. math::

    h_n = \\int_a^b w(x) f_n(x)^2

assume:

.. math::

    P_0 (x) = 1
    \\qquad
    P_{-1} (x) == 0

For the mathematical background, see [golub.welsch-1969-mathcomp]_ and
[abramowitz.stegun-1965]_.

Functions::

  gen_roots_and_weights  -- Generic roots and weights.
  j_roots                -- Jacobi
  js_roots               -- Shifted Jacobi
  la_roots               -- Generalized Laguerre
  h_roots                -- Hermite
  he_roots               -- Hermite (unit-variance)
  cg_roots               -- Ultraspherical (Gegenbauer)
  t_roots                -- Chebyshev of the first kind
  u_roots                -- Chebyshev of the second kind
  c_roots                -- Chebyshev of the first kind ([-2,2] interval)
  s_roots                -- Chebyshev of the second kind ([-2,2] interval)
  ts_roots               -- Shifted Chebyshev of the first kind.
  us_roots               -- Shifted Chebyshev of the second kind.
  p_roots                -- Legendre
  ps_roots               -- Shifted Legendre
  l_roots                -- Laguerre


.. [golub.welsch-1969-mathcomp]
   Golub, Gene H, and John H Welsch. 1969. Calculation of Gauss
   Quadrature Rules. *Mathematics of Computation* 23, 221-230+s1--s10.

.. [abramowitz.stegun-1965]
   Abramowitz, Milton, and Irene A Stegun. (1965) *Handbook of
   Mathematical Functions: with Formulas, Graphs, and Mathematical
   Tables*. Gaithersburg, MD: National Bureau of Standards.
   http://www.math.sfu.ca/~cbm/aands/

"""
#
# Author:  Travis Oliphant 2000
# Updated Sep. 2003 (fixed bugs --- tested to be accurate)

from __future__ import division, print_function, absolute_import

# Scipy imports.
import numpy as np
from numpy import all, any, exp, inf, pi, sqrt
from scipy import linalg

# Local imports.
from . import _ufuncs as cephes
_gam = cephes.gamma

__all__ = ['legendre', 'chebyt', 'chebyu', 'chebyc', 'chebys',
           'jacobi', 'laguerre', 'genlaguerre', 'hermite', 'hermitenorm',
           'gegenbauer', 'sh_legendre', 'sh_chebyt', 'sh_chebyu', 'sh_jacobi',
           'p_roots', 'ps_roots', 'j_roots', 'js_roots', 'l_roots', 'la_roots',
           'he_roots', 'ts_roots', 'us_roots', 's_roots', 't_roots', 'u_roots',
           'c_roots', 'cg_roots', 'h_roots',
           'eval_legendre', 'eval_chebyt', 'eval_chebyu', 'eval_chebyc',
           'eval_chebys', 'eval_jacobi', 'eval_laguerre', 'eval_genlaguerre',
           'eval_hermite', 'eval_hermitenorm', 'eval_gegenbauer',
           'eval_sh_legendre', 'eval_sh_chebyt', 'eval_sh_chebyu',
           'eval_sh_jacobi', 'poch', 'binom']


# For backward compatibility
poch = cephes.poch


class orthopoly1d(np.poly1d):

    def __init__(self, roots, weights=None, hn=1.0, kn=1.0, wfunc=None,
                 limits=None, monic=False, eval_func=None):
        np.poly1d.__init__(self, roots, r=1)
        equiv_weights = [weights[k] / wfunc(roots[k]) for
                         k in range(len(roots))]
        self.__dict__['weights'] = np.array(list(zip(roots,
                                                     weights, equiv_weights)))
        self.__dict__['weight_func'] = wfunc
        self.__dict__['limits'] = limits
        mu = sqrt(hn)
        if monic:
            evf = eval_func
            if evf:
                eval_func = lambda x: evf(x) / kn
            mu = mu / abs(kn)
            kn = 1.0
        self.__dict__['normcoef'] = mu
        self.__dict__['coeffs'] *= kn

        # Note: eval_func will be discarded on arithmetic
        self.__dict__['_eval_func'] = eval_func

    def __call__(self, v):
        if self._eval_func and not isinstance(v, np.poly1d):
            return self._eval_func(v)
        else:
            return np.poly1d.__call__(self, v)

    def _scale(self, p):
        if p == 1.0:
            return
        self.__dict__['coeffs'] *= p
        evf = self.__dict__['_eval_func']
        if evf:
            self.__dict__['_eval_func'] = lambda x: evf(x) * p
        self.__dict__['normcoef'] *= p


def _gen_roots_and_weights(n, mu0, an_func, bn_func, f, df, symmetrize, mu):
    """[x,w] = gen_roots_and_weights(n,an_func,sqrt_bn_func,mu)

    Returns the roots (x) of an nth order orthogonal polynomial,
    and weights (w) to use in appropriate Gaussian quadrature with that
    orthogonal polynomial.

    The polynomials have the recurrence relation
          P_n+1(x) = (x - A_n) P_n(x) - B_n P_n-1(x)

    an_func(n)          should return A_n
    sqrt_bn_func(n)     should return sqrt(B_n)
    mu ( = h_0 )        is the integral of the weight over the orthogonal
                        interval
    """
    k = np.arange(n, dtype='d')
    c = np.zeros((2, n))
    c[0,1:] = bn_func(k[1:])
    c[1,:] = an_func(k)
    x = linalg.eigvals_banded(c, overwrite_a_band=True)

    # improve roots by one application of Newton's method
    y = f(n, x)
    dy = df(n, x)
    x -= y/dy

    fm = f(n-1, x) 
    fm /= np.abs(fm).max()
    dy /= np.abs(dy).max()
    w = 1.0 / (fm * dy)

    if symmetrize:
        w = (w + w[::-1]) / 2
        x = (x - x[::-1]) / 2

    w *= mu0 / w.sum()

    if mu:
        return x, w, mu0
    else:
        return x, w

# Jacobi Polynomials 1               P^(alpha,beta)_n(x)


def j_roots(n, alpha, beta, mu=False):
    """Gauss-Jacobi quadrature

    Computes the sample points and weights for Gauss-Jacobi quadrature. The
    sample points are the roots of the `n`th degree Jacobi polynomial,
    :math:`P^{\\alpha, \\beta}_n(x)`.  These sample points and weights
    correctly integrate polynomials of degree :math:`2*n - 1` or less over the
    interval :math:`[-1, 1]` with weight function
    :math:`f(x) = (1 - x)^{\\alpha} (1 + x)^{\\beta}`.

    Parameters
    ----------
    n : int
        quadrature order 
    alpha : float
        alpha must be > -1
    beta : float
        beta must be > 0
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")
    if alpha <= -1 or beta <= -1:
        raise ValueError("alpha and beta must be greater than -1.")

    if alpha == 0.0 and beta == 0.0:
        return p_roots(m, mu)
    if alpha == beta:
        return cg_roots(m, alpha+0.5, mu)

    mu0 = 2.0**(alpha+beta+1)*cephes.beta(alpha+1, beta+1)
    a = alpha
    b = beta
    if a + b == 0.0:
        an_func = lambda k: np.where(k == 0, (b-a)/(2+a+b), 0.0)
    else:
        an_func = lambda k: np.where(k == 0, (b-a)/(2+a+b),
                  (b*b - a*a) / ((2.0*k+a+b)*(2.0*k+a+b+2)))

    bn_func = lambda k: 2.0 / (2.0*k+a+b)*np.sqrt((k+a)*(k+b) / (2*k+a+b+1)) \
              * np.where(k == 1, 1.0, np.sqrt(k*(k+a+b) / (2.0*k+a+b-1)))

    f = lambda n, x: cephes.eval_jacobi(n, a, b, x)
    df = lambda n, x: 0.5 * (n + a + b + 1) \
                      * cephes.eval_jacobi(n-1, a+1, b+1, x)
    return _gen_roots_and_weights(m, mu0, an_func, bn_func, f, df, False, mu)


def jacobi(n, alpha, beta, monic=False):
    """Returns the nth order Jacobi polynomial, P^(alpha,beta)_n(x)
    orthogonal over [-1,1] with weighting function
    (1-x)**alpha (1+x)**beta with alpha,beta > -1.
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    wfunc = lambda x: (1 - x)**alpha * (1 + x)**beta
    if n == 0:
        return orthopoly1d([], [], 1.0, 1.0, wfunc, (-1, 1), monic,
                           eval_func=np.ones_like)
    x, w, mu = j_roots(n, alpha, beta, mu=True)
    ab1 = alpha + beta + 1.0
    hn = 2**ab1 / (2 * n + ab1) * _gam(n + alpha + 1)
    hn *= _gam(n + beta + 1.0) / _gam(n + 1) / _gam(n + ab1)
    kn = _gam(2 * n + ab1) / 2.0**n / _gam(n + 1) / _gam(n + ab1)
    # here kn = coefficient on x^n term
    p = orthopoly1d(x, w, hn, kn, wfunc, (-1, 1), monic,
                    lambda x: eval_jacobi(n, alpha, beta, x))
    return p

# Jacobi Polynomials shifted         G_n(p,q,x)


def js_roots(n, p1, q1, mu=False):
    """Gauss-Jacobi (shifted) quadrature

    Computes the sample points and weights for Gauss-Jacobi (shifted)
    quadrature. The sample points are the roots of the `n`th degree shifted
    Jacobi polynomial, :math:`G^{p,q}_n(x)`.  These sample points and weights
    correctly integrate polynomials of degree :math:`2*n - 1` or less over the
    interval :math:`[0, 1]` with weight function
    :math:`f(x) = (1 - x)^{p-q} x^{q-1}`

    Parameters
    ----------
    n : int
        quadrature order 
    p1 : float
        (p1 - q1) must be > -1
    q1 : float
        q1 must be > 0
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights
        
    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    if (p1-q1) <= -1 or q1 <= 0:
        raise ValueError("(p - q) must be greater than -1, and q must be greater than 0.")
    xw = j_roots(n, p1-q1, q1-1, mu)
    return ((xw[0] + 1) / 2,) + xw[1:]


def sh_jacobi(n, p, q, monic=False):
    """Returns the nth order Jacobi polynomial, G_n(p,q,x)
    orthogonal over [0,1] with weighting function
    (1-x)**(p-q) (x)**(q-1) with p>q-1 and q > 0.
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    wfunc = lambda x: (1.0 - x)**(p - q) * (x)**(q - 1.)
    if n == 0:
        return orthopoly1d([], [], 1.0, 1.0, wfunc, (-1, 1), monic,
                           eval_func=np.ones_like)
    n1 = n
    x, w, mu0 = js_roots(n1, p, q, mu=True)
    hn = _gam(n + 1) * _gam(n + q) * _gam(n + p) * _gam(n + p - q + 1)
    hn /= (2 * n + p) * (_gam(2 * n + p)**2)
    # kn = 1.0 in standard form so monic is redundant.  Kept for compatibility.
    kn = 1.0
    pp = orthopoly1d(x, w, hn, kn, wfunc=wfunc, limits=(0, 1), monic=monic,
                     eval_func=lambda x: eval_sh_jacobi(n, p, q, x))
    return pp

# Generalized Laguerre               L^(alpha)_n(x)


def la_roots(n, alpha, mu=False):
    """Gauss-generalized Laguerre quadrature

    Computes the sample points and weights for Gauss-generalized Laguerre
    quadrature. The sample points are the roots of the `n`th degree generalized
    Laguerre polynomial, :math:`L^{\\alpha}_n(x)`.  These sample points and
    weights correctly integrate polynomials of degree :math:`2*n - 1` or less
    over the interval :math:`[0, inf]` with weight function
    :math:`f(x) = x^{\\alpha} e^{-x}`.

    Parameters
    ----------
    n : int
        quadrature order 
    alpha : float
        alpha must be > -1
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights
        
    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")
    if alpha < -1:
        raise ValueError("alpha must be greater than -1.")

    mu0 = cephes.gamma(alpha + 1)

    if m == 1:
        x = np.array([alpha+1.0], 'd')
        w = np.array([mu0], 'd')
        if mu:
            return x, w, mu0
        else:
            return x, w

    an_func = lambda k: 2 * k + alpha + 1
    bn_func = lambda k: -np.sqrt(k * (k + alpha))
    f = lambda n, x: cephes.eval_genlaguerre(n, alpha, x)
    df = lambda n, x: (n*cephes.eval_genlaguerre(n, alpha, x)
                     - (n + alpha)*cephes.eval_genlaguerre(n-1, alpha, x))/x
    return _gen_roots_and_weights(m, mu0, an_func, bn_func, f, df, False, mu)


def genlaguerre(n, alpha, monic=False):
    """Returns the nth order generalized (associated) Laguerre polynomial,
    L^(alpha)_n(x), orthogonal over [0,inf) with weighting function
    exp(-x) x**alpha with alpha > -1
    """
    if any(alpha <= -1):
        raise ValueError("alpha must be > -1")
    if n < 0:
        raise ValueError("n must be nonnegative.")

    if n == 0:
        n1 = n + 1
    else:
        n1 = n
    x, w, mu0 = la_roots(n1, alpha, mu=True)
    wfunc = lambda x: exp(-x) * x**alpha
    if n == 0:
        x, w = [], []
    hn = _gam(n + alpha + 1) / _gam(n + 1)
    kn = (-1)**n / _gam(n + 1)
    p = orthopoly1d(x, w, hn, kn, wfunc, (0, inf), monic,
                    lambda x: eval_genlaguerre(n, alpha, x))
    return p

# Laguerre                      L_n(x)


def l_roots(n, mu=False):
    """Gauss-Laguerre quadrature

    Computes the sample points and weights for Gauss-Laguerre quadrature.
    The sample points are the roots of the `n`th degree Laguerre polynomial,
    :math:`L_n(x)`.  These sample points and weights correctly integrate
    polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[0, inf]` with weight function :math:`f(x) = e^{-x}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    numpy.polynomial.laguerre.laggauss
    """
    return la_roots(n, 0.0, mu=mu)


def laguerre(n, monic=False):
    """Return the nth order Laguerre polynoimal, L_n(x), orthogonal over
    [0,inf) with weighting function exp(-x)
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    if n == 0:
        n1 = n + 1
    else:
        n1 = n
    x, w, mu0 = l_roots(n1, mu=True)
    if n == 0:
        x, w = [], []
    hn = 1.0
    kn = (-1)**n / _gam(n + 1)
    p = orthopoly1d(x, w, hn, kn, lambda x: exp(-x), (0, inf), monic,
                    lambda x: eval_laguerre(n, x))
    return p


# Hermite  1                         H_n(x)

def h_roots(n, mu=False):
    """Gauss-Hermite (physicst's) quadrature

    Computes the sample points and weights for Gauss-Hermite quadrature.
    The sample points are the roots of the `n`th degree Hermite polynomial,
    :math:`H_n(x)`.  These sample points and weights correctly integrate
    polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-inf, inf]` with weight function :math:`f(x) = e^{-x^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights
        
    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    numpy.polynomial.hermite.hermgauss
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    mu0 = np.sqrt(np.pi)
    an_func = lambda k: 0.0*k
    bn_func = lambda k: np.sqrt(k/2.0)
    f = cephes.eval_hermite
    df = lambda n, x: 2.0 * n * cephes.eval_hermite(n-1, x)
    return _gen_roots_and_weights(m, mu0, an_func, bn_func, f, df, True, mu)


def hermite(n, monic=False):
    """Return the nth order Hermite polynomial, H_n(x), orthogonal over
    (-inf,inf) with weighting function exp(-x**2)
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    if n == 0:
        n1 = n + 1
    else:
        n1 = n
    x, w, mu0 = h_roots(n1, mu=True)
    wfunc = lambda x: exp(-x * x)
    if n == 0:
        x, w = [], []
    hn = 2**n * _gam(n + 1) * sqrt(pi)
    kn = 2**n
    p = orthopoly1d(x, w, hn, kn, wfunc, (-inf, inf), monic,
                    lambda x: eval_hermite(n, x))
    return p

# Hermite  2                         He_n(x)


def he_roots(n, mu=False):
    """Gauss-Hermite (statistician's) quadrature

    Computes the sample points and weights for Gauss-Hermite quadrature.
    The sample points are the roots of the `n`th degree Hermite polynomial,
    :math:`He_n(x)`.  These sample points and weights correctly integrate
    polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-inf, inf]` with weight function :math:`f(x) = e^{-(x/2)^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights
        
    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    numpy.polynomial.hermite_e.hermegauss
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    mu0 = np.sqrt(np.pi/2.0)
    an_func = lambda k: 0.0*k
    bn_func = lambda k: np.sqrt(k)
    f = cephes.eval_hermitenorm
    df = lambda n, x: n * cephes.eval_hermitenorm(n-1, x)
    return _gen_roots_and_weights(m, mu0, an_func, bn_func, f, df, True, mu)


def hermitenorm(n, monic=False):
    """Return the nth order normalized Hermite polynomial, He_n(x), orthogonal
    over (-inf,inf) with weighting function exp(-(x/2)**2)
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    if n == 0:
        n1 = n + 1
    else:
        n1 = n
    x, w, mu0 = he_roots(n1, mu=True)
    wfunc = lambda x: exp(-x * x / 4.0)
    if n == 0:
        x, w = [], []
    hn = sqrt(2 * pi) * _gam(n + 1)
    kn = 1.0
    p = orthopoly1d(x, w, hn, kn, wfunc=wfunc, limits=(-inf, inf), monic=monic,
                    eval_func=lambda x: eval_hermitenorm(n, x))
    return p

# The remainder of the polynomials can be derived from the ones above.

# Ultraspherical (Gegenbauer)        C^(alpha)_n(x)


def cg_roots(n, alpha, mu=False):
    """Gauss-Gegenbauer quadrature

    Computes the sample points and weights for Gauss-Gegenbauer quadrature.
    The sample points are the roots of the `n`th degree Gegenbauer polynomial,
    :math:`C^{\\alpha}_n(x)`.  These sample points and weights correctly
    integrate polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-1, 1]` with weight function :math:`f(x) = (1-x^2)^{\\alpha-1/2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    alpha : float
        alpha must be > -0.5
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights
        
    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")
    if alpha < -0.5:
        raise ValueError("alpha must be greater than -0.5.")
    elif alpha == 0.0:
        # C(n,0,x) == 0 uniformly, however, as alpha->0, C(n,alpha,x)->T(n,x)
        # strictly, we should just error out here, since the roots are not
        # really defined, but we used to return something useful, so let's 
        # keep doing so.
        return t_roots(n, mu)

    mu0 = np.sqrt(np.pi) * cephes.gamma(alpha + 0.5) / cephes.gamma(alpha + 1)
    an_func = lambda k: 0.0 * k
    bn_func = lambda k: np.sqrt(k * (k + 2 * alpha - 1)
                        / (4 * (k + alpha) * (k + alpha - 1)))
    f = lambda n, x: cephes.eval_gegenbauer(n, alpha, x)
    df = lambda n, x: (-n*x*cephes.eval_gegenbauer(n, alpha, x)
         + (n + 2*alpha - 1)*cephes.eval_gegenbauer(n-1, alpha, x))/(1-x**2)
    return _gen_roots_and_weights(m, mu0, an_func, bn_func, f, df, True, mu)


def gegenbauer(n, alpha, monic=False):
    """Return the nth order Gegenbauer (ultraspherical) polynomial,
    C^(alpha)_n(x), orthogonal over [-1,1] with weighting function
    (1-x**2)**(alpha-1/2) with alpha > -1/2
    """
    base = jacobi(n, alpha - 0.5, alpha - 0.5, monic=monic)
    if monic:
        return base
    #  Abrahmowitz and Stegan 22.5.20
    factor = (_gam(2*alpha + n) * _gam(alpha + 0.5) /
              _gam(2*alpha) / _gam(alpha + 0.5 + n))
    base._scale(factor)
    base.__dict__['_eval_func'] = lambda x: eval_gegenbauer(float(n), alpha, x)
    return base

# Chebyshev of the first kind: T_n(x) =
#     n! sqrt(pi) / _gam(n+1./2)* P^(-1/2,-1/2)_n(x)
# Computed anew.


def t_roots(n, mu=False):
    """Gauss-Chebyshev (first kind) quadrature

    Computes the sample points and weights for Gauss-Chebyshev quadrature.
    The sample points are the roots of the `n`th degree Chebyshev polynomial of
    the first kind, :math:`T_n(x)`.  These sample points and weights correctly
    integrate polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-1, 1]` with weight function :math:`f(x) = 1/\sqrt{1 - x^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights
        
    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    numpy.polynomial.chebyshev.chebgauss
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError('n must be a positive integer.')
    x = np.cos(np.arange(2 * m - 1, 0, -2) * pi / (2 * m))
    w = np.empty_like(x)
    w.fill(pi/m)
    if mu:
        return x, w, pi
    else:
        return x, w


def chebyt(n, monic=False):
    """Return nth order Chebyshev polynomial of first kind, Tn(x).  Orthogonal
    over [-1,1] with weight function (1-x**2)**(-1/2).
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    wfunc = lambda x: 1.0 / sqrt(1 - x * x)
    if n == 0:
        return orthopoly1d([], [], pi, 1.0, wfunc, (-1, 1), monic,
                           lambda x: eval_chebyt(n, x))
    n1 = n
    x, w, mu = t_roots(n1, mu=True)
    hn = pi / 2
    kn = 2**(n - 1)
    p = orthopoly1d(x, w, hn, kn, wfunc, (-1, 1), monic,
                    lambda x: eval_chebyt(n, x))
    return p

# Chebyshev of the second kind
#    U_n(x) = (n+1)! sqrt(pi) / (2*_gam(n+3./2)) * P^(1/2,1/2)_n(x)


def u_roots(n, mu=False):
    """Gauss-Chebyshev (second kind) quadrature

    Computes the sample points and weights for Gauss-Chebyshev quadrature.
    The sample points are the roots of the `n`th degree Chebyshev polynomial of
    the second kind, :math:`U_n(x)`.  These sample points and weights correctly
    integrate polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-1, 1]` with weight function :math:`f(x) = \sqrt{1 - x^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points 
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError('n must be a positive integer.')
    t = np.arange(m, 0, -1) * pi / (m + 1)
    x = np.cos(t)
    w = pi * np.sin(t)**2 / (m + 1)
    if mu:
        return x, w, pi / 2
    else:
        return x, w


def chebyu(n, monic=False):
    """Return nth order Chebyshev polynomial of second kind, Un(x).  Orthogonal
    over [-1,1] with weight function (1-x**2)**(1/2).
    """
    base = jacobi(n, 0.5, 0.5, monic=monic)
    if monic:
        return base
    factor = sqrt(pi) / 2.0 * _gam(n + 2) / _gam(n + 1.5)
    base._scale(factor)
    return base

# Chebyshev of the first kind        C_n(x)


def c_roots(n, mu=False):
    """Gauss-Chebyshev (first kind) quadrature

    Computes the sample points and weights for Gauss-Chebyshev quadrature.
    The sample points are the roots of the `n`th degree Chebyshev polynomial of
    the first kind, :math:`C_n(x)`.  These sample points and weights correctly
    integrate polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-2, 2]` with weight function :math:`f(x) = 1/\sqrt{1 - (x/2)^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points 
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    xw = t_roots(n, mu)
    return (2 * xw[0],) + xw[1:]


def chebyc(n, monic=False):
    """Return nth order Chebyshev polynomial of first kind, Cn(x).  Orthogonal
    over [-2,2] with weight function (1-(x/2)**2)**(-1/2).
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    if n == 0:
        n1 = n + 1
    else:
        n1 = n
    x, w, mu0 = c_roots(n1, mu=True)
    if n == 0:
        x, w = [], []
    hn = 4 * pi * ((n == 0) + 1)
    kn = 1.0
    p = orthopoly1d(x, w, hn, kn,
                    wfunc=lambda x: 1.0 / sqrt(1 - x * x / 4.0),
                    limits=(-2, 2), monic=monic)
    if not monic:
        p._scale(2.0 / p(2))
        p.__dict__['_eval_func'] = lambda x: eval_chebyc(n, x)
    return p

# Chebyshev of the second kind       S_n(x)


def s_roots(n, mu=False):
    """Gauss-Chebyshev (second kind) quadrature

    Computes the sample points and weights for Gauss-Chebyshev quadrature.
    The sample points are the roots of the `n`th degree Chebyshev polynomial of
    the second kind, :math:`S_n(x)`.  These sample points and weights correctly
    integrate polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-2, 2]` with weight function :math:`f(x) = \sqrt{1 - (x/2)^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points 
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    xw = u_roots(n, mu)
    return (2 * xw[0],) + xw[1:]


def chebys(n, monic=False):
    """Return nth order Chebyshev polynomial of second kind, Sn(x).  Orthogonal
    over [-2,2] with weight function (1-(x/)**2)**(1/2).
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    if n == 0:
        n1 = n + 1
    else:
        n1 = n
    x, w, mu0 = s_roots(n1, mu=True)
    if n == 0:
        x, w = [], []
    hn = pi
    kn = 1.0
    p = orthopoly1d(x, w, hn, kn,
                    wfunc=lambda x: sqrt(1 - x * x / 4.0),
                    limits=(-2, 2), monic=monic)
    if not monic:
        factor = (n + 1.0) / p(2)
        p._scale(factor)
        p.__dict__['_eval_func'] = lambda x: eval_chebys(n, x)
    return p

# Shifted Chebyshev of the first kind     T^*_n(x)


def ts_roots(n, mu=False):
    """Gauss-Chebyshev (first kind, shifted) quadrature

    Computes the sample points and weights for Gauss-Chebyshev quadrature.
    The sample points are the roots of the `n`th degree shifted Chebyshev
    polynomial of the first kind, :math:`T_n(x)`.  These sample points and
    weights correctly integrate polynomials of degree :math:`2*n - 1` or less
    over the interval :math:`[0, 1]` with weight function
    :math:`f(x) = 1/\sqrt{x - x^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points 
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    xw = t_roots(n, mu)
    return ((xw[0] + 1) / 2,) + xw[1:]


def sh_chebyt(n, monic=False):
    """Return nth order shifted Chebyshev polynomial of first kind, Tn(x).
    Orthogonal over [0,1] with weight function (x-x**2)**(-1/2).
    """
    base = sh_jacobi(n, 0.0, 0.5, monic=monic)
    if monic:
        return base
    if n > 0:
        factor = 4**n / 2.0
    else:
        factor = 1.0
    base._scale(factor)
    return base


# Shifted Chebyshev of the second kind    U^*_n(x)
def us_roots(n, mu=False):
    """Gauss-Chebyshev (second kind, shifted) quadrature

    Computes the sample points and weights for Gauss-Chebyshev quadrature.
    The sample points are the roots of the `n`th degree shifted Chebyshev
    polynomial of the second kind, :math:`U_n(x)`.  These sample points and
    weights correctly integrate polynomials of degree :math:`2*n - 1` or less
    over the interval :math:`[0, 1]` with weight function
    :math:`f(x) = \sqrt{x - x^2}`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points 
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    xw = u_roots(n, mu)
    return ((xw[0] + 1) / 2,) + xw[1:]


def sh_chebyu(n, monic=False):
    """Return nth order shifted Chebyshev polynomial of second kind, Un(x).
    Orthogonal over [0,1] with weight function (x-x**2)**(1/2).
    """
    base = sh_jacobi(n, 2.0, 1.5, monic=monic)
    if monic:
        return base
    factor = 4**n
    base._scale(factor)
    return base

# Legendre


def p_roots(n, mu=False):
    """Gauss-Legendre quadrature

    Computes the sample points and weights for Gauss-Legendre quadrature.
    The sample points are the roots of the `n`th degree Legendre polynomial
    :math:`P_n(x)`.  These sample points and weights correctly integrate
    polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[-1, 1]` with weight function :math:`f(x) = 1.0`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights
        
    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    numpy.polynomial.legendre.leggauss
    """
    m = int(n)
    if n < 1 or n != m:
        raise ValueError("n must be a positive integer.")

    mu0 = 2.0
    an_func = lambda k: 0.0 * k
    bn_func = lambda k: k * np.sqrt(1.0 / (4 * k * k - 1))
    f = cephes.eval_legendre
    df = lambda n, x: (-n*x*cephes.eval_legendre(n, x)
                     + n*cephes.eval_legendre(n-1, x))/(1-x**2)
    return _gen_roots_and_weights(m, mu0, an_func, bn_func, f, df, True, mu)


def legendre(n, monic=False):
    """
    Legendre polynomial coefficients

    Returns the nth-order Legendre polynomial, P_n(x), orthogonal over
    [-1, 1] with weight function 1.

    Parameters
    ----------
    n
        Order of the polynomial
    monic : bool, optional
        If True, output is a monic polynomial (normalized so the leading
        coefficient is 1).  Default is False.

    Returns
    -------
    P : orthopoly1d
        The Legendre polynomial object

    Examples
    --------
    Generate the 3rd-order Legendre polynomial 1/2*(5x^3 + 0x^2 - 3x + 0):

    >>> legendre(3)
    poly1d([ 2.5,  0. , -1.5, -0. ])

    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    if n == 0:
        n1 = n + 1
    else:
        n1 = n
    x, w, mu0 = p_roots(n1, mu=True)
    if n == 0:
        x, w = [], []
    hn = 2.0 / (2 * n + 1)
    kn = _gam(2 * n + 1) / _gam(n + 1)**2 / 2.0**n
    p = orthopoly1d(x, w, hn, kn, wfunc=lambda x: 1.0, limits=(-1, 1),
                    monic=monic, eval_func=lambda x: eval_legendre(n, x))
    return p

# Shifted Legendre              P^*_n(x)


def ps_roots(n, mu=False):
    """Gauss-Legendre (shifted) quadrature

    Computes the sample points and weights for Gauss-Legendre quadrature.
    The sample points are the roots of the `n`th degree shifted Legendre
    polynomial :math:`P^*_n(x)`.  These sample points and weights correctly
    integrate polynomials of degree :math:`2*n - 1` or less over the interval
    :math:`[0, 1]` with weight function :math:`f(x) = 1.0`.

    Parameters
    ----------
    n : int
        quadrature order 
    mu : boolean
        If True, return the sum of the weights, optional.

    Returns
    -------
    x : ndarray
        Sample points
    w : ndarray
        Weights
    mu : float
        Sum of the weights

    See Also
    --------
    integrate.quadrature
    integrate.fixed_quad
    """
    xw = p_roots(n, mu)
    return ((xw[0] + 1) / 2,) + xw[1:]


def sh_legendre(n, monic=False):
    """Returns the nth order shifted Legendre polynomial, P^*_n(x), orthogonal
    over [0,1] with weighting function 1.
    """
    if n < 0:
        raise ValueError("n must be nonnegative.")

    wfunc = lambda x: 0.0 * x + 1.0
    if n == 0:
        return orthopoly1d([], [], 1.0, 1.0, wfunc, (0, 1), monic,
                           lambda x: eval_sh_legendre(n, x))
    x, w, mu0 = ps_roots(n, mu=True)
    hn = 1.0 / (2 * n + 1.0)
    kn = _gam(2 * n + 1) / _gam(n + 1)**2
    p = orthopoly1d(x, w, hn, kn, wfunc, limits=(0, 1), monic=monic,
                    eval_func=lambda x: eval_sh_legendre(n, x))
    return p

# -----------------------------------------------------------------------------
# Vectorized functions for evaluation
# -----------------------------------------------------------------------------
from ._ufuncs import (binom, eval_jacobi, eval_sh_jacobi, eval_gegenbauer,
                      eval_chebyt, eval_chebyu, eval_chebys, eval_chebyc,
                      eval_sh_chebyt, eval_sh_chebyu, eval_legendre,
                      eval_sh_legendre, eval_genlaguerre, eval_laguerre,
                      eval_hermite, eval_hermitenorm)
