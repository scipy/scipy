#
# Author:  Travis Oliphant, 2002
#

from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.lib.six import xrange
from numpy import (pi, asarray, floor, isscalar, iscomplex, real, imag, sqrt,
                   where, mgrid, cos, sin, exp, place, issubdtype, extract,
                   less, vectorize, inexact, nan, zeros, atleast_1d, sinc)
from ._ufuncs import (ellipkm1, mathieu_a, mathieu_b, iv, jv, gamma, psi, zeta,
                      hankel1, hankel2, yv, kv, gammaln, ndtri, errprint, poch,
                      binom, xlogy)
from . import _ufuncs
import types
from . import specfun
from . import orthogonal
import warnings

__all__ = ['agm', 'ai_zeros', 'assoc_laguerre', 'bei_zeros', 'beip_zeros',
           'ber_zeros', 'bernoulli', 'berp_zeros', 'bessel_diff_formula',
           'bi_zeros', 'clpmn', 'clqmn', 'comb', 'digamma', 'diric', 'ellipk',
           'erf_zeros', 'erfcinv', 'erfinv', 'errprint', 'euler', 'factorial',
           'factorialk', 'factorial2', 'fresnel_zeros',
           'fresnelc_zeros', 'fresnels_zeros', 'gamma', 'gammaln', 'h1vp',
           'h2vp', 'hankel1', 'hankel2', 'hyp0f1', 'iv', 'ivp', 'jn_zeros',
           'jnjnp_zeros', 'jnp_zeros', 'jnyn_zeros', 'jv', 'jvp', 'kei_zeros',
           'keip_zeros', 'kelvin_zeros', 'ker_zeros', 'kerp_zeros', 'kv',
           'kvp', 'lmbda', 'lpmn', 'lpn', 'lqmn', 'lqn', 'mathieu_a',
           'mathieu_b', 'mathieu_even_coef', 'mathieu_odd_coef', 'ndtri',
           'obl_cv_seq', 'pbdn_seq', 'pbdv_seq', 'pbvv_seq', 'perm',
           'polygamma', 'pro_cv_seq', 'psi', 'riccati_jn', 'riccati_yn',
           'sinc', 'sph_in', 'sph_inkn',
           'sph_jn', 'sph_jnyn', 'sph_kn', 'sph_yn', 'y0_zeros', 'y1_zeros',
           'y1p_zeros', 'yn_zeros', 'ynp_zeros', 'yv', 'yvp', 'zeta',
           'entr', 'rel_entr', 'kl_div',
           'SpecialFunctionWarning']


class SpecialFunctionWarning(Warning):
    """Warning that can be issued with ``errprint(True)``"""
    pass
warnings.simplefilter("always", category=SpecialFunctionWarning)


def diric(x,n):
    """Returns the periodic sinc function, also called the Dirichlet function

    The Dirichlet function is defined as::

        diric(x) = sin(x * n/2) / (n * sin(x / 2)),

    where n is a positive integer.

    Parameters
    ----------
    x : array_like
        Input data
    n : int
        Integer defining the periodicity.

    Returns
    -------
    diric : ndarray

    """
    x,n = asarray(x), asarray(n)
    n = asarray(n + (x-x))
    x = asarray(x + (n-n))
    if issubdtype(x.dtype, inexact):
        ytype = x.dtype
    else:
        ytype = float
    y = zeros(x.shape,ytype)

    mask1 = (n <= 0) | (n != floor(n))
    place(y,mask1,nan)

    z = asarray(x / 2.0 / pi)
    mask2 = (1-mask1) & (z == floor(z))
    zsub = extract(mask2,z)
    nsub = extract(mask2,n)
    place(y,mask2,pow(-1,zsub*(nsub-1)))

    mask = (1-mask1) & (1-mask2)
    xsub = extract(mask,x)
    nsub = extract(mask,n)
    place(y,mask,sin(nsub*xsub/2.0)/(nsub*sin(xsub/2.0)))
    return y


def jnjnp_zeros(nt):
    """Compute nt (<=1200) zeros of the Bessel functions Jn and Jn'
    and arange them in order of their magnitudes.

    Returns
    -------
    zo[l-1] : ndarray
        Value of the lth zero of Jn(x) and Jn'(x). Of length `nt`.
    n[l-1] : ndarray
        Order of the Jn(x) or Jn'(x) associated with lth zero. Of length `nt`.
    m[l-1] : ndarray
        Serial number of the zeros of Jn(x) or Jn'(x) associated
        with lth zero. Of length `nt`.
    t[l-1] : ndarray
        0 if lth zero in zo is zero of Jn(x), 1 if it is a zero of Jn'(x). Of
        length `nt`.

    See Also
    --------
    jn_zeros, jnp_zeros : to get separated arrays of zeros.
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt > 1200):
        raise ValueError("Number must be integer <= 1200.")
    nt = int(nt)
    n,m,t,zo = specfun.jdzo(nt)
    return zo[1:nt+1],n[:nt],m[:nt],t[:nt]


def jnyn_zeros(n,nt):
    """Compute nt zeros of the Bessel functions Jn(x), Jn'(x), Yn(x), and
    Yn'(x), respectively. Returns 4 arrays of length nt.

    See jn_zeros, jnp_zeros, yn_zeros, ynp_zeros to get separate arrays.
    """
    if not (isscalar(nt) and isscalar(n)):
        raise ValueError("Arguments must be scalars.")
    if (floor(n) != n) or (floor(nt) != nt):
        raise ValueError("Arguments must be integers.")
    if (nt <= 0):
        raise ValueError("nt > 0")
    return specfun.jyzo(abs(n),nt)


def jn_zeros(n,nt):
    """Compute nt zeros of the Bessel function Jn(x).
    """
    return jnyn_zeros(n,nt)[0]


def jnp_zeros(n,nt):
    """Compute nt zeros of the Bessel function Jn'(x).
    """
    return jnyn_zeros(n,nt)[1]


def yn_zeros(n,nt):
    """Compute nt zeros of the Bessel function Yn(x).
    """
    return jnyn_zeros(n,nt)[2]


def ynp_zeros(n,nt):
    """Compute nt zeros of the Bessel function Yn'(x).
    """
    return jnyn_zeros(n,nt)[3]


def y0_zeros(nt,complex=0):
    """Returns nt (complex or real) zeros of Y0(z), z0, and the value
    of Y0'(z0) = -Y1(z0) at each zero.
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("Arguments must be scalar positive integer.")
    kf = 0
    kc = (complex != 1)
    return specfun.cyzo(nt,kf,kc)


def y1_zeros(nt,complex=0):
    """Returns nt (complex or real) zeros of Y1(z), z1, and the value
    of Y1'(z1) = Y0(z1) at each zero.
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("Arguments must be scalar positive integer.")
    kf = 1
    kc = (complex != 1)
    return specfun.cyzo(nt,kf,kc)


def y1p_zeros(nt,complex=0):
    """Returns nt (complex or real) zeros of Y1'(z), z1', and the value
    of Y1(z1') at each zero.
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("Arguments must be scalar positive integer.")
    kf = 2
    kc = (complex != 1)
    return specfun.cyzo(nt,kf,kc)


def _bessel_diff_formula(v, z, n, L, phase):
    # from AMS55.
    # L(v,z) = J(v,z), Y(v,z), H1(v,z), H2(v,z), phase = -1
    # L(v,z) = I(v,z) or exp(v*pi*i)K(v,z), phase = 1
    # For K, you can pull out the exp((v-k)*pi*i) into the caller
    p = 1.0
    s = L(v-n, z)
    for i in xrange(1, n+1):
        p = phase * (p * (n-i+1)) / i   # = choose(k, i)
        s += p*L(v-n + i*2, z)
    return s / (2.**n)


bessel_diff_formula = np.deprecate(_bessel_diff_formula,
    message="bessel_diff_formula is a private function, do not use it!")


def jvp(v,z,n=1):
    """Return the nth derivative of Jv(z) with respect to z.
    """
    if not isinstance(n,int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return jv(v,z)
    else:
        return _bessel_diff_formula(v, z, n, jv, -1)
#        return (jvp(v-1,z,n-1) - jvp(v+1,z,n-1))/2.0


def yvp(v,z,n=1):
    """Return the nth derivative of Yv(z) with respect to z.
    """
    if not isinstance(n,int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return yv(v,z)
    else:
        return _bessel_diff_formula(v, z, n, yv, -1)
#        return (yvp(v-1,z,n-1) - yvp(v+1,z,n-1))/2.0


def kvp(v,z,n=1):
    """Return the nth derivative of Kv(z) with respect to z.
    """
    if not isinstance(n,int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return kv(v,z)
    else:
        return (-1)**n * _bessel_diff_formula(v, z, n, kv, 1)


def ivp(v,z,n=1):
    """Return the nth derivative of Iv(z) with respect to z.
    """
    if not isinstance(n,int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return iv(v,z)
    else:
        return _bessel_diff_formula(v, z, n, iv, 1)


def h1vp(v,z,n=1):
    """Return the nth derivative of H1v(z) with respect to z.
    """
    if not isinstance(n,int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return hankel1(v,z)
    else:
        return _bessel_diff_formula(v, z, n, hankel1, -1)
#        return (h1vp(v-1,z,n-1) - h1vp(v+1,z,n-1))/2.0


def h2vp(v,z,n=1):
    """Return the nth derivative of H2v(z) with respect to z.
    """
    if not isinstance(n,int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return hankel2(v,z)
    else:
        return _bessel_diff_formula(v, z, n, hankel2, -1)
#        return (h2vp(v-1,z,n-1) - h2vp(v+1,z,n-1))/2.0


def sph_jn(n,z):
    """Compute the spherical Bessel function jn(z) and its derivative for
    all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z):
        nm,jn,jnp,yn,ynp = specfun.csphjy(n1,z)
    else:
        nm,jn,jnp = specfun.sphj(n1,z)
    return jn[:(n+1)], jnp[:(n+1)]


def sph_yn(n,z):
    """Compute the spherical Bessel function yn(z) and its derivative for
    all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z,0):
        nm,jn,jnp,yn,ynp = specfun.csphjy(n1,z)
    else:
        nm,yn,ynp = specfun.sphy(n1,z)
    return yn[:(n+1)], ynp[:(n+1)]


def sph_jnyn(n,z):
    """Compute the spherical Bessel functions, jn(z) and yn(z) and their
    derivatives for all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z,0):
        nm,jn,jnp,yn,ynp = specfun.csphjy(n1,z)
    else:
        nm,yn,ynp = specfun.sphy(n1,z)
        nm,jn,jnp = specfun.sphj(n1,z)
    return jn[:(n+1)],jnp[:(n+1)],yn[:(n+1)],ynp[:(n+1)]


def sph_in(n,z):
    """Compute the spherical Bessel function in(z) and its derivative for
    all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z):
        nm,In,Inp,kn,knp = specfun.csphik(n1,z)
    else:
        nm,In,Inp = specfun.sphi(n1,z)
    return In[:(n+1)], Inp[:(n+1)]


def sph_kn(n,z):
    """Compute the spherical Bessel function kn(z) and its derivative for
    all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z,0):
        nm,In,Inp,kn,knp = specfun.csphik(n1,z)
    else:
        nm,kn,knp = specfun.sphk(n1,z)
    return kn[:(n+1)], knp[:(n+1)]


def sph_inkn(n,z):
    """Compute the spherical Bessel functions, in(z) and kn(z) and their
    derivatives for all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z,0):
        nm,In,Inp,kn,knp = specfun.csphik(n1,z)
    else:
        nm,In,Inp = specfun.sphi(n1,z)
        nm,kn,knp = specfun.sphk(n1,z)
    return In[:(n+1)],Inp[:(n+1)],kn[:(n+1)],knp[:(n+1)]


def riccati_jn(n,x):
    """Compute the Ricatti-Bessel function of the first kind and its
    derivative for all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(x)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n == 0):
        n1 = 1
    else:
        n1 = n
    nm,jn,jnp = specfun.rctj(n1,x)
    return jn[:(n+1)],jnp[:(n+1)]


def riccati_yn(n,x):
    """Compute the Ricatti-Bessel function of the second kind and its
    derivative for all orders up to and including n.
    """
    if not (isscalar(n) and isscalar(x)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n == 0):
        n1 = 1
    else:
        n1 = n
    nm,jn,jnp = specfun.rcty(n1,x)
    return jn[:(n+1)],jnp[:(n+1)]


def erfinv(y):
    """
    Inverse function for erf
    """
    return ndtri((y+1)/2.0)/sqrt(2)


def erfcinv(y):
    """
    Inverse function for erfc
    """
    return -ndtri(0.5*y)/sqrt(2)


def erf_zeros(nt):
    """Compute nt complex zeros of the error function erf(z).
    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.cerzo(nt)


def fresnelc_zeros(nt):
    """Compute nt complex zeros of the cosine Fresnel integral C(z).
    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.fcszo(1,nt)


def fresnels_zeros(nt):
    """Compute nt complex zeros of the sine Fresnel integral S(z).
    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.fcszo(2,nt)


def fresnel_zeros(nt):
    """Compute nt complex zeros of the sine and cosine Fresnel integrals
    S(z) and C(z).
    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.fcszo(2,nt), specfun.fcszo(1,nt)


def hyp0f1(v, z):
    r"""Confluent hypergeometric limit function 0F1.

    Parameters
    ----------
    v, z : array_like
        Input values.

    Returns
    -------
    hyp0f1 : ndarray
        The confluent hypergeometric limit function.

    Notes
    -----
    This function is defined as:

    .. math:: _0F_1(v,z) = \sum_{k=0}^{\inf}\frac{z^k}{(v)_k k!}.

    It's also the limit as q -> infinity of ``1F1(q;v;z/q)``, and satisfies
    the differential equation :math:``f''(z) + vf'(z) = f(z)`.
    """
    v = atleast_1d(v)
    z = atleast_1d(z)
    v, z = np.broadcast_arrays(v, z)
    arg = 2 * sqrt(abs(z))
    old_err = np.seterr(all='ignore')  # for z=0, a<1 and num=inf, next lines
    num = where(z.real >= 0, iv(v - 1, arg), jv(v - 1, arg))
    den = abs(z)**((v - 1.0) / 2)
    num *= gamma(v)
    np.seterr(**old_err)
    num[z == 0] = 1
    den[z == 0] = 1
    return num / den


def assoc_laguerre(x, n, k=0.0):
    """Returns the n-th order generalized (associated) Laguerre polynomial.

    The polynomial :math:`L^(alpha)_n(x)` is orthogonal over ``[0, inf)``,
    with weighting function ``exp(-x) * x**alpha`` with ``alpha > -1``.

    Notes
    -----
    `assoc_laguerre` is a simple wrapper around `eval_genlaguerre`, with
    reversed argument order ``(x, n, k=0.0) --> (n, k, x)``.

    """
    return orthogonal.eval_genlaguerre(n, k, x)

digamma = psi


def polygamma(n, x):
    """Polygamma function which is the nth derivative of the digamma (psi)
    function.

    Parameters
    ----------
    n : array_like of int
        The order of the derivative of `psi`.
    x : array_like
        Where to evaluate the polygamma function.

    Returns
    -------
    polygamma : ndarray
        The result.

    Examples
    --------
    >>> from scipy import special
    >>> x = [2, 3, 25.5]
    >>> special.polygamma(1, x)
    array([ 0.64493407,  0.39493407,  0.03999467])
    >>> special.polygamma(0, x) == special.psi(x)
    array([ True,  True,  True], dtype=bool)

    """
    n, x = asarray(n), asarray(x)
    fac2 = (-1.0)**(n+1) * gamma(n+1.0) * zeta(n+1,x)
    return where(n == 0, psi(x), fac2)


def mathieu_even_coef(m,q):
    """Compute expansion coefficients for even Mathieu functions and
    modified Mathieu functions.
    """
    if not (isscalar(m) and isscalar(q)):
        raise ValueError("m and q must be scalars.")
    if (q < 0):
        raise ValueError("q >=0")
    if (m != floor(m)) or (m < 0):
        raise ValueError("m must be an integer >=0.")

    if (q <= 1):
        qm = 7.5+56.1*sqrt(q)-134.7*q+90.7*sqrt(q)*q
    else:
        qm = 17.0+3.1*sqrt(q)-.126*q+.0037*sqrt(q)*q
    km = int(qm+0.5*m)
    if km > 251:
        print("Warning, too many predicted coefficients.")
    kd = 1
    m = int(floor(m))
    if m % 2:
        kd = 2

    a = mathieu_a(m,q)
    fc = specfun.fcoef(kd,m,q,a)
    return fc[:km]


def mathieu_odd_coef(m,q):
    """Compute expansion coefficients for even Mathieu functions and
    modified Mathieu functions.
    """
    if not (isscalar(m) and isscalar(q)):
        raise ValueError("m and q must be scalars.")
    if (q < 0):
        raise ValueError("q >=0")
    if (m != floor(m)) or (m <= 0):
        raise ValueError("m must be an integer > 0")

    if (q <= 1):
        qm = 7.5+56.1*sqrt(q)-134.7*q+90.7*sqrt(q)*q
    else:
        qm = 17.0+3.1*sqrt(q)-.126*q+.0037*sqrt(q)*q
    km = int(qm+0.5*m)
    if km > 251:
        print("Warning, too many predicted coefficients.")
    kd = 4
    m = int(floor(m))
    if m % 2:
        kd = 3

    b = mathieu_b(m,q)
    fc = specfun.fcoef(kd,m,q,b)
    return fc[:km]


def lpmn(m,n,z):
    """Associated Legendre function of the first kind, Pmn(z)

    Computes the associated Legendre function of the first kind
    of order m and degree n,::

        Pmn(z) = P_n^m(z)

    and its derivative, ``Pmn'(z)``.  Returns two arrays of size
    ``(m+1, n+1)`` containing ``Pmn(z)`` and ``Pmn'(z)`` for all
    orders from ``0..m`` and degrees from ``0..n``.

    This function takes a real argument ``z``. For complex arguments ``z``
    use clpmn instead.

    Parameters
    ----------
    m : int
       ``|m| <= n``; the order of the Legendre function.
    n : int
       where ``n >= 0``; the degree of the Legendre function.  Often
       called ``l`` (lower case L) in descriptions of the associated
       Legendre function
    z : float
        Input value.

    Returns
    -------
    Pmn_z : (m+1, n+1) array
       Values for all orders 0..m and degrees 0..n
    Pmn_d_z : (m+1, n+1) array
       Derivatives for all orders 0..m and degrees 0..n

    See Also
    --------
    clpmn: associated Legendre functions of the first kind for complex z

    Notes
    -----
    In the interval (-1, 1), Ferrer's function of the first kind is
    returned. The phase convention used for the intervals (1, inf)
    and (-inf, -1) is such that the result is always real.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           http://dlmf.nist.gov/14.3

    """
    if not isscalar(m) or (abs(m) > n):
        raise ValueError("m must be <= n.")
    if not isscalar(n) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if not isscalar(z):
        raise ValueError("z must be scalar.")
    if iscomplex(z):
        raise ValueError("Argument must be real. Use clpmn instead.")
    if (m < 0):
        mp = -m
        mf,nf = mgrid[0:mp+1,0:n+1]
        sv = errprint(0)
        if abs(z) < 1:
            # Ferrer function; DLMF 14.9.3
            fixarr = where(mf > nf,0.0,(-1)**mf * gamma(nf-mf+1) / gamma(nf+mf+1))
        else:
            # Match to clpmn; DLMF 14.9.13
            fixarr = where(mf > nf,0.0, gamma(nf-mf+1) / gamma(nf+mf+1))
        sv = errprint(sv)
    else:
        mp = m
    p,pd = specfun.lpmn(mp,n,z)
    if (m < 0):
        p = p * fixarr
        pd = pd * fixarr
    return p,pd


def clpmn(m,n,z,type=3):
    """Associated Legendre function of the first kind, Pmn(z)

    Computes the (associated) Legendre function of the first kind
    of order m and degree n,::

        Pmn(z) = P_n^m(z)

    and its derivative, ``Pmn'(z)``.  Returns two arrays of size
    ``(m+1, n+1)`` containing ``Pmn(z)`` and ``Pmn'(z)`` for all
    orders from ``0..m`` and degrees from ``0..n``.

    Parameters
    ----------
    m : int
       ``|m| <= n``; the order of the Legendre function.
    n : int
       where ``n >= 0``; the degree of the Legendre function.  Often
       called ``l`` (lower case L) in descriptions of the associated
       Legendre function
    z : float or complex
        Input value.
    type : int
       takes values 2 or 3
       2: cut on the real axis |x|>1
       3: cut on the real axis -1<x<1 (default)

    Returns
    -------
    Pmn_z : (m+1, n+1) array
       Values for all orders 0..m and degrees 0..n
    Pmn_d_z : (m+1, n+1) array
       Derivatives for all orders 0..m and degrees 0..n

    See Also
    --------
    lpmn: associated Legendre functions of the first kind for real z

    Notes
    -----
    By default, i.e. for ``type=3``, phase conventions are chosen according
    to [1]_ such that the function is analytic. The cut lies on the interval
    (-1, 1). Approaching the cut from above or below in general yields a phase
    factor with respect to Ferrer's function of the first kind
    (cf. `lpmn`).

    For ``type=2`` a cut at |x|>1 is chosen. Approaching the real values
    on the interval (-1, 1) in the complex plane yields Ferrer's function
    of the first kind.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           http://dlmf.nist.gov/14.21

    """
    if not isscalar(m) or (abs(m) > n):
        raise ValueError("m must be <= n.")
    if not isscalar(n) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if not isscalar(z):
        raise ValueError("z must be scalar.")
    if not(type == 2 or type == 3):
        raise ValueError("type must be either 2 or 3.")
    if (m < 0):
        mp = -m
        mf,nf = mgrid[0:mp+1,0:n+1]
        sv = errprint(0)
        if type == 2:
            fixarr = where(mf > nf,0.0, (-1)**mf * gamma(nf-mf+1) / gamma(nf+mf+1))
        else:
            fixarr = where(mf > nf,0.0,gamma(nf-mf+1) / gamma(nf+mf+1))
        sv = errprint(sv)
    else:
        mp = m
    p,pd = specfun.clpmn(mp,n,real(z),imag(z),type)
    if (m < 0):
        p = p * fixarr
        pd = pd * fixarr
    return p,pd


def lqmn(m,n,z):
    """Associated Legendre function of the second kind, Qmn(z)

    Computes the associated Legendre function of the second kind
    of order m and degree n,::

        Qmn(z) = Q_n^m(z)

    and its derivative, ``Qmn'(z)``.  Returns two arrays of size
    ``(m+1, n+1)`` containing ``Qmn(z)`` and ``Qmn'(z)`` for all
    orders from ``0..m`` and degrees from ``0..n``.

    This function takes a real argument ``z``. For complex arguments ``z``
    use clqmn instead.

    Parameters
    ----------
    m : int
       the order of the Legendre function.
    n : int
       where ``n >= 0``; the degree of the Legendre function.  Often
       called ``l`` (lower case L) in descriptions of the associated
       Legendre function
    z : float
        Input value.

    Returns
    -------
    Pmn_z : (m+1, n+1) array
       Values for all orders 0..m and degrees 0..n
    Pmn_d_z : (m+1, n+1) array
       Derivatives for all orders 0..m and degrees 0..n

    See Also
    --------
    clqmn: associated Legendre functions of the second kind for complex z

    Notes
    -----
    In the interval (-1, 1), Ferrer's function of the second kind is
    returned. The phase convention used for the intervals (1, inf)
    and (-inf, -1) is such that the result is always real.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           http://dlmf.nist.gov/14.3

    """
    if not isscalar(m) or (m < 0):
        raise ValueError("m must be a non-negative integer.")
    if not isscalar(n) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if not isscalar(z):
        raise ValueError("z must be scalar.")
    m = int(m)
    n = int(n)

    # Ensure neither m nor n == 0
    mm = max(1,m)
    nn = max(1,n)

    if iscomplex(z):
        raise ValueError("Argument must be real. Use clqmn instead.")
    else:
        q,qd = specfun.lqmn(mm,nn,z)
    return q[:(m+1),:(n+1)],qd[:(m+1),:(n+1)]


def clqmn(m,n,z,type=3):
    """Associated Legendre function of the second kind, Qmn(z)

    Computes the (associated) Legendre function of the second kind
    of order m and degree n,::

        Qmn(z) = Q_n^m(z)

    and its derivative, ``Qmn'(z)``.  Returns two arrays of size
    ``(m+1, n+1)`` containing ``Qmn(z)`` and ``Qmn'(z)`` for all
    orders from ``0..m`` and degrees from ``0..n``.

    Parameters
    ----------
    m : int
    n : int
       where ``n >= 0``; the degree of the Legendre function.  Often
       called ``l`` (lower case L) in descriptions of the associated
       Legendre function
    z : float or complex
        Input value.
    type : int
       takes values 2 or 3
       2: cut on the real axis |x|>1
       3: cut on the real axis -1<x<1 (default)

    Returns
    -------
    Qmn_z : (m+1, n+1) array
       Values for all orders 0..m and degrees 0..n
    Qmn_d_z : (m+1, n+1) array
       Derivatives for all orders 0..m and degrees 0..n

    See Also
    --------
    lqmn: associated Legendre functions of the second kind for real z

    Notes
    -----
    By default, i.e. for ``type=3``, phase conventions are chosen according
    to [1]_ such that the function is analytic. The cut lies on the interval
    (-1, 1). Approaching the cut from above or below in general yields a phase
    factor with respect to Ferrer's function of the second kind
    (cf. `lqmn`).

    For ``type=2`` a cut at |x|>1 is chosen. Approaching the real values
    on the interval (-1, 1) in the complex plane yields Ferrer's function
    of the first kind.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           http://dlmf.nist.gov/14.21

    """
    if not isscalar(m) or (m < 0):
        raise ValueError("m must be a non-negative integer.")
    if not isscalar(n) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if not isscalar(z):
        raise ValueError("z must be scalar.")
    if not(type == 2 or type == 3):
        raise ValueError("type must be either 2 or 3.")

    # Ensure neither m nor n == 0
    mm = max(1,m)
    nn = max(1,n)

    q,qd = specfun.clqmn(mm,nn,real(z),imag(z),type)
    return q[:(m+1),:(n+1)],qd[:(m+1),:(n+1)]


def bernoulli(n):
    """Return an array of the Bernoulli numbers B0..Bn
    """
    if not isscalar(n) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    n = int(n)
    if (n < 2):
        n1 = 2
    else:
        n1 = n
    return specfun.bernob(int(n1))[:(n+1)]


def euler(n):
    """Return an array of the Euler numbers E0..En (inclusive)
    """
    if not isscalar(n) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    n = int(n)
    if (n < 2):
        n1 = 2
    else:
        n1 = n
    return specfun.eulerb(n1)[:(n+1)]


def lpn(n,z):
    """Compute sequence of Legendre functions of the first kind (polynomials),
    Pn(z) and derivatives for all degrees from 0 to n (inclusive).

    See also special.legendre  for polynomial class.
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z):
        pn,pd = specfun.clpn(n1,z)
    else:
        pn,pd = specfun.lpn(n1,z)
    return pn[:(n+1)],pd[:(n+1)]

## lpni


def lqn(n,z):
    """Compute sequence of Legendre functions of the second kind,
    Qn(z) and derivatives for all degrees from 0 to n (inclusive).
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z):
        qn,qd = specfun.clqn(n1,z)
    else:
        qn,qd = specfun.lqnb(n1,z)
    return qn[:(n+1)],qd[:(n+1)]


def ai_zeros(nt):
    """Compute the zeros of Airy Functions Ai(x) and Ai'(x), a and a'
    respectively, and the associated values of Ai(a') and Ai'(a).

    Returns
    -------
    a[l-1]   -- the lth zero of Ai(x)
    ap[l-1]  -- the lth zero of Ai'(x)
    ai[l-1]  -- Ai(ap[l-1])
    aip[l-1] -- Ai'(a[l-1])
    """
    kf = 1
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be a positive integer scalar.")
    return specfun.airyzo(nt,kf)


def bi_zeros(nt):
    """Compute the zeros of Airy Functions Bi(x) and Bi'(x), b and b'
    respectively, and the associated values of Ai(b') and Ai'(b).

    Returns
    -------
    b[l-1]   -- the lth zero of Bi(x)
    bp[l-1]  -- the lth zero of Bi'(x)
    bi[l-1]  -- Bi(bp[l-1])
    bip[l-1] -- Bi'(b[l-1])
    """
    kf = 2
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be a positive integer scalar.")
    return specfun.airyzo(nt,kf)


def lmbda(v,x):
    """Compute sequence of lambda functions with arbitrary order v
    and their derivatives.  Lv0(x)..Lv(x) are computed with v0=v-int(v).
    """
    if not (isscalar(v) and isscalar(x)):
        raise ValueError("arguments must be scalars.")
    if (v < 0):
        raise ValueError("argument must be > 0.")
    n = int(v)
    v0 = v - n
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    v1 = n1 + v0
    if (v != floor(v)):
        vm, vl, dl = specfun.lamv(v1,x)
    else:
        vm, vl, dl = specfun.lamn(v1,x)
    return vl[:(n+1)], dl[:(n+1)]


def pbdv_seq(v,x):
    """Compute sequence of parabolic cylinder functions Dv(x) and
    their derivatives for Dv0(x)..Dv(x) with v0=v-int(v).
    """
    if not (isscalar(v) and isscalar(x)):
        raise ValueError("arguments must be scalars.")
    n = int(v)
    v0 = v-n
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    v1 = n1 + v0
    dv,dp,pdf,pdd = specfun.pbdv(v1,x)
    return dv[:n1+1],dp[:n1+1]


def pbvv_seq(v,x):
    """Compute sequence of parabolic cylinder functions Dv(x) and
    their derivatives for Dv0(x)..Dv(x) with v0=v-int(v).
    """
    if not (isscalar(v) and isscalar(x)):
        raise ValueError("arguments must be scalars.")
    n = int(v)
    v0 = v-n
    if (n <= 1):
        n1 = 1
    else:
        n1 = n
    v1 = n1 + v0
    dv,dp,pdf,pdd = specfun.pbvv(v1,x)
    return dv[:n1+1],dp[:n1+1]


def pbdn_seq(n,z):
    """Compute sequence of parabolic cylinder functions Dn(z) and
    their derivatives for D0(z)..Dn(z).
    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (floor(n) != n):
        raise ValueError("n must be an integer.")
    if (abs(n) <= 1):
        n1 = 1
    else:
        n1 = n
    cpb,cpd = specfun.cpbdn(n1,z)
    return cpb[:n1+1],cpd[:n1+1]


def ber_zeros(nt):
    """Compute nt zeros of the Kelvin function ber x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,1)


def bei_zeros(nt):
    """Compute nt zeros of the Kelvin function bei x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,2)


def ker_zeros(nt):
    """Compute nt zeros of the Kelvin function ker x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,3)


def kei_zeros(nt):
    """Compute nt zeros of the Kelvin function kei x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,4)


def berp_zeros(nt):
    """Compute nt zeros of the Kelvin function ber' x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,5)


def beip_zeros(nt):
    """Compute nt zeros of the Kelvin function bei' x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,6)


def kerp_zeros(nt):
    """Compute nt zeros of the Kelvin function ker' x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,7)


def keip_zeros(nt):
    """Compute nt zeros of the Kelvin function kei' x
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,8)


def kelvin_zeros(nt):
    """Compute nt zeros of all the Kelvin functions returned in a
    length 8 tuple of arrays of length nt.
    The tuple containse the arrays of zeros of
    (ber, bei, ker, kei, ber', bei', ker', kei')
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt,1), \
           specfun.klvnzo(nt,2), \
           specfun.klvnzo(nt,3), \
           specfun.klvnzo(nt,4), \
           specfun.klvnzo(nt,5), \
           specfun.klvnzo(nt,6), \
           specfun.klvnzo(nt,7), \
           specfun.klvnzo(nt,8)


def pro_cv_seq(m,n,c):
    """Compute a sequence of characteristic values for the prolate
    spheroidal wave functions for mode m and n'=m..n and spheroidal
    parameter c.
    """
    if not (isscalar(m) and isscalar(n) and isscalar(c)):
        raise ValueError("Arguments must be scalars.")
    if (n != floor(n)) or (m != floor(m)):
        raise ValueError("Modes must be integers.")
    if (n-m > 199):
        raise ValueError("Difference between n and m is too large.")
    maxL = n-m+1
    return specfun.segv(m,n,c,1)[1][:maxL]


def obl_cv_seq(m,n,c):
    """Compute a sequence of characteristic values for the oblate
    spheroidal wave functions for mode m and n'=m..n and spheroidal
    parameter c.
    """
    if not (isscalar(m) and isscalar(n) and isscalar(c)):
        raise ValueError("Arguments must be scalars.")
    if (n != floor(n)) or (m != floor(m)):
        raise ValueError("Modes must be integers.")
    if (n-m > 199):
        raise ValueError("Difference between n and m is too large.")
    maxL = n-m+1
    return specfun.segv(m,n,c,-1)[1][:maxL]


def ellipk(m):
    """
    Computes the complete elliptic integral of the first kind.

    This function is defined as

    .. math:: K(m) = \\int_0^{\\pi/2} [1 - m \\sin(t)^2]^{-1/2} dt

    Parameters
    ----------
    m : array_like
        The parameter of the elliptic integral.

    Returns
    -------
    K : array_like
        Value of the elliptic integral.

    Notes
    -----
    For more precision around point m = 1, use `ellipkm1`.

    """
    return ellipkm1(1 - asarray(m))


def agm(a,b):
    """Arithmetic, Geometric Mean

    Start with a_0=a and b_0=b and iteratively compute

    a_{n+1} = (a_n+b_n)/2
    b_{n+1} = sqrt(a_n*b_n)

    until a_n=b_n.   The result is agm(a,b)

    agm(a,b)=agm(b,a)
    agm(a,a) = a
    min(a,b) < agm(a,b) < max(a,b)
    """
    s = a + b + 0.0
    return (pi / 4) * s / ellipkm1(4 * a * b / s ** 2)


def comb(N, k, exact=False, repetition=False):
    """
    The number of combinations of N things taken k at a time.

    This is often expressed as "N choose k".

    Parameters
    ----------
    N : int, ndarray
        Number of things.
    k : int, ndarray
        Number of elements taken.
    exact : bool, optional
        If `exact` is False, then floating point precision is used, otherwise
        exact long integer is computed.
    repetition : bool, optional
        If `repetition` is True, then the number of combinations with
        repetition is computed.

    Returns
    -------
    val : int, ndarray
        The total number of combinations.

    Notes
    -----
    - Array arguments accepted only for exact=False case.
    - If k > N, N < 0, or k < 0, then a 0 is returned.

    Examples
    --------
    >>> k = np.array([3, 4])
    >>> n = np.array([10, 10])
    >>> sc.comb(n, k, exact=False)
    array([ 120.,  210.])
    >>> sc.comb(10, 3, exact=True)
    120L
    >>> sc.comb(10, 3, exact=True, repetition=True)
    220L

    """
    if repetition:
        return comb(N + k - 1, k, exact)
    if exact:
        N = int(N)
        k = int(k)
        if (k > N) or (N < 0) or (k < 0):
            return 0
        val = 1
        for j in xrange(min(k, N-k)):
            val = (val*(N-j))//(j+1)
        return val
    else:
        k,N = asarray(k), asarray(N)
        cond = (k <= N) & (N >= 0) & (k >= 0)
        vals = binom(N, k)
        if isinstance(vals, np.ndarray):
            vals[~cond] = 0
        elif not cond:
            vals = np.float64(0)
        return vals


def perm(N, k, exact=False):
    """
    Permutations of N things taken k at a time, i.e., k-permutations of N.

    It's also known as "partial permutations".

    Parameters
    ----------
    N : int, ndarray
        Number of things.
    k : int, ndarray
        Number of elements taken.
    exact : bool, optional
        If `exact` is False, then floating point precision is used, otherwise
        exact long integer is computed.

    Returns
    -------
    val : int, ndarray
        The number of k-permutations of N.

    Notes
    -----
    - Array arguments accepted only for exact=False case.
    - If k > N, N < 0, or k < 0, then a 0 is returned.

    Examples
    --------
    >>> k = np.array([3, 4])
    >>> n = np.array([10, 10])
    >>> perm(n, k)
    array([  720.,  5040.])
    >>> perm(10, 3, exact=True)
    720

    """
    if exact:
        if (k > N) or (N < 0) or (k < 0):
            return 0
        val = 1
        for i in xrange(N - k + 1, N + 1):
            val *= i
        return val
    else:
        k, N = asarray(k), asarray(N)
        cond = (k <= N) & (N >= 0) & (k >= 0)
        vals = poch(N - k + 1, k)
        if isinstance(vals, np.ndarray):
            vals[~cond] = 0
        elif not cond:
            vals = np.float64(0)
        return vals


def factorial(n,exact=False):
    """
    The factorial function, n! = special.gamma(n+1).

    If exact is 0, then floating point precision is used, otherwise
    exact long integer is computed.

    - Array argument accepted only for exact=False case.
    - If n<0, the return value is 0.

    Parameters
    ----------
    n : int or array_like of ints
        Calculate ``n!``.  Arrays are only supported with `exact` set
        to False.  If ``n < 0``, the return value is 0.
    exact : bool, optional
        The result can be approximated rapidly using the gamma-formula
        above.  If `exact` is set to True, calculate the
        answer exactly using integer arithmetic. Default is False.

    Returns
    -------
    nf : float or int
        Factorial of `n`, as an integer or a float depending on `exact`.

    Examples
    --------
    >>> arr = np.array([3,4,5])
    >>> sc.factorial(arr, exact=False)
    array([   6.,   24.,  120.])
    >>> sc.factorial(5, exact=True)
    120L

    """
    if exact:
        if n < 0:
            return 0
        val = 1
        for k in xrange(1,n+1):
            val *= k
        return val
    else:
        n = asarray(n)
        vals = gamma(n+1)
        return where(n >= 0,vals,0)


def factorial2(n, exact=False):
    """
    Double factorial.

    This is the factorial with every second value skipped, i.e.,
    ``7!! = 7 * 5 * 3 * 1``.  It can be approximated numerically as::

      n!! = special.gamma(n/2+1)*2**((m+1)/2)/sqrt(pi)  n odd
          = 2**(n/2) * (n/2)!                           n even

    Parameters
    ----------
    n : int or array_like
        Calculate ``n!!``.  Arrays are only supported with `exact` set
        to False.  If ``n < 0``, the return value is 0.
    exact : bool, optional
        The result can be approximated rapidly using the gamma-formula
        above (default).  If `exact` is set to True, calculate the
        answer exactly using integer arithmetic.

    Returns
    -------
    nff : float or int
        Double factorial of `n`, as an int or a float depending on
        `exact`.

    Examples
    --------
    >>> factorial2(7, exact=False)
    array(105.00000000000001)
    >>> factorial2(7, exact=True)
    105L

    """
    if exact:
        if n < -1:
            return 0
        if n <= 0:
            return 1
        val = 1
        for k in xrange(n,0,-2):
            val *= k
        return val
    else:
        n = asarray(n)
        vals = zeros(n.shape,'d')
        cond1 = (n % 2) & (n >= -1)
        cond2 = (1-(n % 2)) & (n >= -1)
        oddn = extract(cond1,n)
        evenn = extract(cond2,n)
        nd2o = oddn / 2.0
        nd2e = evenn / 2.0
        place(vals,cond1,gamma(nd2o+1)/sqrt(pi)*pow(2.0,nd2o+0.5))
        place(vals,cond2,gamma(nd2e+1) * pow(2.0,nd2e))
        return vals


def factorialk(n,k,exact=True):
    """
    n(!!...!)  = multifactorial of order k
    k times

    Parameters
    ----------
    n : int
        Calculate multifactorial. If `n` < 0, the return value is 0.
    exact : bool, optional
        If exact is set to True, calculate the answer exactly using
        integer arithmetic.

    Returns
    -------
    val : int
        Multi factorial of `n`.

    Raises
    ------
    NotImplementedError
        Raises when exact is False

    Examples
    --------
    >>> sc.factorialk(5, 1, exact=True)
    120L
    >>> sc.factorialk(5, 3, exact=True)
    10L

    """
    if exact:
        if n < 1-k:
            return 0
        if n <= 0:
            return 1
        val = 1
        for j in xrange(n,0,-k):
            val = val*j
        return val
    else:
        raise NotImplementedError


def entr(x):
    """
    entr(x) = -x*log(x)

    Elementwise function for computing entropy.

    Notes
    -----
    This function is concave.

    """
    return -xlogy(x, x)


def rel_entr(x, y):
    """
    rel_entr(x, y) = x*log(x/y)

    Elementwise function for computing relative entropy.

    Notes
    -----
    This function is jointly convex in (x, y).

    """
    return xlogy(x, x) - xlogy(x, y)


def kl_div(x, y):
    """
    kl_div(x, y) = x*log(x/y) - x + y

    Elementwise function for computing Kullback-Leibler divergence.

    Notes
    -----
    This function is non-negative and is jointly convex in (x, y).
    It is zero exactly when x==y.

    """
    return xlogy(x, x) - xlogy(x, y) - x + y

