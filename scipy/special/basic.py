#
# Author:  Travis Oliphant, 2002
#

from __future__ import division, print_function, absolute_import

import warnings

import numpy as np
from scipy._lib.six import xrange
from numpy import (pi, asarray, floor, isscalar, iscomplex, real, imag, sqrt,
                   where, mgrid, sin, place, issubdtype, extract,
                   less, inexact, nan, zeros, atleast_1d, sinc)
from ._ufuncs import (ellipkm1, mathieu_a, mathieu_b, iv, jv, gamma, psi, zeta,
                      hankel1, hankel2, yv, kv, gammaln, ndtri, errprint, poch,
                      binom)
from . import specfun
from . import orthogonal

__all__ = ['agm', 'ai_zeros', 'assoc_laguerre', 'bei_zeros', 'beip_zeros',
           'ber_zeros', 'bernoulli', 'berp_zeros', 'bessel_diff_formula',
           'bi_zeros', 'clpmn', 'comb', 'digamma', 'diric', 'ellipk',
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
           'SpecialFunctionWarning']


class SpecialFunctionWarning(Warning):
    """Warning that can be issued with ``errprint(True)``"""
    pass
warnings.simplefilter("always", category=SpecialFunctionWarning)


def diric(x, n):
    """Periodic sinc function, also called the Dirichlet function.

    The Dirichlet function is defined as::

        diric(x) = sin(x * n/2) / (n * sin(x / 2)),

    where `n` is a positive integer.

    Parameters
    ----------
    x : array_like
        Input data
    n : int
        Integer defining the periodicity.

    Returns
    -------
    diric : ndarray

    Examples
    --------
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt

    >>> x = np.linspace(-8*np.pi, 8*np.pi, num=201)
    >>> plt.figure(figsize=(8, 8));
    >>> for idx, n in enumerate([2, 3, 4, 9]):
    ...     plt.subplot(2, 2, idx+1)
    ...     plt.plot(x, special.diric(x, n))
    ...     plt.title('diric, n={}'.format(n))
    >>> plt.show()

    The following example demonstrates that `diric` gives the magnitudes
    (modulo the sign and scaling) of the Fourier coefficients of a
    rectangular pulse.

    Suppress output of values that are effectively 0:

    >>> np.set_printoptions(suppress=True)

    Create a signal `x` of length `m` with `k` ones:

    >>> m = 8
    >>> k = 3
    >>> x = np.zeros(m)
    >>> x[:k] = 1

    Use the FFT to compute the Fourier transform of `x`, and
    inspect the magnitudes of the coefficients:

    >>> np.abs(np.fft.fft(x))
    array([ 3.        ,  2.41421356,  1.        ,  0.41421356,  1.        ,
            0.41421356,  1.        ,  2.41421356])

    Now find the same values (up to sign) using `diric`.  We multiply
    by `k` to account for the different scaling conventions of
    `numpy.fft.fft` and `diric`:

    >>> theta = np.linspace(0, 2*np.pi, m, endpoint=False)
    >>> k * special.diric(theta, k)
    array([ 3.        ,  2.41421356,  1.        , -0.41421356, -1.        ,
           -0.41421356,  1.        ,  2.41421356])
    """
    x, n = asarray(x), asarray(n)
    n = asarray(n + (x-x))
    x = asarray(x + (n-n))
    if issubdtype(x.dtype, inexact):
        ytype = x.dtype
    else:
        ytype = float
    y = zeros(x.shape, ytype)

    # empirical minval for 32, 64 or 128 bit float computations
    # where sin(x/2) < minval, result is fixed at +1 or -1
    if np.finfo(ytype).eps < 1e-18:
        minval = 1e-11
    elif np.finfo(ytype).eps < 1e-15:
        minval = 1e-7
    else:
        minval = 1e-3

    mask1 = (n <= 0) | (n != floor(n))
    place(y, mask1, nan)

    x = x / 2
    denom = sin(x)
    mask2 = (1-mask1) & (abs(denom) < minval)
    xsub = extract(mask2, x)
    nsub = extract(mask2, n)
    zsub = xsub / pi
    place(y, mask2, pow(-1, np.round(zsub)*(nsub-1)))

    mask = (1-mask1) & (1-mask2)
    xsub = extract(mask, x)
    nsub = extract(mask, n)
    dsub = extract(mask, denom)
    place(y, mask, sin(nsub*xsub)/(nsub*dsub))
    return y


def jnjnp_zeros(nt):
    """Compute zeros of integer-order Bessel functions Jn and Jn'.

    Results are arranged in order of the magnitudes of the zeros.

    Parameters
    ----------
    nt : int
        Number (<=1200) of zeros to compute

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

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt > 1200):
        raise ValueError("Number must be integer <= 1200.")
    nt = int(nt)
    n, m, t, zo = specfun.jdzo(nt)
    return zo[1:nt+1], n[:nt], m[:nt], t[:nt]


def jnyn_zeros(n, nt):
    """Compute nt zeros of Bessel functions Jn(x), Jn'(x), Yn(x), and Yn'(x).

    Returns 4 arrays of length `nt`, corresponding to the first `nt` zeros of
    Jn(x), Jn'(x), Yn(x), and Yn'(x), respectively.

    Parameters
    ----------
    n : int
        Order of the Bessel functions
    nt : int
        Number (<=1200) of zeros to compute

    See jn_zeros, jnp_zeros, yn_zeros, ynp_zeros to get separate arrays.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(nt) and isscalar(n)):
        raise ValueError("Arguments must be scalars.")
    if (floor(n) != n) or (floor(nt) != nt):
        raise ValueError("Arguments must be integers.")
    if (nt <= 0):
        raise ValueError("nt > 0")
    return specfun.jyzo(abs(n), nt)


def jn_zeros(n, nt):
    """Compute zeros of integer-order Bessel function Jn(x).

    Parameters
    ----------
    n : int
        Order of Bessel function
    nt : int
        Number of zeros to return

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    return jnyn_zeros(n, nt)[0]


def jnp_zeros(n, nt):
    """Compute zeros of integer-order Bessel function derivative Jn'(x).

    Parameters
    ----------
    n : int
        Order of Bessel function
    nt : int
        Number of zeros to return

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    return jnyn_zeros(n, nt)[1]


def yn_zeros(n, nt):
    """Compute zeros of integer-order Bessel function Yn(x).

    Parameters
    ----------
    n : int
        Order of Bessel function
    nt : int
        Number of zeros to return

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    return jnyn_zeros(n, nt)[2]


def ynp_zeros(n, nt):
    """Compute zeros of integer-order Bessel function derivative Yn'(x).

    Parameters
    ----------
    n : int
        Order of Bessel function
    nt : int
        Number of zeros to return

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    return jnyn_zeros(n, nt)[3]


def y0_zeros(nt, complex=False):
    """Compute nt zeros of Bessel function Y0(z), and derivative at each zero.

    The derivatives are given by Y0'(z0) = -Y1(z0) at each zero z0.

    Parameters
    ----------
    nt : int
        Number of zeros to return
    complex : bool, default False
        Set to False to return only the real zeros; set to True to return only
        the complex zeros with negative real part and positive imaginary part.
        Note that the complex conjugates of the latter are also zeros of the
        function, but are not returned by this routine.

    Returns
    -------
    z0n : ndarray
        Location of nth zero of Y0(z)
    y0pz0n : ndarray
        Value of derivative Y0'(z0) for nth zero

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("Arguments must be scalar positive integer.")
    kf = 0
    kc = not complex
    return specfun.cyzo(nt, kf, kc)


def y1_zeros(nt, complex=False):
    """Compute nt zeros of Bessel function Y1(z), and derivative at each zero.

    The derivatives are given by Y1'(z1) = Y0(z1) at each zero z1.

    Parameters
    ----------
    nt : int
        Number of zeros to return
    complex : bool, default False
        Set to False to return only the real zeros; set to True to return only
        the complex zeros with negative real part and positive imaginary part.
        Note that the complex conjugates of the latter are also zeros of the
        function, but are not returned by this routine.

    Returns
    -------
    z1n : ndarray
        Location of nth zero of Y1(z)
    y1pz1n : ndarray
        Value of derivative Y1'(z1) for nth zero

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("Arguments must be scalar positive integer.")
    kf = 1
    kc = not complex
    return specfun.cyzo(nt, kf, kc)


def y1p_zeros(nt, complex=False):
    """Compute nt zeros of Bessel derivative Y1'(z), and value at each zero.

    The values are given by Y1(z1) at each z1 where Y1'(z1)=0.

    Parameters
    ----------
    nt : int
        Number of zeros to return
    complex : bool, default False
        Set to False to return only the real zeros; set to True to return only
        the complex zeros with negative real part and positive imaginary part.
        Note that the complex conjugates of the latter are also zeros of the
        function, but are not returned by this routine.

    Returns
    -------
    z1pn : ndarray
        Location of nth zero of Y1'(z)
    y1z1pn : ndarray
        Value of derivative Y1(z1) for nth zero

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("Arguments must be scalar positive integer.")
    kf = 2
    kc = not complex
    return specfun.cyzo(nt, kf, kc)


def _bessel_diff_formula(v, z, n, L, phase):
    # from AMS55.
    # L(v, z) = J(v, z), Y(v, z), H1(v, z), H2(v, z), phase = -1
    # L(v, z) = I(v, z) or exp(v*pi*i)K(v, z), phase = 1
    # For K, you can pull out the exp((v-k)*pi*i) into the caller
    v = asarray(v)
    p = 1.0
    s = L(v-n, z)
    for i in xrange(1, n+1):
        p = phase * (p * (n-i+1)) / i   # = choose(k, i)
        s += p*L(v-n + i*2, z)
    return s / (2.**n)


bessel_diff_formula = np.deprecate(_bessel_diff_formula,
    message="bessel_diff_formula is a private function, do not use it!")


def jvp(v, z, n=1):
    """Compute nth derivative of Bessel function Jv(z) with respect to `z`.

    Parameters
    ----------
    v : float
        Order of Bessel function
    z : complex
        Argument at which to evaluate the derivative
    n : int, default 1
        Order of derivative

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isinstance(n, int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return jv(v, z)
    else:
        return _bessel_diff_formula(v, z, n, jv, -1)


def yvp(v, z, n=1):
    """Compute nth derivative of Bessel function Yv(z) with respect to `z`.

    Parameters
    ----------
    v : float
        Order of Bessel function
    z : complex
        Argument at which to evaluate the derivative
    n : int, default 1
        Order of derivative

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isinstance(n, int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return yv(v, z)
    else:
        return _bessel_diff_formula(v, z, n, yv, -1)


def kvp(v, z, n=1):
    """Compute nth derivative of real-order modified Bessel function Kv(z)

    Kv(z) is the modified Bessel function of the second kind.
    Derivative is calculated with respect to `z`.

    Parameters
    ----------
    v : array_like of float
        Order of Bessel function
    z : array_like of complex
        Argument at which to evaluate the derivative
    n : int
        Order of derivative.  Default is first derivative.

    Returns
    -------
    out : ndarray
        The results

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 6.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isinstance(n, int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return kv(v, z)
    else:
        return (-1)**n * _bessel_diff_formula(v, z, n, kv, 1)


def ivp(v, z, n=1):
    """Compute nth derivative of modified Bessel function Iv(z) with respect
    to `z`.

    Parameters
    ----------
    v : array_like of float
        Order of Bessel function
    z : array_like of complex
        Argument at which to evaluate the derivative
    n : int, default 1
        Order of derivative

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 6.
           http://jin.ece.illinois.edu/specfunc.html

    Examples
    --------
    Calculate multiple values at order 5:

    >>> from scipy.special import kvp
    >>> kvp(5, (1, 2, 3+5j))
    array([-1849.0354+0.j    ,   -25.7735+0.j    ,    -0.0307+0.0875j])

    Calculate for a single value at multiple orders:

    >>> kvp((4, 4.5, 5), 1)
    array([ -184.0309,  -568.9585, -1849.0354])

    """
    if not isinstance(n, int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return iv(v, z)
    else:
        return _bessel_diff_formula(v, z, n, iv, 1)


def h1vp(v, z, n=1):
    """Compute nth derivative of Hankel function H1v(z) with respect to `z`.

    Parameters
    ----------
    v : float
        Order of Hankel function
    z : complex
        Argument at which to evaluate the derivative
    n : int, default 1
        Order of derivative

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isinstance(n, int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return hankel1(v, z)
    else:
        return _bessel_diff_formula(v, z, n, hankel1, -1)


def h2vp(v, z, n=1):
    """Compute nth derivative of Hankel function H2v(z) with respect to `z`.

    Parameters
    ----------
    v : float
        Order of Hankel function
    z : complex
        Argument at which to evaluate the derivative
    n : int, default 1
        Order of derivative

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 5.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isinstance(n, int) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if n == 0:
        return hankel2(v, z)
    else:
        return _bessel_diff_formula(v, z, n, hankel2, -1)


def sph_jn(n, z):
    """Compute spherical Bessel function jn(z) and derivative.

    This function computes the value and first derivative of jn(z) for all
    orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of jn to compute
    z : complex
        Argument at which to evaluate

    Returns
    -------
    jn : ndarray
        Value of j0(z), ..., jn(z)
    jnp : ndarray
        First derivative j0'(z), ..., jn'(z)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 8.
           http://jin.ece.illinois.edu/specfunc.html

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
        nm, jn, jnp, yn, ynp = specfun.csphjy(n1, z)
    else:
        nm, jn, jnp = specfun.sphj(n1, z)
    return jn[:(n+1)], jnp[:(n+1)]


def sph_yn(n, z):
    """Compute spherical Bessel function yn(z) and derivative.

    This function computes the value and first derivative of yn(z) for all
    orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of yn to compute
    z : complex
        Argument at which to evaluate

    Returns
    -------
    yn : ndarray
        Value of y0(z), ..., yn(z)
    ynp : ndarray
        First derivative y0'(z), ..., yn'(z)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 8.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z, 0):
        nm, jn, jnp, yn, ynp = specfun.csphjy(n1, z)
    else:
        nm, yn, ynp = specfun.sphy(n1, z)
    return yn[:(n+1)], ynp[:(n+1)]


def sph_jnyn(n, z):
    """Compute spherical Bessel functions jn(z) and yn(z) and derivatives.

    This function computes the value and first derivative of jn(z) and yn(z)
    for all orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of jn and yn to compute
    z : complex
        Argument at which to evaluate

    Returns
    -------
    jn : ndarray
        Value of j0(z), ..., jn(z)
    jnp : ndarray
        First derivative j0'(z), ..., jn'(z)
    yn : ndarray
        Value of y0(z), ..., yn(z)
    ynp : ndarray
        First derivative y0'(z), ..., yn'(z)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 8.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z, 0):
        nm, jn, jnp, yn, ynp = specfun.csphjy(n1, z)
    else:
        nm, yn, ynp = specfun.sphy(n1, z)
        nm, jn, jnp = specfun.sphj(n1, z)
    return jn[:(n+1)], jnp[:(n+1)], yn[:(n+1)], ynp[:(n+1)]


def sph_in(n, z):
    """Compute spherical Bessel function in(z) and derivative.

    This function computes the value and first derivative of in(z) for all
    orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of in to compute
    z : complex
        Argument at which to evaluate

    Returns
    -------
    in : ndarray
        Value of i0(z), ..., in(z)
    inp : ndarray
        First derivative i0'(z), ..., in'(z)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 8.
           http://jin.ece.illinois.edu/specfunc.html

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
        nm, In, Inp, kn, knp = specfun.csphik(n1, z)
    else:
        nm, In, Inp = specfun.sphi(n1, z)
    return In[:(n+1)], Inp[:(n+1)]


def sph_kn(n, z):
    """Compute spherical Bessel function kn(z) and derivative.

    This function computes the value and first derivative of kn(z) for all
    orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of kn to compute
    z : complex
        Argument at which to evaluate

    Returns
    -------
    kn : ndarray
        Value of k0(z), ..., kn(z)
    knp : ndarray
        First derivative k0'(z), ..., kn'(z)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 8.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z, 0):
        nm, In, Inp, kn, knp = specfun.csphik(n1, z)
    else:
        nm, kn, knp = specfun.sphk(n1, z)
    return kn[:(n+1)], knp[:(n+1)]


def sph_inkn(n, z):
    """Compute spherical Bessel functions in(z), kn(z), and derivatives.

    This function computes the value and first derivative of in(z) and kn(z)
    for all orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of in and kn to compute
    z : complex
        Argument at which to evaluate

    Returns
    -------
    in : ndarray
        Value of i0(z), ..., in(z)
    inp : ndarray
        First derivative i0'(z), ..., in'(z)
    kn : ndarray
        Value of k0(z), ..., kn(z)
    knp : ndarray
        First derivative k0'(z), ..., kn'(z)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 8.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n < 1):
        n1 = 1
    else:
        n1 = n
    if iscomplex(z) or less(z, 0):
        nm, In, Inp, kn, knp = specfun.csphik(n1, z)
    else:
        nm, In, Inp = specfun.sphi(n1, z)
        nm, kn, knp = specfun.sphk(n1, z)
    return In[:(n+1)], Inp[:(n+1)], kn[:(n+1)], knp[:(n+1)]


def riccati_jn(n, x):
    """Compute Ricatti-Bessel function of the first kind and derivative.

    This function computes the value and first derivative of the function for
    all orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of function to compute
    x : float
        Argument at which to evaluate

    Returns
    -------
    jn : ndarray
        Value of j0(x), ..., jn(x)
    jnp : ndarray
        First derivative j0'(x), ..., jn'(x)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(n) and isscalar(x)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n == 0):
        n1 = 1
    else:
        n1 = n
    nm, jn, jnp = specfun.rctj(n1, x)
    return jn[:(n+1)], jnp[:(n+1)]


def riccati_yn(n, x):
    """Compute Ricatti-Bessel function of the second kind and derivative.

    This function computes the value and first derivative of the function for
    all orders up to and including n.

    Parameters
    ----------
    n : int
        Maximum order of function to compute
    x : float
        Argument at which to evaluate

    Returns
    -------
    yn : ndarray
        Value of y0(x), ..., yn(x)
    ynp : ndarray
        First derivative y0'(x), ..., yn'(x)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(n) and isscalar(x)):
        raise ValueError("arguments must be scalars.")
    if (n != floor(n)) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    if (n == 0):
        n1 = 1
    else:
        n1 = n
    nm, jn, jnp = specfun.rcty(n1, x)
    return jn[:(n+1)], jnp[:(n+1)]


def erfinv(y):
    """Inverse function for erf.
    """
    return ndtri((y+1)/2.0)/sqrt(2)


def erfcinv(y):
    """Inverse function for erfc.
    """
    return -ndtri(0.5*y)/sqrt(2)


def erf_zeros(nt):
    """Compute nt complex zeros of error function erf(z).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.cerzo(nt)


def fresnelc_zeros(nt):
    """Compute nt complex zeros of cosine Fresnel integral C(z).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.fcszo(1, nt)


def fresnels_zeros(nt):
    """Compute nt complex zeros of sine Fresnel integral S(z).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.fcszo(2, nt)


def fresnel_zeros(nt):
    """Compute nt complex zeros of sine and cosine Fresnel integrals S(z) and C(z).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if (floor(nt) != nt) or (nt <= 0) or not isscalar(nt):
        raise ValueError("Argument must be positive scalar integer.")
    return specfun.fcszo(2, nt), specfun.fcszo(1, nt)


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

    .. math:: _0F_1(v, z) = \sum_{k=0}^{\inf}\frac{z^k}{(v)_k k!}.

    It's also the limit as q -> infinity of ``1F1(q;v;z/q)``, and satisfies
    the differential equation :math:`f''(z) + vf'(z) = f(z)`.
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
    """Compute nth-order generalized (associated) Laguerre polynomial.

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
    """Polygamma function n.

    This is the nth derivative of the digamma (psi) function.

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
    fac2 = (-1.0)**(n+1) * gamma(n+1.0) * zeta(n+1, x)
    return where(n == 0, psi(x), fac2)


def mathieu_even_coef(m, q):
    r"""Fourier coefficients for even Mathieu and modified Mathieu functions.

    The Fourier series of the even solutions of the Mathieu differential
    equation are of the form

    .. math:: \mathrm{ce}_{2n}(z, q) = \sum_{k=0}^{\infty} A_{(2n)}^{(2k)} \cos 2kz

    .. math:: \mathrm{ce}_{2n+1}(z, q) = \sum_{k=0}^{\infty} A_{(2n+1)}^{(2k+1)} \cos (2k+1)z

    This function returns the coefficients :math:`A_{(2n)}^{(2k)}` for even
    input m=2n, and the coefficients :math:`A_{(2n+1)}^{(2k+1)}` for odd input
    m=2n+1.

    Parameters
    ----------
    m : int
        Order of Mathieu functions.  Must be non-negative.
    q : float (>=0)
        Parameter of Mathieu functions.  Must be non-negative.

    Returns
    -------
    Ak : ndarray
        Even or odd Fourier coefficients, corresponding to even or odd m.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html
    .. [2] NIST Digital Library of Mathematical Functions
           http://dlmf.nist.gov/28.4#i

    """
    if not (isscalar(m) and isscalar(q)):
        raise ValueError("m and q must be scalars.")
    if (q < 0):
        raise ValueError("q >=0")
    if (m != floor(m)) or (m < 0):
        raise ValueError("m must be an integer >=0.")

    if (q <= 1):
        qm = 7.5 + 56.1*sqrt(q) - 134.7*q + 90.7*sqrt(q)*q
    else:
        qm = 17.0 + 3.1*sqrt(q) - .126*q + .0037*sqrt(q)*q
    km = int(qm + 0.5*m)
    if km > 251:
        print("Warning, too many predicted coefficients.")
    kd = 1
    m = int(floor(m))
    if m % 2:
        kd = 2

    a = mathieu_a(m, q)
    fc = specfun.fcoef(kd, m, q, a)
    return fc[:km]


def mathieu_odd_coef(m, q):
    r"""Fourier coefficients for even Mathieu and modified Mathieu functions.

    The Fourier series of the odd solutions of the Mathieu differential
    equation are of the form

    .. math:: \mathrm{se}_{2n+1}(z, q) = \sum_{k=0}^{\infty} B_{(2n+1)}^{(2k+1)} \sin (2k+1)z

    .. math:: \mathrm{se}_{2n+2}(z, q) = \sum_{k=0}^{\infty} B_{(2n+2)}^{(2k+2)} \sin (2k+2)z

    This function returns the coefficients :math:`B_{(2n+2)}^{(2k+2)}` for even
    input m=2n+2, and the coefficients :math:`B_{(2n+1)}^{(2k+1)}` for odd
    input m=2n+1.

    Parameters
    ----------
    m : int
        Order of Mathieu functions.  Must be non-negative.
    q : float (>=0)
        Parameter of Mathieu functions.  Must be non-negative.

    Returns
    -------
    Bk : ndarray
        Even or odd Fourier coefficients, corresponding to even or odd m.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(m) and isscalar(q)):
        raise ValueError("m and q must be scalars.")
    if (q < 0):
        raise ValueError("q >=0")
    if (m != floor(m)) or (m <= 0):
        raise ValueError("m must be an integer > 0")

    if (q <= 1):
        qm = 7.5 + 56.1*sqrt(q) - 134.7*q + 90.7*sqrt(q)*q
    else:
        qm = 17.0 + 3.1*sqrt(q) - .126*q + .0037*sqrt(q)*q
    km = int(qm + 0.5*m)
    if km > 251:
        print("Warning, too many predicted coefficients.")
    kd = 4
    m = int(floor(m))
    if m % 2:
        kd = 3

    b = mathieu_b(m, q)
    fc = specfun.fcoef(kd, m, q, b)
    return fc[:km]


def lpmn(m, n, z):
    """Associated Legendre function of the first kind, Pmn(z).

    Computes the associated Legendre function of the first kind of order m and
    degree n, ``Pmn(z)`` = :math:`P_n^m(z)`, and its derivative, ``Pmn'(z)``.
    Returns two arrays of size ``(m+1, n+1)`` containing ``Pmn(z)`` and
    ``Pmn'(z)`` for all orders from ``0..m`` and degrees from ``0..n``.

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
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html
    .. [2] NIST Digital Library of Mathematical Functions
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
        mf, nf = mgrid[0:mp+1, 0:n+1]
        sv = errprint(0)
        if abs(z) < 1:
            # Ferrer function; DLMF 14.9.3
            fixarr = where(mf > nf, 0.0,
                           (-1)**mf * gamma(nf-mf+1) / gamma(nf+mf+1))
        else:
            # Match to clpmn; DLMF 14.9.13
            fixarr = where(mf > nf, 0.0, gamma(nf-mf+1) / gamma(nf+mf+1))
        sv = errprint(sv)
    else:
        mp = m
    p, pd = specfun.lpmn(mp, n, z)
    if (m < 0):
        p = p * fixarr
        pd = pd * fixarr
    return p, pd


def clpmn(m, n, z, type=3):
    """Associated Legendre function of the first kind, Pmn(z).

    Computes the associated Legendre function of the first kind of order m and
    degree n, ``Pmn(z)`` = :math:`P_n^m(z)`, and its derivative, ``Pmn'(z)``.
    Returns two arrays of size ``(m+1, n+1)`` containing ``Pmn(z)`` and
    ``Pmn'(z)`` for all orders from ``0..m`` and degrees from ``0..n``.

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
    type : int, optional
       takes values 2 or 3
       2: cut on the real axis ``|x| > 1``
       3: cut on the real axis ``-1 < x < 1`` (default)

    Returns
    -------
    Pmn_z : (m+1, n+1) array
       Values for all orders ``0..m`` and degrees ``0..n``
    Pmn_d_z : (m+1, n+1) array
       Derivatives for all orders ``0..m`` and degrees ``0..n``

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

    For ``type=2`` a cut at ``|x| > 1`` is chosen. Approaching the real values
    on the interval (-1, 1) in the complex plane yields Ferrer's function
    of the first kind.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html
    .. [2] NIST Digital Library of Mathematical Functions
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
        mf, nf = mgrid[0:mp+1, 0:n+1]
        sv = errprint(0)
        if type == 2:
            fixarr = where(mf > nf, 0.0,
                           (-1)**mf * gamma(nf-mf+1) / gamma(nf+mf+1))
        else:
            fixarr = where(mf > nf, 0.0, gamma(nf-mf+1) / gamma(nf+mf+1))
        sv = errprint(sv)
    else:
        mp = m
    p, pd = specfun.clpmn(mp, n, real(z), imag(z), type)
    if (m < 0):
        p = p * fixarr
        pd = pd * fixarr
    return p, pd


def lqmn(m, n, z):
    """Associated Legendre function of the second kind, Qmn(z).

    Computes the associated Legendre function of the second kind of order m and
    degree n, ``Qmn(z)`` = :math:`Q_n^m(z)`, and its derivative, ``Qmn'(z)``.
    Returns two arrays of size ``(m+1, n+1)`` containing ``Qmn(z)`` and
    ``Qmn'(z)`` for all orders from ``0..m`` and degrees from ``0..n``.

    Parameters
    ----------
    m : int
       ``|m| <= n``; the order of the Legendre function.
    n : int
       where ``n >= 0``; the degree of the Legendre function.  Often
       called ``l`` (lower case L) in descriptions of the associated
       Legendre function
    z : complex
        Input value.

    Returns
    -------
    Qmn_z : (m+1, n+1) array
       Values for all orders 0..m and degrees 0..n
    Qmn_d_z : (m+1, n+1) array
       Derivatives for all orders 0..m and degrees 0..n

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

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
    mm = max(1, m)
    nn = max(1, n)

    if iscomplex(z):
        q, qd = specfun.clqmn(mm, nn, z)
    else:
        q, qd = specfun.lqmn(mm, nn, z)
    return q[:(m+1), :(n+1)], qd[:(m+1), :(n+1)]


def bernoulli(n):
    """Bernoulli numbers B0..Bn (inclusive).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

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
    """Euler numbers E0..En (inclusive).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(n) or (n < 0):
        raise ValueError("n must be a non-negative integer.")
    n = int(n)
    if (n < 2):
        n1 = 2
    else:
        n1 = n
    return specfun.eulerb(n1)[:(n+1)]


def lpn(n, z):
    """Legendre functions of the first kind, Pn(z).

    Compute sequence of Legendre functions of the first kind (polynomials),
    Pn(z) and derivatives for all degrees from 0 to n (inclusive).

    See also special.legendre for polynomial class.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

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
        pn, pd = specfun.clpn(n1, z)
    else:
        pn, pd = specfun.lpn(n1, z)
    return pn[:(n+1)], pd[:(n+1)]


def lqn(n, z):
    """Legendre functions of the second kind, Qn(z).

    Compute sequence of Legendre functions of the second kind, Qn(z) and
    derivatives for all degrees from 0 to n (inclusive).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

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
        qn, qd = specfun.clqn(n1, z)
    else:
        qn, qd = specfun.lqnb(n1, z)
    return qn[:(n+1)], qd[:(n+1)]


def ai_zeros(nt):
    """
    Compute `nt` zeros and values of the Airy function Ai and its derivative.

    Computes the first `nt` zeros, `a`, of the Airy function Ai(x);
    first `nt` zeros, `ap`, of the derivative of the Airy function Ai'(x);
    the corresponding values Ai(a');
    and the corresponding values Ai'(a).

    Parameters
    ----------
    nt : int
        Number of zeros to compute

    Returns
    -------
    a : ndarray
        First `nt` zeros of Ai(x)
    ap : ndarray
        First `nt` zeros of Ai'(x)
    ai : ndarray
        Values of Ai(x) evaluated at first `nt` zeros of Ai'(x)
    aip : ndarray
        Values of Ai'(x) evaluated at first `nt` zeros of Ai(x)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    kf = 1
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be a positive integer scalar.")
    return specfun.airyzo(nt, kf)


def bi_zeros(nt):
    """
    Compute `nt` zeros and values of the Airy function Bi and its derivative.

    Computes the first `nt` zeros, b, of the Airy function Bi(x);
    first `nt` zeros, b', of the derivative of the Airy function Bi'(x);
    the corresponding values Bi(b');
    and the corresponding values Bi'(b).

    Parameters
    ----------
    nt : int
        Number of zeros to compute

    Returns
    -------
    b : ndarray
        First `nt` zeros of Bi(x)
    bp : ndarray
        First `nt` zeros of Bi'(x)
    bi : ndarray
        Values of Bi(x) evaluated at first `nt` zeros of Bi'(x)
    bip : ndarray
        Values of Bi'(x) evaluated at first `nt` zeros of Bi(x)

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    kf = 2
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be a positive integer scalar.")
    return specfun.airyzo(nt, kf)


def lmbda(v, x):
    """Jahnke-Emden Lambda function, Lambdav(x).

    Parameters
    ----------
    v : float
        Order of the Lambda function
    x : float
        Value at which to evaluate the function and derivatives

    Returns
    -------
    vl : ndarray
        Values of Lambda_vi(x), for vi=v-int(v), vi=1+v-int(v), ..., vi=v.
    dl : ndarray
        Derivatives Lambda_vi'(x), for vi=v-int(v), vi=1+v-int(v), ..., vi=v.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

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
        vm, vl, dl = specfun.lamv(v1, x)
    else:
        vm, vl, dl = specfun.lamn(v1, x)
    return vl[:(n+1)], dl[:(n+1)]


def pbdv_seq(v, x):
    """Parabolic cylinder functions Dv(x) and derivatives.

    Parameters
    ----------
    v : float
        Order of the parabolic cylinder function
    x : float
        Value at which to evaluate the function and derivatives

    Returns
    -------
    dv : ndarray
        Values of D_vi(x), for vi=v-int(v), vi=1+v-int(v), ..., vi=v.
    dp : ndarray
        Derivatives D_vi'(x), for vi=v-int(v), vi=1+v-int(v), ..., vi=v.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 13.
           http://jin.ece.illinois.edu/specfunc.html

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
    dv, dp, pdf, pdd = specfun.pbdv(v1, x)
    return dv[:n1+1], dp[:n1+1]


def pbvv_seq(v, x):
    """Parabolic cylinder functions Vv(x) and derivatives.

    Parameters
    ----------
    v : float
        Order of the parabolic cylinder function
    x : float
        Value at which to evaluate the function and derivatives

    Returns
    -------
    dv : ndarray
        Values of V_vi(x), for vi=v-int(v), vi=1+v-int(v), ..., vi=v.
    dp : ndarray
        Derivatives V_vi'(x), for vi=v-int(v), vi=1+v-int(v), ..., vi=v.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 13.
           http://jin.ece.illinois.edu/specfunc.html

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
    dv, dp, pdf, pdd = specfun.pbvv(v1, x)
    return dv[:n1+1], dp[:n1+1]


def pbdn_seq(n, z):
    """Parabolic cylinder functions Dn(z) and derivatives.

    Parameters
    ----------
    n : int
        Order of the parabolic cylinder function
    z : complex
        Value at which to evaluate the function and derivatives

    Returns
    -------
    dv : ndarray
        Values of D_i(z), for i=0, ..., i=n.
    dp : ndarray
        Derivatives D_i'(z), for i=0, ..., i=n.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996, chapter 13.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(n) and isscalar(z)):
        raise ValueError("arguments must be scalars.")
    if (floor(n) != n):
        raise ValueError("n must be an integer.")
    if (abs(n) <= 1):
        n1 = 1
    else:
        n1 = n
    cpb, cpd = specfun.cpbdn(n1, z)
    return cpb[:n1+1], cpd[:n1+1]


def ber_zeros(nt):
    """Compute nt zeros of the Kelvin function ber(x).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 1)


def bei_zeros(nt):
    """Compute nt zeros of the Kelvin function bei(x).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 2)


def ker_zeros(nt):
    """Compute nt zeros of the Kelvin function ker(x).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 3)


def kei_zeros(nt):
    """Compute nt zeros of the Kelvin function kei(x).
    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 4)


def berp_zeros(nt):
    """Compute nt zeros of the Kelvin function ber'(x).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 5)


def beip_zeros(nt):
    """Compute nt zeros of the Kelvin function bei'(x).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 6)


def kerp_zeros(nt):
    """Compute nt zeros of the Kelvin function ker'(x).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 7)


def keip_zeros(nt):
    """Compute nt zeros of the Kelvin function kei'(x).

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return specfun.klvnzo(nt, 8)


def kelvin_zeros(nt):
    """Compute nt zeros of all Kelvin functions.

    Returned in a length-8 tuple of arrays of length nt.  The tuple contains
    the arrays of zeros of (ber, bei, ker, kei, ber', bei', ker', kei').

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not isscalar(nt) or (floor(nt) != nt) or (nt <= 0):
        raise ValueError("nt must be positive integer scalar.")
    return (specfun.klvnzo(nt, 1),
            specfun.klvnzo(nt, 2),
            specfun.klvnzo(nt, 3),
            specfun.klvnzo(nt, 4),
            specfun.klvnzo(nt, 5),
            specfun.klvnzo(nt, 6),
            specfun.klvnzo(nt, 7),
            specfun.klvnzo(nt, 8))


def pro_cv_seq(m, n, c):
    """Characteristic values for prolate spheroidal wave functions.

    Compute a sequence of characteristic values for the prolate
    spheroidal wave functions for mode m and n'=m..n and spheroidal
    parameter c.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(m) and isscalar(n) and isscalar(c)):
        raise ValueError("Arguments must be scalars.")
    if (n != floor(n)) or (m != floor(m)):
        raise ValueError("Modes must be integers.")
    if (n-m > 199):
        raise ValueError("Difference between n and m is too large.")
    maxL = n-m+1
    return specfun.segv(m, n, c, 1)[1][:maxL]


def obl_cv_seq(m, n, c):
    """Characteristic values for oblate spheroidal wave functions.

    Compute a sequence of characteristic values for the oblate
    spheroidal wave functions for mode m and n'=m..n and spheroidal
    parameter c.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """
    if not (isscalar(m) and isscalar(n) and isscalar(c)):
        raise ValueError("Arguments must be scalars.")
    if (n != floor(n)) or (m != floor(m)):
        raise ValueError("Modes must be integers.")
    if (n-m > 199):
        raise ValueError("Difference between n and m is too large.")
    maxL = n-m+1
    return specfun.segv(m, n, c, -1)[1][:maxL]


def ellipk(m):
    """Complete elliptic integral of the first kind.

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

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind around m = 1
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind


    """
    return ellipkm1(1 - asarray(m))


def agm(a, b):
    """Arithmetic, Geometric Mean.

    Start with a_0=a and b_0=b and iteratively compute

    a_{n+1} = (a_n+b_n)/2
    b_{n+1} = sqrt(a_n*b_n)

    until a_n=b_n.   The result is agm(a, b)

    agm(a, b)=agm(b, a)
    agm(a, a) = a
    min(a, b) < agm(a, b) < max(a, b)
    """
    s = a + b + 0.0
    return (pi / 4) * s / ellipkm1(4 * a * b / s ** 2)


def comb(N, k, exact=False, repetition=False):
    """The number of combinations of N things taken k at a time.

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
    >>> from scipy.special import comb
    >>> k = np.array([3, 4])
    >>> n = np.array([10, 10])
    >>> comb(n, k, exact=False)
    array([ 120.,  210.])
    >>> comb(10, 3, exact=True)
    120L
    >>> comb(10, 3, exact=True, repetition=True)
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
        k, N = asarray(k), asarray(N)
        cond = (k <= N) & (N >= 0) & (k >= 0)
        vals = binom(N, k)
        if isinstance(vals, np.ndarray):
            vals[~cond] = 0
        elif not cond:
            vals = np.float64(0)
        return vals


def perm(N, k, exact=False):
    """Permutations of N things taken k at a time, i.e., k-permutations of N.

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
    >>> from scipy.special import perm
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


def factorial(n, exact=False):
    """The factorial function, n! = special.gamma(n+1).

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
    >>> from scipy.special import factorial
    >>> arr = np.array([3, 4, 5])
    >>> factorial(arr, exact=False)
    array([   6.,   24.,  120.])
    >>> factorial(5, exact=True)
    120L

    """
    if exact:
        if n < 0:
            return 0
        val = 1
        for k in xrange(1, n+1):
            val *= k
        return val
    else:
        n = asarray(n)
        vals = gamma(n+1)
        return where(n >= 0, vals, 0)


def factorial2(n, exact=False):
    """Double factorial.

    This is the factorial with every second value skipped.  E.g., ``7!! = 7 * 5
    * 3 * 1``.  It can be approximated numerically as::

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
    >>> from scipy.special import factorial2
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
        for k in xrange(n, 0, -2):
            val *= k
        return val
    else:
        n = asarray(n)
        vals = zeros(n.shape, 'd')
        cond1 = (n % 2) & (n >= -1)
        cond2 = (1-(n % 2)) & (n >= -1)
        oddn = extract(cond1, n)
        evenn = extract(cond2, n)
        nd2o = oddn / 2.0
        nd2e = evenn / 2.0
        place(vals, cond1, gamma(nd2o + 1) / sqrt(pi) * pow(2.0, nd2o + 0.5))
        place(vals, cond2, gamma(nd2e + 1) * pow(2.0, nd2e))
        return vals


def factorialk(n, k, exact=True):
    """Multifactorial of n of order k, n(!!...!).

    This is the multifactorial of n skipping k values.  For example,

      factorialk(17, 4) = 17!!!! = 17 * 13 * 9 * 5 * 1

    In particular, for any integer ``n``, we have

      factorialk(n, 1) = factorial(n)

      factorialk(n, 2) = factorial2(n)

    Parameters
    ----------
    n : int
        Calculate multifactorial. If `n` < 0, the return value is 0.
    k : int
        Order of multifactorial.
    exact : bool, optional
        If exact is set to True, calculate the answer exactly using
        integer arithmetic.

    Returns
    -------
    val : int
        Multifactorial of `n`.

    Raises
    ------
    NotImplementedError
        Raises when exact is False

    Examples
    --------
    >>> from scipy.special import factorialk
    >>> factorialk(5, 1, exact=True)
    120L
    >>> factorialk(5, 3, exact=True)
    10L

    """
    if exact:
        if n < 1-k:
            return 0
        if n <= 0:
            return 1
        val = 1
        for j in xrange(n, 0, -k):
            val = val*j
        return val
    else:
        raise NotImplementedError
