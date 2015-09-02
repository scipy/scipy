# Docstrings for generated ufuncs
#
# The syntax is designed to look like the function add_newdoc is being
# called from numpy.lib, but in this file add_newdoc puts the
# docstrings in a dictionary. This dictionary is used in
# generate_ufuncs.py to generate the docstrings for the ufuncs in
# scipy.special at the C level when the ufuncs are created at compile
# time.

from __future__ import division, print_function, absolute_import

docdict = {}


def get(name):
    return docdict.get(name)


def add_newdoc(place, name, doc):
    docdict['.'.join((place, name))] = doc

add_newdoc("scipy.special", "sph_harm",
    r"""
    sph_harm(m, n, theta, phi)

    Compute spherical harmonics.

    .. math:: Y^m_n(\theta,\phi) = \sqrt{\frac{2n+1}{4\pi}\frac{(n-m)!}{(n+m)!}} e^{i m \theta} P^m_n(\cos(\phi))

    Parameters
    ----------
    m : int
       ``|m| <= n``; the order of the harmonic.
    n : int
       where `n` >= 0; the degree of the harmonic.  This is often called
       ``l`` (lower case L) in descriptions of spherical harmonics.
    theta : float
       [0, 2*pi]; the azimuthal (longitudinal) coordinate.
    phi : float
       [0, pi]; the polar (colatitudinal) coordinate.

    Returns
    -------
    y_mn : complex float
       The harmonic :math:`Y^m_n` sampled at `theta` and `phi`

    Notes
    -----
    There are different conventions for the meaning of input arguments
    `theta` and `phi`.  We take `theta` to be the azimuthal angle and
    `phi` to be the polar angle.  It is common to see the opposite
    convention - that is `theta` as the polar angle and `phi` as the
    azimuthal angle.

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 14.30. http://dlmf.nist.gov/14.30
    """)

add_newdoc("scipy.special", "_ellip_harm",
    """
    Internal function, use `ellip_harm` instead.
    """)

add_newdoc("scipy.special", "_ellip_norm",
    """
    Internal function, use `ellip_norm` instead.
    """)

add_newdoc("scipy.special", "_lambertw",
    """
    Internal function, use `lambertw` instead.
    """)

add_newdoc("scipy.special", "airy",
    """
    airy(z)

    Airy functions and their derivatives.

    Parameters
    ----------
    z : float or complex
        Argument.

    Returns
    -------
    Ai, Aip, Bi, Bip
        Airy functions Ai and Bi, and their derivatives Aip and Bip

    Notes
    -----
    The Airy functions Ai and Bi are two independent solutions of y''(x) = x y.

    """)

add_newdoc("scipy.special", "airye",
    """
    airye(z)

    Exponentially scaled Airy functions and their derivatives.

    Scaling::

        eAi  = Ai  * exp(2.0/3.0*z*sqrt(z))
        eAip = Aip * exp(2.0/3.0*z*sqrt(z))
        eBi  = Bi  * exp(-abs((2.0/3.0*z*sqrt(z)).real))
        eBip = Bip * exp(-abs((2.0/3.0*z*sqrt(z)).real))

    Parameters
    ----------
    z : float or complex
        Argument.

    Returns
    -------
    eAi, eAip, eBi, eBip
        Airy functions Ai and Bi, and their derivatives Aip and Bip

    """)

add_newdoc("scipy.special", "bdtr",
    """
    bdtr(k, n, p)

    Binomial distribution cumulative distribution function.

    Sum of the terms 0 through `k` of the Binomial probability density.

    ::

        y = sum(nCj p**j (1-p)**(n-j), j=0..k)

    Parameters
    ----------
    k, n : int
        Terms to include
    p : float
        Probability

    Returns
    -------
    y : float
        Sum of terms

    """)

add_newdoc("scipy.special", "bdtrc",
    """
    bdtrc(k, n, p)

    Binomial distribution survival function.

    Sum of the terms k+1 through `n` of the Binomial probability density

    ::

        y = sum(nCj p**j (1-p)**(n-j), j=k+1..n)

    Parameters
    ----------
    k, n : int
        Terms to include
    p : float
        Probability

    Returns
    -------
    y : float
        Sum of terms

    """)

add_newdoc("scipy.special", "bdtri",
    """
    bdtri(k, n, y)

    Inverse function to `bdtr` vs. `p`

    Finds probability `p` such that for the cumulative binomial
    probability ``bdtr(k, n, p) == y``.
    """)

add_newdoc("scipy.special", "bdtrik",
    """
    bdtrik(y, n, p)

    Inverse function to `bdtr` vs `k`
    """)

add_newdoc("scipy.special", "bdtrin",
    """
    bdtrin(k, y, p)

    Inverse function to `bdtr` vs `n`
    """)

add_newdoc("scipy.special", "binom",
    """
    binom(n, k)

    Binomial coefficient
    """)

add_newdoc("scipy.special", "btdtria",
    """
    btdtria(p, b, x)

    Inverse of `btdtr` vs `a`

    """)

add_newdoc("scipy.special", "btdtrib",
    """
    btdtria(a, p, x)

    Inverse of `btdtr` vs `b`

    """)

add_newdoc("scipy.special", "bei",
    """
    bei(x)

    Kelvin function bei
    """)

add_newdoc("scipy.special", "beip",
    """
    beip(x)

    Derivative of the Kelvin function `bei`
    """)

add_newdoc("scipy.special", "ber",
    """
    ber(x)

    Kelvin function ber.
    """)

add_newdoc("scipy.special", "berp",
    """
    berp(x)

    Derivative of the Kelvin function `ber`
    """)

add_newdoc("scipy.special", "besselpoly",
    r"""
    besselpoly(a, lmb, nu)

    Weighted integral of a Bessel function.

    .. math::

       \int_0^1 x^\lambda J_\nu(2 a x) \, dx

    where :math:`J_\nu` is a Bessel function and :math:`\lambda=lmb`,
    :math:`\nu=nu`.

    """)

add_newdoc("scipy.special", "beta",
    """
    beta(a, b)

    Beta function.

    ::

        beta(a, b) =  gamma(a) * gamma(b) / gamma(a+b)
    """)

add_newdoc("scipy.special", "betainc",
    """
    betainc(a, b, x)

    Incomplete beta integral.

    Compute the incomplete beta integral of the arguments, evaluated
    from zero to `x`::

        gamma(a+b) / (gamma(a)*gamma(b)) * integral(t**(a-1) (1-t)**(b-1), t=0..x).

    Notes
    -----
    The incomplete beta is also sometimes defined without the terms
    in gamma, in which case the above definition is the so-called regularized
    incomplete beta. Under this definition, you can get the incomplete beta by
    multiplying the result of the scipy function by beta(a, b).

    """)

add_newdoc("scipy.special", "betaincinv",
    """
    betaincinv(a, b, y)

    Inverse function to beta integral.

    Compute `x` such that betainc(a, b, x) = y.
    """)

add_newdoc("scipy.special", "betaln",
    """
    betaln(a, b)

    Natural logarithm of absolute value of beta function.

    Computes ``ln(abs(beta(x)))``.
    """)

add_newdoc("scipy.special", "boxcox",
    """
    boxcox(x, lmbda)

    Compute the Box-Cox transformation.

    The Box-Cox transformation is::

        y = (x**lmbda - 1) / lmbda  if lmbda != 0
            log(x)                  if lmbda == 0

    Returns `nan` if ``x < 0``.
    Returns `-inf` if ``x == 0`` and ``lmbda < 0``.

    Parameters
    ----------
    x : array_like
        Data to be transformed.
    lmbda : array_like
        Power parameter of the Box-Cox transform.

    Returns
    -------
    y : array
        Transformed data.

    Notes
    -----

    .. versionadded:: 0.14.0

    Examples
    --------
    >>> from scipy.special import boxcox
    >>> boxcox([1, 4, 10], 2.5)
    array([   0.        ,   12.4       ,  126.09110641])
    >>> boxcox(2, [0, 1, 2])
    array([ 0.69314718,  1.        ,  1.5       ])
    """)

add_newdoc("scipy.special", "boxcox1p",
    """
    boxcox1p(x, lmbda)

    Compute the Box-Cox transformation of 1 + `x`.

    The Box-Cox transformation computed by `boxcox1p` is::

        y = ((1+x)**lmbda - 1) / lmbda  if lmbda != 0
            log(1+x)                    if lmbda == 0

    Returns `nan` if ``x < -1``.
    Returns `-inf` if ``x == -1`` and ``lmbda < 0``.

    Parameters
    ----------
    x : array_like
        Data to be transformed.
    lmbda : array_like
        Power parameter of the Box-Cox transform.

    Returns
    -------
    y : array
        Transformed data.

    Notes
    -----

    .. versionadded:: 0.14.0

    Examples
    --------
    >>> from scipy.special import boxcox1p
    >>> boxcox1p(1e-4, [0, 0.5, 1])
    array([  9.99950003e-05,   9.99975001e-05,   1.00000000e-04])
    >>> boxcox1p([0.01, 0.1], 0.25)
    array([ 0.00996272,  0.09645476])
    """)

add_newdoc("scipy.special", "inv_boxcox",
    """
    inv_boxcox(y, lmbda)

    Compute the inverse of the Box-Cox transformation.

    Find ``x`` such that::

        y = (x**lmbda - 1) / lmbda  if lmbda != 0
            log(x)                  if lmbda == 0

    Parameters
    ----------
    y : array_like
        Data to be transformed.
    lmbda : array_like
        Power parameter of the Box-Cox transform.

    Returns
    -------
    x : array
        Transformed data.

    Notes
    -----

    .. versionadded:: 0.16.0

    Examples
    --------
    >>> from scipy.special import boxcox, inv_boxcox
    >>> y = boxcox([1, 4, 10], 2.5)
    >>> inv_boxcox(y, 2.5)
    array([1., 4., 10.])
    """)

add_newdoc("scipy.special", "inv_boxcox1p",
    """
    inv_boxcox1p(y, lmbda)

    Compute the inverse of the Box-Cox transformation.

    Find ``x`` such that::

        y = ((1+x)**lmbda - 1) / lmbda  if lmbda != 0
            log(1+x)                    if lmbda == 0

    Parameters
    ----------
    y : array_like
        Data to be transformed.
    lmbda : array_like
        Power parameter of the Box-Cox transform.

    Returns
    -------
    x : array
        Transformed data.

    Notes
    -----

    .. versionadded:: 0.16.0

    Examples
    --------
    >>> from scipy.special import boxcox1p, inv_boxcox1p
    >>> y = boxcox1p([1, 4, 10], 2.5)
    >>> inv_boxcox1p(y, 2.5)
    array([1., 4., 10.])
    """)

add_newdoc("scipy.special", "btdtr",
    """
    btdtr(a, b, x)

    Cumulative beta distribution.

    Returns the area from zero to `x` under the beta density function::

        gamma(a+b)/(gamma(a)*gamma(b)))*integral(t**(a-1) (1-t)**(b-1), t=0..x)

    See Also
    --------
    betainc
    """)

add_newdoc("scipy.special", "btdtri",
    """
    btdtri(a, b, p)

    p-th quantile of the beta distribution.

    This is effectively the inverse of `btdtr` returning the value of `x` for which
    ``btdtr(a, b, x) = p``

    See Also
    --------
    betaincinv
    """)

add_newdoc("scipy.special", "cbrt",
    """
    cbrt(x)

    Cube root of `x`
    """)

add_newdoc("scipy.special", "chdtr",
    """
    chdtr(v, x)

    Chi square cumulative distribution function

    Returns the area under the left hand tail (from 0 to `x`) of the Chi
    square probability density function with `v` degrees of freedom::

        1/(2**(v/2) * gamma(v/2)) * integral(t**(v/2-1) * exp(-t/2), t=0..x)
    """)

add_newdoc("scipy.special", "chdtrc",
    """
    chdtrc(v, x)

    Chi square survival function

    Returns the area under the right hand tail (from `x` to
    infinity) of the Chi square probability density function with `v`
    degrees of freedom::

        1/(2**(v/2) * gamma(v/2)) * integral(t**(v/2-1) * exp(-t/2), t=x..inf)
    """)

add_newdoc("scipy.special", "chdtri",
    """
    chdtri(v, p)

    Inverse to `chdtrc`

    Returns the argument x such that ``chdtrc(v, x) == p``.
    """)

add_newdoc("scipy.special", "chdtriv",
    """
    chdtri(p, x)

    Inverse to `chdtr` vs `v`

    Returns the argument v such that ``chdtr(v, x) == p``.
    """)

add_newdoc("scipy.special", "chndtr",
    """
    chndtr(x, df, nc)

    Non-central chi square cumulative distribution function

    """)

add_newdoc("scipy.special", "chndtrix",
    """
    chndtrix(p, df, nc)

    Inverse to `chndtr` vs `x`
    """)

add_newdoc("scipy.special", "chndtridf",
    """
    chndtridf(x, p, nc)

    Inverse to `chndtr` vs `df`
    """)

add_newdoc("scipy.special", "chndtrinc",
    """
    chndtrinc(x, df, p)

    Inverse to `chndtr` vs `nc`
    """)

add_newdoc("scipy.special", "cosdg",
    """
    cosdg(x)

    Cosine of the angle `x` given in degrees.
    """)

add_newdoc("scipy.special", "cosm1",
    """
    cosm1(x)

    cos(x) - 1 for use when `x` is near zero.
    """)

add_newdoc("scipy.special", "cotdg",
    """
    cotdg(x)

    Cotangent of the angle `x` given in degrees.
    """)

add_newdoc("scipy.special", "dawsn",
    """
    dawsn(x)

    Dawson's integral.

    Computes::

        exp(-x**2) * integral(exp(t**2), t=0..x).

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva
    """)

add_newdoc("scipy.special", "ellipe",
    """
    ellipe(m)

    Complete elliptic integral of the second kind

    This function is defined as

    .. math:: E(m) = \\int_0^{\\pi/2} [1 - m \\sin(t)^2]^{1/2} dt

    Parameters
    ----------
    m : array_like
        Defines the parameter of the elliptic integral.

    Returns
    -------
    E : ndarray
        Value of the elliptic integral.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipeinc : Incomplete elliptic integral of the second kind
    """)

add_newdoc("scipy.special", "ellipeinc",
    """
    ellipeinc(phi, m)

    Incomplete elliptic integral of the second kind

    This function is defined as

    .. math:: E(\\phi, m) = \\int_0^{\\phi} [1 - m \\sin(t)^2]^{1/2} dt

    Parameters
    ----------
    phi : array_like
        amplitude of the elliptic integral.

    m : array_like
        parameter of the elliptic integral.

    Returns
    -------
    E : ndarray
        Value of the elliptic integral.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    """)

add_newdoc("scipy.special", "ellipj",
    """
    ellipj(u, m)

    Jacobian elliptic functions

    Calculates the Jacobian elliptic functions of parameter `m` between
    0 and 1, and real `u`.

    Parameters
    ----------
    m, u
        Parameters

    Returns
    -------
    sn, cn, dn, ph
        The returned functions::

            sn(u|m), cn(u|m), dn(u|m)

        The value ``ph`` is such that if ``u = ellik(ph, m)``,
        then ``sn(u|m) = sin(ph)`` and ``cn(u|m) = cos(ph)``.

    """)

add_newdoc("scipy.special", "ellipkm1",
    """
    ellipkm1(p)

    Complete elliptic integral of the first kind around `m` = 1

    This function is defined as

    .. math:: K(p) = \\int_0^{\\pi/2} [1 - m \\sin(t)^2]^{-1/2} dt

    where `m = 1 - p`.

    Parameters
    ----------
    p : array_like
        Defines the parameter of the elliptic integral as `m` = 1 - p.

    Returns
    -------
    K : ndarray
        Value of the elliptic integral.

    See Also
    --------
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind
    """)

add_newdoc("scipy.special", "ellipkinc",
    """
    ellipkinc(phi, m)

    Incomplete elliptic integral of the first kind

    This function is defined as

    .. math:: K(\\phi, m) = \\int_0^{\\phi} [1 - m \\sin(t)^2]^{-1/2} dt

    Parameters
    ----------
    phi : array_like
        amplitude of the elliptic integral

    m : array_like
        parameter of the elliptic integral

    Returns
    -------
    K : ndarray
        Value of the elliptic integral

    Notes
    -----
    This function is also called ``F(phi, m)``.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind

    """)

add_newdoc("scipy.special", "entr",
    r"""
    entr(x)

    Elementwise function for computing entropy.

    .. math:: \text{entr}(x) = \begin{cases} - x \log(x) & x > 0  \\ 0 & x = 0 \\ -\infty & \text{otherwise} \end{cases}

    Parameters
    ----------
    x : ndarray
        Input array.

    Returns
    -------
    res : ndarray
        The value of the elementwise entropy function at the given points `x`.

    See Also
    --------
    kl_div, rel_entr

    Notes
    -----
    This function is concave.

    .. versionadded:: 0.14.0

    """)

add_newdoc("scipy.special", "erf",
    """
    erf(z)

    Returns the error function of complex argument.

    It is defined as ``2/sqrt(pi)*integral(exp(-t**2), t=0..z)``.

    Parameters
    ----------
    x : ndarray
        Input array.

    Returns
    -------
    res : ndarray
        The values of the error function at the given points `x`.

    See Also
    --------
    erfc, erfinv, erfcinv

    Notes
    -----
    The cumulative of the unit normal distribution is given by
    ``Phi(z) = 1/2[1 + erf(z/sqrt(2))]``.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Error_function
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover,
        1972. http://www.math.sfu.ca/~cbm/aands/page_297.htm
    .. [3] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "erfc",
    """
    erfc(x)

    Complementary error function, 1 - erf(x).

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "erfi",
    """
    erfi(z)

    Imaginary error function, -i erf(i z).

    Notes
    -----

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "erfcx",
    """
    erfcx(x)

    Scaled complementary error function, exp(x^2) erfc(x).

    Notes
    -----

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    """)

add_newdoc("scipy.special", "eval_jacobi",
    """
    eval_jacobi(n, alpha, beta, x, out=None)

    Evaluate Jacobi polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_jacobi",
    """
    eval_sh_jacobi(n, p, q, x, out=None)

    Evaluate shifted Jacobi polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_gegenbauer",
    """
    eval_gegenbauer(n, alpha, x, out=None)

    Evaluate Gegenbauer polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_chebyt",
    """
    eval_chebyt(n, x, out=None)

    Evaluate Chebyshev T polynomial at a point.

    This routine is numerically stable for `x` in ``[-1, 1]`` at least
    up to order ``10000``.
    """)

add_newdoc("scipy.special", "eval_chebyu",
    """
    eval_chebyu(n, x, out=None)

    Evaluate Chebyshev U polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_chebys",
    """
    eval_chebys(n, x, out=None)

    Evaluate Chebyshev S polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_chebyc",
    """
    eval_chebyc(n, x, out=None)

    Evaluate Chebyshev C polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_chebyt",
    """
    eval_sh_chebyt(n, x, out=None)

    Evaluate shifted Chebyshev T polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_chebyu",
    """
    eval_sh_chebyu(n, x, out=None)

    Evaluate shifted Chebyshev U polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_legendre",
    """
    eval_legendre(n, x, out=None)

    Evaluate Legendre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_sh_legendre",
    """
    eval_sh_legendre(n, x, out=None)

    Evaluate shifted Legendre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_genlaguerre",
    """
    eval_genlaguerre(n, alpha, x, out=None)

    Evaluate generalized Laguerre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_laguerre",
     """
    eval_laguerre(n, x, out=None)

    Evaluate Laguerre polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_hermite",
    """
    eval_hermite(n, x, out=None)

    Evaluate Hermite polynomial at a point.
    """)

add_newdoc("scipy.special", "eval_hermitenorm",
    """
    eval_hermitenorm(n, x, out=None)

    Evaluate normalized Hermite polynomial at a point.
    """)

add_newdoc("scipy.special", "exp1",
    """
    exp1(z)

    Exponential integral E_1 of complex argument z

    ::

        integral(exp(-z*t)/t, t=1..inf).
    """)

add_newdoc("scipy.special", "exp10",
    """
    exp10(x)

    10**x
    """)

add_newdoc("scipy.special", "exp2",
    """
    exp2(x)

    2**x
    """)

add_newdoc("scipy.special", "expi",
    """
    expi(x)

    Exponential integral Ei

    Defined as::

        integral(exp(t)/t, t=-inf..x)

    See `expn` for a different exponential integral.
    """)

add_newdoc('scipy.special', 'expit',
    """
    expit(x)

    Expit ufunc for ndarrays.

    The expit function, also known as the logistic function, is defined as
    expit(x) = 1/(1+exp(-x)). It is the inverse of the logit function.

    Parameters
    ----------
    x : ndarray
        The ndarray to apply expit to element-wise.

    Returns
    -------
    out : ndarray
        An ndarray of the same shape as x. Its entries
        are expit of the corresponding entry of x.

    Notes
    -----
    As a ufunc expit takes a number of optional
    keyword arguments. For more information
    see `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_

    .. versionadded:: 0.10.0

    """)

add_newdoc("scipy.special", "expm1",
    """
    expm1(x)

    exp(x) - 1 for use when `x` is near zero.
    """)

add_newdoc("scipy.special", "expn",
    """
    expn(n, x)

    Exponential integral E_n

    Returns the exponential integral for integer `n` and non-negative `x` and
    `n`::

        integral(exp(-x*t) / t**n, t=1..inf).
    """)

add_newdoc("scipy.special", "exprel",
    r"""
    exprel(x)

    Relative error exponential, (exp(x)-1)/x, for use when `x` is near zero.

    Parameters
    ----------
    x : ndarray
        Input array.

    Returns
    -------
    res : ndarray
        Output array.

    See Also
    --------
    expm1

    .. versionadded:: 0.17.0
    """)

add_newdoc("scipy.special", "fdtr",
    """
    fdtr(dfn, dfd, x)

    F cumulative distribution function

    Returns the area from zero to `x` under the F density function (also
    known as Snedcor's density or the variance ratio density).  This
    is the density of X = (unum/dfn)/(uden/dfd), where unum and uden
    are random variables having Chi square distributions with dfn and
    dfd degrees of freedom, respectively.
    """)

add_newdoc("scipy.special", "fdtrc",
    """
    fdtrc(dfn, dfd, x)

    F survival function

    Returns the complemented F distribution function.
    """)

add_newdoc("scipy.special", "fdtri",
    """
    fdtri(dfn, dfd, p)

    Inverse to `fdtr` vs x

    Finds the F density argument `x` such that ``fdtr(dfn, dfd, x) == p``.
    """)

add_newdoc("scipy.special", "fdtridfd",
    """
    fdtridfd(dfn, p, x)

    Inverse to `fdtr` vs dfd

    Finds the F density argument dfd such that ``fdtr(dfn, dfd, x) == p``.
    """)

add_newdoc("scipy.special", "fdtridfn",
    """
    fdtridfn(p, dfd, x)

    Inverse to `fdtr` vs dfn

    finds the F density argument dfn such that ``fdtr(dfn, dfd, x) == p``.
    """)

add_newdoc("scipy.special", "fresnel",
    """
    fresnel(z)

    Fresnel sin and cos integrals

    Defined as::

        ssa = integral(sin(pi/2 * t**2), t=0..z)
        csa = integral(cos(pi/2 * t**2), t=0..z)

    Parameters
    ----------
    z : float or complex array_like
        Argument

    Returns
    -------
    ssa, csa
        Fresnel sin and cos integral values

    """)

add_newdoc("scipy.special", "gamma",
    """
    gamma(z)

    Gamma function

    The gamma function is often referred to as the generalized
    factorial since ``z*gamma(z) = gamma(z+1)`` and ``gamma(n+1) =
    n!`` for natural number *n*.
    """)

add_newdoc("scipy.special", "gammainc",
    """
    gammainc(a, x)

    Incomplete gamma function

    Defined as::

        1 / gamma(a) * integral(exp(-t) * t**(a-1), t=0..x)

    `a` must be positive and `x` must be >= 0.
    """)

add_newdoc("scipy.special", "gammaincc",
    """
    gammaincc(a, x)

    Complemented incomplete gamma integral

    Defined as::

        1 / gamma(a) * integral(exp(-t) * t**(a-1), t=x..inf) = 1 - gammainc(a, x)

    `a` must be positive and `x` must be >= 0.
    """)

add_newdoc("scipy.special", "gammainccinv",
    """
    gammainccinv(a, y)

    Inverse to `gammaincc`

    Returns `x` such that ``gammaincc(a, x) == y``.
    """)

add_newdoc("scipy.special", "gammaincinv",
    """
    gammaincinv(a, y)

    Inverse to `gammainc`

    Returns `x` such that ``gammainc(a, x) = y``.
    """)

add_newdoc("scipy.special", "gammaln",
    """
    gammaln(z)

    Logarithm of absolute value of gamma function

    Defined as::

        ln(abs(gamma(z)))

    See Also
    --------
    gammasgn
    """)

add_newdoc("scipy.special", "gammasgn",
    """
    gammasgn(x)

    Sign of the gamma function.

    See Also
    --------
    gammaln
    """)

add_newdoc("scipy.special", "gdtr",
    """
    gdtr(a, b, x)

    Gamma distribution cumulative density function.

    Returns the integral from zero to `x` of the gamma probability
    density function::

        a**b / gamma(b) * integral(t**(b-1) exp(-at), t=0..x).

    The arguments `a` and `b` are used differently here than in other
    definitions.
    """)

add_newdoc("scipy.special", "gdtrc",
    """
    gdtrc(a, b, x)

    Gamma distribution survival function.

    Integral from `x` to infinity of the gamma probability density
    function.

    See Also
    --------
    gdtr, gdtri
    """)

add_newdoc("scipy.special", "gdtria",
    """
    gdtria(p, b, x, out=None)

    Inverse of `gdtr` vs a.

    Returns the inverse with respect to the parameter `a` of ``p =
    gdtr(a, b, x)``, the cumulative distribution function of the gamma
    distribution.

    Parameters
    ----------
    p : array_like
        Probability values.
    b : array_like
        `b` parameter values of `gdtr(a, b, x)`.  `b` is the "shape" parameter
        of the gamma distribution.
    x : array_like
        Nonnegative real values, from the domain of the gamma distribution.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`.  `out` is then the
        array returned by the function.

    Returns
    -------
    a : ndarray
        Values of the `a` parameter such that `p = gdtr(a, b, x)`.  `1/a`
        is the "scale" parameter of the gamma distribution.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtrib : Inverse with respect to `b` of `gdtr(a, b, x)`.
    gdtrix : Inverse with respect to `x` of `gdtr(a, b, x)`.

    Examples
    --------
    First evaluate `gdtr`.

    >>> from scipy.special import gdtr, gdtria
    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtria(p, 3.4, 5.6)
    1.2
    """)

add_newdoc("scipy.special", "gdtrib",
    """
    gdtrib(a, p, x, out=None)

    Inverse of `gdtr` vs b.

    Returns the inverse with respect to the parameter `b` of ``p =
    gdtr(a, b, x)``, the cumulative distribution function of the gamma
    distribution.

    Parameters
    ----------
    a : array_like
        `a` parameter values of `gdtr(a, b, x)`. `1/a` is the "scale"
        parameter of the gamma distribution.
    p : array_like
        Probability values.
    x : array_like
        Nonnegative real values, from the domain of the gamma distribution.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`.  `out` is then the
        array returned by the function.

    Returns
    -------
    b : ndarray
        Values of the `b` parameter such that `p = gdtr(a, b, x)`.  `b` is
        the "shape" parameter of the gamma distribution.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtria : Inverse with respect to `a` of `gdtr(a, b, x)`.
    gdtrix : Inverse with respect to `x` of `gdtr(a, b, x)`.

    Examples
    --------
    First evaluate `gdtr`.

    >>> from scipy.special import gdtr, gdtrib
    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtrib(1.2, p, 5.6)
    3.3999999999723882
    """)

add_newdoc("scipy.special", "gdtrix",
    """
    gdtrix(a, b, p, out=None)

    Inverse of `gdtr` vs x.

    Returns the inverse with respect to the parameter `x` of ``p =
    gdtr(a, b, x)``, the cumulative distribution function of the gamma
    distribution. This is also known as the p'th quantile of the
    distribution.

    Parameters
    ----------
    a : array_like
        `a` parameter values of `gdtr(a, b, x)`.  `1/a` is the "scale"
        parameter of the gamma distribution.
    b : array_like
        `b` parameter values of `gdtr(a, b, x)`.  `b` is the "shape" parameter
        of the gamma distribution.
    p : array_like
        Probability values.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`.  `out` is then the
        array returned by the function.

    Returns
    -------
    x : ndarray
        Values of the `x` parameter such that `p = gdtr(a, b, x)`.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtria : Inverse with respect to `a` of `gdtr(a, b, x)`.
    gdtrib : Inverse with respect to `b` of `gdtr(a, b, x)`.

    Examples
    --------
    First evaluate `gdtr`.

    >>> from scipy.special import gdtr, gdtrix
    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtrix(1.2, 3.4, p)
    5.5999999999999996
    """)

add_newdoc("scipy.special", "hankel1",
    """
    hankel1(v, z)

    Hankel function of the first kind

    Parameters
    ----------
    v : float
        Order
    z : float or complex
        Argument

    """)

add_newdoc("scipy.special", "hankel1e",
    """
    hankel1e(v, z)

    Exponentially scaled Hankel function of the first kind

    Defined as::

        hankel1e(v, z) = hankel1(v, z) * exp(-1j * z)

    Parameters
    ----------
    v : float
        Order
    z : complex
        Argument
    """)

add_newdoc("scipy.special", "hankel2",
    """
    hankel2(v, z)

    Hankel function of the second kind

    Parameters
    ----------
    v : float
        Order
    z : complex
        Argument

    """)

add_newdoc("scipy.special", "hankel2e",
    """
    hankel2e(v, z)

    Exponentially scaled Hankel function of the second kind

    Defined as::

        hankel1e(v, z) = hankel1(v, z) * exp(1j * z)

    Parameters
    ----------
    v : float
        Order
    z : complex
        Argument
    """)

add_newdoc("scipy.special", "huber",
    r"""
    huber(delta, r)

    Huber loss function.

    .. math:: \text{huber}(\delta, r) = \begin{cases} \infty & \delta < 0  \\ \frac{1}{2}r^2 & 0 \le \delta, | r | \le \delta \\ \delta ( |r| - \frac{1}{2}\delta ) & \text{otherwise} \end{cases}

    Parameters
    ----------
    delta : ndarray
        Input array, indicating the quadratic vs. linear loss changepoint.
    r : ndarray
        Input array, possibly representing residuals.

    Returns
    -------
    res : ndarray
        The computed Huber loss function values.

    Notes
    -----
    This function is convex in r.

    .. versionadded:: 0.15.0

    """)

add_newdoc("scipy.special", "hyp1f1",
    """
    hyp1f1(a, b, x)

    Confluent hypergeometric function 1F1(a, b; x)
    """)

add_newdoc("scipy.special", "hyp1f2",
    """
    hyp1f2(a, b, c, x)

    Hypergeometric function 1F2 and error estimate

    Returns
    -------
    y
        Value of the function
    err
        Error estimate
    """)

add_newdoc("scipy.special", "hyp2f0",
    """
    hyp2f0(a, b, x, type)

    Hypergeometric function 2F0 in y and an error estimate

    The parameter `type` determines a convergence factor and can be
    either 1 or 2.

    Returns
    -------
    y
        Value of the function
    err
        Error estimate
    """)

add_newdoc("scipy.special", "hyp2f1",
    """
    hyp2f1(a, b, c, z)

    Gauss hypergeometric function 2F1(a, b; c; z).
    """)

add_newdoc("scipy.special", "hyp3f0",
    """
    hyp3f0(a, b, c, x)

    Hypergeometric function 3F0 in y and an error estimate

    Returns
    -------
    y
        Value of the function
    err
        Error estimate
    """)

add_newdoc("scipy.special", "hyperu",
    """
    hyperu(a, b, x)

    Confluent hypergeometric function U(a, b, x) of the second kind
    """)

add_newdoc("scipy.special", "i0",
    """
    i0(x)

    Modified Bessel function of order 0
    """)

add_newdoc("scipy.special", "i0e",
    """
    i0e(x)

    Exponentially scaled modified Bessel function of order 0.

    Defined as::

        i0e(x) = exp(-abs(x)) * i0(x).
    """)

add_newdoc("scipy.special", "i1",
    """
    i1(x)

    Modified Bessel function of order 1
    """)

add_newdoc("scipy.special", "i1e",
    """
    i1e(x)

    Exponentially scaled modified Bessel function of order 1.

    Defined as::

        i1e(x) = exp(-abs(x)) * i1(x)
    """)

add_newdoc("scipy.special", "it2i0k0",
    """
    it2i0k0(x)

    Integrals related to modified Bessel functions of order 0

    Returns
    -------
    ii0
        ``integral((i0(t)-1)/t, t=0..x)``
    ik0
        ``int(k0(t)/t, t=x..inf)``
    """)

add_newdoc("scipy.special", "it2j0y0",
    """
    it2j0y0(x)

    Integrals related to Bessel functions of order 0

    Returns
    -------
    ij0
        ``integral((1-j0(t))/t, t=0..x)``
    iy0
        ``integral(y0(t)/t, t=x..inf)``
    """)

add_newdoc("scipy.special", "it2struve0",
    """
    it2struve0(x)

    Integral related to Struve function of order 0

    Returns
    -------
    i
        ``integral(H0(t)/t, t=x..inf)``
    """)

add_newdoc("scipy.special", "itairy",
    """
    itairy(x)

    Integrals of Airy functions

    Calculates the integral of Airy functions from 0 to `x`

    Returns
    -------
    Apt, Bpt
        Integrals for positive arguments
    Ant, Bnt
        Integrals for negative arguments

    """)

add_newdoc("scipy.special", "iti0k0",
    """
    iti0k0(x)

    Integrals of modified Bessel functions of order 0

    Returns simple integrals from 0 to `x` of the zeroth order modified
    Bessel functions `i0` and `k0`.

    Returns
    -------
    ii0, ik0
    """)

add_newdoc("scipy.special", "itj0y0",
    """
    itj0y0(x)

    Integrals of Bessel functions of order 0

    Returns simple integrals from 0 to `x` of the zeroth order Bessel
    functions `j0` and `y0`.

    Returns
    -------
    ij0, iy0
    """)

add_newdoc("scipy.special", "itmodstruve0",
    """
    itmodstruve0(x)

    Integral of the modified Struve function of order 0

    Returns
    -------
    i
        ``integral(L0(t), t=0..x)``
    """)

add_newdoc("scipy.special", "itstruve0",
    """
    itstruve0(x)

    Integral of the Struve function of order 0

    Returns
    -------
    i
        ``integral(H0(t), t=0..x)``
    """)

add_newdoc("scipy.special", "iv",
    """
    iv(v, z)

    Modified Bessel function of the first kind  of real order

    Parameters
    ----------
    v
        Order. If `z` is of real type and negative, `v` must be integer valued.
    z
        Argument.

    """)

add_newdoc("scipy.special", "ive",
    """
    ive(v, z)

    Exponentially scaled modified Bessel function of the first kind

    Defined as::

        ive(v, z) = iv(v, z) * exp(-abs(z.real))
    """)

add_newdoc("scipy.special", "j0",
    """
    j0(x)

    Bessel function the first kind of order 0
    """)

add_newdoc("scipy.special", "j1",
    """
    j1(x)

    Bessel function of the first kind of order 1
    """)

add_newdoc("scipy.special", "jn",
    """
    jn(n, x)

    Bessel function of the first kind of integer order `n`.

    Notes
    -----
    `jn` is an alias of `jv`.

    """)

add_newdoc("scipy.special", "jv",
    """
    jv(v, z)

    Bessel function of the first kind of real order `v`
    """)

add_newdoc("scipy.special", "jve",
    """
    jve(v, z)

    Exponentially scaled Bessel function of order `v`

    Defined as::

        jve(v, z) = jv(v, z) * exp(-abs(z.imag))
    """)

add_newdoc("scipy.special", "k0",
    """
    k0(x)

    Modified Bessel function K of order 0

    Modified Bessel function of the second kind (sometimes called the
    third kind) of order 0.
    """)

add_newdoc("scipy.special", "k0e",
    """
    k0e(x)

    Exponentially scaled modified Bessel function K of order 0

    Defined as::

        k0e(x) = exp(x) * k0(x).
    """)

add_newdoc("scipy.special", "k1",
    """
    i1(x)

    Modified Bessel function of the first kind of order 1
    """)

add_newdoc("scipy.special", "k1e",
    """
    k1e(x)

    Exponentially scaled modified Bessel function K of order 1

    Defined as::

        k1e(x) = exp(x) * k1(x)
    """)

add_newdoc("scipy.special", "kei",
    """
    kei(x)

    Kelvin function ker
    """)

add_newdoc("scipy.special", "keip",
    """
    keip(x)

    Derivative of the Kelvin function kei
    """)

add_newdoc("scipy.special", "kelvin",
    """
    kelvin(x)

    Kelvin functions as complex numbers

    Returns
    -------
    Be, Ke, Bep, Kep
        The tuple (Be, Ke, Bep, Kep) contains complex numbers
        representing the real and imaginary Kelvin functions and their
        derivatives evaluated at `x`.  For example, kelvin(x)[0].real =
        ber x and kelvin(x)[0].imag = bei x with similar relationships
        for ker and kei.
    """)

add_newdoc("scipy.special", "ker",
    """
    ker(x)

    Kelvin function ker
    """)

add_newdoc("scipy.special", "kerp",
    """
    kerp(x)

    Derivative of the Kelvin function ker
    """)

add_newdoc("scipy.special", "kl_div",
    r"""
    kl_div(x, y)

    Elementwise function for computing Kullback-Leibler divergence.

    .. math:: \mathrm{kl\_div}(x, y) = \begin{cases} x \log(x / y) - x + y & x > 0, y > 0 \\ y & x = 0, y \ge 0 \\ \infty & \text{otherwise} \end{cases}

    Parameters
    ----------
    x : ndarray
        First input array.
    y : ndarray
        Second input array.

    Returns
    -------
    res : ndarray
        Output array.

    See Also
    --------
    entr, rel_entr

    Notes
    -----
    This function is non-negative and is jointly convex in `x` and `y`.

    .. versionadded:: 0.14.0

    """)

add_newdoc("scipy.special", "kn",
    """
    kn(n, x)

    Modified Bessel function of the second kind of integer order `n`

    Returns the modified Bessel function of the second kind for integer order
    `n` at real `z`.

    These are also sometimes called functions of the third kind, Basset
    functions, or Macdonald functions.

    Parameters
    ----------
    n : array_like of int
        Order of Bessel functions (floats will truncate with a warning)
    z : array_like of float
        Argument at which to evaluate the Bessel functions

    Returns
    -------
    out : ndarray
        The results

    Examples
    --------
    Plot the function of several orders for real input:

    >>> from scipy.special import kn
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(0, 5, 1000)
    >>> for N in range(6):
    ...     plt.plot(x, kn(N, x), label='$K_{}(x)$'.format(N))
    >>> plt.ylim(0, 10)
    >>> plt.legend()
    >>> plt.title(r'Modified Bessel function of the second kind $K_n(x)$')
    >>> plt.show()

    Calculate for a single value at multiple orders:

    >>> print(kn([4, 5, 6], 1))
    [   44.23241425   360.96060181  3653.83837891]
    """)

add_newdoc("scipy.special", "kolmogi",
    """
    kolmogi(p)

    Inverse function to kolmogorov

    Returns y such that ``kolmogorov(y) == p``.
    """)

add_newdoc("scipy.special", "kolmogorov",
    """
    kolmogorov(y)

    Complementary cumulative distribution function of Kolmogorov distribution

    Returns the complementary cumulative distribution function of
    Kolmogorov's limiting distribution (Kn* for large n) of a
    two-sided test for equality between an empirical and a theoretical
    distribution. It is equal to the (limit as n->infinity of the)
    probability that sqrt(n) * max absolute deviation > y.
    """)

add_newdoc("scipy.special", "kv",
    """
    kv(v, z)

    Modified Bessel function of the second kind of real order `v`

    Returns the modified Bessel function of the second kind for real order
    `v` at complex `z`.

    These are also sometimes called functions of the third kind, Basset
    functions, or Macdonald functions.

    Parameters
    ----------
    v : array_like of float
        Order of Bessel functions
    z : array_like of complex
        Argument at which to evaluate the Bessel functions

    Returns
    -------
    out : ndarray
        The results

    Examples
    --------
    Plot the function of several orders for real input:

    >>> from scipy.special import kv
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(0, 5, 1000)
    >>> for N in np.linspace(0, 6, 5):
    ...     plt.plot(x, kv(N, x), label='$K_{{{}}}(x)$'.format(N))
    >>> plt.ylim(0, 10)
    >>> plt.legend()
    >>> plt.title(r'Modified Bessel function of the second kind $K_\nu(x)$')
    >>> plt.show()

    Calculate for a single value at multiple orders:

    >>> print(kv([4, 4.5, 5], 1+2j))
    [ 0.19920848+2.38918114j  2.34926017+3.59995073j  7.28267680+3.81037734j]

    """)

add_newdoc("scipy.special", "kve",
    """
    kve(v, z)

    Exponentially scaled modified Bessel function of the second kind.

    Returns the exponentially scaled, modified Bessel function of the
    second kind (sometimes called the third kind) for real order `v` at
    complex `z`::

        kve(v, z) = kv(v, z) * exp(z)
    """)

add_newdoc("scipy.special", "log1p",
    """
    log1p(x)

    Calculates log(1+x) for use when `x` is near zero
    """)

add_newdoc('scipy.special', 'logit',
    """
    logit(x)

    Logit ufunc for ndarrays.

    The logit function is defined as logit(p) = log(p/(1-p)).
    Note that logit(0) = -inf, logit(1) = inf, and logit(p)
    for p<0 or p>1 yields nan.

    Parameters
    ----------
    x : ndarray
        The ndarray to apply logit to element-wise.

    Returns
    -------
    out : ndarray
        An ndarray of the same shape as x. Its entries
        are logit of the corresponding entry of x.

    Notes
    -----
    As a ufunc logit takes a number of optional
    keyword arguments. For more information
    see `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_

    .. versionadded:: 0.10.0

    """)

add_newdoc("scipy.special", "lpmv",
    """
    lpmv(m, v, x)

    Associated legendre function of integer order.

    Parameters
    ----------
    m : int
        Order
    v : real
        Degree. Must be ``v>-m-1`` or ``v<m``
    x : complex
        Argument. Must be ``|x| <= 1``.

    """)

add_newdoc("scipy.special", "mathieu_a",
    """
    mathieu_a(m, q)

    Characteristic value of even Mathieu functions

    Returns the characteristic value for the even solution,
    ``ce_m(z, q)``, of Mathieu's equation.
    """)

add_newdoc("scipy.special", "mathieu_b",
    """
    mathieu_b(m, q)

    Characteristic value of odd Mathieu functions

    Returns the characteristic value for the odd solution,
    ``se_m(z, q)``, of Mathieu's equation.
    """)

add_newdoc("scipy.special", "mathieu_cem",
    """
    mathieu_cem(m, q, x)

    Even Mathieu function and its derivative

    Returns the even Mathieu function, ``ce_m(x, q)``, of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` of ce_m(x, q)

    Parameters
    ----------
    m
        Order of the function
    q
        Parameter of the function
    x
        Argument of the function, *given in degrees, not radians*

    Returns
    -------
    y
        Value of the function
    yp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "mathieu_modcem1",
    """
    mathieu_modcem1(m, q, x)

    Even modified Mathieu function of the first kind and its derivative

    Evaluates the even modified Mathieu function of the first kind,
    ``Mc1m(x, q)``, and its derivative at `x` for order `m` and parameter
    `q`.

    Returns
    -------
    y
        Value of the function
    yp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "mathieu_modcem2",
    """
    mathieu_modcem2(m, q, x)

    Even modified Mathieu function of the second kind and its derivative

    Evaluates the even modified Mathieu function of the second kind,
    Mc2m(x, q), and its derivative at `x` (given in degrees) for order `m`
    and parameter `q`.

    Returns
    -------
    y
        Value of the function
    yp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "mathieu_modsem1",
    """
    mathieu_modsem1(m, q, x)

    Odd modified Mathieu function of the first kind and its derivative

    Evaluates the odd modified Mathieu function of the first kind,
    Ms1m(x, q), and its derivative at `x` (given in degrees) for order `m`
    and parameter `q`.

    Returns
    -------
    y
        Value of the function
    yp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "mathieu_modsem2",
    """
    mathieu_modsem2(m, q, x)

    Odd modified Mathieu function of the second kind and its derivative

    Evaluates the odd modified Mathieu function of the second kind,
    Ms2m(x, q), and its derivative at `x` (given in degrees) for order `m`
    and parameter q.

    Returns
    -------
    y
        Value of the function
    yp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "mathieu_sem",
    """
    mathieu_sem(m, q, x)

    Odd Mathieu function and its derivative

    Returns the odd Mathieu function, se_m(x, q), of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` of se_m(x, q).

    Parameters
    ----------
    m
        Order of the function
    q
        Parameter of the function
    x
        Argument of the function, *given in degrees, not radians*.

    Returns
    -------
    y
        Value of the function
    yp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "modfresnelm",
    """
    modfresnelm(x)

    Modified Fresnel negative integrals

    Returns
    -------
    fm
        Integral ``F_-(x)``: ``integral(exp(-1j*t*t), t=x..inf)``
    km
        Integral ``K_-(x)``: ``1/sqrt(pi)*exp(1j*(x*x+pi/4))*fp``
    """)

add_newdoc("scipy.special", "modfresnelp",
    """
    modfresnelp(x)

    Modified Fresnel positive integrals

    Returns
    -------
    fp
        Integral ``F_+(x)``: ``integral(exp(1j*t*t), t=x..inf)``
    kp
        Integral ``K_+(x)``: ``1/sqrt(pi)*exp(-1j*(x*x+pi/4))*fp``
    """)

add_newdoc("scipy.special", "modstruve",
    """
    modstruve(v, x)

    Modified Struve function

    Returns the modified Struve function Lv(x) of order `v` at `x`, `x` must
    be positive unless `v` is an integer.
    """)

add_newdoc("scipy.special", "nbdtr",
    """
    nbdtr(k, n, p)

    Negative binomial cumulative distribution function

    Returns the sum of the terms 0 through `k` of the negative binomial
    distribution::

        sum((n+j-1)Cj p**n (1-p)**j, j=0..k).

    In a sequence of Bernoulli trials this is the probability that k
    or fewer failures precede the nth success.
    """)

add_newdoc("scipy.special", "nbdtrc",
    """
    nbdtrc(k, n, p)

    Negative binomial survival function

    Returns the sum of the terms k+1 to infinity of the negative
    binomial distribution.
    """)

add_newdoc("scipy.special", "nbdtri",
    """
    nbdtri(k, n, y)

    Inverse of `nbdtr` vs `p`

    Finds the argument p such that ``nbdtr(k, n, p) = y``.
    """)

add_newdoc("scipy.special", "nbdtrik",
    """
    nbdtrik(y, n, p)

    Inverse of `nbdtr` vs `k`

    Finds the argument k such that ``nbdtr(k, n, p) = y``.
    """)

add_newdoc("scipy.special", "nbdtrin",
    """
    nbdtrin(k, y, p)

    Inverse of `nbdtr` vs `n`

    Finds the argument `n` such that ``nbdtr(k, n, p) = y``.
    """)

add_newdoc("scipy.special", "ncfdtr",
    """
    ncfdtr(dfn, dfd, nc, f)

    Cumulative distribution function of the non-central F distribution.

    Parameters
    ----------
    dfn : array_like
        Degrees of freedom of the numerator sum of squares.  Range (0, inf).
    dfd : array_like
        Degrees of freedom of the denominator sum of squares.  Range (0, inf).
    nc : array_like
        Noncentrality parameter.  Should be in range (0, 1e4).
    f : array_like
        Quantiles, i.e. the upper limit of integration.

    Returns
    -------
    cdf : float or ndarray
        The calculated CDF.  If all inputs are scalar, the return will be a
        float.  Otherwise it will be an array.

    See Also
    --------
    ncdfdtri : Inverse CDF (iCDF) of the non-central F distribution.
    ncdfdtridfd : Calculate dfd, given CDF and iCDF values.
    ncdfdtridfn : Calculate dfn, given CDF and iCDF values.
    ncdfdtrinc : Calculate noncentrality parameter, given CDF, iCDF, dfn, dfd.

    Examples
    --------
    >>> from scipy import special
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    Plot the CDF of the non-central F distribution, for nc=0.  Compare with the
    F-distribution from scipy.stats:

    >>> x = np.linspace(-1, 8, num=500)
    >>> dfn = 3
    >>> dfd = 2
    >>> ncf_stats = stats.f.cdf(x, dfn, dfd)
    >>> ncf_special = special.ncfdtr(dfn, dfd, 0, x)

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, ncf_stats, 'b-', lw=3)
    >>> ax.plot(x, ncf_special, 'r-')
    >>> plt.show()

    """)

add_newdoc("scipy.special", "ncfdtri",
    """
    ncfdtri(p, dfn, dfd, nc)

    Inverse cumulative distribution function of the non-central F distribution.

    See `ncfdtr` for more details.

    """)

add_newdoc("scipy.special", "ncfdtridfd",
    """
    ncfdtridfd(p, f, dfn, nc)

    Calculate degrees of freedom (denominator) for the noncentral F-distribution.

    See `ncfdtr` for more details.

    """)

add_newdoc("scipy.special", "ncfdtridfn",
    """
    ncfdtridfn(p, f, dfd, nc)

    Calculate degrees of freedom (numerator) for the noncentral F-distribution.

    See `ncfdtr` for more details.

    """)

add_newdoc("scipy.special", "ncfdtrinc",
    """
    ncfdtrinc(p, f, dfn, dfd)

    Calculate non-centrality parameter for non-central F distribution.

    See `ncfdtr` for more details.

    """)

add_newdoc("scipy.special", "nctdtr",
    """
    nctdtr(df, nc, t)

    Cumulative distribution function of the non-central `t` distribution.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution.  Should be in range (0, inf).
    nc : array_like
        Noncentrality parameter.  Should be in range (-1e6, 1e6).
    t : array_like
        Quantiles, i.e. the upper limit of integration.

    Returns
    -------
    cdf : float or ndarray
        The calculated CDF.  If all inputs are scalar, the return will be a
        float.  Otherwise it will be an array.

    See Also
    --------
    nctdtrit : Inverse CDF (iCDF) of the non-central t distribution.
    nctdtridf : Calculate degrees of freedom, given CDF and iCDF values.
    nctdtrinc : Calculate non-centrality parameter, given CDF iCDF values.

    Examples
    --------
    >>> from scipy import special
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    Plot the CDF of the non-central t distribution, for nc=0.  Compare with the
    t-distribution from scipy.stats:

    >>> x = np.linspace(-5, 5, num=500)
    >>> df = 3
    >>> nct_stats = stats.t.cdf(x, df)
    >>> nct_special = special.nctdtr(df, 0, x)

    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.plot(x, nct_stats, 'b-', lw=3)
    >>> ax.plot(x, nct_special, 'r-')
    >>> plt.show()

    """)

add_newdoc("scipy.special", "nctdtridf",
    """
    nctdtridf(p, nc, t)

    Calculate degrees of freedom for non-central t distribution.

    See `nctdtr` for more details.

    Parameters
    ----------
    p : array_like
        CDF values, in range (0, 1].
    nc : array_like
        Noncentrality parameter.  Should be in range (-1e6, 1e6).
    t : array_like
        Quantiles, i.e. the upper limit of integration.

    """)

add_newdoc("scipy.special", "nctdtrinc",
    """
    nctdtrinc(df, p, t)

    Calculate non-centrality parameter for non-central t distribution.

    See `nctdtr` for more details.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution.  Should be in range (0, inf).
    p : array_like
        CDF values, in range (0, 1].
    t : array_like
        Quantiles, i.e. the upper limit of integration.

    """)

add_newdoc("scipy.special", "nctdtrit",
    """
    nctdtrit(df, nc, p)

    Inverse cumulative distribution function of the non-central t distribution.

    See `nctdtr` for more details.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution.  Should be in range (0, inf).
    nc : array_like
        Noncentrality parameter.  Should be in range (-1e6, 1e6).
    p : array_like
        CDF values, in range (0, 1].

    """)

add_newdoc("scipy.special", "ndtr",
    """
    ndtr(x)

    Gaussian cumulative distribution function

    Returns the area under the standard Gaussian probability
    density function, integrated from minus infinity to `x`::

        1/sqrt(2*pi) * integral(exp(-t**2 / 2), t=-inf..x)

    """)

add_newdoc("scipy.special", "nrdtrimn",
    """
    nrdtrimn(p, x, std)

    Calculate mean of normal distribution given other params.

    Parameters
    ----------
    p : array_like
        CDF values, in range (0, 1].
    x : array_like
        Quantiles, i.e. the upper limit of integration.
    std : array_like
        Standard deviation.

    Returns
    -------
    mn : float or ndarray
        The mean of the normal distribution.

    See Also
    --------
    nrdtrimn, ndtr

    """)

add_newdoc("scipy.special", "nrdtrisd",
    """
    nrdtrisd(p, x, mn)

    Calculate standard deviation of normal distribution given other params.

    Parameters
    ----------
    p : array_like
        CDF values, in range (0, 1].
    x : array_like
        Quantiles, i.e. the upper limit of integration.
    mn : float or ndarray
        The mean of the normal distribution.

    Returns
    -------
    std : array_like
        Standard deviation.

    See Also
    --------
    nrdtristd, ndtr

    """)

add_newdoc("scipy.special", "log_ndtr",
    """
    log_ndtr(x)

    Logarithm of Gaussian cumulative distribution function

    Returns the log of the area under the standard Gaussian probability
    density function, integrated from minus infinity to `x`::

        log(1/sqrt(2*pi) * integral(exp(-t**2 / 2), t=-inf..x))
    """)

add_newdoc("scipy.special", "ndtri",
    """
    ndtri(y)

    Inverse of `ndtr` vs x

    Returns the argument x for which the area under the Gaussian
    probability density function (integrated from minus infinity to `x`)
    is equal to y.
    """)

add_newdoc("scipy.special", "obl_ang1",
    """
    obl_ang1(m, n, c, x)

    Oblate spheroidal angular function of the first kind and its derivative

    Computes the oblate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "obl_ang1_cv",
    """
    obl_ang1_cv(m, n, c, cv, x)

    Oblate spheroidal angular function obl_ang1 for precomputed characteristic value

    Computes the oblate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "obl_cv",
    """
    obl_cv(m, n, c)

    Characteristic value of oblate spheroidal function

    Computes the characteristic value of oblate spheroidal wave
    functions of order `m`, `n` (n>=m) and spheroidal parameter `c`.
    """)

add_newdoc("scipy.special", "obl_rad1",
    """
    obl_rad1(m, n, c, x)

    Oblate spheroidal radial function of the first kind and its derivative

    Computes the oblate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "obl_rad1_cv",
    """
    obl_rad1_cv(m, n, c, cv, x)

    Oblate spheroidal radial function obl_rad1 for precomputed characteristic value

    Computes the oblate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "obl_rad2",
    """
    obl_rad2(m, n, c, x)

    Oblate spheroidal radial function of the second kind and its derivative.

    Computes the oblate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "obl_rad2_cv",
    """
    obl_rad2_cv(m, n, c, cv, x)

    Oblate spheroidal radial function obl_rad2 for precomputed characteristic value

    Computes the oblate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pbdv",
    """
    pbdv(v, x)

    Parabolic cylinder function D

    Returns (d, dp) the parabolic cylinder function Dv(x) in d and the
    derivative, Dv'(x) in dp.

    Returns
    -------
    d
        Value of the function
    dp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pbvv",
    """
    pbvv(v, x)

    Parabolic cylinder function V

    Returns the parabolic cylinder function Vv(x) in v and the
    derivative, Vv'(x) in vp.

    Returns
    -------
    v
        Value of the function
    vp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pbwa",
    """
    pbwa(a, x)

    Parabolic cylinder function W

    Returns the parabolic cylinder function W(a, x) in w and the
    derivative, W'(a, x) in wp.

    .. warning::

       May not be accurate for large (>5) arguments in a and/or x.

    Returns
    -------
    w
        Value of the function
    wp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pdtr",
    """
    pdtr(k, m)

    Poisson cumulative distribution function

    Returns the sum of the first `k` terms of the Poisson distribution:
    sum(exp(-m) * m**j / j!, j=0..k) = gammaincc( k+1, m).  Arguments
    must both be positive and `k` an integer.
    """)

add_newdoc("scipy.special", "pdtrc",
    """
    pdtrc(k, m)

    Poisson survival function

    Returns the sum of the terms from k+1 to infinity of the Poisson
    distribution: sum(exp(-m) * m**j / j!, j=k+1..inf) = gammainc(
    k+1, m).  Arguments must both be positive and `k` an integer.
    """)

add_newdoc("scipy.special", "pdtri",
    """
    pdtri(k, y)

    Inverse to `pdtr` vs m

    Returns the Poisson variable `m` such that the sum from 0 to `k` of
    the Poisson density is equal to the given probability `y`:
    calculated by gammaincinv(k+1, y). `k` must be a nonnegative
    integer and `y` between 0 and 1.
    """)

add_newdoc("scipy.special", "pdtrik",
    """
    pdtrik(p, m)

    Inverse to `pdtr` vs k

    Returns the quantile k such that ``pdtr(k, m) = p``
    """)

add_newdoc("scipy.special", "poch",
    """
    poch(z, m)

    Rising factorial (z)_m

    The Pochhammer symbol (rising factorial), is defined as::

        (z)_m = gamma(z + m) / gamma(z)

    For positive integer `m` it reads::

        (z)_m = z * (z + 1) * ... * (z + m - 1)
    """)

add_newdoc("scipy.special", "pro_ang1",
    """
    pro_ang1(m, n, c, x)

    Prolate spheroidal angular function of the first kind and its derivative

    Computes the prolate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pro_ang1_cv",
    """
    pro_ang1_cv(m, n, c, cv, x)

    Prolate spheroidal angular function pro_ang1 for precomputed characteristic value

    Computes the prolate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pro_cv",
    """
    pro_cv(m, n, c)

    Characteristic value of prolate spheroidal function

    Computes the characteristic value of prolate spheroidal wave
    functions of order `m`, `n` (n>=m) and spheroidal parameter `c`.
    """)

add_newdoc("scipy.special", "pro_rad1",
    """
    pro_rad1(m, n, c, x)

    Prolate spheroidal radial function of the first kind and its derivative

    Computes the prolate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pro_rad1_cv",
    """
    pro_rad1_cv(m, n, c, cv, x)

    Prolate spheroidal radial function pro_rad1 for precomputed characteristic value

    Computes the prolate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pro_rad2",
    """
    pro_rad2(m, n, c, x)

    Prolate spheroidal radial function of the secon kind and its derivative

    Computes the prolate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pro_rad2_cv",
    """
    pro_rad2_cv(m, n, c, cv, x)

    Prolate spheroidal radial function pro_rad2 for precomputed characteristic value

    Computes the prolate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Returns
    -------
    s
        Value of the function
    sp
        Value of the derivative vs x
    """)

add_newdoc("scipy.special", "pseudo_huber",
    r"""
    pseudo_huber(delta, r)

    Pseudo-Huber loss function.

    .. math:: \mathrm{pseudo\_huber}(\delta, r) = \delta^2 \left( \sqrt{ 1 + \left( \frac{r}{\delta} \right)^2 } - 1 \right)

    Parameters
    ----------
    delta : ndarray
        Input array, indicating the soft quadratic vs. linear loss changepoint.
    r : ndarray
        Input array, possibly representing residuals.

    Returns
    -------
    res : ndarray
        The computed Pseudo-Huber loss function values.

    Notes
    -----
    This function is convex in :math:`r`.

    .. versionadded:: 0.15.0

    """)

add_newdoc("scipy.special", "psi",
    """
    psi(z)

    Digamma function

    The derivative of the logarithm of the gamma function evaluated at
    `z` (also called the digamma function).
    """)

add_newdoc("scipy.special", "radian",
    """
    radian(d, m, s)

    Convert from degrees to radians

    Returns the angle given in (d)egrees, (m)inutes, and (s)econds in
    radians.
    """)

add_newdoc("scipy.special", "rel_entr",
    r"""
    rel_entr(x, y)

    Elementwise function for computing relative entropy.

    .. math:: \mathrm{rel\_entr}(x, y) = \begin{cases} x \log(x / y) & x > 0, y > 0 \\ 0 & x = 0, y \ge 0 \\ \infty & \text{otherwise} \end{cases}

    Parameters
    ----------
    x : ndarray
        First input array.
    y : ndarray
        Second input array.

    Returns
    -------
    res : ndarray
        Output array.

    See Also
    --------
    entr, kl_div

    Notes
    -----
    This function is jointly convex in x and y.

    .. versionadded:: 0.14.0

    """)

add_newdoc("scipy.special", "rgamma",
    """
    rgamma(z)

    Gamma function inverted

    Returns ``1/gamma(x)``
    """)

add_newdoc("scipy.special", "round",
    """
    round(x)

    Round to nearest integer

    Returns the nearest integer to `x` as a double precision floating
    point result.  If `x` ends in 0.5 exactly, the nearest even integer
    is chosen.
    """)

add_newdoc("scipy.special", "shichi",
    """
    shichi(x)

    Hyperbolic sine and cosine integrals

    Returns
    -------
    shi
        ``integral(sinh(t)/t, t=0..x)``
    chi
        ``eul + ln x + integral((cosh(t)-1)/t, t=0..x)``
        where ``eul`` is Euler's constant.
    """)

add_newdoc("scipy.special", "sici",
    """
    sici(x)

    Sine and cosine integrals

    Returns
    -------
    si
        ``integral(sin(t)/t, t=0..x)``
    ci
        ``eul + ln x + integral((cos(t) - 1)/t, t=0..x)``
        where ``eul`` is Euler's constant.
    """)

add_newdoc("scipy.special", "sindg",
    """
    sindg(x)

    Sine of angle given in degrees
    """)

add_newdoc("scipy.special", "smirnov",
    """
    smirnov(n, e)

    Kolmogorov-Smirnov complementary cumulative distribution function

    Returns the exact Kolmogorov-Smirnov complementary cumulative
    distribution function (Dn+ or Dn-) for a one-sided test of
    equality between an empirical and a theoretical distribution. It
    is equal to the probability that the maximum difference between a
    theoretical distribution and an empirical one based on `n` samples
    is greater than e.
    """)

add_newdoc("scipy.special", "smirnovi",
    """
    smirnovi(n, y)

    Inverse to `smirnov`

    Returns ``e`` such that ``smirnov(n, e) = y``.
    """)

add_newdoc("scipy.special", "spence",
    """
    spence(x)

    Dilogarithm integral

    Returns the dilogarithm integral::

        -integral(log t / (t-1), t=1..x)
    """)

add_newdoc("scipy.special", "stdtr",
    """
    stdtr(df, t)

    Student t distribution cumulative density function

    Returns the integral from minus infinity to t of the Student t
    distribution with df > 0 degrees of freedom::

       gamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2)) *
       integral((1+x**2/df)**(-df/2-1/2), x=-inf..t)

    """)

add_newdoc("scipy.special", "stdtridf",
    """
    stdtridf(p, t)

    Inverse of `stdtr` vs df

    Returns the argument df such that stdtr(df, t) is equal to `p`.
    """)

add_newdoc("scipy.special", "stdtrit",
    """
    stdtrit(df, p)

    Inverse of `stdtr` vs `t`

    Returns the argument `t` such that stdtr(df, t) is equal to `p`.
    """)

add_newdoc("scipy.special", "struve",
    """
    struve(v, x)

    Struve function

    Computes the struve function Hv(x) of order `v` at `x`, `x` must be
    positive unless `v` is an integer.
    """)

add_newdoc("scipy.special", "tandg",
    """
    tandg(x)

    Tangent of angle x given in degrees.
    """)

add_newdoc("scipy.special", "tklmbda",
    """
    tklmbda(x, lmbda)

    Tukey-Lambda cumulative distribution function

    """)

add_newdoc("scipy.special", "wofz",
    """
    wofz(z)

    Faddeeva function

    Returns the value of the Faddeeva function for complex argument::

        exp(-z**2)*erfc(-i*z)

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva
    """)

add_newdoc("scipy.special", "xlogy",
    """
    xlogy(x, y)

    Compute ``x*log(y)`` so that the result is 0 if ``x = 0``.

    Parameters
    ----------
    x : array_like
        Multiplier
    y : array_like
        Argument

    Returns
    -------
    z : array_like
        Computed x*log(y)

    Notes
    -----

    .. versionadded:: 0.13.0

    """)

add_newdoc("scipy.special", "xlog1py",
    """
    xlog1py(x, y)

    Compute ``x*log1p(y)`` so that the result is 0 if ``x = 0``.

    Parameters
    ----------
    x : array_like
        Multiplier
    y : array_like
        Argument

    Returns
    -------
    z : array_like
        Computed x*log1p(y)

    Notes
    -----

    .. versionadded:: 0.13.0

    """)

add_newdoc("scipy.special", "y0",
    """
    y0(x)

    Bessel function of the second kind of order 0

    Returns the Bessel function of the second kind of order 0 at `x`.
    """)

add_newdoc("scipy.special", "y1",
    """
    y1(x)

    Bessel function of the second kind of order 1

    Returns the Bessel function of the second kind of order 1 at `x`.
    """)

add_newdoc("scipy.special", "yn",
    """
    yn(n, x)

    Bessel function of the second kind of integer order

    Returns the Bessel function of the second kind of integer order `n`
    at `x`.
    """)

add_newdoc("scipy.special", "yv",
    """
    yv(v, z)

    Bessel function of the second kind of real order

    Returns the Bessel function of the second kind of real order `v` at
    complex `z`.
    """)

add_newdoc("scipy.special", "yve",
    """
    yve(v, z)

    Exponentially scaled Bessel function of the second kind of real order

    Returns the exponentially scaled Bessel function of the second
    kind of real order `v` at complex `z`::

        yve(v, z) = yv(v, z) * exp(-abs(z.imag))

    """)

add_newdoc("scipy.special", "zeta",
    """
    zeta(x, q)

    Hurwitz zeta function

    The Riemann zeta function of two arguments (also known as the
    Hurwitz zeta function).

    This function is defined as

    .. math:: \\zeta(x, q) = \\sum_{k=0}^{\\infty} 1 / (k+q)^x,

    where ``x > 1`` and ``q > 0``.

    See also
    --------
    zetac

    """)

add_newdoc("scipy.special", "zetac",
    """
    zetac(x)

    Riemann zeta function minus 1.

    This function is defined as

    .. math:: \\zeta(x) = \\sum_{k=2}^{\\infty} 1 / k^x,

    where ``x > 1``.

    See Also
    --------
    zeta

    """)

add_newdoc("scipy.special", "_struve_asymp_large_z",
    """
    _struve_asymp_large_z(v, z, is_h)

    Internal function for testing `struve` & `modstruve`

    Evaluates using asymptotic expansion

    Returns
    -------
    v, err
    """)

add_newdoc("scipy.special", "_struve_power_series",
    """
    _struve_power_series(v, z, is_h)

    Internal function for testing `struve` & `modstruve`

    Evaluates using power series

    Returns
    -------
    v, err
    """)

add_newdoc("scipy.special", "_struve_bessel_series",
    """
    _struve_bessel_series(v, z, is_h)

    Internal function for testing `struve` & `modstruve`

    Evaluates using Bessel function series

    Returns
    -------
    v, err
    """)
