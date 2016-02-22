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
    r"""
    airy(z)

    Airy functions and their derivatives.

    Parameters
    ----------
    z : array_like
        Real or complex argument.

    Returns
    -------
    Ai, Aip, Bi, Bip : ndarrays
        Airy functions Ai and Bi, and their derivatives Aip and Bip.

    Notes
    -----
    The Airy functions Ai and Bi are two independent solutions of

    .. math:: y''(x) = x y(x).

    For real `z` in [-10, 10], the computation is carried out by calling
    the Cephes [1]_ `airy` routine, which uses power series summation
    for small `z` and rational minimax approximations for large `z`.

    Outside this range, the AMOS [2]_ `zairy` and `zbiry` routines are
    employed.  They are computed using power series for :math:`|z| < 1` and
    the following relations to modified Bessel functions for larger `z`
    (where :math:`t \equiv 2 z^{3/2}/3`):

    .. math::

        Ai(z) = \frac{1}{\pi \sqrt{3}} K_{1/3}(t)

        Ai'(z) = -\frac{z}{\pi \sqrt{3}} K_{2/3}(t)

        Bi(z) = \sqrt{\frac{z}{3}} \left(I_{-1/3}(t) + I_{1/3}(t) \right)

        Bi'(z) = \frac{z}{\sqrt{3}} \left(I_{-2/3}(t) + I_{2/3}(t)\right)

    See also
    --------
    airye : exponentially scaled Airy functions.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    .. [2] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/.org/amos/
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
    z : array_like
        Real or complex argument.

    Returns
    -------
    eAi, eAip, eBi, eBip : array_like
        Airy functions Ai and Bi, and their derivatives Aip and Bip

    Notes
    -----
    Wrapper for the AMOS [1]_ routines `zairy` and `zbiry`.

    See also
    --------
    airy

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "bdtr",
    r"""
    bdtr(k, n, p)

    Binomial distribution cumulative distribution function.

    Sum of the terms 0 through `k` of the Binomial probability density.

    .. math::
        \mathrm{bdtr}(k, n, p) = \sum_{j=0}^k {{n}\choose{j}} p^j (1-p)^{n-j}

    Parameters
    ----------
    k : array_like
        Number of successes (int).
    n : array_like
        Number of events (int).
    p : array_like
        Probability of success in a single event (float).

    Returns
    -------
    y : ndarray
        Probability of `k` or fewer successes in `n` independent events with
        success probabilities of `p`.

    Notes
    -----
    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{bdtr}(k, n, p) = I_{1 - p}(n - k, k + 1).

    Wrapper for the Cephes [1]_ routine `bdtr`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

    """)

add_newdoc("scipy.special", "bdtrc",
    r"""
    bdtrc(k, n, p)

    Binomial distribution survival function.

    Sum of the terms `k + 1` through `n` of the binomial probability density,

    .. math::
        \mathrm{bdtrc}(k, n, p) = \sum_{j=k+1}^n {{n}\choose{j}} p^j (1-p)^{n-j}

    Parameters
    ----------
    k : array_like
        Number of successes (int).
    n : array_like
        Number of events (int)
    p : array_like
        Probability of success in a single event.

    Returns
    -------
    y : ndarray
        Probability of `k + 1` or more successes in `n` independent events
        with success probabilities of `p`.

    See also
    --------
    bdtr
    betainc

    Notes
    -----
    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{bdtrc}(k, n, p) = I_{p}(k + 1, n - k).

    Wrapper for the Cephes [1]_ routine `bdtrc`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

    """)

add_newdoc("scipy.special", "bdtri",
    """
    bdtri(k, n, y)

    Inverse function to `bdtr` with respect to `p`.

    Finds the event probability `p` such that the sum of the terms 0 through
    `k` of the binomial probability density is equal to the given cumulative
    probability `y`.

    Parameters
    ----------
    k : array_like
        Number of successes (float).
    n : array_like
        Number of events (float)
    y : array_like
        Cumulative probability (probability of `k` or fewer successes in `n`
        events).

    Returns
    -------
    p : ndarray
        The event probability such that `bdtr(k, n, p) = y`.

    See also
    --------
    bdtr
    betaincinv

    Notes
    -----
    The computation is carried out using the inverse beta integral function
    and the relation,::

        1 - p = betaincinv(n - k, k + 1, y).

    Wrapper for the Cephes [1]_ routine `bdtri`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "bdtrik",
    """
    bdtrik(y, n, p)

    Inverse function to `bdtr` with respect to `k`.

    Finds the number of successes `k` such that the sum of the terms 0 through
    `k` of the Binomial probability density for `n` events with probability
    `p` is equal to the given cumulative probability `y`.

    Parameters
    ----------
    y : array_like
        Cumulative probability (probability of `k` or fewer successes in `n`
        events).
    n : array_like
        Number of events (float).
    p : array_like
        Success probability (float).

    Returns
    -------
    k : ndarray
        The number of successes `k` such that `bdtr(k, n, p) = y`.

    See also
    --------
    bdtr

    Notes
    -----
    Formula 26.5.24 of [1]_ is used to reduce the binomial distribution to the
    cumulative incomplete beta distribution.

    Computation of `k` involves a seach for a value that produces the desired
    value of `y`.  The search relies on the monotinicity of `y` with `k`.

    Wrapper for the CDFLIB [2]_ Fortran routine `cdfbin`.

    References
    ----------
    .. [1] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [2] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.

    """)

add_newdoc("scipy.special", "bdtrin",
    """
    bdtrin(k, y, p)

    Inverse function to `bdtr` with respect to `n`.

    Finds the number of events `n` such that the sum of the terms 0 through
    `k` of the Binomial probability density for events with probability `p` is
    equal to the given cumulative probability `y`.

    Parameters
    ----------
    k : array_like
        Number of successes (float).
    y : array_like
        Cumulative probability (probability of `k` or fewer successes in `n`
        events).
    p : array_like
        Success probability (float).

    Returns
    -------
    n : ndarray
        The number of events `n` such that `bdtr(k, n, p) = y`.

    See also
    --------
    bdtr

    Notes
    -----
    Formula 26.5.24 of [1]_ is used to reduce the binomial distribution to the
    cumulative incomplete beta distribution.

    Computation of `n` involves a seach for a value that produces the desired
    value of `y`.  The search relies on the monotinicity of `y` with `n`.

    Wrapper for the CDFLIB [2]_ Fortran routine `cdfbin`.

    References
    ----------
    .. [1] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [2] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    """)

add_newdoc("scipy.special", "binom",
    """
    binom(n, k)

    Binomial coefficient
    """)

add_newdoc("scipy.special", "btdtria",
    r"""
    btdtria(p, b, x)

    Inverse of `btdtr` with respect to `a`.

    This is the inverse of the beta cumulative distribution function, `btdtr`,
    considered as a function of `a`, returning the value of `a` for which
    `btdtr(a, b, x) = p`, or

    .. math::
        p = \int_0^x \frac{\Gamma(a + b)}{\Gamma(a)\Gamma(b)} t^{a-1} (1-t)^{b-1}\,dt

    Parameters
    ----------
    p : array_like
        Cumulative probability, in [0, 1].
    b : array_like
        Shape parameter (`b` > 0).
    x : array_like
        The quantile, in [0, 1].

    Returns
    -------
    a : ndarray
        The value of the shape parameter `a` such that `btdtr(a, b, x) = p`.

    See Also
    --------
    btdtr : Cumulative density function of the beta distribution.
    btdtri : Inverse with respect to `x`.
    btdtrib : Inverse with respect to `b`.

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdfbet`.

    The cumulative distribution function `p` is computed using a routine by
    DiDinato and Morris [2]_.  Computation of `a` involves a seach for a value
    that produces the desired value of `p`.  The search relies on the
    monotinicity of `p` with `a`.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] DiDinato, A. R. and Morris, A. H.,
           Algorithm 708: Significant Digit Computation of the Incomplete Beta
           Function Ratios. ACM Trans. Math. Softw. 18 (1993), 360-373.

    """)

add_newdoc("scipy.special", "btdtrib",
    r"""
    btdtria(a, p, x)

    Inverse of `btdtr` with respect to `b`.

    This is the inverse of the beta cumulative distribution function, `btdtr`,
    considered as a function of `b`, returning the value of `b` for which
    `btdtr(a, b, x) = p`, or

    .. math::
        p = \int_0^x \frac{\Gamma(a + b)}{\Gamma(a)\Gamma(b)} t^{a-1} (1-t)^{b-1}\,dt

    Parameters
    ----------
    a : array_like
        Shape parameter (`a` > 0).
    p : array_like
        Cumulative probability, in [0, 1].
    x : array_like
        The quantile, in [0, 1].

    Returns
    -------
    b : ndarray
        The value of the shape parameter `b` such that `btdtr(a, b, x) = p`.

    See Also
    --------
    btdtr : Cumulative density function of the beta distribution.
    btdtri : Inverse with respect to `x`.
    btdtria : Inverse with respect to `a`.

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdfbet`.

    The cumulative distribution function `p` is computed using a routine by
    DiDinato and Morris [2]_.  Computation of `b` involves a seach for a value
    that produces the desired value of `p`.  The search relies on the
    monotinicity of `p` with `b`.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] DiDinato, A. R. and Morris, A. H.,
           Algorithm 708: Significant Digit Computation of the Incomplete Beta
           Function Ratios. ACM Trans. Math. Softw. 18 (1993), 360-373.


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

    Computes ``ln(abs(beta(a, b)))``.
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
    r"""
    btdtr(a, b, x)

    Cumulative density function of the beta distribution.

    Returns the integral from zero to `x` of the beta probability density
    function,

    .. math::
        I = \int_0^x \frac{\Gamma(a + b)}{\Gamma(a)\Gamma(b)} t^{a-1} (1-t)^{b-1}\,dt

    where :math:`\Gamma` is the gamma function.

    Parameters
    ----------
    a : array_like
        Shape parameter (a > 0).
    b : array_like
        Shape parameter (b > 0).
    x : array_like
        Upper limit of integration, in [0, 1].

    Returns
    -------
    I : ndarray
        Cumulative density function of the beta distribution with parameters
        `a` and `b` at `x`.

    See Also
    --------
    betainc

    Notes
    -----
    This function is identical to the incomplete beta integral function
    `betainc`.

    Wrapper for the Cephes [1]_ routine `btdtr`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

    """)

add_newdoc("scipy.special", "btdtri",
    r"""
    btdtri(a, b, p)

    The `p`-th quantile of the beta distribution.

    This function is the inverse of the beta cumulative distribution function,
    `btdtr`, returning the value of `x` for which `btdtr(a, b, x) = p`, or

    .. math::
        p = \int_0^x \frac{\Gamma(a + b)}{\Gamma(a)\Gamma(b)} t^{a-1} (1-t)^{b-1}\,dt

    Parameters
    ----------
    a : array_like
        Shape parameter (`a` > 0).
    b : array_like
        Shape parameter (`b` > 0).
    p : array_like
        Cumulative probability, in [0, 1].

    Returns
    -------
    x : ndarray
        The quantile corresponding to `p`.

    See Also
    --------
    betaincinv
    btdtr

    Notes
    -----
    The value of `x` is found by interval halving or Newton iterations.

    Wrapper for the Cephes [1]_ routine `incbi`, which solves the equivalent
    problem of finding the inverse of the incomplete beta integral.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

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

    See Also
    --------
    wofz, erf, erfc, erfcx, erfi

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-15, 15, num=1000)
    >>> plt.plot(x, special.dawsn(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$dawsn(x)$')
    >>> plt.show()

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

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellpe`.

    For `m > 0` the computation uses the approximation,

    .. math:: E(m) \\approx P(1-m) - (1-m) \\log(1-m) Q(1-m),

    where :math:`P` and :math:`Q` are tenth-order polynomials.  For
    `m < 0`, the relation

    .. math:: E(m) = E(m/(m - 1)) \\sqrt(1-m)

    is used.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipeinc : Incomplete elliptic integral of the second kind

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
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

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellie`.

    Computation uses arithmetic-geometric means algorithm.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "ellipj",
    """
    ellipj(u, m)

    Jacobian elliptic functions

    Calculates the Jacobian elliptic functions of parameter `m` between
    0 and 1, and real argument `u`.

    Parameters
    ----------
    m : array_like
        Parameter.
    u : array_like
        Argument.

    Returns
    -------
    sn, cn, dn, ph : ndarrays
        The returned functions::

            sn(u|m), cn(u|m), dn(u|m)

        The value `ph` is such that if `u = ellipk(ph, m)`,
        then `sn(u|m) = sin(ph)` and `cn(u|m) = cos(ph)`.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellpj`.

    These functions are periodic, with quarter-period on the real axis
    equal to the complete elliptic integral `ellipk(m)`.

    Relation to incomplete elliptic integral: If `u = ellipk(phi,m)`, then
    `sn(u|m) = sin(phi)`, and `cn(u|m) = cos(phi)`.  The `phi` is called
    the amplitude of `u`.

    Computation is by means of the arithmetic-geometric mean algorithm,
    except when `m` is within 1e-9 of 0 or 1.  In the latter case with `m`
    close to 1, the approximation applies only for `phi < pi/2`.

    See also
    --------
    ellipk : Complete elliptic integral of the first kind.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
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
        Defines the parameter of the elliptic integral as `m = 1 - p`.

    Returns
    -------
    K : ndarray
        Value of the elliptic integral.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellpk`.

    For `p <= 1`, computation uses the approximation,

    .. math:: K(p) \\approx P(p) - \\log(p) Q(p),

    where :math:`P` and :math:`Q` are tenth-order polynomials.  The
    argument `p` is used internally rather than `m` so that the logarithmic
    singularity at `m = 1` will be shifted to the origin; this preserves
    maximum accuracy.  For `p > 1`, the identity

    .. math:: K(p) = K(1/p)/\\sqrt(p)

    is used.

    See Also
    --------
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "ellipkinc",
    """
    ellipkinc(phi, m)

    Incomplete elliptic integral of the first kind

    This function is defined as

    .. math:: K(\\phi, m) = \\int_0^{\\phi} [1 - m \\sin(t)^2]^{-1/2} dt

    This function is also called `F(phi, m)`.

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
    Wrapper for the Cephes [1]_ routine `ellik`.  The computation is
    carried out using the arithmetic-geometric mean algorithm.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
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
    erfc, erfinv, erfcinv, wofz, erfcx, erfi

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

    Examples
    --------
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.erf(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$erf(x)$')
    >>> plt.show()

    """)

add_newdoc("scipy.special", "erfc",
    """
    erfc(x)

    Complementary error function, ``1 - erf(x)``.

    See Also
    --------
    erf, erfi, erfcx, dawsn, wofz

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.erfc(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$erfc(x)$')
    >>> plt.show()

    """)

add_newdoc("scipy.special", "erfi",
    """
    erfi(z)

    Imaginary error function, ``-i erf(i z)``.

    See Also
    --------
    erf, erfc, erfcx, dawsn, wofz

    Notes
    -----

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.erfi(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$erfi(x)$')
    >>> plt.show()

    """)

add_newdoc("scipy.special", "erfcx",
    """
    erfcx(x)

    Scaled complementary error function, ``exp(x**2) * erfc(x)``.

    See Also
    --------
    erf, erfc, erfi, dawsn, wofz

    Notes
    -----

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.erfcx(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$erfcx(x)$')
    >>> plt.show()

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
    see `ufuncs <https://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_

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
    r"""
    fdtr(dfn, dfd, x)

    F cumulative distribution function.

    Returns the value of the cumulative density function of the
    F-distribution, also known as Snedecor's F-distribution or the
    Fisher-Snedecor distribution.

    The F-distribution with parameters :math:`d_n` and :math:`d_d` is the
    distribution of the random variable,

    .. math::
        X = \frac{U_n/d_n}{U_d/d_d},

    where :math:`U_n` and :math:`U_d` are random variables distributed
    :math:`\chi^2`, with :math:`d_n` and :math:`d_d` degrees of freedom,
    respectively.

    Parameters
    ----------
    dfn : array_like
        First parameter (positive float).
    dfd : array_like
        Second parameter (positive float).
    x : array_like
        Argument (nonnegative float).

    Returns
    -------
    y : ndarray
        The CDF of the F-distribution with parameters `dfn` and `dfd` at `x`.

    Notes
    -----
    The regularized incomplete beta function is used, according to the
    formula,

    .. math::
        F(d_n, d_d; x) = I_{xd_n/(d_d + xd_n)}(d_n/2, d_d/2).

    Wrapper for the Cephes [1]_ routine `fdtr`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

    """)

add_newdoc("scipy.special", "fdtrc",
    r"""
    fdtrc(dfn, dfd, x)

    F survival function.

    Returns the complemented F-distribution function (the integral of the
    density from `x` to infinity).

    Parameters
    ----------
    dfn : array_like
        First parameter (positive float).
    dfd : array_like
        Second parameter (positive float).
    x : array_like
        Argument (nonnegative float).

    Returns
    -------
    y : ndarray
        The complemented F-distribution function with parameters `dfn` and
        `dfd` at `x`.

    See also
    --------
    fdtr

    Notes
    -----
    The regularized incomplete beta function is used, according to the
    formula,

    .. math::
        F(d_n, d_d; x) = I_{d_d/(d_d + xd_n)}(d_d/2, d_n/2).

    Wrapper for the Cephes [1]_ routine `fdtrc`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "fdtri",
    r"""
    fdtri(dfn, dfd, p)

    The `p`-th quantile of the F-distribution.

    This function is the inverse of the F-distribution CDF, `fdtr`, returning
    the `x` such that `fdtr(dfn, dfd, x) = p`.

    Parameters
    ----------
    dfn : array_like
        First parameter (positive float).
    dfd : array_like
        Second parameter (positive float).
    p : array_like
        Cumulative probability, in [0, 1].

    Returns
    -------
    x : ndarray
        The quantile corresponding to `p`.

    Notes
    -----
    The computation is carried out using the relation to the inverse
    regularized beta function, :math:`I^{-1}_x(a, b)`.  Let
    :math:`z = I^{-1}_p(d_d/2, d_n/2).`  Then,

    .. math::
        x = \frac{d_d (1 - z)}{d_n z}.

    If `p` is such that :math:`x < 0.5`, the following relation is used
    instead for improved stability: let
    :math:`z' = I^{-1}_{1 - p}(d_n/2, d_d/2).` Then,

    .. math::
        x = \frac{d_d z'}{d_n (1 - z')}.

    Wrapper for the Cephes [1]_ routine `fdtri`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

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

    Performs a logarithmic transformation of the
    values of the gamma function in one of two
    ways, depending on the input `z`:

    1) `z` is not complex (i.e. `z` is a purely
       real number *or* it is array_like and
       contains purely real elements)

    The natural logarithm of the absolute value of
    gamma(z) is computed. Thus, it is defined as:

        ln(abs(gamma(z)))

    2) `z` is complex (i.e. `z` is a complex
       number *or* it is array_like and contains
       at least one complex element)

    The natural logarithm of gamma(z) is computed.
    Thus, it is defined as:

        ln((gamma(z))

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
    r"""
    gdtr(a, b, x)

    Gamma distribution cumulative density function.

    Returns the integral from zero to `x` of the gamma probability density
    function,

    .. math::

        F = \int_0^x \frac{a^b}{\Gamma(b)} t^{b-1} e^{-at}\,dt,

    where :math:`\Gamma` is the gamma function.

    Parameters
    ----------
    a : array_like
        The rate parameter of the gamma distribution, sometimes denoted
        :math:`\beta` (float).  It is also the reciprocal of the scale
        parameter :math:`\theta`.
    b : array_like
        The shape parameter of the gamma distribution, sometimes denoted
        :math:`\alpha` (float).
    x : array_like
        The quantile (upper limit of integration; float).

    See also
    --------
    gdtrc : 1 - CDF of the gamma distribution.

    Returns
    -------
    F : ndarray
        The CDF of the gamma distribution with parameters `a` and `b`
        evaluated at `x`.

    Notes
    -----
    The evaluation is carried out using the relation to the incomplete gamma
    integral (regularized gamma function).

    Wrapper for the Cephes [1]_ routine `gdtr`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

    """)

add_newdoc("scipy.special", "gdtrc",
    r"""
    gdtrc(a, b, x)

    Gamma distribution survival function.

    Integral from `x` to infinity of the gamma probability density function,

    .. math::

        F = \int_x^\infty \frac{a^b}{\Gamma(b)} t^{b-1} e^{-at}\,dt,

    where :math:`\Gamma` is the gamma function.

    Parameters
    ----------
    a : array_like
        The rate parameter of the gamma distribution, sometimes denoted
        :math:`\beta` (float).  It is also the reciprocal of the scale
        parameter :math:`\theta`.
    b : array_like
        The shape parameter of the gamma distribution, sometimes denoted
        :math:`\alpha` (float).
    x : array_like
        The quantile (lower limit of integration; float).

    Returns
    -------
    F : ndarray
        The survival function of the gamma distribution with parameters `a`
        and `b` evaluated at `x`.

    See Also
    --------
    gdtr, gdtri

    Notes
    -----
    The evaluation is carried out using the relation to the incomplete gamma
    integral (regularized gamma function).

    Wrapper for the Cephes [1]_ routine `gdtrc`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

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

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdfgam`.

    The cumulative distribution function `p` is computed using a routine by
    DiDinato and Morris [2]_.  Computation of `a` involves a seach for a value
    that produces the desired value of `p`.  The search relies on the
    monotinicity of `p` with `a`.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] DiDinato, A. R. and Morris, A. H.,
           Computation of the incomplete gamma function ratios and their
           inverse.  ACM Trans. Math. Softw. 12 (1986), 377-393.

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

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdfgam`.

    The cumulative distribution function `p` is computed using a routine by
    DiDinato and Morris [2]_.  Computation of `b` involves a seach for a value
    that produces the desired value of `p`.  The search relies on the
    monotinicity of `p` with `b`.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] DiDinato, A. R. and Morris, A. H.,
           Computation of the incomplete gamma function ratios and their
           inverse.  ACM Trans. Math. Softw. 12 (1986), 377-393.

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

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdfgam`.

    The cumulative distribution function `p` is computed using a routine by
    DiDinato and Morris [2]_.  Computation of `x` involves a seach for a value
    that produces the desired value of `p`.  The search relies on the
    monotinicity of `p` with `x`.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] DiDinato, A. R. and Morris, A. H.,
           Computation of the incomplete gamma function ratios and their
           inverse.  ACM Trans. Math. Softw. 12 (1986), 377-393.

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
    r"""
    hankel1(v, z)

    Hankel function of the first kind

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    out : Values of the Hankel function of the first kind.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(1)}_v(z) = \frac{2}{\imath\pi} \exp(-\imath \pi v/2) K_v(z \exp(-\imath\pi/2))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(1)}_{-v}(z) = H^{(1)}_v(z) \exp(\imath\pi v)

    is used.

    See also
    --------
    hankel1e : this function with leading exponential behavior stripped off.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "hankel1e",
    r"""
    hankel1e(v, z)

    Exponentially scaled Hankel function of the first kind

    Defined as::

        hankel1e(v, z) = hankel1(v, z) * exp(-1j * z)

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    out : Values of the exponentially scaled Hankel function.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(1)}_v(z) = \frac{2}{\imath\pi} \exp(-\imath \pi v/2) K_v(z \exp(-\imath\pi/2))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(1)}_{-v}(z) = H^{(1)}_v(z) \exp(\imath\pi v)

    is used.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "hankel2",
    r"""
    hankel2(v, z)

    Hankel function of the second kind

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    out : Values of the Hankel function of the second kind.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(2)}_v(z) = -\frac{2}{\imath\pi} \exp(\imath \pi v/2) K_v(z \exp(\imath\pi/2))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(2)}_{-v}(z) = H^{(2)}_v(z) \exp(-\imath\pi v)

    is used.

    See also
    --------
    hankel2e : this function with leading exponential behavior stripped off.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "hankel2e",
    r"""
    hankel2e(v, z)

    Exponentially scaled Hankel function of the second kind

    Defined as::

        hankel2e(v, z) = hankel2(v, z) * exp(1j * z)

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    out : Values of the exponentially scaled Hankel function of the second kind.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(2)}_v(z) = -\frac{2}{\imath\pi} \exp(\frac{\imath \pi v}{2}) K_v(z exp(\frac{\imath\pi}{2}))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(2)}_{-v}(z) = H^{(2)}_v(z) \exp(-\imath\pi v)

    is used.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

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

add_newdoc("scipy.special", "hyp0f1", 
    r"""
    hyp0f1(v, x)

    Confluent hypergeometric limit function 0F1.

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

    .. math:: _0F_1(v, z) = \sum_{k=0}^{\infty}\frac{z^k}{(v)_k k!}.

    It's also the limit as :math:`q \to \infty` of :math:`_1F_1(q; v; z/q)`,
    and satisfies the differential equation :math:`f''(z) + vf'(z) = f(z)`.
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
    r"""
    i0(x)

    Modified Bessel function of order 0.

    Defined as,

    .. math::
        I_0(x) = \sum_{k=0}^\infty \frac{(x^2/4)^k}{(k!)^2} = J_0(\imath x),

    where :math:`J_0` is the Bessel function of the first kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float)

    Returns
    -------
    I : ndarray
        Value of the modified Bessel function of order 0 at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `i0`.

    See also
    --------
    iv
    i0e

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "i0e",
    """
    i0e(x)

    Exponentially scaled modified Bessel function of order 0.

    Defined as::

        i0e(x) = exp(-abs(x)) * i0(x).

    Parameters
    ----------
    x : array_like
        Argument (float)

    Returns
    -------
    I : ndarray
        Value of the exponentially scaled modified Bessel function of order 0
        at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval.  The
    polynomial expansions used are the same as those in `i0`, but
    they are not multiplied by the dominant exponential factor.

    This function is a wrapper for the Cephes [1]_ routine `i0e`.

    See also
    --------
    iv
    i0

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "i1",
    r"""
    i1(x)

    Modified Bessel function of order 1.

    Defined as,

    .. math::
        I_1(x) = \frac{1}{2}x \sum_{k=0}^\infty \frac{(x^2/4)^k}{k! (k + 1)!}
               = -\imath J_1(\imath x),

    where :math:`J_1` is the Bessel function of the first kind of order 1.

    Parameters
    ----------
    x : array_like
        Argument (float)

    Returns
    -------
    I : ndarray
        Value of the modified Bessel function of order 1 at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `i1`.

    See also
    --------
    iv
    i1e

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "i1e",
    """
    i1e(x)

    Exponentially scaled modified Bessel function of order 1.

    Defined as::

        i1e(x) = exp(-abs(x)) * i1(x)

    Parameters
    ----------
    x : array_like
        Argument (float)

    Returns
    -------
    I : ndarray
        Value of the exponentially scaled modified Bessel function of order 1
        at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval. The
    polynomial expansions used are the same as those in `i1`, but
    they are not multiplied by the dominant exponential factor.

    This function is a wrapper for the Cephes [1]_ routine `i1e`.

    See also
    --------
    iv
    i1

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
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
    r"""
    it2struve0(x)

    Integral related to the Struve function of order 0.

    Returns the integral,

    .. math::
        \int_x^\infty \frac{H_0(t)}{t}\,dt

    where :math:`H_0` is the Struve function of order 0.

    Parameters
    ----------
    x : array_like
        Lower limit of integration.

    Returns
    -------
    I : ndarray
        The value of the integral.

    See also
    --------
    struve

    Notes
    -----
    Wrapper for a Fortran routine created by Shanjie Zhang and Jianming
    Jin [1]_.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html
    """)

add_newdoc("scipy.special", "itairy",
    """
    itairy(x)

    Integrals of Airy functions

    Calculates the integrals of Airy functions from 0 to `x`.

    Parameters
    ----------

    x: array_like
        Upper limit of integration (float).

    Returns
    -------
    Apt
        Integral of Ai(t) from 0 to x.
    Bpt
        Integral of Bi(t) from 0 to x.
    Ant
        Integral of Ai(-t) from 0 to x.
    Bnt
        Integral of Bi(-t) from 0 to x.

    Notes
    -----

    Wrapper for a Fortran routine created by Shanjie Zhang and Jianming
    Jin [1]_.

    References
    ----------

    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html
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
    r"""
    itmodstruve0(x)

    Integral of the modified Struve function of order 0.

    .. math::
        I = \int_0^x L_0(t)\,dt

    Parameters
    ----------
    x : array_like
        Upper limit of integration (float).

    Returns
    -------
    I : ndarray
        The integral of :math:`L_0` from 0 to `x`.

    Notes
    -----
    Wrapper for a Fortran routine created by Shanjie Zhang and Jianming
    Jin [1]_.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """)

add_newdoc("scipy.special", "itstruve0",
    r"""
    itstruve0(x)

    Integral of the Struve function of order 0.

    .. math::
        I = \int_0^x H_0(t)\,dt

    Parameters
    ----------
    x : array_like
        Upper limit of integration (float).

    Returns
    -------
    I : ndarray
        The integral of :math:`H_0` from 0 to `x`.

    See also
    --------
    struve

    Notes
    -----
    Wrapper for a Fortran routine created by Shanjie Zhang and Jianming
    Jin [1]_.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           http://jin.ece.illinois.edu/specfunc.html

    """)

add_newdoc("scipy.special", "iv",
    r"""
    iv(v, z)

    Modified Bessel function of the first kind of real order.

    Parameters
    ----------
    v : array_like
        Order. If `z` is of real type and negative, `v` must be integer
        valued.
    z : array_like of float or complex
        Argument.

    Returns
    -------
    out : ndarray
        Values of the modified Bessel function.

    Notes
    -----
    For real `z` and :math:`v \in [-50, 50]`, the evaluation is carried out
    using Temme's method [1]_.  For larger orders, uniform asymptotic
    expansions are applied.

    For complex `z` and positive `v`, the AMOS [2]_ `zbesi` routine is
    called. It uses a power series for small `z`, the asymptitic expansion
    for large `abs(z)`, the Miller algorithm normalized by the Wronskian
    and a Neumann series for intermediate magnitudes, and the uniform
    asymptitic expansions for :math:`I_v(z)` and :math:`J_v(z)` for large
    orders.  Backward recurrence is used to generate sequences or reduce
    orders when necessary.

    The calculations above are done in the right half plane and continued
    into the left half plane by the formula,

    .. math:: I_v(z \exp(\pm\imath\pi)) = \exp(\pm\pi v) I_v(z)

    (valid when the real part of `z` is positive).  For negative `v`, the
    formula

    .. math:: I_{-v}(z) = I_v(z) + \frac{2}{\pi} \sin(\pi v) K_v(z)

    is used, where :math:`K_v(z)` is the modified Bessel function of the
    second kind, evaluated using the AMOS routine `zbesk`.

    See also
    --------
    kve : This function with leading exponential behavior stripped off.

    References
    ----------
    .. [1] Temme, Journal of Computational Physics, vol 21, 343 (1976)
    .. [2] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "ive",
    r"""
    ive(v, z)

    Exponentially scaled modified Bessel function of the first kind

    Defined as::

        ive(v, z) = iv(v, z) * exp(-abs(z.real))

    Parameters
    ----------
    v : array_like of float
        Order.
    z : array_like of float or complex
        Argument.

    Returns
    -------
    out : ndarray
        Values of the exponentially scaled modified Bessel function.

    Notes
    -----
    For positive `v`, the AMOS [1]_ `zbesi` routine is called. It uses a
    power series for small `z`, the asymptitic expansion for large
    `abs(z)`, the Miller algorithm normalized by the Wronskian and a
    Neumann series for intermediate magnitudes, and the uniform asymptitic
    expansions for :math:`I_v(z)` and :math:`J_v(z)` for large orders.
    Backward recurrence is used to generate sequences or reduce orders when
    necessary.

    The calculations above are done in the right half plane and continued
    into the left half plane by the formula,

    .. math:: I_v(z \exp(\pm\imath\pi)) = \exp(\pm\pi v) I_v(z)

    (valid when the real part of `z` is positive).  For negative `v`, the
    formula

    .. math:: I_{-v}(z) = I_v(z) + \frac{2}{\pi} \sin(\pi v) K_v(z)

    is used, where :math:`K_v(z)` is the modified Bessel function of the
    second kind, evaluated using the AMOS routine `zbesk`.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "j0",
    r"""
    j0(x)

    Bessel function of the first kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float).

    Returns
    -------
    J : ndarray
        Value of the Bessel function of the first kind of order 0 at `x`.

    Notes
    -----
    The domain is divided into the intervals [0, 5] and (5, infinity). In the
    first interval the following rational approximation is used:

    .. math::

        J_0(x) \approx (w - r_1^2)(w - r_2^2) \frac{P_3(w)}{Q_8(w)},

    where :math:`w = x^2` and :math:`r_1`, :math:`r_2` are the zeros of
    :math:`J_0`, and :math:`P_3` and :math:`Q_8` are polynomials of degrees 3
    and 8, respectively.

    In the second interval, the Hankel asymptotic expansion is employed with
    two rational functions of degree 6/6 and 7/7.

    This function is a wrapper for the Cephes [1]_ routine `j0`.

    See also
    --------
    jv : Bessel function of real order and complex argument.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "j1",
    """
    j1(x)

    Bessel function of the first kind of order 1.

    Parameters
    ----------
    x : array_like
        Argument (float).

    Returns
    -------
    J : ndarray
        Value of the Bessel function of the first kind of order 1 at `x`.

    Notes
    -----
    The domain is divided into the intervals [0, 8] and (8, infinity). In the
    first interval a 24 term Chebyshev expansion is used. In the second, the
    asymptotic trigonometric representation is employed using two rational
    functions of degree 5/5.

    This function is a wrapper for the Cephes [1]_ routine `j1`.

    See also
    --------
    jv

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

    """)

add_newdoc("scipy.special", "jn",
    """
    jn(n, x)

    Bessel function of the first kind of integer order and real argument.

    Notes
    -----
    `jn` is an alias of `jv`.

    See also
    --------
    jv

    """)

add_newdoc("scipy.special", "jv",
    r"""
    jv(v, z)

    Bessel function of the first kind of real order and complex argument.

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    J : ndarray
        Value of the Bessel function, :math:`J_v(z)`.

    Notes
    -----
    For positive `v` values, the computation is carried out using the AMOS
    [1]_ `zbesj` routine, which exploits the connection to the modified
    Bessel function :math:`I_v`,

    .. math::
        J_v(z) = \exp(n\pi\imath/2) I_v(-\imath z)\qquad (\Im z > 0)

        J_v(z) = \exp(-n\pi\imath/2) I_v(\imath z)\qquad (\Im z < 0)

    For negative `v` values the formula,

    .. math:: J_{-v}(z) = J_v(z) \cos(\pi v) - Y_v(z) \sin(\pi v)

    is used, where :math:`Y_v(z)` is the Bessel function of the second
    kind, computed using the AMOS routine `zbesy`.  Note that the second
    term is exactly zero for integer `v`; to improve accuracy the second
    term is explicitly omitted for `v` values such that `v = floor(v)`.

    See also
    --------
    jve : :math:`J_v` with leading exponential behavior stripped off.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "jve",
    r"""
    jve(v, z)

    Exponentially scaled Bessel function of order `v`.

    Defined as::

        jve(v, z) = jv(v, z) * exp(-abs(z.imag))

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    J : ndarray
        Value of the exponentially scaled Bessel function.

    Notes
    -----
    For positive `v` values, the computation is carried out using the AMOS
    [1]_ `zbesj` routine, which exploits the connection to the modified
    Bessel function :math:`I_v`,

    .. math::
        J_v(z) = \exp(n\pi\imath/2) I_v(-\imath z)\qquad (\Im z > 0)

        J_v(z) = \exp(-n\pi\imath/2) I_v(\imath z)\qquad (\Im z < 0)

    For negative `v` values the formula,

    .. math:: J_{-v}(z) = J_v(z) \cos(\pi v) - Y_v(z) \sin(\pi v)

    is used, where :math:`Y_v(z)` is the Bessel function of the second
    kind, computed using the AMOS routine `zbesy`.  Note that the second
    term is exactly zero for integer `v`; to improve accuracy the second
    term is explicitly omitted for `v` values such that `v = floor(v)`.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    """)

add_newdoc("scipy.special", "k0",
    r"""
    k0(x)

    Modified Bessel function of the second kind of order 0, :math:`K_0`.

    This function is also sometimes referred to as the modified Bessel
    function of the third kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float).

    Returns
    -------
    K : ndarray
        Value of the modified Bessel function :math:`K_0` at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k0`.

    See also
    --------
    kv
    k0e

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "k0e",
    """
    k0e(x)

    Exponentially scaled modified Bessel function K of order 0

    Defined as::

        k0e(x) = exp(x) * k0(x).

    Parameters
    ----------
    x : array_like
        Argument (float)

    Returns
    -------
    K : ndarray
        Value of the exponentially scaled modified Bessel function K of order
        0 at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k0e`.

    See also
    --------
    kv
    k0

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "k1",
    """
    k1(x)

    Modified Bessel function of the second kind of order 1, :math:`K_1(x)`.

    Parameters
    ----------
    x : array_like
        Argument (float)

    Returns
    -------
    K : ndarray
        Value of the modified Bessel function K of order 1 at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k1`.

    See also
    --------
    kv
    k1e

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "k1e",
    """
    k1e(x)

    Exponentially scaled modified Bessel function K of order 1

    Defined as::

        k1e(x) = exp(x) * k1(x)

    Parameters
    ----------
    x : array_like
        Argument (float)

    Returns
    -------
    K : ndarray
        Value of the exponentially scaled modified Bessel function K of order
        1 at `x`.

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k1e`.

    See also
    --------
    kv
    k1

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
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
    r"""
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

    Notes
    -----
    Wrapper for AMOS [1]_ routine `zbesk`.  For a discussion of the
    algorithm used, see [2]_ and the references therein.

    See Also
    --------
    kv : Same function, but accepts real order and complex argument
    kvp : Derivative of this function

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    .. [2] Donald E. Amos, "Algorithm 644: A portable package for Bessel
           functions of a complex argument and nonnegative order", ACM
           TOMS Vol. 12 Issue 3, Sept. 1986, p. 265

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

    >>> kn([4, 5, 6], 1)
    array([   44.23241585,   360.9605896 ,  3653.83831186])
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
    r"""
    kv(v, z)

    Modified Bessel function of the second kind of real order `v`

    Returns the modified Bessel function of the second kind for real order
    `v` at complex `z`.

    These are also sometimes called functions of the third kind, Basset
    functions, or Macdonald functions.  They are defined as those solutions
    of the modified Bessel equation for which,

    .. math::
        K_v(x) \sim \sqrt{\pi/(2x)} \exp(-x)

    as :math:`x \to \infty` [3]_.

    Parameters
    ----------
    v : array_like of float
        Order of Bessel functions
    z : array_like of complex
        Argument at which to evaluate the Bessel functions

    Returns
    -------
    out : ndarray
        The results. Note that input must be of complex type to get complex
        output, e.g. ``kv(3, -2+0j)`` instead of ``kv(3, -2)``.

    Notes
    -----
    Wrapper for AMOS [1]_ routine `zbesk`.  For a discussion of the
    algorithm used, see [2]_ and the references therein.

    See Also
    --------
    kve : This function with leading exponential behavior stripped off.
    kvp : Derivative of this function

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    .. [2] Donald E. Amos, "Algorithm 644: A portable package for Bessel
           functions of a complex argument and nonnegative order", ACM
           TOMS Vol. 12 Issue 3, Sept. 1986, p. 265
    .. [3] NIST Digital Library of Mathematical Functions,
           Eq. 10.25.E3. http://dlmf.nist.gov/10.25.E3

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

    >>> kv([4, 4.5, 5], 1+2j)
    array([ 0.1992+2.3892j,  2.3493+3.6j   ,  7.2827+3.8104j])

    """)

add_newdoc("scipy.special", "kve",
    r"""
    kve(v, z)

    Exponentially scaled modified Bessel function of the second kind.

    Returns the exponentially scaled, modified Bessel function of the
    second kind (sometimes called the third kind) for real order `v` at
    complex `z`::

        kve(v, z) = kv(v, z) * exp(z)

    Parameters
    ----------
    v : array_like of float
        Order of Bessel functions
    z : array_like of complex
        Argument at which to evaluate the Bessel functions

    Returns
    -------
    out : ndarray
        The exponentially scaled modified Bessel function of the second kind.

    Notes
    -----
    Wrapper for AMOS [1]_ routine `zbesk`.  For a discussion of the
    algorithm used, see [2]_ and the references therein.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    .. [2] Donald E. Amos, "Algorithm 644: A portable package for Bessel
           functions of a complex argument and nonnegative order", ACM
           TOMS Vol. 12 Issue 3, Sept. 1986, p. 265
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
    see `ufuncs <https://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_

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
    v : float
        Degree.
    x : float
        Argument. Must be ``|x| <= 1``.

    Returns
    -------
    res : float
        The value of the function.

    See Also
    --------
    lpmn : Similar, but computes values for all orders 0..m and degrees 0..n.
    clpmn : Similar to `lpmn` but allows a complex argument.

    Notes
    -----
    It is possible to extend the domain of this function to all
    complex m, v, x, but this is not yet implemented.

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
    r"""
    modstruve(v, x)

    Modified Struve function.

    Return the value of the modified Struve function of order `v` at `x`.  The
    modified Struve function is defined as,

    .. math::
        L_v(x) = -\imath \exp(-\pi\imath v/2) H_v(x),

    where :math:`H_v` is the Struve function.

    Parameters
    ----------
    v : array_like
        Order of the modified Struve function (float).
    x : array_like
        Argument of the Struve function (float; must be positive unless `v` is
        an integer).

    Returns
    -------
    L : ndarray
        Value of the modified Struve function of order `v` at `x`.

    Notes
    -----
    Three methods discussed in [1]_ are used to evaluate the function:

    - power series
    - expansion in Bessel functions (if :math:`|z| < |v| + 20`)
    - asymptotic large-z expansion (if :math:`z \geq 0.7v + 12`)

    Rounding errors are estimated based on the largest terms in the sums, and
    the result associated with the smallest error is returned.

    See also
    --------
    struve

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           http://dlmf.nist.gov/11
    """)

add_newdoc("scipy.special", "nbdtr",
    r"""
    nbdtr(k, n, p)

    Negative binomial cumulative distribution function.

    Returns the sum of the terms 0 through `k` of the negative binomial
    distribution probability mass function,

    .. math::

        F = \sum_{j=0}^k {{n + j - 1}\choose{j}} p^n (1 - p)^j.

    In a sequence of Bernoulli trials with individual success probabilities
    `p`, this is the probability that `k` or fewer failures precede the nth
    success.

    Parameters
    ----------
    k : array_like
        The maximum number of allowed failures (nonnegative int).
    n : array_like
        The target number of successes (positive int).
    p : array_like
        Probability of success in a single event (float).

    Returns
    -------
    F : ndarray
        The probability of `k` or fewer failures before `n` successes in a
        sequence of events with individual success probability `p`.

    See also
    --------
    nbdtrc

    Notes
    -----
    If floating point values are passed for `k` or `n`, they will be truncated
    to integers.

    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{nbdtr}(k, n, p) = I_{p}(n, k + 1).

    Wrapper for the Cephes [1]_ routine `nbdtr`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html

    """)

add_newdoc("scipy.special", "nbdtrc",
    r"""
    nbdtrc(k, n, p)

    Negative binomial survival function.

    Returns the sum of the terms `k + 1` to infinity of the negative binomial
    distribution probability mass function,

    .. math::

        F = \sum_{j=k + 1}^\infty {{n + j - 1}\choose{j}} p^n (1 - p)^j.

    In a sequence of Bernoulli trials with individual success probabilities
    `p`, this is the probability that more than `k` failures precede the nth
    success.

    Parameters
    ----------
    k : array_like
        The maximum number of allowed failures (nonnegative int).
    n : array_like
        The target number of successes (positive int).
    p : array_like
        Probability of success in a single event (float).

    Returns
    -------
    F : ndarray
        The probability of `k + 1` or more failures before `n` successes in a
        sequence of events with individual success probability `p`.

    Notes
    -----
    If floating point values are passed for `k` or `n`, they will be truncated
    to integers.

    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{nbdtrc}(k, n, p) = I_{1 - p}(k + 1, n).

    Wrapper for the Cephes [1]_ routine `nbdtrc`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "nbdtri",
    """
    nbdtri(k, n, y)

    Inverse of `nbdtr` vs `p`.

    Returns the inverse with respect to the parameter `p` of
    `y = nbdtr(k, n, p)`, the negative binomial cumulative distribution
    function.

    Parameters
    ----------
    k : array_like
        The maximum number of allowed failures (nonnegative int).
    n : array_like
        The target number of successes (positive int).
    y : array_like
        The probability of `k` or fewer failures before `n` successes (float).

    Returns
    -------
    p : ndarray
        Probability of success in a single event (float) such that
        `nbdtr(k, n, p) = y`.

    See also
    --------
    nbdtr : Cumulative distribution function of the negative binomial.
    nbdtrik : Inverse with respect to `k` of `nbdtr(k, n, p)`.
    nbdtrin : Inverse with respect to `n` of `nbdtr(k, n, p)`.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `nbdtri`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    
    """)

add_newdoc("scipy.special", "nbdtrik",
    r"""
    nbdtrik(y, n, p)

    Inverse of `nbdtr` vs `k`.

    Returns the inverse with respect to the parameter `k` of
    `y = nbdtr(k, n, p)`, the negative binomial cumulative distribution
    function.

    Parameters
    ----------
    y : array_like
        The probability of `k` or fewer failures before `n` successes (float).
    n : array_like
        The target number of successes (positive int).
    p : array_like
        Probability of success in a single event (float).
    
    Returns
    -------
    k : ndarray
        The maximum number of allowed failures such that `nbdtr(k, n, p) = y`.
    
    See also
    --------
    nbdtr : Cumulative distribution function of the negative binomial.
    nbdtri : Inverse with respect to `p` of `nbdtr(k, n, p)`.
    nbdtrin : Inverse with respect to `n` of `nbdtr(k, n, p)`.

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdfnbn`.

    Formula 26.5.26 of [2]_,

    .. math::
        \sum_{j=k + 1}^\infty {{n + j - 1}\choose{j}} p^n (1 - p)^j = I_{1 - p}(k + 1, n),

    is used to reduce calculation of the cumulative distribution function to
    that of a regularized incomplete beta :math:`I`.

    Computation of `k` involves a seach for a value that produces the desired
    value of `y`.  The search relies on the monotinicity of `y` with `k`.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.

    """)

add_newdoc("scipy.special", "nbdtrin",
    r"""
    nbdtrin(k, y, p)

    Inverse of `nbdtr` vs `n`.

    Returns the inverse with respect to the parameter `n` of
    `y = nbdtr(k, n, p)`, the negative binomial cumulative distribution
    function.

    Parameters
    ----------
    k : array_like
        The maximum number of allowed failures (nonnegative int).
    y : array_like
        The probability of `k` or fewer failures before `n` successes (float).
    p : array_like
        Probability of success in a single event (float).
    
    Returns
    -------
    n : ndarray
        The number of successes `n` such that `nbdtr(k, n, p) = y`.
    
    See also
    --------
    nbdtr : Cumulative distribution function of the negative binomial.
    nbdtri : Inverse with respect to `p` of `nbdtr(k, n, p)`.
    nbdtrik : Inverse with respect to `k` of `nbdtr(k, n, p)`.

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdfnbn`.

    Formula 26.5.26 of [2]_,

    .. math::
        \sum_{j=k + 1}^\infty {{n + j - 1}\choose{j}} p^n (1 - p)^j = I_{1 - p}(k + 1, n),

    is used to reduce calculation of the cumulative distribution function to
    that of a regularized incomplete beta :math:`I`.

    Computation of `n` involves a seach for a value that produces the desired
    value of `y`.  The search relies on the monotinicity of `y` with `n`.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.

    """)

add_newdoc("scipy.special", "ncfdtr",
    r"""
    ncfdtr(dfn, dfd, nc, f)

    Cumulative distribution function of the non-central F distribution.

    The non-central F describes the distribution of,

    .. math::
        Z = \frac{X/d_n}{Y/d_d}

    where :math:`X` and :math:`Y` are independently distributed, with
    :math:`X` distributed non-central :math:`\chi^2` with noncentrality
    parameter `nc` and :math:`d_n` degrees of freedom, and :math:`Y`
    distributed :math:`\chi^2` with :math:`d_d` degrees of freedom.

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

    Notes
    -----
    Wrapper for the CDFLIB [1]_ Fortran routine `cdffnc`.

    The cumulative distribution function is computed using Formula 26.6.20 of
    [2]_:

    .. math::
        F(d_n, d_d, n_c, f) = \sum_{j=0}^\infty e^{-n_c/2} \frac{(n_c/2)^j}{j!} I_{x}(\frac{d_n}{2} + j, \frac{d_d}{2}),

    where :math:`I` is the regularized incomplete beta function, and
    :math:`x = f d_n/(f d_n + d_d)`.

    The computation time required for this routine is proportional to the
    noncentrality parameter `nc`.  Very large values of this parameter can
    consume immense computer resources.  This is why the search range is
    bounded by 10,000.

    References
    ----------
    .. [1] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.

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

    Notes
    -----
    The value of the cumulative noncentral F distribution is not necessarily
    monotone in either degrees of freedom.  There thus may be two values that
    provide a given CDF value.  This routine assumes monotonicity and will
    find an arbitrary one of the two values.

    """)

add_newdoc("scipy.special", "ncfdtridfn",
    """
    ncfdtridfn(p, f, dfd, nc)

    Calculate degrees of freedom (numerator) for the noncentral F-distribution.

    See `ncfdtr` for more details.

    Notes
    -----
    The value of the cumulative noncentral F distribution is not necessarily
    monotone in either degrees of freedom.  There thus may be two values that
    provide a given CDF value.  This routine assumes monotonicity and will
    find an arbitrary one of the two values.

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
    r"""
    struve(v, x)

    Struve function.

    Return the value of the Struve function of order `v` at `x`.  The Struve
    function is defined as,

    .. math::
        H_v(x) = (z/2)^{v + 1} \sum_{n=0}^\infty \frac{(-1)^n (z/2)^{2n}}{\Gamma(n + \frac{3}{2}) \Gamma(n + v + \frac{3}{2})},

    where :math:`\Gamma` is the gamma function.

    Parameters
    ----------
    v : array_like
        Order of the Struve function (float).
    x : array_like
        Argument of the Struve function (float; must be positive unless `v` is
        an integer).

    Returns
    -------
    H : ndarray
        Value of the Struve function of order `v` at `x`.

    Notes
    -----
    Three methods discussed in [1]_ are used to evaluate the Struve function:

    - power series
    - expansion in Bessel functions (if :math:`|z| < |v| + 20`)
    - asymptotic large-z expansion (if :math:`z \geq 0.7v + 12`)

    Rounding errors are estimated based on the largest terms in the sums, and
    the result associated with the smallest error is returned.

    See also
    --------
    modstruve

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           http://dlmf.nist.gov/11

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

        exp(-z**2) * erfc(-i*z)

    See Also
    --------
    dawsn, erf, erfc, erfcx, erfi

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.wofz(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$wofz(x)$')
    >>> plt.show()

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
    r"""
    y0(x)

    Bessel function of the second kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float).

    Returns
    -------
    Y : ndarray
        Value of the Bessel function of the second kind of order 0 at `x`.

    Notes
    -----
    
    The domain is divided into the intervals [0, 5] and (5, infinity). In the
    first interval a rational approximation :math:`R(x)` is employed to
    compute,
    
    .. math::
        
        Y_0(x) = R(x) + \frac{2 \log(x) J_0(x)}{\pi},

    where :math:`J_0` is the Bessel function of the first kind of order 0.

    In the second interval, the Hankel asymptotic expansion is employed with
    two rational functions of degree 6/6 and 7/7.

    This function is a wrapper for the Cephes [1]_ routine `y0`.

    See also
    --------
    j0
    yv

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "y1",
    """
    y1(x)

    Bessel function of the second kind of order 1.

    Parameters
    ----------
    x : array_like
        Argument (float).

    Returns
    -------
    Y : ndarray
        Value of the Bessel function of the second kind of order 1 at `x`.

    Notes
    -----

    The domain is divided into the intervals [0, 8] and (8, infinity). In the
    first interval a 25 term Chebyshev expansion is used, and computing
    :math:`J_1` (the Bessel function of the first kind) is required. In the
    second, the asymptotic trigonometric representation is employed using two
    rational functions of degree 5/5.

    This function is a wrapper for the Cephes [1]_ routine `y1`.

    See also
    --------
    j1
    yn
    yv
    
    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "yn",
    r"""
    yn(n, x)

    Bessel function of the second kind of integer order and real argument.

    Parameters
    ----------
    n : array_like
        Order (integer).
    z : array_like
        Argument (float).

    Returns
    -------
    Y : ndarray
        Value of the Bessel function, :math:`Y_n(x)`.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `yn`.

    The function is evaluated by forward recurrence on `n`, starting with
    values computed by the Cephes routines `y0` and `y1`. If `n = 0` or 1,
    the routine for `y0` or `y1` is called directly.

    See also
    --------
    yv : For real order and real or complex argument.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/index.html
    """)

add_newdoc("scipy.special", "yv",
    r"""
    yv(v, z)

    Bessel function of the second kind of real order and complex argument.

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    Y : ndarray
        Value of the Bessel function of the second kind, :math:`Y_v(x)`.

    Notes
    -----
    For positive `v` values, the computation is carried out using the
    AMOS [1]_ `zbesy` routine, which exploits the connection to the Hankel
    Bessel functions :math:`H_v^{(1)}` and :math:`H_v^{(2)}`,

    .. math:: Y_v(z) = \frac{1}{2\imath} (H_v^{(1)} - H_v^{(2)}).

    For negative `v` values the formula,

    .. math:: Y_{-v}(z) = Y_v(z) \cos(\pi v) + J_v(z) \sin(\pi v)

    is used, where :math:`J_v(z)` is the Bessel function of the first kind,
    computed using the AMOS routine `zbesj`.  Note that the second term is
    exactly zero for integer `v`; to improve accuracy the second term is
    explicitly omitted for `v` values such that `v = floor(v)`.

    See also
    --------
    yve : :math:`Y_v` with leading exponential behavior stripped off.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    """)

add_newdoc("scipy.special", "yve",
    r"""
    yve(v, z)

    Exponentially scaled Bessel function of the second kind of real order.

    Returns the exponentially scaled Bessel function of the second
    kind of real order `v` at complex `z`::

        yve(v, z) = yv(v, z) * exp(-abs(z.imag))

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).

    Returns
    -------
    Y : ndarray
        Value of the exponentially scaled Bessel function.

    Notes
    -----
    For positive `v` values, the computation is carried out using the
    AMOS [1]_ `zbesy` routine, which exploits the connection to the Hankel
    Bessel functions :math:`H_v^{(1)}` and :math:`H_v^{(2)}`,

    .. math:: Y_v(z) = \frac{1}{2\imath} (H_v^{(1)} - H_v^{(2)}).

    For negative `v` values the formula,

    .. math:: Y_{-v}(z) = Y_v(z) \cos(\pi v) + J_v(z) \sin(\pi v)

    is used, where :math:`J_v(z)` is the Bessel function of the first kind,
    computed using the AMOS routine `zbesj`.  Note that the second term is
    exactly zero for integer `v`; to improve accuracy the second term is
    explicitly omitted for `v` values such that `v = floor(v)`.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
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

add_newdoc("scipy.special", "_spherical_jn",
    """
    Internal function, use `spherical_jn` instead.
    """)

add_newdoc("scipy.special", "_spherical_jn_d",
    """
    Internal function, use `spherical_jn` instead.
    """)

add_newdoc("scipy.special", "_spherical_yn",
    """
    Internal function, use `spherical_yn` instead.
    """)

add_newdoc("scipy.special", "_spherical_yn_d",
    """
    Internal function, use `spherical_yn` instead.
    """)

add_newdoc("scipy.special", "_spherical_in",
    """
    Internal function, use `spherical_in` instead.
    """)

add_newdoc("scipy.special", "_spherical_in_d",
    """
    Internal function, use `spherical_in` instead.
    """)

add_newdoc("scipy.special", "_spherical_kn",
    """
    Internal function, use `spherical_kn` instead.
    """)

add_newdoc("scipy.special", "_spherical_kn_d",
    """
    Internal function, use `spherical_kn` instead.
    """)
