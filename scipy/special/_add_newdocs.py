# Docstrings for generated ufuncs
#
# The syntax is designed to look like the function add_newdoc is being
# called from numpy.lib, but in this file add_newdoc puts the
# docstrings in a dictionary. This dictionary is used in
# _generate_pyx.py to generate the docstrings for the ufuncs in
# scipy.special at the C level when the ufuncs are created at compile
# time.

docdict: dict[str, str] = {}


def get(name):
    return docdict.get(name)


def add_newdoc(name, doc):
    docdict[name] = doc


add_newdoc("_sf_error_test_function",
    """
    Private function; do not use.
    """)


add_newdoc("_cosine_cdf",
    """
    _cosine_cdf(x)

    Cumulative distribution function (CDF) of the cosine distribution::

                 {             0,              x < -pi
        cdf(x) = { (pi + x + sin(x))/(2*pi),   -pi <= x <= pi
                 {             1,              x > pi

    Parameters
    ----------
    x : array_like
        `x` must contain real numbers.

    Returns
    -------
    scalar or ndarray
        The cosine distribution CDF evaluated at `x`.

    """)

add_newdoc("_cosine_invcdf",
    """
    _cosine_invcdf(p)

    Inverse of the cumulative distribution function (CDF) of the cosine
    distribution.

    The CDF of the cosine distribution is::

        cdf(x) = (pi + x + sin(x))/(2*pi)

    This function computes the inverse of cdf(x).

    Parameters
    ----------
    p : array_like
        `p` must contain real numbers in the interval ``0 <= p <= 1``.
        `nan` is returned for values of `p` outside the interval [0, 1].

    Returns
    -------
    scalar or ndarray
        The inverse of the cosine distribution CDF evaluated at `p`.

    """)

add_newdoc("_ellip_harm",
    """
    Internal function, use `ellip_harm` instead.
    """)

add_newdoc("wrightomega",
    r"""
    wrightomega(z, out=None)

    Wright Omega function.

    Defined as the solution to

    .. math::

        \omega + \log(\omega) = z

    where :math:`\log` is the principal branch of the complex logarithm.

    Parameters
    ----------
    z : array_like
        Points at which to evaluate the Wright Omega function
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    omega : scalar or ndarray
        Values of the Wright Omega function

    See Also
    --------
    lambertw : The Lambert W function

    Notes
    -----
    .. versionadded:: 0.19.0

    The function can also be defined as

    .. math::

        \omega(z) = W_{K(z)}(e^z)

    where :math:`K(z) = \lceil (\Im(z) - \pi)/(2\pi) \rceil` is the
    unwinding number and :math:`W` is the Lambert W function.

    The implementation here is taken from [1]_.

    References
    ----------
    .. [1] Lawrence, Corless, and Jeffrey, "Algorithm 917: Complex
           Double-Precision Evaluation of the Wright :math:`\omega`
           Function." ACM Transactions on Mathematical Software,
           2012. :doi:`10.1145/2168773.2168779`.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import wrightomega, lambertw

    >>> wrightomega([-2, -1, 0, 1, 2])
    array([0.12002824, 0.27846454, 0.56714329, 1.        , 1.5571456 ])

    Complex input:

    >>> wrightomega(3 + 5j)
    (1.5804428632097158+3.8213626783287937j)

    Verify that ``wrightomega(z)`` satisfies ``w + log(w) = z``:

    >>> w = -5 + 4j
    >>> wrightomega(w + np.log(w))
    (-5+4j)

    Verify the connection to ``lambertw``:

    >>> z = 0.5 + 3j
    >>> wrightomega(z)
    (0.0966015889280649+1.4937828458191993j)
    >>> lambertw(np.exp(z))
    (0.09660158892806493+1.4937828458191993j)

    >>> z = 0.5 + 4j
    >>> wrightomega(z)
    (-0.3362123489037213+2.282986001579032j)
    >>> lambertw(np.exp(z), k=1)
    (-0.33621234890372115+2.282986001579032j)
    """)

add_newdoc("bdtr",
    r"""
    bdtr(k, n, p, out=None)

    Binomial distribution cumulative distribution function.

    Sum of the terms 0 through `floor(k)` of the Binomial probability density.

    .. math::
        \mathrm{bdtr}(k, n, p) =
        \sum_{j=0}^{\lfloor k \rfloor} {{n}\choose{j}} p^j (1-p)^{n-j}

    Parameters
    ----------
    k : array_like
        Number of successes (double), rounded down to the nearest integer.
    n : array_like
        Number of events (int).
    p : array_like
        Probability of success in a single event (float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    y : scalar or ndarray
        Probability of `floor(k)` or fewer successes in `n` independent events with
        success probabilities of `p`.

    Notes
    -----
    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{bdtr}(k, n, p) =
        I_{1 - p}(n - \lfloor k \rfloor, \lfloor k \rfloor + 1).

    Wrapper for the Cephes [1]_ routine `bdtr`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    """)

add_newdoc("bdtrc",
    r"""
    bdtrc(k, n, p, out=None)

    Binomial distribution survival function.

    Sum of the terms `floor(k) + 1` through `n` of the binomial probability
    density,

    .. math::
        \mathrm{bdtrc}(k, n, p) =
        \sum_{j=\lfloor k \rfloor +1}^n {{n}\choose{j}} p^j (1-p)^{n-j}

    Parameters
    ----------
    k : array_like
        Number of successes (double), rounded down to nearest integer.
    n : array_like
        Number of events (int)
    p : array_like
        Probability of success in a single event.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    y : scalar or ndarray
        Probability of `floor(k) + 1` or more successes in `n` independent
        events with success probabilities of `p`.

    See Also
    --------
    bdtr
    betainc

    Notes
    -----
    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{bdtrc}(k, n, p) = I_{p}(\lfloor k \rfloor + 1, n - \lfloor k \rfloor).

    Wrapper for the Cephes [1]_ routine `bdtrc`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    """)

add_newdoc("bdtri",
    r"""
    bdtri(k, n, y, out=None)

    Inverse function to `bdtr` with respect to `p`.

    Finds the event probability `p` such that the sum of the terms 0 through
    `k` of the binomial probability density is equal to the given cumulative
    probability `y`.

    Parameters
    ----------
    k : array_like
        Number of successes (float), rounded down to the nearest integer.
    n : array_like
        Number of events (float)
    y : array_like
        Cumulative probability (probability of `k` or fewer successes in `n`
        events).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    p : scalar or ndarray
        The event probability such that `bdtr(\lfloor k \rfloor, n, p) = y`.

    See Also
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
           https://netlib.org/cephes/
    """)

add_newdoc("bdtrik",
    """
    bdtrik(y, n, p, out=None)

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    k : scalar or ndarray
        The number of successes `k` such that `bdtr(k, n, p) = y`.

    See Also
    --------
    bdtr

    Notes
    -----
    Formula 26.5.24 of [1]_ (or equivalently [2]_) is used to reduce the binomial
    distribution to the cumulative incomplete beta distribution.

    Computation of `k` involves a search for a value that produces the desired
    value of `y`. The search relies on the monotonicity of `y` with `k`.

    Wrapper for the CDFLIB [3]_ Fortran routine `cdfbin`.

    References
    ----------
    .. [1] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [2] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17.5#E5
    .. [3] Barry Brown, James Lovato, and Kathy Russell,
           CDFLIB: Library of Fortran Routines for Cumulative Distribution
           Functions, Inverses, and Other Parameters.

    """)

add_newdoc("bdtrin",
    r"""
    bdtrin(k, y, p, out=None)

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    n : scalar or ndarray
        The number of events `n` such that `bdtr(k, n, p) = y`.

    See Also
    --------
    bdtr

    Notes
    -----
    This function uses the `find_minimum_number_of_trials` method of the
    `binomial_distribution` class of the Boost.Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    How often do we have to flip a fair coin to have at least a 90% chance
    of getting 10 heads? `bdtrin` answers this question:

    >>> import scipy.special as sc
    >>> k = 10  # number of times we want heads
    >>> p = 0.5  # probability of flipping heads
    >>> y = 0.9  # cumulative probability
    >>> result = sc.bdtrin(k, y, p)
    >>> result
    15.90442928275109

    To verify, compute the cumulative probability of getting 10 or fewer
    successes in 16 trials with probability 0.5 using the binomial
    distribution from `scipy.stats`. Since `bdtrin` returns a non-integer
    number of trials, we round up to the next integer:

    >>> from scipy.stats import Binomial
    >>> Binomial(n=16, p=p).cdf(k)
    0.8949432373046875

    """)

add_newdoc("btdtria",
    r"""
    btdtria(p, b, x, out=None)

    Inverse of `betainc` with respect to `a`.

    This is the inverse of the beta cumulative distribution function, `betainc`,
    considered as a function of `a`, returning the value of `a` for which
    `betainc(a, b, x) = p`, or

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    a : scalar or ndarray
        The value of the shape parameter `a` such that `betainc(a, b, x) = p`.

    See Also
    --------
    betainc : Regularized incomplete beta function
    betaincinv : Inverse of the regularized incomplete beta function
    btdtrib : Inverse of the beta cumulative distribution function, with respect to `b`.

    Notes
    -----
    This function wraps the ``ibeta_inva`` routine from the
    Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import scipy.special as sc

    This function is the inverse of `betainc` for fixed
    values of :math:`b` and :math:`x`.

    >>> a, b, x = 1.2, 3.1, 0.2
    >>> y = sc.betainc(a, b, x)
    >>> sc.btdtria(y, b, x)
    1.2

    """)

add_newdoc("btdtrib",
    r"""
    btdtria(a, p, x, out=None)

    Inverse of `betainc` with respect to `b`.

    This is the inverse of the beta cumulative distribution function, `betainc`,
    considered as a function of `b`, returning the value of `b` for which
    `betainc(a, b, x) = p`, or

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    b : scalar or ndarray
        The value of the shape parameter `b` such that `betainc(a, b, x) = p`.

    See Also
    --------
    betainc : Regularized incomplete beta function
    betaincinv : Inverse of the regularized incomplete beta function with
                 respect to `x`.
    btdtria : Inverse of the beta cumulative distribution function, with respect to `a`.

    Notes
    -----
    Wrapper for the `ibeta_invb` routine from the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import scipy.special as sc
    >>> a, b, x = 1.2, 3.1, 0.2
    >>> y = sc.betainc(a, b, x)

    `btdtrib` is the inverse of `betainc` for fixed values of :math:`a` and
    :math:`x`:

    >>> sc.btdtrib(a, y, x)
    3.1

    """)

add_newdoc(
    "betainc",
    r"""
    betainc(a, b, x, out=None)

    Regularized incomplete beta function.

    Computes the regularized incomplete beta function, defined as [1]_:

    .. math::

        I_x(a, b) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \int_0^x
        t^{a-1}(1-t)^{b-1}dt,

    for :math:`0 \leq x \leq 1`.

    This function is the cumulative distribution function for the beta
    distribution; its range is [0, 1].

    Parameters
    ----------
    a, b : array_like
           Positive, real-valued parameters.
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    scalar or ndarray
        Value of the regularized incomplete beta function.

    See Also
    --------
    beta : beta function
    betaincinv : inverse of the regularized incomplete beta function
    betaincc : complement of the regularized incomplete beta function
    scipy.stats.beta : beta distribution

    Notes
    -----
    The term *regularized* in the name of this function refers to the
    scaling of the function by the gamma function terms shown in the
    formula.  When not qualified as *regularized*, the name *incomplete
    beta function* often refers to just the integral expression,
    without the gamma terms.  One can use the function `beta` from
    `scipy.special` to get this "nonregularized" incomplete beta
    function by multiplying the result of ``betainc(a, b, x)`` by
    ``beta(a, b)``.

    ``betainc(a, b, x)`` is treated as a two parameter family of functions
    of a single variable `x`, rather than as a function of three variables.
    This impacts only the limiting cases ``a = 0``, ``b = 0``, ``a = inf``,
    ``b = inf``.

    In general

    .. math::

        \lim_{(a, b) \rightarrow (a_0, b_0)} \mathrm{betainc}(a, b, x)

    is treated as a pointwise limit in ``x``. Thus for example,
    ``betainc(0, b, 0)`` equals ``0`` for ``b > 0``, although it would be
    indeterminate when considering the simultaneous limit ``(a, x) -> (0+, 0+)``.

    This function wraps the ``ibeta`` routine from the
    Boost Math C++ library [2]_.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17
    .. [2] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------

    Let :math:`B(a, b)` be the `beta` function.

    >>> import scipy.special as sc

    The coefficient in terms of `gamma` is equal to
    :math:`1/B(a, b)`. Also, when :math:`x=1`
    the integral is equal to :math:`B(a, b)`.
    Therefore, :math:`I_{x=1}(a, b) = 1` for any :math:`a, b`.

    >>> sc.betainc(0.2, 3.5, 1.0)
    1.0

    It satisfies
    :math:`I_x(a, b) = x^a F(a, 1-b, a+1, x)/ (aB(a, b))`,
    where :math:`F` is the hypergeometric function `hyp2f1`:

    >>> a, b, x = 1.4, 3.1, 0.5
    >>> x**a * sc.hyp2f1(a, 1 - b, a + 1, x)/(a * sc.beta(a, b))
    0.8148904036225295
    >>> sc.betainc(a, b, x)
    0.8148904036225296

    This functions satisfies the relationship
    :math:`I_x(a, b) = 1 - I_{1-x}(b, a)`:

    >>> sc.betainc(2.2, 3.1, 0.4)
    0.49339638807619446
    >>> 1 - sc.betainc(3.1, 2.2, 1 - 0.4)
    0.49339638807619446

    """)


add_newdoc(
    "betaincc",
    r"""
    betaincc(a, b, x, out=None)

    Complement of the regularized incomplete beta function.

    Computes the complement of the regularized incomplete beta function,
    defined as [1]_:

    .. math::

        \bar{I}_x(a, b) = 1 - I_x(a, b)
                        = 1 - \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \int_0^x
                                  t^{a-1}(1-t)^{b-1}dt,

    for :math:`0 \leq x \leq 1`.

    Parameters
    ----------
    a, b : array_like
           Positive, real-valued parameters.
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    scalar or ndarray
        Value of the complement of the regularized incomplete beta function.

    See Also
    --------
    betainc : regularized incomplete beta function
    betaincinv : inverse of the regularized incomplete beta function
    betainccinv :
        inverse of the complement of the regularized incomplete beta function
    beta : beta function
    scipy.stats.beta : beta distribution

    Notes
    -----
    .. versionadded:: 1.11.0

    Like `betainc`, ``betaincc(a, b, x)`` is treated as a two parameter
    family of functions of a single variable `x`, rather than as a function of
    three variables. See the `betainc` docstring for more info on how this
    impacts limiting cases.

    This function wraps the ``ibetac`` routine from the
    Boost Math C++ library [2]_.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17
    .. [2] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import betaincc, betainc

    The naive calculation ``1 - betainc(a, b, x)`` loses precision when
    the values of ``betainc(a, b, x)`` are close to 1:

    >>> 1 - betainc(0.5, 8, [0.9, 0.99, 0.999])
    array([2.0574632e-09, 0.0000000e+00, 0.0000000e+00])

    By using ``betaincc``, we get the correct values:

    >>> betaincc(0.5, 8, [0.9, 0.99, 0.999])
    array([2.05746321e-09, 1.97259354e-17, 1.96467954e-25])

    """)

add_newdoc(
    "betaincinv",
    r"""
    betaincinv(a, b, y, out=None)

    Inverse of the regularized incomplete beta function.

    Computes :math:`x` such that:

    .. math::

        y = I_x(a, b) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}
        \int_0^x t^{a-1}(1-t)^{b-1}dt,

    where :math:`I_x` is the normalized incomplete beta function `betainc`
    and :math:`\Gamma` is the `gamma` function [1]_.

    Parameters
    ----------
    a, b : array_like
        Positive, real-valued parameters
    y : array_like
        Real-valued input
    out : ndarray, optional
        Optional output array for function values

    Returns
    -------
    scalar or ndarray
        Value of the inverse of the regularized incomplete beta function

    See Also
    --------
    betainc : regularized incomplete beta function
    gamma : gamma function

    Notes
    -----
    This function wraps the ``ibeta_inv`` routine from the
    Boost Math C++ library [2]_.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17
    .. [2] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import scipy.special as sc

    This function is the inverse of `betainc` for fixed
    values of :math:`a` and :math:`b`.

    >>> a, b = 1.2, 3.1
    >>> y = sc.betainc(a, b, 0.2)
    >>> sc.betaincinv(a, b, y)
    0.2
    >>>
    >>> a, b = 7.5, 0.4
    >>> x = sc.betaincinv(a, b, 0.5)
    >>> sc.betainc(a, b, x)
    0.5

    """)


add_newdoc(
    "betainccinv",
    r"""
    betainccinv(a, b, y, out=None)

    Inverse of the complemented regularized incomplete beta function.

    Computes :math:`x` such that:

    .. math::

        y = 1 - I_x(a, b) = 1 - \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}
        \int_0^x t^{a-1}(1-t)^{b-1}dt,

    where :math:`I_x` is the normalized incomplete beta function `betainc`
    and :math:`\Gamma` is the `gamma` function [1]_.

    Parameters
    ----------
    a, b : array_like
        Positive, real-valued parameters
    y : array_like
        Real-valued input
    out : ndarray, optional
        Optional output array for function values

    Returns
    -------
    scalar or ndarray
        Value of the inverse of the regularized incomplete beta function

    See Also
    --------
    betainc : regularized incomplete beta function
    betaincc : complement of the regularized incomplete beta function

    Notes
    -----
    .. versionadded:: 1.11.0

    This function wraps the ``ibetac_inv`` routine from the
    Boost Math C++ library [2]_.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17
    .. [2] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import betainccinv, betaincc

    This function is the inverse of `betaincc` for fixed
    values of :math:`a` and :math:`b`.

    >>> a, b = 1.2, 3.1
    >>> y = betaincc(a, b, 0.2)
    >>> betainccinv(a, b, y)
    0.2

    >>> a, b = 7, 2.5
    >>> x = betainccinv(a, b, 0.875)
    >>> betaincc(a, b, x)
    0.875

    """)

add_newdoc("chdtriv",
    """
    chdtriv(p, x, out=None)

    Inverse to `chdtr` with respect to `v`.

    Returns `v` such that ``chdtr(v, x) == p``.

    Parameters
    ----------
    p : array_like
        Probability that the Chi square random variable is less than
        or equal to `x`.
    x : array_like
        Nonnegative input.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Degrees of freedom.

    See Also
    --------
    chdtr, chdtrc, chdtri

    Notes
    -----
    This function wraps routines from the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.
    .. [2] Chi-Square distribution,
        https://www.itl.nist.gov/div898/handbook/eda/section3/eda3666.htm

    Examples
    --------
    >>> import scipy.special as sc

    It inverts `chdtr`.

    >>> p, x = 0.5, 1
    >>> sc.chdtr(sc.chdtriv(p, x), x)
    0.5000000000000003
    >>> v = 1
    >>> sc.chdtriv(sc.chdtr(v, x), x)
    1.0

    """)

add_newdoc("chndtr",
    r"""
    chndtr(x, df, nc, out=None)

    Non-central chi-squared cumulative distribution function.

    The cumulative distribution function is given by

    .. math::

        F_{\nu,\lambda}(x)
        = \sum_{j=0}^{\infty}
          e^{-\lambda / 2}
          \frac{(\lambda / 2)^j}{j!}
          F_{\chi^2_{\nu + 2j}}(x),

    where :math:`\nu > 0` is the degrees of freedom (``df``), :math:`\lambda \geq 0`
    is the non-centrality parameter (``nc``), and :math:`F_{\chi^2_{\nu + 2j}}` is the
    CDF of the central chi-squared distribution with :math:`\nu + 2j` degrees of
    freedom.

    Parameters
    ----------
    x : array_like
        Upper bound of the integral; must satisfy ``x >= 0``.
    df : array_like
        Degrees of freedom; must satisfy ``df > 0``.
    nc : array_like
        Non-centrality parameter; must satisfy ``nc >= 0``.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    cdf : scalar or ndarray
        Value of the non-central chi-squared cumulative distribution function.

    See Also
    --------
    chndtrix: Noncentral Chi Squared distribution quantile
    chndtridf: Inverse of `chndtr` with respect to `df`
    chndtrinc: Inverse of `chndtr` with respect to `nc`
    scipy.stats.ncx2: Non-central chi-squared distribution

    Notes
    -----
    The noncentral chi squared distribution is also available in
    `scipy.stats.ncx2`. ``scipy.stats.ncx2.cdf`` is equivalent to `chndtr`.

    This function wraps routines from the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    Compute the noncentral chi squared distribution CDF at one point.

    >>> x = 4.0
    >>> df = 1.0
    >>> nc = 5.0
    >>> sc.chndtr(x, df, nc)
    0.40667858759710945

    Plot the noncentral chi squared distribution CDF for different parameters.

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(0, 40, 1000)
    >>> plt.plot(x, sc.chndtr(x, 1, 5), label=r"$df=1,\ nc=5$")
    >>> plt.plot(x, sc.chndtr(x, 5, 10), label=r"$df=5,\ nc=10$")
    >>> plt.legend()
    >>> plt.show()

    """)

add_newdoc("chndtrix",
    """
    chndtrix(p, df, nc, out=None)

    Inverse to `chndtr` vs `x`.

    Calculated using a search to find a value for `x` that produces the
    desired value of `p`.

    Parameters
    ----------
    p : array_like
        Probability; must satisfy ``0 <= p < 1``
    df : array_like
        Degrees of freedom; must satisfy ``df > 0``
    nc : array_like
        Non-centrality parameter; must satisfy ``nc >= 0``
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    x : scalar or ndarray
        Value so that the probability a non-central Chi square random variable
        with `df` degrees of freedom and non-centrality, `nc`, is greater than
        `x` equals `p`.

    See Also
    --------
    chndtr : Noncentral chi-squared distribution CDF
    chndtridf : inverse of `chndtr` with respect to `cdf`
    chndtrinc : inverse of `chndtr` with respect to `nc`
    scipy.stats.ncx2 : Non-central chi-squared distribution

    Notes
    -----
    The noncentral chi squared distribution is also available in
    `scipy.stats.ncx2`. ``scipy.stats.ncx2.ppf`` is equivalent to `chndtrix`.

    This function wraps routines from the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import chndtrix, chndtr

    Compute the noncentral chi squared distribution CDF at one point.
    >>> x, df, nc = 3, 5, 10
    >>> p = chndtr(x, df, nc)

    `chndtrix` is the inverse of `chndtr` with respect to `x`:

    >>> chndtrix(p, df, nc)
    3.0

    """)

add_newdoc("chndtridf",
    """
    chndtridf(x, p, nc, out=None)

    Inverse to `chndtr` vs `df`.

    Calculated using a search to find a value for `df` that produces the
    desired value of `p`.

    Parameters
    ----------
    x : array_like
        Upper bound of the integral; must satisfy ``x >= 0``
    p : array_like
        Probability; must satisfy ``0 <= p < 1``
    nc : array_like
        Non-centrality parameter; must satisfy ``nc >= 0``
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    df : scalar or ndarray
        Degrees of freedom

    See Also
    --------
    chndtr : Noncentral chi-squared distribution CDF
    chndtrix : inverse of `chndtr` with respect to `x`
    chndtrinc : inverse of `chndtr` with respect to `nc`
    scipy.stats.ncx2 : Non-central chi-squared distribution

    Notes
    -----
    The noncentral chi squared distribution is also available in
    `scipy.stats.ncx2`.

    This function wraps routines from the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import chndtridf, chndtr

    Compute the noncentral chi squared distribution CDF at one point.

    >>> x, df, nc = 3, 5, 10
    >>> p = chndtr(x, df, nc)

    `chndtridf` is the inverse of `chndtr` with respect to `df`:

    >>> chndtridf(x, p, nc)
    5.0

    """)

add_newdoc("chndtrinc",
    """
    chndtrinc(x, df, p, out=None)

    Inverse of `chndtr` with respect to `nc`.

    Finds the non-centrality parameter `nc` such that

    .. math::

        \\operatorname{chndtr}(x, df, nc) = p.

    Parameters
    ----------
    x : array_like
        Upper bound of the integral; must satisfy ``x >= 0``.
    df : array_like
        Degrees of freedom; must satisfy ``df > 0``.
    p : array_like
        Probability; must satisfy ``0 <= p < 1``.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    nc : scalar or ndarray
        Non-centrality parameter.

    See Also
    --------
    chndtr : Noncentral chi-squared distribution CDF
    chndtridf : Inverse of `chndtr` with respect to `df`
    chndtrinc : Inverse of `chndtr` with respect to `nc`
    scipy.stats.ncx2 : Non-central chi-squared distribution

    Notes
    -----
    The noncentral chi squared distribution is also available in
    `scipy.stats.ncx2`.

    This function wraps routines from the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import chndtrinc, chndtr

    Compute the noncentral chi squared distribution CDF at one point.

    >>> x, df, nc = 3, 5, 10
    >>> p = chndtr(x, df, nc)

    `chndtrinc` is the inverse of `chndtr` with respect to `nc`:

    >>> chndtrinc(x, df, p)
    10.0

    """)

add_newdoc(
    "elliprc",
    r"""
    elliprc(x, y, out=None)

    Degenerate symmetric elliptic integral.

    The function RC is defined as [1]_

    .. math::

        R_{\mathrm{C}}(x, y) =
           \frac{1}{2} \int_0^{+\infty} (t + x)^{-1/2} (t + y)^{-1} dt
           = R_{\mathrm{F}}(x, y, y)

    Parameters
    ----------
    x, y : array_like
        Real or complex input parameters. `x` can be any number in the
        complex plane cut along the negative real axis. `y` must be non-zero.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    R : scalar or ndarray
        Value of the integral. If `y` is real and negative, the Cauchy
        principal value is returned. If both of `x` and `y` are real, the
        return value is real. Otherwise, the return value is complex.

    See Also
    --------
    elliprf : Completely-symmetric elliptic integral of the first kind.
    elliprd : Symmetric elliptic integral of the second kind.
    elliprg : Completely-symmetric elliptic integral of the second kind.
    elliprj : Symmetric elliptic integral of the third kind.

    Notes
    -----
    RC is a degenerate case of the symmetric integral RF: ``elliprc(x, y) ==
    elliprf(x, y, y)``. It is an elementary function rather than an elliptic
    integral.

    The code implements Carlson's algorithm based on the duplication theorems
    and series expansion up to the 7th order. [2]_

    .. versionadded:: 1.8.0

    References
    ----------
    .. [1] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
           Functions," NIST, US Dept. of Commerce.
           https://dlmf.nist.gov/19.16.E6
    .. [2] B. C. Carlson, "Numerical computation of real or complex elliptic
           integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
           :doi:`10.1007/BF02198293`. https://arxiv.org/abs/math/9409227

    Examples
    --------
    Basic homogeneity property:

    >>> import numpy as np
    >>> from scipy.special import elliprc

    >>> x = 1.2 + 3.4j
    >>> y = 5.
    >>> scale = 0.3 + 0.4j
    >>> elliprc(scale*x, scale*y)
    (0.5484493976710874-0.4169557678995833j)

    >>> elliprc(x, y)/np.sqrt(scale)
    (0.5484493976710874-0.41695576789958333j)

    When the two arguments coincide, the integral is particularly
    simple:

    >>> x = 1.2 + 3.4j
    >>> elliprc(x, x)
    (0.4299173120614631-0.3041729818745595j)

    >>> 1/np.sqrt(x)
    (0.4299173120614631-0.30417298187455954j)

    Another simple case: the first argument vanishes:

    >>> y = 1.2 + 3.4j
    >>> elliprc(0, y)
    (0.6753125346116815-0.47779380263880866j)

    >>> np.pi/2/np.sqrt(y)
    (0.6753125346116815-0.4777938026388088j)

    When `x` and `y` are both positive, we can express
    :math:`R_C(x,y)` in terms of more elementary functions.  For the
    case :math:`0 \le x < y`,

    >>> x = 3.2
    >>> y = 6.
    >>> elliprc(x, y)
    0.44942991498453444

    >>> np.arctan(np.sqrt((y-x)/x))/np.sqrt(y-x)
    0.44942991498453433

    And for the case :math:`0 \le y < x`,

    >>> x = 6.
    >>> y = 3.2
    >>> elliprc(x,y)
    0.4989837501576147

    >>> np.log((np.sqrt(x)+np.sqrt(x-y))/np.sqrt(y))/np.sqrt(x-y)
    0.49898375015761476

    """)

add_newdoc(
    "elliprd",
    r"""
    elliprd(x, y, z, out=None)

    Symmetric elliptic integral of the second kind.

    The function RD is defined as [1]_

    .. math::

        R_{\mathrm{D}}(x, y, z) =
           \frac{3}{2} \int_0^{+\infty} [(t + x) (t + y)]^{-1/2} (t + z)^{-3/2}
           dt

    Parameters
    ----------
    x, y, z : array_like
        Real or complex input parameters. `x` or `y` can be any number in the
        complex plane cut along the negative real axis, but at most one of them
        can be zero, while `z` must be non-zero.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    R : scalar or ndarray
        Value of the integral. If all of `x`, `y`, and `z` are real, the
        return value is real. Otherwise, the return value is complex.

    See Also
    --------
    elliprc : Degenerate symmetric elliptic integral.
    elliprf : Completely-symmetric elliptic integral of the first kind.
    elliprg : Completely-symmetric elliptic integral of the second kind.
    elliprj : Symmetric elliptic integral of the third kind.

    Notes
    -----
    RD is a degenerate case of the elliptic integral RJ: ``elliprd(x, y, z) ==
    elliprj(x, y, z, z)``.

    The code implements Carlson's algorithm based on the duplication theorems
    and series expansion up to the 7th order. [2]_

    .. versionadded:: 1.8.0

    References
    ----------
    .. [1] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
           Functions," NIST, US Dept. of Commerce.
           https://dlmf.nist.gov/19.16.E5
    .. [2] B. C. Carlson, "Numerical computation of real or complex elliptic
           integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
           :doi:`10.1007/BF02198293`. https://arxiv.org/abs/math/9409227

    Examples
    --------
    Basic homogeneity property:

    >>> import numpy as np
    >>> from scipy.special import elliprd

    >>> x = 1.2 + 3.4j
    >>> y = 5.
    >>> z = 6.
    >>> scale = 0.3 + 0.4j
    >>> elliprd(scale*x, scale*y, scale*z)
    (-0.03703043835680379-0.24500934665683802j)

    >>> elliprd(x, y, z)*np.power(scale, -1.5)
    (-0.0370304383568038-0.24500934665683805j)

    All three arguments coincide:

    >>> x = 1.2 + 3.4j
    >>> elliprd(x, x, x)
    (-0.03986825876151896-0.14051741840449586j)

    >>> np.power(x, -1.5)
    (-0.03986825876151894-0.14051741840449583j)

    The so-called "second lemniscate constant":

    >>> elliprd(0, 2, 1)/3
    0.5990701173677961

    >>> from scipy.special import gamma
    >>> gamma(0.75)**2/np.sqrt(2*np.pi)
    0.5990701173677959

    """)

add_newdoc(
    "elliprf",
    r"""
    elliprf(x, y, z, out=None)

    Completely-symmetric elliptic integral of the first kind.

    The function RF is defined as [1]_

    .. math::

        R_{\mathrm{F}}(x, y, z) =
           \frac{1}{2} \int_0^{+\infty} [(t + x) (t + y) (t + z)]^{-1/2} dt

    Parameters
    ----------
    x, y, z : array_like
        Real or complex input parameters. `x`, `y`, or `z` can be any number in
        the complex plane cut along the negative real axis, but at most one of
        them can be zero.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    R : scalar or ndarray
        Value of the integral. If all of `x`, `y`, and `z` are real, the return
        value is real. Otherwise, the return value is complex.

    See Also
    --------
    elliprc : Degenerate symmetric integral.
    elliprd : Symmetric elliptic integral of the second kind.
    elliprg : Completely-symmetric elliptic integral of the second kind.
    elliprj : Symmetric elliptic integral of the third kind.

    Notes
    -----
    The code implements Carlson's algorithm based on the duplication theorems
    and series expansion up to the 7th order (cf.:
    https://dlmf.nist.gov/19.36.i) and the AGM algorithm for the complete
    integral. [2]_

    .. versionadded:: 1.8.0

    References
    ----------
    .. [1] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
           Functions," NIST, US Dept. of Commerce.
           https://dlmf.nist.gov/19.16.E1
    .. [2] B. C. Carlson, "Numerical computation of real or complex elliptic
           integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
           :doi:`10.1007/BF02198293`. https://arxiv.org/abs/math/9409227

    Examples
    --------
    Basic homogeneity property:

    >>> import numpy as np
    >>> from scipy.special import elliprf

    >>> x = 1.2 + 3.4j
    >>> y = 5.
    >>> z = 6.
    >>> scale = 0.3 + 0.4j
    >>> elliprf(scale*x, scale*y, scale*z)
    (0.5328051227278146-0.4008623567957094j)

    >>> elliprf(x, y, z)/np.sqrt(scale)
    (0.5328051227278147-0.4008623567957095j)

    All three arguments coincide:

    >>> x = 1.2 + 3.4j
    >>> elliprf(x, x, x)
    (0.42991731206146316-0.30417298187455954j)

    >>> 1/np.sqrt(x)
    (0.4299173120614631-0.30417298187455954j)

    The so-called "first lemniscate constant":

    >>> elliprf(0, 1, 2)
    1.3110287771460598

    >>> from scipy.special import gamma
    >>> gamma(0.25)**2/(4*np.sqrt(2*np.pi))
    1.3110287771460598

    """)

add_newdoc(
    "elliprg",
    r"""
    elliprg(x, y, z, out=None)

    Completely-symmetric elliptic integral of the second kind.

    The function RG is defined as [1]_

    .. math::

        R_{\mathrm{G}}(x, y, z) =
           \frac{1}{4} \int_0^{+\infty} [(t + x) (t + y) (t + z)]^{-1/2}
           \left(\frac{x}{t + x} + \frac{y}{t + y} + \frac{z}{t + z}\right) t
           dt

    Parameters
    ----------
    x, y, z : array_like
        Real or complex input parameters. `x`, `y`, or `z` can be any number in
        the complex plane cut along the negative real axis.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    R : scalar or ndarray
        Value of the integral. If all of `x`, `y`, and `z` are real, the return
        value is real. Otherwise, the return value is complex.

    See Also
    --------
    elliprc : Degenerate symmetric integral.
    elliprd : Symmetric elliptic integral of the second kind.
    elliprf : Completely-symmetric elliptic integral of the first kind.
    elliprj : Symmetric elliptic integral of the third kind.

    Notes
    -----
    The implementation uses the relation [1]_

    .. math::

        2 R_{\mathrm{G}}(x, y, z) =
           z R_{\mathrm{F}}(x, y, z) -
           \frac{1}{3} (x - z) (y - z) R_{\mathrm{D}}(x, y, z) +
           \sqrt{\frac{x y}{z}}

    and the symmetry of `x`, `y`, `z` when at least one non-zero parameter can
    be chosen as the pivot. When one of the arguments is close to zero, the AGM
    method is applied instead. Other special cases are computed following Ref.
    [2]_

    .. versionadded:: 1.8.0

    References
    ----------
    .. [1] B. C. Carlson, "Numerical computation of real or complex elliptic
           integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
           :doi:`10.1007/BF02198293`. https://arxiv.org/abs/math/9409227
    .. [2] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
           Functions," NIST, US Dept. of Commerce.
           https://dlmf.nist.gov/19.16.E1
           https://dlmf.nist.gov/19.20.ii

    Examples
    --------
    Basic homogeneity property:

    >>> import numpy as np
    >>> from scipy.special import elliprg

    >>> x = 1.2 + 3.4j
    >>> y = 5.
    >>> z = 6.
    >>> scale = 0.3 + 0.4j
    >>> elliprg(scale*x, scale*y, scale*z)
    (1.195936862005246+0.8470988320464167j)

    >>> elliprg(x, y, z)*np.sqrt(scale)
    (1.195936862005246+0.8470988320464165j)

    Simplifications:

    >>> elliprg(0, y, y)
    1.756203682760182

    >>> 0.25*np.pi*np.sqrt(y)
    1.7562036827601817

    >>> elliprg(0, 0, z)
    1.224744871391589

    >>> 0.5*np.sqrt(z)
    1.224744871391589

    The surface area of a triaxial ellipsoid with semiaxes ``a``, ``b``, and
    ``c`` is given by

    .. math::

        S = 4 \pi a b c R_{\mathrm{G}}(1 / a^2, 1 / b^2, 1 / c^2).

    >>> def ellipsoid_area(a, b, c):
    ...     r = 4.0 * np.pi * a * b * c
    ...     return r * elliprg(1.0 / (a * a), 1.0 / (b * b), 1.0 / (c * c))
    >>> print(ellipsoid_area(1, 3, 5))
    108.62688289491807
    """)

add_newdoc(
    "elliprj",
    r"""
    elliprj(x, y, z, p, out=None)

    Symmetric elliptic integral of the third kind.

    The function RJ is defined as [1]_

    .. math::

        R_{\mathrm{J}}(x, y, z, p) =
           \frac{3}{2} \int_0^{+\infty} [(t + x) (t + y) (t + z)]^{-1/2}
           (t + p)^{-1} dt

    .. warning::
        This function should be considered experimental when the inputs are
        unbalanced.  Check correctness with another independent implementation.

    Parameters
    ----------
    x, y, z, p : array_like
        Real or complex input parameters. `x`, `y`, or `z` are numbers in
        the complex plane cut along the negative real axis (subject to further
        constraints, see Notes), and at most one of them can be zero. `p` must
        be non-zero.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    R : scalar or ndarray
        Value of the integral. If all of `x`, `y`, `z`, and `p` are real, the
        return value is real. Otherwise, the return value is complex.

        If `p` is real and negative, while `x`, `y`, and `z` are real,
        non-negative, and at most one of them is zero, the Cauchy principal
        value is returned. [1]_ [2]_

    See Also
    --------
    elliprc : Degenerate symmetric integral.
    elliprd : Symmetric elliptic integral of the second kind.
    elliprf : Completely-symmetric elliptic integral of the first kind.
    elliprg : Completely-symmetric elliptic integral of the second kind.

    Notes
    -----
    The code implements Carlson's algorithm based on the duplication theorems
    and series expansion up to the 7th order. [3]_ The algorithm is slightly
    different from its earlier incarnation as it appears in [1]_, in that the
    call to `elliprc` (or ``atan``/``atanh``, see [4]_) is no longer needed in
    the inner loop. Asymptotic approximations are used where arguments differ
    widely in the order of magnitude. [5]_

    The input values are subject to certain sufficient but not necessary
    constraints when input arguments are complex. Notably, ``x``, ``y``, and
    ``z`` must have non-negative real parts, unless two of them are
    non-negative and complex-conjugates to each other while the other is a real
    non-negative number. [1]_ If the inputs do not satisfy the sufficient
    condition described in Ref. [1]_ they are rejected outright with the output
    set to NaN.

    In the case where one of ``x``, ``y``, and ``z`` is equal to ``p``, the
    function ``elliprd`` should be preferred because of its less restrictive
    domain.

    .. versionadded:: 1.8.0

    References
    ----------
    .. [1] B. C. Carlson, "Numerical computation of real or complex elliptic
           integrals," Numer. Algorithm, vol. 10, no. 1, pp. 13-26, 1995.
           :doi:`10.1007/BF02198293`. https://arxiv.org/abs/math/9409227
    .. [2] B. C. Carlson, ed., Chapter 19 in "Digital Library of Mathematical
           Functions," NIST, US Dept. of Commerce.
           https://dlmf.nist.gov/19.20.iii
    .. [3] B. C. Carlson, J. FitzSimmons, "Reduction Theorems for Elliptic
           Integrands with the Square Root of Two Quadratic Factors," J.
           Comput. Appl. Math., vol. 118, nos. 1-2, pp. 71-85, 2000.
           :doi:`10.1016/S0377-0427(00)00282-X`.
    .. [4] F. Johansson, "Numerical Evaluation of Elliptic Functions, Elliptic
           Integrals and Modular Forms," in J. Blumlein, C. Schneider, P.
           Paule, eds., "Elliptic Integrals, Elliptic Functions and Modular
           Forms in Quantum Field Theory," pp. 269-293, 2019 (Cham,
           Switzerland: Springer Nature Switzerland).
           :doi:`10.1007/978-3-030-04480-0`. https://arxiv.org/abs/1806.06725
    .. [5] B. C. Carlson, J. L. Gustafson, "Asymptotic Approximations for
           Symmetric Elliptic Integrals," SIAM J. Math. Anls., vol. 25, no. 2,
           pp. 288-303, 1994. :doi:`10.1137/S0036141092228477`.
           https://arxiv.org/abs/math/9310223

    Examples
    --------
    Basic homogeneity property:

    >>> import numpy as np
    >>> from scipy.special import elliprj

    >>> x = 1.2 + 3.4j
    >>> y = 5.
    >>> z = 6.
    >>> p = 7.
    >>> scale = 0.3 - 0.4j
    >>> elliprj(scale*x, scale*y, scale*z, scale*p)
    (0.10834905565679157+0.19694950747103812j)

    >>> elliprj(x, y, z, p)*np.power(scale, -1.5)
    (0.10834905565679556+0.19694950747103854j)

    Reduction to simpler elliptic integral:

    >>> elliprj(x, y, z, z)
    (0.08288462362195129-0.028376809745123258j)

    >>> from scipy.special import elliprd
    >>> elliprd(x, y, z)
    (0.08288462362195136-0.028376809745123296j)

    All arguments coincide:

    >>> elliprj(x, x, x, x)
    (-0.03986825876151896-0.14051741840449586j)

    >>> np.power(x, -1.5)
    (-0.03986825876151894-0.14051741840449583j)

    """)

add_newdoc(
    "erfinv",
    """
    erfinv(y, out=None)

    Inverse of the error function.

    Computes the inverse of the error function.

    In the complex domain, there is no unique complex number w satisfying
    erf(w)=z. This indicates a true inverse function would be multivalued.
    When the domain restricts to the real, -1 < x < 1, there is a unique real
    number satisfying erf(erfinv(x)) = x.

    Parameters
    ----------
    y : ndarray
        Argument at which to evaluate. Domain: [-1, 1]
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    erfinv : scalar or ndarray
        The inverse of erf of y, element-wise

    See Also
    --------
    erf : Error function of a complex argument
    erfc : Complementary error function, ``1 - erf(x)``
    erfcinv : Inverse of the complementary error function

    Notes
    -----
    This function wraps the ``erf_inv`` routine from the
    Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.special import erfinv, erf

    >>> erfinv(0.5)
    0.4769362762044699

    >>> y = np.linspace(-1.0, 1.0, num=9)
    >>> x = erfinv(y)
    >>> x
    array([       -inf, -0.81341985, -0.47693628, -0.22531206,  0.        ,
            0.22531206,  0.47693628,  0.81341985,         inf])

    Verify that ``erf(erfinv(y))`` is ``y``.

    >>> erf(x)
    array([-1.  , -0.75, -0.5 , -0.25,  0.  ,  0.25,  0.5 ,  0.75,  1.  ])

    Plot the function:

    >>> y = np.linspace(-1, 1, 200)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(y, erfinv(y))
    >>> ax.grid(True)
    >>> ax.set_xlabel('y')
    >>> ax.set_title('erfinv(y)')
    >>> plt.show()

    """)

add_newdoc("eval_jacobi",
    r"""
    eval_jacobi(n, alpha, beta, x, out=None)

    Evaluate Jacobi polynomial at a point.

    The Jacobi polynomials can be defined via the Gauss hypergeometric
    function :math:`{}_2F_1` as

    .. math::

        P_n^{(\alpha, \beta)}(x) = \frac{(\alpha + 1)_n}{\Gamma(n + 1)}
          {}_2F_1(-n, 1 + \alpha + \beta + n; \alpha + 1; (1 - x)/2)

    where :math:`(\cdot)_n` is the Pochhammer symbol; see `poch`. When
    :math:`n` is an integer the result is a polynomial of degree
    :math:`n`. See 22.5.42 in [AS]_ or [DLMF]_ for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer the result is
        determined via the relation to the Gauss hypergeometric
        function.
    alpha : array_like
        Parameter.
    beta : array_like
        Parameter.
    x : array_like
        Points at which to evaluate the polynomial.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    P : scalar or ndarray
        Values of the Jacobi polynomial.

    See Also
    --------
    roots_jacobi : roots and quadrature weights of Jacobi polynomials
    jacobi : Jacobi polynomial object
    hyp2f1 : Gauss hypergeometric function

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.E7

    """)

add_newdoc("eval_sh_jacobi",
    r"""
    eval_sh_jacobi(n, p, q, x, out=None)

    Evaluate shifted Jacobi polynomial at a point.

    Defined by

    .. math::

        G_n^{(p, q)}(x)
          = \binom{2n + p - 1}{n}^{-1} P_n^{(p - q, q - 1)}(2x - 1),

    where :math:`P_n^{(\cdot, \cdot)}` is the n-th Jacobi polynomial.
    See 22.5.2 in [AS]_ (or equivalently [DLMF]_)  for details.

    Parameters
    ----------
    n : int
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to `binom` and `eval_jacobi`.
    p : float
        Parameter
    q : float
        Parameter
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    G : scalar or ndarray
        Values of the shifted Jacobi polynomial.

    See Also
    --------
    roots_sh_jacobi : roots and quadrature weights of shifted Jacobi
                      polynomials
    sh_jacobi : shifted Jacobi polynomial object
    eval_jacobi : evaluate Jacobi polynomials

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.1.E2

    """)

add_newdoc("eval_gegenbauer",
    r"""
    eval_gegenbauer(n, alpha, x, out=None)

    Evaluate Gegenbauer (ultraspherical) polynomial at a point.

    The Gegenbauer polynomials can be defined via the Gauss
    hypergeometric function :math:`{}_2F_1` as

    .. math::

        C_n^{(\alpha)}(x) = \frac{(2\alpha)_n}{\Gamma(n + 1)}
          {}_2F_1(-n, 2\alpha + n; \alpha + 1/2; (1 - x)/2).

    where :math:`(\cdot)_n` is the Pochhammer symbol; see `poch`. When
    :math:`n` is an integer the result is a polynomial of degree
    :math:`n`. See 22.5.46 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to the Gauss hypergeometric
        function.
    alpha : array_like
        Parameter.
    x : array_like
        Points at which to evaluate the Gegenbauer polynomial.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    C : scalar or ndarray
        Values of the Gegenbauer polynomial.

    See Also
    --------
    roots_gegenbauer : roots and quadrature weights of Gegenbauer
                       polynomials
    gegenbauer : Gegenbauer polynomial object
    hyp2f1 : Gauss hypergeometric function

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.E9

    """)

add_newdoc("eval_chebyt",
    r"""
    eval_chebyt(n, x, out=None)

    Evaluate Chebyshev polynomial of the first kind at a point.

    The Chebyshev polynomials of the first kind can be defined via the
    Gauss hypergeometric function :math:`{}_2F_1` as

    .. math::

        T_n(x) = {}_2F_1(n, -n; 1/2; (1 - x)/2).

    When :math:`n` is an integer the result is a polynomial of degree
    :math:`n`. See 22.5.47 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to the Gauss hypergeometric
        function.
    x : array_like
        Points at which to evaluate the Chebyshev polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    T : scalar or ndarray
        Values of the Chebyshev polynomial

    See Also
    --------
    roots_chebyt : roots and quadrature weights of Chebyshev
                   polynomials of the first kind
    chebyu : Chebychev polynomial object
    eval_chebyu : evaluate Chebyshev polynomials of the second kind
    hyp2f1 : Gauss hypergeometric function
    numpy.polynomial.chebyshev.Chebyshev : Chebyshev series

    Notes
    -----
    This routine is numerically stable for `x` in ``[-1, 1]`` at least
    up to order ``10000``.

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.E11_2

    """)

add_newdoc("eval_chebyu",
    r"""
    eval_chebyu(n, x, out=None)

    Evaluate Chebyshev polynomial of the second kind at a point.

    The Chebyshev polynomials of the second kind can be defined via
    the Gauss hypergeometric function :math:`{}_2F_1` as

    .. math::

        U_n(x) = (n + 1) {}_2F_1(-n, n + 2; 3/2; (1 - x)/2).

    When :math:`n` is an integer the result is a polynomial of degree
    :math:`n`. See 22.5.48 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to the Gauss hypergeometric
        function.
    x : array_like
        Points at which to evaluate the Chebyshev polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    U : scalar or ndarray
        Values of the Chebyshev polynomial

    See Also
    --------
    roots_chebyu : roots and quadrature weights of Chebyshev
                   polynomials of the second kind
    chebyu : Chebyshev polynomial object
    eval_chebyt : evaluate Chebyshev polynomials of the first kind
    hyp2f1 : Gauss hypergeometric function

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.E11_4

    """)

add_newdoc("eval_chebys",
    r"""
    eval_chebys(n, x, out=None)

    Evaluate Chebyshev polynomial of the second kind on [-2, 2] at a
    point.

    These polynomials are defined as

    .. math::

        S_n(x) = U_n(x/2)

    where :math:`U_n` is a Chebyshev polynomial of the second kind.
    See 22.5.13 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to `eval_chebyu`.
    x : array_like
        Points at which to evaluate the Chebyshev polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    S : scalar or ndarray
        Values of the Chebyshev polynomial

    See Also
    --------
    roots_chebys : roots and quadrature weights of Chebyshev
                   polynomials of the second kind on [-2, 2]
    chebys : Chebyshev polynomial object
    eval_chebyu : evaluate Chebyshev polynomials of the second kind

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.1.E3

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    They are a scaled version of the Chebyshev polynomials of the
    second kind.

    >>> x = np.linspace(-2, 2, 6)
    >>> sc.eval_chebys(3, x)
    array([-4.   ,  0.672,  0.736, -0.736, -0.672,  4.   ])
    >>> sc.eval_chebyu(3, x / 2)
    array([-4.   ,  0.672,  0.736, -0.736, -0.672,  4.   ])

    """)

add_newdoc("eval_chebyc",
    r"""
    eval_chebyc(n, x, out=None)

    Evaluate Chebyshev polynomial of the first kind on [-2, 2] at a
    point.

    These polynomials are defined as

    .. math::

        C_n(x) = 2 T_n(x/2)

    where :math:`T_n` is a Chebyshev polynomial of the first kind. See
    22.5.11 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to `eval_chebyt`.
    x : array_like
        Points at which to evaluate the Chebyshev polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    C : scalar or ndarray
        Values of the Chebyshev polynomial

    See Also
    --------
    roots_chebyc : roots and quadrature weights of Chebyshev
                   polynomials of the first kind on [-2, 2]
    chebyc : Chebyshev polynomial object
    numpy.polynomial.chebyshev.Chebyshev : Chebyshev series
    eval_chebyt : evaluate Chebycshev polynomials of the first kind

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.1.E3

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    They are a scaled version of the Chebyshev polynomials of the
    first kind.

    >>> x = np.linspace(-2, 2, 6)
    >>> sc.eval_chebyc(3, x)
    array([-2.   ,  1.872,  1.136, -1.136, -1.872,  2.   ])
    >>> 2 * sc.eval_chebyt(3, x / 2)
    array([-2.   ,  1.872,  1.136, -1.136, -1.872,  2.   ])

    """)

add_newdoc("eval_sh_chebyt",
    r"""
    eval_sh_chebyt(n, x, out=None)

    Evaluate shifted Chebyshev polynomial of the first kind at a
    point.

    These polynomials are defined as

    .. math::

        T_n^*(x) = T_n(2x - 1)

    where :math:`T_n` is a Chebyshev polynomial of the first kind. See
    22.5.14 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to `eval_chebyt`.
    x : array_like
        Points at which to evaluate the shifted Chebyshev polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    T : scalar or ndarray
        Values of the shifted Chebyshev polynomial

    See Also
    --------
    roots_sh_chebyt : roots and quadrature weights of shifted
                      Chebyshev polynomials of the first kind
    sh_chebyt : shifted Chebyshev polynomial object
    eval_chebyt : evaluate Chebyshev polynomials of the first kind
    numpy.polynomial.chebyshev.Chebyshev : Chebyshev series

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.7.E7

    """)

add_newdoc("eval_sh_chebyu",
    r"""
    eval_sh_chebyu(n, x, out=None)

    Evaluate shifted Chebyshev polynomial of the second kind at a
    point.

    These polynomials are defined as

    .. math::

        U_n^*(x) = U_n(2x - 1)

    where :math:`U_n` is a Chebyshev polynomial of the first kind. See
    22.5.15 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to `eval_chebyu`.
    x : array_like
        Points at which to evaluate the shifted Chebyshev polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    U : scalar or ndarray
        Values of the shifted Chebyshev polynomial

    See Also
    --------
    roots_sh_chebyu : roots and quadrature weights of shifted
                      Chebychev polynomials of the second kind
    sh_chebyu : shifted Chebyshev polynomial object
    eval_chebyu : evaluate Chebyshev polynomials of the second kind

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.7.E8

    """)

add_newdoc("eval_legendre",
    r"""
    eval_legendre(n, x, out=None)

    Evaluate Legendre polynomial at a point.

    The Legendre polynomials can be defined via the Gauss
    hypergeometric function :math:`{}_2F_1` as

    .. math::

        P_n(x) = {}_2F_1(-n, n + 1; 1; (1 - x)/2).

    When :math:`n` is an integer the result is a polynomial of degree
    :math:`n`. See 22.5.49 in [AS]_ (or equivalently [DLMF]_) for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to the Gauss hypergeometric
        function.
    x : array_like
        Points at which to evaluate the Legendre polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    P : scalar or ndarray
        Values of the Legendre polynomial

    See Also
    --------
    roots_legendre : roots and quadrature weights of Legendre
                     polynomials
    legendre : Legendre polynomial object
    hyp2f1 : Gauss hypergeometric function
    numpy.polynomial.legendre.Legendre : Legendre series

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/15.9.E7

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import eval_legendre

    Evaluate the zero-order Legendre polynomial at x = 0

    >>> eval_legendre(0, 0)
    1.0

    Evaluate the first-order Legendre polynomial between -1 and 1

    >>> X = np.linspace(-1, 1, 5)  # Domain of Legendre polynomials
    >>> eval_legendre(1, X)
    array([-1. , -0.5,  0. ,  0.5,  1. ])

    Evaluate Legendre polynomials of order 0 through 4 at x = 0

    >>> N = range(0, 5)
    >>> eval_legendre(N, 0)
    array([ 1.   ,  0.   , -0.5  ,  0.   ,  0.375])

    Plot Legendre polynomials of order 0 through 4

    >>> X = np.linspace(-1, 1)

    >>> import matplotlib.pyplot as plt
    >>> for n in range(0, 5):
    ...     y = eval_legendre(n, X)
    ...     plt.plot(X, y, label=r'$P_{}(x)$'.format(n))

    >>> plt.title("Legendre Polynomials")
    >>> plt.xlabel("x")
    >>> plt.ylabel(r'$P_n(x)$')
    >>> plt.legend(loc='lower right')
    >>> plt.show()

    """)

add_newdoc("eval_sh_legendre",
    r"""
    eval_sh_legendre(n, x, out=None)

    Evaluate shifted Legendre polynomial at a point.

    These polynomials are defined as

    .. math::

        P_n^*(x) = P_n(2x - 1)

    where :math:`P_n` is a Legendre polynomial. See 2.2.11 in [AS]_
    or [DLMF]_ for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the value is
        determined via the relation to `eval_legendre`.
    x : array_like
        Points at which to evaluate the shifted Legendre polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    P : scalar or ndarray
        Values of the shifted Legendre polynomial

    See Also
    --------
    roots_sh_legendre : roots and quadrature weights of shifted
                        Legendre polynomials
    sh_legendre : shifted Legendre polynomial object
    eval_legendre : evaluate Legendre polynomials
    numpy.polynomial.legendre.Legendre : Legendre series

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.7.E10

    """)

add_newdoc("eval_genlaguerre",
    r"""
    eval_genlaguerre(n, alpha, x, out=None)

    Evaluate generalized Laguerre polynomial at a point.

    The generalized Laguerre polynomials can be defined via the
    confluent hypergeometric function :math:`{}_1F_1` as

    .. math::

        L_n^{(\alpha)}(x) = \binom{n + \alpha}{n}
          {}_1F_1(-n, \alpha + 1, x).

    When :math:`n` is an integer the result is a polynomial of degree
    :math:`n`. See 22.5.54 in [AS]_ or [DLMF]_ for details. The Laguerre
    polynomials are the special case where :math:`\alpha = 0`.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer, the result is
        determined via the relation to the confluent hypergeometric
        function.
    alpha : array_like
        Parameter; must have ``alpha > -1``
    x : array_like
        Points at which to evaluate the generalized Laguerre
        polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    L : scalar or ndarray
        Values of the generalized Laguerre polynomial

    See Also
    --------
    roots_genlaguerre : roots and quadrature weights of generalized
                        Laguerre polynomials
    genlaguerre : generalized Laguerre polynomial object
    hyp1f1 : confluent hypergeometric function
    eval_laguerre : evaluate Laguerre polynomials

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.E12

    """)

add_newdoc("eval_laguerre",
    r"""
    eval_laguerre(n, x, out=None)

    Evaluate Laguerre polynomial at a point.

    The Laguerre polynomials can be defined via the confluent
    hypergeometric function :math:`{}_1F_1` as

    .. math::

        L_n(x) = {}_1F_1(-n, 1, x).

    See 22.5.16 and 22.5.54 in [AS]_ (or equivalently [DLMF1]_ and [DLMF2]_)
    for details. When :math:`n` is an integer the result is a polynomial
    of degree :math:`n`.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial. If not an integer the result is
        determined via the relation to the confluent hypergeometric
        function.
    x : array_like
        Points at which to evaluate the Laguerre polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    L : scalar or ndarray
        Values of the Laguerre polynomial

    See Also
    --------
    roots_laguerre : roots and quadrature weights of Laguerre
                     polynomials
    laguerre : Laguerre polynomial object
    numpy.polynomial.laguerre.Laguerre : Laguerre series
    eval_genlaguerre : evaluate generalized Laguerre polynomials

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF1] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.1#I1.ix7.p1
    .. [DLMF2] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.E12

     """)

add_newdoc("eval_hermite",
    r"""
    eval_hermite(n, x, out=None)

    Evaluate physicist's Hermite polynomial at a point.

    Defined by

    .. math::

        H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2};

    :math:`H_n` is a polynomial of degree :math:`n`. See 22.11.7 in
    [AS]_ or [DLMF]_ for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial
    x : array_like
        Points at which to evaluate the Hermite polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    H : scalar or ndarray
        Values of the Hermite polynomial

    See Also
    --------
    roots_hermite : roots and quadrature weights of physicist's
                    Hermite polynomials
    hermite : physicist's Hermite polynomial object
    numpy.polynomial.hermite.Hermite : Physicist's Hermite series
    eval_hermitenorm : evaluate Probabilist's Hermite polynomials

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.T1

    """)

add_newdoc("eval_hermitenorm",
    r"""
    eval_hermitenorm(n, x, out=None)

    Evaluate probabilist's (normalized) Hermite polynomial at a
    point.

    Defined by

    .. math::

        He_n(x) = (-1)^n e^{x^2/2} \frac{d^n}{dx^n} e^{-x^2/2};

    :math:`He_n` is a polynomial of degree :math:`n`. See 22.11.8 in
    [AS]_ or [DLMF]_ for details.

    Parameters
    ----------
    n : array_like
        Degree of the polynomial
    x : array_like
        Points at which to evaluate the Hermite polynomial
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    He : scalar or ndarray
        Values of the Hermite polynomial

    See Also
    --------
    roots_hermitenorm : roots and quadrature weights of probabilist's
                        Hermite polynomials
    hermitenorm : probabilist's Hermite polynomial object
    numpy.polynomial.hermite_e.HermiteE : Probabilist's Hermite series
    eval_hermite : evaluate physicist's Hermite polynomials

    References
    ----------
    .. [AS] Milton Abramowitz and Irene A. Stegun, eds.
        Handbook of Mathematical Functions with Formulas,
        Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [DLMF] NIST Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/18.5.T1

    """)

add_newdoc("expn",
    r"""
    expn(n, x, out=None)

    Generalized exponential integral En.

    For integer :math:`n \geq 0` and real :math:`x \geq 0` the
    generalized exponential integral is defined as [DLMF]_

    .. math::

        E_n(x) = x^{n - 1} \int_x^\infty \frac{e^{-t}}{t^n} dt.

    Parameters
    ----------
    n : array_like
        Non-negative integers
    x : array_like
        Real argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the generalized exponential integral

    See Also
    --------
    exp1 : special case of :math:`E_n` for :math:`n = 1`
    expi : related to :math:`E_n` when :math:`n = 1`

    References
    ----------
    .. [DLMF] Digital Library of Mathematical Functions, 8.19.2
              https://dlmf.nist.gov/8.19#E2

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    Its domain is nonnegative n and x.

    >>> sc.expn(-1, 1.0), sc.expn(1, -1.0)
    (nan, nan)

    It has a pole at ``x = 0`` for ``n = 1, 2``; for larger ``n`` it
    is equal to ``1 / (n - 1)``.

    >>> sc.expn([0, 1, 2, 3, 4], 0)
    array([       inf,        inf, 1.        , 0.5       , 0.33333333])

    For n equal to 0 it reduces to ``exp(-x) / x``.

    >>> x = np.array([1, 2, 3, 4])
    >>> sc.expn(0, x)
    array([0.36787944, 0.06766764, 0.01659569, 0.00457891])
    >>> np.exp(-x) / x
    array([0.36787944, 0.06766764, 0.01659569, 0.00457891])

    For n equal to 1 it reduces to `exp1`.

    >>> sc.expn(1, x)
    array([0.21938393, 0.04890051, 0.01304838, 0.00377935])
    >>> sc.exp1(x)
    array([0.21938393, 0.04890051, 0.01304838, 0.00377935])

    """)

add_newdoc("fdtr",
    r"""
    fdtr(dfn, dfd, x, out=None)

    F cumulative distribution function.

    Returns the value of the cumulative distribution function of the
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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    y : scalar or ndarray
        The CDF of the F-distribution with parameters `dfn` and `dfd` at `x`.

    See Also
    --------
    fdtrc : F distribution survival function
    fdtri : F distribution inverse cumulative distribution
    scipy.stats.f : F distribution

    Notes
    -----
    The regularized incomplete beta function is used, according to the
    formula,

    .. math::
        F(d_n, d_d; x) = I_{xd_n/(d_d + xd_n)}(d_n/2, d_d/2).

    Wrapper for a routine from the Boost Math C++ library [1]_. The
    F distribution is also available as `scipy.stats.f`. Calling
    `fdtr` directly can improve performance compared to the ``cdf``
    method of `scipy.stats.f` (see last example below).

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.


    Examples
    --------
    Calculate the function for ``dfn=1`` and ``dfd=2`` at ``x=1``.

    >>> import numpy as np
    >>> from scipy.special import fdtr
    >>> fdtr(1, 2, 1)
    0.5773502691896258

    Calculate the function at several points by providing a NumPy array for
    `x`.

    >>> x = np.array([0.5, 2., 3.])
    >>> fdtr(1, 2, x)
    array([0.4472136 , 0.70710678, 0.77459667])

    Plot the function for several parameter sets.

    >>> import matplotlib.pyplot as plt
    >>> dfn_parameters = [1, 5, 10, 50]
    >>> dfd_parameters = [1, 1, 2, 3]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(dfn_parameters, dfd_parameters,
    ...                            linestyles))
    >>> x = np.linspace(0, 30, 1000)
    >>> fig, ax = plt.subplots()
    >>> for parameter_set in parameters_list:
    ...     dfn, dfd, style = parameter_set
    ...     fdtr_vals = fdtr(dfn, dfd, x)
    ...     ax.plot(x, fdtr_vals, label=rf"$d_n={dfn},\, d_d={dfd}$",
    ...             ls=style)
    >>> ax.legend()
    >>> ax.set_xlabel("$x$")
    >>> ax.set_title("F distribution cumulative distribution function")
    >>> plt.show()

    The F distribution is also available as `scipy.stats.f`. Using `fdtr`
    directly can be much faster than calling the ``cdf`` method of
    `scipy.stats.f`, especially for small arrays or individual values.
    To get the same results one must use the following parametrization:
    ``stats.f(dfn, dfd).cdf(x)=fdtr(dfn, dfd, x)``.

    >>> from scipy.stats import f
    >>> dfn, dfd = 1, 2
    >>> x = 1
    >>> fdtr_res = fdtr(dfn, dfd, x)  # this will often be faster than below
    >>> f_dist_res = f(dfn, dfd).cdf(x)
    >>> fdtr_res == f_dist_res  # test that results are equal
    True
    """)

add_newdoc("fdtrc",
    r"""
    fdtrc(dfn, dfd, x, out=None)

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    y : scalar or ndarray
        The complemented F-distribution function with parameters `dfn` and
        `dfd` at `x`.

    See Also
    --------
    fdtr : F distribution cumulative distribution function
    fdtri : F distribution inverse cumulative distribution function
    scipy.stats.f : F distribution

    Notes
    -----
    The regularized incomplete beta function is used, according to the
    formula,

    .. math::
        F(d_n, d_d; x) = I_{d_d/(d_d + xd_n)}(d_d/2, d_n/2).

    Wrapper for a routine from the Boost Math C++ library [1]_. The
    F distribution is also available as `scipy.stats.f`. Calling
    `fdtrc` directly can improve performance compared to the ``sf``
    method of `scipy.stats.f` (see last example below).

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    Calculate the function for ``dfn=1`` and ``dfd=2`` at ``x=1``.

    >>> import numpy as np
    >>> from scipy.special import fdtrc
    >>> fdtrc(1, 2, 1)
    0.42264973081037427

    Calculate the function at several points by providing a NumPy array for
    `x`.

    >>> x = np.array([0.5, 2., 3.])
    >>> fdtrc(1, 2, x)
    array([0.5527864 , 0.29289322, 0.22540333])

    Plot the function for several parameter sets.

    >>> import matplotlib.pyplot as plt
    >>> dfn_parameters = [1, 5, 10, 50]
    >>> dfd_parameters = [1, 1, 2, 3]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(dfn_parameters, dfd_parameters,
    ...                            linestyles))
    >>> x = np.linspace(0, 30, 1000)
    >>> fig, ax = plt.subplots()
    >>> for parameter_set in parameters_list:
    ...     dfn, dfd, style = parameter_set
    ...     fdtrc_vals = fdtrc(dfn, dfd, x)
    ...     ax.plot(x, fdtrc_vals, label=rf"$d_n={dfn},\, d_d={dfd}$",
    ...             ls=style)
    >>> ax.legend()
    >>> ax.set_xlabel("$x$")
    >>> ax.set_title("F distribution survival function")
    >>> plt.show()

    The F distribution is also available as `scipy.stats.f`. Using `fdtrc`
    directly can be much faster than calling the ``sf`` method of
    `scipy.stats.f`, especially for small arrays or individual values.
    To get the same results one must use the following parametrization:
    ``stats.f(dfn, dfd).sf(x)=fdtrc(dfn, dfd, x)``.

    >>> from scipy.stats import f
    >>> dfn, dfd = 1, 2
    >>> x = 1
    >>> fdtrc_res = fdtrc(dfn, dfd, x)  # this will often be faster than below
    >>> f_dist_res = f(dfn, dfd).sf(x)
    >>> f_dist_res == fdtrc_res  # test that results are equal
    True
    """)

add_newdoc("fdtri",
    r"""
    fdtri(dfn, dfd, p, out=None)

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    x : scalar or ndarray
        The quantile corresponding to `p`.

    See Also
    --------
    fdtr : F distribution cumulative distribution function
    fdtrc : F distribution survival function
    scipy.stats.f : F distribution

    Notes
    -----
    Wrapper for a routine from the Boost Math C++ library [1]_. The
    F distribution is also available as `scipy.stats.f`. Calling
    `fdtri` directly can improve performance compared to the ``ppf``
    method of `scipy.stats.f` (see last example below).

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    `fdtri` represents the inverse of the F distribution CDF which is
    available as `fdtr`. Here, we calculate the CDF for ``df1=1``, ``df2=2``
    at ``x=3``. `fdtri` then returns ``3`` given the same values for `df1`,
    `df2` and the computed CDF value.

    >>> import numpy as np
    >>> from scipy.special import fdtri, fdtr
    >>> df1, df2 = 1, 2
    >>> x = 3
    >>> cdf_value =  fdtr(df1, df2, x)
    >>> fdtri(df1, df2, cdf_value)
    3.000000000000006

    Calculate the function at several points by providing a NumPy array for
    `x`.

    >>> x = np.array([0.1, 0.4, 0.7])
    >>> fdtri(1, 2, x)
    array([0.02020202, 0.38095238, 1.92156863])

    Plot the function for several parameter sets.

    >>> import matplotlib.pyplot as plt
    >>> dfn_parameters = [50, 10, 1, 50]
    >>> dfd_parameters = [0.5, 1, 1, 5]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(dfn_parameters, dfd_parameters,
    ...                            linestyles))
    >>> x = np.linspace(0, 1, 1000)
    >>> fig, ax = plt.subplots()
    >>> for parameter_set in parameters_list:
    ...     dfn, dfd, style = parameter_set
    ...     fdtri_vals = fdtri(dfn, dfd, x)
    ...     ax.plot(x, fdtri_vals, label=rf"$d_n={dfn},\, d_d={dfd}$",
    ...             ls=style)
    >>> ax.legend()
    >>> ax.set_xlabel("$x$")
    >>> title = "F distribution inverse cumulative distribution function"
    >>> ax.set_title(title)
    >>> ax.set_ylim(0, 30)
    >>> plt.show()

    The F distribution is also available as `scipy.stats.f`. Using `fdtri`
    directly can be much faster than calling the ``ppf`` method of
    `scipy.stats.f`, especially for small arrays or individual values.
    To get the same results one must use the following parametrization:
    ``stats.f(dfn, dfd).ppf(x)=fdtri(dfn, dfd, x)``.

    >>> from scipy.stats import f
    >>> dfn, dfd = 1, 2
    >>> x = 0.7
    >>> fdtri_res = fdtri(dfn, dfd, x)  # this will often be faster than below
    >>> f_dist_res = f(dfn, dfd).ppf(x)
    >>> f_dist_res == fdtri_res  # test that results are equal
    True
    """)

add_newdoc("fdtridfd",
    """
    fdtridfd(dfn, p, x, out=None)

    Inverse to `fdtr` vs dfd.

    Finds the F density argument dfd such that ``fdtr(dfn, dfd, x) == p``.

    Parameters
    ----------
    dfn : array_like
        First parameter (positive float).
    p : array_like
        Cumulative probability, in [0, 1].
    x : array_like
        Argument (nonnegative float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    dfd : scalar or ndarray
        `dfd` such that ``fdtr(dfn, dfd, x) == p``.

    See Also
    --------
    fdtr : F distribution cumulative distribution function
    fdtrc : F distribution survival function
    fdtri : F distribution quantile function
    scipy.stats.f : F distribution

    Examples
    --------
    Compute the F distribution cumulative distribution function for one
    parameter set.

    >>> from scipy.special import fdtridfd, fdtr
    >>> dfn, dfd, x = 10, 5, 2
    >>> cdf_value = fdtr(dfn, dfd, x)
    >>> cdf_value
    0.7700248806501017

    Verify that `fdtridfd` recovers the original value for `dfd`:

    >>> fdtridfd(dfn, cdf_value, x)
    5.0
    """)

'''
commented out as fdtridfn seems to have bugs and is not in functions.json
see: https://github.com/scipy/scipy/pull/15622#discussion_r811440983

add_newdoc(
    "fdtridfn",
    """
    fdtridfn(p, dfd, x, out=None)

    Inverse to `fdtr` vs dfn

    finds the F density argument dfn such that ``fdtr(dfn, dfd, x) == p``.


    Parameters
    ----------
    p : array_like
        Cumulative probability, in [0, 1].
    dfd : array_like
        Second parameter (positive float).
    x : array_like
        Argument (nonnegative float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    dfn : scalar or ndarray
        `dfn` such that ``fdtr(dfn, dfd, x) == p``.

    See Also
    --------
    fdtr, fdtrc, fdtri, fdtridfd


    """)
'''

add_newdoc("hyp0f1",
    r"""
    hyp0f1(v, z, out=None)

    Confluent hypergeometric limit function 0F1.

    Parameters
    ----------
    v : array_like
        Real-valued parameter
    z : array_like
        Real- or complex-valued argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The confluent hypergeometric limit function

    Notes
    -----
    This function is defined as:

    .. math:: _0F_1(v, z) = \sum_{k=0}^{\infty}\frac{z^k}{(v)_k k!}.

    It's also the limit as :math:`q \to \infty` of :math:`_1F_1(q; v; z/q)`,
    and satisfies the differential equation :math:`f''(z) + vf'(z) =
    f(z)`. See [1]_ for more information.

    References
    ----------
    .. [1] Wolfram MathWorld, "Confluent Hypergeometric Limit Function",
           https://mathworld.wolfram.com/ConfluentHypergeometricLimitFunction.html

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is one when `z` is zero.

    >>> sc.hyp0f1(1, 0)
    1.0

    It is the limit of the confluent hypergeometric function as `q`
    goes to infinity.

    >>> q = np.array([1, 10, 100, 1000])
    >>> v = 1
    >>> z = 1
    >>> sc.hyp1f1(q, v, z / q)
    array([2.71828183, 2.31481985, 2.28303778, 2.27992985])
    >>> sc.hyp0f1(v, z)
    2.2795853023360673

    It is related to Bessel functions.

    >>> n = 1
    >>> x = np.linspace(0, 1, 5)
    >>> sc.jv(n, x)
    array([0.        , 0.12402598, 0.24226846, 0.3492436 , 0.44005059])
    >>> (0.5 * x)**n / sc.factorial(n) * sc.hyp0f1(n + 1, -0.25 * x**2)
    array([0.        , 0.12402598, 0.24226846, 0.3492436 , 0.44005059])

    """)

add_newdoc("hyp1f1",
    r"""
    hyp1f1(a, b, x, out=None)

    Confluent hypergeometric function 1F1.

    The confluent hypergeometric function is defined by the series

    .. math::

       {}_1F_1(a; b; x) = \sum_{k = 0}^\infty \frac{(a)_k}{(b)_k k!} x^k.

    See [DLMF]_ for more details. Here :math:`(\cdot)_k` is the
    Pochhammer symbol; see `poch`.

    Parameters
    ----------
    a, b : array_like
        Real parameters
    x : array_like
        Real or complex argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the confluent hypergeometric function

    See Also
    --------
    hyperu : another confluent hypergeometric function
    hyp0f1 : confluent hypergeometric limit function
    hyp2f1 : Gaussian hypergeometric function

    Notes
    -----
    For real values, this function uses the ``hyp1f1`` routine from the C++ Boost
    library [2]_, for complex values a C translation of the specfun
    Fortran library [3]_.

    References
    ----------
    .. [DLMF] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/13.2#E2
    .. [2] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.
    .. [3] S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996.

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is one when `x` is zero:

    >>> sc.hyp1f1(0.5, 0.5, 0)
    1.0

    It is singular when `b` is a nonpositive integer.

    >>> sc.hyp1f1(0.5, -1, 0)
    inf

    It is a polynomial when `a` is a nonpositive integer.

    >>> a, b, x = -1, 0.5, np.array([1.0, 2.0, 3.0, 4.0])
    >>> sc.hyp1f1(a, b, x)
    array([-1., -3., -5., -7.])
    >>> 1 + (a / b) * x
    array([-1., -3., -5., -7.])

    It reduces to the exponential function when ``a = b``.

    >>> sc.hyp1f1(2, 2, [1, 2, 3, 4])
    array([ 2.71828183,  7.3890561 , 20.08553692, 54.59815003])
    >>> np.exp([1, 2, 3, 4])
    array([ 2.71828183,  7.3890561 , 20.08553692, 54.59815003])

    """)

add_newdoc("kn",
    r"""
    kn(n, x, out=None)

    Modified Bessel function of the second kind of integer order `n`.

    Returns the modified Bessel function of the second kind for integer order
    `n` at real `z`.

    These are also sometimes called functions of the third kind, Basset
    functions, or Macdonald functions.

    Parameters
    ----------
    n : array_like of int
        Order of Bessel functions (floats will truncate with a warning)
    x : array_like of float
        Argument at which to evaluate the Bessel functions
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Value of the Modified Bessel function of the second kind,
        :math:`K_n(x)`.

    See Also
    --------
    kv : Same function, but accepts real order and complex argument
    kvp : Derivative of this function

    Notes
    -----
    Wrapper for AMOS [1]_ routine `zbesk`.  For a discussion of the
    algorithm used, see [2]_ and the references therein.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           https://netlib.org/amos/
    .. [2] Donald E. Amos, "Algorithm 644: A portable package for Bessel
           functions of a complex argument and nonnegative order", ACM
           TOMS Vol. 12 Issue 3, Sept. 1986, p. 265.

    Examples
    --------
    Plot the function of several orders for real input:

    >>> import numpy as np
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

add_newdoc(
    "_landau_pdf",
    """
    _landau_pdf(x, loc, scale)

    Probability density function of the Landau distribution.

    Parameters
    ----------
    x : array_like
        Real-valued argument
    loc : array_like
        Real-valued distribution location
    scale : array_like
        Positive, real-valued distribution scale

    Returns
    -------
    scalar or ndarray
    """)

add_newdoc(
    "_landau_cdf",
    """
    _landau_cdf(x, loc, scale)

    Cumulative distribution function of the Landau distribution.

    Parameters
    ----------
    x : array_like
        Real-valued argument
    loc : array_like
        Real-valued distribution location
    scale : array_like
        Positive, real-valued distribution scale

    Returns
    -------
    scalar or ndarray
    """)

add_newdoc(
    "_landau_sf",
    """
    _landau_sf(x, loc, scale)

    Survival function of the Landau distribution.

    Parameters
    ----------
    x : array_like
        Real-valued argument
    loc : array_like
        Real-valued distribution location
    scale : array_like
        Positive, real-valued distribution scale

    Returns
    -------
    scalar or ndarray
    """)

add_newdoc(
    "_landau_ppf",
    """
    _landau_ppf(p, loc, scale)

    Percent point function of the Landau distribution.

    Parameters
    ----------
    p : array_like
        Real-valued argument between 0 and 1
    loc : array_like
        Real-valued distribution location
    scale : array_like
        Positive, real-valued distribution scale

    Returns
    -------
    scalar or ndarray
    """)

add_newdoc(
    "_landau_isf",
    """
    _landau_isf(p, loc, scale)

    Inverse survival function of the Landau distribution.

    Parameters
    ----------
    p : array_like
        Real-valued argument between 0 and 1
    loc : array_like
        Real-valued distribution location
    scale : array_like
        Positive, real-valued distribution scale

    Returns
    -------
    scalar or ndarray
    """)

add_newdoc("log_gammainc",
    r"""
    log_gammainc(a, x, out=None)

    Logarithm of the regularized lower incomplete gamma function.

    Defined as

    .. math::

        \log P(a, x) = \log \frac{1}{\Gamma(a)} \int_0^x t^{a-1} e^{-t} dt

    for :math:`a > 0` and :math:`x \geq 0`. This function is more accurate
    than computing ``log(gammainc(a, x))`` directly.

    Parameters
    ----------
    a : array_like
        Positive real parameter.
    x : array_like
        Nonneg real argument.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    scalar or ndarray
        Values of the log of the regularized lower incomplete gamma function.

    See Also
    --------
    gammainc : regularized lower incomplete gamma function
    gammaincc : regularized upper incomplete gamma function
    log_gammaincc : log of the regularized upper incomplete gamma function

    Notes
    -----
    This function wraps the ``lgamma_p`` routine from the
    Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    This function is useful when the naive computation ``log(gammainc(a, x))``
    underflows for small values.

    >>> import numpy as np
    >>> from scipy.special import log_gammainc, gammainc

    >>> with np.errstate(divide='ignore'):
    ...    print(np.log(gammainc(500, 30)))
    -inf

    >>> log_gammainc(500, 30)
    -940.6700276993504

    """)

add_newdoc("log_gammaincc",
    r"""
    log_gammaincc(a, x, out=None)

    Logarithm of the regularized upper incomplete gamma function.

    Defined as

    .. math::

        \log Q(a, x) = \log \frac{1}{\Gamma(a)} \int_x^{\infty} t^{a-1} e^{-t} dt

    for :math:`a > 0` and :math:`x \geq 0`. This function is more accurate
    than computing ``log(gammaincc(a, x))`` directly when the survival
    function value is very small.

    Parameters
    ----------
    a : array_like
        Positive real parameter.
    x : array_like
        Nonnegative real argument.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    scalar or ndarray
        Values of the log of the regularized upper incomplete gamma function.

    See Also
    --------
    gammainc : regularized lower incomplete gamma function
    gammaincc : regularized upper incomplete gamma function
    log_gammainc : log of the regularized lower incomplete gamma function

    Notes
    -----
    This function wraps the ``lgamma_q`` routine from the
    Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import log_gammaincc, gammaincc

    For very small survival function values, ``log(gammaincc(a, x))``
    underflows to ``-inf`` while ``log_gammaincc`` retains precision:

    >>> with np.errstate(divide='ignore'):
    ...     print(np.log(gammaincc(10, 1000.0)))
    -inf

    >>> log_gammaincc(10, 1000.0)
    -950.622998370156

    """)

add_newdoc("nbdtr",
    r"""
    nbdtr(k, n, p, out=None)

    Negative binomial cumulative distribution function.

    Returns the sum of the terms 0 through `k` of the negative binomial
    distribution probability mass function,

    .. math::

        F(k) = \sum_{j=0}^k {{n + j - 1}\choose{j}} p^n (1 - p)^j.

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
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The probability of `k` or fewer failures before `n` successes in a
        sequence of events with individual success probability `p`.

    See Also
    --------
    nbdtrc : Negative binomial survival function
    nbdtrik : Negative binomial quantile function
    scipy.stats.nbinom : Negative binomial distribution

    Notes
    -----
    If floating point values are passed for `k` or `n`, they will be truncated
    to integers.

    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{nbdtr}(k, n, p) = I_{p}(n, k + 1).

    Wrapper for the Cephes [1]_ routine `nbdtr`.

    The negative binomial distribution is also available as
    `scipy.stats.nbinom`. Using `nbdtr` directly can improve performance
    compared to the ``cdf`` method of `scipy.stats.nbinom` (see last example).

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    Examples
    --------
    Compute the function for ``k=10`` and ``n=5`` at ``p=0.5``.

    >>> import numpy as np
    >>> from scipy.special import nbdtr
    >>> nbdtr(10, 5, 0.5)
    0.940765380859375

    Compute the function for ``n=10`` and ``p=0.5`` at several points by
    providing a NumPy array or list for `k`.

    >>> nbdtr([5, 10, 15], 10, 0.5)
    array([0.15087891, 0.58809853, 0.88523853])

    Plot the function for four different parameter sets.

    >>> import matplotlib.pyplot as plt
    >>> k = np.arange(130)
    >>> n_parameters = [20, 20, 20, 80]
    >>> p_parameters = [0.2, 0.5, 0.8, 0.5]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(p_parameters, n_parameters,
    ...                            linestyles))
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> for parameter_set in parameters_list:
    ...     p, n, style = parameter_set
    ...     nbdtr_vals = nbdtr(k, n, p)
    ...     ax.plot(k, nbdtr_vals, label=rf"$n={n},\, p={p}$",
    ...             ls=style)
    >>> ax.legend()
    >>> ax.set_xlabel("$k$")
    >>> ax.set_title("Negative binomial cumulative distribution function")
    >>> plt.show()

    The negative binomial distribution is also available as
    `scipy.stats.nbinom`. Using `nbdtr` directly can be much faster than
    calling the ``cdf`` method of `scipy.stats.nbinom`, especially for small
    arrays or individual values. To get the same results one must use the
    following parametrization: ``nbinom(n, p).cdf(k)=nbdtr(k, n, p)``.

    >>> from scipy.stats import nbinom
    >>> k, n, p = 5, 3, 0.5
    >>> nbdtr_res = nbdtr(k, n, p)  # this will often be faster than below
    >>> stats_res = nbinom(n, p).cdf(k)
    >>> stats_res, nbdtr_res  # test that results are equal
    (0.85546875, 0.85546875)

    `nbdtr` can evaluate different parameter sets by providing arrays with
    shapes compatible for broadcasting for `k`, `n` and `p`. Here we compute
    the function for three different `k` at four locations `p`, resulting in
    a 3x4 array.

    >>> k = np.array([[5], [10], [15]])
    >>> p = np.array([0.3, 0.5, 0.7, 0.9])
    >>> k.shape, p.shape
    ((3, 1), (4,))

    >>> nbdtr(k, 5, p)
    array([[0.15026833, 0.62304687, 0.95265101, 0.9998531 ],
           [0.48450894, 0.94076538, 0.99932777, 0.99999999],
           [0.76249222, 0.99409103, 0.99999445, 1.        ]])
    """)

add_newdoc("nbdtrc",
    r"""
    nbdtrc(k, n, p, out=None)

    Negative binomial survival function.

    Returns the sum of the terms `k + 1` to infinity of the negative binomial
    distribution probability mass function,

    .. math::

        S(k) = \sum_{j=k + 1}^\infty {{n + j - 1}\choose{j}} p^n (1 - p)^j.

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
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The probability of `k + 1` or more failures before `n` successes in a
        sequence of events with individual success probability `p`.

    See Also
    --------
    nbdtr : Negative binomial cumulative distribution function
    nbdtrik : Negative binomial percentile function
    scipy.stats.nbinom : Negative binomial distribution

    Notes
    -----
    If floating point values are passed for `k` or `n`, they will be truncated
    to integers.

    The terms are not summed directly; instead the regularized incomplete beta
    function is employed, according to the formula,

    .. math::
        \mathrm{nbdtrc}(k, n, p) = I_{1 - p}(k + 1, n).

    Wrapper for the Cephes [1]_ routine `nbdtrc`.

    The negative binomial distribution is also available as
    `scipy.stats.nbinom`. Using `nbdtrc` directly can improve performance
    compared to the ``sf`` method of `scipy.stats.nbinom` (see last example).

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    Examples
    --------
    Compute the function for ``k=10`` and ``n=5`` at ``p=0.5``.

    >>> import numpy as np
    >>> from scipy.special import nbdtrc
    >>> nbdtrc(10, 5, 0.5)
    0.059234619140624986

    Compute the function for ``n=10`` and ``p=0.5`` at several points by
    providing a NumPy array or list for `k`.

    >>> nbdtrc([5, 10, 15], 10, 0.5)
    array([0.84912109, 0.41190147, 0.11476147])

    Plot the function for four different parameter sets.

    >>> import matplotlib.pyplot as plt
    >>> k = np.arange(130)
    >>> n_parameters = [20, 20, 20, 80]
    >>> p_parameters = [0.2, 0.5, 0.8, 0.5]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(p_parameters, n_parameters,
    ...                            linestyles))
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> for parameter_set in parameters_list:
    ...     p, n, style = parameter_set
    ...     nbdtrc_vals = nbdtrc(k, n, p)
    ...     ax.plot(k, nbdtrc_vals, label=rf"$n={n},\, p={p}$",
    ...             ls=style)
    >>> ax.legend()
    >>> ax.set_xlabel("$k$")
    >>> ax.set_title("Negative binomial distribution survival function")
    >>> plt.show()

    The negative binomial distribution is also available as
    `scipy.stats.nbinom`. Using `nbdtrc` directly can be much faster than
    calling the ``sf`` method of `scipy.stats.nbinom`, especially for small
    arrays or individual values. To get the same results one must use the
    following parametrization: ``nbinom(n, p).sf(k)=nbdtrc(k, n, p)``.

    >>> from scipy.stats import nbinom
    >>> k, n, p = 3, 5, 0.5
    >>> nbdtr_res = nbdtrc(k, n, p)  # this will often be faster than below
    >>> stats_res = nbinom(n, p).sf(k)
    >>> stats_res, nbdtr_res  # test that results are equal
    (0.6367187499999999, 0.6367187499999999)

    `nbdtrc` can evaluate different parameter sets by providing arrays with
    shapes compatible for broadcasting for `k`, `n` and `p`. Here we compute
    the function for three different `k` at four locations `p`, resulting in
    a 3x4 array.

    >>> k = np.array([[5], [10], [15]])
    >>> p = np.array([0.3, 0.5, 0.7, 0.9])
    >>> k.shape, p.shape
    ((3, 1), (4,))

    >>> nbdtrc(k, 5, p)
    array([[8.49731667e-01, 3.76953125e-01, 4.73489874e-02, 1.46902600e-04],
           [5.15491059e-01, 5.92346191e-02, 6.72234070e-04, 9.29610100e-09],
           [2.37507779e-01, 5.90896606e-03, 5.55025308e-06, 3.26346760e-13]])
    """)

add_newdoc(
    "nbdtri",
    r"""
    nbdtri(k, n, y, out=None)

    Returns the inverse with respect to the parameter `p` of
    ``y = nbdtr(k, n, p)``, the negative binomial cumulative distribution
    function.

    Parameters
    ----------
    k : array_like
        The maximum number of allowed failures (nonnegative int).
    n : array_like
        The target number of successes (positive int).
    y : array_like
        The probability of `k` or fewer failures before `n` successes (float).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    p : scalar or ndarray
        Probability of success in a single event (float) such that
        `nbdtr(k, n, p) = y`.

    See Also
    --------
    nbdtr : Cumulative distribution function of the negative binomial.
    nbdtrc : Negative binomial survival function.
    scipy.stats.nbinom : negative binomial distribution.
    nbdtrik : Inverse with respect to `k` of `nbdtr(k, n, p)`.
    nbdtrin : Inverse with respect to `n` of `nbdtr(k, n, p)`.
    scipy.stats.nbinom : Negative binomial distribution

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `nbdtri`.

    The negative binomial distribution is also available as
    `scipy.stats.nbinom`. Using `nbdtri` directly can improve performance
    compared to the ``ppf`` method of `scipy.stats.nbinom`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    Examples
    --------
    `nbdtri` is the inverse of `nbdtr` with respect to `p`.
    Up to floating point errors the following holds:
    ``nbdtri(k, n, nbdtr(k, n, p))=p``.

    >>> import numpy as np
    >>> from scipy.special import nbdtri, nbdtr
    >>> k, n, p = 5, 10, 0.2
    >>> cdf_val = nbdtr(k, n, p)
    >>> nbdtri(k, n, cdf_val)
    0.20000000000000004

    Compute the function for ``k=10`` and ``n=5`` at several points by
    providing a NumPy array or list for `y`.

    >>> y = np.array([0.1, 0.4, 0.8])
    >>> nbdtri(3, 5, y)
    array([0.34462319, 0.51653095, 0.69677416])

    Plot the function for three different parameter sets.

    >>> import matplotlib.pyplot as plt
    >>> n_parameters = [5, 20, 30, 30]
    >>> k_parameters = [20, 20, 60, 80]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(n_parameters, k_parameters, linestyles))
    >>> cdf_vals = np.linspace(0, 1, 1000)
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> for parameter_set in parameters_list:
    ...     n, k, style = parameter_set
    ...     nbdtri_vals = nbdtri(k, n, cdf_vals)
    ...     ax.plot(cdf_vals, nbdtri_vals, label=rf"$k={k},\ n={n}$",
    ...             ls=style)
    >>> ax.legend()
    >>> ax.set_ylabel("$p$")
    >>> ax.set_xlabel("$CDF$")
    >>> title = "nbdtri: inverse of negative binomial CDF with respect to $p$"
    >>> ax.set_title(title)
    >>> plt.show()

    `nbdtri` can evaluate different parameter sets by providing arrays with
    shapes compatible for broadcasting for `k`, `n` and `p`. Here we compute
    the function for three different `k` at four locations `p`, resulting in
    a 3x4 array.

    >>> k = np.array([[5], [10], [15]])
    >>> y = np.array([0.3, 0.5, 0.7, 0.9])
    >>> k.shape, y.shape
    ((3, 1), (4,))

    >>> nbdtri(k, 5, y)
    array([[0.37258157, 0.45169416, 0.53249956, 0.64578407],
           [0.24588501, 0.30451981, 0.36778453, 0.46397088],
           [0.18362101, 0.22966758, 0.28054743, 0.36066188]])
    """)

add_newdoc("nbdtrik",
    r"""
    nbdtrik(y, n, p, out=None)

    Negative binomial percentile function.

    Returns the inverse with respect to the parameter `k` of
    ``y = nbdtr(k, n, p)``, the negative binomial cumulative distribution
    function.

    Parameters
    ----------
    y : array_like
        The probability of `k` or fewer failures before `n` successes (float).
    n : array_like
        The target number of successes (positive int).
    p : array_like
        Probability of success in a single event (float).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    k : scalar or ndarray
        The maximum number of allowed failures such that `nbdtr(k, n, p) = y`.

    See Also
    --------
    nbdtr : Cumulative distribution function of the negative binomial.
    nbdtrc : Survival function of the negative binomial.
    nbdtri : Inverse with respect to `p` of `nbdtr(k, n, p)`.
    nbdtrin : Inverse with respect to `n` of `nbdtr(k, n, p)`.
    scipy.stats.nbinom : Negative binomial distribution

    Notes
    -----
    This function wraps routines from the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    Compute the negative binomial cumulative distribution function for an
    exemplary parameter set.

    >>> import numpy as np
    >>> from scipy.special import nbdtr, nbdtrik
    >>> k, n, p = 5, 2, 0.5
    >>> cdf_value = nbdtr(k, n, p)
    >>> cdf_value
    0.9375

    Verify that `nbdtrik` recovers the original value for `k`.

    >>> nbdtrik(cdf_value, n, p)
    5.0

    Plot the function for different parameter sets.

    >>> import matplotlib.pyplot as plt
    >>> p_parameters = [0.2, 0.5, 0.7, 0.5]
    >>> n_parameters = [30, 30, 30, 80]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(p_parameters, n_parameters, linestyles))
    >>> cdf_vals = np.linspace(0, 1, 1000)
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> for parameter_set in parameters_list:
    ...     p, n, style = parameter_set
    ...     nbdtrik_vals = nbdtrik(cdf_vals, n, p)
    ...     ax.plot(cdf_vals, nbdtrik_vals, label=rf"$n={n},\ p={p}$",
    ...             ls=style)
    >>> ax.legend()
    >>> ax.set_ylabel("$k$")
    >>> ax.set_xlabel("$CDF$")
    >>> ax.set_title("Negative binomial percentile function")
    >>> plt.show()

    The negative binomial distribution is also available as
    `scipy.stats.nbinom`. The percentile function  method ``ppf``
    returns the result of `nbdtrik` rounded up to integers:

    >>> from scipy.stats import nbinom
    >>> q, n, p = 0.6, 5, 0.5
    >>> nbinom.ppf(q, n, p), nbdtrik(q, n, p)
    (5.0, 4.800428460273882)

    """)

add_newdoc("nbdtrin",
    r"""
    nbdtrin(k, y, p, out=None)

    Inverse of `nbdtr` vs `n`.

    Returns the inverse with respect to the parameter `n` of
    ``y = nbdtr(k, n, p)``, the negative binomial cumulative distribution
    function.

    Parameters
    ----------
    k : array_like
        The maximum number of allowed failures (nonnegative int).
    y : array_like
        The probability of `k` or fewer failures before `n` successes (float).
    p : array_like
        Probability of success in a single event (float).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    n : scalar or ndarray
        The number of successes `n` such that `nbdtr(k, n, p) = y`.

    See Also
    --------
    nbdtr : Cumulative distribution function of the negative binomial.
    nbdtri : Inverse with respect to `p` of `nbdtr(k, n, p)`.
    nbdtrik : Inverse with respect to `k` of `nbdtr(k, n, p)`.

    Notes
    -----
    This function wraps routines from the Boost Math C++ library [1]_.

    Formula 26.5.26 of [2]_ or [3]_,

    .. math::
        \sum_{j=k + 1}^\infty {{n + j - 1}
        \choose{j}} p^n (1 - p)^j = I_{1 - p}(k + 1, n),

    is used to reduce calculation of the cumulative distribution function to
    that of a regularized incomplete beta :math:`I`.

    Computation of `n` involves a search for a value that produces the desired
    value of `y`.  The search relies on the monotonicity of `y` with `n`.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [3] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17.E24

    Examples
    --------
    Compute the negative binomial cumulative distribution function for an
    exemplary parameter set.

    >>> from scipy.special import nbdtr, nbdtrin
    >>> k, n, p = 5, 2, 0.5
    >>> cdf_value = nbdtr(k, n, p)
    >>> cdf_value
    0.9375

    Verify that `nbdtrin` recovers the original value for `n` up to floating
    point accuracy.

    >>> nbdtrin(k, cdf_value, p)
    2.0
    """)

add_newdoc("ncfdtr",
    r"""
    ncfdtr(dfn, dfd, nc, f, out=None)

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
        Noncentrality parameter.  Range [0, inf).
    f : array_like
        Quantiles, i.e. the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    cdf : scalar or ndarray
        The calculated CDF.  If all inputs are scalar, the return will be a
        float.  Otherwise it will be an array.

    See Also
    --------
    ncfdtri : Quantile function; inverse of `ncfdtr` with respect to `f`.
    ncfdtridfd : Inverse of `ncfdtr` with respect to `dfd`.
    ncfdtridfn : Inverse of `ncfdtr` with respect to `dfn`.
    ncfdtrinc : Inverse of `ncfdtr` with respect to `nc`.
    scipy.stats.ncf : Non-central F distribution.

    Notes
    -----
    This function calculates the CDF of the non-central f distribution using
    the Boost Math C++ library [1]_.

    The cumulative distribution function is computed using Formula 26.6.20 of
    [2]_:

    .. math::
        F(d_n, d_d, n_c, f) = \sum_{j=0}^\infty e^{-n_c/2}
        \frac{(n_c/2)^j}{j!} I_{x}(\frac{d_n}{2} + j, \frac{d_d}{2}),

    where :math:`I` is the regularized incomplete beta function, and
    :math:`x = f d_n/(f d_n + d_d)`.

    Note that argument order of `ncfdtr` is different from that of the
    similar ``cdf`` method of `scipy.stats.ncf`: `f` is the last
    parameter of `ncfdtr` but the first parameter of ``scipy.stats.ncf.cdf``.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.

    Examples
    --------
    >>> import numpy as np
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

add_newdoc("ncfdtri",
    """
    ncfdtri(dfn, dfd, nc, p, out=None)

    Inverse with respect to `f` of the CDF of the non-central F distribution.

    See `ncfdtr` for more details.

    Parameters
    ----------
    dfn : array_like
        Degrees of freedom of the numerator sum of squares.  Range (0, inf).
    dfd : array_like
        Degrees of freedom of the denominator sum of squares.  Range (0, inf).
    nc : array_like
        Noncentrality parameter.  Range [0, inf).
    p : array_like
        Value of the cumulative distribution function.  Must be in the
        range [0, 1].
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    f : scalar or ndarray
        Quantiles, i.e., the upper limit of integration.

    See Also
    --------
    ncfdtr : CDF of the non-central F distribution.
    ncfdtridfd : Inverse of `ncfdtr` with respect to `dfd`.
    ncfdtridfn : Inverse of `ncfdtr` with respect to `dfn`.
    ncfdtrinc : Inverse of `ncfdtr` with respect to `nc`.
    scipy.stats.ncf : Non-central F distribution.

    Notes
    -----
    This function calculates the Quantile of the non-central f distribution
    using the Boost Math C++ library [1]_.

    Note that argument order of `ncfdtri` is different from that of the
    similar ``ppf`` method of `scipy.stats.ncf`. `p` is the last parameter
    of `ncfdtri` but the first parameter of ``scipy.stats.ncf.ppf``.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import ncfdtr, ncfdtri

    Compute the CDF for several values of `f`:

    >>> f = [0.5, 1, 1.5]
    >>> p = ncfdtr(2, 3, 1.5, f)
    >>> p
    array([ 0.20782291,  0.36107392,  0.47345752])

    Compute the inverse.  We recover the values of `f`, as expected:

    >>> ncfdtri(2, 3, 1.5, p)
    array([ 0.5,  1. ,  1.5])

    """)

add_newdoc("ncfdtridfd",
    """
    ncfdtridfd(dfn, p, nc, f, out=None)

    Calculate degrees of freedom (denominator) for the noncentral F-distribution.

    This is the inverse with respect to `dfd` of `ncfdtr`.
    See `ncfdtr` for more details.

    Parameters
    ----------
    dfn : array_like
        Degrees of freedom of the numerator sum of squares.  Range (0, inf).
    p : array_like
        Value of the cumulative distribution function.  Must be in the
        range [0, 1].
    nc : array_like
        Noncentrality parameter.  Should be in range (0, 1e4).
    f : array_like
        Quantiles, i.e., the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    dfd : scalar or ndarray
        Degrees of freedom of the denominator sum of squares.

    See Also
    --------
    ncfdtr : CDF of the non-central F distribution.
    ncfdtri : Quantile function; inverse of `ncfdtr` with respect to `f`.
    ncfdtridfn : Inverse of `ncfdtr` with respect to `dfn`.
    ncfdtrinc : Inverse of `ncfdtr` with respect to `nc`.

    Notes
    -----
    The value of the cumulative noncentral F distribution is not necessarily
    monotone in either degrees of freedom. There thus may be two values that
    provide a given CDF value. This routine assumes monotonicity and will
    find an arbitrary one of the two values.

    Examples
    --------
    >>> from scipy.special import ncfdtr, ncfdtridfd

    Compute the CDF for several values of `dfd`:

    >>> dfd = [1, 2, 3]
    >>> p = ncfdtr(2, dfd, 0.25, 15)
    >>> p
    array([ 0.8097138 ,  0.93020416,  0.96787852])

    Compute the inverse.  We recover the values of `dfd`, as expected:

    >>> ncfdtridfd(2, p, 0.25, 15)
    array([ 1.,  2.,  3.])

    """)

add_newdoc("ncfdtridfn",
    """
    ncfdtridfn(p, dfd, nc, f, out=None)

    Calculate degrees of freedom (numerator) for the noncentral F-distribution.

    This is the inverse with respect to `dfn` of `ncfdtr`.
    See `ncfdtr` for more details.

    Parameters
    ----------
    p : array_like
        Value of the cumulative distribution function. Must be in the
        range [0, 1].
    dfd : array_like
        Degrees of freedom of the denominator sum of squares. Range (0, inf).
    nc : array_like
        Noncentrality parameter.  Should be in range (0, 1e4).
    f : float
        Quantiles, i.e., the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    dfn : scalar or ndarray
        Degrees of freedom of the numerator sum of squares.

    See Also
    --------
    ncfdtr : CDF of the non-central F distribution.
    ncfdtri : Quantile function; inverse of `ncfdtr` with respect to `f`.
    ncfdtridfd : Inverse of `ncfdtr` with respect to `dfd`.
    ncfdtrinc : Inverse of `ncfdtr` with respect to `nc`.

    Notes
    -----
    The value of the cumulative noncentral F distribution is not necessarily
    monotone in either degrees of freedom. There thus may be two values that
    provide a given CDF value. This routine assumes monotonicity and will
    find an arbitrary one of the two values.

    Examples
    --------
    >>> from scipy.special import ncfdtr, ncfdtridfn

    Compute the CDF for several values of `dfn`:

    >>> dfn = [1, 2, 3]
    >>> p = ncfdtr(dfn, 2, 0.25, 15)
    >>> p
    array([ 0.92562363,  0.93020416,  0.93188394])

    Compute the inverse. We recover the values of `dfn`, as expected:

    >>> ncfdtridfn(p, 2, 0.25, 15)
    array([ 1.,  2.,  3.])

    """)

add_newdoc("ncfdtrinc",
    """
    ncfdtrinc(dfn, dfd, p, f, out=None)

    Calculate non-centrality parameter for non-central F distribution.

    This is the inverse with respect to `nc` of `ncfdtr`.
    See `ncfdtr` for more details.

    Parameters
    ----------
    dfn : array_like
        Degrees of freedom of the numerator sum of squares. Range (0, inf).
    dfd : array_like
        Degrees of freedom of the denominator sum of squares. Range (0, inf).
    p : array_like
        Value of the cumulative distribution function. Must be in the
        range [0, 1].
    f : array_like
        Quantiles, i.e., the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    nc : scalar or ndarray
        Noncentrality parameter.

    See Also
    --------
    ncfdtr : CDF of the non-central F distribution.
    ncfdtri : Quantile function; inverse of `ncfdtr` with respect to `f`.
    ncfdtridfd : Inverse of `ncfdtr` with respect to `dfd`.
    ncfdtridfn : Inverse of `ncfdtr` with respect to `dfn`.

    Notes
    -----
    This function calculates the non-centrality parameter of the
    non-central F distribution given a probability, quantile, and degrees
    of freedom using the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import ncfdtr, ncfdtrinc

    Compute the CDF for several values of `nc`:

    >>> nc = [0.5, 1.5, 2.0]
    >>> p = ncfdtr(2, 3, nc, 15)
    >>> p
    array([ 0.96309246,  0.94327955,  0.93304098])

    Compute the inverse. We recover the values of `nc`, as expected:

    >>> ncfdtrinc(2, 3, p, 15)
    array([ 0.5,  1.5,  2. ])

    """)

add_newdoc("nctdtr",
    """
    nctdtr(df, nc, t, out=None)

    Cumulative distribution function of the non-central `t` distribution.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution. Should be in range (0, inf).
    nc : array_like
        Noncentrality parameter.
    t : array_like
        Quantiles, i.e., the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    cdf : scalar or ndarray
        The calculated CDF. If all inputs are scalar, the return will be a
        float. Otherwise, it will be an array.

    See Also
    --------
    nctdtrit : Inverse CDF (iCDF) of the non-central t distribution.
    nctdtridf : Calculate degrees of freedom, given CDF and iCDF values.
    nctdtrinc : Calculate non-centrality parameter, given CDF iCDF values.

    Notes
    -----
    This function calculates the CDF of the non-central t distribution using
    the Boost Math C++ library [1]_.

    Note that the argument order of `nctdtr` is different from that of the
    similar ``cdf`` method of `scipy.stats.nct`: `t` is the last
    parameter of `nctdtr` but the first parameter of ``scipy.stats.nct.cdf``.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import special
    >>> from scipy import stats
    >>> import matplotlib.pyplot as plt

    Plot the CDF of the non-central t distribution, for nc=0. Compare with the
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

add_newdoc("nctdtridf",
    """
    nctdtridf(p, nc, t, out=None)

    Calculate degrees of freedom for non-central t distribution.

    See `nctdtr` for more details.

    Parameters
    ----------
    p : array_like
        CDF values, in range (0, 1].
    nc : array_like
        Noncentrality parameter. Should be in range (-1e6, 1e6).
    t : array_like
        Quantiles, i.e., the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    df : scalar or ndarray
        The degrees of freedom. If all inputs are scalar, the return will be a
        float. Otherwise, it will be an array.

    See Also
    --------
    nctdtr :  CDF of the non-central `t` distribution.
    nctdtrit : Inverse CDF (iCDF) of the non-central t distribution.
    nctdtrinc : Calculate non-centrality parameter, given CDF iCDF values.

    Notes
    -----
    This function calculates the degrees of freedom of the non-central t
    distribution given a probability, quantile and non-centrality parameter
    using the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import nctdtr, nctdtridf

    Compute the CDF for several values of `df`:

    >>> df = [1, 2, 3]
    >>> p = nctdtr(df, 0.25, 1)
    >>> p
    array([0.67491974, 0.716464  , 0.73349456])

    Compute the inverse. We recover the values of `df`, as expected:

    >>> nctdtridf(p, 0.25, 1)
    array([1., 2., 3.])

    """)

add_newdoc("nctdtrinc",
    """
    nctdtrinc(df, p, t, out=None)

    Calculate non-centrality parameter for non-central t distribution.

    See `nctdtr` for more details.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution. Should be in range (0, inf).
    p : array_like
        CDF values, in range (0, 1].
    t : array_like
        Quantiles, i.e., the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    nc : scalar or ndarray
        Noncentrality parameter

    See Also
    --------
    nctdtr :  CDF of the non-central `t` distribution.
    nctdtrit : Inverse CDF (iCDF) of the non-central t distribution.
    nctdtridf : Calculate degrees of freedom, given CDF and iCDF values.

    Notes
    -----
    This function calculates the non-centrality parameter of the
    non-central t distribution given a probability, quantile and
    degrees of freedom using the Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import nctdtr, nctdtrinc

    Compute the CDF for several values of `nc`:

    >>> nc = [0.5, 1.5, 2.5]
    >>> p = nctdtr(3, nc, 1.5)
    >>> p
    array([0.77569497, 0.45524533, 0.1668691 ])

    Compute the inverse. We recover the values of `nc`, as expected:

    >>> nctdtrinc(3, p, 1.5)
    array([0.5, 1.5, 2.5])

    """)

add_newdoc("nctdtrit",
    """
    nctdtrit(df, nc, p, out=None)

    Inverse cumulative distribution function of the non-central t distribution.

    See `nctdtr` for more details.

    Parameters
    ----------
    df : array_like
        Degrees of freedom of the distribution. Should be in range (0, inf).
    nc : array_like
        Noncentrality parameter.
    p : array_like
        CDF values, in range (0, 1].
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    t : scalar or ndarray
        Quantiles

    See Also
    --------
    nctdtr :  CDF of the non-central `t` distribution.
    nctdtridf : Calculate degrees of freedom, given CDF and iCDF values.
    nctdtrinc : Calculate non-centrality parameter, given CDF iCDF values.

    Notes
    -----
    This function calculates the quantile of the non-central t distribution using
    the Boost Math C++ library [1]_.

    Note that the argument order of `nctdtrit` is different from that of the
    similar ``ppf`` method of `scipy.stats.nct`: `t` is the last
    parameter of `nctdtrit` but the first parameter of ``scipy.stats.nct.ppf``.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> from scipy.special import nctdtr, nctdtrit

    Compute the CDF for several values of `t`:

    >>> t = [0.5, 1, 1.5]
    >>> p = nctdtr(3, 1, t)
    >>> p
    array([0.29811049, 0.46922687, 0.6257559 ])

    Compute the inverse. We recover the values of `t`, as expected:

    >>> nctdtrit(3, 1, p)
    array([0.5, 1. , 1.5])

    """)

add_newdoc("pdtri",
    """
    pdtri(k, y, out=None)

    Inverse of `pdtr` with respect to `m`.

    Returns the Poisson variable `m` such that the sum from 0 to `k` of
    the Poisson density is equal to the given probability `y`:
    calculated by ``gammainccinv(k + 1, y)``. `k` must be a nonnegative
    integer and `y` between 0 and 1.

    Parameters
    ----------
    k : array_like
        Number of occurrences (nonnegative, real).
    y : array_like
        Probability.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the shape parameter `m` such that ``pdtr(k, m) = y``.

    See Also
    --------
    pdtr : Poisson cumulative distribution function
    pdtrc : Poisson survival function
    pdtrik : Inverse of `pdtr` with respect to `k`

    Examples
    --------
    >>> import scipy.special as sc

    Compute the CDF for several values of `m`:

    >>> k = 1
    >>> m = [0.5, 1, 1.5]
    >>> p = sc.pdtr(k, m)
    >>> p
    array([0.90979599, 0.73575888, 0.5578254 ])

    Invert the CDF with respect to the Poisson mean. We recover the values
    of `m`, as expected:

    >>> sc.pdtri(k, p)
    array([0.5, 1. , 1.5])

    Verify the relation with `gammainccinv`:

    >>> sc.gammainccinv(k + 1, p)
    array([0.5, 1. , 1.5])


    """)

add_newdoc("pdtrik",
    """
    pdtrik(p, m, out=None)

    Inverse of `pdtr` with respect to `k`.

    Parameters
    ----------
    p : array_like
        Probability
    m : array_like
        Shape parameter (nonnegative, real)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The number of occurrences `k` such that ``pdtr(k, m) = p``

    See Also
    --------
    pdtr : Poisson cumulative distribution function
    pdtrc : Poisson survival function
    pdtri : inverse of `pdtr` with respect to `m`

    Notes
    -----
    This function relies on the ``gamma_q_inva`` function from the Boost
    Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import scipy.special as sc

    Compute the CDF for several values of `k`:

    >>> k = [1, 2, 3]
    >>> p = sc.pdtr(k, 2)
    >>> p
    array([0.40600585, 0.67667642, 0.85712346])

    Compute the inverse. We recover the values of `k`, as expected:

    >>> sc.pdtrik(p, 2)
    array([1., 2., 3.])

    """)

add_newdoc("powm1", """
    powm1(x, y, out=None)

    Computes ``x**y - 1``.

    This function is useful when `y` is near 0, or when `x` is near 1.

    The function is implemented for real types only (unlike ``numpy.power``,
    which accepts complex inputs).

    Parameters
    ----------
    x : array_like
        The base. Must be a real type (i.e. integer or float, not complex).
    y : array_like
        The exponent. Must be a real type (i.e. integer or float, not complex).

    Returns
    -------
    array_like
        Result of the calculation

    Notes
    -----
    .. versionadded:: 1.10.0

    The underlying code is implemented for single precision and double
    precision floats only.  Unlike `numpy.power`, integer inputs to
    `powm1` are converted to floating point, and complex inputs are
    not accepted.

    Note the following edge cases:

    * ``powm1(x, 0)`` returns 0 for any ``x``, including 0, ``inf``
      and ``nan``.
    * ``powm1(1, y)`` returns 0 for any ``y``, including ``nan``
      and ``inf``.

    This function wraps the ``powm1`` routine from the
    Boost Math C++ library [1]_.

    References
    ----------
    .. [1] The Boost Developers. "Boost C++ Libraries". https://www.boost.org/.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import powm1

    >>> x = np.array([1.2, 10.0, 0.9999999975])
    >>> y = np.array([1e-9, 1e-11, 0.1875])
    >>> powm1(x, y)
    array([ 1.82321557e-10,  2.30258509e-11, -4.68749998e-10])

    It can be verified that the relative errors in those results
    are less than 2.5e-16.

    Compare that to the result of ``x**y - 1``, where the
    relative errors are all larger than 8e-8:

    >>> x**y - 1
    array([ 1.82321491e-10,  2.30258035e-11, -4.68750039e-10])

    """)


add_newdoc("shichi",
    r"""
    shichi(x, out=None)

    Hyperbolic sine and cosine integrals.

    The hyperbolic sine integral is

    .. math::

      \int_0^x \frac{\sinh{t}}{t}dt

    and the hyperbolic cosine integral is

    .. math::

      \gamma + \log(x) + \int_0^x \frac{\cosh{t} - 1}{t} dt

    where :math:`\gamma` is Euler's constant and :math:`\log` is the
    principal branch of the logarithm [1]_ (see also [2]_).

    Parameters
    ----------
    x : array_like
        Real or complex points at which to compute the hyperbolic sine
        and cosine integrals.
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    si : scalar or ndarray
        Hyperbolic sine integral at ``x``
    ci : scalar or ndarray
        Hyperbolic cosine integral at ``x``

    See Also
    --------
    sici : Sine and cosine integrals.
    exp1 : Exponential integral E1.
    expi : Exponential integral Ei.

    Notes
    -----
    For real arguments with ``x < 0``, ``chi`` is the real part of the
    hyperbolic cosine integral. For such points ``chi(x)`` and ``chi(x
    + 0j)`` differ by a factor of ``1j*pi``.

    For real arguments the function is computed by calling Cephes'
    [3]_ *shichi* routine. For complex arguments the algorithm is based
    on Mpmath's [4]_ *shi* and *chi* routines.

    References
    ----------
    .. [1] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
           (See Section 5.2.)
    .. [2] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/6.2.E15 and https://dlmf.nist.gov/6.2.E16
    .. [3] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/
    .. [4] Fredrik Johansson and others.
           "mpmath: a Python library for arbitrary-precision floating-point
           arithmetic" (Version 0.19) https://mpmath.org/

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.special import shichi, sici

    `shichi` accepts real or complex input:

    >>> shichi(0.5)
    (0.5069967498196671, -0.05277684495649357)
    >>> shichi(0.5 + 2.5j)
    ((0.11772029666668238+1.831091777729851j),
     (0.29912435887648825+1.7395351121166562j))

    The hyperbolic sine and cosine integrals Shi(z) and Chi(z) are
    related to the sine and cosine integrals Si(z) and Ci(z) by

    * Shi(z) = -i*Si(i*z)
    * Chi(z) = Ci(-i*z) + i*pi/2

    >>> z = 0.25 + 5j
    >>> shi, chi = shichi(z)
    >>> shi, -1j*sici(1j*z)[0]            # Should be the same.
    ((-0.04834719325101729+1.5469354086921228j),
     (-0.04834719325101729+1.5469354086921228j))
    >>> chi, sici(-1j*z)[1] + 1j*np.pi/2  # Should be the same.
    ((-0.19568708973868087+1.556276312103824j),
     (-0.19568708973868087+1.556276312103824j))

    Plot the functions evaluated on the real axis:

    >>> xp = np.geomspace(1e-8, 4.0, 250)
    >>> x = np.concatenate((-xp[::-1], xp))
    >>> shi, chi = shichi(x)

    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, shi, label='Shi(x)')
    >>> ax.plot(x, chi, '--', label='Chi(x)')
    >>> ax.set_xlabel('x')
    >>> ax.set_title('Hyperbolic Sine and Cosine Integrals')
    >>> ax.legend(shadow=True, framealpha=1, loc='lower right')
    >>> ax.grid(True)
    >>> plt.show()

    """)

add_newdoc("sici",
    r"""
    sici(x, out=None)

    Sine and cosine integrals.

    The sine integral is

    .. math::

      \int_0^x \frac{\sin{t}}{t}dt

    and the cosine integral is

    .. math::

      \gamma + \log(x) + \int_0^x \frac{\cos{t} - 1}{t}dt

    where :math:`\gamma` is Euler's constant and :math:`\log` is the
    principal branch of the logarithm [1]_ (see also [2]_).

    Parameters
    ----------
    x : array_like
        Real or complex points at which to compute the sine and cosine
        integrals.
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    si : scalar or ndarray
        Sine integral at ``x``
    ci : scalar or ndarray
        Cosine integral at ``x``

    See Also
    --------
    shichi : Hyperbolic sine and cosine integrals.
    exp1 : Exponential integral E1.
    expi : Exponential integral Ei.

    Notes
    -----
    For real arguments with ``x < 0``, ``ci`` is the real part of the
    cosine integral. For such points ``ci(x)`` and ``ci(x + 0j)``
    differ by a factor of ``1j*pi``.

    For real arguments the function is computed by calling Cephes'
    [3]_ *sici* routine. For complex arguments the algorithm is based
    on Mpmath's [4]_ *si* and *ci* routines.

    References
    ----------
    .. [1] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
           (See Section 5.2.)
    .. [2] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/6.2.E9, https://dlmf.nist.gov/6.2.E12,
           and https://dlmf.nist.gov/6.2.E13
    .. [3] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/
    .. [4] Fredrik Johansson and others.
           "mpmath: a Python library for arbitrary-precision floating-point
           arithmetic" (Version 0.19) https://mpmath.org/

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.special import sici, exp1

    `sici` accepts real or complex input:

    >>> sici(2.5)
    (1.7785201734438267, 0.2858711963653835)
    >>> sici(2.5 + 3j)
    ((4.505735874563953+0.06863305018999577j),
    (0.0793644206906966-2.935510262937543j))

    For z in the right half plane, the sine and cosine integrals are
    related to the exponential integral E1 (implemented in SciPy as
    `scipy.special.exp1`) by

    * Si(z) = (E1(i*z) - E1(-i*z))/2i + pi/2
    * Ci(z) = -(E1(i*z) + E1(-i*z))/2

    See [1]_ (equations 5.2.21 and 5.2.23).

    We can verify these relations:

    >>> z = 2 - 3j
    >>> sici(z)
    ((4.54751388956229-1.3991965806460565j),
    (1.408292501520851+2.9836177420296055j))

    >>> (exp1(1j*z) - exp1(-1j*z))/2j + np.pi/2  # Same as sine integral
    (4.54751388956229-1.3991965806460565j)

    >>> -(exp1(1j*z) + exp1(-1j*z))/2            # Same as cosine integral
    (1.408292501520851+2.9836177420296055j)

    Plot the functions evaluated on the real axis; the dotted horizontal
    lines are at pi/2 and -pi/2:

    >>> x = np.linspace(-16, 16, 150)
    >>> si, ci = sici(x)

    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, si, label='Si(x)')
    >>> ax.plot(x, ci, '--', label='Ci(x)')
    >>> ax.legend(shadow=True, framealpha=1, loc='upper left')
    >>> ax.set_xlabel('x')
    >>> ax.set_title('Sine and Cosine Integrals')
    >>> ax.axhline(np.pi/2, linestyle=':', alpha=0.5, color='k')
    >>> ax.axhline(-np.pi/2, linestyle=':', alpha=0.5, color='k')
    >>> ax.grid(True)
    >>> plt.show()

    """)

add_newdoc("smirnov",
    r"""
    smirnov(n, d, out=None)

    Kolmogorov-Smirnov complementary cumulative distribution function.

    Returns the exact Kolmogorov-Smirnov complementary cumulative
    distribution function,(aka the Survival Function) of Dn+ (or Dn-)
    for a one-sided test of equality between an empirical and a
    theoretical distribution. It is equal to the probability that the
    maximum difference between a theoretical distribution and an empirical
    one based on `n` samples is greater than d.

    Parameters
    ----------
    n : int
      Number of samples
    d : float array_like
      Deviation between the Empirical CDF (ECDF) and the target CDF.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The value(s) of smirnov(n, d), Prob(Dn+ >= d) (Also Prob(Dn- >= d))

    See Also
    --------
    smirnovi : The Inverse Survival Function for the distribution
    scipy.stats.ksone : Provides the functionality as a continuous distribution
    kolmogorov, kolmogi : Functions for the two-sided distribution

    Notes
    -----
    `smirnov` is used by `stats.kstest` in the application of the
    Kolmogorov-Smirnov Goodness of Fit test. For historical reasons this
    function is exposed in `scpy.special`, but the recommended way to achieve
    the most accurate CDF/SF/PDF/PPF/ISF computations is to use the
    `stats.ksone` distribution.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import smirnov
    >>> from scipy.stats import norm

    Show the probability of a gap at least as big as 0, 0.5 and 1.0 for a
    sample of size 5.

    >>> smirnov(5, [0, 0.5, 1.0])
    array([ 1.   ,  0.056,  0.   ])

    Compare a sample of size 5 against N(0, 1), the standard normal
    distribution with mean 0 and standard deviation 1.

    `x` is the sample.

    >>> x = np.array([-1.392, -0.135, 0.114, 0.190, 1.82])

    >>> target = norm(0, 1)
    >>> cdfs = target.cdf(x)
    >>> cdfs
    array([0.0819612 , 0.44630594, 0.5453811 , 0.57534543, 0.9656205 ])

    Construct the empirical CDF and the K-S statistics (Dn+, Dn-, Dn).

    >>> n = len(x)
    >>> ecdfs = np.arange(n+1, dtype=float)/n
    >>> cols = np.column_stack([x, ecdfs[1:], cdfs, cdfs - ecdfs[:n],
    ...                        ecdfs[1:] - cdfs])
    >>> with np.printoptions(precision=3):
    ...    print(cols)
    [[-1.392  0.2    0.082  0.082  0.118]
     [-0.135  0.4    0.446  0.246 -0.046]
     [ 0.114  0.6    0.545  0.145  0.055]
     [ 0.19   0.8    0.575 -0.025  0.225]
     [ 1.82   1.     0.966  0.166  0.034]]
    >>> gaps = cols[:, -2:]
    >>> Dnpm = np.max(gaps, axis=0)
    >>> print(f'Dn-={Dnpm[0]:f}, Dn+={Dnpm[1]:f}')
    Dn-=0.246306, Dn+=0.224655
    >>> probs = smirnov(n, Dnpm)
    >>> print(f'For a sample of size {n} drawn from N(0, 1):',
    ...       f' Smirnov n={n}: Prob(Dn- >= {Dnpm[0]:f}) = {probs[0]:.4f}',
    ...       f' Smirnov n={n}: Prob(Dn+ >= {Dnpm[1]:f}) = {probs[1]:.4f}',
    ...       sep='\n')
    For a sample of size 5 drawn from N(0, 1):
     Smirnov n=5: Prob(Dn- >= 0.246306) = 0.4711
     Smirnov n=5: Prob(Dn+ >= 0.224655) = 0.5245

    Plot the empirical CDF and the standard normal CDF.

    >>> import matplotlib.pyplot as plt
    >>> plt.step(np.concatenate(([-2.5], x, [2.5])),
    ...          np.concatenate((ecdfs, [1])),
    ...          where='post', label='Empirical CDF')
    >>> xx = np.linspace(-2.5, 2.5, 100)
    >>> plt.plot(xx, target.cdf(xx), '--', label='CDF for N(0, 1)')

    Add vertical lines marking Dn+ and Dn-.

    >>> iminus, iplus = np.argmax(gaps, axis=0)
    >>> plt.vlines([x[iminus]], ecdfs[iminus], cdfs[iminus], color='r',
    ...            alpha=0.5, lw=4)
    >>> plt.vlines([x[iplus]], cdfs[iplus], ecdfs[iplus+1], color='m',
    ...            alpha=0.5, lw=4)

    >>> plt.grid(True)
    >>> plt.legend(framealpha=1, shadow=True)
    >>> plt.show()
    """)

add_newdoc("smirnovi",
    """
    smirnovi(n, p, out=None)

    Inverse to `smirnov`.

    Returns `d` such that ``smirnov(n, d) == p``, the critical value
    corresponding to `p`.

    Parameters
    ----------
    n : int
      Number of samples
    p : float array_like
        Probability
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The value(s) of smirnovi(n, p), the critical values.

    See Also
    --------
    smirnov : The Survival Function (SF) for the distribution
    scipy.stats.ksone : Provides the functionality as a continuous distribution
    kolmogorov, kolmogi : Functions for the two-sided distribution
    scipy.stats.kstwobign : Two-sided Kolmogorov-Smirnov distribution, large n

    Notes
    -----
    `smirnov` is used by `stats.kstest` in the application of the
    Kolmogorov-Smirnov Goodness of Fit test. For historical reasons this
    function is exposed in `scpy.special`, but the recommended way to achieve
    the most accurate CDF/SF/PDF/PPF/ISF computations is to use the
    `stats.ksone` distribution.

    Examples
    --------
    >>> from scipy.special import smirnovi, smirnov

    >>> n = 24
    >>> deviations = [0.1, 0.2, 0.3]

    Use `smirnov` to compute the complementary CDF of the Smirnov
    distribution for the given number of samples and deviations.

    >>> p = smirnov(n, deviations)
    >>> p
    array([0.58105083, 0.12826832, 0.01032231])

    The inverse function ``smirnovi(n, p)`` returns ``deviations``.

    >>> smirnovi(n, p)
    array([0.1, 0.2, 0.3])

    """)

add_newdoc("_smirnovc",
    """
    _smirnovc(n, d)
     Internal function, do not use.
    """)

add_newdoc("_smirnovci",
    """
     Internal function, do not use.
    """)

add_newdoc("_smirnovp",
    """
    _smirnovp(n, p)
     Internal function, do not use.
    """)

add_newdoc(
    "stdtr",
    r"""
    stdtr(df, t, out=None)

    Student t distribution cumulative distribution function.

    Returns the integral:

    .. math::
        \frac{\Gamma((df+1)/2)}{\sqrt{\pi df} \Gamma(df/2)}
        \int_{-\infty}^t (1+x^2/df)^{-(df+1)/2}\, dx

    Parameters
    ----------
    df : array_like
        Degrees of freedom
    t : array_like
        Upper bound of the integral
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Value of the Student t CDF at t

    See Also
    --------
    stdtridf : inverse of stdtr with respect to `df`
    stdtrit : inverse of stdtr with respect to `t`
    scipy.stats.t : student t distribution

    Notes
    -----
    The student t distribution is also available as `scipy.stats.t`.
    Calling `stdtr` directly can improve performance compared to the
    ``cdf`` method of `scipy.stats.t` (see last example below).

    The function is computed using the Boost Math library [1]_, which
    relies on the incomplete beta function.

    References
    ----------
    .. [1] Boost C++ Libraries, https://www.boost.org/

    Examples
    --------
    Calculate the function for ``df=3`` at ``t=1``.

    >>> import numpy as np
    >>> from scipy.special import stdtr
    >>> import matplotlib.pyplot as plt
    >>> stdtr(3, 1)
    0.8044988905221148

    Plot the function for three different degrees of freedom.

    >>> x = np.linspace(-10, 10, 1000)
    >>> fig, ax = plt.subplots()
    >>> parameters = [(1, "solid"), (3, "dashed"), (10, "dotted")]
    >>> for (df, linestyle) in parameters:
    ...     ax.plot(x, stdtr(df, x), ls=linestyle, label=f"$df={df}$")
    >>> ax.legend()
    >>> ax.set_title("Student t distribution cumulative distribution function")
    >>> plt.show()

    The function can be computed for several degrees of freedom at the same
    time by providing a NumPy array or list for `df`:

    >>> stdtr([1, 2, 3], 1)
    array([0.75      , 0.78867513, 0.80449889])

    It is possible to calculate the function at several points for several
    different degrees of freedom simultaneously by providing arrays for `df`
    and `t` with shapes compatible for broadcasting. Compute `stdtr` at
    4 points for 3 degrees of freedom resulting in an array of shape 3x4.

    >>> dfs = np.array([[1], [2], [3]])
    >>> t = np.array([2, 4, 6, 8])
    >>> dfs.shape, t.shape
    ((3, 1), (4,))

    >>> stdtr(dfs, t)
    array([[0.85241638, 0.92202087, 0.94743154, 0.96041658],
           [0.90824829, 0.97140452, 0.98666426, 0.99236596],
           [0.93033702, 0.98599577, 0.99536364, 0.99796171]])

    The t distribution is also available as `scipy.stats.t`. Calling `stdtr`
    directly can be much faster than calling the ``cdf`` method of
    `scipy.stats.t`. To get the same results, one must use the following
    parametrization: ``scipy.stats.t(df).cdf(x) = stdtr(df, x)``.

    >>> from scipy.stats import t
    >>> df, x = 3, 1
    >>> stdtr_result = stdtr(df, x)  # this can be faster than below
    >>> stats_result = t(df).cdf(x)
    >>> stats_result == stdtr_result  # test that results are equal
    True
    """)

add_newdoc("stdtridf",
    """
    stdtridf(p, t, out=None)

    Inverse of `stdtr` vs df.

    Returns the argument df such that stdtr(df, t) is equal to `p`.

    Parameters
    ----------
    p : array_like
        Probability
    t : array_like
        Upper bound of the integral
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    df : scalar or ndarray
        Value of `df` such that ``stdtr(df, t) == p``

    See Also
    --------
    stdtr : Student t CDF
    stdtrit : inverse of stdtr with respect to `t`
    scipy.stats.t : Student t distribution

    Examples
    --------
    Compute the student t cumulative distribution function for one
    parameter set.

    >>> from scipy.special import stdtr, stdtridf
    >>> df, x = 5, 2
    >>> cdf_value = stdtr(df, x)
    >>> cdf_value
    0.9490302605850709

    Verify that `stdtridf` recovers the original value for `df` given
    the CDF value and `x`.

    >>> stdtridf(cdf_value, x)
    5.0
    """)

add_newdoc("stdtrit",
    """
    stdtrit(df, p, out=None)

    The `p`-th quantile of the student t distribution.

    This function is the inverse of the student t distribution cumulative
    distribution function (CDF), returning `t` such that `stdtr(df, t) = p`.

    Returns the argument `t` such that stdtr(df, t) is equal to `p`.

    Parameters
    ----------
    df : array_like
        Degrees of freedom
    p : array_like
        Probability
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    t : scalar or ndarray
        Value of `t` such that ``stdtr(df, t) == p``

    See Also
    --------
    stdtr : Student t CDF
    stdtridf : inverse of stdtr with respect to `df`
    scipy.stats.t : Student t distribution

    Notes
    -----
    The student t distribution is also available as `scipy.stats.t`. Calling
    `stdtrit` directly can improve performance compared to the ``ppf``
    method of `scipy.stats.t` (see last example below).

    The function is computed using the Boost Math library [1]_, which
    relies on the incomplete beta function.

    References
    ----------
    .. [1] Boost C++ Libraries, https://www.boost.org/

    Examples
    --------
    `stdtrit` represents the inverse of the student t distribution CDF which
    is available as `stdtr`. Here, we calculate the CDF for ``df`` at
    ``x=1``. `stdtrit` then returns ``1`` up to floating point errors
    given the same value for `df` and the computed CDF value.

    >>> import numpy as np
    >>> from scipy.special import stdtr, stdtrit
    >>> import matplotlib.pyplot as plt
    >>> df = 3
    >>> x = 1
    >>> cdf_value = stdtr(df, x)
    >>> stdtrit(df, cdf_value)
    0.9999999994418539

    Plot the function for three different degrees of freedom.

    >>> x = np.linspace(0, 1, 1000)
    >>> parameters = [(1, "solid"), (2, "dashed"), (5, "dotted")]
    >>> fig, ax = plt.subplots()
    >>> for (df, linestyle) in parameters:
    ...     ax.plot(x, stdtrit(df, x), ls=linestyle, label=f"$df={df}$")
    >>> ax.legend()
    >>> ax.set_ylim(-10, 10)
    >>> ax.set_title("Student t distribution quantile function")
    >>> plt.show()

    The function can be computed for several degrees of freedom at the same
    time by providing a NumPy array or list for `df`:

    >>> stdtrit([1, 2, 3], 0.7)
    array([0.72654253, 0.6172134 , 0.58438973])

    It is possible to calculate the function at several points for several
    different degrees of freedom simultaneously by providing arrays for `df`
    and `p` with shapes compatible for broadcasting. Compute `stdtrit` at
    4 points for 3 degrees of freedom resulting in an array of shape 3x4.

    >>> dfs = np.array([[1], [2], [3]])
    >>> p = np.array([0.2, 0.4, 0.7, 0.8])
    >>> dfs.shape, p.shape
    ((3, 1), (4,))

    >>> stdtrit(dfs, p)
    array([[-1.37638192, -0.3249197 ,  0.72654253,  1.37638192],
           [-1.06066017, -0.28867513,  0.6172134 ,  1.06066017],
           [-0.97847231, -0.27667066,  0.58438973,  0.97847231]])

    The t distribution is also available as `scipy.stats.t`. Calling `stdtrit`
    directly can be much faster than calling the ``ppf`` method of
    `scipy.stats.t`. To get the same results, one must use the following
    parametrization: ``scipy.stats.t(df).ppf(x) = stdtrit(df, x)``.

    >>> from scipy.stats import t
    >>> df, x = 3, 0.5
    >>> stdtrit_result = stdtrit(df, x)  # this can be faster than below
    >>> stats_result = t(df).ppf(x)
    >>> stats_result == stdtrit_result  # test that results are equal
    True
    """)

add_newdoc("yn",
    r"""
    yn(n, x, out=None)

    Bessel function of the second kind of integer order and real argument.

    Parameters
    ----------
    n : array_like
        Order (integer).
    x : array_like
        Argument (float).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    Y : scalar or ndarray
        Value of the Bessel function, :math:`Y_n(x)`.

    See Also
    --------
    yv : For real order and real or complex argument.
    y0: faster implementation of this function for order 0
    y1: faster implementation of this function for order 1

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `yn`.

    The function is evaluated by forward recurrence on `n`, starting with
    values computed by the Cephes routines `y0` and `y1`. If ``n = 0`` or 1,
    the routine for `y0` or `y1` is called directly.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    Examples
    --------
    Evaluate the function of order 0 at one point.

    >>> from scipy.special import yn
    >>> yn(0, 1.)
    0.08825696421567697

    Evaluate the function at one point for different orders.

    >>> yn(0, 1.), yn(1, 1.), yn(2, 1.)
    (0.08825696421567697, -0.7812128213002888, -1.6506826068162546)

    The evaluation for different orders can be carried out in one call by
    providing a list or NumPy array as argument for the `v` parameter:

    >>> yn([0, 1, 2], 1.)
    array([ 0.08825696, -0.78121282, -1.65068261])

    Evaluate the function at several points for order 0 by providing an
    array for `z`.

    >>> import numpy as np
    >>> points = np.array([0.5, 3., 8.])
    >>> yn(0, points)
    array([-0.44451873,  0.37685001,  0.22352149])

    If `z` is an array, the order parameter `v` must be broadcastable to
    the correct shape if different orders shall be computed in one call.
    To calculate the orders 0 and 1 for a 1D array:

    >>> orders = np.array([[0], [1]])
    >>> orders.shape
    (2, 1)

    >>> yn(orders, points)
    array([[-0.44451873,  0.37685001,  0.22352149],
           [-1.47147239,  0.32467442, -0.15806046]])

    Plot the functions of order 0 to 3 from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, yn(i, x), label=f'$Y_{i!r}$')
    >>> ax.set_ylim(-3, 1)
    >>> ax.legend()
    >>> plt.show()
    """)

add_newdoc("_stirling2_inexact",
    r"""
    Internal function, do not use.
    """)

add_newdoc(
    "_beta_pdf",
    r"""
    _beta_pdf(x, a, b)

    Probability density function of beta distribution.

    Parameters
    ----------
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration
    a, b : array_like
           Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_beta_ppf",
    r"""
    _beta_ppf(x, a, b)

    Percent point function of beta distribution.

    Parameters
    ----------
    x : array_like
        Real-valued such that :math:`0 \leq x \leq 1`,
        the upper limit of integration
    a, b : array_like
           Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_invgauss_ppf",
    """
    _invgauss_ppf(x, mu)

    Percent point function of inverse gaussian distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    mu : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_invgauss_isf",
    """
    _invgauss_isf(x, mu, s)

    Inverse survival function of inverse gaussian distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    mu : array_like
        Positive, real-valued parameters
    s : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_cauchy_ppf",
    """
    _cauchy_ppf(p, loc, scale)

    Percent point function (i.e. quantile) of the Cauchy distribution.

    Parameters
    ----------
    p : array_like
        Probabilities
    loc : array_like
        Location parameter of the distribution.
    scale : array_like
        Scale parameter of the distribution.

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_cauchy_isf",
    """
    _cauchy_isf(p, loc, scale)

    Inverse survival function of the Cauchy distribution.

    Parameters
    ----------
    p : array_like
        Probabilities
    loc : array_like
        Location parameter of the distribution.
    scale : array_like
        Scale parameter of the distribution.

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_pdf",
    """
    _ncx2_pdf(x, k, l)

    Probability density function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_cdf",
    """
    _ncx2_cdf(x, k, l)

    Cumulative density function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_ppf",
    """
    _ncx2_ppf(x, k, l)

    Percent point function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_sf",
    """
    _ncx2_sf(x, k, l)

    Survival function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncx2_isf",
    """
    _ncx2_isf(x, k, l)

    Inverse survival function of Non-central chi-squared distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    k, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_pdf",
    """
    _ncf_pdf(x, v1, v2, l)

    Probability density function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_cdf",
    """
    _ncf_cdf(x, v1, v2, l)

    Cumulative density function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_ppf",
    """
    _ncf_ppf(x, v1, v2, l)

    Percent point function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_sf",
    """
    _ncf_sf(x, v1, v2, l)

    Survival function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_isf",
    """
    _ncf_isf(x, v1, v2, l)

    Inverse survival function of noncentral F-distribution.

    Parameters
    ----------
    x : array_like
        Positive real-valued
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_mean",
    """
    _ncf_mean(v1, v2, l)

    Mean of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_variance",
    """
    _ncf_variance(v1, v2, l)

    Variance of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_skewness",
    """
    _ncf_skewness(v1, v2, l)

    Skewness of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_ncf_kurtosis_excess",
    """
    _ncf_kurtosis_excess(v1, v2, l)

    Kurtosis excess of noncentral F-distribution.

    Parameters
    ----------
    v1, v2, l : array_like
        Positive, real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_cdf",
    """
    _nct_cdf(x, v, l)

    Cumulative density function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_pdf",
    """
    _nct_pdf(x, v, l)

    Probability density function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)


add_newdoc(
    "_nct_ppf",
    """
    _nct_ppf(x, v, l)

    Percent point function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_sf",
    """
    _nct_sf(x, v, l)

    Survival function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_isf",
    """
    _nct_isf(x, v, l)

    Inverse survival function of noncentral t-distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_mean",
    """
    _nct_mean(v, l)

    Mean of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_variance",
    """
    _nct_variance(v, l)

    Variance of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_skewness",
    """
    _nct_skewness(v, l)

    Skewness of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nct_kurtosis_excess",
    """
    _nct_kurtosis_excess(v, l)

    Kurtosis excess of noncentral t-distribution.

    Parameters
    ----------
    v : array_like
        Positive, real-valued parameters
    l : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_skewnorm_cdf",
    """
    _skewnorm_cdf(x, l, sc, sh)

    Cumulative density function of skewnorm distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    l : array_like
        Real-valued parameters
    sc : array_like
        Positive, Real-valued parameters
    sh : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_skewnorm_ppf",
    """
    _skewnorm_ppf(x, l, sc, sh)

    Percent point function of skewnorm distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    l : array_like
        Real-valued parameters
    sc : array_like
        Positive, Real-valued parameters
    sh : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_skewnorm_isf",
    """
    _skewnorm_isf(x, l, sc, sh)

    Inverse survival function of skewnorm distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    l : array_like
        Real-valued parameters
    sc : array_like
        Positive, Real-valued parameters
    sh : array_like
        Real-valued parameters

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_pmf",
    """
    _binom_pmf(x, n, p)

    Probability mass function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_cdf",
    """
    _binom_cdf(x, n, p)

    Cumulative density function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_ppf",
    """
    _binom_ppf(x, n, p)

    Percent point function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_sf",
    """
    _binom_sf(x, n, p)

    Survival function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_binom_isf",
    """
    _binom_isf(x, n, p)

    Inverse survival function of binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    n : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_pmf",
    """
    _nbinom_pmf(x, r, p)

    Probability mass function of negative binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_cdf",
    """
    _nbinom_cdf(x, r, p)

    Cumulative density function of negative binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_ppf",
    """
    _nbinom_ppf(x, r, p)

    Percent point function of negative binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_sf",
    """
    _nbinom_sf(x, r, p)

    Survival function of negative binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_isf",
    """
    _nbinom_isf(x, r, p)

    Inverse survival function of negative binomial distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_mean",
    """
    _nbinom_mean(r, p)

    Mean of negative binomial distribution.

    Parameters
    ----------
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_variance",
    """
    _nbinom_variance(r, p)

    Variance of negative binomial distribution.

    Parameters
    ----------
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_skewness",
    """
    _nbinom_skewness(r, p)

    Skewness of negative binomial distribution.

    Parameters
    ----------
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_nbinom_kurtosis_excess",
    """
    _nbinom_kurtosis_excess(r, p)

    Kurtosis excess of negative binomial distribution.

    Parameters
    ----------
    r : array_like
        Positive, integer-valued parameter
    p : array_like
        Positive, real-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_hypergeom_pmf",
    """
    _hypergeom_pmf(x, r, N, M)

    Probability mass function of hypergeometric distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r, N, M : array_like
        Positive, integer-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_hypergeom_cdf",
    """
    _hypergeom_cdf(x, r, N, M)

    Cumulative density function of hypergeometric distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r, N, M : array_like
        Positive, integer-valued parameter

    Returns
    -------
    scalar or ndarray
    """)

add_newdoc(
    "_hypergeom_sf",
    """
    _hypergeom_sf(x, r, N, M)

    Survival function of hypergeometric distribution.

    Parameters
    ----------
    x : array_like
        Real-valued
    r, N, M : array_like
        Positive, integer-valued parameter

    Returns
    -------
    scalar or ndarray
    """)

add_newdoc(
    "_hypergeom_mean",
    """
    _hypergeom_mean(r, N, M)

    Mean of hypergeometric distribution.

    Parameters
    ----------
    r, N, M : array_like
        Positive, integer-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_hypergeom_variance",
    """
    _hypergeom_variance(r, N, M)

    Variance of hypergeometric distribution.

    Parameters
    ----------
    r, N, M : array_like
        Positive, integer-valued parameter

    Returns
    -------
    scalar or ndarray

    """)

add_newdoc(
    "_hypergeom_skewness",
    """
    _hypergeom_skewness(r, N, M)

    Skewness of hypergeometric distribution.

    Parameters
    ----------
    r, N, M : array_like
        Positive, integer-valued parameter

    Returns
    -------
    scalar or ndarray

    """)
