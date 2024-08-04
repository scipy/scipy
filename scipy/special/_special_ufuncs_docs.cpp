const char *_cospi_doc = R"(
    Internal function, do not use.
    )";

const char *_sinpi_doc = R"(
    Internal function, do not use.
    )";

const char *_zeta_doc = R"(
    _zeta(x, q)

    Internal function, Hurwitz zeta.

    )";

const char *airy_doc = R"(
    airy(z, out=None)

    Airy functions and their derivatives.

    Parameters
    ----------
    z : array_like
        Real or complex argument.
    out : tuple of ndarray, optional
        Optional output arrays for the function values

    Returns
    -------
    Ai, Aip, Bi, Bip : 4-tuple of scalar or ndarray
        Airy functions Ai and Bi, and their derivatives Aip and Bip.

    See Also
    --------
    airye : exponentially scaled Airy functions.

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

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/
    .. [2] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    Compute the Airy functions on the interval [-15, 5].

    >>> import numpy as np
    >>> from scipy import special
    >>> x = np.linspace(-15, 5, 201)
    >>> ai, aip, bi, bip = special.airy(x)

    Plot Ai(x) and Bi(x).

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, ai, 'r', label='Ai(x)')
    >>> plt.plot(x, bi, 'b--', label='Bi(x)')
    >>> plt.ylim(-0.5, 1.0)
    >>> plt.grid()
    >>> plt.legend(loc='upper left')
    >>> plt.show()

    )";

const char *airye_doc = R"(
    airye(z, out=None)

    Exponentially scaled Airy functions and their derivatives.

    Scaling::

        eAi  = Ai  * exp(2.0/3.0*z*sqrt(z))
        eAip = Aip * exp(2.0/3.0*z*sqrt(z))
        eBi  = Bi  * exp(-abs(2.0/3.0*(z*sqrt(z)).real))
        eBip = Bip * exp(-abs(2.0/3.0*(z*sqrt(z)).real))

    Parameters
    ----------
    z : array_like
        Real or complex argument.
    out : tuple of ndarray, optional
        Optional output arrays for the function values

    Returns
    -------
    eAi, eAip, eBi, eBip : 4-tuple of scalar or ndarray
        Exponentially scaled Airy functions eAi and eBi, and their derivatives
        eAip and eBip

    See Also
    --------
    airy

    Notes
    -----
    Wrapper for the AMOS [1]_ routines `zairy` and `zbiry`.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    We can compute exponentially scaled Airy functions and their derivatives:

    >>> import numpy as np
    >>> from scipy.special import airye
    >>> import matplotlib.pyplot as plt
    >>> z = np.linspace(0, 50, 500)
    >>> eAi, eAip, eBi, eBip = airye(z)
    >>> f, ax = plt.subplots(2, 1, sharex=True)
    >>> for ind, data in enumerate([[eAi, eAip, ["eAi", "eAip"]],
    ...                             [eBi, eBip, ["eBi", "eBip"]]]):
    ...     ax[ind].plot(z, data[0], "-r", z, data[1], "-b")
    ...     ax[ind].legend(data[2])
    ...     ax[ind].grid(True)
    >>> plt.show()

    We can compute these using usual non-scaled Airy functions by:

    >>> from scipy.special import airy
    >>> Ai, Aip, Bi, Bip = airy(z)
    >>> np.allclose(eAi, Ai * np.exp(2.0 / 3.0 * z * np.sqrt(z)))
    True
    >>> np.allclose(eAip, Aip * np.exp(2.0 / 3.0 * z * np.sqrt(z)))
    True
    >>> np.allclose(eBi, Bi * np.exp(-abs(np.real(2.0 / 3.0 * z * np.sqrt(z)))))
    True
    >>> np.allclose(eBip, Bip * np.exp(-abs(np.real(2.0 / 3.0 * z * np.sqrt(z)))))
    True

    Comparing non-scaled and exponentially scaled ones, the usual non-scaled
    function quickly underflows for large values, whereas the exponentially
    scaled function does not.

    >>> airy(200)
    (0.0, 0.0, nan, nan)
    >>> airye(200)
    (0.07501041684381093, -1.0609012305109042, 0.15003188417418148, 2.1215836725571093)

    )";

const char *bei_doc = R"(
    bei(x, out=None)

    Kelvin function bei.

    Defined as

    .. math::

        \mathrm{bei}(x) = \Im[J_0(x e^{3 \pi i / 4})]

    where :math:`J_0` is the Bessel function of the first kind of
    order zero (see `jv`). See [dlmf]_ for more details.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the Kelvin function.

    See Also
    --------
    ber : the corresponding real part
    beip : the derivative of bei
    jv : Bessel function of the first kind

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10.61

    Examples
    --------
    It can be expressed using Bessel functions.

    >>> import numpy as np
    >>> import scipy.special as sc
    >>> x = np.array([1.0, 2.0, 3.0, 4.0])
    >>> sc.jv(0, x * np.exp(3 * np.pi * 1j / 4)).imag
    array([0.24956604, 0.97229163, 1.93758679, 2.29269032])
    >>> sc.bei(x)
    array([0.24956604, 0.97229163, 1.93758679, 2.29269032])

    )";

const char *beip_doc = R"(
    beip(x, out=None)

    Derivative of the Kelvin function bei.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        The values of the derivative of bei.

    See Also
    --------
    bei

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10#PT5

    )";

const char *ber_doc = R"(
    ber(x, out=None)

    Kelvin function ber.

    Defined as

    .. math::

        \mathrm{ber}(x) = \Re[J_0(x e^{3 \pi i / 4})]

    where :math:`J_0` is the Bessel function of the first kind of
    order zero (see `jv`). See [dlmf]_ for more details.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the Kelvin function.

    See Also
    --------
    bei : the corresponding real part
    berp : the derivative of bei
    jv : Bessel function of the first kind

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10.61

    Examples
    --------
    It can be expressed using Bessel functions.

    >>> import numpy as np
    >>> import scipy.special as sc
    >>> x = np.array([1.0, 2.0, 3.0, 4.0])
    >>> sc.jv(0, x * np.exp(3 * np.pi * 1j / 4)).real
    array([ 0.98438178,  0.75173418, -0.22138025, -2.56341656])
    >>> sc.ber(x)
    array([ 0.98438178,  0.75173418, -0.22138025, -2.56341656])

    )";

const char *berp_doc = R"(
    berp(x, out=None)

    Derivative of the Kelvin function ber.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        The values of the derivative of ber.

    See Also
    --------
    ber

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10#PT5

    )";

const char *binom_doc = R"(
    binom(x, y, out=None)

    Binomial coefficient considered as a function of two real variables.

    For real arguments, the binomial coefficient is defined as

    .. math::

        \binom{x}{y} = \frac{\Gamma(x + 1)}{\Gamma(y + 1)\Gamma(x - y + 1)} =
            \frac{1}{(x + 1)\mathrm{B}(x - y + 1, y + 1)}

    Where :math:`\Gamma` is the Gamma function (`gamma`) and :math:`\mathrm{B}`
    is the Beta function (`beta`) [1]_.

    Parameters
    ----------
    x, y: array_like
       Real arguments to :math:`\binom{x}{y}`.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Value of binomial coefficient.

    See Also
    --------
    comb : The number of combinations of N things taken k at a time.

    Notes
    -----
    The Gamma function has poles at non-positive integers and tends to either
    positive or negative infinity depending on the direction on the real line
    from which a pole is approached. When considered as a function of two real
    variables, :math:`\binom{x}{y}` is thus undefined when `x` is a negative
    integer.  `binom` returns ``nan`` when ``x`` is a negative integer. This
    is the case even when ``x`` is a negative integer and ``y`` an integer,
    contrary to the usual convention for defining :math:`\binom{n}{k}` when it
    is considered as a function of two integer variables.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Binomial_coefficient

    Examples
    --------
    The following examples illustrate the ways in which `binom` differs from
    the function `comb`.

    >>> from scipy.special import binom, comb

    When ``exact=False`` and ``x`` and ``y`` are both positive, `comb` calls
    `binom` internally.

    >>> x, y = 3, 2
    >>> (binom(x, y), comb(x, y), comb(x, y, exact=True))
    (3.0, 3.0, 3)

    For larger values, `comb` with ``exact=True`` no longer agrees
    with `binom`.

    >>> x, y = 43, 23
    >>> (binom(x, y), comb(x, y), comb(x, y, exact=True))
    (960566918219.9999, 960566918219.9999, 960566918220)

    `binom` returns ``nan`` when ``x`` is a negative integer, but is otherwise
    defined for negative arguments. `comb` returns 0 whenever one of ``x`` or
    ``y`` is negative or ``x`` is less than ``y``.

    >>> x, y = -3, 2
    >>> (binom(x, y), comb(x, y))
    (nan, 0.0)

    >>> x, y = -3.1, 2.2
    >>> (binom(x, y), comb(x, y))
    (18.714147876804432, 0.0)

    >>> x, y = 2.2, 3.1
    >>> (binom(x, y), comb(x, y))
    (0.037399983365134115, 0.0)
    )";

const char *exp1_doc = R"(
    exp1(z, out=None)

    Exponential integral E1.

    For complex :math:`z \ne 0` the exponential integral can be defined as
    [1]_

    .. math::

       E_1(z) = \int_z^\infty \frac{e^{-t}}{t} dt,

    where the path of the integral does not cross the negative real
    axis or pass through the origin.

    Parameters
    ----------
    z: array_like
        Real or complex argument.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the exponential integral E1

    See Also
    --------
    expi : exponential integral :math:`Ei`
    expn : generalization of :math:`E_1`

    Notes
    -----
    For :math:`x > 0` it is related to the exponential integral
    :math:`Ei` (see `expi`) via the relation

    .. math::

       E_1(x) = -Ei(-x).

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 6.2.1
           https://dlmf.nist.gov/6.2#E1

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It has a pole at 0.

    >>> sc.exp1(0)
    inf

    It has a branch cut on the negative real axis.

    >>> sc.exp1(-1)
    nan
    >>> sc.exp1(complex(-1, 0))
    (-1.8951178163559368-3.141592653589793j)
    >>> sc.exp1(complex(-1, -0.0))
    (-1.8951178163559368+3.141592653589793j)

    It approaches 0 along the positive real axis.

    >>> sc.exp1([1, 10, 100, 1000])
    array([2.19383934e-01, 4.15696893e-06, 3.68359776e-46, 0.00000000e+00])

    It is related to `expi`.

    >>> x = np.array([1, 2, 3, 4])
    >>> sc.exp1(x)
    array([0.21938393, 0.04890051, 0.01304838, 0.00377935])
    >>> -sc.expi(-x)
    array([0.21938393, 0.04890051, 0.01304838, 0.00377935])

    )";

const char *expi_doc = R"(
    expi(x, out=None)

    Exponential integral Ei.

    For real :math:`x`, the exponential integral is defined as [1]_

    .. math::

        Ei(x) = \int_{-\infty}^x \frac{e^t}{t} dt.

    For :math:`x > 0` the integral is understood as a Cauchy principal
    value.

    It is extended to the complex plane by analytic continuation of
    the function on the interval :math:`(0, \infty)`. The complex
    variant has a branch cut on the negative real axis.

    Parameters
    ----------
    x : array_like
        Real or complex valued argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the exponential integral

    See Also
    --------
    exp1 : Exponential integral :math:`E_1`
    expn : Generalized exponential integral :math:`E_n`

    Notes
    -----
    The exponential integrals :math:`E_1` and :math:`Ei` satisfy the
    relation

    .. math::

        E_1(x) = -Ei(-x)

    for :math:`x > 0`.

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 6.2.5
           https://dlmf.nist.gov/6.2#E5

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is related to `exp1`.

    >>> x = np.array([1, 2, 3, 4])
    >>> -sc.expi(-x)
    array([0.21938393, 0.04890051, 0.01304838, 0.00377935])
    >>> sc.exp1(x)
    array([0.21938393, 0.04890051, 0.01304838, 0.00377935])

    The complex variant has a branch cut on the negative real axis.

    >>> sc.expi(-1 + 1e-12j)
    (-0.21938393439552062+3.1415926535894254j)
    >>> sc.expi(-1 - 1e-12j)
    (-0.21938393439552062-3.1415926535894254j)

    As the complex variant approaches the branch cut, the real parts
    approach the value of the real variant.

    >>> sc.expi(-1)
    -0.21938393439552062

    The SciPy implementation returns the real variant for complex
    values on the branch cut.

    >>> sc.expi(complex(-1, 0.0))
    (-0.21938393439552062-0j)
    >>> sc.expi(complex(-1, -0.0))
    (-0.21938393439552062-0j)

    )";

const char *expit_doc = R"(
    expit(x, out=None)

    Expit (a.k.a. logistic sigmoid) ufunc for ndarrays.

    The expit function, also known as the logistic sigmoid function, is
    defined as ``expit(x) = 1/(1+exp(-x))``.  It is the inverse of the
    logit function.

    Parameters
    ----------
    x : ndarray
        The ndarray to apply expit to element-wise.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        An ndarray of the same shape as x. Its entries
        are `expit` of the corresponding entry of x.

    See Also
    --------
    logit

    Notes
    -----
    As a ufunc expit takes a number of optional
    keyword arguments. For more information
    see `ufuncs <https://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_

    .. versionadded:: 0.10.0

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import expit, logit

    >>> expit([-np.inf, -1.5, 0, 1.5, np.inf])
    array([ 0.        ,  0.18242552,  0.5       ,  0.81757448,  1.        ])

    `logit` is the inverse of `expit`:

    >>> logit(expit([-2.5, 0, 3.1, 5.0]))
    array([-2.5,  0. ,  3.1,  5. ])

    Plot expit(x) for x in [-6, 6]:

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-6, 6, 121)
    >>> y = expit(x)
    >>> plt.plot(x, y)
    >>> plt.grid()
    >>> plt.xlim(-6, 6)
    >>> plt.xlabel('x')
    >>> plt.title('expit(x)')
    >>> plt.show()
    )";

const char *exprel_doc = R"(
    exprel(x, out=None)

    Relative error exponential, ``(exp(x) - 1)/x``.

    When `x` is near zero, ``exp(x)`` is near 1, so the numerical calculation
    of ``exp(x) - 1`` can suffer from catastrophic loss of precision.
    ``exprel(x)`` is implemented to avoid the loss of precision that occurs when
    `x` is near zero.

    Parameters
    ----------
    x : ndarray
        Input array.  `x` must contain real numbers.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        ``(exp(x) - 1)/x``, computed element-wise.

    See Also
    --------
    expm1

    Notes
    -----
    .. versionadded:: 0.17.0

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import exprel

    >>> exprel(0.01)
    1.0050167084168056
    >>> exprel([-0.25, -0.1, 0, 0.1, 0.25])
    array([ 0.88479687,  0.95162582,  1.        ,  1.05170918,  1.13610167])

    Compare ``exprel(5e-9)`` to the naive calculation.  The exact value
    is ``1.00000000250000000416...``.

    >>> exprel(5e-9)
    1.0000000025

    >>> (np.exp(5e-9) - 1)/5e-9
    0.99999999392252903
    )";

const char *gamma_doc = R"(
    gamma(z, out=None)

    gamma function.

    The gamma function is defined as

    .. math::

       \Gamma(z) = \int_0^\infty t^{z-1} e^{-t} dt

    for :math:`\Re(z) > 0` and is extended to the rest of the complex
    plane by analytic continuation. See [dlmf]_ for more details.

    Parameters
    ----------
    z : array_like
        Real or complex valued argument
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the gamma function

    Notes
    -----
    The gamma function is often referred to as the generalized
    factorial since :math:`\Gamma(n + 1) = n!` for natural numbers
    :math:`n`. More generally it satisfies the recurrence relation
    :math:`\Gamma(z + 1) = z \cdot \Gamma(z)` for complex :math:`z`,
    which, combined with the fact that :math:`\Gamma(1) = 1`, implies
    the above identity for :math:`z = n`.

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/5.2#E1

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import gamma, factorial

    >>> gamma([0, 0.5, 1, 5])
    array([         inf,   1.77245385,   1.        ,  24.        ])

    >>> z = 2.5 + 1j
    >>> gamma(z)
    (0.77476210455108352+0.70763120437959293j)
    >>> gamma(z+1), z*gamma(z)  # Recurrence property
    ((1.2292740569981171+2.5438401155000685j),
     (1.2292740569981158+2.5438401155000658j))

    >>> gamma(0.5)**2  # gamma(0.5) = sqrt(pi)
    3.1415926535897927

    Plot gamma(x) for real x

    >>> x = np.linspace(-3.5, 5.5, 2251)
    >>> y = gamma(x)

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, y, 'b', alpha=0.6, label='gamma(x)')
    >>> k = np.arange(1, 7)
    >>> plt.plot(k, factorial(k-1), 'k*', alpha=0.6,
    ...          label='(x-1)!, x = 1, 2, ...')
    >>> plt.xlim(-3.5, 5.5)
    >>> plt.ylim(-10, 25)
    >>> plt.grid()
    >>> plt.xlabel('x')
    >>> plt.legend(loc='lower right')
    >>> plt.show()

    )";

const char *gammaln_doc = R"(
    gammaln(x, out=None)

    Logarithm of the absolute value of the gamma function.

    Defined as

    .. math::

       \ln(\lvert\Gamma(x)\rvert)

    where :math:`\Gamma` is the gamma function. For more details on
    the gamma function, see [dlmf]_.

    Parameters
    ----------
    x : array_like
        Real argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the log of the absolute value of gamma

    See Also
    --------
    gammasgn : sign of the gamma function
    loggamma : principal branch of the logarithm of the gamma function

    Notes
    -----
    It is the same function as the Python standard library function
    :func:`math.lgamma`.

    When used in conjunction with `gammasgn`, this function is useful
    for working in logspace on the real axis without having to deal
    with complex numbers via the relation ``exp(gammaln(x)) =
    gammasgn(x) * gamma(x)``.

    For complex-valued log-gamma, use `loggamma` instead of `gammaln`.

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/5

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It has two positive zeros.

    >>> sc.gammaln([1, 2])
    array([0., 0.])

    It has poles at nonpositive integers.

    >>> sc.gammaln([0, -1, -2, -3, -4])
    array([inf, inf, inf, inf, inf])

    It asymptotically approaches ``x * log(x)`` (Stirling's formula).

    >>> x = np.array([1e10, 1e20, 1e40, 1e80])
    >>> sc.gammaln(x)
    array([2.20258509e+11, 4.50517019e+21, 9.11034037e+41, 1.83206807e+82])
    >>> x * np.log(x)
    array([2.30258509e+11, 4.60517019e+21, 9.21034037e+41, 1.84206807e+82])

    )";

const char *hankel1_doc = R"(
    hankel1(v, z, out=None)

    Hankel function of the first kind

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the Hankel function of the first kind.

    See Also
    --------
    hankel1e : ndarray
        This function with leading exponential behavior stripped off.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(1)}_v(z) =
              \frac{2}{\imath\pi} \exp(-\imath \pi v/2) K_v(z \exp(-\imath\pi/2))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(1)}_{-v}(z) = H^{(1)}_v(z) \exp(\imath\pi v)

    is used.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    )";

const char *hankel1e_doc = R"(
    hankel1e(v, z, out=None)

    Exponentially scaled Hankel function of the first kind

    Defined as::

        hankel1e(v, z) = hankel1(v, z) * exp(-1j * z)

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the exponentially scaled Hankel function.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(1)}_v(z) =
              \frac{2}{\imath\pi} \exp(-\imath \pi v/2) K_v(z \exp(-\imath\pi/2))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(1)}_{-v}(z) = H^{(1)}_v(z) \exp(\imath\pi v)

    is used.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    )";

const char *hankel2_doc = R"(
    hankel2(v, z, out=None)

    Hankel function of the second kind

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the Hankel function of the second kind.

    See Also
    --------
    hankel2e : this function with leading exponential behavior stripped off.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(2)}_v(z) =
              -\frac{2}{\imath\pi} \exp(\imath \pi v/2) K_v(z \exp(\imath\pi/2))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(2)}_{-v}(z) = H^{(2)}_v(z) \exp(-\imath\pi v)

    is used.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/
    )";

const char *hankel2e_doc = R"(
    hankel2e(v, z, out=None)

    Exponentially scaled Hankel function of the second kind

    Defined as::

        hankel2e(v, z) = hankel2(v, z) * exp(1j * z)

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the exponentially scaled Hankel function of the second kind.

    Notes
    -----
    A wrapper for the AMOS [1]_ routine `zbesh`, which carries out the
    computation using the relation,

    .. math:: H^{(2)}_v(z) = -\frac{2}{\imath\pi}
              \exp(\frac{\imath \pi v}{2}) K_v(z exp(\frac{\imath\pi}{2}))

    where :math:`K_v` is the modified Bessel function of the second kind.
    For negative orders, the relation

    .. math:: H^{(2)}_{-v}(z) = H^{(2)}_v(z) \exp(-\imath\pi v)

    is used.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    )";

const char *hyp2f1_doc = R"(
    hyp2f1(a, b, c, z, out=None)

    Gauss hypergeometric function 2F1(a, b; c; z)

    Parameters
    ----------
    a, b, c : array_like
        Arguments, should be real-valued.
    z : array_like
        Argument, real or complex.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    hyp2f1 : scalar or ndarray
        The values of the gaussian hypergeometric function.

    See Also
    --------
    hyp0f1 : confluent hypergeometric limit function.
    hyp1f1 : Kummer's (confluent hypergeometric) function.

    Notes
    -----
    This function is defined for :math:`|z| < 1` as

    .. math::

       \mathrm{hyp2f1}(a, b, c, z) = \sum_{n=0}^\infty
       \frac{(a)_n (b)_n}{(c)_n}\frac{z^n}{n!},

    and defined on the rest of the complex z-plane by analytic
    continuation [1]_.
    Here :math:`(\cdot)_n` is the Pochhammer symbol; see `poch`. When
    :math:`n` is an integer the result is a polynomial of degree :math:`n`.

    The implementation for complex values of ``z`` is described in [2]_,
    except for ``z`` in the region defined by

    .. math::

         0.9 <= \left|z\right| < 1.1,
         \left|1 - z\right| >= 0.9,
         \mathrm{real}(z) >= 0

    in which the implementation follows [4]_.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/15.2
    .. [2] S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996
    .. [3] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/
    .. [4] J.L. Lopez and N.M. Temme, "New series expansions of the Gauss
           hypergeometric function", Adv Comput Math 39, 349-365 (2013).
           https://doi.org/10.1007/s10444-012-9283-y

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It has poles when `c` is a negative integer.

    >>> sc.hyp2f1(1, 1, -2, 1)
    inf

    It is a polynomial when `a` or `b` is a negative integer.

    >>> a, b, c = -1, 1, 1.5
    >>> z = np.linspace(0, 1, 5)
    >>> sc.hyp2f1(a, b, c, z)
    array([1.        , 0.83333333, 0.66666667, 0.5       , 0.33333333])
    >>> 1 + a * b * z / c
    array([1.        , 0.83333333, 0.66666667, 0.5       , 0.33333333])

    It is symmetric in `a` and `b`.

    >>> a = np.linspace(0, 1, 5)
    >>> b = np.linspace(0, 1, 5)
    >>> sc.hyp2f1(a, b, 1, 0.5)
    array([1.        , 1.03997334, 1.1803406 , 1.47074441, 2.        ])
    >>> sc.hyp2f1(b, a, 1, 0.5)
    array([1.        , 1.03997334, 1.1803406 , 1.47074441, 2.        ])

    It contains many other functions as special cases.

    >>> z = 0.5
    >>> sc.hyp2f1(1, 1, 2, z)
    1.3862943611198901
    >>> -np.log(1 - z) / z
    1.3862943611198906

    >>> sc.hyp2f1(0.5, 1, 1.5, z**2)
    1.098612288668109
    >>> np.log((1 + z) / (1 - z)) / (2 * z)
    1.0986122886681098

    >>> sc.hyp2f1(0.5, 1, 1.5, -z**2)
    0.9272952180016117
    >>> np.arctan(z) / z
    0.9272952180016122
    )";

const char *it2i0k0_doc = R"(
    it2i0k0(x, out=None)

    Integrals related to modified Bessel functions of order 0.

    Computes the integrals

    .. math::

        \int_0^x \frac{I_0(t) - 1}{t} dt \\
        \int_x^\infty \frac{K_0(t)}{t} dt.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the integrals.
    out : tuple of ndarrays, optional
        Optional output arrays for the function results.

    Returns
    -------
    ii0 : scalar or ndarray
        The integral for `i0`
    ik0 : scalar or ndarray
        The integral for `k0`

    References
    ----------
    .. [1] S. Zhang and J.M. Jin, "Computation of Special Functions",
           Wiley 1996

    Examples
    --------
    Evaluate the functions at one point.

    >>> from scipy.special import it2i0k0
    >>> int_i, int_k = it2i0k0(1.)
    >>> int_i, int_k
    (0.12897944249456852, 0.2085182909001295)

    Evaluate the functions at several points.

    >>> import numpy as np
    >>> points = np.array([0.5, 1.5, 3.])
    >>> int_i, int_k = it2i0k0(points)
    >>> int_i, int_k
    (array([0.03149527, 0.30187149, 1.50012461]),
     array([0.66575102, 0.0823715 , 0.00823631]))

    Plot the functions from 0 to 5.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 5., 1000)
    >>> int_i, int_k = it2i0k0(x)
    >>> ax.plot(x, int_i, label=r"$\int_0^x \frac{I_0(t)-1}{t}\,dt$")
    >>> ax.plot(x, int_k, label=r"$\int_x^{\infty} \frac{K_0(t)}{t}\,dt$")
    >>> ax.legend()
    >>> ax.set_ylim(0, 10)
    >>> plt.show()
    )";

const char *it2j0y0_doc = R"(
    it2j0y0(x, out=None)

    Integrals related to Bessel functions of the first kind of order 0.

    Computes the integrals

    .. math::

        \int_0^x \frac{1 - J_0(t)}{t} dt \\
        \int_x^\infty \frac{Y_0(t)}{t} dt.

    For more on :math:`J_0` and :math:`Y_0` see `j0` and `y0`.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the integrals.
    out : tuple of ndarrays, optional
        Optional output arrays for the function results.

    Returns
    -------
    ij0 : scalar or ndarray
        The integral for `j0`
    iy0 : scalar or ndarray
        The integral for `y0`

    References
    ----------
    .. [1] S. Zhang and J.M. Jin, "Computation of Special Functions",
           Wiley 1996

    Examples
    --------
    Evaluate the functions at one point.

    >>> from scipy.special import it2j0y0
    >>> int_j, int_y = it2j0y0(1.)
    >>> int_j, int_y
    (0.12116524699506871, 0.39527290169929336)

    Evaluate the functions at several points.

    >>> import numpy as np
    >>> points = np.array([0.5, 1.5, 3.])
    >>> int_j, int_y = it2j0y0(points)
    >>> int_j, int_y
    (array([0.03100699, 0.26227724, 0.85614669]),
     array([ 0.26968854,  0.29769696, -0.02987272]))

    Plot the functions from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> int_j, int_y = it2j0y0(x)
    >>> ax.plot(x, int_j, label=r"$\int_0^x \frac{1-J_0(t)}{t}\,dt$")
    >>> ax.plot(x, int_y, label=r"$\int_x^{\infty} \frac{Y_0(t)}{t}\,dt$")
    >>> ax.legend()
    >>> ax.set_ylim(-2.5, 2.5)
    >>> plt.show()
    )";

const char *it2struve0_doc = R"(
    it2struve0(x, out=None)

    Integral related to the Struve function of order 0.

    Returns the integral,

    .. math::
        \int_x^\infty \frac{H_0(t)}{t}\,dt

    where :math:`H_0` is the Struve function of order 0.

    Parameters
    ----------
    x : array_like
        Lower limit of integration.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    I : scalar or ndarray
        The value of the integral.

    See Also
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
           https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html

    Examples
    --------
    Evaluate the function at one point.

    >>> import numpy as np
    >>> from scipy.special import it2struve0
    >>> it2struve0(1.)
    0.9571973506383524

    Evaluate the function at several points by supplying
    an array for `x`.

    >>> points = np.array([1., 2., 3.5])
    >>> it2struve0(points)
    array([0.95719735, 0.46909296, 0.10366042])

    Plot the function from -10 to 10.

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-10., 10., 1000)
    >>> it2struve0_values = it2struve0(x)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, it2struve0_values)
    >>> ax.set_xlabel(r'$x$')
    >>> ax.set_ylabel(r'$\int_x^{\infty}\frac{H_0(t)}{t}\,dt$')
    >>> plt.show()
    )";

const char *itairy_doc = R"(
    itairy(x, out=None)

    Integrals of Airy functions

    Calculates the integrals of Airy functions from 0 to `x`.

    Parameters
    ----------

    x : array_like
        Upper limit of integration (float).
    out : tuple of ndarray, optional
        Optional output arrays for the function values

    Returns
    -------
    Apt : scalar or ndarray
        Integral of Ai(t) from 0 to x.
    Bpt : scalar or ndarray
        Integral of Bi(t) from 0 to x.
    Ant : scalar or ndarray
        Integral of Ai(-t) from 0 to x.
    Bnt : scalar or ndarray
        Integral of Bi(-t) from 0 to x.

    Notes
    -----

    Wrapper for a Fortran routine created by Shanjie Zhang and Jianming
    Jin [1]_.

    References
    ----------

    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html

    Examples
    --------
    Compute the functions at ``x=1.``.

    >>> import numpy as np
    >>> from scipy.special import itairy
    >>> import matplotlib.pyplot as plt
    >>> apt, bpt, ant, bnt = itairy(1.)
    >>> apt, bpt, ant, bnt
    (0.23631734191710949,
     0.8727691167380077,
     0.46567398346706845,
     0.3730050096342943)

    Compute the functions at several points by providing a NumPy array for `x`.

    >>> x = np.array([1., 1.5, 2.5, 5])
    >>> apt, bpt, ant, bnt = itairy(x)
    >>> apt, bpt, ant, bnt
    (array([0.23631734, 0.28678675, 0.324638  , 0.33328759]),
     array([  0.87276912,   1.62470809,   5.20906691, 321.47831857]),
     array([0.46567398, 0.72232876, 0.93187776, 0.7178822 ]),
     array([ 0.37300501,  0.35038814, -0.02812939,  0.15873094]))

    Plot the functions from -10 to 10.

    >>> x = np.linspace(-10, 10, 500)
    >>> apt, bpt, ant, bnt = itairy(x)
    >>> fig, ax = plt.subplots(figsize=(6, 5))
    >>> ax.plot(x, apt, label=r"$\int_0^x\, Ai(t)\, dt$")
    >>> ax.plot(x, bpt, ls="dashed", label=r"$\int_0^x\, Bi(t)\, dt$")
    >>> ax.plot(x, ant, ls="dashdot", label=r"$\int_0^x\, Ai(-t)\, dt$")
    >>> ax.plot(x, bnt, ls="dotted", label=r"$\int_0^x\, Bi(-t)\, dt$")
    >>> ax.set_ylim(-2, 1.5)
    >>> ax.legend(loc="lower right")
    >>> plt.show()
    )";

const char *iti0k0_doc = R"(
    iti0k0(x, out=None)

    Integrals of modified Bessel functions of order 0.

    Computes the integrals

    .. math::

        \int_0^x I_0(t) dt \\
        \int_0^x K_0(t) dt.

    For more on :math:`I_0` and :math:`K_0` see `i0` and `k0`.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the integrals.
    out : tuple of ndarrays, optional
        Optional output arrays for the function results.

    Returns
    -------
    ii0 : scalar or ndarray
        The integral for `i0`
    ik0 : scalar or ndarray
        The integral for `k0`

    References
    ----------
    .. [1] S. Zhang and J.M. Jin, "Computation of Special Functions",
           Wiley 1996

    Examples
    --------
    Evaluate the functions at one point.

    >>> from scipy.special import iti0k0
    >>> int_i, int_k = iti0k0(1.)
    >>> int_i, int_k
    (1.0865210970235892, 1.2425098486237771)

    Evaluate the functions at several points.

    >>> import numpy as np
    >>> points = np.array([0., 1.5, 3.])
    >>> int_i, int_k = iti0k0(points)
    >>> int_i, int_k
    (array([0.        , 1.80606937, 6.16096149]),
     array([0.        , 1.39458246, 1.53994809]))

    Plot the functions from 0 to 5.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 5., 1000)
    >>> int_i, int_k = iti0k0(x)
    >>> ax.plot(x, int_i, label=r"$\int_0^x I_0(t)\,dt$")
    >>> ax.plot(x, int_k, label=r"$\int_0^x K_0(t)\,dt$")
    >>> ax.legend()
    >>> plt.show()
    )";

const char *itj0y0_doc = R"(
    itj0y0(x, out=None)

    Integrals of Bessel functions of the first kind of order 0.

    Computes the integrals

    .. math::

        \int_0^x J_0(t) dt \\
        \int_0^x Y_0(t) dt.

    For more on :math:`J_0` and :math:`Y_0` see `j0` and `y0`.

    Parameters
    ----------
    x : array_like
        Values at which to evaluate the integrals.
    out : tuple of ndarrays, optional
        Optional output arrays for the function results.

    Returns
    -------
    ij0 : scalar or ndarray
        The integral of `j0`
    iy0 : scalar or ndarray
        The integral of `y0`

    References
    ----------
    .. [1] S. Zhang and J.M. Jin, "Computation of Special Functions",
           Wiley 1996

    Examples
    --------
    Evaluate the functions at one point.

    >>> from scipy.special import itj0y0
    >>> int_j, int_y = itj0y0(1.)
    >>> int_j, int_y
    (0.9197304100897596, -0.637069376607422)

    Evaluate the functions at several points.

    >>> import numpy as np
    >>> points = np.array([0., 1.5, 3.])
    >>> int_j, int_y = itj0y0(points)
    >>> int_j, int_y
    (array([0.        , 1.24144951, 1.38756725]),
     array([ 0.        , -0.51175903,  0.19765826]))

    Plot the functions from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> int_j, int_y = itj0y0(x)
    >>> ax.plot(x, int_j, label=r"$\int_0^x J_0(t)\,dt$")
    >>> ax.plot(x, int_y, label=r"$\int_0^x Y_0(t)\,dt$")
    >>> ax.legend()
    >>> plt.show()

    )";

const char *itmodstruve0_doc = R"(
    itmodstruve0(x, out=None)

    Integral of the modified Struve function of order 0.

    .. math::
        I = \int_0^x L_0(t)\,dt

    Parameters
    ----------
    x : array_like
        Upper limit of integration (float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    I : scalar or ndarray
        The integral of :math:`L_0` from 0 to `x`.

    See Also
    --------
    modstruve: Modified Struve function which is integrated by this function

    Notes
    -----
    Wrapper for a Fortran routine created by Shanjie Zhang and Jianming
    Jin [1]_.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html

    Examples
    --------
    Evaluate the function at one point.

    >>> import numpy as np
    >>> from scipy.special import itmodstruve0
    >>> itmodstruve0(1.)
    0.3364726286440384

    Evaluate the function at several points by supplying
    an array for `x`.

    >>> points = np.array([1., 2., 3.5])
    >>> itmodstruve0(points)
    array([0.33647263, 1.588285  , 7.60382578])

    Plot the function from -10 to 10.

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-10., 10., 1000)
    >>> itmodstruve0_values = itmodstruve0(x)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, itmodstruve0_values)
    >>> ax.set_xlabel(r'$x$')
    >>> ax.set_ylabel(r'$\int_0^xL_0(t)\,dt$')
    >>> plt.show()
    )";

const char *itstruve0_doc = R"(
    itstruve0(x, out=None)

    Integral of the Struve function of order 0.

    .. math::
        I = \int_0^x H_0(t)\,dt

    Parameters
    ----------
    x : array_like
        Upper limit of integration (float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    I : scalar or ndarray
        The integral of :math:`H_0` from 0 to `x`.

    See Also
    --------
    struve: Function which is integrated by this function

    Notes
    -----
    Wrapper for a Fortran routine created by Shanjie Zhang and Jianming
    Jin [1]_.

    References
    ----------
    .. [1] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html

    Examples
    --------
    Evaluate the function at one point.

    >>> import numpy as np
    >>> from scipy.special import itstruve0
    >>> itstruve0(1.)
    0.30109042670805547

    Evaluate the function at several points by supplying
    an array for `x`.

    >>> points = np.array([1., 2., 3.5])
    >>> itstruve0(points)
    array([0.30109043, 1.01870116, 1.96804581])

    Plot the function from -20 to 20.

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-20., 20., 1000)
    >>> istruve0_values = itstruve0(x)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, istruve0_values)
    >>> ax.set_xlabel(r'$x$')
    >>> ax.set_ylabel(r'$\int_0^{x}H_0(t)\,dt$')
    >>> plt.show()
    )";

const char *iv_doc = R"(
    iv(v, z, out=None)

    Modified Bessel function of the first kind of real order.

    Parameters
    ----------
    v : array_like
        Order. If `z` is of real type and negative, `v` must be integer
        valued.
    z : array_like of float or complex
        Argument.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the modified Bessel function.

    See Also
    --------
    ive : This function with leading exponential behavior stripped off.
    i0 : Faster version of this function for order 0.
    i1 : Faster version of this function for order 1.

    Notes
    -----
    For real `z` and :math:`v \in [-50, 50]`, the evaluation is carried out
    using Temme's method [1]_.  For larger orders, uniform asymptotic
    expansions are applied.

    For complex `z` and positive `v`, the AMOS [2]_ `zbesi` routine is
    called. It uses a power series for small `z`, the asymptotic expansion
    for large `abs(z)`, the Miller algorithm normalized by the Wronskian
    and a Neumann series for intermediate magnitudes, and the uniform
    asymptotic expansions for :math:`I_v(z)` and :math:`J_v(z)` for large
    orders. Backward recurrence is used to generate sequences or reduce
    orders when necessary.

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
    .. [1] Temme, Journal of Computational Physics, vol 21, 343 (1976)
    .. [2] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    Evaluate the function of order 0 at one point.

    >>> from scipy.special import iv
    >>> iv(0, 1.)
    1.2660658777520084

    Evaluate the function at one point for different orders.

    >>> iv(0, 1.), iv(1, 1.), iv(1.5, 1.)
    (1.2660658777520084, 0.565159103992485, 0.2935253263474798)

    The evaluation for different orders can be carried out in one call by
    providing a list or NumPy array as argument for the `v` parameter:

    >>> iv([0, 1, 1.5], 1.)
    array([1.26606588, 0.5651591 , 0.29352533])

    Evaluate the function at several points for order 0 by providing an
    array for `z`.

    >>> import numpy as np
    >>> points = np.array([-2., 0., 3.])
    >>> iv(0, points)
    array([2.2795853 , 1.        , 4.88079259])

    If `z` is an array, the order parameter `v` must be broadcastable to
    the correct shape if different orders shall be computed in one call.
    To calculate the orders 0 and 1 for an 1D array:

    >>> orders = np.array([[0], [1]])
    >>> orders.shape
    (2, 1)

    >>> iv(orders, points)
    array([[ 2.2795853 ,  1.        ,  4.88079259],
           [-1.59063685,  0.        ,  3.95337022]])

    Plot the functions of order 0 to 3 from -5 to 5.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-5., 5., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, iv(i, x), label=f'$I_{i!r}$')
    >>> ax.legend()
    >>> plt.show()

    )";

const char *iv_ratio_doc = R"(
    _iv_ratio(v, x, out=None)

    Internal function, do not use.

    Return `iv(v, x) / iv(v-1, x)` for `v >= 1` and `x >= 0`.
    )";

const char *ive_doc = R"(
    ive(v, z, out=None)

    Exponentially scaled modified Bessel function of the first kind.

    Defined as::

        ive(v, z) = iv(v, z) * exp(-abs(z.real))

    For imaginary numbers without a real part, returns the unscaled
    Bessel function of the first kind `iv`.

    Parameters
    ----------
    v : array_like of float
        Order.
    z : array_like of float or complex
        Argument.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the exponentially scaled modified Bessel function.

    See Also
    --------
    iv: Modified Bessel function of the first kind
    i0e: Faster implementation of this function for order 0
    i1e: Faster implementation of this function for order 1

    Notes
    -----
    For positive `v`, the AMOS [1]_ `zbesi` routine is called. It uses a
    power series for small `z`, the asymptotic expansion for large
    `abs(z)`, the Miller algorithm normalized by the Wronskian and a
    Neumann series for intermediate magnitudes, and the uniform asymptotic
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

    `ive` is useful for large arguments `z`: for these, `iv` easily overflows,
    while `ive` does not due to the exponential scaling.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    In the following example `iv` returns infinity whereas `ive` still returns
    a finite number.

    >>> from scipy.special import iv, ive
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> iv(3, 1000.), ive(3, 1000.)
    (inf, 0.01256056218254712)

    Evaluate the function at one point for different orders by
    providing a list or NumPy array as argument for the `v` parameter:

    >>> ive([0, 1, 1.5], 1.)
    array([0.46575961, 0.20791042, 0.10798193])

    Evaluate the function at several points for order 0 by providing an
    array for `z`.

    >>> points = np.array([-2., 0., 3.])
    >>> ive(0, points)
    array([0.30850832, 1.        , 0.24300035])

    Evaluate the function at several points for different orders by
    providing arrays for both `v` for `z`. Both arrays have to be
    broadcastable to the correct shape. To calculate the orders 0, 1
    and 2 for a 1D array of points:

    >>> ive([[0], [1], [2]], points)
    array([[ 0.30850832,  1.        ,  0.24300035],
           [-0.21526929,  0.        ,  0.19682671],
           [ 0.09323903,  0.        ,  0.11178255]])

    Plot the functions of order 0 to 3 from -5 to 5.

    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-5., 5., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, ive(i, x), label=fr'$I_{i!r}(z)\cdot e^{{-|z|}}$')
    >>> ax.legend()
    >>> ax.set_xlabel(r"$z$")
    >>> plt.show()
    )";

const char *jv_doc = R"(
    jv(v, z, out=None)

    Bessel function of the first kind of real order and complex argument.

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    J : scalar or ndarray
        Value of the Bessel function, :math:`J_v(z)`.

    See Also
    --------
    jve : :math:`J_v` with leading exponential behavior stripped off.
    spherical_jn : spherical Bessel functions.
    j0 : faster version of this function for order 0.
    j1 : faster version of this function for order 1.

    Notes
    -----
    For positive `v` values, the computation is carried out using the AMOS
    [1]_ `zbesj` routine, which exploits the connection to the modified
    Bessel function :math:`I_v`,

    .. math::
        J_v(z) = \exp(v\pi\imath/2) I_v(-\imath z)\qquad (\Im z > 0)

        J_v(z) = \exp(-v\pi\imath/2) I_v(\imath z)\qquad (\Im z < 0)

    For negative `v` values the formula,

    .. math:: J_{-v}(z) = J_v(z) \cos(\pi v) - Y_v(z) \sin(\pi v)

    is used, where :math:`Y_v(z)` is the Bessel function of the second
    kind, computed using the AMOS routine `zbesy`.  Note that the second
    term is exactly zero for integer `v`; to improve accuracy the second
    term is explicitly omitted for `v` values such that `v = floor(v)`.

    Not to be confused with the spherical Bessel functions (see `spherical_jn`).

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    Evaluate the function of order 0 at one point.

    >>> from scipy.special import jv
    >>> jv(0, 1.)
    0.7651976865579666

    Evaluate the function at one point for different orders.

    >>> jv(0, 1.), jv(1, 1.), jv(1.5, 1.)
    (0.7651976865579666, 0.44005058574493355, 0.24029783912342725)

    The evaluation for different orders can be carried out in one call by
    providing a list or NumPy array as argument for the `v` parameter:

    >>> jv([0, 1, 1.5], 1.)
    array([0.76519769, 0.44005059, 0.24029784])

    Evaluate the function at several points for order 0 by providing an
    array for `z`.

    >>> import numpy as np
    >>> points = np.array([-2., 0., 3.])
    >>> jv(0, points)
    array([ 0.22389078,  1.        , -0.26005195])

    If `z` is an array, the order parameter `v` must be broadcastable to
    the correct shape if different orders shall be computed in one call.
    To calculate the orders 0 and 1 for an 1D array:

    >>> orders = np.array([[0], [1]])
    >>> orders.shape
    (2, 1)

    >>> jv(orders, points)
    array([[ 0.22389078,  1.        , -0.26005195],
           [-0.57672481,  0.        ,  0.33905896]])

    Plot the functions of order 0 to 3 from -10 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-10., 10., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, jv(i, x), label=f'$J_{i!r}$')
    >>> ax.legend()
    >>> plt.show()

    )";

const char *jve_doc = R"(
    jve(v, z, out=None)

    Exponentially scaled Bessel function of the first kind of order `v`.

    Defined as::

        jve(v, z) = jv(v, z) * exp(-abs(z.imag))

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    J : scalar or ndarray
        Value of the exponentially scaled Bessel function.

    See Also
    --------
    jv: Unscaled Bessel function of the first kind

    Notes
    -----
    For positive `v` values, the computation is carried out using the AMOS
    [1]_ `zbesj` routine, which exploits the connection to the modified
    Bessel function :math:`I_v`,

    .. math::
        J_v(z) = \exp(v\pi\imath/2) I_v(-\imath z)\qquad (\Im z > 0)

        J_v(z) = \exp(-v\pi\imath/2) I_v(\imath z)\qquad (\Im z < 0)

    For negative `v` values the formula,

    .. math:: J_{-v}(z) = J_v(z) \cos(\pi v) - Y_v(z) \sin(\pi v)

    is used, where :math:`Y_v(z)` is the Bessel function of the second
    kind, computed using the AMOS routine `zbesy`.  Note that the second
    term is exactly zero for integer `v`; to improve accuracy the second
    term is explicitly omitted for `v` values such that `v = floor(v)`.

    Exponentially scaled Bessel functions are useful for large arguments `z`:
    for these, the unscaled Bessel functions can easily under-or overflow.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    Compare the output of `jv` and `jve` for large complex arguments for `z`
    by computing their values for order ``v=1`` at ``z=1000j``. We see that
    `jv` overflows but `jve` returns a finite number:

    >>> import numpy as np
    >>> from scipy.special import jv, jve
    >>> v = 1
    >>> z = 1000j
    >>> jv(v, z), jve(v, z)
    ((inf+infj), (7.721967686709077e-19+0.012610930256928629j))

    For real arguments for `z`, `jve` returns the same as `jv`.

    >>> v, z = 1, 1000
    >>> jv(v, z), jve(v, z)
    (0.004728311907089523, 0.004728311907089523)

    The function can be evaluated for several orders at the same time by
    providing a list or NumPy array for `v`:

    >>> jve([1, 3, 5], 1j)
    array([1.27304208e-17+2.07910415e-01j, -4.99352086e-19-8.15530777e-03j,
           6.11480940e-21+9.98657141e-05j])

    In the same way, the function can be evaluated at several points in one
    call by providing a list or NumPy array for `z`:

    >>> jve(1, np.array([1j, 2j, 3j]))
    array([1.27308412e-17+0.20791042j, 1.31814423e-17+0.21526929j,
           1.20521602e-17+0.19682671j])

    It is also possible to evaluate several orders at several points
    at the same time by providing arrays for `v` and `z` with
    compatible shapes for broadcasting. Compute `jve` for two different orders
    `v` and three points `z` resulting in a 2x3 array.

    >>> v = np.array([[1], [3]])
    >>> z = np.array([1j, 2j, 3j])
    >>> v.shape, z.shape
    ((2, 1), (3,))

    >>> jve(v, z)
    array([[1.27304208e-17+0.20791042j,  1.31810070e-17+0.21526929j,
            1.20517622e-17+0.19682671j],
           [-4.99352086e-19-0.00815531j, -1.76289571e-18-0.02879122j,
            -2.92578784e-18-0.04778332j]])
    )";

const char *kei_doc = R"(
    kei(x, out=None)

    Kelvin function kei.

    Defined as

    .. math::

        \mathrm{kei}(x) = \Im[K_0(x e^{\pi i / 4})]

    where :math:`K_0` is the modified Bessel function of the second
    kind (see `kv`). See [dlmf]_ for more details.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the Kelvin function.

    See Also
    --------
    ker : the corresponding real part
    keip : the derivative of kei
    kv : modified Bessel function of the second kind

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10.61

    Examples
    --------
    It can be expressed using the modified Bessel function of the
    second kind.

    >>> import numpy as np
    >>> import scipy.special as sc
    >>> x = np.array([1.0, 2.0, 3.0, 4.0])
    >>> sc.kv(0, x * np.exp(np.pi * 1j / 4)).imag
    array([-0.49499464, -0.20240007, -0.05112188,  0.0021984 ])
    >>> sc.kei(x)
    array([-0.49499464, -0.20240007, -0.05112188,  0.0021984 ])

    )";

const char *keip_doc = R"(
    keip(x, out=None)

    Derivative of the Kelvin function kei.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        The values of the derivative of kei.

    See Also
    --------
    kei

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10#PT5

    )";

const char *kelvin_doc = R"(
    kelvin(x, out=None)

    Kelvin functions as complex numbers

    Parameters
    ----------
    x : array_like
        Argument
    out : tuple of ndarray, optional
        Optional output arrays for the function values

    Returns
    -------
    Be, Ke, Bep, Kep : 4-tuple of scalar or ndarray
        The tuple (Be, Ke, Bep, Kep) contains complex numbers
        representing the real and imaginary Kelvin functions and their
        derivatives evaluated at `x`.  For example, kelvin(x)[0].real =
        ber x and kelvin(x)[0].imag = bei x with similar relationships
        for ker and kei.
    )";

const char *ker_doc = R"(
    ker(x, out=None)

    Kelvin function ker.

    Defined as

    .. math::

        \mathrm{ker}(x) = \Re[K_0(x e^{\pi i / 4})]

    Where :math:`K_0` is the modified Bessel function of the second
    kind (see `kv`). See [dlmf]_ for more details.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the Kelvin function.

    See Also
    --------
    kei : the corresponding imaginary part
    kerp : the derivative of ker
    kv : modified Bessel function of the second kind

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10.61

    Examples
    --------
    It can be expressed using the modified Bessel function of the
    second kind.

    >>> import numpy as np
    >>> import scipy.special as sc
    >>> x = np.array([1.0, 2.0, 3.0, 4.0])
    >>> sc.kv(0, x * np.exp(np.pi * 1j / 4)).real
    array([ 0.28670621, -0.04166451, -0.06702923, -0.03617885])
    >>> sc.ker(x)
    array([ 0.28670621, -0.04166451, -0.06702923, -0.03617885])

    )";

const char *kerp_doc = R"(
    kerp(x, out=None)

    Derivative of the Kelvin function ker.

    Parameters
    ----------
    x : array_like
        Real argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the derivative of ker.

    See Also
    --------
    ker

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10#PT5

    )";

const char *kv_doc = R"(
    kv(v, z, out=None)

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
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The results. Note that input must be of complex type to get complex
        output, e.g. ``kv(3, -2+0j)`` instead of ``kv(3, -2)``.

    See Also
    --------
    kve : This function with leading exponential behavior stripped off.
    kvp : Derivative of this function

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
    .. [3] NIST Digital Library of Mathematical Functions,
           Eq. 10.25.E3. https://dlmf.nist.gov/10.25.E3

    Examples
    --------
    Plot the function of several orders for real input:

    >>> import numpy as np
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

    )";

const char *kve_doc = R"(
    kve(v, z, out=None)

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
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The exponentially scaled modified Bessel function of the second kind.

    See Also
    --------
    kv : This function without exponential scaling.
    k0e : Faster version of this function for order 0.
    k1e : Faster version of this function for order 1.

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

    Examples
    --------
    In the following example `kv` returns 0 whereas `kve` still returns
    a useful finite number.

    >>> import numpy as np
    >>> from scipy.special import kv, kve
    >>> import matplotlib.pyplot as plt
    >>> kv(3, 1000.), kve(3, 1000.)
    (0.0, 0.03980696128440973)

    Evaluate the function at one point for different orders by
    providing a list or NumPy array as argument for the `v` parameter:

    >>> kve([0, 1, 1.5], 1.)
    array([1.14446308, 1.63615349, 2.50662827])

    Evaluate the function at several points for order 0 by providing an
    array for `z`.

    >>> points = np.array([1., 3., 10.])
    >>> kve(0, points)
    array([1.14446308, 0.6977616 , 0.39163193])

    Evaluate the function at several points for different orders by
    providing arrays for both `v` for `z`. Both arrays have to be
    broadcastable to the correct shape. To calculate the orders 0, 1
    and 2 for a 1D array of points:

    >>> kve([[0], [1], [2]], points)
    array([[1.14446308, 0.6977616 , 0.39163193],
           [1.63615349, 0.80656348, 0.41076657],
           [4.41677005, 1.23547058, 0.47378525]])

    Plot the functions of order 0 to 3 from 0 to 5.

    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 5., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, kve(i, x), label=fr'$K_{i!r}(z)\cdot e^z$')
    >>> ax.legend()
    >>> ax.set_xlabel(r"$z$")
    >>> ax.set_ylim(0, 4)
    >>> ax.set_xlim(0, 5)
    >>> plt.show()
    )";

const char *lambertw_doc = R"(
    lambertw(z, k=0, tol=1e-8)

    Lambert W function.

    The Lambert W function `W(z)` is defined as the inverse function
    of ``w * exp(w)``. In other words, the value of ``W(z)`` is
    such that ``z = W(z) * exp(W(z))`` for any complex number
    ``z``.

    The Lambert W function is a multivalued function with infinitely
    many branches. Each branch gives a separate solution of the
    equation ``z = w exp(w)``. Here, the branches are indexed by the
    integer `k`.

    Parameters
    ----------
    z : array_like
        Input argument.
    k : int, optional
        Branch index.
    tol : float, optional
        Evaluation tolerance.

    Returns
    -------
    w : array
        `w` will have the same shape as `z`.

    See Also
    --------
    wrightomega : the Wright Omega function

    Notes
    -----
    All branches are supported by `lambertw`:

    * ``lambertw(z)`` gives the principal solution (branch 0)
    * ``lambertw(z, k)`` gives the solution on branch `k`

    The Lambert W function has two partially real branches: the
    principal branch (`k = 0`) is real for real ``z > -1/e``, and the
    ``k = -1`` branch is real for ``-1/e < z < 0``. All branches except
    ``k = 0`` have a logarithmic singularity at ``z = 0``.

    **Possible issues**

    The evaluation can become inaccurate very close to the branch point
    at ``-1/e``. In some corner cases, `lambertw` might currently
    fail to converge, or can end up on the wrong branch.

    **Algorithm**

    Halley's iteration is used to invert ``w * exp(w)``, using a first-order
    asymptotic approximation (O(log(w)) or `O(w)`) as the initial estimate.

    The definition, implementation and choice of branches is based on [2]_.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Lambert_W_function
    .. [2] Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
       (1996) 329-359.
       https://cs.uwaterloo.ca/research/tr/1993/03/W.pdf

    Examples
    --------
    The Lambert W function is the inverse of ``w exp(w)``:

    >>> import numpy as np
    >>> from scipy.special import lambertw
    >>> w = lambertw(1)
    >>> w
    (0.56714329040978384+0j)
    >>> w * np.exp(w)
    (1.0+0j)

    Any branch gives a valid inverse:

    >>> w = lambertw(1, k=3)
    >>> w
    (-2.8535817554090377+17.113535539412148j)
    >>> w*np.exp(w)
    (1.0000000000000002+1.609823385706477e-15j)

    **Applications to equation-solving**

    The Lambert W function may be used to solve various kinds of
    equations.  We give two examples here.

    First, the function can be used to solve implicit equations of the
    form

        :math:`x = a + b e^{c x}`

    for :math:`x`.  We assume :math:`c` is not zero.  After a little
    algebra, the equation may be written

        :math:`z e^z = -b c e^{a c}`

    where :math:`z = c (a - x)`.  :math:`z` may then be expressed using
    the Lambert W function

        :math:`z = W(-b c e^{a c})`

    giving

        :math:`x = a - W(-b c e^{a c})/c`

    For example,

    >>> a = 3
    >>> b = 2
    >>> c = -0.5

    The solution to :math:`x = a + b e^{c x}` is:

    >>> x = a - lambertw(-b*c*np.exp(a*c))/c
    >>> x
    (3.3707498368978794+0j)

    Verify that it solves the equation:

    >>> a + b*np.exp(c*x)
    (3.37074983689788+0j)

    The Lambert W function may also be used find the value of the infinite
    power tower :math:`z^{z^{z^{\ldots}}}`:

    >>> def tower(z, n):
    ...     if n == 0:
    ...         return z
    ...     return z ** tower(z, n-1)
    ...
    >>> tower(0.5, 100)
    0.641185744504986
    >>> -lambertw(-np.log(0.5)) / np.log(0.5)
    (0.64118574450498589+0j)
    )";

const char *loggamma_doc = R"(
    loggamma(z, out=None)

    Principal branch of the logarithm of the gamma function.

    Defined to be :math:`\log(\Gamma(x))` for :math:`x > 0` and
    extended to the complex plane by analytic continuation. The
    function has a single branch cut on the negative real axis.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    z : array_like
        Values in the complex plane at which to compute ``loggamma``
    out : ndarray, optional
        Output array for computed values of ``loggamma``

    Returns
    -------
    loggamma : scalar or ndarray
        Values of ``loggamma`` at z.

    See Also
    --------
    gammaln : logarithm of the absolute value of the gamma function
    gammasgn : sign of the gamma function

    Notes
    -----
    It is not generally true that :math:`\log\Gamma(z) =
    \log(\Gamma(z))`, though the real parts of the functions do
    agree. The benefit of not defining `loggamma` as
    :math:`\log(\Gamma(z))` is that the latter function has a
    complicated branch cut structure whereas `loggamma` is analytic
    except for on the negative real axis.

    The identities

    .. math::
      \exp(\log\Gamma(z)) &= \Gamma(z) \\
      \log\Gamma(z + 1) &= \log(z) + \log\Gamma(z)

    make `loggamma` useful for working in complex logspace.

    On the real line `loggamma` is related to `gammaln` via
    ``exp(loggamma(x + 0j)) = gammasgn(x)*exp(gammaln(x))``, up to
    rounding error.

    The implementation here is based on [hare1997]_.

    References
    ----------
    .. [hare1997] D.E.G. Hare,
      *Computing the Principal Branch of log-Gamma*,
      Journal of Algorithms, Volume 25, Issue 2, November 1997, pages 221-236.
    )";

const char *logit_doc = R"(
    """
    logit(x, out=None)

    Logit ufunc for ndarrays.

    The logit function is defined as logit(p) = log(p/(1-p)).
    Note that logit(0) = -inf, logit(1) = inf, and logit(p)
    for p<0 or p>1 yields nan.

    Parameters
    ----------
    x : ndarray
        The ndarray to apply logit to element-wise.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        An ndarray of the same shape as x. Its entries
        are logit of the corresponding entry of x.

    See Also
    --------
    expit

    Notes
    -----
    As a ufunc logit takes a number of optional
    keyword arguments. For more information
    see `ufuncs <https://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_

    .. versionadded:: 0.10.0

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import logit, expit

    >>> logit([0, 0.25, 0.5, 0.75, 1])
    array([       -inf, -1.09861229,  0.        ,  1.09861229,         inf])

    `expit` is the inverse of `logit`:

    >>> expit(logit([0.1, 0.75, 0.999]))
    array([ 0.1  ,  0.75 ,  0.999])

    Plot logit(x) for x in [0, 1]:

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(0, 1, 501)
    >>> y = logit(x)
    >>> plt.plot(x, y)
    >>> plt.grid()
    >>> plt.ylim(-6, 6)
    >>> plt.xlabel('x')
    >>> plt.title('logit(x)')
    >>> plt.show()
    )";

const char *log_expit_doc = R"(
    log_expit(x, out=None)

    Logarithm of the logistic sigmoid function.

    The SciPy implementation of the logistic sigmoid function is
    `scipy.special.expit`, so this function is called ``log_expit``.

    The function is mathematically equivalent to ``log(expit(x))``, but
    is formulated to avoid loss of precision for inputs with large
    (positive or negative) magnitude.

    Parameters
    ----------
    x : array_like
        The values to apply ``log_expit`` to element-wise.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    out : scalar or ndarray
        The computed values, an ndarray of the same shape as ``x``.

    See Also
    --------
    expit

    Notes
    -----
    As a ufunc, ``log_expit`` takes a number of optional keyword arguments.
    For more information see
    `ufuncs <https://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_

    .. versionadded:: 1.8.0

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import log_expit, expit

    >>> log_expit([-3.0, 0.25, 2.5, 5.0])
    array([-3.04858735, -0.57593942, -0.07888973, -0.00671535])

    Large negative values:

    >>> log_expit([-100, -500, -1000])
    array([ -100.,  -500., -1000.])

    Note that ``expit(-1000)`` returns 0, so the naive implementation
    ``log(expit(-1000))`` return ``-inf``.

    Large positive values:

    >>> log_expit([29, 120, 400])
    array([-2.54366565e-013, -7.66764807e-053, -1.91516960e-174])

    Compare that to the naive implementation:

    >>> np.log(expit([29, 120, 400]))
    array([-2.54463117e-13,  0.00000000e+00,  0.00000000e+00])

    The first value is accurate to only 3 digits, and the larger inputs
    lose all precision and return 0.
    )";

const char *log_wright_bessel_doc = R"(
    log_wright_bessel(a, b, x, out=None)

    Natural logarithm of Wright's generalized Bessel function, see `wright_bessel`.
    This function comes in handy in particular for large values of x.

    Parameters
    ----------
    a : array_like of float
        a >= 0
    b : array_like of float
        b >= 0
    x : array_like of float
        x >= 0
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Value of the logarithm of Wright's generalized Bessel function

    Notes
    -----
    Due to the complexity of the function with its three parameters, only
    non-negative arguments are implemented.

    .. versionadded:: 1.14.0

    Examples
    --------
    >>> from scipy.special import log_wright_bessel
    >>> a, b, x = 1.5, 1.1, 2.5
    >>> log_wright_bessel(a, b, x)
    1.1947654935299217
    )";

const char *mathieu_a_doc = R"(
    mathieu_a(m, q, out=None)

    Characteristic value of even Mathieu functions

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Characteristic value for the even solution, ``ce_m(z, q)``, of
        Mathieu's equation.

    See Also
    --------
    mathieu_b, mathieu_cem, mathieu_sem

    )";

const char *mathieu_b_doc = R"(
    mathieu_b(m, q, out=None)

    Characteristic value of odd Mathieu functions

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Characteristic value for the odd solution, ``se_m(z, q)``, of Mathieu's
        equation.

    See Also
    --------
    mathieu_a, mathieu_cem, mathieu_sem

    )";

const char *mathieu_cem_doc = R"(
    mathieu_cem(m, q, x, out=None)

    Even Mathieu function and its derivative

    Returns the even Mathieu function, ``ce_m(x, q)``, of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` of ce_m(x, q)

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    mathieu_a, mathieu_b, mathieu_sem

    )";

const char *mathieu_modcem1_doc = R"(
    mathieu_modcem1(m, q, x, out=None)

    Even modified Mathieu function of the first kind and its derivative

    Evaluates the even modified Mathieu function of the first kind,
    ``Mc1m(x, q)``, and its derivative at `x` for order `m` and parameter
    `q`.

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    mathieu_modsem1

    )";

const char *mathieu_modcem2_doc = R"(
    mathieu_modcem2(m, q, x, out=None)

    Even modified Mathieu function of the second kind and its derivative

    Evaluates the even modified Mathieu function of the second kind,
    Mc2m(x, q), and its derivative at `x` (given in degrees) for order `m`
    and parameter `q`.

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    mathieu_modsem2

    )";

const char *mathieu_modsem1_doc = R"(
    mathieu_modsem1(m, q, x, out=None)

    Odd modified Mathieu function of the first kind and its derivative

    Evaluates the odd modified Mathieu function of the first kind,
    Ms1m(x, q), and its derivative at `x` (given in degrees) for order `m`
    and parameter `q`.

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    mathieu_modcem1

    )";

const char *mathieu_modsem2_doc = R"(
    mathieu_modsem2(m, q, x, out=None)

    Odd modified Mathieu function of the second kind and its derivative

    Evaluates the odd modified Mathieu function of the second kind,
    Ms2m(x, q), and its derivative at `x` (given in degrees) for order `m`
    and parameter q.

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    mathieu_modcem2

    )";

const char *mathieu_sem_doc = R"(
    mathieu_sem(m, q, x, out=None)

    Odd Mathieu function and its derivative

    Returns the odd Mathieu function, se_m(x, q), of order `m` and
    parameter `q` evaluated at `x` (given in degrees).  Also returns the
    derivative with respect to `x` of se_m(x, q).

    Parameters
    ----------
    m : array_like
        Order of the function
    q : array_like
        Parameter of the function
    x : array_like
        Argument of the function, *given in degrees, not radians*.
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    y : scalar or ndarray
        Value of the function
    yp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    mathieu_a, mathieu_b, mathieu_cem

    )";

const char *modfresnelm_doc = R"(
    modfresnelm(x, out=None)

    Modified Fresnel negative integrals

    Parameters
    ----------
    x : array_like
        Function argument
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    fm : scalar or ndarray
        Integral ``F_-(x)``: ``integral(exp(-1j*t*t), t=x..inf)``
    km : scalar or ndarray
        Integral ``K_-(x)``: ``1/sqrt(pi)*exp(1j*(x*x+pi/4))*fp``

    See Also
    --------
    modfresnelp

    )";

const char *modfresnelp_doc = R"(
    modfresnelp(x, out=None)

    Modified Fresnel positive integrals

    Parameters
    ----------
    x : array_like
        Function argument
    out : tuple of ndarray, optional
        Optional output arrays for the function results

    Returns
    -------
    fp : scalar or ndarray
        Integral ``F_+(x)``: ``integral(exp(1j*t*t), t=x..inf)``
    kp : scalar or ndarray
        Integral ``K_+(x)``: ``1/sqrt(pi)*exp(-1j*(x*x+pi/4))*fp``

    See Also
    --------
    modfresnelm

    )";

const char *obl_ang1_doc = R"(
    obl_ang1(m, n, c, x, out=None)

    Oblate spheroidal angular function of the first kind and its derivative

    Computes the oblate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Parameter x (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    obl_ang1_cv

    )";

const char *obl_ang1_cv_doc = R"(
    obl_ang1_cv(m, n, c, cv, x, out=None)

    Oblate spheroidal angular function obl_ang1 for precomputed characteristic value

    Computes the oblate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    cv : array_like
        Characteristic value
    x : array_like
        Parameter x (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    obl_ang1

    )";

const char *obl_cv_doc = R"(
    obl_cv(m, n, c, out=None)

    Characteristic value of oblate spheroidal function

    Computes the characteristic value of oblate spheroidal wave
    functions of order `m`, `n` (n>=m) and spheroidal parameter `c`.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    cv : scalar or ndarray
        Characteristic value

    )";

const char *obl_rad1_doc = R"(
    obl_rad1(m, n, c, x, out=None)

    Oblate spheroidal radial function of the first kind and its derivative

    Computes the oblate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Parameter x (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    obl_rad1_cv

    )";

const char *obl_rad1_cv_doc = R"(
    obl_rad1_cv(m, n, c, cv, x, out=None)

    Oblate spheroidal radial function obl_rad1 for precomputed characteristic value

    Computes the oblate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    cv : array_like
        Characteristic value
    x : array_like
        Parameter x (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    obl_rad1

    )";

const char *obl_rad2_doc = R"(
    obl_rad2(m, n, c, x, out=None)

    Oblate spheroidal radial function of the second kind and its derivative.

    Computes the oblate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Parameter x (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    obl_rad2_cv

    )";

const char *obl_rad2_cv_doc = R"(
    obl_rad2_cv(m, n, c, cv, x, out=None)

    Oblate spheroidal radial function obl_rad2 for precomputed characteristic value

    Computes the oblate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    cv : array_like
        Characteristic value
    x : array_like
        Parameter x (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x

    See Also
    --------
    obl_rad2
    )";

const char *pbdv_doc = R"(
    pbdv(v, x, out=None)

    Parabolic cylinder function D

    Returns (d, dp) the parabolic cylinder function Dv(x) in d and the
    derivative, Dv'(x) in dp.

    Parameters
    ----------
    v : array_like
        Real parameter
    x : array_like
        Real argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    d : scalar or ndarray
        Value of the function
    dp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *pbvv_doc = R"(
    pbvv(v, x, out=None)

    Parabolic cylinder function V

    Returns the parabolic cylinder function Vv(x) in v and the
    derivative, Vv'(x) in vp.

    Parameters
    ----------
    v : array_like
        Real parameter
    x : array_like
        Real argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    v : scalar or ndarray
        Value of the function
    vp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *pbwa_doc = R"(
    pbwa(a, x, out=None)

    Parabolic cylinder function W.

    The function is a particular solution to the differential equation

    .. math::

        y'' + \left(\frac{1}{4}x^2 - a\right)y = 0,

    for a full definition see section 12.14 in [1]_.

    Parameters
    ----------
    a : array_like
        Real parameter
    x : array_like
        Real argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    w : scalar or ndarray
        Value of the function
    wp : scalar or ndarray
        Value of the derivative in x

    Notes
    -----
    The function is a wrapper for a Fortran routine by Zhang and Jin
    [2]_. The implementation is accurate only for ``|a|, |x| < 5`` and
    returns NaN outside that range.

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 14.30.
           https://dlmf.nist.gov/14.30
    .. [2] Zhang, Shanjie and Jin, Jianming. "Computation of Special
           Functions", John Wiley and Sons, 1996.
           https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
    )";

const char *pro_ang1_doc = R"(
    pro_ang1(m, n, c, x, out=None)

    Prolate spheroidal angular function of the first kind and its derivative

    Computes the prolate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Real parameter (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *pro_ang1_cv_doc = R"(
    pro_ang1_cv(m, n, c, cv, x, out=None)

    Prolate spheroidal angular function pro_ang1 for precomputed characteristic value

    Computes the prolate spheroidal angular function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    cv : array_like
        Characteristic value
    x : array_like
        Real parameter (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *pro_cv_doc = R"(
    pro_cv(m, n, c, out=None)

    Characteristic value of prolate spheroidal function

    Computes the characteristic value of prolate spheroidal wave
    functions of order `m`, `n` (n>=m) and spheroidal parameter `c`.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    cv : scalar or ndarray
        Characteristic value
    )";

const char *pro_rad1_doc = R"(
    pro_rad1(m, n, c, x, out=None)

    Prolate spheroidal radial function of the first kind and its derivative

    Computes the prolate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Real parameter (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *pro_rad1_cv_doc = R"(
    pro_rad1_cv(m, n, c, cv, x, out=None)

    Prolate spheroidal radial function pro_rad1 for precomputed characteristic value

    Computes the prolate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    cv : array_like
        Characteristic value
    x : array_like
        Real parameter (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *pro_rad2_doc = R"(
    pro_rad2(m, n, c, x, out=None)

    Prolate spheroidal radial function of the second kind and its derivative

    Computes the prolate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    cv : array_like
        Characteristic value
    x : array_like
        Real parameter (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *pro_rad2_cv_doc = R"(
    pro_rad2_cv(m, n, c, cv, x, out=None)

    Prolate spheroidal radial function pro_rad2 for precomputed characteristic value

    Computes the prolate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``|x| < 1.0``. Requires
    pre-computed characteristic value.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    cv : array_like
        Characteristic value
    x : array_like
        Real parameter (``|x| < 1.0``)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Value of the function
    sp : scalar or ndarray
        Value of the derivative vs x
    )";

const char *psi_doc = R"(
    psi(z, out=None)

    The digamma function.

    The logarithmic derivative of the gamma function evaluated at ``z``.

    Parameters
    ----------
    z : array_like
        Real or complex argument.
    out : ndarray, optional
        Array for the computed values of ``psi``.

    Returns
    -------
    digamma : scalar or ndarray
        Computed values of ``psi``.

    Notes
    -----
    For large values not close to the negative real axis, ``psi`` is
    computed using the asymptotic series (5.11.2) from [1]_. For small
    arguments not close to the negative real axis, the recurrence
    relation (5.5.2) from [1]_ is used until the argument is large
    enough to use the asymptotic series. For values close to the
    negative real axis, the reflection formula (5.5.4) from [1]_ is
    used first. Note that ``psi`` has a family of zeros on the
    negative real axis which occur between the poles at nonpositive
    integers. Around the zeros the reflection formula suffers from
    cancellation and the implementation loses precision. The sole
    positive zero and the first negative zero, however, are handled
    separately by precomputing series expansions using [2]_, so the
    function should maintain full accuracy around the origin.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/5
    .. [2] Fredrik Johansson and others.
           "mpmath: a Python library for arbitrary-precision floating-point arithmetic"
           (Version 0.19) http://mpmath.org/

    Examples
    --------
    >>> from scipy.special import psi
    >>> z = 3 + 4j
    >>> psi(z)
    (1.55035981733341+1.0105022091860445j)

    Verify psi(z) = psi(z + 1) - 1/z:

    >>> psi(z + 1) - 1/z
    (1.55035981733341+1.0105022091860445j)
    )";

const char *rgamma_doc = R"(
    rgamma(z, out=None)

    Reciprocal of the gamma function.

    Defined as :math:`1 / \Gamma(z)`, where :math:`\Gamma` is the
    gamma function. For more on the gamma function see `gamma`.

    Parameters
    ----------
    z : array_like
        Real or complex valued input
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Function results

    See Also
    --------
    gamma, gammaln, loggamma

    Notes
    -----
    The gamma function has no zeros and has simple poles at
    nonpositive integers, so `rgamma` is an entire function with zeros
    at the nonpositive integers. See the discussion in [dlmf]_ for
    more details.

    References
    ----------
    .. [dlmf] Nist, Digital Library of Mathematical functions,
        https://dlmf.nist.gov/5.2#i

    Examples
    --------
    >>> import scipy.special as sc

    It is the reciprocal of the gamma function.

    >>> sc.rgamma([1, 2, 3, 4])
    array([1.        , 1.        , 0.5       , 0.16666667])
    >>> 1 / sc.gamma([1, 2, 3, 4])
    array([1.        , 1.        , 0.5       , 0.16666667])

    It is zero at nonpositive integers.

    >>> sc.rgamma([0, -1, -2, -3])
    array([0., 0., 0., 0.])

    It rapidly underflows to zero along the positive real axis.

    >>> sc.rgamma([10, 100, 179])
    array([2.75573192e-006, 1.07151029e-156, 0.00000000e+000])
    )";

const char *scaled_exp1_doc = R"(
    _scaled_exp1(x, out=None):

    Compute the scaled exponential integral.

    This is a private function, subject to change or removal with no
    deprecation.

    This function computes F(x), where F is the factor remaining in E_1(x)
    when exp(-x)/x is factored out.  That is,::

        E_1(x) = exp(-x)/x * F(x)

    or

        F(x) = x * exp(x) * E_1(x)

    The function is defined for real x >= 0.  For x < 0, nan is returned.

    F has the properties:

    * F(0) = 0
    * F(x) is increasing on [0, inf).
    * The limit as x goes to infinity of F(x) is 1.

    Parameters
    ----------
    x: array_like
        The input values. Must be real.  The implementation is limited to
        double precision floating point, so other types will be cast to
        to double precision.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the scaled exponential integral.

    See Also
    --------
    exp1 : exponential integral E_1

    Examples
    --------
    >>> from scipy.special import _scaled_exp1
    >>> _scaled_exp1([0, 0.1, 1, 10, 100])
    )";

const char *spherical_jn_doc = R"(
    Internal function, use `spherical_jn` instead.
    )";

const char *spherical_jn_d_doc = R"(
    Internal function, use `spherical_jn` instead.
    )";

const char *spherical_yn_doc = R"(
    Internal function, use `spherical_yn` instead.
    )";

const char *spherical_yn_d_doc = R"(
    Internal function, use `spherical_yn` instead.
    )";

const char *spherical_in_doc = R"(
    Internal function, use `spherical_in` instead.
    )";

const char *spherical_in_d_doc = R"(
    Internal function, use `spherical_in` instead.
    )";

const char *spherical_kn_doc = R"(
    Internal function, use `spherical_kn` instead.
    )";

const char *spherical_kn_d_doc = R"(
    Internal function, use `spherical_kn` instead.
    )";

const char *sph_harm_doc = R"(
    sph_harm(m, n, theta, phi, out=None)

    Compute spherical harmonics.

    The spherical harmonics are defined as

    .. math::

        Y^m_n(\theta,\phi) = \sqrt{\frac{2n+1}{4\pi} \frac{(n-m)!}{(n+m)!}}
          e^{i m \theta} P^m_n(\cos(\phi))

    where :math:`P_n^m` are the associated Legendre functions; see `lpmv`.

    .. deprecated:: 1.15.0
        This function is deprecated and will be removed in a future version.
        Use `scipy.special.sph_harm_y` instead.

    Parameters
    ----------
    m : array_like
        Order of the harmonic (int); must have ``|m| <= n``.
    n : array_like
       Degree of the harmonic (int); must have ``n >= 0``. This is
       often denoted by ``l`` (lower case L) in descriptions of
       spherical harmonics.
    theta : array_like
       Azimuthal (longitudinal) coordinate; must be in ``[0, 2*pi]``.
    phi : array_like
       Polar (colatitudinal) coordinate; must be in ``[0, pi]``.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    y_mn : complex scalar or ndarray
       The harmonic :math:`Y^m_n` sampled at ``theta`` and ``phi``.

    Notes
    -----
    There are different conventions for the meanings of the input
    arguments ``theta`` and ``phi``. In SciPy ``theta`` is the
    azimuthal angle and ``phi`` is the polar angle. It is common to
    see the opposite convention, that is, ``theta`` as the polar angle
    and ``phi`` as the azimuthal angle.

    Note that SciPy's spherical harmonics include the Condon-Shortley
    phase [2]_ because it is part of `lpmv`.

    With SciPy's conventions, the first several spherical harmonics
    are

    .. math::

        Y_0^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{1}{\pi}} \\
        Y_1^{-1}(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                    e^{-i\theta} \sin(\phi) \\
        Y_1^0(\theta, \phi) &= \frac{1}{2} \sqrt{\frac{3}{\pi}}
                                 \cos(\phi) \\
        Y_1^1(\theta, \phi) &= -\frac{1}{2} \sqrt{\frac{3}{2\pi}}
                                 e^{i\theta} \sin(\phi).

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 14.30.
           https://dlmf.nist.gov/14.30
    .. [2] https://en.wikipedia.org/wiki/Spherical_harmonics#Condon.E2.80.93Shortley_phase
    )";

const char *wright_bessel_doc = R"(
    wright_bessel(a, b, x, out=None)

    Wright's generalized Bessel function.

    Wright's generalized Bessel function is an entire function and defined as

    .. math:: \Phi(a, b; x) = \sum_{k=0}^\infty \frac{x^k}{k! \Gamma(a k + b)}

    See Also [1].

    Parameters
    ----------
    a : array_like of float
        a >= 0
    b : array_like of float
        b >= 0
    x : array_like of float
        x >= 0
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Value of the Wright's generalized Bessel function

    Notes
    -----
    Due to the complexity of the function with its three parameters, only
    non-negative arguments are implemented.

    .. versionadded:: 1.7.0

    References
    ----------
    .. [1] Digital Library of Mathematical Functions, 10.46.
           https://dlmf.nist.gov/10.46.E1

    Examples
    --------
    >>> from scipy.special import wright_bessel
    >>> a, b, x = 1.5, 1.1, 2.5
    >>> wright_bessel(a, b-1, x)
    4.5314465939443025

    Now, let us verify the relation

    .. math:: \Phi(a, b-1; x) = a x \Phi(a, b+a; x) + (b-1) \Phi(a, b; x)

    >>> a * x * wright_bessel(a, b+a, x) + (b-1) * wright_bessel(a, b, x)
    4.5314465939443025
    )";

const char *yv_doc = R"(
    yv(v, z, out=None)

    Bessel function of the second kind of real order and complex argument.

    Parameters
    ----------
    v : array_like
        Order (float).
    z : array_like
        Argument (float or complex).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    Y : scalar or ndarray
        Value of the Bessel function of the second kind, :math:`Y_v(x)`.

    See Also
    --------
    yve : :math:`Y_v` with leading exponential behavior stripped off.
    y0: faster implementation of this function for order 0
    y1: faster implementation of this function for order 1

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

    Examples
    --------
    Evaluate the function of order 0 at one point.

    >>> from scipy.special import yv
    >>> yv(0, 1.)
    0.088256964215677

    Evaluate the function at one point for different orders.

    >>> yv(0, 1.), yv(1, 1.), yv(1.5, 1.)
    (0.088256964215677, -0.7812128213002889, -1.102495575160179)

    The evaluation for different orders can be carried out in one call by
    providing a list or NumPy array as argument for the `v` parameter:

    >>> yv([0, 1, 1.5], 1.)
    array([ 0.08825696, -0.78121282, -1.10249558])

    Evaluate the function at several points for order 0 by providing an
    array for `z`.

    >>> import numpy as np
    >>> points = np.array([0.5, 3., 8.])
    >>> yv(0, points)
    array([-0.44451873,  0.37685001,  0.22352149])

    If `z` is an array, the order parameter `v` must be broadcastable to
    the correct shape if different orders shall be computed in one call.
    To calculate the orders 0 and 1 for an 1D array:

    >>> orders = np.array([[0], [1]])
    >>> orders.shape
    (2, 1)

    >>> yv(orders, points)
    array([[-0.44451873,  0.37685001,  0.22352149],
           [-1.47147239,  0.32467442, -0.15806046]])

    Plot the functions of order 0 to 3 from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, yv(i, x), label=f'$Y_{i!r}$')
    >>> ax.set_ylim(-3, 1)
    >>> ax.legend()
    >>> plt.show()

    )";

const char *yve_doc = R"(
    yve(v, z, out=None)

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
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    Y : scalar or ndarray
        Value of the exponentially scaled Bessel function.

    See Also
    --------
    yv: Unscaled Bessel function of the second kind of real order.

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

    Exponentially scaled Bessel functions are useful for large `z`:
    for these, the unscaled Bessel functions can easily under-or overflow.

    References
    ----------
    .. [1] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    Compare the output of `yv` and `yve` for large complex arguments for `z`
    by computing their values for order ``v=1`` at ``z=1000j``. We see that
    `yv` returns nan but `yve` returns a finite number:

    >>> import numpy as np
    >>> from scipy.special import yv, yve
    >>> v = 1
    >>> z = 1000j
    >>> yv(v, z), yve(v, z)
    ((nan+nanj), (-0.012610930256928629+7.721967686709076e-19j))

    For real arguments for `z`, `yve` returns the same as `yv` up to
    floating point errors.

    >>> v, z = 1, 1000
    >>> yv(v, z), yve(v, z)
    (-0.02478433129235178, -0.02478433129235179)

    The function can be evaluated for several orders at the same time by
    providing a list or NumPy array for `v`:

    >>> yve([1, 2, 3], 1j)
    array([-0.20791042+0.14096627j,  0.38053618-0.04993878j,
           0.00815531-1.66311097j])

    In the same way, the function can be evaluated at several points in one
    call by providing a list or NumPy array for `z`:

    >>> yve(1, np.array([1j, 2j, 3j]))
    array([-0.20791042+0.14096627j, -0.21526929+0.01205044j,
           -0.19682671+0.00127278j])

    It is also possible to evaluate several orders at several points
    at the same time by providing arrays for `v` and `z` with
    broadcasting compatible shapes. Compute `yve` for two different orders
    `v` and three points `z` resulting in a 2x3 array.

    >>> v = np.array([[1], [2]])
    >>> z = np.array([3j, 4j, 5j])
    >>> v.shape, z.shape
    ((2, 1), (3,))

    >>> yve(v, z)
    array([[-1.96826713e-01+1.27277544e-03j, -1.78750840e-01+1.45558819e-04j,
            -1.63972267e-01+1.73494110e-05j],
           [1.94960056e-03-1.11782545e-01j,  2.02902325e-04-1.17626501e-01j,
            2.27727687e-05-1.17951906e-01j]])
    )";
