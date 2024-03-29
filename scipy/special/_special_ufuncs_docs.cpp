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