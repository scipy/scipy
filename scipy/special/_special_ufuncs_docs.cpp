const char *_cospi_doc = R"(
    Internal function, do not use.
    )";

const char *_bivariate_normal_sf_doc = R"(
    Internal function, do not use.
    )";

const char *_cosine_cdf_doc = R"(
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
    )";

const char *_cosine_invcdf_doc = R"(
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
    )";

const char *_igam_fac_doc = R"(
    Internal function, do not use.
    )";

const char *_kolmogc_doc = R"(
    Internal function, do not use.
    )";

const char *_kolmogci_doc = R"(
    Internal function, do not use.
    )";

const char *_kolmogp_doc = R"(
    Internal function, do not use.
    )";

const char *_lgam1p_doc = R"(
    Internal function, do not use.
    )";

const char *_stirling2_inexact_doc = R"(
    Internal function, do not use.
    )";

const char *hyp0f1_doc = R"(
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

    )";

const char *agm_doc = R"(
    agm(a, b, out=None)

    Compute the arithmetic-geometric mean of `a` and `b`.

    Start with :math:`a_0 = a` and :math:`b_0 = b` and iteratively compute

    .. math::

        a_{n+1} = \frac{a_n + b_n}{2}, \quad
        b_{n+1} = \sqrt{a_n b_n}.

    :math:`a_n` and :math:`b_n` converge to the same limit as :math:`n`
    increases; their common limit is :math:`\operatorname{agm}(a, b)`.

    Parameters
    ----------
    a, b : array_like
        Real values only. If the values are both negative, the result
        is negative. If one value is negative and the other is positive,
        `nan` is returned.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        The arithmetic-geometric mean of :math:`a` and :math:`b`.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import agm
    >>> a, b = 24.0, 6.0
    >>> agm(a, b)
    13.458171481725614

    Compare that result to the iteration:

    >>> while a != b:
    ...     a, b = (a + b)/2, np.sqrt(a*b)
    ...     print("a = %19.16f  b=%19.16f" % (a, b))
    ...
    a = 15.0000000000000000  b=12.0000000000000000
    a = 13.5000000000000000  b=13.4164078649987388
    a = 13.4582039324993694  b=13.4581390309909850
    a = 13.4581714817451772  b=13.4581714817060547
    a = 13.4581714817256159  b=13.4581714817256159

    When array-like arguments are given, broadcasting applies:

    >>> a = np.array([[1.5], [3], [6]])  # a has shape (3, 1).
    >>> b = np.array([6, 12, 24, 48])    # b has shape (4,).
    >>> agm(a, b)
    array([[  3.36454287,   5.42363427,   9.05798751,  15.53650756],
           [  4.37037309,   6.72908574,  10.84726853,  18.11597502],
           [  6.        ,   8.74074619,  13.45817148,  21.69453707]])
    )";

const char *entr_doc = R"(
    entr(x, out=None)

    Elementwise function for computing entropy.

    .. math:: \text{entr}(x) = \begin{cases} - x \log(x) & x > 0  \\ 0 & x = 0
              \\ -\infty & \text{otherwise} \end{cases}

    Parameters
    ----------
    x : ndarray
        Input array.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    res : scalar or ndarray
        The value of the elementwise entropy function at the given points `x`.

    See Also
    --------
    kl_div, rel_entr, scipy.stats.entropy

    Notes
    -----
    .. versionadded:: 0.15.0

    This function is concave.

    The origin of this function is in convex programming; see [1]_.
    Given a probability distribution :math:`p_1, \ldots, p_n`,
    the definition of entropy in the context of *information theory* is

    .. math::

        \sum_{i = 1}^n \mathrm{entr}(p_i).

    To compute the latter quantity, use `scipy.stats.entropy`.

    References
    ----------
    .. [1] Boyd, Stephen and Lieven Vandenberghe. *Convex optimization*.
           Cambridge University Press, 2004.
           :doi:`10.1017/CBO9780511804441`.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import entr

    Calculate the entropy (in nats) of a 3-outcome probability distribution

    >>> p = np.array([0.2, 0.5, 0.3])
    >>> entr(p)
    array([0.32188758, 0.34657359, 0.36119184])
    >>> entr(p).sum()
    1.0296530140645737

    )";

const char *chdtr_doc = R"(
    chdtr(v, x, out=None)

    Chi square cumulative distribution function.

    Returns the area under the left tail (from 0 to `x`) of the Chi
    square probability density function with `v` degrees of freedom:

    .. math::

        \frac{1}{2^{v/2} \Gamma(v/2)} \int_0^x t^{v/2 - 1} e^{-t/2} dt

    Here :math:`\Gamma` is the Gamma function; see `gamma`. This
    integral can be expressed in terms of the regularized lower
    incomplete gamma function `gammainc` as
    ``gammainc(v / 2, x / 2)``. [1]_

    Parameters
    ----------
    v : array_like
        Degrees of freedom.
    x : array_like
        Upper bound of the integral.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the cumulative distribution function.

    See Also
    --------
    chdtrc, chdtri, chdtriv, gammainc

    References
    ----------
    .. [1] Chi-Square distribution,
        https://www.itl.nist.gov/div898/handbook/eda/section3/eda3666.htm

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It can be expressed in terms of the regularized lower incomplete
    gamma function.

    >>> v = 1
    >>> x = np.arange(4)
    >>> sc.chdtr(v, x)
    array([0.        , 0.68268949, 0.84270079, 0.91673548])
    >>> sc.gammainc(v / 2, x / 2)
    array([0.        , 0.68268949, 0.84270079, 0.91673548])

    )";

const char *chdtrc_doc = R"(
    chdtrc(v, x, out=None)

    Chi square survival function.

    Returns the area under the right hand tail (from `x` to infinity)
    of the Chi square probability density function with `v` degrees of
    freedom:

    .. math::

        \frac{1}{2^{v/2} \Gamma(v/2)} \int_x^\infty t^{v/2 - 1} e^{-t/2} dt

    Here :math:`\Gamma` is the Gamma function; see `gamma`. This
    integral can be expressed in terms of the regularized upper
    incomplete gamma function `gammaincc` as
    ``gammaincc(v / 2, x / 2)``. [1]_

    Parameters
    ----------
    v : array_like
        Degrees of freedom.
    x : array_like
        Lower bound of the integral.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the survival function.

    See Also
    --------
    chdtr, chdtri, chdtriv, gammaincc

    References
    ----------
    .. [1] Chi-Square distribution,
        https://www.itl.nist.gov/div898/handbook/eda/section3/eda3666.htm

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It can be expressed in terms of the regularized upper incomplete
    gamma function.

    >>> v = 1
    >>> x = np.arange(4)
    >>> sc.chdtrc(v, x)
    array([1.        , 0.31731051, 0.15729921, 0.08326452])
    >>> sc.gammaincc(v / 2, x / 2)
    array([1.        , 0.31731051, 0.15729921, 0.08326452])

    )";

const char *chdtri_doc = R"(
    chdtri(v, p, out=None)

    Inverse to `chdtrc` with respect to `x`.

    Returns `x` such that ``chdtrc(v, x) == p``.

    Parameters
    ----------
    v : array_like
        Degrees of freedom.
    p : array_like
        Probability.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    x : scalar or ndarray
        Value so that the probability a Chi square random variable
        with `v` degrees of freedom is greater than `x` equals `p`.

    See Also
    --------
    chdtrc, chdtr, chdtriv

    References
    ----------
    .. [1] Chi-Square distribution,
        https://www.itl.nist.gov/div898/handbook/eda/section3/eda3666.htm

    Examples
    --------
    >>> import scipy.special as sc

    It inverts `chdtrc`.

    >>> v, p = 1, 0.3
    >>> sc.chdtrc(v, sc.chdtri(v, p))
    0.3
    >>> x = 1
    >>> sc.chdtri(v, sc.chdtrc(v, x))
    1.0

    )";

const char *erfcinv_doc = R"(
    erfcinv(y, out=None)

    Inverse of the complementary error function.

    Computes the inverse of the complementary error function.

    In the complex domain, there is no unique complex number :math:`w` satisfying
    :math:`\\operatorname{erfc}(w) = z`. This indicates a true inverse function
    would be multivalued.
    When the domain restricts to the real interval :math:`0 < x < 2`, there is
    a unique real number satisfying

    .. math::

        \\operatorname{erfc}(\\operatorname{erfcinv}(x)) = x

    It is related to the inverse of the error function by

    .. math::

        \\operatorname{erfcinv}(1 - x) = \\operatorname{erfinv}(x)

    Parameters
    ----------
    y : ndarray
        Argument at which to evaluate. Domain: :math:`[0, 2]`
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    erfcinv : scalar or ndarray
        The inverse of :math:`\\operatorname{erfc}` of :math:`y`, element-wise

    See Also
    --------
    erf : Error function
    erfc : Complementary error function
    erfinv : Inverse of the error function

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.special import erfcinv

    >>> erfcinv(0.5)
    0.4769362762044699

    >>> y = np.linspace(0.0, 2.0, num=11)
    >>> erfcinv(y)
    array([        inf,  0.9061938 ,  0.59511608,  0.37080716,  0.17914345,
           -0.        , -0.17914345, -0.37080716, -0.59511608, -0.9061938 ,
                  -inf])

    Plot the function:

    >>> y = np.linspace(0, 2, 200)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(y, erfcinv(y))
    >>> ax.grid(True)
    >>> ax.set_xlabel('y')
    >>> ax.set_title('erfcinv(y)')
    >>> plt.show()

    )";

const char *lpmv_doc = R"(
    lpmv(m, v, x, out=None)

    Associated Legendre function of integer order and real degree.

    Defined as

    .. math::

        P_v^m(x) = (-1)^m (1 - x^2)^{m/2} \frac{d^m}{dx^m} (P_v(x))

    where

    .. math::

        P_v(x) = \sum_{k = 0}^\infty \frac{(-v)_k (v + 1)_k}{(k!)^2}
                \left(\frac{1 - x}{2}\right)^k

    is the Legendre function of the first kind. Here :math:`(\cdot)_k`
    is the Pochhammer symbol; see `poch`.

    Parameters
    ----------
    m : array_like
        Order (int or float). If passed a float not equal to an
        integer the function returns NaN.
    v : array_like
        Degree (float).
    x : array_like
        Argument (float). Must have ``|x| <= 1``.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    pmv : scalar or ndarray
        Value of the associated Legendre function.

    Notes
    -----
    Note that this implementation includes the Condon-Shortley phase.

    References
    ----------
    .. [1] S. Zhang and J.M. Jin, "Computation of Special Functions", Wiley 1996.

    Examples
    --------
    Compute a single associated Legendre function value:

    >>> from scipy.special import lpmv
    >>> lpmv(1, 2.5, 0.3)
    -0.1292299

    Plot the associated Legendre functions for orders ``m = 1`` and
    ``m = 2`` and degrees ``v = 0, ..., 5``:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> nus = np.arange(0, 6)
    >>> xs = np.linspace(-1+1e-5, 1-1e-5, 101)

    Order 1.

    >>> fig, ax = plt.subplots()
    >>> ax.plot(xs, lpmv(1, nus, xs[:, np.newaxis]), "-")
    >>> ax.legend([rf"$\nu={nu}$" for nu in nus])
    >>> ax.set_xlabel(r"$x$")
    >>> ax.set_ylabel(r"$P_\nu^1(x)$")
    >>> plt.show()

    Order 2.

    >>> fig, ax = plt.subplots()
    >>> ax.plot(xs, lpmv(2, nus, xs[:, np.newaxis]), "-")
    >>> ax.legend([rf"$\nu={nu}$" for nu in nus])
    >>> ax.set_xlabel(r"$x$")
    >>> ax.set_ylabel(r"$P_\nu^2(x)$")
    >>> plt.show()
    )";

const char *huber_doc = R"(
    huber(delta, r, out=None)

    Huber loss function.

    .. math:: \text{huber}(\delta, r) = \begin{cases} \infty & \delta < 0  \\
              \frac{1}{2}r^2 & 0 \le \delta, | r | \le \delta \\
              \delta ( |r| - \frac{1}{2}\delta ) & \text{otherwise} \end{cases}

    Parameters
    ----------
    delta : ndarray
        Input array, indicating the quadratic vs. linear loss changepoint.
    r : ndarray
        Input array, possibly representing residuals.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        The computed Huber loss function values.

    See Also
    --------
    pseudo_huber : smooth approximation of this function

    Notes
    -----
    `huber` is useful as a loss function in robust statistics or machine
    learning to reduce the influence of outliers as compared to the common
    squared error loss, residuals with a magnitude higher than `delta` are
    not squared [1]_.

    Typically, `r` represents residuals, the difference
    between a model prediction and data. Then, for :math:`|r|\leq\delta`,
    `huber` resembles the squared error and for :math:`|r|>\delta` the
    absolute error. This way, the Huber loss often achieves
    a fast convergence in model fitting for small residuals like the squared
    error loss function and still reduces the influence of outliers
    (:math:`|r|>\delta`) like the absolute error loss. As :math:`\delta` is
    the cutoff between squared and absolute error regimes, it has
    to be tuned carefully for each problem. `huber` is also
    convex, making it suitable for gradient based optimization.

    .. versionadded:: 0.15.0

    References
    ----------
    .. [1] Peter Huber. "Robust Estimation of a Location Parameter",
           1964. Annals of Statistics. 53 (1): 73 - 101.

    Examples
    --------
    Import all necessary modules.

    >>> import numpy as np
    >>> from scipy.special import huber
    >>> import matplotlib.pyplot as plt

    Compute the function for ``delta=1`` at ``r=2``

    >>> huber(1., 2.)
    1.5

    Compute the function for different `delta` by providing a NumPy array or
    list for `delta`.

    >>> huber([1., 3., 5.], 4.)
    array([3.5, 7.5, 8. ])

    Compute the function at different points by providing a NumPy array or
    list for `r`.

    >>> huber(2., np.array([1., 1.5, 3.]))
    array([0.5  , 1.125, 4.   ])

    The function can be calculated for different `delta` and `r` by
    providing arrays for both with compatible shapes for broadcasting.

    >>> r = np.array([1., 2.5, 8., 10.])
    >>> deltas = np.array([[1.], [5.], [9.]])
    >>> print(r.shape, deltas.shape)
    (4,) (3, 1)

    >>> huber(deltas, r)
    array([[ 0.5  ,  2.   ,  7.5  ,  9.5  ],
           [ 0.5  ,  3.125, 27.5  , 37.5  ],
           [ 0.5  ,  3.125, 32.   , 49.5  ]])

    Plot the function for different `delta`.

    >>> x = np.linspace(-4, 4, 500)
    >>> deltas = [1, 2, 3]
    >>> linestyles = ["dashed", "dotted", "dashdot"]
    >>> fig, ax = plt.subplots()
    >>> combined_plot_parameters = list(zip(deltas, linestyles))
    >>> for delta, style in combined_plot_parameters:
    ...     ax.plot(x, huber(delta, x), label=fr"$\delta={delta}$", ls=style)
    >>> ax.legend(loc="upper center")
    >>> ax.set_xlabel("$x$")
    >>> ax.set_title(r"Huber loss function $h_{\delta}(x)$")
    >>> ax.set_xlim(-4, 4)
    >>> ax.set_ylim(0, 8)
    >>> plt.show()
    )";

const char *kl_div_doc = R"(
    kl_div(x, y, out=None)

    Elementwise function for computing Kullback-Leibler divergence.

    .. math::

        \mathrm{kl\_div}(x, y) =
          \begin{cases}
            x \log(x / y) - x + y & x > 0, y > 0 \\
            y & x = 0, y \ge 0 \\
            \infty & \text{otherwise}
          \end{cases}

    Parameters
    ----------
    x, y : array_like
        Real arguments
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the Kullback-Liebler divergence.

    See Also
    --------
    entr, rel_entr, scipy.stats.entropy

    Notes
    -----
    .. versionadded:: 0.15.0

    This function is non-negative and is jointly convex in `x` and `y`.

    The origin of this function is in convex programming; see [1]_ for
    details. This is why the function contains the extra :math:`-x
    + y` terms over what might be expected from the Kullback-Leibler
    divergence. For a version of the function without the extra terms,
    see `rel_entr`.

    References
    ----------
    .. [1] Boyd, Stephen and Lieven Vandenberghe. *Convex optimization*.
           Cambridge University Press, 2004.
           :doi:`10.1017/CBO9780511804441`.

    )";

const char *pseudo_huber_doc = R"(
    pseudo_huber(delta, r, out=None)

    Pseudo-Huber loss function.

    .. math:: \mathrm{pseudo\_huber}(\delta, r) =
              \delta^2 \left( \sqrt{ 1 + \left( \frac{r}{\delta} \right)^2 } - 1 \right)

    Parameters
    ----------
    delta : array_like
        Input array, indicating the soft quadratic vs. linear loss changepoint.
    r : array_like
        Input array, possibly representing residuals.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    res : scalar or ndarray
        The computed Pseudo-Huber loss function values.

    See Also
    --------
    huber: Similar function which this function approximates

    Notes
    -----
    Like `huber`, `pseudo_huber` often serves as a robust loss function
    in statistics or machine learning to reduce the influence of outliers.
    Unlike `huber`, `pseudo_huber` is smooth.

    Typically, `r` represents residuals, the difference
    between a model prediction and data. Then, for :math:`|r|\leq\delta`,
    `pseudo_huber` resembles the squared error and for :math:`|r|>\delta` the
    absolute error. This way, the Pseudo-Huber loss often achieves
    a fast convergence in model fitting for small residuals like the squared
    error loss function and still reduces the influence of outliers
    (:math:`|r|>\delta`) like the absolute error loss. As :math:`\delta` is
    the cutoff between squared and absolute error regimes, it has
    to be tuned carefully for each problem. `pseudo_huber` is also
    convex, making it suitable for gradient based optimization. [1]_ [2]_

    .. versionadded:: 0.15.0

    References
    ----------
    .. [1] Hartley, Zisserman, "Multiple View Geometry in Computer Vision".
           2003. Cambridge University Press. p. 619
    .. [2] Charbonnier et al. "Deterministic edge-preserving regularization
           in computed imaging". 1997. IEEE Trans. Image Processing.
           6 (2): 298 - 311.

    Examples
    --------
    Import all necessary modules.

    >>> import numpy as np
    >>> from scipy.special import pseudo_huber, huber
    >>> import matplotlib.pyplot as plt

    Calculate the function for ``delta=1`` at ``r=2``.

    >>> pseudo_huber(1., 2.)
    1.2360679774997898

    Calculate the function at ``r=2`` for different `delta` by providing
    a list or NumPy array for `delta`.

    >>> pseudo_huber([1., 2., 4.], 3.)
    array([2.16227766, 3.21110255, 4.        ])

    Calculate the function for ``delta=1`` at several points by providing
    a list or NumPy array for `r`.

    >>> pseudo_huber(2., np.array([1., 1.5, 3., 4.]))
    array([0.47213595, 1.        , 3.21110255, 4.94427191])

    The function can be calculated for different `delta` and `r` by
    providing arrays for both with compatible shapes for broadcasting.

    >>> r = np.array([1., 2.5, 8., 10.])
    >>> deltas = np.array([[1.], [5.], [9.]])
    >>> print(r.shape, deltas.shape)
    (4,) (3, 1)

    >>> pseudo_huber(deltas, r)
    array([[ 0.41421356,  1.6925824 ,  7.06225775,  9.04987562],
           [ 0.49509757,  2.95084972, 22.16990566, 30.90169944],
           [ 0.49846624,  3.06693762, 27.37435121, 40.08261642]])

    Plot the function for different `delta`.

    >>> x = np.linspace(-4, 4, 500)
    >>> deltas = [1, 2, 3]
    >>> linestyles = ["dashed", "dotted", "dashdot"]
    >>> fig, ax = plt.subplots()
    >>> combined_plot_parameters = list(zip(deltas, linestyles))
    >>> for delta, style in combined_plot_parameters:
    ...     ax.plot(x, pseudo_huber(delta, x), label=rf"$\delta={delta}$",
    ...             ls=style)
    >>> ax.legend(loc="upper center")
    >>> ax.set_xlabel("$x$")
    >>> ax.set_title(r"Pseudo-Huber loss function $h_{\delta}(x)$")
    >>> ax.set_xlim(-4, 4)
    >>> ax.set_ylim(0, 8)
    >>> plt.show()

    Finally, illustrate the difference between `huber` and `pseudo_huber` by
    plotting them and their gradients with respect to `r`. The plot shows
    that `pseudo_huber` is continuously differentiable while `huber` is not
    at the points :math:`\pm\delta`.

    >>> def huber_grad(delta, x):
    ...     grad = np.copy(x)
    ...     linear_area = np.argwhere(np.abs(x) > delta)
    ...     grad[linear_area]=delta*np.sign(x[linear_area])
    ...     return grad
    >>> def pseudo_huber_grad(delta, x):
    ...     return x* (1+(x/delta)**2)**(-0.5)
    >>> x=np.linspace(-3, 3, 500)
    >>> delta = 1.
    >>> fig, ax = plt.subplots(figsize=(7, 7))
    >>> ax.plot(x, huber(delta, x), label="Huber", ls="dashed")
    >>> ax.plot(x, huber_grad(delta, x), label="Huber Gradient", ls="dashdot")
    >>> ax.plot(x, pseudo_huber(delta, x), label="Pseudo-Huber", ls="dotted")
    >>> ax.plot(x, pseudo_huber_grad(delta, x), label="Pseudo-Huber Gradient",
    ...         ls="solid")
    >>> ax.legend(loc="upper center")
    >>> plt.show()
    )";

const char *rel_entr_doc = R"(
    rel_entr(x, y, out=None)

    Elementwise function for computing relative entropy.

    .. math::

        \mathrm{rel\_entr}(x, y) =
            \begin{cases}
                x \log(x / y) & x > 0, y > 0 \\
                0 & x = 0, y \ge 0 \\
                \infty & \text{otherwise}
            \end{cases}

    Parameters
    ----------
    x, y : array_like
        Input arrays
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Relative entropy of the inputs

    See Also
    --------
    entr, kl_div, scipy.stats.entropy

    Notes
    -----
    .. versionadded:: 0.15.0

    This function is jointly convex in x and y.

    The origin of this function is in convex programming; see
    [1]_. Given two discrete probability distributions :math:`p_1,
    \ldots, p_n` and :math:`q_1, \ldots, q_n`, the definition of relative
    entropy in the context of *information theory* is

    .. math::

        \sum_{i = 1}^n \mathrm{rel\_entr}(p_i, q_i).

    To compute the latter quantity, use `scipy.stats.entropy`.

    See [2]_ for details.

    References
    ----------
    .. [1] Boyd, Stephen and Lieven Vandenberghe. *Convex optimization*.
           Cambridge University Press, 2004.
           :doi:`10.1017/CBO9780511804441`.
    .. [2] Kullback-Leibler divergence,
           https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence

    )";

const char *gdtr_doc = R"(
    gdtr(a, b, x, out=None)

    Gamma distribution cumulative distribution function.

    Returns the integral from zero to `x` of the gamma probability density
    function,

    .. math::

        F(x) = \int_0^x \frac{a^b}{\Gamma(b)} t^{b-1} e^{-at}\,dt,

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        The CDF of the gamma distribution with parameters `a` and `b`
        evaluated at `x`.

    See Also
    --------
    gdtrc : 1 - CDF of the gamma distribution.
    scipy.stats.gamma: Gamma distribution

    Notes
    -----
    The evaluation is carried out using the relation to the incomplete gamma
    integral (regularized gamma function).

    Wrapper for the Cephes [1]_ routine `gdtr`. Calling `gdtr` directly can
    improve performance compared to the ``cdf`` method of `scipy.stats.gamma`
    (see last example below).

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    Examples
    --------
    Compute the function for ``a=1``, ``b=2`` at ``x=5``.

    >>> import numpy as np
    >>> from scipy.special import gdtr
    >>> import matplotlib.pyplot as plt
    >>> gdtr(1., 2., 5.)
    0.9595723180054873

    Compute the function for ``a=1`` and ``b=2`` at several points by
    providing a NumPy array for `x`.

    >>> xvalues = np.array([1., 2., 3., 4])
    >>> gdtr(1., 1., xvalues)
    array([0.63212056, 0.86466472, 0.95021293, 0.98168436])

    `gdtr` can evaluate different parameter sets by providing arrays with
    broadcasting compatible shapes for `a`, `b` and `x`. Here we compute the
    function for three different `a` at four positions `x` and ``b=3``,
    resulting in a 3x4 array.

    >>> a = np.array([[0.5], [1.5], [2.5]])
    >>> x = np.array([1., 2., 3., 4])
    >>> a.shape, x.shape
    ((3, 1), (4,))

    >>> gdtr(a, 3., x)
    array([[0.01438768, 0.0803014 , 0.19115317, 0.32332358],
           [0.19115317, 0.57680992, 0.82642193, 0.9380312 ],
           [0.45618688, 0.87534798, 0.97974328, 0.9972306 ]])

    Plot the function for four different parameter sets.

    >>> a_parameters = [0.3, 1, 2, 6]
    >>> b_parameters = [2, 10, 15, 20]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(a_parameters, b_parameters, linestyles))
    >>> x = np.linspace(0, 30, 1000)
    >>> fig, ax = plt.subplots()
    >>> for parameter_set in parameters_list:
    ...     a, b, style = parameter_set
    ...     gdtr_vals = gdtr(a, b, x)
    ...     ax.plot(x, gdtr_vals, label=fr"$a= {a},\, b={b}$", ls=style)
    >>> ax.legend()
    >>> ax.set_xlabel("$x$")
    >>> ax.set_title("Gamma distribution cumulative distribution function")
    >>> plt.show()

    The gamma distribution is also available as `scipy.stats.gamma`. Using
    `gdtr` directly can be much faster than calling the ``cdf`` method of
    `scipy.stats.gamma`, especially for small arrays or individual values.
    To get the same results one must use the following parametrization:
    ``stats.gamma(b, scale=1/a).cdf(x)=gdtr(a, b, x)``.

    >>> from scipy.stats import gamma
    >>> a = 2.
    >>> b = 3
    >>> x = 1.
    >>> gdtr_result = gdtr(a, b, x)  # this will often be faster than below
    >>> gamma_dist_result = gamma(b, scale=1/a).cdf(x)
    >>> gdtr_result == gamma_dist_result  # test that results are equal
    True
    )";

const char *gdtrc_doc = R"(
    gdtrc(a, b, x, out=None)

    Gamma distribution survival function.

    Integral from `x` to infinity of the gamma probability density function,

    .. math::

        S(x) = \int_x^\infty \frac{a^b}{\Gamma(b)} t^{b-1} e^{-at}\,dt,

    where :math:`\Gamma` is the gamma function.

    Parameters
    ----------
    a : array_like
        The rate parameter of the gamma distribution, sometimes denoted
        :math:`\beta` (float). It is also the reciprocal of the scale
        parameter :math:`\theta`.
    b : array_like
        The shape parameter of the gamma distribution, sometimes denoted
        :math:`\alpha` (float).
    x : array_like
        The quantile (lower limit of integration; float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        The survival function of the gamma distribution with parameters `a`
        and `b` evaluated at `x`.

    See Also
    --------
    gdtr: Gamma distribution cumulative distribution function
    scipy.stats.gamma: Gamma distribution
    gdtrix

    Notes
    -----
    The evaluation is carried out using the relation to the incomplete gamma
    integral (regularized gamma function).

    Wrapper for the Cephes [1]_ routine `gdtrc`. Calling `gdtrc` directly can
    improve performance compared to the ``sf`` method of `scipy.stats.gamma`
    (see last example below).

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/

    Examples
    --------
    Compute the function for ``a=1`` and ``b=2`` at ``x=5``.

    >>> import numpy as np
    >>> from scipy.special import gdtrc
    >>> import matplotlib.pyplot as plt
    >>> gdtrc(1., 2., 5.)
    0.04042768199451279

    Compute the function for ``a=1``, ``b=2`` at several points by providing
    a NumPy array for `x`.

    >>> xvalues = np.array([1., 2., 3., 4])
    >>> gdtrc(1., 1., xvalues)
    array([0.36787944, 0.13533528, 0.04978707, 0.01831564])

    `gdtrc` can evaluate different parameter sets by providing arrays with
    broadcasting compatible shapes for `a`, `b` and `x`. Here we compute the
    function for three different `a` at four positions `x` and ``b=3``,
    resulting in a 3x4 array.

    >>> a = np.array([[0.5], [1.5], [2.5]])
    >>> x = np.array([1., 2., 3., 4])
    >>> a.shape, x.shape
    ((3, 1), (4,))

    >>> gdtrc(a, 3., x)
    array([[0.98561232, 0.9196986 , 0.80884683, 0.67667642],
           [0.80884683, 0.42319008, 0.17357807, 0.0619688 ],
           [0.54381312, 0.12465202, 0.02025672, 0.0027694 ]])

    Plot the function for four different parameter sets.

    >>> a_parameters = [0.3, 1, 2, 6]
    >>> b_parameters = [2, 10, 15, 20]
    >>> linestyles = ['solid', 'dashed', 'dotted', 'dashdot']
    >>> parameters_list = list(zip(a_parameters, b_parameters, linestyles))
    >>> x = np.linspace(0, 30, 1000)
    >>> fig, ax = plt.subplots()
    >>> for parameter_set in parameters_list:
    ...     a, b, style = parameter_set
    ...     gdtrc_vals = gdtrc(a, b, x)
    ...     ax.plot(x, gdtrc_vals, label=fr"$a= {a},\, b={b}$", ls=style)
    >>> ax.legend()
    >>> ax.set_xlabel("$x$")
    >>> ax.set_title("Gamma distribution survival function")
    >>> plt.show()

    The gamma distribution is also available as `scipy.stats.gamma`.
    Using `gdtrc` directly can be much faster than calling the ``sf`` method
    of `scipy.stats.gamma`, especially for small arrays or individual
    values. To get the same results one must use the following parametrization:
    ``stats.gamma(b, scale=1/a).sf(x)=gdtrc(a, b, x)``.

    >>> from scipy.stats import gamma
    >>> a = 2
    >>> b = 3
    >>> x = 1.
    >>> gdtrc_result = gdtrc(a, b, x)  # this will often be faster than below
    >>> gamma_dist_result = gamma(b, scale=1/a).sf(x)
    >>> gdtrc_result == gamma_dist_result  # test that results are equal
    True
    )";

const char *owens_t_doc = R"(
    owens_t(h, a, out=None)

    Owen's T Function.

    The function T(h, a) gives the probability of the event
    (X > h and 0 < Y < a * X) where X and Y are independent
    standard normal random variables.

    Parameters
    ----------
    h : array_like
        Input value.
    a : array_like
        Input value.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    t: scalar or ndarray
        Probability of the event (X > h and 0 < Y < a * X),
        where X and Y are independent standard normal random variables.

    References
    ----------
    .. [1] M. Patefield and D. Tandy, "Fast and accurate calculation of
           Owen's T Function", Statistical Software vol. 5, pp. 1-25, 2000.

    Examples
    --------
    >>> from scipy import special
    >>> a = 3.5
    >>> h = 0.78
    >>> special.owens_t(h, a)
    0.10877216734852274
    )";

const char *pdtr_doc = R"(
    pdtr(k, m, out=None)

    Poisson cumulative distribution function.

    Defined as the probability that a Poisson-distributed random
    variable with event rate :math:`m` is less than or equal to
    :math:`k`. More concretely, this works out to be [1]_

    .. math::

       \exp(-m) \sum_{j = 0}^{\lfloor{k}\rfloor} \frac{m^j}{j!}.

    Parameters
    ----------
    k : array_like
        Number of occurrences (nonnegative, real)
    m : array_like
        Shape parameter (nonnegative, real)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the Poisson cumulative distribution function

    See Also
    --------
    pdtrc : Poisson survival function
    pdtrik : inverse of `pdtr` with respect to `k`
    pdtri : inverse of `pdtr` with respect to `m`

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Poisson_distribution

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is a cumulative distribution function, so it converges to 1
    monotonically as `k` goes to infinity.

    >>> sc.pdtr([1, 10, 100, np.inf], 1)
    array([0.73575888, 0.99999999, 1.        , 1.        ])

    It is discontinuous at integers and constant between integers.

    >>> sc.pdtr([1, 1.5, 1.9, 2], 1)
    array([0.73575888, 0.73575888, 0.73575888, 0.9196986 ])

    )";

const char *pdtrc_doc = R"(
    pdtrc(k, m, out=None)

    Poisson survival function.

    Returns the sum of the terms from k+1 to infinity of the Poisson
    distribution: sum(exp(-m) * m**j / j!, j=k+1..inf) = gammainc(
    k+1, m). Arguments must both be non-negative doubles.

    Parameters
    ----------
    k : array_like
        Number of occurrences (nonnegative, real)
    m : array_like
        Shape parameter (nonnegative, real)
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the Poisson survival function

    See Also
    --------
    pdtr : Poisson cumulative distribution function
    pdtrik : inverse of `pdtr` with respect to `k`
    pdtri : inverse of `pdtr` with respect to `m`

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is a survival function, so it decreases to 0
    monotonically as `k` goes to infinity.

    >>> k = np.array([1, 10, 100, np.inf])
    >>> sc.pdtrc(k, 1)
    array([2.64241118e-001, 1.00477664e-008, 3.94147589e-161, 0.00000000e+000])

    It can be expressed in terms of the lower incomplete gamma
    function `gammainc`.

    >>> sc.gammainc(k + 1, 1)
    array([2.64241118e-001, 1.00477664e-008, 3.94147589e-161, 0.00000000e+000])

    )";

const char *poch_doc = R"(
    poch(z, m, out=None)

    Pochhammer symbol.

    The Pochhammer symbol (rising factorial) is defined as

    .. math::

        (z)_m = \frac{\Gamma(z + m)}{\Gamma(z)}

    For positive integer `m` it reads

    .. math::

        (z)_m = z (z + 1) ... (z + m - 1)

    See [DLMF]_ for more details.

    Parameters
    ----------
    z, m : array_like
        Real-valued arguments.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The value of the function.

    References
    ----------
    .. [DLMF] Nist, Digital Library of Mathematical Functions
        https://dlmf.nist.gov/5.2#iii

    Examples
    --------
    >>> import scipy.special as sc

    It is 1 when m is 0.

    >>> sc.poch([1, 2, 3, 4], 0)
    array([1., 1., 1., 1.])

    For z equal to 1 it reduces to the factorial function.

    >>> sc.poch(1, 5)
    120.0
    >>> 1 * 2 * 3 * 4 * 5
    120

    It can be expressed in terms of the gamma function.

    >>> z, m = 3.7, 2.1
    >>> sc.poch(z, m)
    20.529581933776953
    >>> sc.gamma(z + m) / sc.gamma(z)
    20.52958193377696

    )";

const char *round_doc = R"(
    round(x, out=None)

    Round to the nearest integer.

    Returns the nearest integer to `x`.  If `x` ends in 0.5 exactly,
    the nearest even integer is chosen.

    Parameters
    ----------
    x : array_like
        Real valued input.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        The nearest integers to the elements of `x`. The result is of
        floating type, not integer type.

    Examples
    --------
    >>> import scipy.special as sc

    It rounds to even.

    >>> sc.round([0.5, 1.5])
    array([0., 2.])

    )";

const char *gdtria_doc = R"(
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
        `b` parameter values of `gdtr(a, b, x)`. `b` is the "shape" parameter
        of the gamma distribution.
    x : array_like
        Nonnegative real values, from the domain of the gamma distribution.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`.  `out` is then the
        array returned by the function.

    Returns
    -------
    a : scalar or ndarray
        Values of the `a` parameter such that ``p = gdtr(a, b, x)`.  ``1/a``
        is the "scale" parameter of the gamma distribution.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtrib : Inverse with respect to `b` of `gdtr(a, b, x)`.
    gdtrix : Inverse with respect to `x` of `gdtr(a, b, x)`.
    gammaincinv : Inverse of the incomplete regularized gamma function.

    Notes
    -----
    `gdtria` is implemented in terms of the incomplete gamma inverse as
    ``gdtria(p, b, x) = gammaincinv(b, p)/x``.

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
    )";

const char *gdtrib_doc = R"(
    gdtrib(a, p, x, out=None)

    Inverse of `gdtr` vs b.

    Returns the inverse with respect to the parameter `b` of ``p =
    gdtr(a, b, x)``, the cumulative distribution function of the gamma
    distribution.

    Parameters
    ----------
    a : array_like
        `a` parameter values of ``gdtr(a, b, x)`. ``1/a`` is the "scale"
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
    b : scalar or ndarray
        Values of the `b` parameter such that `p = gdtr(a, b, x)`.  `b` is
        the "shape" parameter of the gamma distribution.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtria : Inverse with respect to `a` of `gdtr(a, b, x)`.
    gdtrix : Inverse with respect to `x` of `gdtr(a, b, x)`.

    Notes
    -----

    The cumulative distribution function `p` is computed using the Cephes [1]_
    routines `igam` and `igamc`. Computation of `b` involves a search for a value
    that produces the desired value of `p` using Chandrupatla's bracketing
    root finding algorithm [2]_.

    Note that there are some edge cases where `gdtrib` is extended by taking
    limits where they are uniquely defined. In particular
    ``x == 0`` with ``p > 0`` and ``p == 0`` with ``x > 0``.
    For these edge cases, a numerical result will be returned for
    ``gdtrib(a, p, x)`` even though ``gdtr(a, gdtrib(a, p, x), x)`` is
    undefined.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           https://netlib.org/cephes/
    .. [2] Chandrupatla, Tirupathi R.
           "A new hybrid quadratic/bisection algorithm for finding the zero of a
           nonlinear function without using derivatives".
           Advances in Engineering Software, 28(3), 145-149.
           :doi:`10.1016/s0965-9978(96)00051-8`.

    Examples
    --------
    First evaluate `gdtr`.

    >>> from scipy.special import gdtr, gdtrib
    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtrib(1.2, p, 5.6)
    3.3999999999999995
    )";

const char *gdtrix_doc = R"(
    gdtrix(a, b, p, out=None)

    Inverse of `gdtr` vs x.

    Returns the inverse with respect to the parameter `x` of ``p =
    gdtr(a, b, x)``, the cumulative distribution function of the gamma
    distribution. This is also known as the pth quantile of the
    distribution.

    Parameters
    ----------
    a : array_like
        `a` parameter values of ``gdtr(a, b, x)``. ``1/a`` is the "scale"
        parameter of the gamma distribution.
    b : array_like
        `b` parameter values of ``gdtr(a, b, x)``. `b` is the "shape" parameter
        of the gamma distribution.
    p : array_like
        Probability values.
    out : ndarray, optional
        If a fourth argument is given, it must be a numpy.ndarray whose size
        matches the broadcast result of `a`, `b` and `x`. `out` is then the
        array returned by the function.

    Returns
    -------
    x : scalar or ndarray
        Values of the `x` parameter such that `p = gdtr(a, b, x)`.

    See Also
    --------
    gdtr : CDF of the gamma distribution.
    gdtria : Inverse with respect to `a` of ``gdtr(a, b, x)``.
    gdtrib : Inverse with respect to `b` of ``gdtr(a, b, x)``.
    gammaincinv : Inverse of the incomplete regularized gamma function.

    Notes
    -----
    `gdtrix` is implemented in terms of the incomplete gamma inverse as
    ``gdtrix(a, b, p) = gammaincinv(b, p)/a``.

    Examples
    --------
    First evaluate `gdtr`.

    >>> from scipy.special import gdtr, gdtrix
    >>> p = gdtr(1.2, 3.4, 5.6)
    >>> print(p)
    0.94378087442

    Verify the inverse.

    >>> gdtrix(1.2, 3.4, p)
    5.6
    )";

const char *kolmogorov_doc = R"(
    kolmogorov(x, out=None)

    Complementary cumulative distribution function of the Kolmogorov distribution.

    Returns the survival function of Kolmogorov's limiting distribution,
    i.e., the distribution of :math:`\sqrt{n}\, D_n` as :math:`n \to \infty`,
    where :math:`D_n` is the Kolmogorov--Smirnov statistic:

    .. math::

       D_n = \sup_x \left| F_n(x) - F(x) \right|

    with :math:`F_n` the empirical cumulative distribution function (ECDF) and
    :math:`F` the theoretical CDF. Specifically, this function computes:

    .. math::

       \mathbb{P}(\sqrt{n}\, D_n > x)

    Parameters
    ----------
    x : array_like
        Absolute deviation between the Empirical CDF (ECDF) and the target CDF,
        multiplied by :math:`\sqrt{n}`.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        The probability that the test statistic exceeds `x`, in the range
        :math:`[0, 1]`.

    See Also
    --------
    kolmogi : Inverse survival function of the Kolmogorov distribution.
    scipy.stats.kstwobign : Provides the functionality as a continuous distribution.
    smirnov, smirnovi : Functions for the one-sided distribution.

    Notes
    -----
    `kolmogorov` is used by `scipy.stats.kstest` in the application of the
    Kolmogorov-Smirnov Goodness of Fit test. For historical reasons this
    function is exposed in `scipy.special`, but the recommended way to achieve
    the most accurate CDF/SF/PDF/PPF/ISF computations is to use the
    `scipy.stats.kstwobign` distribution.

    References
    ----------
    .. [1] Marsaglia, G., Tsang, W. W., & Wang, J. (2003). "Evaluating
       Kolmogorov's distribution." *Journal of Statistical Software*, 8(18), 1-4.
    .. [2] "Kolmogorov-Smirnov test", Wikipedia,
       https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test

    Examples
    --------
    Show the probability of a gap at least as big as 0, 0.5 and 1.0.

    >>> import numpy as np
    >>> from scipy.special import kolmogorov
    >>> from scipy.stats import kstwobign
    >>> kolmogorov([0, 0.5, 1.0])
    array([ 1.        ,  0.96394524,  0.26999967])

    Compare a sample of size 1000 drawn from a Laplace(0, 1) distribution against
    the target distribution, a Normal(0, 1) distribution.

    >>> from scipy.stats import norm, laplace
    >>> rng = np.random.default_rng()
    >>> n = 1000
    >>> lap01 = laplace(0, 1)
    >>> x = np.sort(lap01.rvs(n, random_state=rng))
    >>> np.mean(x), np.std(x)
    (np.float64(-0.05841730131499543), np.float64(1.3968109101997568))

    Construct the Empirical CDF and the K-S statistic Dn.

    >>> target = norm(0,1)  # Normal mean 0, stddev 1
    >>> cdfs = target.cdf(x)
    >>> ecdfs = np.arange(n+1, dtype=float)/n
    >>> gaps = np.column_stack([cdfs - ecdfs[:n], ecdfs[1:] - cdfs])
    >>> Dn = np.max(gaps)
    >>> Kn = np.sqrt(n) * Dn
    >>> print(f'Dn={Dn:f}, sqrt(n)*Dn={Kn:f}')
    Dn=0.043363, sqrt(n)*Dn=1.371265
    >>> print(chr(10).join(['For a Laplace sample tested against N(0, 1):',
    ...   f' the approximate Kolmogorov probability that sqrt(n)*Dn>={Kn:f} '
    ...   f'is {kolmogorov(Kn):f}',
    ...   f' the approximate Kolmogorov probability that sqrt(n)*Dn<={Kn:f} '
    ...   f'is {kstwobign.cdf(Kn):f}']))
    For a Laplace sample tested against N(0, 1):
     the approximate Kolmogorov probability that sqrt(n)*Dn>=1.371265 is 0.046533
     the approximate Kolmogorov probability that sqrt(n)*Dn<=1.371265 is 0.953467

    Plot the Empirical CDF against the target N(0, 1) CDF.

    >>> import matplotlib.pyplot as plt
    >>> plt.step(np.concatenate([[-3], x]), ecdfs, where='post', label='Empirical CDF')
    >>> x3 = np.linspace(-3, 3, 100)
    >>> plt.plot(x3, target.cdf(x3), label='CDF for N(0, 1)')
    >>> plt.ylim([0, 1]); plt.grid(True); plt.legend();
    >>> # Add vertical lines marking Dn+ and Dn-
    >>> iminus, iplus = np.argmax(gaps, axis=0)
    >>> plt.vlines([x[iminus]], ecdfs[iminus], cdfs[iminus],
    ...            color='r', linestyle='dashed', lw=4)
    >>> plt.vlines([x[iplus]], cdfs[iplus], ecdfs[iplus+1],
    ...            color='r', linestyle='dashed', lw=4)
    >>> plt.show()
    )";

const char *kolmogi_doc = R"(
    kolmogi(p, out=None)

    Inverse Survival Function of Kolmogorov distribution.

    It is the inverse function to `kolmogorov`.
    Returns y such that ``kolmogorov(y) == p``.

    Parameters
    ----------
    p : float array_like
        Probability
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The value(s) of kolmogi(p)

    See Also
    --------
    kolmogorov : The Survival Function for the distribution
    scipy.stats.kstwobign : Provides the functionality as a continuous distribution
    smirnov, smirnovi : Functions for the one-sided distribution

    Notes
    -----
    For historical reasons this function is exposed in `scipy.special`, but the
    recommended way to achieve the most accurate CDF/SF/PDF/PPF/ISF computations is to
    use the `scipy.stats.kstwobign` distribution.

    Examples
    --------
    >>> from scipy.special import kolmogi
    >>> kolmogi([0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    array([        inf,  1.22384787,  1.01918472,  0.82757356,  0.67644769,
            0.57117327,  0.        ])

    )";

const char *nrdtrimn_doc = R"(
    nrdtrimn(p, std, x, out=None)

    Calculate mean of normal distribution given other params.

    Parameters
    ----------
    p : array_like
        CDF values, in range (0, 1].
    std : array_like
        Standard deviation.
    x : array_like
        Quantiles, i.e. the upper limit of integration.
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    mn : scalar or ndarray
        The mean of the normal distribution.

    See Also
    --------
    scipy.stats.norm : Normal distribution
    ndtr : Standard normal cumulative probability distribution
    ndtri : Inverse of standard normal CDF with respect to quantile
    nrdtrisd : Inverse of normal distribution CDF with respect to
               standard deviation

    Examples
    --------
    `nrdtrimn` can be used to recover the mean of a normal distribution
    if we know the CDF value `p` for a given quantile `x` and the
    standard deviation `std`. First, we calculate
    the normal distribution CDF for an exemplary parameter set.

    >>> from scipy.stats import norm
    >>> mean = 3.
    >>> std = 2.
    >>> x = 6.
    >>> p = norm.cdf(x, loc=mean, scale=std)
    >>> p
    0.9331927987311419

    Verify that `nrdtrimn` returns the original value for `mean`.

    >>> from scipy.special import nrdtrimn
    >>> nrdtrimn(p, std, x)
    3.0000000000000004
    )";

const char *nrdtrisd_doc = R"(
    nrdtrisd(mn, p, x, out=None)

    Calculate standard deviation of normal distribution given other params.

    Parameters
    ----------
    mn : scalar or ndarray
        The mean of the normal distribution.
    p : array_like
        CDF values, in range (0, 1].
    x : array_like
        Quantiles, i.e. the upper limit of integration.

    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    std : scalar or ndarray
        Standard deviation.

    See Also
    --------
    scipy.stats.norm : Normal distribution
    ndtr : Standard normal cumulative probability distribution
    ndtri : Inverse of standard normal CDF with respect to quantile
    nrdtrimn : Inverse of normal distribution CDF with respect to
               mean

    Examples
    --------
    `nrdtrisd` can be used to recover the standard deviation of a normal
    distribution if we know the CDF value `p` for a given quantile `x` and
    the mean `mn`. First, we calculate the normal distribution CDF for an
    exemplary parameter set.

    >>> from scipy.stats import norm
    >>> mean = 3.
    >>> std = 2.
    >>> x = 6.
    >>> p = norm.cdf(x, loc=mean, scale=std)
    >>> p
    0.9331927987311419

    Verify that `nrdtrisd` returns the original value for `std`.

    >>> from scipy.special import nrdtrisd
    >>> nrdtrisd(mean, p, x)
    2.0000000000000004
    )";

const char *_gen_harmonic_doc = R"(
    _gen_harmonic(n, a)

    Internal private function.

    Compute sum_{i=1}^{n} i**-a for 1 <= m <= n.

    This is the generalized harmonic number.

    nan is returned if n < 1.

    When n is type double, if it is nan, nan is returned.
    Otherwise when n is double it is assumed to be an integer
    value that is between 0 and 2**53.  This is not checked by
    the function.

    This function is used in scipy.stats.zipfian.
    )";

const char *_normalized_gen_harmonic_doc = R"(
    _normalized_gen_harmonic(j, k, n, a)

    Internal private function.

    Compute (sum_{i=j}^{k} i**-a)/(sum_{i=1}^{n} i**-a) for 1 <= j <= k <= n.

    When j, k and n are type double, nan is returned if any are nan.
    Otherwise when the type is double it is assumed that i, j and k have integer
    values that are between 0 and 2**53.  This is not checked by the function.
    Failure to ensure this condition could result in invalid results and
    possibly an infinite loop in the underlying C++ code.

    This function is used in scipy.stats.zipfian.
    )";

const char *_von_mises_cdf_doc = R"(
    _von_mises_cdf(kappa, x, out=None)

    Internal private function.

    Compute the CDF of the von Mises distribution with concentration
    ``kappa``, extended periodically over ``x``.
    )";

const char *besselpoly_doc = R"(
    besselpoly(a, lmb, nu, out=None)

    Weighted integral of the Bessel function of the first kind.

    Computes

    .. math::

       \int_0^1 x^\lambda J_\nu(2 a x) \, dx

    where :math:`J_\nu` is a Bessel function and :math:`\lambda=lmb`,
    :math:`\nu=nu`.

    Parameters
    ----------
    a : array_like
        Scale factor inside the Bessel function.
    lmb : array_like
        Power of `x`
    nu : array_like
        Order of the Bessel function.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Value of the integral.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Evaluate the function for one parameter set.

    >>> from scipy.special import besselpoly
    >>> besselpoly(1, 1, 1)
    0.24449718372863877

    Evaluate the function for different scale factors.

    >>> import numpy as np
    >>> factors = np.array([0., 3., 6.])
    >>> besselpoly(factors, 1, 1)
    array([ 0.        , -0.00549029,  0.00140174])

    Plot the function for varying powers, orders and scales.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> powers = np.linspace(0, 10, 100)
    >>> orders = [1, 2, 3]
    >>> scales = [1, 2]
    >>> all_combinations = [(order, scale) for order in orders
    ...                     for scale in scales]
    >>> for order, scale in all_combinations:
    ...     ax.plot(powers, besselpoly(scale, powers, order),
    ...             label=rf"$\nu={order}, a={scale}$")
    >>> ax.legend()
    >>> ax.set_xlabel(r"$\lambda$")
    >>> ax.set_ylabel(r"$\int_0^1 x^{\lambda} J_{\nu}(2ax)\,dx$")
    >>> plt.show()
    )";

const char *beta_doc = R"(
    beta(a, b, out=None)

    Beta function.

    This function is defined in [1]_ as

    .. math::

        B(a, b) = \int_0^1 t^{a-1}(1-t)^{b-1}dt
                = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)},

    where :math:`\Gamma` is the gamma function.

    Parameters
    ----------
    a, b : array_like
        Real-valued arguments
    out : ndarray, optional
        Optional output array for the function result

    Returns
    -------
    scalar or ndarray
        Value of the beta function

    See Also
    --------
    gamma : the gamma function
    betainc :  the regularized incomplete beta function
    betaln : the natural logarithm of the absolute
             value of the beta function

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions,
           Eq. 5.12.1. https://dlmf.nist.gov/5.12

    Examples
    --------
    >>> import scipy.special as sc

    The beta function relates to the gamma function by the
    definition given above:

    >>> sc.beta(2, 3)
    0.08333333333333333
    >>> sc.gamma(2)*sc.gamma(3)/sc.gamma(2 + 3)
    0.08333333333333333

    As this relationship demonstrates, the beta function
    is symmetric:

    >>> sc.beta(1.7, 2.4)
    0.16567527689031739
    >>> sc.beta(2.4, 1.7)
    0.16567527689031739

    This function satisfies :math:`B(1, b) = 1/b`:

    >>> sc.beta(1, 4)
    0.25
    )";

const char *betaln_doc = R"(
    betaln(a, b, out=None)

    Natural logarithm of absolute value of beta function.

    Computes ``ln(abs(beta(a, b)))``.

    Parameters
    ----------
    a, b : array_like
        Positive, real-valued parameters
    out : ndarray, optional
        Optional output array for function values

    Returns
    -------
    scalar or ndarray
        Value of the betaln function

    See Also
    --------
    gamma : the gamma function
    betainc :  the regularized incomplete beta function
    beta : the beta function

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import betaln, beta

    Verify that, for moderate values of ``a`` and ``b``, ``betaln(a, b)``
    is the same as ``log(beta(a, b))``:

    >>> betaln(3, 4)
    -4.0943445622221

    >>> np.log(beta(3, 4))
    -4.0943445622221

    In the following ``beta(a, b)`` underflows to 0, so we can't compute
    the logarithm of the actual value.

    >>> a = 400
    >>> b = 900
    >>> beta(a, b)
    0.0

    We can compute the logarithm of ``beta(a, b)`` by using `betaln`:

    >>> betaln(a, b)
    -804.3069951764146
    )";

const char *cbrt_doc = R"(
    cbrt(x, out=None)

    Element-wise cube root of `x`.

    Parameters
    ----------
    x : array_like
        `x` must contain real numbers.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        The cube root of each value in `x`.

    Examples
    --------
    >>> from scipy.special import cbrt

    >>> cbrt(8)
    2.0
    >>> cbrt([-8, -3, 0.125, 1.331])
    array([-2.        , -1.44224957,  0.5       ,  1.1       ])
    )";

const char *cosdg_doc = R"(
    cosdg(x, out=None)

    Cosine of the angle `x` given in degrees.

    Parameters
    ----------
    x : array_like
        Angle, given in degrees.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Cosine of the input.

    See Also
    --------
    sindg, tandg, cotdg

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is more accurate than using cosine directly.

    >>> x = 90 + 180 * np.arange(3)
    >>> sc.cosdg(x)
    array([-0.,  0., -0.])
    >>> np.cos(x * np.pi / 180)
    array([ 6.1232340e-17, -1.8369702e-16,  3.0616170e-16])
    )";

const char *cosm1_doc = R"(
    cosm1(x, out=None)

    Compute ``cos(x) - 1``, especially when `x` is near zero.

    Parameters
    ----------
    x : array_like
        Real valued argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of ``cos(x) - 1``.

    See Also
    --------
    expm1, log1p

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is more accurate than computing ``cos(x) - 1`` directly for
    ``x`` around 0.

    >>> x = 1e-30
    >>> np.cos(x) - 1
    0.0
    >>> sc.cosm1(x)
    -5.0000000000000005e-61
    )";

const char *_sinpi_doc = R"(
    Internal function, do not use.
    )";

const char *sindg_doc = R"(
    sindg(x, out=None)

    Sine of the angle `x` given in degrees.

    Parameters
    ----------
    x : array_like
        Angle, given in degrees.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Sine at the input.

    See Also
    --------
    cosdg, tandg, cotdg

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is more accurate than using sine directly.

    >>> x = 180 * np.arange(3)
    >>> sc.sindg(x)
    array([ 0., -0.,  0.])
    >>> np.sin(x * np.pi / 180)
    array([ 0.0000000e+00,  1.2246468e-16, -2.4492936e-16])
    )";

const char *_zeta_doc = R"(
    _zeta(x, q)

    Internal function, Hurwitz zeta.

    )";

const char *zetac_doc = R"(
    zetac(x, out=None)

    Riemann zeta function minus 1.

    This function is defined as

    .. math::

        \zeta(x) - 1 = \sum_{k=2}^{\infty} \frac{1}{k^x}

    where ``x > 1``.  For ``x < 1`` the analytic continuation is
    computed. For more information on the Riemann zeta function, see
    [dlmf]_.

    Parameters
    ----------
    x : array_like of float
        Values at which to compute zeta(x) - 1 (must be real).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of zeta(x) - 1.

    See Also
    --------
    zeta

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/25

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import zetac, zeta

    Some special values:

    >>> zetac(2), np.pi**2/6 - 1
    (0.64493406684822641, 0.6449340668482264)

    >>> zetac(-1), -1.0/12 - 1
    (-1.0833333333333333, -1.0833333333333333)

    Compare ``zetac(x)`` to ``zeta(x) - 1`` for large `x`:

    >>> zetac(60), zeta(60) - 1
    (8.673617380119933e-19, 0.0)
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
    The Airy functions :math:`\operatorname{Ai}` and :math:`\operatorname{Bi}` are two
    independent solutions of

    .. math:: y''(z) = z y(z).

    For real :math:`z` in :math:`[-10, 10]`, the computation is carried out by calling
    the Cephes [1]_ `airy` routine, which uses power series summation
    for small :math:`z` and rational minimax approximations for large :math:`z`.

    Outside this range, the AMOS [2]_ `zairy` and `zbiry` routines are
    employed.  They are computed using power series for :math:`|z| < 1` and
    the following relations to modified Bessel functions for larger :math:`z`
    (where :math:`t \equiv 2 z^{3/2}/3`):

    .. math::

        \operatorname{Ai}(z) = \frac{1}{\pi}\sqrt{\frac{z}{3}} \, K_{1/3}(t)

        \operatorname{Ai}'(z) = -\frac{z}{\pi \sqrt{3}} \, K_{2/3}(t)

        \operatorname{Bi}(z) = \sqrt{\frac{z}{3}} \left(I_{-1/3}(t) + I_{1/3}(t)\right)

        \operatorname{Bi}'(z) = \frac{z}{\sqrt{3}} \left(I_{-2/3}(t) + I_{2/3}(t)\right)

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/
    .. [2] Donald E. Amos, "AMOS, A Portable Package for Bessel Functions
           of a Complex Argument and Nonnegative Order",
           http://netlib.org/amos/

    Examples
    --------
    Compute the Airy functions on the interval :math:`[-15, 5]`.

    >>> import numpy as np
    >>> from scipy import special
    >>> x = np.linspace(-15, 5, 201)
    >>> ai, aip, bi, bip = special.airy(x)

    Plot :math:`\operatorname{Ai}(x)` and :math:`\operatorname{Bi}(x)`.

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
    bei : the corresponding imaginary part
    berp : the derivative of ber
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
    ber : Kelvin function ber

    References
    ----------
    .. [dlmf] NIST, Digital Library of Mathematical Functions,
        https://dlmf.nist.gov/10#PT5

    )";

const char *bdtrik_doc = R"(
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

    Implemented using the Boost library.

    References
    ----------
    .. [1] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [2] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/8.17.5#E5
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
    x, y : array_like
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

const char *boxcox_doc = R"(
    boxcox(x, lmbda, out=None)

    Compute the Box-Cox transformation.

    The Box-Cox transformation is

    .. math::

        y = \begin{cases}
            (x^\lambda - 1) / \lambda & \text{if } \lambda \neq 0 \\
            \log(x) & \text{if } \lambda = 0
        \end{cases}

    Returns ``nan`` if :math:`x < 0`.
    Returns ``-inf`` if :math:`x = 0` and :math:`\lambda \leq 0`.

    Parameters
    ----------
    x : array_like
        Data to be transformed.
    lmbda : array_like
        Power parameter :math:`\lambda` of the Box-Cox transform.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    y : scalar or ndarray
        Transformed data.

    See Also
    --------
    boxcox1p : Box-Cox transformation of ``1 + x``.
    inv_boxcox : Inverse of the Box-Cox transformation.

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
    )";

const char *boxcox1p_doc = R"(
    boxcox1p(x, lmbda, out=None)

    Compute the Box-Cox transformation of :math:`1 + x`.

    The Box-Cox transformation computed by `boxcox1p` is

    .. math::

        y = \begin{cases}
            ((1+x)^\lambda - 1) / \lambda & \text{if } \lambda \neq 0 \\
            \log(1+x) & \text{if } \lambda = 0
        \end{cases}

    Returns ``nan`` if :math:`x < -1`.
    Returns ``-inf`` if :math:`x = -1` and :math:`\lambda \leq 0`.

    Parameters
    ----------
    x : array_like
        Data to be transformed.
    lmbda : array_like
        Power parameter :math:`\lambda` of the Box-Cox transform.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    y : scalar or ndarray
        Transformed data.

    See Also
    --------
    boxcox : Box-Cox transformation.
    inv_boxcox1p : Inverse of the Box-Cox transformation of ``1 + x``.

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
    )";

const char *inv_boxcox_doc = R"(
    inv_boxcox(y, lmbda, out=None)

    Compute the inverse of the Box-Cox transformation.

    Find :math:`x` such that

    .. math::

        y = \begin{cases}
            (x^\lambda - 1) / \lambda & \text{if } \lambda \neq 0 \\
            \log(x) & \text{if } \lambda = 0
        \end{cases}

    Parameters
    ----------
    y : array_like
        Transformed data (input to the inverse transform).
    lmbda : array_like
        Power parameter :math:`\lambda` of the Box-Cox transform.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    x : scalar or ndarray
        Original data (inverse Box-Cox transform of `y`).

    See Also
    --------
    boxcox : Box-Cox transformation.
    inv_boxcox1p : Inverse of the Box-Cox transformation of ``1 + x``.

    Notes
    -----

    .. versionadded:: 0.16.0

    Examples
    --------
    >>> from scipy.special import boxcox, inv_boxcox
    >>> y = boxcox([1, 4, 10], 2.5)
    >>> inv_boxcox(y, 2.5)
    array([1., 4., 10.])
    )";

const char *inv_boxcox1p_doc = R"(
    inv_boxcox1p(y, lmbda, out=None)

    Compute the inverse of the Box-Cox transformation of :math:`1 + x`.

    Find :math:`x` such that

    .. math::

        y = \begin{cases}
            ((1+x)^\lambda - 1) / \lambda & \text{if } \lambda \neq 0 \\
            \log(1+x) & \text{if } \lambda = 0
        \end{cases}

    Parameters
    ----------
    y : array_like
        Transformed data (input to the inverse transform).
    lmbda : array_like
        Power parameter :math:`\lambda` of the Box-Cox transform.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    x : scalar or ndarray
        Original data (inverse Box-Cox transform of `y`).

    See Also
    --------
    boxcox1p : Box-Cox transformation of ``1 + x``.
    inv_boxcox : Inverse of the Box-Cox transformation.

    Notes
    -----

    .. versionadded:: 0.16.0

    Examples
    --------
    >>> from scipy.special import boxcox1p, inv_boxcox1p
    >>> y = boxcox1p([1, 4, 10], 2.5)
    >>> inv_boxcox1p(y, 2.5)
    array([1., 4., 10.])
    )";

const char *cotdg_doc = R"(
    cotdg(x, out=None)

    Cotangent of the angle `x` given in degrees.

    Parameters
    ----------
    x : array_like
        Angle, given in degrees.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Cotangent at the input.

    See Also
    --------
    sindg, cosdg, tandg

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is more accurate than using cotangent directly.

    >>> x = 90 + 180 * np.arange(3)
    >>> sc.cotdg(x)
    array([0., 0., 0.])
    >>> 1 / np.tan(x * np.pi / 180)
    array([6.1232340e-17, 1.8369702e-16, 3.0616170e-16])
    )";

const char *ellipe_doc = R"(
    ellipe(m, out=None)

    Complete elliptic integral of the second kind.

    This function is defined as

    .. math:: E(m) = \int_0^{\pi/2} [1 - m \sin(t)^2]^{1/2} dt

    Parameters
    ----------
    m : array_like
        Defines the parameter of the elliptic integral.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    E : scalar or ndarray
        Value of the elliptic integral.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipeinc : Incomplete elliptic integral of the second kind
    elliprd : Symmetric elliptic integral of the second kind.
    elliprg : Completely-symmetric elliptic integral of the second kind.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellpe`.

    For ``m > 0`` the computation uses the approximation,

    .. math:: E(m) \approx P(1-m) - (1-m) \log(1-m) Q(1-m),

    where :math:`P` and :math:`Q` are tenth-order polynomials.  For
    ``m < 0``, the relation

    .. math:: E(m) = E(m/(m - 1)) \sqrt(1-m)

    is used.

    The parameterization in terms of :math:`m` follows that of section
    17.2 in [2]_. Other parameterizations in terms of the
    complementary parameter :math:`1 - m`, modular angle
    :math:`\sin^2(\alpha) = m`, or modulus :math:`k^2 = m` are also
    used, so be careful that you choose the correct parameter.

    The Legendre E integral is related to Carlson's symmetric R_D or R_G
    functions in multiple ways [3]_. For example,

    .. math:: E(m) = 2 R_G(0, 1-k^2, 1) .

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [3] NIST Digital Library of Mathematical
           Functions. http://dlmf.nist.gov/, Release 1.0.28 of
           2020-09-15. See Sec. 19.25(i) https://dlmf.nist.gov/19.25#i

    Examples
    --------
    This function is used in finding the circumference of an
    ellipse with semi-major axis `a` and semi-minor axis `b`.

    >>> import numpy as np
    >>> from scipy import special

    >>> a = 3.5
    >>> b = 2.1
    >>> e_sq = 1.0 - b**2/a**2  # eccentricity squared

    Then the circumference is found using the following:

    >>> C = 4*a*special.ellipe(e_sq)  # circumference formula
    >>> C
    17.868899204378693

    When `a` and `b` are the same (meaning eccentricity is 0),
    this reduces to the circumference of a circle.

    >>> 4*a*special.ellipe(0.0)  # formula for ellipse with a = b
    21.991148575128552
    >>> 2*np.pi*a  # formula for circle of radius a
    21.991148575128552
    )";

const char *ellipeinc_doc = R"(
    ellipeinc(phi, m, out=None)

    Incomplete elliptic integral of the second kind.

    This function is defined as

    .. math:: E(\phi, m) = \int_0^{\phi} [1 - m \sin(t)^2]^{1/2} dt

    Parameters
    ----------
    phi : array_like
        amplitude of the elliptic integral.
    m : array_like
        parameter of the elliptic integral.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    E : scalar or ndarray
        Value of the elliptic integral.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    elliprd : Symmetric elliptic integral of the second kind.
    elliprf : Completely-symmetric elliptic integral of the first kind.
    elliprg : Completely-symmetric elliptic integral of the second kind.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellie`.

    Computation uses arithmetic-geometric means algorithm.

    The parameterization in terms of :math:`m` follows that of section
    17.2 in [2]_. Other parameterizations in terms of the
    complementary parameter :math:`1 - m`, modular angle
    :math:`\sin^2(\alpha) = m`, or modulus :math:`k^2 = m` are also
    used, so be careful that you choose the correct parameter.

    The Legendre E incomplete integral can be related to combinations
    of Carlson's symmetric integrals R_D, R_F, and R_G in multiple
    ways [3]_. For example, with :math:`c = \csc^2\phi`,

    .. math::
      E(\phi, m) = R_F(c-1, c-k^2, c)
        - \frac{1}{3} k^2 R_D(c-1, c-k^2, c) .

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [3] NIST Digital Library of Mathematical
           Functions. http://dlmf.nist.gov/, Release 1.0.28 of
           2020-09-15. See Sec. 19.25(i) https://dlmf.nist.gov/19.25#i

    Examples
    --------
    The elliptic integral of the second kind can be used to find the circumference of an
    ellipse with semi-major axis ``a`` and semi-minor axis ``b``.

    >>> import numpy as np
    >>> from scipy.special import ellipeinc
    >>> a, b = 3.5, 2.1
    >>> e = np.sqrt(1.0 - b**2/a**2)  # eccentricity
    >>> 4*a*ellipeinc(np.pi/2, e**2)
    np.float64(17.86889920437869)
    )";

const char *ellipj_doc = R"(
    ellipj(u, m, out=None)

    Jacobi elliptic functions.

    Calculates the Jacobi elliptic functions of parameter `m` less than or equal to
    1, and real argument `u`.

    Parameters
    ----------
    u : array_like
        Argument.
    m : array_like
        Parameter.
    out : tuple of ndarray, optional
        Optional output arrays for the function values

    Returns
    -------
    sn, cn, dn, ph : 4-tuple of scalar or ndarray
        The returned functions::

            sn(u|m), cn(u|m), dn(u|m)

        The value `ph` is such that if ``u = ellipkinc(ph, m)``,
        then ``sn(u|m) = sin(ph)`` and ``cn(u|m) = cos(ph)``.

    See Also
    --------
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind

    Notes
    -----
    Wrapper for the Cephes [1]_ routine ``ellpj``.

    These functions are periodic, with quarter-period on the real axis
    equal to the complete elliptic integral ``ellipk(m)``.

    Relation to incomplete elliptic integral: If ``u = ellipkinc(phi,m)``, then
    ``sn(u|m) = sin(phi)``, and ``cn(u|m) = cos(phi)``. The ``phi`` is called
    the amplitude of `u`.

    Computation is by means of the arithmetic-geometric mean algorithm,
    except when `m` is within 1e-9 of 0 or 1. In the latter case with `m`
    close to 1, the approximation applies only for ``phi < pi/2``.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    The elliptic sine sn(u|m) interpolates between the sine function and the
    hyperbolic tangent when m changes from 0 to 1.

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from scipy.special import ellipj
    >>> u = np.linspace(0, 2*np.pi, 100)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(u, np.sin(u), '--', label='sin(u)')
    >>> for m in (0.2, 0.8, 0.99):
    ...     ax.plot(u, ellipj(u, m)[0], label=f'sn(u|{m})')
    >>> ax.plot(u, np.tanh(u), '--', label='tanh(u)')
    >>> ax.set_xlabel('u')
    >>> ax.legend(loc='lower left')
    >>> plt.show()

    Like for sine and cosine, the squares of elliptic sine and elliptic cosine
    add up to one.

    >>> u = np.linspace(0, 5, 11)
    >>> sn, cn, _, _ = ellipj(u, 0.7)
    >>> sn**2 + cn**2
    array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
    )";

const char *ellipkm1_doc = R"(
    ellipkm1(p, out=None)

    Complete elliptic integral of the first kind around `m` = 1.

    This function is defined as

    .. math:: K(p) = \int_0^{\pi/2} [1 - m \sin(t)^2]^{-1/2} dt

    where `m = 1 - p`.

    Parameters
    ----------
    p : array_like
        Defines the parameter of the elliptic integral as `m = 1 - p`.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    K : scalar or ndarray
        Value of the elliptic integral.

    See Also
    --------
    ellipk : Complete elliptic integral of the first kind
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind
    elliprf : Completely-symmetric elliptic integral of the first kind.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellpk`.

    For ``p <= 1``, computation uses the approximation,

    .. math:: K(p) \approx P(p) - \log(p) Q(p)

    where :math:`P` and :math:`Q` are tenth-order polynomials.  The
    argument `p` is used internally rather than `m` so that the logarithmic
    singularity at ``m = 1`` will be shifted to the origin; this preserves
    maximum accuracy.  For ``p > 1``, the identity

    .. math:: K(p) = K(1/p)/\sqrt{p}

    is used.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/
    .. [2] NIST Digital Library of Mathematical Functions,
           Eq. 19.12.1. https://dlmf.nist.gov/19.12.E1

    Examples
    --------
    >>> from scipy.special import ellipk, ellipkm1
    >>> p = 1e-10
    >>> m = 1-p
    >>> ellipk(m), ellipkm1(p)
    (np.float64(12.899219785017415), np.float64(12.8992198263876))

    In order to decide which one of the two results is closer to the correct
    one, one can use the asymptotic expansion including the next-to-leading
    order [2]_.

    >>> from math import log, sqrt
    >>> log(4/sqrt(p)) + 0.25*p*(log(4/sqrt(p))-1)
    12.8992198263876
    
    We can conclude that for such small values of :math:`p`, `ellipkm1` yields
    the better result. For even smaller values of :math:`p`, the difference
    becomes more apparent.

    >>> p = 1e-15
    >>> m = 1-p
    >>> ellipk(m), ellipkm1(p)
    (np.float64(18.656082357290334), np.float64(18.655682558575236))
    >>> log(4/sqrt(p)) + 0.25*p*(log(4/sqrt(p))-1)
    18.655682558575236

    For even smaller values of :math:`p`, the finite spacing between float
    numbers becomes relevant. Then `ellipkm1` needs to be used in any case.
    )";

const char *ellipk_doc = R"(
    ellipk(m, out=None)

    Complete elliptic integral of the first kind.

    This function is defined as

    .. math:: K(m) = \int_0^{\pi/2} [1 - m \sin(t)^2]^{-1/2} dt

    Parameters
    ----------
    m : array_like
        The parameter of the elliptic integral.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    K : scalar or ndarray
        Value of the elliptic integral.

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind around m = 1
    ellipkinc : Incomplete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind
    elliprf : Completely-symmetric elliptic integral of the first kind.

    Notes
    -----
    For more precision around point m = 1, use `ellipkm1`, which this
    function calls.

    The parameterization in terms of :math:`m` follows that of section
    17.2 in [1]_. Other parameterizations in terms of the
    complementary parameter :math:`1 - m`, modular angle
    :math:`\sin^2(\alpha) = m`, or modulus :math:`k^2 = m` are also
    used, so be careful that you choose the correct parameter.

    The Legendre K integral is related to Carlson's symmetric R_F
    function by [2]_:

    .. math:: K(m) = R_F(0, 1-k^2, 1) .

    References
    ----------
    .. [1] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [2] NIST Digital Library of Mathematical
           Functions. http://dlmf.nist.gov/, Release 1.0.28 of
           2020-09-15. See Sec. 19.25(i) https://dlmf.nist.gov/19.25#i

    Examples
    --------
    The period :math:`T` of a simple pendulum increases with growing oscillation
    amplitude and can be expressed through the complete elliptic integral of the
    first kind. The plot displays the ratio of the period of the pendulum and
    the period :math:`T_0` for small amplitude as a function of the maximum
    angle :math:`\alpha` reached by the pendulum.

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from scipy.special import ellipk
    >>> alpha = np.linspace(0, np.pi, 100)
    >>> m = np.sin(alpha/2)**2
    >>> fig, ax = plt.subplots()
    >>> ax.plot(alpha, ellipk(m)*2/np.pi)
    >>> ax.set_xlabel(r'$\alpha$')
    >>> ax.set_ylabel('$T/T_0$')
    >>> plt.show()
    )";

const char *ellipkinc_doc = R"(
    ellipkinc(phi, m, out=None)

    Incomplete elliptic integral of the first kind.

    This function is defined as

    .. math:: K(\phi, m) = \int_0^{\phi} [1 - m \sin(t)^2]^{-1/2} dt

    This function is also called :math:`F(\phi, m)`.

    Parameters
    ----------
    phi : array_like
        amplitude of the elliptic integral
    m : array_like
        parameter of the elliptic integral
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    K : scalar or ndarray
        Value of the elliptic integral

    See Also
    --------
    ellipkm1 : Complete elliptic integral of the first kind, near `m` = 1
    ellipk : Complete elliptic integral of the first kind
    ellipe : Complete elliptic integral of the second kind
    ellipeinc : Incomplete elliptic integral of the second kind
    elliprf : Completely-symmetric elliptic integral of the first kind.

    Notes
    -----
    Wrapper for the Cephes [1]_ routine `ellik`.  The computation is
    carried out using the arithmetic-geometric mean algorithm.

    The parameterization in terms of :math:`m` follows that of section
    17.2 in [2]_. Other parameterizations in terms of the
    complementary parameter :math:`1 - m`, modular angle
    :math:`\sin^2(\alpha) = m`, or modulus :math:`k^2 = m` are also
    used, so be careful that you choose the correct parameter.

    The Legendre K incomplete integral (or F integral) is related to
    Carlson's symmetric R_F function [3]_.
    Setting :math:`c = \csc^2\phi`,

    .. math:: F(\phi, m) = R_F(c-1, c-k^2, c) .

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [3] NIST Digital Library of Mathematical
           Functions. http://dlmf.nist.gov/, Release 1.0.28 of
           2020-09-15. See Sec. 19.25(i) https://dlmf.nist.gov/19.25#i

    Examples
    --------
    >>> from scipy.special import ellipkinc
    >>> phi = 0.3
    >>> m = 0.8
    >>> u = ellipkinc(phi, m)
    >>> u
    np.float64(0.30365239221539364)

    The result should be consistent with the known relations
    for the Jacobi elliptic functions: ``sn(u|m) = sin(phi)``
    and ``cn(u|m) = cos(phi)``.

    >>> from math import cos, pi, sin
    >>> from scipy.special import ellipj
    >>> sn, cn, _, _ = ellipj(u, m)
    >>> sn, sin(phi)
    (np.float64(0.2955202066613395), 0.29552020666133955)
    >>> cn, cos(phi)
    (np.float64(0.9553364891256061), 0.955336489125606)

    For :math:`\phi=\pi/2`, the incomplete elliptic integral should
    equal the complete elliptic integral.

    >>> from scipy.special import ellipk
    >>> ellipkinc(pi/2, m), ellipk(m)
    (np.float64(2.257205326820854), np.float64(2.257205326820854))
    )";

const char *xlogy_doc = R"(
    xlogy(x, y, out=None)

    Compute ``x*log(y)`` so that the result is 0 if ``x = 0``.

    Parameters
    ----------
    x : array_like
        Multiplier
    y : array_like
        Argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    z : scalar or ndarray
        Computed x*log(y)

    Notes
    -----
    The log function used in the computation is the natural log.

    .. versionadded:: 0.13.0

    Examples
    --------
    We can use this function to calculate the binary logistic loss also
    known as the binary cross entropy. This loss function is used for
    binary classification problems and is defined as:

    .. math::
        L = \frac{1}{n} \sum_{i=0}^n -[y_i*\log({y_{pred}}_i) + (1-y_i)*\log(1-{y_{pred}}_i)]

    We can define the parameters `x` and `y` as y and y_pred respectively.
    y is the array of the actual labels which over here can be either 0 or 1.
    y_pred is the array of the predicted probabilities with respect to
    the positive class (1).

    >>> import numpy as np
    >>> from scipy.special import xlogy
    >>> y = np.array([0, 1, 0, 1, 1, 0])
    >>> y_pred = np.array([0.3, 0.8, 0.4, 0.7, 0.9, 0.2])
    >>> n = len(y)
    >>> loss = -(xlogy(y, y_pred) + xlogy(1 - y, 1 - y_pred)).sum()
    >>> loss /= n
    >>> loss
    0.29597052165495025

    A lower loss is usually better as it indicates that the predictions are
    similar to the actual labels. In this example since our predicted
    probabilities are close to the actual labels, we get an overall loss
    that is reasonably low and appropriate.
    )";

const char *xlog1py_doc = R"(
    xlog1py(x, y, out=None)

    Compute ``x*log1p(y)`` so that the result is 0 if ``x = 0``.

    Parameters
    ----------
    x : array_like
        Multiplier
    y : array_like
        Argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    z : scalar or ndarray
        Computed x*log1p(y)

    Notes
    -----

    .. versionadded:: 0.13.0

    Examples
    --------
    This example shows how the function can be used to calculate the log of
    the probability mass function for a geometric discrete random variable.
    The probability mass function of the geometric distribution is defined
    as follows:

    .. math:: f(k) = (1-p)^{k-1} p

    where :math:`p` is the probability of a single success
    and :math:`1-p` is the probability of a single failure
    and :math:`k` is the number of trials to get the first success.

    >>> import numpy as np
    >>> from scipy.special import xlog1py
    >>> p = 0.5
    >>> k = 100
    >>> _pmf = np.power(1 - p, k - 1) * p
    >>> _pmf
    7.888609052210118e-31

    If we take k as a relatively large number the value of the probability
    mass function can become very low. In such cases taking the log of the
    pmf would be more suitable as the log function can change the values
    to a scale that is more appropriate to work with.

    >>> _log_pmf = xlog1py(k - 1, -p) + np.log(p)
    >>> _log_pmf
    -69.31471805599453

    We can confirm that we get a value close to the original pmf value by
    taking the exponential of the log pmf.

    >>> _orig_pmf = np.exp(_log_pmf)
    >>> np.isclose(_pmf, _orig_pmf)
    True
    )";

const char *_log1mexp_doc = R"(
    Internal function, do not use.
    )";

const char *_log1pmx_doc = R"(
    Internal function, do not use.
    )";

const char *log1p_doc = R"(
    log1p(x, out=None)

    Calculates log(1 + x) for use when `x` is near zero.

    Parameters
    ----------
    x : array_like
        Real or complex valued input.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of ``log(1 + x)``.

    See Also
    --------
    expm1, cosm1

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is more accurate than using ``log(1 + x)`` directly for ``x``
    near 0. Note that in the below example ``1 + 1e-17 == 1`` to
    double precision.

    >>> sc.log1p(1e-17)
    1e-17
    >>> np.log(1 + 1e-17)
    0.0
    )";

const char *expm1_doc = R"(
    expm1(x, out=None)

    Compute ``exp(x) - 1``.

    When `x` is near zero, ``exp(x)`` is near 1, so the numerical calculation
    of ``exp(x) - 1`` can suffer from catastrophic loss of precision.
    ``expm1(x)`` is implemented to avoid the loss of precision that occurs when
    `x` is near zero.

    Parameters
    ----------
    x : array_like
        `x` must contain real numbers.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        ``exp(x) - 1`` computed element-wise.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import expm1

    >>> expm1(1.0)
    1.7182818284590451
    >>> expm1([-0.2, -0.1, 0, 0.1, 0.2])
    array([-0.18126925, -0.09516258,  0.        ,  0.10517092,  0.22140276])

    The exact value of ``exp(7.5e-13) - 1`` is::

        7.5000000000028125000000007031250000001318...*10**-13.

    Here is what ``expm1(7.5e-13)`` gives:

    >>> expm1(7.5e-13)
    7.5000000000028135e-13

    Compare that to ``exp(7.5e-13) - 1``, where the subtraction results in
    a "catastrophic" loss of precision:

    >>> np.exp(7.5e-13) - 1
    7.5006667543675576e-13
    )";

const char *exp2_doc = R"(
    exp2(x, out=None)

    Compute ``2**x`` element-wise.

    Parameters
    ----------
    x : array_like
        `x` must contain real numbers.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        ``2**x``, computed element-wise.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import exp2

    >>> exp2(3)
    8.0
    >>> x = np.array([[-1, -0.5, 0], [0.5, 1, 1.5]])
    >>> exp2(x)
    array([[ 0.5       ,  0.70710678,  1.        ],
           [ 1.41421356,  2.        ,  2.82842712]])
    )";

const char *exp10_doc = R"(
    exp10(x, out=None)

    Compute ``10**x`` element-wise.

    Parameters
    ----------
    x : array_like
        `x` must contain real numbers.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        ``10**x``, computed element-wise.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import exp10

    >>> exp10(3)
    1000.0
    >>> x = np.array([[-1, -0.5, 0], [0.5, 1, 1.5]])
    >>> exp10(x)
    array([[  0.1       ,   0.31622777,   1.        ],
           [  3.16227766,  10.        ,  31.6227766 ]])
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
    z : array_like
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

const char *erf_doc = R"(
    erf(z, out=None)

    Error function of real or complex argument.

    .. math::

        \operatorname{erf}(z) = \frac{2}{\sqrt{\pi}} \int_0^z e^{-t^2} dt

    Parameters
    ----------
    z : ndarray
        Input array.
    out : ndarray, optional
        Optional output array for the function values.

    Returns
    -------
    res : scalar or ndarray
        The values of the error function at the given points `z`.

    See Also
    --------
    erfc, erfcx, erfi, erfinv, erfcinv, wofz

    Notes
    -----
    The cumulative distribution function (CDF) of the standard normal distribution can 
    be expressed in terms of the error function as

    .. math::

        \Phi(z) = \frac{1}{2}
        \left[1 + \operatorname{erf} \left(\frac{z}{\sqrt{2}}\right)\right]

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Error_function
    .. [2] Milton Abramowitz and Irene A. Stegun, eds.
           Handbook of Mathematical Functions with Formulas,
           Graphs, and Mathematical Tables. New York: Dover, 1972.
    .. [3] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> z = np.linspace(-3, 3)
    >>> plt.plot(z, special.erf(z))
    >>> plt.xlabel('$z$')
    >>> plt.ylabel('$erf(z)$')
    >>> plt.show()
    )";

const char *erfc_doc = R"(
    erfc(x, out=None)

    Complementary error function.

    The complementary error function is defined as

    .. math::

        \operatorname{erfc}(x) = 1 - \operatorname{erf}(x)

    Parameters
    ----------
    x : array_like
        Real or complex valued argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the complementary error function

    See Also
    --------
    erf, erfi, erfcx, dawsn, wofz

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.erfc(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$erfc(x)$')
    >>> plt.show()
    )";

const char *erfi_doc = R"(
    erfi(z, out=None)

    Imaginary error function.

    The imaginary error function is defined as

    .. math::

        \operatorname{erfi}(z) = -i \operatorname{erf}(i z)

    Parameters
    ----------
    z : array_like
        Real or complex valued argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the imaginary error function

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
    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.erfi(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$erfi(x)$')
    >>> plt.show()
    )";

const char *erfcx_doc = R"(
    erfcx(x, out=None)

    Scaled complementary error function.

    The scaled complementary error function is defined as

    .. math::

        \operatorname{erfcx}(x) = e^{x^2} \operatorname{erfc}(x)

    Parameters
    ----------
    x : array_like
        Real or complex valued argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Values of the scaled complementary error function


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
    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-3, 3)
    >>> plt.plot(x, special.erfcx(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$erfcx(x)$')
    >>> plt.show()
    )";

const char *expit_doc = R"(
    expit(x, out=None)

    Expit (also known as logistic sigmoid) ufunc for ndarrays.

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

const char *dawsn_doc = R"(
    dawsn(x, out=None)

    Dawson's integral.

    Computes::

        exp(-x**2) * integral(exp(t**2), t=0..x).

    Parameters
    ----------
    x : array_like
        Function parameter.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    y : scalar or ndarray
        Value of the integral.

    See Also
    --------
    wofz, erf, erfc, erfcx, erfi

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-15, 15, num=1000)
    >>> plt.plot(x, special.dawsn(x))
    >>> plt.xlabel('$x$')
    >>> plt.ylabel('$dawsn(x)$')
    >>> plt.show()
    )";

const char *fresnel_doc = R"(
    fresnel(z, out=None)

    Fresnel integrals.

    The Fresnel integrals are defined as

    .. math::

       S(z) &= \int_0^z \sin(\pi t^2 /2) dt, \\
       C(z) &= \int_0^z \cos(\pi t^2 /2) dt.

    See [dlmf]_ for details.

    Parameters
    ----------
    z : array_like
        Real or complex valued argument
    out : 2-tuple of ndarrays, optional
        Optional output arrays for the function results

    Returns
    -------
    S, C : 2-tuple of scalar or ndarray
        Values of the Fresnel integrals

    See Also
    --------
    fresnel_zeros : zeros of the Fresnel integrals

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/7.2#iii

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import fresnel, erf

    First, we verify the following limits:

    .. math::

         \lim_{z \to \infty} S(z) = \lim_{z \to \infty} C(z) = 1/2.

    >>> S, C = fresnel([0.1, 1, 10, 100, np.inf])
    >>> S
    array([0.00052359, 0.43825915, 0.46816998, 0.4968169 , 0.5       ])
    >>> C
    array([0.09999753, 0.7798934 , 0.49989869, 0.4999999 , 0.5       ])

    Next, we verify the following relation to the error function:

    .. math::

        C(z) + i S(z) = \frac{1 + i}{2}
        \operatorname{erf}\left(\frac{\sqrt{\pi}(1 - i)}{2} z\right).


    >>> z = np.linspace(-10, 10)
    >>> S, C = fresnel(z)
    >>> RHS = (1 + 1j)/2 * erf((np.sqrt(np.pi)*(1 - 1j))/2*z)
    >>> np.allclose(C + 1j*S, RHS)
    True

    Finally, we plot :math:`C(z)` against :math:`S(z)` to get the Euler or Cornu spiral.

    >>> import matplotlib.pyplot as plt
    >>> z = np.linspace(-10, 10, num=1000)
    >>> S, C = fresnel(z)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(C, S)
    >>> ax.set_aspect('equal')
    >>> ax.set_xlabel('$C(z)$')
    >>> ax.set_ylabel('$S(z)$')
    >>> plt.show()
    )";

const char *gamma_doc = R"(
    gamma(z, out=None)

    Compute the gamma function.

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

    The gamma function has poles at non-negative integers and the sign
    of infinity as z approaches each pole depends upon the direction in
    which the pole is approached. For this reason, the consistent thing
    is for gamma(z) to return NaN at negative integers, and to return
    -inf when x = -0.0 and +inf when x = 0.0, using the signbit of zero
    to signify the direction in which the origin is being approached. This
    is for instance what is recommended for the gamma function in annex F
    entry 9.5.4 of the Iso C 99 standard [isoc99]_.

    Prior to SciPy version 1.15, ``scipy.special.gamma(z)`` returned ``+inf``
    at each pole. This was fixed in version 1.15, but with the following
    consequence. Expressions where gamma appears in the denominator
    such as

    ``gamma(u) * gamma(v) / (gamma(w) * gamma(x))``

    no longer evaluate to 0 if the numerator is well defined but there is a
    pole in the denominator. Instead such expressions evaluate to NaN. We
    recommend instead using the function `rgamma` for the reciprocal gamma
    function in such cases. The above expression could for instance be written
    as

    ``gamma(u) * gamma(v) * (rgamma(w) * rgamma(x))``

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/5.2#E1
    .. [isoc99] https://www.open-std.org/jtc1/sc22/wg14/www/docs/n1256.pdf

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

const char *gammainc_doc = R"(
    gammainc(a, x, out=None)

    Regularized lower incomplete gamma function.

    It is defined as

    .. math::

        P(a, x) = \frac{1}{\Gamma(a)} \int_0^x t^{a - 1}e^{-t} dt

    for :math:`a > 0` and :math:`x \geq 0`. See [dlmf]_ for details.

    Parameters
    ----------
    a : array_like
        Positive parameter
    x : array_like
        Nonnegative argument
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the lower incomplete gamma function

    See Also
    --------
    gammaincc : regularized upper incomplete gamma function
    gammaincinv : inverse of the regularized lower incomplete gamma function
    gammainccinv : inverse of the regularized upper incomplete gamma function

    Notes
    -----
    The function satisfies the relation ``gammainc(a, x) +
    gammaincc(a, x) = 1`` where `gammaincc` is the regularized upper
    incomplete gamma function.

    The implementation largely follows that of [boost]_.

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical functions
              https://dlmf.nist.gov/8.2#E4
    .. [boost] Maddock et. al., "Incomplete Gamma Functions",
       https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html

    Examples
    --------
    >>> import scipy.special as sc

    It is the CDF of the gamma distribution, so it starts at 0 and
    monotonically increases to 1.

    >>> sc.gammainc(0.5, [0, 1, 10, 100])
    array([0.        , 0.84270079, 0.99999226, 1.        ])

    It is equal to one minus the upper incomplete gamma function.

    >>> a, x = 0.5, 0.4
    >>> sc.gammainc(a, x)
    0.6289066304773024
    >>> 1 - sc.gammaincc(a, x)
    0.6289066304773024
    )";

const char *gammaincc_doc = R"(
    gammaincc(a, x, out=None)

    Regularized upper incomplete gamma function.

    It is defined as

    .. math::

        Q(a, x) = \frac{1}{\Gamma(a)} \int_x^\infty t^{a - 1}e^{-t} dt

    for :math:`a > 0` and :math:`x \geq 0`. See [dlmf]_ for details.

    Parameters
    ----------
    a : array_like
        Positive parameter
    x : array_like
        Nonnegative argument
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the upper incomplete gamma function

    See Also
    --------
    gammainc : regularized lower incomplete gamma function
    gammaincinv : inverse of the regularized lower incomplete gamma function
    gammainccinv : inverse of the regularized upper incomplete gamma function

    Notes
    -----
    The function satisfies the relation ``gammainc(a, x) +
    gammaincc(a, x) = 1`` where `gammainc` is the regularized lower
    incomplete gamma function.

    The implementation largely follows that of [boost]_.

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical functions
              https://dlmf.nist.gov/8.2#E4
    .. [boost] Maddock et. al., "Incomplete Gamma Functions",
       https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html

    Examples
    --------
    >>> import scipy.special as sc

    It is the survival function of the gamma distribution, so it
    starts at 1 and monotonically decreases to 0.

    >>> sc.gammaincc(0.5, [0, 1, 10, 100, 1000])
    array([1.00000000e+00, 1.57299207e-01, 7.74421643e-06, 2.08848758e-45,
           0.00000000e+00])

    It is equal to one minus the lower incomplete gamma function.

    >>> a, x = 0.5, 0.4
    >>> sc.gammaincc(a, x)
    0.37109336952269756
    >>> 1 - sc.gammainc(a, x)
    0.37109336952269756
    )";

const char *gammainccinv_doc = R"(
    gammainccinv(a, y, out=None)

    Inverse of the regularized upper incomplete gamma function.

    Given an input :math:`y` between 0 and 1, returns :math:`x` such
    that :math:`y = Q(a, x)`. Here :math:`Q` is the regularized upper
    incomplete gamma function; see `gammaincc`. This is well-defined
    because the upper incomplete gamma function is monotonic as can
    be seen from its definition in [dlmf]_.

    Parameters
    ----------
    a : array_like
        Positive parameter
    y : array_like
        Argument between 0 and 1, inclusive
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the inverse of the upper incomplete gamma function

    See Also
    --------
    gammaincc : regularized upper incomplete gamma function
    gammainc : regularized lower incomplete gamma function
    gammaincinv : inverse of the regularized lower incomplete gamma function

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/8.2#E4

    Examples
    --------
    >>> import scipy.special as sc

    It starts at infinity and monotonically decreases to 0.

    >>> sc.gammainccinv(0.5, [0, 0.1, 0.5, 1])
    array([       inf, 1.35277173, 0.22746821, 0.        ])

    It inverts the upper incomplete gamma function.

    >>> a, x = 0.5, [0, 0.1, 0.5, 1]
    >>> sc.gammaincc(a, sc.gammainccinv(a, x))
    array([0. , 0.1, 0.5, 1. ])

    >>> a, x = 0.5, [0, 10, 50]
    >>> sc.gammainccinv(a, sc.gammaincc(a, x))
    array([ 0., 10., 50.])
    )";

const char *gammaincinv_doc = R"(
    gammaincinv(a, y, out=None)

    Inverse to the regularized lower incomplete gamma function.

    Given an input :math:`y` between 0 and 1, returns :math:`x` such
    that :math:`y = P(a, x)`. Here :math:`P` is the regularized lower
    incomplete gamma function; see `gammainc`. This is well-defined
    because the lower incomplete gamma function is monotonic as can be
    seen from its definition in [dlmf]_.

    Parameters
    ----------
    a : array_like
        Positive parameter
    y : array_like
        Parameter between 0 and 1, inclusive
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of the inverse of the lower incomplete gamma function

    See Also
    --------
    gammainc : regularized lower incomplete gamma function
    gammaincc : regularized upper incomplete gamma function
    gammainccinv : inverse of the regularized upper incomplete gamma function

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/8.2#E4

    Examples
    --------
    >>> import scipy.special as sc

    It starts at 0 and monotonically increases to infinity.

    >>> sc.gammaincinv(0.5, [0, 0.1 ,0.5, 1])
    array([0.        , 0.00789539, 0.22746821,        inf])

    It inverts the lower incomplete gamma function.

    >>> a, x = 0.5, [0, 0.1, 0.5, 1]
    >>> sc.gammainc(a, sc.gammaincinv(a, x))
    array([0. , 0.1, 0.5, 1. ])

    >>> a, x = 0.5, [0, 10, 25]
    >>> sc.gammaincinv(a, sc.gammainc(a, x))
    array([ 0.        , 10.        , 25.00001465])
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

const char *gammasgn_doc = R"(
    gammasgn(x, out=None)

    Sign of the gamma function.

    It is defined as

    .. math::

       \text{gammasgn}(x) =
       \begin{cases}
         +1 & \Gamma(x) > 0 \\
         -1 & \Gamma(x) < 0
       \end{cases}

    where :math:`\Gamma` is the gamma function; see `gamma`. This
    definition is complete since the gamma function is never zero;
    see the discussion after [dlmf]_.

    Parameters
    ----------
    x : array_like
        Real argument
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Sign of the gamma function

    See Also
    --------
    gamma : the gamma function
    gammaln : log of the absolute value of the gamma function
    loggamma : analytic continuation of the log of the gamma function

    Notes
    -----
    The gamma function can be computed as ``gammasgn(x) *
    np.exp(gammaln(x))``.

    References
    ----------
    .. [dlmf] NIST Digital Library of Mathematical Functions
              https://dlmf.nist.gov/5.2#E1

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is 1 for ``x > 0``.

    >>> sc.gammasgn([1, 2, 3, 4])
    array([1., 1., 1., 1.])

    It alternates between -1 and 1 for negative integers.

    >>> sc.gammasgn([-0.5, -1.5, -2.5, -3.5])
    array([-1.,  1., -1.,  1.])

    It can be used to compute the gamma function.

    >>> x = [1.5, 0.5, -0.5, -1.5]
    >>> sc.gammasgn(x) * np.exp(sc.gammaln(x))
    array([ 0.88622693,  1.77245385, -3.5449077 ,  2.3632718 ])
    >>> sc.gamma(x)
    array([ 0.88622693,  1.77245385, -3.5449077 ,  2.3632718 ])
    )";

const char *hankel1_doc = R"(
    hankel1(v, z, out=None)

    Hankel function of the first kind.

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
    hankel1e : This function with leading exponential behavior stripped off.

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

    Examples
    --------
    For the inhomogeneous Helmholtz equation in :math:`\mathbb{R}^2` subject to radiation
    boundary conditions, the Green's function is given by

    .. math::

        G(\mathbf{x}, \mathbf{x}^\prime) =
        \frac{1}{4i} H^{(1)}_0(k|\mathbf{x} - \mathbf{x^\prime}|)

    where :math:`k` is the wavenumber and :math:`H^{(1)}_0` is the Hankel function
    of the first kind and of order zero. In the following example, we will solve the
    Helmholtz equation with two Dirac sources.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.special import hankel1
    >>> k = 10  # Wavenumber
    >>> sources = [0.5, -0.5]  # Location of two point sources
    >>> x, y = np.linspace(-3, 3, 300), np.linspace(-3, 3, 300)
    >>> Z = np.add.outer(1j*y, x)
    >>> U = np.zeros_like(Z)
    >>> for sz in sources:
    ...     r = np.abs(Z - sz)
    ...     U += (1j/4)*hankel1(0, k*r)

    Finally, we will plot the real part of the solution.

    >>> fig, ax = plt.subplots()
    >>> ax.pcolormesh(np.real(Z), np.imag(Z), np.real(U))
    >>> plt.show()
    )";

const char *hankel1e_doc = R"(
    hankel1e(v, z, out=None)

    Exponentially scaled Hankel function of the first kind.

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

    Hankel function of the second kind.

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

    Exponentially scaled Hankel function of the second kind.

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

    Gauss hypergeometric function 2F1(a, b; c; z).

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

const char *hyperu_doc = R"(
    hyperu(a, b, x, out=None)

    Confluent hypergeometric function U.

    It is defined as the solution to the equation

    .. math::

       x \frac{d^2w}{dx^2} + (b - x) \frac{dw}{dx} - aw = 0

    which satisfies the property

    .. math::

       U(a, b, x) \sim x^{-a}

    as :math:`x \to \infty`. See [DLMF]_ for more details.

    Parameters
    ----------
    a, b : array_like
        Real-valued parameters
    x : array_like
        Real-valued argument
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Values of `U`

    References
    ----------
    .. [DLMF] NIST Digital Library of Mathematics Functions
              https://dlmf.nist.gov/13.2#E6

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It has a branch cut along the negative `x` axis.

    >>> x = np.linspace(-0.1, -10, 5)
    >>> sc.hyperu(1, 1, x)
    array([nan, nan, nan, nan, nan])

    It approaches zero as `x` goes to infinity.

    >>> x = np.array([1, 10, 100])
    >>> sc.hyperu(1, 1, x)
    array([0.59634736, 0.09156333, 0.00990194])

    It satisfies Kummer's transformation.

    >>> a, b, x = 2, 1, 1
    >>> sc.hyperu(a, b, x)
    0.1926947246463881
    >>> x**(1 - b) * sc.hyperu(a - b + 1, 2 - b, x)
    0.1926947246463881
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

    Integrals of Airy functions.

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

const char *i0_doc = R"(
    i0(x, out=None)

    Modified Bessel function of order 0.

    Defined as,

    .. math::
        I_0(x) = \sum_{k=0}^\infty \frac{(x^2/4)^k}{(k!)^2} = J_0(\imath x),

    where :math:`J_0` is the Bessel function of the first kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float)
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    I : scalar or ndarray
        Value of the modified Bessel function of order 0 at `x`.

    See Also
    --------
    iv: Modified Bessel function of any order
    i0e: Exponentially scaled modified Bessel function of order 0

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `i0`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import i0
    >>> i0(1.)
    1.2660658777520082

    Calculate at several points:

    >>> import numpy as np
    >>> i0(np.array([-2., 0., 3.5]))
    array([2.2795853 , 1.        , 7.37820343])

    Plot the function from -10 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-10., 10., 1000)
    >>> y = i0(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *i0e_doc = R"(
    i0e(x, out=None)

    Exponentially scaled modified Bessel function of order 0.

    Defined as::

        i0e(x) = exp(-abs(x)) * i0(x).

    Parameters
    ----------
    x : array_like
        Argument (float)
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    I : scalar or ndarray
        Value of the exponentially scaled modified Bessel function of order 0
        at `x`.

    See Also
    --------
    iv: Modified Bessel function of the first kind
    i0: Modified Bessel function of order 0

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval. The
    polynomial expansions used are the same as those in `i0`, but
    they are not multiplied by the dominant exponential factor.

    This function is a wrapper for the Cephes [1]_ routine `i0e`. `i0e`
    is useful for large arguments `x`: for these, `i0` quickly overflows.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    In the following example `i0` returns infinity whereas `i0e` still returns
    a finite number.

    >>> from scipy.special import i0, i0e
    >>> i0(1000.), i0e(1000.)
    (inf, 0.012617240455891257)

    Calculate the function at several points by providing a NumPy array or
    list for `x`:

    >>> import numpy as np
    >>> i0e(np.array([-2., 0., 3.]))
    array([0.30850832, 1.        , 0.24300035])

    Plot the function from -10 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-10., 10., 1000)
    >>> y = i0e(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *i1_doc = R"(
    i1(x, out=None)

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
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    I : scalar or ndarray
        Value of the modified Bessel function of order 1 at `x`.

    See Also
    --------
    iv: Modified Bessel function of the first kind
    i1e: Exponentially scaled modified Bessel function of order 1

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `i1`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import i1
    >>> i1(1.)
    0.5651591039924851

    Calculate the function at several points:

    >>> import numpy as np
    >>> i1(np.array([-2., 0., 6.]))
    array([-1.59063685,  0.        , 61.34193678])

    Plot the function between -10 and 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-10., 10., 1000)
    >>> y = i1(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *i1e_doc = R"(
    i1e(x, out=None)

    Exponentially scaled modified Bessel function of order 1.

    Defined as::

        i1e(x) = exp(-abs(x)) * i1(x)

    Parameters
    ----------
    x : array_like
        Argument (float)
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    I : scalar or ndarray
        Value of the exponentially scaled modified Bessel function of order 1
        at `x`.

    See Also
    --------
    iv: Modified Bessel function of the first kind
    i1: Modified Bessel function of order 1

    Notes
    -----
    The range is partitioned into the two intervals [0, 8] and (8, infinity).
    Chebyshev polynomial expansions are employed in each interval. The
    polynomial expansions used are the same as those in `i1`, but
    they are not multiplied by the dominant exponential factor.

    This function is a wrapper for the Cephes [1]_ routine `i1e`. `i1e`
    is useful for large arguments `x`: for these, `i1` quickly overflows.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    In the following example `i1` returns infinity whereas `i1e` still returns
    a finite number.

    >>> from scipy.special import i1, i1e
    >>> i1(1000.), i1e(1000.)
    (inf, 0.01261093025692863)

    Calculate the function at several points by providing a NumPy array or
    list for `x`:

    >>> import numpy as np
    >>> i1e(np.array([-2., 0., 6.]))
    array([-0.21526929,  0.        ,  0.15205146])

    Plot the function between -10 and 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-10., 10., 1000)
    >>> y = i1e(x)
    >>> ax.plot(x, y)
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
    To calculate the orders 0 and 1 for a 1D array:

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

    Internal function.

    Return `iv(v, x) / iv(v-1, x)` for `v >= 0.5` and `x >= 0`, where `iv`
    is the modified Bessel function of the first kind.

    Parameters
    ----------
    v : array_like of float
        Order.  Must be `>= 0.5`.  May be `+inf` if `x` is finite.
    x : array_like of float
        Argument.  Must be `>= 0`.  May be `+inf` if `v` is finite.
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        Ratio between modified Bessel function of the first kind of adjacent
        orders.  The returned value is between `0` and `1`, inclusive.

        If either `v` or `x` is `nan`, `nan` is returned.  Otherwise, the
        special values are:

        - If `v < 0.5` or `x < 0`, set "domain" error and return `nan`.
        - If `v >= 0.5` and `x == 0`, return `x`.
        - If `v >= 0.5` and `x == +inf`, return `1.0`.
        - If `v == +inf` and `0 < x < +inf`, return `0.0`.
        - If `v == +inf` and `x == +inf`, set "domain" error and return `nan`.

    See Also
    --------
    iv : modified Bessel function of the first kind

    Notes
    -----
    The function is computed using the _Perron continued fraction_ of [1]_.
    The continued fraction is evaluated using the "series method" of [2]_.
    Kahan summation is used to evaluate the series.

    The accuracy is tested numerically with 600,000 trials.  The peak
    relative error is `3.4e-16`; the RMSE is `0.9e-16`.

    Reference
    ---------
    .. [1] Gautschi, W. and Slavik, J. (1978). "On the computation of
           modified Bessel function ratios." Mathematics of Computation,
           32(143):865-875.

    .. [2] Gautschi, W. (1967). “Computational Aspects of Three-Term
           Recurrence Relations.” SIAM Review, 9(1):24-82.

    )";

const char *iv_ratio_c_doc = R"(
    _iv_ratio_c(v, x, out=None)

    Internal function.

    Return `1 - iv(v, x) / iv(v-1, x)` for `v >= 0.5` and `x >= 0`, where
    `iv` is the modified Bessel function of the first kind.

    Notes
    -----
    See `_iv_ratio` for details about the parameters, return value, and
    algorithm.

    The accuracy is tested numerically with 600,000 trials.  The peak
    relative error is `9.0e-16`; the RMSE is `1.5e-16`.

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

const char *j0_doc = R"(
    j0(x, out=None)

    Bessel function of the first kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    J : scalar or ndarray
        Value of the Bessel function of the first kind of order 0 at `x`.

    See Also
    --------
    jv : Bessel function of real order and complex argument.
    spherical_jn : spherical Bessel functions.

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
    It should not be confused with the spherical Bessel functions (see
    `spherical_jn`).

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import j0
    >>> j0(1.)
    0.7651976865579665

    Calculate the function at several points:

    >>> import numpy as np
    >>> j0(np.array([-2., 0., 4.]))
    array([ 0.22389078,  1.        , -0.39714981])

    Plot the function from -20 to 20.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-20., 20., 1000)
    >>> y = j0(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *j1_doc = R"(
    j1(x, out=None)

    Bessel function of the first kind of order 1.

    Parameters
    ----------
    x : array_like
        Argument (float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    J : scalar or ndarray
        Value of the Bessel function of the first kind of order 1 at `x`.

    See Also
    --------
    jv: Bessel function of the first kind
    spherical_jn: spherical Bessel functions.

    Notes
    -----
    The domain is divided into the intervals [0, 8] and (8, infinity). In the
    first interval a 24 term Chebyshev expansion is used. In the second, the
    asymptotic trigonometric representation is employed using two rational
    functions of degree 5/5.

    This function is a wrapper for the Cephes [1]_ routine `j1`.
    It should not be confused with the spherical Bessel functions (see
    `spherical_jn`).

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import j1
    >>> j1(1.)
    0.44005058574493355

    Calculate the function at several points:

    >>> import numpy as np
    >>> j1(np.array([-2., 0., 4.]))
    array([-0.57672481,  0.        , -0.06604333])

    Plot the function from -20 to 20.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-20., 20., 1000)
    >>> y = j1(x)
    >>> ax.plot(x, y)
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
    To calculate the orders 0 and 1 for a 1D array:

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

    Kelvin functions as complex numbers.

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

const char *k0_doc = R"(
    k0(x, out=None)

    Modified Bessel function of the second kind of order 0, :math:`K_0`.

    This function is also sometimes referred to as the modified Bessel
    function of the third kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float).
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    K : scalar or ndarray
        Value of the modified Bessel function :math:`K_0` at `x`.

    See Also
    --------
    kv: Modified Bessel function of the second kind of any order
    k0e: Exponentially scaled modified Bessel function of the second kind

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k0`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import k0
    >>> k0(1.)
    0.42102443824070823

    Calculate the function at several points:

    >>> import numpy as np
    >>> k0(np.array([0.5, 2., 3.]))
    array([0.92441907, 0.11389387, 0.0347395 ])

    Plot the function from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> y = k0(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *k0e_doc = R"(
    k0e(x, out=None)

    Exponentially scaled modified Bessel function K of order 0.

    Defined as::

        k0e(x) = exp(x) * k0(x).

    Parameters
    ----------
    x : array_like
        Argument (float)
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    K : scalar or ndarray
        Value of the exponentially scaled modified Bessel function K of order
        0 at `x`.

    See Also
    --------
    kv: Modified Bessel function of the second kind of any order
    k0: Modified Bessel function of the second kind

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k0e`. `k0e` is
    useful for large arguments: for these, `k0` easily underflows.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    In the following example `k0` returns 0 whereas `k0e` still returns a
    useful finite number:

    >>> from scipy.special import k0, k0e
    >>> k0(1000.), k0e(1000)
    (0., 0.03962832160075422)

    Calculate the function at several points by providing a NumPy array or
    list for `x`:

    >>> import numpy as np
    >>> k0e(np.array([0.5, 2., 3.]))
    array([1.52410939, 0.84156822, 0.6977616 ])

    Plot the function from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> y = k0e(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *k1_doc = R"(
    k1(x, out=None)

    Modified Bessel function of the second kind of order 1, :math:`K_1(x)`.

    Parameters
    ----------
    x : array_like
        Argument (float)
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    K : scalar or ndarray
        Value of the modified Bessel function K of order 1 at `x`.

    See Also
    --------
    kv: Modified Bessel function of the second kind of any order
    k1e: Exponentially scaled modified Bessel function K of order 1

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k1`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import k1
    >>> k1(1.)
    0.6019072301972346

    Calculate the function at several points:

    >>> import numpy as np
    >>> k1(np.array([0.5, 2., 3.]))
    array([1.65644112, 0.13986588, 0.04015643])

    Plot the function from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> y = k1(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *k1e_doc = R"(
    k1e(x, out=None)

    Exponentially scaled modified Bessel function K of order 1.

    Defined as::

        k1e(x) = exp(x) * k1(x)

    Parameters
    ----------
    x : array_like
        Argument (float)
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    K : scalar or ndarray
        Value of the exponentially scaled modified Bessel function K of order
        1 at `x`.

    See Also
    --------
    kv: Modified Bessel function of the second kind of any order
    k1: Modified Bessel function of the second kind of order 1

    Notes
    -----
    The range is partitioned into the two intervals [0, 2] and (2, infinity).
    Chebyshev polynomial expansions are employed in each interval.

    This function is a wrapper for the Cephes [1]_ routine `k1e`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    In the following example `k1` returns 0 whereas `k1e` still returns a
    useful floating point number.

    >>> from scipy.special import k1, k1e
    >>> k1(1000.), k1e(1000.)
    (0., 0.03964813081296021)

    Calculate the function at several points by providing a NumPy array or
    list for `x`:

    >>> import numpy as np
    >>> k1e(np.array([0.5, 2., 3.]))
    array([2.73100971, 1.03347685, 0.80656348])

    Plot the function from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> y = k1e(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *kv_doc = R"(
    kv(v, z, out=None)

    Modified Bessel function of the second kind of real order `v`.

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

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import loggamma, gamma

    >>> z = 1.5 + 2j
    >>> loggamma(z)
    np.complex128(-1.4991963725850939+0.7332806816909994j)

    Verify :math:`\exp(\log \Gamma(z)) = \Gamma(z)`:

    >>> np.exp(loggamma(z))
    np.complex128(0.165915108938991+0.14946347326641998j)
    >>> gamma(z)
    np.complex128(0.165915108938991+0.14946347326641998j)

    Verify the recurrence :math:`\log \Gamma(z+1) = \log(z) + \log \Gamma(z)`:

    >>> loggamma(z + 1)
    np.complex128(-0.5829056407109388+1.6605758996926108j)
    >>> np.log(z) + loggamma(z)
    np.complex128(-0.5829056407109388+1.6605758996926117j)
    )";

const char *logit_doc = R"(
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

const char *log_ndtr_doc = R"(
    log_ndtr(x, out=None)

    Logarithm of Gaussian cumulative distribution function.

    Returns the log of the area under the standard Gaussian probability
    density function, integrated from minus infinity to `x`::

        log(1/sqrt(2*pi) * integral(exp(-t**2 / 2), t=-inf..x))

    Parameters
    ----------
    x : array_like, real or complex
        Argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The value of the log of the normal CDF evaluated at `x`

    See Also
    --------
    erf
    erfc
    scipy.stats.norm
    ndtr

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import log_ndtr, ndtr

    The benefit of ``log_ndtr(x)`` over the naive implementation
    ``np.log(ndtr(x))`` is most evident with moderate to large positive
    values of ``x``:

    >>> x = np.array([6, 7, 9, 12, 15, 25])
    >>> log_ndtr(x)
    array([-9.86587646e-010, -1.27981254e-012, -1.12858841e-019,
           -1.77648211e-033, -3.67096620e-051, -3.05669671e-138])

    The results of the naive calculation for the moderate ``x`` values
    have only 5 or 6 correct significant digits. For values of ``x``
    greater than approximately 8.3, the naive expression returns 0:

    >>> np.log(ndtr(x))
    array([-9.86587701e-10, -1.27986510e-12,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00])
    )";

const char *ndtri_doc = R"(
    ndtri(p, out=None)

    Inverse of `ndtr`.

    Returns the quantile `x` such that the cumulative distribution function of the
    standard normal distribution evaluated at `x` equals `p`, that is, ``ndtr(x) == p``.

    Parameters
    ----------
    p : array_like
        Probability values.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    x : scalar or ndarray
        Quantile(s) corresponding to the probabilitie(s) in `p`.

    See Also
    --------
    ndtr : Standard normal cumulative distribution function
    ndtri_exp : Inverse of log_ndtr

    Examples
    --------
    `ndtri` is the percentile (quantile) function of the standard normal distribution,
    i.e., the inverse of the cumulative distribution function `ndtr`.

    First, compute a cumulative distribution value:

    >>> import numpy as np
    >>> from scipy.special import ndtri, ndtr
    >>> cdf_val = ndtr(2)
    >>> cdf_val
    0.9772498680518208

    Verify that `ndtri` yields the original value for `x` up to floating point errors.

    >>> ndtri(cdf_val)
    2.0000000000000004

    Plot the percentile function over a range of probabilities.

    >>> import matplotlib.pyplot as plt
    >>> p = np.linspace(1e-3, 1 - 1e-3, 201)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(p, ndtri(p))
    >>> ax.set_title("Standard normal percentile function")
    >>> plt.show()
    )";

const char *ndtri_exp_doc = R"(
    ndtri_exp(y, out=None)

    Inverse of `log_ndtr` vs x. Allows for greater precision than
    `ndtri` composed with `numpy.exp` for very small values of y and for
    y close to 0.

    Parameters
    ----------
    y : array_like of float
        Function argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Inverse of the log CDF of the standard normal distribution, evaluated
        at y.

    See Also
    --------
    log_ndtr : log of the standard normal cumulative distribution function
    ndtr : standard normal cumulative distribution function
    ndtri : standard normal percentile function

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    `ndtri_exp` agrees with the naive implementation when the latter does
    not suffer from underflow.

    >>> sc.ndtri_exp(-1)
    -0.33747496376420244
    >>> sc.ndtri(np.exp(-1))
    -0.33747496376420244

    For extreme values of y, the naive approach fails

    >>> sc.ndtri(np.exp(-800))
    -inf
    >>> sc.ndtri(np.exp(-1e-20))
    inf

    whereas `ndtri_exp` is still able to compute the result to high precision.

    >>> sc.ndtri_exp(-800)
    -39.88469483825668
    >>> sc.ndtri_exp(-1e-20)
    9.262340089798409
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

    Characteristic value of even Mathieu functions.

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

    Characteristic value of odd Mathieu functions.

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

    Even Mathieu function and its derivative.

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

    Notes
    -----
    The even Mathieu functions are the solutions to Mathieu's differential equation

    .. math::

        \frac{d^2y}{dx^2} + (a_m - 2q \cos(2x))y = 0

    for which the characteristic number :math:`a_m` (calculated with `mathieu_a`)
    results in an odd, periodic solution :math:`y(x)` with period 180 degrees
    (for even :math:`m`) or 360 degrees (for odd :math:`m`).

    References
    ----------
    .. [1] 'Mathieu function'. *Wikipedia*.
           https://en.wikipedia.org/wiki/Mathieu_function

    Examples
    --------
    Plot even Mathieu functions of orders ``2`` and ``4``.

    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> m = np.asarray([2, 4])
    >>> q = 50
    >>> x = np.linspace(-180, 180, 300)[:, np.newaxis]
    >>> y, _ = special.mathieu_cem(m, q, x)
    >>> plt.plot(x, y)
    >>> plt.xlabel('x (degrees)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 2', 'm = 4'))

    Because the orders ``2`` and
    ``4`` are even, the period of each function is 180 degrees.

    )";

const char *mathieu_modcem1_doc = R"(
    mathieu_modcem1(m, q, x, out=None)

    Even modified Mathieu function of the first kind and its derivative.

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

    Even modified Mathieu function of the second kind and its derivative.

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

    Odd modified Mathieu function of the first kind and its derivative.

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

    Odd modified Mathieu function of the second kind and its derivative.

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

    Odd Mathieu function and its derivative.

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

    Notes
    -----
    Odd Mathieu functions are the solutions to Mathieu's differential equation

    .. math::

        \frac{d^2y}{dx^2} + (b_m - 2q \cos(2x))y = 0

    for which the characteristic number :math:`b_m` (calculated with `mathieu_b`)
    results in an odd, periodic solution :math:`y(x)` with period 180 degrees
    (for even :math:`m`) or 360 degrees (for odd :math:`m`).

    References
    ----------
    .. [1] 'Mathieu function'. *Wikipedia*.
           https://en.wikipedia.org/wiki/Mathieu_function

    Examples
    --------
    Plot odd Mathieu functions of orders ``2`` and ``4``.

    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt
    >>> m = np.asarray([2, 4])
    >>> q = 50
    >>> x = np.linspace(-180, 180, 300)[:, np.newaxis]
    >>> y, _ = special.mathieu_sem(m, q, x)
    >>> plt.plot(x, y)
    >>> plt.xlabel('x (degrees)')
    >>> plt.ylabel('y')
    >>> plt.legend(('m = 2', 'm = 4'))

    Because the orders ``2`` and
    ``4`` are even, the period of each function is 180 degrees.

    )";

const char *modfresnelm_doc = R"(
    modfresnelm(x, out=None)

    Modified Fresnel negative integrals.

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
        Integral ``K_-(x)``: ``1/sqrt(pi)*exp(1j*(x*x+pi/4))*fm``

    See Also
    --------
    modfresnelp

    )";

const char *modfresnelp_doc = R"(
    modfresnelp(x, out=None)

    Modified Fresnel positive integrals.

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

    Oblate spheroidal angular function of the first kind and its derivative.

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

const char *ndtr_doc = R"(
    ndtr(x, out=None)

    Cumulative distribution of the standard normal distribution.

    Returns the area under the standard Gaussian probability
    density function, integrated from minus infinity to `x`

    .. math::

       \frac{1}{\sqrt{2\pi}} \int_{-\infty}^x \exp(-t^2/2) dt

    Parameters
    ----------
    x : array_like, real or complex
        Argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        The value of the normal CDF evaluated at `x`

    See Also
    --------
    log_ndtr : Logarithm of ndtr
    ndtri : Inverse of ndtr, standard normal percentile function
    erf : Error function
    erfc : 1 - erf
    scipy.stats.norm : Normal distribution

    Examples
    --------
    Evaluate `ndtr` at one point.

    >>> import numpy as np
    >>> from scipy.special import ndtr
    >>> ndtr(0.5)
    0.6914624612740131

    Evaluate the function at several points by providing a NumPy array
    or list for `x`.

    >>> ndtr([0, 0.5, 2])
    array([0.5       , 0.69146246, 0.97724987])

    Plot the function.

    >>> import matplotlib.pyplot as plt
    >>> x = np.linspace(-5, 5, 100)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, ndtr(x))
    >>> ax.set_title(r"Standard normal cumulative distribution function $\Phi$")
    >>> plt.show()
    )";

const char *obl_ang1_cv_doc = R"(
    obl_ang1_cv(m, n, c, cv, x, out=None)

    Oblate spheroidal angular function obl_ang1 for precomputed characteristic value.

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

    Characteristic value of oblate spheroidal function.

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

    Oblate spheroidal radial function of the first kind and its derivative.

    Computes the oblate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``x >= 0.0``.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Parameter x (``x >= 0.0``)
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

    Oblate spheroidal radial function obl_rad1 for precomputed characteristic value.

    Computes the oblate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``x >= 0.0``. Requires
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
        Parameter x (``x >= 0.0``)
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
    and n>=m, spheroidal parameter `c` and ``x >= 0.0``.

    Parameters
    ----------
    m : array_like
        Mode parameter m (nonnegative)
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Parameter x (``x >= 0.0``)
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

    Oblate spheroidal radial function obl_rad2 for precomputed characteristic value.

    Computes the oblate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``x >= 0.0``. Requires
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
        Parameter x (``x >= 0.0``)
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

    Parabolic cylinder function D.

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

    Parabolic cylinder function V.

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

    Prolate spheroidal angular function of the first kind and its derivative.

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

    Prolate spheroidal angular function pro_ang1 for precomputed characteristic value.

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

    Characteristic value of prolate spheroidal function.

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

    Prolate spheroidal radial function of the first kind and its derivative.

    Computes the prolate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``x > 1.0``.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Real parameter (``x > 1.0``)
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

    Prolate spheroidal radial function pro_rad1 for precomputed characteristic value.

    Computes the prolate spheroidal radial function of the first kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``x > 1.0``. Requires
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
        Real parameter (``x > 1.0``)
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

    Prolate spheroidal radial function of the second kind and its derivative.

    Computes the prolate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``x > 1.0``.

    Parameters
    ----------
    m : array_like
        Nonnegative mode parameter m
    n : array_like
        Mode parameter n (>= m)
    c : array_like
        Spheroidal parameter
    x : array_like
        Real parameter (``x > 1.0``)
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

    Prolate spheroidal radial function pro_rad2 for precomputed characteristic value.

    Computes the prolate spheroidal radial function of the second kind
    and its derivative (with respect to `x`) for mode parameters m>=0
    and n>=m, spheroidal parameter `c` and ``x > 1.0``. Requires
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
        Real parameter (``x > 1.0``)
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

    .. math:: \psi(z) = \frac{d}{dz}\log\Gamma(z) = \frac{\Gamma'(z)}{\Gamma(z)}

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

const char *digammainv_doc = R"(
    digammainv(y, out=None)

    Inverse of the digamma function.

    Parameters
    ----------
    y : array_like of float
        Real-valued argument.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Positive real solution of :math:`\psi(x) = y`.

    See Also
    --------
    psi : The digamma function.
    digamma : Alias for `psi`.

    Notes
    -----
    .. versionadded:: 1.19.0

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.special import digamma, digammainv
    >>> x = np.logspace(-2, 2, 5)
    >>> np.allclose(digammainv(digamma(x)), x)
    True

    Plot the digamma function and its inverse:

    >>> x = np.linspace(0.3, 5, 200)
    >>> y = digamma(x)
    >>> t = np.linspace(-3, 4, 200)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, y, label=r"$\psi(x)$")
    >>> ax.plot(y, digammainv(y), label=r"$\psi^{-1}(y)$")
    >>> ax.plot(t, t, "--", label=r"$y = x$")
    >>> ax.set_xlabel("input")
    >>> ax.legend()
    >>> plt.show()
    )";

const char *radian_doc = R"(
    radian(d, m, s, out=None)

    Convert from degrees to radians.

    Returns the angle given in (d)egrees, (m)inutes, and (s)econds in
    radians.

    Parameters
    ----------
    d : array_like
        Degrees, can be real-valued.
    m : array_like
        Minutes, can be real-valued.
    s : array_like
        Seconds, can be real-valued.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Values of the inputs in radians.

    Examples
    --------
    >>> import scipy.special as sc

    There are many ways to specify an angle.

    >>> sc.radian(90, 0, 0)
    1.5707963267948966
    >>> sc.radian(0, 60 * 90, 0)
    1.5707963267948966
    >>> sc.radian(0, 0, 60**2 * 90)
    1.5707963267948966

    The inputs can be real-valued.

    >>> sc.radian(1.5, 0, 0)
    0.02617993877991494
    >>> sc.radian(1, 30, 0)
    0.02617993877991494
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

const char *_riemann_zeta_doc = R"(
    Internal function, use `zeta` instead.
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

const char *tandg_doc = R"(
    tandg(x, out=None)

    Tangent of angle `x` given in degrees.

    Parameters
    ----------
    x : array_like
        Angle, given in degrees.
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    scalar or ndarray
        Tangent at the input.

    See Also
    --------
    sindg, cosdg, cotdg

    Examples
    --------
    >>> import numpy as np
    >>> import scipy.special as sc

    It is more accurate than using tangent directly.

    >>> x = 180 * np.arange(3)
    >>> sc.tandg(x)
    array([0., 0., 0.])
    >>> np.tan(x * np.pi / 180)
    array([ 0.0000000e+00, -1.2246468e-16, -2.4492936e-16])
    )";

const char *tklmbda_doc = R"(
    tklmbda(x, lmbda, out=None)

    Cumulative distribution function of the Tukey lambda distribution.

    Parameters
    ----------
    x, lmbda : array_like
        Parameters
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    cdf : scalar or ndarray
        Value of the Tukey lambda CDF

    See Also
    --------
    scipy.stats.tukeylambda : Tukey lambda distribution

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.special import tklmbda, expit

    Compute the cumulative distribution function (CDF) of the Tukey lambda
    distribution at several ``x`` values for `lmbda` = -1.5.

    >>> x = np.linspace(-2, 2, 9)
    >>> x
    array([-2. , -1.5, -1. , -0.5,  0. ,  0.5,  1. ,  1.5,  2. ])
    >>> tklmbda(x, -1.5)
    array([0.34688734, 0.3786554 , 0.41528805, 0.45629737, 0.5       ,
           0.54370263, 0.58471195, 0.6213446 , 0.65311266])

    When `lmbda` is 0, the function is the logistic sigmoid function,
    which is implemented in `scipy.special` as `expit`.

    >>> tklmbda(x, 0)
    array([0.11920292, 0.18242552, 0.26894142, 0.37754067, 0.5       ,
           0.62245933, 0.73105858, 0.81757448, 0.88079708])
    >>> expit(x)
    array([0.11920292, 0.18242552, 0.26894142, 0.37754067, 0.5       ,
           0.62245933, 0.73105858, 0.81757448, 0.88079708])

    When `lmbda` is 1, the Tukey lambda distribution is uniform on the
    interval [-1, 1], so the CDF increases linearly.

    >>> t = np.linspace(-1, 1, 9)
    >>> tklmbda(t, 1)
    array([0.   , 0.125, 0.25 , 0.375, 0.5  , 0.625, 0.75 , 0.875, 1.   ])

    In the following, we generate plots for several values of `lmbda`.

    The first figure shows graphs for `lmbda` <= 0.

    >>> styles = ['-', '-.', '--', ':']
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-12, 12, 500)
    >>> for k, lmbda in enumerate([-1.0, -0.5, 0.0]):
    ...     y = tklmbda(x, lmbda)
    ...     ax.plot(x, y, styles[k], label=rf'$\lambda$ = {lmbda:-4.1f}')

    >>> ax.set_title(r'tklmbda(x, $\lambda$)')
    >>> ax.set_label('x')
    >>> ax.legend(framealpha=1, shadow=True)
    >>> ax.grid(True)

    The second figure shows graphs for `lmbda` > 0.  The dots in the
    graphs show the bounds of the support of the distribution.

    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-4.2, 4.2, 500)
    >>> lmbdas = [0.25, 0.5, 1.0, 1.5]
    >>> for k, lmbda in enumerate(lmbdas):
    ...     y = tklmbda(x, lmbda)
    ...     ax.plot(x, y, styles[k], label=fr'$\lambda$ = {lmbda}')

    >>> ax.set_prop_cycle(None)
    >>> for lmbda in lmbdas:
    ...     ax.plot([-1/lmbda, 1/lmbda], [0, 1], '.', ms=8)

    >>> ax.set_title(r'tklmbda(x, $\lambda$)')
    >>> ax.set_xlabel('x')
    >>> ax.legend(framealpha=1, shadow=True)
    >>> ax.grid(True)

    >>> plt.tight_layout()
    >>> plt.show()

    The CDF of the Tukey lambda distribution is also implemented as the
    ``cdf`` method of `scipy.stats.tukeylambda`.  In the following,
    ``tukeylambda.cdf(x, -0.5)`` and ``tklmbda(x, -0.5)`` compute the
    same values:

    >>> from scipy.stats import tukeylambda
    >>> x = np.linspace(-2, 2, 9)

    >>> tukeylambda.cdf(x, -0.5)
    array([0.21995157, 0.27093858, 0.33541677, 0.41328161, 0.5       ,
           0.58671839, 0.66458323, 0.72906142, 0.78004843])

    >>> tklmbda(x, -0.5)
    array([0.21995157, 0.27093858, 0.33541677, 0.41328161, 0.5       ,
           0.58671839, 0.66458323, 0.72906142, 0.78004843])

    The implementation in ``tukeylambda`` also provides location and scale
    parameters, and other methods such as ``pdf()`` (the probability
    density function) and ``ppf()`` (the inverse of the CDF), so for
    working with the Tukey lambda distribution, ``tukeylambda`` is more
    generally useful.  The primary advantage of ``tklmbda`` is that it is
    significantly faster than ``tukeylambda.cdf``.
    )";

const char *spence_doc = R"(
    spence(z, out=None)

    Spence's function, also known as the dilogarithm.

    It is defined to be

    .. math::
      \int_1^z \frac{\log(t)}{1 - t}dt

    for complex :math:`z`, where the contour of integration is taken
    to avoid the branch cut of the logarithm. Spence's function is
    analytic everywhere except the negative real axis where it has a
    branch cut.

    Parameters
    ----------
    z : array_like
        Points at which to evaluate Spence's function
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    s : scalar or ndarray
        Computed values of Spence's function

    Notes
    -----
    There is a different convention which defines Spence's function by
    the integral

    .. math::
      -\int_0^z \frac{\log(1 - t)}{t}dt,

    which is equivalent to ``spence(1 - z)``.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.special import spence
    >>> import matplotlib.pyplot as plt

    The function is defined for complex inputs:

    >>> spence([1-1j, 1.5+2j, 3j, -10-5j])
    array([-0.20561676+0.91596559j, -0.86766909-1.39560134j,
           -0.59422064-2.49129918j, -1.14044398+6.80075924j])

    For complex inputs on the branch cut, which is the negative real axis,
    the function returns the limit for ``z`` with positive imaginary part.
    For example, in the following, note the sign change of the imaginary
    part of the output for ``z = -2`` and ``z = -2 - 1e-8j``:

    >>> spence([-2 + 1e-8j, -2, -2 - 1e-8j])
    array([2.32018041-3.45139229j, 2.32018042-3.4513923j ,
           2.32018041+3.45139229j])

    The function returns ``nan`` for real inputs on the branch cut:

    >>> spence(-1.5)
    nan

    Verify some particular values: ``spence(0) = pi**2/6``,
    ``spence(1) = 0`` and ``spence(2) = -pi**2/12``.

    >>> spence([0, 1, 2])
    array([ 1.64493407,  0.        , -0.82246703])
    >>> np.pi**2/6, -np.pi**2/12
    (1.6449340668482264, -0.8224670334241132)

    Verify the identity::

        spence(z) + spence(1 - z) = pi**2/6 - log(z)*log(1 - z)

    >>> z = 3 + 4j
    >>> spence(z) + spence(1 - z)
    (-2.6523186143876067+1.8853470951513935j)
    >>> np.pi**2/6 - np.log(z)*np.log(1 - z)
    (-2.652318614387606+1.885347095151394j)

    Plot the function for positive real input.

    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0, 6, 400)
    >>> ax.plot(x, spence(x))
    >>> ax.grid()
    >>> ax.set_xlabel('x')
    >>> ax.set_title('spence(x)')
    >>> plt.show()
)";

const char *struve_h_doc = R"(
    struve(v, x, out=None)

    Struve function.

    Return the value of the Struve function of order `v` at `x`.  The Struve
    function is defined as,

    .. math::
        H_v(x) = (x/2)^{v + 1} \sum_{n=0}^\infty
        \frac{(-1)^n (x/2)^{2n}}{\Gamma(n + \frac{3}{2}) \Gamma(n + v + \frac{3}{2})},

    where :math:`\Gamma` is the gamma function.

    Parameters
    ----------
    v : array_like
        Order of the Struve function (float).
    x : array_like
        Argument of the Struve function (float; must be positive unless `v` is
        an integer).
    out : ndarray, optional
        Optional output array for the function results.

    Returns
    -------
    H : scalar or ndarray
        Value of the Struve function of order `v` at `x`.

    See Also
    --------
    modstruve : Modified Struve function

    Notes
    -----
    Three methods discussed in [1]_ are used to evaluate the Struve function:

    - power series
    - expansion in Bessel functions (if :math:`|x| < |v| + 20`)
    - asymptotic large-x expansion (if :math:`x \geq 0.7v + 12`)

    Rounding errors are estimated based on the largest terms in the sums, and
    the result associated with the smallest error is returned.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/11

    Examples
    --------
    Calculate the Struve function of order 1 at 2.

    >>> import numpy as np
    >>> from scipy.special import struve
    >>> import matplotlib.pyplot as plt
    >>> struve(1, 2.)
    0.6467637282835622

    Calculate the Struve function at 2 for orders 1, 2 and 3 by providing
    a list for the order parameter `v`.

    >>> struve([1, 2, 3], 2.)
    array([0.64676373, 0.28031806, 0.08363767])

    Calculate the Struve function of order 1 for several points by providing
    an array for `x`.

    >>> points = np.array([2., 5., 8.])
    >>> struve(1, points)
    array([0.64676373, 0.80781195, 0.48811605])

    Compute the Struve function for several orders at several points by
    providing arrays for `v` and `x`. The arrays have to be broadcastable
    to the correct shapes.

    >>> orders = np.array([[1], [2], [3]])
    >>> points.shape, orders.shape
    ((3,), (3, 1))

    >>> struve(orders, points)
    array([[0.64676373, 0.80781195, 0.48811605],
           [0.28031806, 1.56937455, 1.51769363],
           [0.08363767, 1.50872065, 2.98697513]])

    Plot the Struve functions of order 0 to 3 from -10 to 10.

    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-10., 10., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, struve(i, x), label=f'$H_{i!r}$')
    >>> ax.legend(ncol=2)
    >>> ax.set_xlim(-10, 10)
    >>> ax.set_title(r"Struve functions $H_{\nu}$")
    >>> plt.show()
    )";

const char *struve_l_doc = R"(
    modstruve(v, x, out=None)

    Modified Struve function.

    Return the value of the modified Struve function of order `v` at `x`.  The
    modified Struve function is defined as,

    .. math::
        L_v(x) = -\imath \exp(-\pi\imath v/2) H_v(\imath x),

    where :math:`H_v` is the Struve function.

    Parameters
    ----------
    v : array_like
        Order of the modified Struve function (float).
    x : array_like
        Argument of the Struve function (float; must be positive unless `v` is
        an integer).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    L : scalar or ndarray
        Value of the modified Struve function of order `v` at `x`.

    See Also
    --------
    struve

    Notes
    -----
    Three methods discussed in [1]_ are used to evaluate the function:

    - power series
    - expansion in Bessel functions (if :math:`|x| < |v| + 20`)
    - asymptotic large-x expansion (if :math:`x \geq 0.7v + 12`)

    Rounding errors are estimated based on the largest terms in the sums, and
    the result associated with the smallest error is returned.

    References
    ----------
    .. [1] NIST Digital Library of Mathematical Functions
           https://dlmf.nist.gov/11

    Examples
    --------
    Calculate the modified Struve function of order 1 at 2.

    >>> import numpy as np
    >>> from scipy.special import modstruve
    >>> import matplotlib.pyplot as plt
    >>> modstruve(1, 2.)
    1.102759787367716

    Calculate the modified Struve function at 2 for orders 1, 2 and 3 by
    providing a list for the order parameter `v`.

    >>> modstruve([1, 2, 3], 2.)
    array([1.10275979, 0.41026079, 0.11247294])

    Calculate the modified Struve function of order 1 for several points
    by providing an array for `x`.

    >>> points = np.array([2., 5., 8.])
    >>> modstruve(1, points)
    array([  1.10275979,  23.72821578, 399.24709139])

    Compute the modified Struve function for several orders at several
    points by providing arrays for `v` and `x`. The arrays have to be
    broadcastable to the correct shapes.

    >>> orders = np.array([[1], [2], [3]])
    >>> points.shape, orders.shape
    ((3,), (3, 1))

    >>> modstruve(orders, points)
    array([[1.10275979e+00, 2.37282158e+01, 3.99247091e+02],
           [4.10260789e-01, 1.65535979e+01, 3.25973609e+02],
           [1.12472937e-01, 9.42430454e+00, 2.33544042e+02]])

    Plot the modified Struve functions of order 0 to 3 from -5 to 5.

    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(-5., 5., 1000)
    >>> for i in range(4):
    ...     ax.plot(x, modstruve(i, x), label=f'$L_{i!r}$')
    >>> ax.legend(ncol=2)
    >>> ax.set_xlim(-5, 5)
    >>> ax.set_title(r"Modified Struve functions $L_{\nu}$")
    >>> plt.show()
    )";

const char *voigt_profile_doc = R"(
    voigt_profile(x, sigma, gamma, out=None)

    Voigt profile.

    The Voigt profile is a convolution of a 1-D Normal distribution with
    standard deviation ``sigma`` and a 1-D Cauchy distribution with half-width at
    half-maximum ``gamma``.

    If ``sigma = 0``, PDF of Cauchy distribution is returned.
    Conversely, if ``gamma = 0``, PDF of Normal distribution is returned.
    If ``sigma = gamma = 0``, the return value is ``Inf`` for ``x = 0``,
    and ``0`` for all other ``x``.

    Parameters
    ----------
    x : array_like
        Real argument
    sigma : array_like
        The standard deviation of the Normal distribution part
    gamma : array_like
        The half-width at half-maximum of the Cauchy distribution part
    out : ndarray, optional
        Optional output array for the function values

    Returns
    -------
    scalar or ndarray
        The Voigt profile at the given arguments

    See Also
    --------
    wofz : Faddeeva function

    Notes
    -----
    It can be expressed in terms of Faddeeva function

    .. math:: V(x; \sigma, \gamma) = \frac{Re[w(z)]}{\sigma\sqrt{2\pi}},
    .. math:: z = \frac{x + i\gamma}{\sqrt{2}\sigma}

    where :math:`w(z)` is the Faddeeva function.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Voigt_profile

    Examples
    --------
    Calculate the function at point 2 for ``sigma=1`` and ``gamma=1``.

    >>> from scipy.special import voigt_profile
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> voigt_profile(2, 1., 1.)
    0.09071519942627544

    Calculate the function at several points by providing a NumPy array
    for `x`.

    >>> values = np.array([-2., 0., 5])
    >>> voigt_profile(values, 1., 1.)
    array([0.0907152 , 0.20870928, 0.01388492])

    Plot the function for different parameter sets.

    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> x = np.linspace(-10, 10, 500)
    >>> parameters_list = [(1.5, 0., "solid"), (1.3, 0.5, "dashed"),
    ...                    (0., 1.8, "dotted"), (1., 1., "dashdot")]
    >>> for params in parameters_list:
    ...     sigma, gamma, linestyle = params
    ...     voigt = voigt_profile(x, sigma, gamma)
    ...     ax.plot(x, voigt, label=rf"$\sigma={sigma},\, \gamma={gamma}$",
    ...             ls=linestyle)
    >>> ax.legend()
    >>> plt.show()

    Verify visually that the Voigt profile indeed arises as the convolution
    of a normal and a Cauchy distribution.

    >>> from scipy.signal import convolve
    >>> x, dx = np.linspace(-10, 10, 500, retstep=True)
    >>> def gaussian(x, sigma):
    ...     return np.exp(-0.5 * x**2/sigma**2)/(sigma * np.sqrt(2*np.pi))
    >>> def cauchy(x, gamma):
    ...     return gamma/(np.pi * (np.square(x)+gamma**2))
    >>> sigma = 2
    >>> gamma = 1
    >>> gauss_profile = gaussian(x, sigma)
    >>> cauchy_profile = cauchy(x, gamma)
    >>> convolved = dx * convolve(cauchy_profile, gauss_profile, mode="same")
    >>> voigt = voigt_profile(x, sigma, gamma)
    >>> fig, ax = plt.subplots(figsize=(8, 8))
    >>> ax.plot(x, gauss_profile, label="Gauss: $G$", c='b')
    >>> ax.plot(x, cauchy_profile, label="Cauchy: $C$", c='y', ls="dashed")
    >>> xx = 0.5*(x[1:] + x[:-1])  # midpoints
    >>> ax.plot(xx, convolved[1:], label="Convolution: $G * C$", ls='dashdot',
    ...         c='k')
    >>> ax.plot(x, voigt, label="Voigt", ls='dotted', c='r')
    >>> ax.legend()
    >>> plt.show()
    )";

const char *wofz_doc = R"(
    wofz(z, out=None)

    Faddeeva function.

    Returns the value of the Faddeeva function for complex argument:

    .. math::

        w(z) = e^{-z^2} \operatorname{erfc}(-i z)

    Parameters
    ----------
    z : array_like
        complex argument
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    scalar or ndarray
        Value of the Faddeeva function

    See Also
    --------
    dawsn, erf, erfc, erfcx, erfi

    References
    ----------
    .. [1] Steven G. Johnson, Faddeeva W function implementation.
       http://ab-initio.mit.edu/Faddeeva

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import special
    >>> import matplotlib.pyplot as plt

    >>> x = np.linspace(-3, 3)
    >>> z = special.wofz(x)

    >>> plt.plot(x, z.real, label='wofz(x).real')
    >>> plt.plot(x, z.imag, label='wofz(x).imag')
    >>> plt.xlabel('$x$')
    >>> plt.legend(framealpha=1, shadow=True)
    >>> plt.grid(alpha=0.25)
    >>> plt.show()
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

const char *wrightomega_doc = R"(
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
    )";

const char *y0_doc = R"(
    y0(x, out=None)

    Bessel function of the second kind of order 0.

    Parameters
    ----------
    x : array_like
        Argument (float).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    Y : scalar or ndarray
        Value of the Bessel function of the second kind of order 0 at `x`.

    See Also
    --------
    j0: Bessel function of the first kind of order 0
    yv: Bessel function of the first kind

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

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import y0
    >>> y0(1.)
    0.08825696421567697

    Calculate at several points:

    >>> import numpy as np
    >>> y0(np.array([0.5, 2., 3.]))
    array([-0.44451873,  0.51037567,  0.37685001])

    Plot the function from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> y = y0(x)
    >>> ax.plot(x, y)
    >>> plt.show()
    )";

const char *y1_doc = R"(
    y1(x, out=None)

    Bessel function of the second kind of order 1.

    Parameters
    ----------
    x : array_like
        Argument (float).
    out : ndarray, optional
        Optional output array for the function results

    Returns
    -------
    Y : scalar or ndarray
        Value of the Bessel function of the second kind of order 1 at `x`.

    See Also
    --------
    j1: Bessel function of the first kind of order 1
    yn: Bessel function of the second kind
    yv: Bessel function of the second kind

    Notes
    -----
    The domain is divided into the intervals [0, 8] and (8, infinity). In the
    first interval a 25 term Chebyshev expansion is used, and computing
    :math:`J_1` (the Bessel function of the first kind) is required. In the
    second, the asymptotic trigonometric representation is employed using two
    rational functions of degree 5/5.

    This function is a wrapper for the Cephes [1]_ routine `y1`.

    References
    ----------
    .. [1] Cephes Mathematical Functions Library,
           http://www.netlib.org/cephes/

    Examples
    --------
    Calculate the function at one point:

    >>> from scipy.special import y1
    >>> y1(1.)
    -0.7812128213002888

    Calculate at several points:

    >>> import numpy as np
    >>> y1(np.array([0.5, 2., 3.]))
    array([-1.47147239, -0.10703243,  0.32467442])

    Plot the function from 0 to 10.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x = np.linspace(0., 10., 1000)
    >>> y = y1(x)
    >>> ax.plot(x, y)
    >>> plt.show()
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
    To calculate the orders 0 and 1 for a 1D array:

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
