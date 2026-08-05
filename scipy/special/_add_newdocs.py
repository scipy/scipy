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


add_newdoc("_ellip_harm",
    """
    Internal function, use `ellip_harm` instead.
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

    Examples
    --------
    We have a coin for which the probability of showing heads when flipped
    is 0.525. The coin is flipped 16 times.  What is the probability that
    the number of heads is less than or equal to 5?

    Let X be the number of heads.  We want to find the probability that
    X <= `k`, given `k` is 5, the total number of flips `n` is 16 and
    the probability `p` of heads for a single flip is 0.525.
    This is what ``bdtr(k, n, p)`` computes:

    >>> from scipy.special import bdtr

    >>> bdtr(5, 16, 0.525)
    np.float64(0.07281293895810999)

    The following plot shows the graph of ``bdtr(k, 16, 0.525)``:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

    >>> n = 16
    >>> p = 0.525
    >>> k = np.arange(n + 1)
    >>> plt.plot(k, bdtr(k, n, p), 'o')
    >>> plt.grid(True)
    >>> plt.xlabel('k')
    >>> plt.title(f"bdtr(k, {n}, {p})")
    >>> plt.show()
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

    Examples
    --------
    We have a coin for which the probability of showing heads when flipped
    is 0.525. The coin is flipped 16 times.  What is the probability that
    the number of heads is greater than 10?

    Let X be the number of heads.  We want to find the probability that
    X > `k`, given `k` is 10, the total number of flips `n` is 16 and the
    probability `p` of heads in a single flip is 0.525.  This is what
    ``bdtrc(k, n, p)`` computes:

    >>> from scipy.special import bdtrc

    >>> bdtrc(10, 16, 0.525)
    np.float64(0.14643065267555583)

    The following plot shows the graph ``bdtrc(k, 16, 0.525)``:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt

    >>> n = 16
    >>> p = 0.525
    >>> k = np.arange(n + 1)
    >>> plt.plot(k, bdtrc(k, n, p), 'o')
    >>> plt.grid(True)
    >>> plt.xlabel('k')
    >>> plt.title(f"bdtrc(k, {n}, {p})")
    >>> plt.show()
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

    Examples
    --------
    An "unfair" coin is to be created that has probability `p` of showing
    heads when flipped.  What is the value of `p` that will ensure that
    the probability of getting heads at most once in 4 tosses is 0.5?

    Let X be the number of heads. We want to find `p` such that the probability
    that X <= `k` is `y`, where `k` is 1, the total number of flips `n` is 4,
    and the cumulative probability `y` is 0.5. This is what ``bdtri(k, n, y)``
    computes:

    >>> from scipy.special import bdtri, bdtr

    >>> p = bdtri(1, 4, 0.5)
    >>> p
    np.float64(0.3857275681323896)

    Verify the result:

    >>> bdtr(1, 4, p)  # Should be 0.5.
    np.float64(0.5)
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
    x : array_like
        Points at which to evaluate the polynomial.
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
    `smirnov` is used by `scipy.stats.kstest` in the application of the
    Kolmogorov-Smirnov Goodness of Fit test. For historical reasons this
    function is exposed in `scipy.special`, but the recommended way to achieve
    the most accurate CDF/SF/PDF/PPF/ISF computations is to use the
    `scipy.stats.ksone` distribution.

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
    For historical reasons this function is exposed in `scipy.special`, but the
    recommended way to achieve the most accurate CDF/SF/PDF/PPF/ISF computations is to
    use the `scipy.stats.ksone` distribution.

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
