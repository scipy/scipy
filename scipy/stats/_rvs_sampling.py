# -*- coding: utf-8 -*-
import numpy as np
from scipy._lib._util import check_random_state
from scipy.interpolate import CubicHermiteSpline


def rvs_ratio_uniforms(pdf, umax, vmin, vmax, size=1, c=0, random_state=None):
    """
    Generate random samples from a probability density function using the
    ratio-of-uniforms method.

    Parameters
    ----------
    pdf : callable
        A function with signature `pdf(x)` that is proportional to the
        probability density function of the distribution.
    umax : float
        The upper bound of the bounding rectangle in the u-direction.
    vmin : float
        The lower bound of the bounding rectangle in the v-direction.
    vmax : float
        The upper bound of the bounding rectangle in the v-direction.
    size : int or tuple of ints, optional
        Defining number of random variates (default is 1).
    c : float, optional.
        Shift parameter of ratio-of-uniforms method, see Notes. Default is 0.
    random_state : {None, int, `numpy.random.Generator`,
                    `numpy.random.RandomState`}, optional

        If `seed` is None (or `np.random`), the `numpy.random.RandomState`
        singleton is used.
        If `seed` is an int, a new ``RandomState`` instance is used,
        seeded with `seed`.
        If `seed` is already a ``Generator`` or ``RandomState`` instance then
        that instance is used.

    Returns
    -------
    rvs : ndarray
        The random variates distributed according to the probability
        distribution defined by the pdf.

    Notes
    -----
    Given a univariate probability density function `pdf` and a constant `c`,
    define the set ``A = {(u, v) : 0 < u <= sqrt(pdf(v/u + c))}``.
    If `(U, V)` is a random vector uniformly distributed over `A`,
    then `V/U + c` follows a distribution according to `pdf`.

    The above result (see [1]_, [2]_) can be used to sample random variables
    using only the pdf, i.e. no inversion of the cdf is required. Typical
    choices of `c` are zero or the mode of `pdf`. The set `A` is a subset of
    the rectangle ``R = [0, umax] x [vmin, vmax]`` where

    - ``umax = sup sqrt(pdf(x))``
    - ``vmin = inf (x - c) sqrt(pdf(x))``
    - ``vmax = sup (x - c) sqrt(pdf(x))``

    In particular, these values are finite if `pdf` is bounded and
    ``x**2 * pdf(x)`` is bounded (i.e. subquadratic tails).
    One can generate `(U, V)` uniformly on `R` and return
    `V/U + c` if `(U, V)` are also in `A` which can be directly
    verified.

    The algorithm is not changed if one replaces `pdf` by k * `pdf` for any
    constant k > 0. Thus, it is often convenient to work with a function
    that is proportional to the probability density function by dropping
    unneccessary normalization factors.

    Intuitively, the method works well if `A` fills up most of the
    enclosing rectangle such that the probability is high that `(U, V)`
    lies in `A` whenever it lies in `R` as the number of required
    iterations becomes too large otherwise. To be more precise, note that
    the expected number of iterations to draw `(U, V)` uniformly
    distributed on `R` such that `(U, V)` is also in `A` is given by
    the ratio ``area(R) / area(A) = 2 * umax * (vmax - vmin) / area(pdf)``,
    where `area(pdf)` is the integral of `pdf` (which is equal to one if the
    probability density function is used but can take on other values if a
    function proportional to the density is used). The equality holds since
    the area of `A` is equal to 0.5 * area(pdf) (Theorem 7.1 in [1]_).
    If the sampling fails to generate a single random variate after 50000
    iterations (i.e. not a single draw is in `A`), an exception is raised.

    If the bounding rectangle is not correctly specified (i.e. if it does not
    contain `A`), the algorithm samples from a distribution different from
    the one given by `pdf`. It is therefore recommended to perform a
    test such as `~scipy.stats.kstest` as a check.

    References
    ----------
    .. [1] L. Devroye, "Non-Uniform Random Variate Generation",
       Springer-Verlag, 1986.

    .. [2] W. Hoermann and J. Leydold, "Generating generalized inverse Gaussian
       random variates", Statistics and Computing, 24(4), p. 547--557, 2014.

    .. [3] A.J. Kinderman and J.F. Monahan, "Computer Generation of Random
       Variables Using the Ratio of Uniform Deviates",
       ACM Transactions on Mathematical Software, 3(3), p. 257--260, 1977.

    Examples
    --------
    >>> from scipy import stats
    >>> rng = np.random.default_rng()

    Simulate normally distributed random variables. It is easy to compute the
    bounding rectangle explicitly in that case. For simplicity, we drop the
    normalization factor of the density.

    >>> f = lambda x: np.exp(-x**2 / 2)
    >>> v_bound = np.sqrt(f(np.sqrt(2))) * np.sqrt(2)
    >>> umax, vmin, vmax = np.sqrt(f(0)), -v_bound, v_bound
    >>> rvs = stats.rvs_ratio_uniforms(f, umax, vmin, vmax, size=2500,
    ...                                random_state=rng)

    The K-S test confirms that the random variates are indeed normally
    distributed (normality is not rejected at 5% significance level):

    >>> stats.kstest(rvs, 'norm')[1]
    0.250634764150542

    The exponential distribution provides another example where the bounding
    rectangle can be determined explicitly.

    >>> rvs = stats.rvs_ratio_uniforms(lambda x: np.exp(-x), umax=1,
    ...                                vmin=0, vmax=2*np.exp(-1), size=1000,
    ...                                random_state=rng)
    >>> stats.kstest(rvs, 'expon')[1]
    0.21121052054580314

    """
    if vmin >= vmax:
        raise ValueError("vmin must be smaller than vmax.")

    if umax <= 0:
        raise ValueError("umax must be positive.")

    size1d = tuple(np.atleast_1d(size))
    N = np.prod(size1d)  # number of rvs needed, reshape upon return

    # start sampling using ratio of uniforms method
    rng = check_random_state(random_state)
    x = np.zeros(N)
    simulated, i = 0, 1

    # loop until N rvs have been generated: expected runtime is finite.
    # to avoid infinite loop, raise exception if not a single rv has been
    # generated after 50000 tries. even if the expected numer of iterations
    # is 1000, the probability of this event is (1-1/1000)**50000
    # which is of order 10e-22
    while simulated < N:
        k = N - simulated
        # simulate uniform rvs on [0, umax] and [vmin, vmax]
        u1 = umax * rng.uniform(size=k)
        v1 = rng.uniform(vmin, vmax, size=k)
        # apply rejection method
        rvs = v1 / u1 + c
        accept = (u1**2 <= pdf(rvs))
        num_accept = np.sum(accept)
        if num_accept > 0:
            x[simulated:(simulated + num_accept)] = rvs[accept]
            simulated += num_accept

        if (simulated == 0) and (i*N >= 50000):
            msg = ("Not a single random variate could be generated in {} "
                   "attempts. The ratio of uniforms method does not appear "
                   "to work for the provided parameters. Please check the "
                   "pdf and the bounds.".format(i*N))
            raise RuntimeError(msg)
        i += 1

    return np.reshape(x, size1d)


class FastNumericalInverse():
    r"""
    Object representing a Fast Numerical Inverse of a probability distribution

    The initializer of `FastNumericalInverse` accepts `dist`, a frozen instance
    of `scipy.stats.rv_continuous`, and provides an object with methods
    that approximate `dist.ppf` and `dist.rvs`. For some distributions,
    these methods may be much faster than those of `dist` itself.

    Parameters
    ----------
    dist : scipy.stats.rv_frozen
        Frozen distribution for which a fast numerical inverse is desired
    tol : float, optional
        u-error tolerance. The default is 1e-12.
    max_intervals : int, optional
        Maximum number of intervals in the cubic Hermite Spline used to
        approximate the percent point function. The default is 100000.

    Attributes
    -------
    intervals : int
        The number of intervals of the interpolant
    midpoint_error : float
        The maximum u-error at an interpolant interval midpoint

    Notes
    -----
    `FastNumericalInverse` approximates the inverse of a continuous
    statistical distribution's CDF with a cubic Hermite spline.

    As described in [1]_, it begins by evaluating the distribution's PDF and
    CDF at a mesh of quantiles ``x`` within the distribution's support.
    It uses the results to fit a cubic Hermite spline ``H`` such that
    ``H(p) == x``, where ``p`` is the array of percentiles corresponding
    with the quantiles ``x``. Therefore, the spline approximates the inverse
    of the distribution's CDF to machine precision at the percentiles ``p``,
    but typically, the spline will not be as accurate at the midpoints between
    the percentile points::

        p_mid = (p[:-1] + p[1:])/2

    so the mesh of quantiles is refined as needed to reduce the maximum
    "u-error"::

        u_error = np.max(np.abs(dist.cdf(H(p_mid)) - p_mid))

    below the specified tolerance `tol`. Refinement stops when the required
    tolerance is achieved or when the number of mesh intervals after the next
    refinement could exceed the maximum allowed number `max_intervals`.

    The approximation will not be accurate for the extreme tails of the
    distribution (e.g. ``H(1e-17)``, ``H(1)``).

    References
    ----------
    .. [1] Hörmann, Wolfgang, and Josef Leydold. "Continuous random variate
           generation by fast numerical inversion." ACM Transactions on
           Modeling and Computer Simulation (TOMACS) 13.4 (2003): 347-362.

    Examples
    --------
    For some distributions, ``dist.ppf`` are ``dist.rvs`` are quite slow.

    >>> import numpy as np
    >>> from scipy.stats import genexpon
    >>> from timeit import timeit
    >>> time_once = lambda f: f"{timeit(f, number=1)*1000:.6} ms"
    >>> dist = genexpon(9, 16, 3)  # freeze the distribution
    >>> p = np.linspace(0.01, 0.99, 99)  # percentiles from 1% to 99%
    >>> time_once(lambda: dist.ppf(p))
    '154.565 ms'  # may vary

    >>> time_once(lambda: dist.rvs(size=100))
    '148.979 ms'  # may vary

    The `FastNumericalInverse` has a method that approximates ``dist.ppf``.

    >>> from scipy.stats import FastNumericalInverse
    >>> fni = FastNumericalInverse(dist)
    >>> np.allclose(fni.ppf(p), dist.ppf(p))
    True

    In some cases, it is faster to both generate the fast numerical inverse
    and use it than to call ``dist.ppf``.

    >>> def time_me():
    ...     fni = FastNumericalInverse(dist)
    ...     fni.ppf(p)
    >>> time_once(time_me)
    '11.9222 ms'  # may vary

    After generating the fast numerical inverse, subsequent calls to its
    methods are much faster.
    >>> time_once(lambda: fni.ppf(p))
    '0.0819 ms'  # may vary

    The fast numerical inverse can also be used to generate random variates
    using inverse transform sampling.

    >>> time_once(lambda: fni.rvs(size=100))
    '0.0911 ms'  # may vary

    Depending on the implementation of the distribution's random sampling
    method, the random variates generated may be nearly identical, given
    the same random state.

    >>> rvs1 = dist.rvs(size=100, random_state=0)
    >>> rvs2 = fni.rvs(size=100, random_state=0)
    >>> np.allclose(rvs1, rvs2)
    True

    """
    def __init__(self, dist, tol=1e-12, max_intervals=100000):
        H, eu, intervals = _fast_numerical_inverse(dist, tol, max_intervals)
        self.H = H
        self.midpoint_error = eu
        self.intervals = intervals

    def ppf(self, q):
        r'''
        Approximate percent point function (inverse of `cdf`) at q of the
        given RV.

        Parameters
        ----------
        q : array_like
            lower tail probability

        Returns
        -------
        x : array_like
            quantile corresponding to the lower tail probability q.

        '''
        q = np.asarray(q)  # no harm; self.H always returns an array
        result = np.zeros_like(q, dtype=np.float64)
        i = (q >= 0) & (q <= 1)
        result[i] = self.H(q[i])
        result[~i] = np.nan
        return result

    def rvs(self, size=1, random_state=None):
        """
        Random variates of the given RV.

        Parameters
        ----------
        size : int or tuple of ints, optional
            Defining number of random variates (Default is 1).
        random_state : (see below), optional
            This parameter defines the object to use for drawing random
            variates.
            If `random_state` is `None` the `~np.random.RandomState` singleton
            is used.
            If `random_state` is an int, a new ``RandomState`` instance is
            used, seeded with `random_state`.
            If `random_state` is already a ``RandomState`` or ``Generator``
            instance, then that object is used.
            Default is None.

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of given `size`.
        """
        rng = check_random_state(random_state)
        return self.ppf(rng.random(size=size))


def _fni_input_validation(dist, tol, max_intervals):
    if int(max_intervals) != max_intervals or max_intervals <= 1:
        raise ValueError("`max_intervals' must be an integer greater than 1.")

    tol = float(tol)  # if there's an exception, raise it now

    return dist, tol, max_intervals


def _fast_numerical_inverse(dist, tol=1e-12, max_intervals=100000):
    """
    Generate a fast approximate PPF (inverse CDF) of a probability distribution

    `_fast_numerical_inverse` accepts `dist`, a frozen instance of
    `scipy.stats.rv_continuous`, and returns a object `fni` with methods
    that approximate `dist.ppf` and `dist.rvs`. For some distributions,
    these methods may be faster than those of `dist` itself.

    Parameters
    ----------
    dist : scipy.stats.rv_frozen
        Frozen distribution for which fast approximate PPF is desired
    tol : float, optional
        u-error tolerance. The default is 1e-12.
    max_intervals : int, optional
        Maximum number of intervals in the cubic Hermite Spline used to
        approximate the percent point function. The default is 100000.

    Returns
    -------
    H: scipy.interpolate.CubicHermiteSpline
        Interpolant of the distributions's PPF
    intervals : int
        The number of intervals of the interpolant
    midpoint_error : float
        The maximum u-error at an interpolant interval midpoint

    """
    dist, tol, max_intervals = _fni_input_validation(dist, tol, max_intervals)

    # [1] Section 2.1: "For distributions with unbounded domain, we have to
    # chop off its tails at [a] and [b] such that F(a) and 1-F(b) are small
    # compared to the maximal tolerated approximation error."
    p = np.array([dist.ppf(tol/10), dist.isf(tol/10)])  # initial interval

    # [1] Section 2.3: "We then halve this interval recursively until
    # |u[i+1]-u[i]| is smaller than some threshold value, for example, 0.05."
    while p.size-1 <= np.ceil(max_intervals/2):
        u = dist.cdf(p)
        i = np.nonzero(np.diff(u) > 0.05)[0]
        if not i.size:
            break

        p_mid = (p[i] + p[i+1])/2
        p = np.sort(np.concatenate((p, p_mid)))

    # [1] Section 2.3: "Now we continue with checking the error estimate in
    # each of the intervals and continue with splitting them until [it] is
    # smaller than a given error bound."
    while p.size-1 <= max_intervals:
        # [1] Equation 4-8
        u = dist.cdf(p)
        f = dist.pdf(p)
        H = CubicHermiteSpline(u, p, 1/f)

        # [1] Equation 12
        u_mid = (u[:-1] + u[1:])/2
        eu = np.abs(dist.cdf(H(u_mid)) - u_mid)

        i = np.nonzero(eu > tol)[0]
        if not i.size:
            break

        p_mid = (p[i] + p[i+1])/2
        p = np.sort(np.concatenate((p, p_mid)))

    # todo: add test for monotonicity [1] Section 2.4
    # todo: deal with vanishing density [1] Section 2.5
    return H, eu, p.size-1
