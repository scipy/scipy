from __future__ import division, print_function, absolute_import

import numpy as np
import warnings
from scipy._lib._util import check_random_state


def rvs_ratio_uniforms(pdf, umax, vmin, vmax, size=1, c=0, random_state=None):
    """
    Generate random samples from a probability density function using the
    ratio-of-uniforms method.

    Parameters
    ----------
    pdf : callable
        A function with signature `pdf(x)` that is the probability
        density function of the distribution.
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
    random_state : int or np.random.RandomState instance, optional
        If already a RandomState instance, use it.
        If seed is an int, return a new RandomState instance seeded with seed.
        If None, use np.random.RandomState. Default is None.

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

    Intuitively, the method works well if `A` fills up most of the
    enclosing rectangle such that the probability is high that `(U, V)`
    lies in `A` whenever it lies in `R` as the number of required
    iterations becomes too large otherwise. To be more precise, note that
    the expected number of iterations to draw `(U, V)` uniformly
    distributed on `R` such that `(U, V)` is also in `A` is given by
    the ratio ``area(R) / area(A) = 2 * umax * (vmax - vmin)``, using the fact
    that the area of `A` is equal to 1/2 (Theorem 7.1 in [1]_). A warning
    is displayed if this ratio is larger than 20. Moreover, if the sampling
    fails to generate a single random variate after 50000 iterations (i.e.
    not a single draw is in `A`), an exception is raised.

    If the bounding rectangle is not correctly specified (i.e. if it does not
    contain `A`), the algorithm samples from a distribution different from
    the one given by `pdf`. It is therefore recommended to perform a
    test such as `stats.kstest` as a check.

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

    Simulate normally distributed random variables. It is easy to compute the
    bounding rectangle explicitly in that case.

    >>> f = stats.norm.pdf
    >>> v_bound = np.sqrt(f(np.sqrt(2))) * np.sqrt(2)
    >>> umax, vmin, vmax = np.sqrt(f(0)), -v_bound, v_bound
    >>> np.random.seed(12345)
    >>> rvs = stats.rvs_ratio_uniforms(f, umax, vmin, vmax, size=2500)

    The K-S test confirms that the random variates are indeed normally
    distributed (normality is not rejected at 5% significance level):

    >>> stats.kstest(rvs, 'norm')[1]
    0.3420173467307603

    The exponential distribution provides another example where the bounding
    rectangle can be determined explicitly.

    >>> np.random.seed(12345)
    >>> rvs = stats.rvs_ratio_uniforms(lambda x: np.exp(-x), umax=1,
    ...                                vmin=0, vmax=2*np.exp(-1), size=1000)
    >>> stats.kstest(rvs, 'expon')[1]
    0.928454552559516

    Sometimes it can be helpful to use a non-zero shift parameter `c`, see e.g.
    [2]_ above in the case of the generalized inverse Gaussian distribution.

    """

    if vmin >= vmax:
        raise ValueError("vmin must be smaller than vmax.")

    if umax <= 0:
        raise ValueError("umax must be positive.")

    exp_iter = 2 * (vmax - vmin) * umax  # rejection constant (see [1])
    if exp_iter > 20:
        msg = ("The expected number of iterations to generate a single random "
               "number from the desired distribution is larger than {}, "
               "potentially causing bad performance.".format(int(exp_iter)))
        warnings.warn(msg, RuntimeWarning)

    size1d = tuple(np.atleast_1d(size))
    N = np.prod(size1d)  # number of rvs needed, reshape upon return

    # start sampling using ratio of uniforms method
    rng = check_random_state(random_state)
    x = np.zeros(N)
    simulated, i = 0, 1

    # loop until N rvs have been generated: expected runtime is finite
    # to avoid infinite loop, raise exception if not a single rv has been
    # generated after 50000 tries. even if exp_iter = 1000, probability of
    # this event is (1-1/1000)**50000 which is of order 10e-22
    while True:
        k = N - simulated
        # simulate uniform rvs on [0, umax] and [vmin, vmax]
        u1 = umax * rng.random_sample(size=k)
        v1 = vmin + (vmax - vmin) * rng.random_sample(size=k)
        # apply rejection method
        rvs = v1 / u1 + c
        accept = (u1**2 <= pdf(rvs))
        num_accept = np.sum(accept)
        if num_accept > 0:
            take = min(num_accept, N - simulated)
            x[simulated:(simulated + take)] = rvs[accept][0:take]
            simulated += take
        if simulated >= N:
            return np.reshape(x, size1d)
        if (simulated == 0) and (i*N >= 50000):
            msg = ("Not a single random variate could be generated in {} "
                   "attempts. The ratio of uniforms method does not appear "
                   "to work for the provided parameters. Please check the "
                   "pdf and the bounds.".format(i*N))
            raise RuntimeError(msg)
        i += 1
