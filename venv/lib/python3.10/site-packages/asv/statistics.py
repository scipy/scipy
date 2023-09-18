# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Author: Pauli Virtanen, 2016

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import math
import operator

from .util import inf, nan, is_na


def compute_stats(samples, number):
    """
    Statistical analysis of measured samples.

    Parameters
    ----------
    samples : list of float
        List of total times (y) of benchmarks.
    number : int
        Repeat number for each sample.

    Returns
    -------
    beta_hat : float
        Estimated time per iteration
    stats : dict
        Information on statistics of the estimator.

    """

    if len(samples) < 1:
        return None, None

    Y = list(samples)

    # Median and quantiles
    y_50, ci_50 = quantile_ci(Y, 0.5, alpha_min=0.99)
    y_25 = quantile(Y, 0.25)
    y_75 = quantile(Y, 0.75)

    # If nonparametric CI estimation didn't give an estimate,
    # use the credible interval of a bayesian posterior distribution.
    a, b = ci_50
    if (math.isinf(a) or math.isinf(b)) and len(Y) > 1:
        # Compute posterior distribution for location, assuming
        # exponential noise. The MLE is equal to the median.
        c = LaplacePosterior(Y)

        # Use the CI from that distribution to extend beyond sample
        # bounds
        if math.isinf(a):
            a = min(c.ppf(0.01/2), min(Y))
        if math.isinf(b):
            b = max(c.ppf(1 - 0.01/2), max(Y))

        ci_50 = (a, b)

    # Produce results
    mean = sum(Y) / len(Y)
    var = sum((yp - mean)**2 for yp in Y) / len(Y)   # mle
    std = math.sqrt(var)

    result = y_50

    stats = {'ci_99': list(ci_50),
             'q_25': y_25,
             'q_75': y_75,
             'min': min(Y),
             'max': max(Y),
             'mean': mean,
             'std': std,
             'repeat': len(Y),
             'number': number}

    return result, stats


def get_err(result, stats):
    """
    Return an 'error measure' suitable for informing the user
    about the spread of the measurement results.
    """
    a, b = stats['q_25'], stats['q_75']
    return (b - a)/2


def get_weight(stats):
    """
    Return a data point weight for the result.
    """
    if stats is None or 'ci_99' not in stats:
        return None

    try:
        return 2 / abs(stats['ci_99'][1] - stats['ci_99'][0])
    except ZeroDivisionError:
        return None


def is_different(samples_a, samples_b, stats_a, stats_b, p_threshold=0.002):
    """Check whether the samples are statistically different.

    If sample data is not provided, or the sample is too small, falls
    back to a pessimistic CI-based check. If it returns True, then the
    difference is statistically significant. If it returns False, it
    might or might not be statistically significant.

    Parameters
    ----------
    samples_a, samples_b
        Input samples
    stats_a, stats_b
        Input stats data

    """

    if samples_a is not None and samples_b is not None:
        # Raw data present: Mann-Whitney U test, but only if there's
        # enough data so that the test can return True
        a = [x for x in samples_a if not is_na(x)]
        b = [x for x in samples_b if not is_na(x)]

        p_min = 1 / binom(len(a) + len(b), min(len(a), len(b)))
        if p_min < p_threshold:
            _, p = mann_whitney_u(a, b)
            return p < p_threshold

    # If confidence intervals overlap, reject.
    # Corresponds to a test with ill-specified threshold p-value,
    # which generally can be significantly smaller than p <= 0.01
    # depending on the actual data. For normal test (known variance),
    # 0.00027 <= p <= 0.01.
    ci_a = stats_a['ci_99']
    ci_b = stats_b['ci_99']

    if ci_a[1] >= ci_b[0] and ci_a[0] <= ci_b[1]:
        return False

    return True


def quantile_ci(x, q, alpha_min=0.01):
    """
    Compute a quantile and a confidence interval.

    Assumes independence, but otherwise nonparametric.

    Parameters
    ----------
    x : list of float
        Samples
    q : float
        Quantile to compute, in [0,1].
    alpha_min : float, optional
        Limit for coverage.
        The result has coverage >= 1 - alpha_min.

    Returns
    -------
    m : float
        Quantile of x
    ci : tuple of floats
        Confidence interval (a, b), of coverage >= alpha_min.

    """

    y = sorted(x)
    n = len(y)

    alpha_min = min(alpha_min, 1 - alpha_min)
    pa = alpha_min / 2
    pb = 1 - pa

    a = -inf
    b = inf

    # It's known that
    #
    # Pr[X_{(r)} < m < X_{(s)}] = Pr[r <= K <= s-1], K ~ Bin(n,p)
    #
    # where cdf(m) = p defines the quantile.
    #
    # Simplest median CI follows by picking r,s such that
    #
    #   F(r;n,q) <= alpha/2
    #   F(s;n,q) >= 1 - alpha/2
    #
    #   F(k;n,q) = sum(binom_pmf(n, j, q) for j in range(k))
    #
    # Then (y[r-1], y[s-1]) is a CI.
    # If no such r or s exists, replace by +-inf.

    F = 0
    for k, yp in enumerate(y):
        F += binom_pmf(n, k, q)
        # F = F(k+1;n,q)

        if F <= pa:
            a = yp

        if F >= pb:
            b = yp
            break

    m = quantile(y, q)
    return m, (a, b)


def quantile(x, q):
    """
    Compute quantile/percentile of the data

    Parameters
    ----------
    x : list of float
        Data set
    q : float
        Quantile to compute, 0 <= q <= 1

    """
    if not 0 <= q <= 1:
        raise ValueError("Invalid quantile")

    y = sorted(x)
    n = len(y)

    z = (n - 1) * q
    j = int(math.floor(z))
    z -= j

    if j == n - 1:
        m = y[-1]
    else:
        m = (1 - z)*y[j] + z*y[j+1]

    return m


_mann_whitney_u_memo = {}

def mann_whitney_u(x, y, method='auto'):
    """
    Mann-Whitney U test

    Ties are handled conservatively, returning the least significant
    tie breaking.

    Parameters
    ----------
    x, y : list of float
        Samples to test
    method : {'auto', 'exact', 'normal'}
        Whether to compute p-value exactly of via normal approximation.
        The option 'auto' switches to approximation for sample size > 20.

    Returns
    -------
    u : int
        U-statistic
    p : float
        p-value for two-sided alternative

    References
    ----------
    .. [1] Mann & Whitney, Ann. Math. Statist. 18, 50 (1947).
    .. [2] Gibbons & Chakraborti, "Nonparametric statistical inference". (2003)

    """
    memo = _mann_whitney_u_memo
    if len(memo) > 100000:
        memo.clear()

    m = len(x)
    n = len(y)

    if method == 'auto':
        if max(m, n) > 20:
            method = 'normal'
        else:
            method = 'exact'

    u, ties = mann_whitney_u_u(x, y)

    # Conservative tie breaking
    if u <= m*n//2 and u + ties >= m*n//2:
        ties = m*n//2 - u

    ux1 = min(u, m*n - u)
    ux2 = min(u + ties, m*n - (u + ties))

    if ux1 >= ux2:
        ux = ux1
    else:
        u = u + ties
        ux = ux2

    # Get p-value
    if method == 'exact':
        p1 = mann_whitney_u_cdf(m, n, ux, memo)
        p2 = 1.0 - mann_whitney_u_cdf(m, n, max(m*n//2, m*n - ux - 1), memo)
        p = p1 + p2
    elif method == 'normal':
        N = m + n
        var = m*n*(N + 1) / 12
        z = (ux - m*n/2) / math.sqrt(var)
        cdf = 0.5 * math.erfc(-z / math.sqrt(2))
        p = 2 * cdf
    else:
        raise ValueError("Unknown method {!r}".format(method))

    return u, p


def mann_whitney_u_u(x, y):
    u = 0
    ties = 0
    for xx in x:
        for yy in y:
            if xx > yy:
                u += 1
            elif xx == yy:
                ties += 1
    return u, ties


def mann_whitney_u_cdf(m, n, u, memo=None):
    if memo is None:
        memo = {}
    cdf = 0
    for uu in range(u + 1):
         cdf += mann_whitney_u_pmf(m, n, uu, memo)
    return cdf


def mann_whitney_u_pmf(m, n, u, memo=None):
    if memo is None:
        memo = {}
    return mann_whitney_u_r(m, n, u, memo) / binom(m + n, m)


def mann_whitney_u_r(m, n, u, memo=None):
    """
    Number of orderings in Mann-Whitney U test.

    The PMF of U for samples of sizes (m, n) is given by
    p(u) = r(m, n, u) / binom(m + n, m).

    References
    ----------
    .. [1] Mann & Whitney, Ann. Math. Statist. 18, 50 (1947).
    """
    if u < 0:
        value = 0
    elif m == 0 or n == 0:
        value = 1 if u == 0 else 0
    else:
        # Don't bother figuring out table construction, memoization
        # sorts it out
        if memo is None:
            memo = {}
        key = (m, n, u)
        value = memo.get(key)
        if value is not None:
            return value

        value = (mann_whitney_u_r(m, n - 1, u, memo)
                 + mann_whitney_u_r(m - 1, n, u - n, memo))

        memo[key] = value
    return value


def binom_pmf(n, k, p):
    """Binomial pmf = (n choose k) p**k (1 - p)**(n - k)"""
    if not (0 <= k <= n):
        return 0
    if p == 0:
        return 1.0 * (k == 0)
    elif p == 1.0:
        return 1.0 * (k == n)

    logp = math.log(p)
    log1mp = math.log(1 - p)
    return math.exp(lgamma(1 + n) - lgamma(1 + n - k) - lgamma(1 + k)
                        + k*logp + (n - k)*log1mp)


_BERNOULLI = [1.0, -0.5, 0.166666666667, 0.0, -0.0333333333333, 0.0, 0.0238095238095]


def lgamma(x):
    """
    Log gamma function. Only implemented at integers.
    """

    if x <= 0:
        raise ValueError("Domain error")

    if x > 100:
        # DLMF 5.11.1
        r = 0.5 * math.log(2*math.pi) + (x - 0.5) * math.log(x) - x
        for k in range(1, len(_BERNOULLI)//2 + 1):
            r += _BERNOULLI[2*k] / (2*k*(2*k - 1) * x**(2*k - 1))
        return r

    # Fall back to math.factorial
    int_x = int(x)
    err_int = abs(x - int_x)

    if err_int < 1e-12 * abs(x):
        return math.log(math.factorial(int_x - 1))

    # Would need full implementation
    return nan


def binom(n, k):
    """
    Binomial coefficient (n over k)
    """
    n = operator.index(n)
    k = operator.index(k)
    if not 0 <= k <= n:
        return 0
    m = n + 1
    num = 1
    den = 1
    for j in range(1, min(k, n - k) + 1):
        num *= m - j
        den *= j
    return num // den


class LaplacePosterior(object):
    """
    Univariate distribution::

        p(beta|y) = N [sum(|y_j - beta|)]**(-nu-1)

    where N is the normalization factor.

    Parameters
    ----------
    y : list of float
        Samples
    nu : float, optional
        Degrees of freedom. Default: len(y)-1

    Notes
    -----
    This is the posterior distribution in the Bayesian model assuming
    Laplace distributed noise::

        p(y|beta,sigma) = N exp(- sum_j (1/sigma) |y_j - beta|)

        p(sigma) ~ 1/sigma

        nu = len(y) - 1

    The MLE for beta is median(y).

    Note that the same approach applied to a Gaussian model::

        p(y|beta,sigma) = N exp(- sum_j 1/(2 sigma^2) (y_j - beta)^2)

    results to::

        p(beta|y) = N T(t, m-1);   t = (beta - mean(y)) / (sstd(y) / sqrt(m))

    where ``T(t, nu)`` is the Student t-distribution pdf, which then gives
    the standard textbook formulas for the mean.

    """

    def __init__(self, y, nu=None):
        if len(y) == 0:
            raise ValueError("empty input")

        if nu is None:
            self.nu = len(y) - 1
        else:
            self.nu = nu

        # Sort input
        y = sorted(y)

        # Get location and scale so that data is centered at MLE, and
        # the unnormalized PDF at MLE has amplitude ~ 1/nu.
        #
        # Proper scaling of inputs is important to avoid overflows
        # when computing the unnormalized CDF integrals below.
        self.mle = quantile(y, 0.5)
        self._y_scale = sum(abs(yp - self.mle) for yp in y)
        self._y_scale *= self.nu**(1/(self.nu + 1))

        # Shift and scale
        if self._y_scale != 0:
            self.y = [(yp - self.mle)/self._y_scale for yp in y]
        else:
            self.y = [0 for yp in y]

        self._cdf_norm = None
        self._cdf_memo = {}

    def _cdf_unnorm(self, beta):
        """
        Unnormalized CDF of this distribution::

            cdf_unnorm(b) = int_{-oo}^{b} 1/(sum_j |y - b'|)**(m+1) db'

        """
        if beta != beta:
            return beta

        for k, y in enumerate(self.y):
            if y > beta:
                k0 = k
                break
        else:
            k0 = len(self.y)

        cdf = 0

        nu = self.nu

        # Save some work by memoizing intermediate results
        if k0 - 1 in self._cdf_memo:
            k_start = k0
            cdf = self._cdf_memo[k0 - 1]
        else:
            k_start = 0
            cdf = 0

        # Do the integral piecewise, resolving the absolute values
        for k in range(k_start, k0 + 1):
            c = 2*k - len(self.y)
            y = sum(self.y[k:]) - sum(self.y[:k])

            if k == 0:
                a = -inf
            else:
                a = self.y[k-1]

            if k == k0:
                b = beta
            else:
                b = self.y[k]

            if c == 0:
                term = (b - a) / y**(nu+1)
            else:
                term = 1/(nu*c) * ((a*c + y)**(-nu) - (b*c + y)**(-nu))

            cdf += max(0, term)  # avoid rounding error

            if k != k0:
                self._cdf_memo[k] = cdf

        if beta == inf:
            self._cdf_memo[len(self.y)] = cdf

        return cdf

    def _ppf_unnorm(self, cdfx):
        """
        Inverse function for _cdf_unnorm
        """
        # Find interval
        for k in range(len(self.y) + 1):
            if cdfx <= self._cdf_memo[k]:
                break

        # Invert on interval
        c = 2*k - len(self.y)
        y = sum(self.y[k:]) - sum(self.y[:k])
        nu = self.nu
        if k == 0:
            term = cdfx
        else:
            a = self.y[k-1]
            term = cdfx - self._cdf_memo[k-1]

        if k == 0:
            z = -nu*c*term
            if z > 0:
                beta = (z**(-1/nu) - y) / c
            else:
                beta = -inf
        elif c == 0:
            beta = a + term * y**(nu+1)
        else:
            z = (a*c + y)**(-nu) - nu*c*term
            if z > 0:
                beta = (z**(-1/nu) - y)/c
            else:
                beta = inf

        if k < len(self.y):
            beta = min(beta, self.y[k])

        return beta

    def pdf(self, beta):
        """
        Probability distribution function
        """
        return math.exp(self.logpdf(beta))

    def logpdf(self, beta):
        """
        Logarithm of probability distribution function
        """
        if self._y_scale == 0:
            return inf if beta == self.mle else -inf

        beta = (beta - self.mle) / self._y_scale

        if self._cdf_norm is None:
            self._cdf_norm = self._cdf_unnorm(inf)

        ws = sum(abs(yp - beta) for yp in self.y)
        m = self.nu
        return -(m+1)*math.log(ws) - math.log(self._cdf_norm) - math.log(self._y_scale)

    def cdf(self, beta):
        """
        Cumulative probability distribution function
        """
        if self._y_scale == 0:
            return 1.0*(beta > self.mle)

        beta = (beta - self.mle) / self._y_scale

        if self._cdf_norm is None:
            self._cdf_norm = self._cdf_unnorm(inf)
        return self._cdf_unnorm(beta) / self._cdf_norm

    def ppf(self, cdf):
        """
        Percent point function (inverse function for cdf)
        """
        if cdf < 0 or cdf > 1.0:
            return nan

        if self._y_scale == 0:
            return self.mle

        if self._cdf_norm is None:
            self._cdf_norm = self._cdf_unnorm(inf)

        cdfx = min(cdf * self._cdf_norm, self._cdf_norm)
        beta = self._ppf_unnorm(cdfx)
        return beta * self._y_scale + self.mle

