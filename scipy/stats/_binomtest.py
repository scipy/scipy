from collections import namedtuple
import operator
from math import sqrt
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtri
from ._discrete_distns import binom


ConfidenceInterval = namedtuple('ConfidenceInterval', ['low', 'high'])


class BinomTestResult:

    def __init__(self, k, n, alternative, pvalue, proportion_estimate):
        self.k = k
        self.n = n
        self.alternative = alternative
        self.proportion_estimate = proportion_estimate
        self.pvalue = pvalue

    def __repr__(self):
        s = ("BinomTestResult("
             f"k={self.k}, "
             f"n={self.n}, "
             f"alternative={self.alternative!r}, "
             f"proportion_estimate={self.proportion_estimate}, "
             f"pvalue={self.pvalue})")
        return s

    def proportion_ci(self, confidence_level=0.95, method='exact'):
        """
        Compute the confidence interval for the estimated proportion.

        Parameters
        ----------
        confidence_level : float, optional
            Confidence level for the computed confidence interval
            of the estimated proportion. Default is 0.95.
        method : {'exact', 'wilson', 'wilsoncc'}, optional
            Selects the method used to compute the confidence interval
            for the estimate of the proportion:

                ``'exact'`` :
                    Use the Clopper-Pearson exact method
                ``'wilson'`` :
                    Wilson's method, without continuity correction
                ``'wilsoncc'`` :
                    Wilson's method, with continuity correction

            Default is ``'exact'``.

        Returns
        -------
        ci : A namedtuple containing the lower and upper limits
            of the confidence interval.
        """
        if method not in ('exact', 'wilson', 'wilsoncc'):
            raise ValueError("method must be one of 'exact', 'wilson' or "
                             "'wilsoncc'.")
        if not (0 <= confidence_level <= 1):
            raise ValueError('confidence_level must be in the interval [0, 1].')
        if method == 'exact':
            low, high = _binom_exact_conf_int(self.k, self.n,
                                              confidence_level,
                                              self.alternative)
        else:
            # method is 'wilson' or 'wilsoncc'
            low, high = _binom_wilson_conf_int(self.k, self.n,
                                               confidence_level,
                                               self.alternative,
                                               correction=method == 'wilsoncc')
        return ConfidenceInterval(low=low, high=high)


def _binom_exact_conf_int(k, n, confidence_level, alternative):
    """
    Compute the estimate and confidence interval for the binomial test.

    Returns proportion, prop_low, prop_high
    """
    if alternative == 'two-sided':
        alpha = (1 - confidence_level) / 2
        if k == 0:
            plow = 0.0
        else:
            plow = brentq(lambda p: binom.sf(k-1, n, p) - alpha, 0, 1)
        if k == n:
            phigh = 1.0
        else:
            phigh = brentq(lambda p: binom.cdf(k, n, p) - alpha, 0, 1)
    elif alternative == 'less':
        alpha = 1 - confidence_level
        plow = 0.0
        if k == n:
            phigh = 1.0
        else:
            phigh = brentq(lambda p: binom.cdf(k, n, p) - alpha, 0, 1)
    elif alternative == 'greater':
        alpha = 1 - confidence_level
        if k == 0:
            plow = 0.0
        else:
            plow = brentq(lambda p: binom.sf(k-1, n, p) - alpha, 0, 1)
        phigh = 1.0
    return plow, phigh


def _binom_wilson_conf_int(k, n, confidence_level, alternative, correction):
    # This function assumes that the arguments have already been validated.
    # In particular, `alternative` must be one of 'two-sided', 'less' or
    # 'greater'.
    p = k / n
    if alternative == 'two-sided':
        z = ndtri(0.5 + 0.5*confidence_level)
    else:
        z = ndtri(confidence_level)

    t = 1 + z**2/n
    r = (p + z**2/(2*n)) / t

    if correction:
        if alternative == 'less' or k == 0:
            lo = 0.0
        else:
            dlo = ((z * sqrt(z**2 - 1/n + 4*n*p*(1 - p) + (4*p - 2)) + 1) /
                   (2*n*t))
            lo = r - dlo
        if alternative == 'greater' or k == n:
            hi = 1.0
        else:
            dhi = ((z * sqrt(z**2 - 1/n + 4*n*p*(1 - p) - (4*p - 2)) + 1) /
                   (2*n*t))
            hi = r + dhi
    else:
        d = z/t * sqrt(p*(1-p)/n + (z/(2*n))**2)
        if alternative == 'less' or k == 0:
            lo = 0.0
        else:
            lo = r - d
        if alternative == 'greater' or k == n:
            hi = 1.0
        else:
            hi = r + d

    return lo, hi


def _validate_nonneg_int(k, name):
    try:
        k = operator.index(k)
    except TypeError:
        raise TypeError(f'{name} must be an integer.')
    if k < 0:
        raise ValueError(f'{name} must be nonnegative')
    return k


def binomtest(k, *, n=None, p=0.5, alternative='two-sided'):
    """
    Perform a test that the probability of success is p.

    This is a test of the null hypothesis that the probability of success
    in a Bernoulli experiment is `p`.

    Parameters
    ----------
    k : int
        The number of successes.
    n : int
        The number of trials.
    p : float, optional
        The hypothesized probability of success.  ``0 <= p <= 1``. The
        default value is ``p = 0.5``.
    alternative : {'two-sided', 'greater', 'less'}, optional
        Indicates the alternative hypothesis. The default value is
        'two-sided'.

    Returns
    -------
    result : BinomTestResult instance
        The return value is an object with the following attributes:

        * `k` : int
        * `n` : int
        * `alternative` : str
        * `pvalue` : float
        * `proportion_estimate` : float

        The object has the following methods:

        * `proportion_ci(confidence_level=0.95, method='exact')`
            Compute the confidence interval for ``proportion_estimate``.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Binomial_test

    Examples
    --------
    >>> from scipy.stats import binomtest

    A car manufacturer claims that no more than 10% of their cars are unsafe.
    15 cars are inspected for safety, 3 were found to be unsafe. Test the
    manufacturer's claim:

    >>> result = binomtest(3, n=15, p=0.1, alternative='greater')
    >>> result.pvalue
    0.18406106910639114

    The null hypothesis cannot be rejected at the 5% level of significance
    because the returned p-value is greater than the critical value of 5%.

    The estimated proportion is simply ``3/15``:

    >>> result.proportion_estimate
    0.2

    We can use the `proportion_ci()` method of the result to compute the
    confidence interval of the estimate:

    >>> result.proportion_ci(confidence_level=0.95)
    ConfidenceInterval(low=0.056846867590246826, high=1.0)

    """
    k = _validate_nonneg_int(k, 'k')
    n = _validate_nonneg_int(n, 'n')
    if k > n:
        raise ValueError('k must not be greater than n.')

    if (p > 1.0) or (p < 0.0):
        raise ValueError("p must be in range [0,1]")

    if alternative not in ('two-sided', 'less', 'greater'):
        raise ValueError("alternative not recognized; \n"
                         "must be 'two-sided', 'less' or 'greater'")

    if alternative == 'less':
        pval = binom.cdf(k, n, p)
    elif alternative == 'greater':
        pval = binom.sf(k-1, n, p)
    else:
        # alternative is 'two-sided'
        d = binom.pmf(k, n, p)
        rerr = 1 + 1e-7
        if k == p * n:
            # special case as shortcut, would also be handled by `else` below
            pval = 1.
        elif k < p * n:
            i = np.arange(np.ceil(p * n), n+1)
            y = np.sum(binom.pmf(i, n, p) <= d*rerr, axis=0)
            pval = binom.cdf(k, n, p) + binom.sf(n - y, n, p)
        else:
            i = np.arange(np.floor(p*n) + 1)
            y = np.sum(binom.pmf(i, n, p) <= d*rerr, axis=0)
            pval = binom.cdf(y-1, n, p) + binom.sf(k-1, n, p)

        pval = min(1.0, pval)

    result = BinomTestResult(k=k, n=n, alternative=alternative,
                             proportion_estimate=k/n, pvalue=pval)
    return result
