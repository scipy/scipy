
from functools import lru_cache
from dataclasses import dataclass

import numpy as np

from scipy.optimize import brentq
from scipy._lib._bunch import _make_tuple_bunch
from .distributions import hypergeom


def _sample_odds_ratio(table):
    """
    Given a table [[a, b], [c, d]], compute a*d/(b*c).

    Return nan if the numerator and denominator are 0.
    Return inf if just the denominator is 0.
    """
    # table must be a 2x2 numpy array.
    if table[1, 0] > 0 and table[0, 1] > 0:
        oddsratio = table[0, 0] * table[1, 1] / (table[1, 0] * table[0, 1])
    elif table[0, 0] == 0 or table[1, 1] == 0:
        oddsratio = np.nan
    else:
        oddsratio = np.inf
    return oddsratio


def _hypergeom_support_bounds(total, ngood, nsample):
    """
    Returns the inclusive bounds of the support.
    """
    nbad = total - ngood
    return max(0, nsample - nbad), min(nsample, ngood)


@lru_cache()
def _hypergeom_support_and_logpmf(total, ngood, nsample):
    """
    The support and the log of the PMF of the hypergeometric distribution.
    """
    lo, hi = _hypergeom_support_bounds(total, ngood, nsample)
    support = np.arange(lo, hi + 1)
    logpmf = hypergeom._logpmf(support, total, ngood, nsample)
    return support, logpmf


def _nc_hypergeom_support_and_pmf(lognc, total, ngood, nsample):
    """
    The support and the PMF of the noncentral hypergeometric distribution.

    The first argument is the log of the noncentrality parameter.
    """
    # Get the support and log of the PMF for the corresponding
    # hypergeometric distribution.
    support, hg_logpmf = _hypergeom_support_and_logpmf(total, ngood, nsample)

    # Calculate the PMF using the ideas from
    #   Bruce Barret (2017), A note on exact calculation of the non central
    #   hypergeometric distribution, Communications in Statistics - Theory
    #   and Methods, 46:13, 6737-6741, DOI: 10.1080/03610926.2015.1134573
    #
    # Create the log of the PMF from the log of the PMF of the hypergeometric
    # distribution.
    nchg_logpmf = hg_logpmf + lognc * support

    # Exponentiation after subtracting the max to get the "quasi" PMF (not
    # correctly normalized yet).
    nchg_qpmf = np.exp(nchg_logpmf - nchg_logpmf.max())

    # Normalize to the actual PMF of the noncentral hypergeometric dist.
    nchg_pmf = nchg_qpmf / nchg_qpmf.sum()

    return support, nchg_pmf


def _nc_hypergeom_mean(lognc, total, ngood, nsample):
    """
    Expected value of Fisher's noncentral hypergeometric distribution.

    lognc is the logarithm of the noncentrality parameter.  The remaining
    parameters are the same as those of the hypergeometric distribution.
    """
    # Explicit check for the edge cases.  (XXX These checks might not be
    # necessary, but double-check that we never call this function with +/- inf
    # before removing.  Maybe we'll get these values here inadvertently because
    # of overflow or underflow.)
    if lognc == -np.inf:
        return _hypergeom_support_bounds(total, ngood, nsample)[0]
    if lognc == np.inf:
        return _hypergeom_support_bounds(total, ngood, nsample)[1]

    support, nchg_pmf = _nc_hypergeom_support_and_pmf(lognc, total,
                                                      ngood, nsample)

    # Compute the expected value of the distribution.
    mean = support.dot(nchg_pmf)
    return mean


def _solve(func):
    """
    Solve func(lognc) = 0.  func must be an increasing function.
    """
    # We could just as well call the variable `x` instead of `lognc`, but we
    # always call this function with functions for which lognc (the log of
    # the noncentrality parameter) is the variable for which we are solving.
    # It also helps to remind that the linear search for a bracketing interval
    # is actually working with logarithms of the variable of interest, so
    # changing lognc by 5 corresponds to changing nc by a multiplicative
    # factor of exp(5) ~= 148.  While the strategy used here would be terrible
    # as a general purpose solver, here it works well.
    lognc = 0
    value = func(lognc)
    if value == 0:
        return lognc

    # Find a bracketing interval.
    if value > 0:
        lognc -= 5
        while func(lognc) > 0:
            lognc -= 5
        lo = lognc
        hi = lognc + 5
    else:
        lognc += 5
        while func(lognc) < 0:
            lognc += 5
        lo = lognc - 5
        hi = lognc

    # lo and hi bracket the solution for lognc.
    lognc = brentq(func, lo, hi, xtol=1e-13)
    return lognc


def _nc_hypergeom_mean_inverse(x, total, ngood, nsample):
    """
    For the given total, ngood, and nsample, find the noncentrality
    parameter of Fisher's noncentral hypergeometric distribution whose
    mean is x.
    """
    lognc = _solve(lambda lognc: _nc_hypergeom_mean(lognc, total, ngood,
                                                    nsample) - x)
    return np.exp(lognc)


def _nc_hypergeom_cdf(x, lognc, total, ngood, nsample):
    """
    CDF of the noncentral hypergeometric distribution.
    """
    lo, hi = _hypergeom_support_bounds(total, ngood, nsample)
    if lognc == -np.inf:
        return 1.0*(x >= lo)
    elif lognc == np.inf:
        return 1.0*(x >= hi)

    support, nchg_pmf = _nc_hypergeom_support_and_pmf(lognc, total,
                                                      ngood, nsample)
    return nchg_pmf[support <= x].sum()


def _nc_hypergeom_sf(x, lognc, total, ngood, nsample):
    """
    Survival function of the noncentral hypergeometric distribution.
    """
    lo, hi = _hypergeom_support_bounds(total, ngood, nsample)
    if lognc == -np.inf:
        return 1.0*(x < lo)
    elif lognc == np.inf:
        return 1.0*(x < hi)

    support, nchg_pmf = _nc_hypergeom_support_and_pmf(lognc, total,
                                                      ngood, nsample)
    return nchg_pmf[support > x].sum()


def _hypergeom_params_from_table(table):
    x = table[0, 0]
    total = table.sum()
    ngood = table[0].sum()
    nsample = table[:, 0].sum()
    return x, total, ngood, nsample


def _ci_upper(table, alpha):
    """
    Compute the upper end of the confidence interval.
    """
    if _sample_odds_ratio(table) == np.inf:
        return np.inf

    x, total, ngood, nsample = _hypergeom_params_from_table(table)

    # _nc_hypergeom_cdf is a decreasing function of lognc, so we negate
    # it in the lambda expression.
    lognc = _solve(lambda lognc: -_nc_hypergeom_cdf(x, lognc, total,
                                                    ngood, nsample) + alpha)
    return np.exp(lognc)


def _ci_lower(table, alpha):
    """
    Compute the lower end of the confidence interval.
    """
    if _sample_odds_ratio(table) == 0:
        return 0

    x, total, ngood, nsample = _hypergeom_params_from_table(table)

    lognc = _solve(lambda lognc: _nc_hypergeom_sf(x - 1, lognc, total,
                                                  ngood, nsample) - alpha)
    return np.exp(lognc)


def _conditional_oddsratio(table):
    """
    Conditional MLE of the odds ratio for the 2x2 contingency table.
    """
    x, total, ngood, nsample = _hypergeom_params_from_table(table)
    lo, hi = _hypergeom_support_bounds(total, ngood, nsample)

    # Check if x is at one of the extremes of the support.  If so, we know
    # the odds ratio is either 0 or inf.
    if x == lo:
        # x is at the low end of the support.
        return 0
    if x == hi:
        # x is at the high end of the support.
        return np.inf

    nc = _nc_hypergeom_mean_inverse(x, total, ngood, nsample)
    return nc


def _conditional_oddsratio_ci(table, confidence_level=0.95,
                              alternative='two-sided'):
    """
    Conditional exact confidence interval for the odds ratio.

    This function implements Cornfield's "exact confidence limits",
    as explained in section 2 of the paper:

        J. Cornfield (1956), A statistical problem arising from
        retrospective studies. In Neyman, J. (ed.), Proceedings of
        the Third Berkeley Symposium on Mathematical Statistics and
        Probability 4, pp. 135-148.

    """
    if alternative == 'two-sided':
        alpha = 0.5*(1 - confidence_level)
        lower = _ci_lower(table, alpha)
        upper = _ci_upper(table, alpha)
    elif alternative == 'less':
        lower = 0.0
        upper = _ci_upper(table, 1 - confidence_level)
    else:
        # alternative == 'greater'
        lower = _ci_lower(table, 1 - confidence_level)
        upper = np.inf

    return lower, upper


@dataclass
class ConfidenceInterval:
    low: float
    high: float


FisherExactBaseResult = _make_tuple_bunch('FisherExactBaseResult',
                                          ['sample_odds_ratio', 'pvalue'],
                                          ['table', 'alternative',
                                           'conditional_odds_ratio'])


class FisherExactResult(FisherExactBaseResult):

    def conditional_odds_ratio_ci(self, confidence_level=0.95):
        """
        Confidence interval for the conditional odds ratio.

        The limits of the confidence interval are the conditional "exact
        confidence limits" as defined in section 2 of Cornfield [2]_, and
        originally described by Fisher [1]_.  The conditional odds ratio
        and confidence interval are also discussed in Section 4.1.2 of the
        text by Sahai and Khurshid [3]_.

        Parameters
        ----------
        confidence_level: float
            Desired confidence level for the confidence interval.
            The value must be given as a fraction between 0 and 1.
            Default is 0.95 (meaning 95%).

        Returns
        -------
        ci : ``ConfidenceInterval`` instance
            The confidence interval, represented as an object with
            attributes ``low`` and ``high``.

        References
        ----------
        .. [1] R. A. Fisher (1935), The logic of inductive inference,
               Journal of the Royal Statistical Society, Vol. 98, No. 1,
               pp. 39-82.
        .. [2] J. Cornfield (1956), A statistical problem arising from
               retrospective studies. In Neyman, J. (ed.), Proceedings of
               the Third Berkeley Symposium on Mathematical Statistics
               and Probability 4, pp. 135-148.
        .. [3] H. Sahai and A. Khurshid (1996), Statistics in Epidemiology:
               Methods, Techniques, and Applications, CRC Press LLC, Boca
               Raton, Florida.
        """
        if confidence_level < 0 or confidence_level > 1:
            raise ValueError('confidence_level must be between 0 and 1')

        table = self.table
        if 0 in table.sum(axis=0) or 0 in table.sum(axis=1):
            # If both values in a row or column are zero, the p-value is 1,
            # the odds ratio is NaN and the confidence interval is (0, inf).
            ci = (0, np.inf)
        else:
            ci = _conditional_oddsratio_ci(table,
                                           confidence_level=confidence_level,
                                           alternative=self.alternative)
        return ConfidenceInterval(low=ci[0], high=ci[1])


def fisher_exact(table, alternative='two-sided'):
    r"""
    Perform a Fisher exact test on a 2x2 contingency table.

    Parameters
    ----------
    table : array_like of ints
        A 2x2 contingency table.  Elements must be non-negative integers.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available (default is 'two-sided'):

        * 'two-sided'
        * 'less': one-sided
        * 'greater': one-sided

    Returns
    -------
    result : object
        The returned object has three computed attributes:

        * ``sample_odds_ratio``, float, is
          ``table[0, 0]*table[1, 1]/(table[0, 1]*table[1, 0])``.
          This is the prior odds ratio and not a posterior estimate.
        * ``pvalue``, float, is the probability of obtaining a
          distribution at least as extreme as the one that was
          actually observed, assuming that the null hypothesis is
          true.
        * ``conditional_odds_ratio``, float, is the conditional
          maximum likelihood estimate for the odds ratio.   It is
          the noncentrality parameter of Fisher's noncentral
          hypergeometric distribution with the same hypergeometric
          parameters as ``table`` and whose mean is ``table[0, 0]``.

        The object also has these parameters that were passed in to
        the function:

        * ``table``
        * ``alternative``

        The object has the method ``conditional_odds_ratio_ci``
        to compute the confidence interval of the conditional
        odds ratio.

    See Also
    --------
    chi2_contingency : Chi-square test of independence of variables in a
        contingency table.

    Notes
    -----
    For tables with large numbers, the (inexact) chi-square test implemented
    in the function `chi2_contingency` can also be used.

    For backwards compatibility, the object returned by the function
    also acts like a sequence of length 2 that holds just the
    ``sample_odds_ratio`` and the ``pvalue`` fields.  The following
    code to unpack the two values,

        oddsratio, pvalue = fisher_exact(table)

    still works, but that loses the other parameters attached to the
    return value.  It is recommended to use the return value as an
    object with attributes, e.g.:

        result = fisher_exact(table)
        print(result.pvalue, result.conditional_odds_ratio)

    References
    ----------
    .. [1] R. A. Fisher (1935), The logic of inductive inference, Journal
           of the Royal Statistical Society, Vol. 98, No. 1, pp. 39-82.
    .. [2] J. Cornfield (1956), A statistical problem arising from
           retrospective studies. In Neyman, J. (ed.), Proceedings of
           the Third Berkeley Symposium on Mathematical Statistics and
           Probability 4, pp. 135-148.
    .. [3] H. Sahai and A. Khurshid (1996), Statistics in Epidemiology:
           Methods, Techniques, and Applications, CRC Press LLC, Boca
           Raton, Florida.

    Examples
    --------
    Say we spend a few days counting whales and sharks in the Atlantic and
    Indian oceans. In the Atlantic ocean we find 8 whales and 1 shark, in
    the Indian ocean 2 whales and 5 sharks. Then our contingency table is::

                Atlantic  Indian
        whales     8        2
        sharks     1        5

    We use this table to find the p-value:

    >>> from scipy.stats import fisher_exact
    >>> result = fisher_exact([[8, 2], [1, 5]])
    >>> result.pvalue
    0.03496503496503495

    Under the null hypothesis, the probability that we would observe this
    or an even more imbalanced ratio by chance is about 3.5%.  A commonly
    used significance level is 5%--if we adopt that, we can therefore
    conclude that our observed imbalance is statistically significant;
    whales prefer the Atlantic while sharks prefer the Indian ocean.

    In epidemiology, individuals are classified as "exposed" or
    "unexposed" to some factor or treatment. If the occurrence of some
    illness is under study, those who have the illness are often
    classifed as "cases", and those without it are "noncases".  The
    counts of the occurrences of these classes gives a contingency
    table::

                    exposed    unexposed
        cases          a           b
        noncases       c           d

    The sample odds ratio may be written ``(a/c) / (b/d)``.  ``a/c`` can
    be interpreted as the odds of a case occurring in the exposed group,
    and ``b/d`` as the odds of a case occurring in the unexposed group.
    The sample odds ratio is the ratio of these odds.  If the odds ratio
    is greater than 1, it suggests that there is a positive association
    between being exposed and being a case.

    Interchanging the rows or columns of the contingency table inverts
    the odds ratio, so it is import to understand the meaning of labels
    given to the rows and columns of the table when interpreting the
    odds ratio.  This is in contrast to the p-value, which is invariant
    under interchanging the rows or columns.


    Consider a hypothetical example where it is hypothesized that
    exposure to a certain chemical is assocated with increased occurrence
    of a certain disease.  Suppose we have the following table for a
    collection of 410 people::

                  exposed   unexposed
        cases         7         15
        noncases     58        472

    The question we ask is "Is exposure to the chemical associated with
    increased risk of the disease?"

    Compute the test, and first check the p-value.

    >>> test = fisher_exact([[7, 15], [58, 472]])
    >>> test.pvalue
    0.009208708293019454

    The p-value is less than 0.01, which suggests there is an
    association.

    >>> test.sample_odds_ratio
    3.7977011494252872
    >>> test.conditional_odds_ratio
    3.783668770554967

    The sample odds ratio is 3.80 and the conditional odds ratio
    is 3.78, again suggesting an association.  We can also check
    the 95% confidence interval of the conditional odds ratio.

    >>> test.conditional_odds_ratio_ci(confidence_level=0.95)
    ConfidenceInterval(low=1.2514829132265854, high=10.363493716701255)

    So the 95% confidence interval for the conditional odds ratio
    is approximately (1.25, 10.4), which is strong evidence that the
    odds ratio is greater than 1.
    """
    c = np.asarray(table)

    if c.shape != (2, 2):
        raise ValueError(f"Invalid shape {c.shape}. The input `table` must be "
                         "of shape (2, 2).")

    if not np.issubdtype(c.dtype, np.integer):
        raise ValueError("`table` must be an array of integers, but got "
                         f"type {c.dtype}")
    c = c.astype(np.int64)

    if np.any(c < 0):
        raise ValueError("All values in `table` must be nonnegative.")

    if 0 in c.sum(axis=0) or 0 in c.sum(axis=1):
        # If both values in a row or column are zero, the p-value is 1 and
        # the odds ratio is NaN.
        result = FisherExactResult(table=c, alternative=alternative,
                                   sample_odds_ratio=np.nan, pvalue=1.0,
                                   conditional_odds_ratio=np.nan)
        return result

    n1 = c[0, 0] + c[0, 1]
    n2 = c[1, 0] + c[1, 1]
    n = c[0, 0] + c[1, 0]

    def binary_search(n, n1, n2, side):
        """Binary search for where to begin halves in two-sided test."""
        if side == "upper":
            minval = mode
            maxval = n
        else:
            minval = 0
            maxval = mode
        guess = -1
        while maxval - minval > 1:
            if maxval == minval + 1 and guess == minval:
                guess = maxval
            else:
                guess = (maxval + minval) // 2
            pguess = hypergeom.pmf(guess, n1 + n2, n1, n)
            if side == "upper":
                ng = guess - 1
            else:
                ng = guess + 1
            if pguess <= pexact < hypergeom.pmf(ng, n1 + n2, n1, n):
                break
            elif pguess < pexact:
                maxval = guess
            else:
                minval = guess
        if guess == -1:
            guess = minval
        if side == "upper":
            while guess > 0 and (hypergeom.pmf(guess, n1 + n2, n1, n) <
                                 pexact * epsilon):
                guess -= 1
            while hypergeom.pmf(guess, n1 + n2, n1, n) > pexact / epsilon:
                guess += 1
        else:
            while hypergeom.pmf(guess, n1 + n2, n1, n) < pexact * epsilon:
                guess += 1
            while guess > 0 and (hypergeom.pmf(guess, n1 + n2, n1, n) >
                                 pexact / epsilon):
                guess -= 1
        return guess

    if alternative == 'less':
        pvalue = hypergeom.cdf(c[0, 0], n1 + n2, n1, n)
    elif alternative == 'greater':
        # Same formula as the 'less' case, but with the second column.
        pvalue = hypergeom.cdf(c[0, 1], n1 + n2, n1, c[0, 1] + c[1, 1])
    elif alternative == 'two-sided':
        mode = int((n + 1) * (n1 + 1) / (n1 + n2 + 2))
        pexact = hypergeom.pmf(c[0, 0], n1 + n2, n1, n)
        pmode = hypergeom.pmf(mode, n1 + n2, n1, n)

        epsilon = 1 - 1e-4
        if np.abs(pexact - pmode) / np.maximum(pexact, pmode) <= 1 - epsilon:
            pvalue = 1.0
        elif c[0, 0] < mode:
            plower = hypergeom.cdf(c[0, 0], n1 + n2, n1, n)
            if hypergeom.pmf(n, n1 + n2, n1, n) > pexact / epsilon:
                pvalue = plower
            else:
                guess = binary_search(n, n1, n2, "upper")
                pvalue = plower + hypergeom.sf(guess - 1, n1 + n2, n1, n)
        else:
            pupper = hypergeom.sf(c[0, 0] - 1, n1 + n2, n1, n)
            if hypergeom.pmf(0, n1 + n2, n1, n) > pexact / epsilon:
                pvalue = pupper
            else:
                guess = binary_search(n, n1, n2, "lower")
                pvalue = pupper + hypergeom.cdf(guess, n1 + n2, n1, n)
    else:
        msg = "`alternative` must be one of {'two-sided', 'less', 'greater'}"
        raise ValueError(msg)

    pvalue = min(pvalue, 1.0)

    oddsratio = _sample_odds_ratio(c)
    cond_odds = _conditional_oddsratio(c)
    result = FisherExactResult(table=c, alternative=alternative,
                               sample_odds_ratio=oddsratio, pvalue=pvalue,
                               conditional_odds_ratio=cond_odds)
    return result
