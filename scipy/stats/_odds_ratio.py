
from dataclasses import dataclass

import numpy as np

from scipy.special import ndtr, ndtri
from scipy.optimize import brentq
from .distributions import nchypergeom_fisher, hypergeom
from .stats import fisher_exact


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


def _solve(func):
    """
    Solve func(nc) = 0.  func must be an increasing function.
    """
    # We could just as well call the variable `x` instead of `nc`, but we
    # always call this function with functions for which nc (the noncentrality
    # parameter) is the variable for which we are solving.
    nc = 1.0
    value = func(nc)
    if value == 0:
        return nc

    # Multiplicative factor by which to increase or decrease nc when
    # searching for a bracketing interval.
    factor = 2.0
    # Find a bracketing interval.
    if value > 0:
        nc /= factor
        while func(nc) > 0:
            nc /= factor
        lo = nc
        hi = factor*nc
    else:
        nc *= factor
        while func(nc) < 0:
            nc *= factor
        lo = nc/factor
        hi = nc

    # lo and hi bracket the solution for nc.
    nc = brentq(func, lo, hi, xtol=1e-13)
    return nc


def _nc_hypergeom_mean_inverse(x, total, ngood, nsample):
    """
    For the given total, ngood, and nsample, find the noncentrality
    parameter of Fisher's noncentral hypergeometric distribution whose
    mean is x.
    """
    nc = _solve(lambda nc: nchypergeom_fisher.mean(total, ngood, nsample,
                                                   nc) - x)
    return nc


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

    # nchypergeom_fisher.cdf is a decreasing function of nc, so we negate
    # it in the lambda expression.
    nc = _solve(lambda nc: -nchypergeom_fisher.cdf(x, total, ngood, nsample,
                                                   nc) + alpha)
    return nc


def _ci_lower(table, alpha):
    """
    Compute the lower end of the confidence interval.
    """
    if _sample_odds_ratio(table) == 0:
        return 0

    x, total, ngood, nsample = _hypergeom_params_from_table(table)

    nc = _solve(lambda nc: nchypergeom_fisher.sf(x - 1, total, ngood, nsample,
                                                 nc) - alpha)
    return nc


def _conditional_oddsratio(table):
    """
    Conditional MLE of the odds ratio for the 2x2 contingency table.
    """
    x, total, ngood, nsample = _hypergeom_params_from_table(table)
    lo, hi = hypergeom.support(total, ngood, nsample)

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


def _sample_odds_ratio_ci(table, confidence_level=0.95,
                          alternative='two-sided'):
    oddsratio = _sample_odds_ratio(table)
    log_or = np.log(oddsratio)
    se = np.sqrt((1/table).sum())
    if alternative == 'less':
        z = ndtri(confidence_level)
        loglow = -np.inf
        loghigh = log_or + z*se
    elif alternative == 'greater':
        z = ndtri(confidence_level)
        loglow = log_or - z*se
        loghigh = np.inf
    else:
        # alternative is 'two-sided'
        z = ndtri(0.5*confidence_level + 0.5)
        loglow = log_or - z*se
        loghigh = log_or + z*se

    return np.exp(loglow), np.exp(loghigh)


@dataclass
class ConfidenceInterval:
    low: float
    high: float


@dataclass
class OddsRatioResult:
    """
    Result of `scipy.stats.contingency.odds_ratio`.

    Attributes
    ----------
    table : numpy.ndarray
        The table that was passed to `odds_ratio`.
    kind : str
        The `kind` that was passed to `odds_ratio`. This will be
        either ``'conditional'`` or ``'sample'``.
    alternative : str
        The `alternative` that was passed to `odds_ratio`.  This will
        be ``'two-sided'``, ``'less'`` or ``'greater'``.
    odds_ratio : float
        The computed odds ratio.
    pvalue : float
        The p-value of the estimate of the odds ratio.

    Methods
    -------
    odds_ratio_ci
    """

    table: np.ndarray
    kind: str
    alternative: str
    odds_ratio: float
    pvalue: float

    def odds_ratio_ci(self, confidence_level=0.95):
        """
        Confidence interval for the odds ratio.

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

        Notes
        -----
        When `kind` is ``'conditional'``, the limits of the confidence
        interval are the conditional "exact confidence limits" as described
        by Fisher [1]_. The conditional odds ratio and confidence interval are
        also discussed in Section 4.1.2 of the text by Sahai and Khurshid [2]_.

        When `kind` is ``'sample'``, the confidence interval is computed
        under the assumption that the logarithm of the odds ratio is normally
        distributed with standard error given by::

            se = sqrt(1/a + 1/b + 1/c + 1/d)

        where ``a``, ``b``, ``c`` and ``d`` are the elements of the
        contingency table.  (See, for example, [2]_, section 3.1.3.2,
        or [3]_, section 2.3.3).

        References
        ----------
        .. [1] R. A. Fisher (1935), The logic of inductive inference,
               Journal of the Royal Statistical Society, Vol. 98, No. 1,
               pp. 39-82.
        .. [2] H. Sahai and A. Khurshid (1996), Statistics in Epidemiology:
               Methods, Techniques, and Applications, CRC Press LLC, Boca
               Raton, Florida.
        .. [3] Alan Agresti, An Introduction to Categorical Data Analyis
               (second edition), Wiley, Hoboken, NJ, USA (2007).
        """
        if self.kind == 'conditional':
            ci = self._conditional_odds_ratio_ci(confidence_level)
        else:
            ci = self._sample_odds_ratio_ci(confidence_level)
        return ci

    def _conditional_odds_ratio_ci(self, confidence_level=0.95):
        """
        Confidence interval for the conditional odds ratio.
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

    def _sample_odds_ratio_ci(self, confidence_level=0.95):
        """
        Confidence interval for the sample odds ratio.
        """
        if confidence_level < 0 or confidence_level > 1:
            raise ValueError('confidence_level must be between 0 and 1')

        table = self.table
        if 0 in table.sum(axis=0) or 0 in table.sum(axis=1):
            # If both values in a row or column are zero, the p-value is 1,
            # the odds ratio is NaN and the confidence interval is (0, inf).
            ci = (0, np.inf)
        else:
            ci = _sample_odds_ratio_ci(table,
                                       confidence_level=confidence_level,
                                       alternative=self.alternative)
        return ConfidenceInterval(low=ci[0], high=ci[1])


def odds_ratio(table, kind='conditional', alternative='two-sided'):
    r"""
    Compute the odds ratio for a 2x2 contingency table.

    Parameters
    ----------
    table : array_like of ints
        A 2x2 contingency table.  Elements must be non-negative integers.
    kind : str, optional
        Which kind of odds ratio to compute, either the sample
        odds ratio (``kind='sample'``) or the conditional odds ratio
        (``kind='conditional'``).  Default is ``'conditional'``.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available (default is 'two-sided'):

        * 'two-sided'
        * 'less': one-sided
        * 'greater': one-sided

    Returns
    -------
    result : `~scipy.stats._result_classes.OddsRatioResult` instance
        The returned object has two computed attributes:

        odds_ratio : float
            * If `kind` is ``'sample'``, this is
              ``table[0, 0]*table[1, 1]/(table[0, 1]*table[1, 0])``.
              This is the prior odds ratio and not a posterior estimate.
            * If `kind` is ``'conditional'``, this is the conditional
              maximum likelihood estimate for the odds ratio. It is
              the noncentrality parameter of Fisher's noncentral
              hypergeometric distribution with the same hypergeometric
              parameters as `table` and whose mean is ``table[0, 0]``.
        pvalue : float
            The p-value associated with the computed odds ratio.

            * If `kind` is ``'sample'``, the p-value is based on the
              normal approximation to the distribution of the log of
              the sample odds ratio.
            * If `kind` is ``'conditional'``, the p-value is computed
              by `scipy.stats.fisher_exact`.

        The object also stores the input arguments `table`, `kind`
        and `alternative` as attributes.

        The object has the method `odds_ratio_ci` that computes
        the confidence interval of the odds ratio.

    Notes
    -----
    The conditional odds ratio was discussed by Fisher (see "Example 1"
    of [1]_).  Texts that cover the odds ratio include [2]_ and [3]_.

    References
    ----------
    .. [1] R. A. Fisher (1935), The logic of inductive inference,
           Journal of the Royal Statistical Society, Vol. 98, No. 1,
           pp. 39-82.
    .. [2] Breslow NE, Day NE (1980). Statistical methods in cancer research.
           Volume I - The analysis of case-control studies. IARC Sci Publ.
           (32):5-338. PMID: 7216345. (See section 4.2.)
    .. [3] H. Sahai and A. Khurshid (1996), Statistics in Epidemiology:
           Methods, Techniques, and Applications, CRC Press LLC, Boca
           Raton, Florida.

    Examples
    --------
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
    odds ratio.

    Consider a hypothetical example where it is hypothesized that
    exposure to a certain chemical is assocated with increased occurrence
    of a certain disease.  Suppose we have the following table for a
    collection of 410 people::

                  exposed   unexposed
        cases         7         15
        noncases     58        472

    The question we ask is "Is exposure to the chemical associated with
    increased risk of the disease?"

    Compute the odds ratio:

    >>> from scipy.stats.contingency import odds_ratio
    >>> test = odds_ratio([[7, 15], [58, 472]])
    >>> test.odds_ratio
    3.7836687705553564

    For this sample, the odds of getting the disease for those who have
    been exposed to the chemical are almost 3.8 times that of those who
    have not been exposed.

    Check the p-value:

    >>> test.pvalue
    0.009208708293019454

    The null hypothesis for this test is that the odds ratio is 1.
    The p-value tells us that the probability of getting such an
    extreme odds ratio under the null hypothesis is 0.0092.

    We can compute the 95% confidence interval for the odds ratio:

    >>> test.odds_ratio_ci(confidence_level=0.95)
    ConfidenceInterval(low=1.251482913226682, high=10.363493716701287)

    The 95% confidence interval for the conditional odds ratio is
    approximately (1.25, 10.4), which is evidence that the odds ratio
    is greater than 1.
    """
    if kind not in ['conditional', 'sample']:
        raise ValueError("`kind` must be 'conditional' or 'sample'.")
    if alternative not in ['two-sided', 'less', 'greater']:
        raise ValueError("`alternative` must be 'two-sided', 'less' or "
                         "'greater'.")

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
        # If both values in a row or column are zero, the p-value is NaN and
        # the odds ratio is NaN.
        result = OddsRatioResult(table=c, kind=kind, alternative=alternative,
                                 odds_ratio=np.nan, pvalue=np.nan)
        return result

    if kind == 'sample':
        oddsratio = _sample_odds_ratio(c)
        log_or = np.log(oddsratio)
        se = np.sqrt((1/c).sum())
        if alternative == 'two-sided':
            pvalue = 2*ndtr(-abs(log_or)/se)
        elif alternative == 'less':
            pvalue = ndtr(log_or/se)
        else:
            pvalue = ndtr(-log_or/se)
    else:
        # kind is 'conditional'
        oddsratio = _conditional_oddsratio(c)
        # We can use fisher_exact to compute the p-value.
        pvalue = fisher_exact(c, alternative=alternative)[1]

    result = OddsRatioResult(table=c, kind=kind, alternative=alternative,
                             odds_ratio=oddsratio, pvalue=pvalue)
    return result
