import math
import numpy as np
from functools import partial
from scipy.stats._resampling import permutation_test


def _bws_input_validation(x, y, alternative, variant):
    ''' Input validation and standardization for bws test'''
    x, y = np.atleast_1d(x), np.atleast_1d(y)
    if np.isnan(x).any() or np.isnan(y).any():
        raise ValueError('`x` and `y` must not contain NaNs.')
    if np.size(x) == 0 or np.size(y) == 0:
        raise ValueError('`x` and `y` must be of nonzero size.')

    alternatives = {'two-sided', 'less', 'greater'}
    alternative = alternative.lower()
    if alternative not in alternatives:
        raise ValueError(f'`alternative` must be one of {alternatives}.')

    variants = {1, 2, 3, 4, 5}
    if variant not in variants:
        raise ValueError(f'`variant` must be one of {variants}.')

    return x, y, alternative, variant


def _bws_statistic(x, y, alternative, variant):
    '''Compute the BWS test statistic for two independent samples'''
    Ri, Hj = np.sort(x), np.sort(y)
    n, m = Ri.shape[-1], Hj.shape[-1]
    i, j = np.arange(1, n+1), np.arange(1, m+1)

    if alternative == 'two-sided':
        match variant:
            case 1:
                Bx_num = (Ri - (m + n + 1) / (n + 1) * i) ** 2
                By_num = (Hj - (m + n + 1) / (m + 1) * j) ** 2
                Bx_den = i / (n + 1) * (1 - i / (n + 1)) * \
                    m * (m + n + 1) / (n + 2)
                By_den = j / (m + 1) * (1 - j / (m + 1)) * \
                    n * (m + n + 1) / (m + 2)
            case 2:
                Bx_num = (Ri - (m + n) / n * i) ** 2
                By_num = (Hj - (m + n) / m * j) ** 2
                Bx_den = i / (n + 1) * (1 - i / (n + 1)) * m * (m + n) / n
                By_den = j / (m + 1) * (1 - j / (m + 1)) * n * (m + n) / m
            case 3:
                Bx_num = (Ri - (m + n + 1) / (n + 1) * i) ** 2
                By_num = (Hj - (m + n + 1) / (m + 1) * j) ** 2
                Bx_den = (i / (n + 1) * (1 - i / (n + 1))
                          * m * (m + n + 1) / (n + 2)) ** 2
                By_den = (j / (m + 1) * (1 - j / (m + 1))
                          * n * (m + n + 1) / (m + 2)) ** 2
            case 4:
                Bx_num = abs(Ri - (m + n + 1) / (n + 1) * i)
                By_num = abs(Hj - (m + n + 1) / (m + 1) * j)
                Bx_den = (i / (n + 1) * (1 - i / (n + 1))
                          * m * (m + n + 1) / (n + 2)) ** 2
                By_den = (j / (m + 1) * (1 - j / (m + 1))
                          * n * (m + n + 1) / (m + 2)) ** 2
            case 3:
                Bx_num = (Ri - (m + n + 1) / (n + 1) * i) ** 2
                By_num = (Hj - (m + n + 1) / (m + 1) * j) ** 2
                Bx_den = math.log(i / (n + 1) * (1 - i / (n + 1))
                                  * m * (m + n + 1) / (n + 2))
                By_den = math.log(j / (m + 1) * (1 - j / (m + 1))
                                  * n * (m + n + 1) / (m + 2))

        Bx = 1 / n * np.sum(Bx_num / Bx_den)
        By = 1 / m * np.sum(By_num / By_den)

        B = (Bx + By) * 0.5
    else:
        match variant:
            case 1:
                Bx_num = (Ri - (m + n + 1) / (n + 1) * i) * \
                    abs(Ri - (m + n + 1) / (n + 1) * i)
                By_num = (Hj - (m + n + 1) / (m + 1) * j) * \
                    abs(Hj - (m + n + 1) / (m + 1) * j)
                Bx_den = i / (n + 1) * (1 - i / (n + 1)) * \
                    m * (m + n + 1) / (n + 2)
                By_den = j / (m + 1) * (1 - j / (m + 1)) * \
                    n * (m + n + 1) / (m + 2)
            case 2:
                Bx_num = (Ri - (m + n) / n * i) * abs(Ri - (m + n) / n * i)
                By_num = (Hj - (m + n) / m * j) * abs(Hj - (m + n) / m * j)
                Bx_den = i / (n + 1) * (1 - i / (n + 1)) * m * (m + n) / n
                By_den = j / (m + 1) * (1 - j / (m + 1)) * n * (m + n) / m
            case 3:
                Bx_num = (Ri - (m + n + 1) / (n + 1) * i) * \
                    abs(Ri - (m + n + 1) / (n + 1) * i)
                By_num = (Hj - (m + n + 1) / (m + 1) * j) * \
                    abs(Hj - (m + n + 1) / (m + 1) * j)
                Bx_den = (i / (n + 1) * (1 - i / (n + 1))
                          * m * (m + n + 1) / (n + 2)) ** 2
                By_den = (j / (m + 1) * (1 - j / (m + 1))
                          * n * (m + n + 1) / (m + 2)) ** 2
            case 4:
                Bx_num = abs(Ri - (m + n + 1) / (n + 1) * i)
                By_num = abs(Hj - (m + n + 1) / (m + 1) * j)
                Bx_den = (i / (n + 1) * (1 - i / (n + 1))
                          * m * (m + n + 1) / (n + 2)) ** 2
                By_den = (j / (m + 1) * (1 - j / (m + 1))
                          * n * (m + n + 1) / (m + 2)) ** 2
            case 5:
                Bx_num = (Ri - (m + n + 1) / (n + 1) * i) * \
                    abs(Ri - (m + n + 1) / (n + 1) * i)
                By_num = (Hj - (m + n + 1) / (m + 1) * j) * \
                    abs(Hj - (m + n + 1) / (m + 1) * j)
                Bx_den = math.log(i / (n + 1) * (1 - i / (n + 1))
                                  * m * (m + n + 1) / (n + 2))
                By_den = math.log(j / (m + 1) * (1 - j / (m + 1))
                                  * n * (m + n + 1) / (m + 2))

        Bx = 1 / n * np.sum(Bx_num / Bx_den)
        By = 1 / m * np.sum(By_num / By_den)

        if alternative == 'greater':
            B = (Bx - By) * 0.5
        else:
            B = (By - Bx) * 0.5

    return B


def bws_test(x, y, *, alternative="two-sided", variant=1, n_resamples=9999, random_state=None):
    r'''Perform the Baumgartner-Weiss-Schindler test on two independent samples.

    The Baumgartner-Weiss-Schindler (BWS) test is a nonparametric test of 
    the null hypothesis that the distribution underlying sample `x` 
    is the same as the distribution underlying sample `y`. Unlike 
    the Kolmogorov-Smirnov, Wilcoxon, and Cramer-Von Mises tests, 
    the BWS test weights the integral by the variance of the difference
    in CDFs, emphasizing the tails of the distributions, which increases
    the power of the test for a lot of applications.

    Parameters
    ----------
    x, y : array-like
        1-d arrays of samples.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        Let *F(u)* and *G(u)* be the cumulative distribution functions of the
        distributions underlying `x` and `y`, respectively. Then the following
        alternative hypotheses are available:

        * 'two-sided': the distributions are not equal, i.e. *F(u) ≠ G(u)* for
          at least one *u*.
        * 'less': the distribution underlying `x` is stochastically less
          than the distribution underlying `y`, i.e. *F(u) >= G(u)* for all *u*.
        * 'greater': the distribution underlying `x` is stochastically greater
          than the distribution underlying `y`, i.e. *F(u) <= G(u)* for all *u*.

        Under a more restrictive set of assumptions, the alternative hypotheses
        can be expressed in terms of the locations of the distributions;
        see [3] section 5.1.
    variant : {1, 2, 3, 4, 5}, optional
        Selects the method used to calculate the :math:`B` statistic.
        Default is 'naive'. The following options are available.

        * ``1``: Murakami's `B1` statistic, from his 2012 paper [2]_.
        * ``2``: Murakami's `B2` statistic, from his 2012 paper [2]_, which use the same method as [1]_.
        * ``3``: Murakami's `B3` statistic, from his 2012 paper [2]_.
        * ``4``: Murakami's `B4` statistic, from his 2012 paper [2]_.
        * ``5``: Murakami's `B5` statistic, from his 2012 paper [2]_.
    n_resamples : int or np.inf, default: 9999
        Number of random permutations (resamples) used to approximate the null
        distribution. If greater than or equal to the number of distinct
        permutations, the exact null distribution will be computed.
        Note that the number of distinct permutations grows very rapidly with
        the sizes of samples, so exact tests are feasible only for very small
        data sets.
    random_state : {None, int, `numpy.random.Generator`,
                    `numpy.random.RandomState`}, optional

        Pseudorandom number generator state used to generate permutations.

        If `random_state` is ``None`` (default), the
        `numpy.random.RandomState` singleton is used.
        If `random_state` is an int, a new ``RandomState`` instance is used,
        seeded with `random_state`.
        If `random_state` is already a ``Generator`` or ``RandomState``
        instance then that instance is used.

    Returns
    -------
    res : PermutationTestResult
    An object with attributes:

    statistic : float or ndarray
        The observed test statistic of the data.
    pvalue : float or ndarray
        The p-value for the given alternative.
    null_distribution : ndarray
        The values of the test statistic generated under the null hypothesis.

    See also
    --------
    scipy.stats.wilcoxon, scipy.stats.mannwhitneyu, scipy.stats.ttest_ind

    References
    ----------
    .. [1] Neuhäuser, M. (2005). Exact Tests Based on the Baumgartner-Wei-Schindler 
           StatisticA Survey. Statistical Papers, 46(1), 1-29.
    .. [2] Murakami, H. (2012). Modified Baumgartner statistics for the two-sample
           and multisample problems: a numerical comparison. Journal of Statistical
           Computation and Simulation, 82(5), 711-728.
    .. [3] Fay, M. P., & Proschan, M. A. (2010). Wilcoxon-Mann-Whitney or t-test? 
           On assumptions for hypothesis tests and multiple interpretations of 
           decision rules. Statistics surveys, 4, 1.

    Examples
    --------
    We follow the example of table 3 in [1]_: Fourteen children were divided randomly into two groups, 
    an experimental group and a control group.

    >>> control_group_ranks = np.array([1, 2, 3, 4, 6, 7, 8])
    >>> experimental_group_ranks = np.array([5, 9, 10, 11, 12, 13, 14])

    We use the BWS test to assess whether there is a 
    statistically significant difference in the methods of training preschool children.
    The null hypothesis is that the distribution of first method is the same as 
    the distribution of second one. We decide that a confidence level of 95% is 
    required to reject the null hypothesis in favor of the alternative that 
    the distributions are different.
    Since the number of samples is very small and there are no ties in the
    data, we can compare the observed test statistic against the *exact*
    distribution of the test statistic under the null hypothesis.

    >>> from scipy.stats import bws_test
    >>> B, p = bws_test(control_group_ranks, experimental_group_ranks)
    >>> print(B)
    5.132167152575315

    This agrees with :math:`B = 5.132` reported in [1]_. The *p*-value produced
    by `bws_test` also agrees with :math:`p = 0.0029` reported in [1]_.

    >>> print(p)
    0.002913752913752914

    It normally takes a long time to calculate the asymptotic *p*-value when the passing
    arrays are large, so we can approximate *p*-value using a resampling permutation test.
    Here we use `n_resamples=10000` to perform a random permutation test with 10000 samples.
    Note that *p*-value varies depending on the number of samples.

    >>> B, p = bws_test(control_group_ranks, experimental_group_ranks, n_resamples=10000)
    >>> print(p)
    0.003275267196716418

    Under this assumption, the *p*-value would be low enough to reject the
    null hypothesis in favor of the alternative.
    '''

    x, y, alternative, variant = _bws_input_validation(
        x, y, alternative, variant)
    bws_statistic = partial(
        _bws_statistic, alternative=alternative, variant=variant)
    res = permutation_test((x, y), bws_statistic, alternative='greater',
                           n_resamples=n_resamples, random_state=random_state)

    return res
