from collections import namedtuple
from itertools import permutations
from functools import lru_cache
import os
import numpy as np
from ._continuous_distns import norm
import scipy.stats


Page_L_Result = namedtuple('Page_L_Result',
                           ('statistic', 'pvalue', 'method'))


def pagel(data, ranked=True, predicted_ranks=None, method='auto', n_s=3000,
          ties=None):
    r"""
    Compute Page's L, a measure of trend in observations between treatments.

    Page's L is useful when:

        * there are :math:`n \geq 3` treatments,
        * :math:`m \geq 2` subjects are observed for each treatment, and
        * the observations are hypothesized to have a particular order.

    Specifically, the test considers the null hypothesis that

    .. math::

        m_1 = m_2 = m_3 \cdots = m_n,

    where :math:`m_j` is the mean of the observed quantity under treatment
    :math:`j`, against the alternative hypothesis that

    .. math::

        m_1 > m_2 > m_3 > \cdots > m_n.

    As noted by [4]_, Page's :math:`L` test has greater statistical power than
    the Friedman test against the alternative that there is a difference in
    trend, as Friedman's test only considers a difference in the means of the
    observations without considering their order. Whereas Spearman :math:`\rho`
    considers the correlation between the ranked observations of two variables
    (e.g. the airspeed velocity of a swallow vs. the weight of the coconut it
    carries), Page's :math:`L` is concerned with a trend in an observation
    (e.g. the airspeed velocity of a swallow) across several distinct
    treatments (e.g. carrying each of five coconuts of different weight) even
    as the observation is repeated with multiple subjects (e.g. one European
    swallow and one African swallow).

    Parameters
    ----------
    data : array-like
        A :math:`m \times n` array; the element in row :math:`i` and
        column :math:`j` is the observation corresponding with subject
        :math:`i` and treatment :math:`j`. By default, the data
        are assumed to be ranked, and the columns are assumed to
        be arranged in order of decreasing predicted mean.

    ranked : boolean, default = ``True``
        By default, the data are assumed to be ranked. If ``False``,
        the ``data`` will be ranked with `scipy.stats.rankdata` along
        ``axis=1``.

    predicted_ranks : array-like, optional
        The predicted ranks of the column means. If not specified,
        the columns are assumed to be arranged in order of decreasing
        predicted mean, so the `predicted_ranks` are
        :math:`[n, n-1, \dots, 2, 1]`.

    method : {'auto', 'asymptotic', 'exact', 'mc'}, optional
        Selects the method used to calculate the *p*-value. The following
        options are available.

            * 'auto': selects between 'exact', 'mc', and 'asymptotic' to
              achieve reasonably accurate results in reasonable time
            * 'asymptotic': compares the standardized test statistic against
              the chi-square distribution with one degree of freedom
            * 'exact': computes the exact *p*-value by comparing the observed
              :math:`L` statistic against those realized by all possible
              permutations of ranks (under the null hypothesis that each
              permutation is equally likely)
            * 'mc': computes the *p*-value by comparing the observed :math:`L`
              statistic against those from randomly generated rank data

    n_s : int, default = 2000
        If ``method='mc'``, `n_s` is the number of Monte Carlo samples
        used to approximate a *p*-value; ignored otherwise.

    ties : boolean, optional
        Whether to adjust the calculation for the possibility of ties. If
        not specified, the calculations are adjusted for the possibility of
        ties only if ties are present in at least one row of the data. Ignored
        when the `'exact'` method is used.

    Returns
    -------
    statistic : float
        Page's :math:`L` test statistic.
    pvalue : float
        The associated *p*-value
    method : {'asymptotic', 'exact', 'mc'}
        The method used to compute the *p*-value

    See Also
    --------
    rankdata, friedmanchisquare, spearmanr

    Notes
    -----
    As noted in [1]_, "the :math:`n` 'treatments' could just as well represent
    :math:`n` objects or events or performances or persons or trials ranked."
    Similarly, the :math:`m` 'subjects' could equally stand for :math:`m`
    "groupings by ability or some other control variable, or judges doing
    the ranking, or random replications of some other sort."

    The procedure for calculating the :math:`L` statistic, adapted from
    [1]_, is:

        1. "Predetermine with careful logic the appropriate hypotheses
           concerning the predicted ording of the experimental results.
           If no reasonable basis for ordering any treatments is known, the
           :math:`L` test is not appropriate."
        2. "As in other experiments, determine at what level of confidence
           you will reject the null hypothesis that there is no agreement of
           experimental results with the monotonic hypothesis."
        3. "Cast the experimental material into a two-way table of :math:`n`
           columns (treatments, objects ranked, conditions) and :math:`m`
           rows (subjects, replication groups, levels of control variables)."
        4. "When experimental observations are recorded, rank them across each
           row", e.g. ``ranks = scipy.stats.rankdata(data, axis=1)``.
        5. "Add the ranks in each column", e.g.
           ``colsums = np.sum(ranks, axis=0)``.
        6. "Multiply each sum of ranks by the predicted rank for that same
           column", e.g. ``products = predicted_ranks * colsums``.
        7. "Sum all such products", e.g. ``L = products.sum()``.

    [1]_ continues by suggesting use of the standardized statistic

    .. math::

        \chi_L^2 = \frac{\left[12L-3mn(n+1)^2\right]^2}{mn^2(n^2-1)(n+1)}

    "which is distributed approximately as chi-square with 1 degree of
    freedom. The ordinary use of :math:`\chi^2` tables would be
    equivalent to a two-sided test of agreement. If a one-sided test
    is desired, *as will almost always be the case*, the probability
    discovered in the chi-square table should be *halved*."

    However, this standardized statistic does not distinguish between the
    observed values being well correlated with the predicted ranks and being
    _anti_-correlated with the predicted ranks. Instead, we follow [2]_
    and calculate the standardized statistic

    .. math::

        \Lambda = \frac{L - E_0}{\sqrt{V_0}},

    where :math:`E_0 = \frac{1}{4} mn(n+1)^2` and
    :math:`V_0 = \frac{1}{144} mn^2(n+1)(n^2-1)`, "which is asymptotically
    normal under the null hypothesis". The asymptotic calculation is not
    adjusted for the possibility of ties, which, according to [2]_, does not
    violate the significance level of the result.

    The *p*-value for ``method='exact'`` is generated by comparing the observed
    value of :math:`L` against the :math:`L` values generated for all
    :math:`(n!)^m` possible permutations of ranks. The calculation is performed
    using the recursive method of [5]. For ``method='mc'``,
    the observed value of :math:`L`  is compared against ``n_s``
    randomly sampled rank permutations. If the calculation is to be adjusted
    for the possibility of ties, then the observed :math:`L` is compared
    against the :math:`L` values generated for all  :math:`(m^m)^n` possible
    tables of ranks (`exact`) or a randomly generated subset of them (`mc`).


    References
    ----------
    .. [1] Ellis Batten Page, "Ordered hypotheses for multiple treatments:
       a significant test for linear ranks", *Journal of the American
       Statistical Association* 58(301), p. 216--230, 1963.

    .. [2] Markus Neuhauser, *Nonparametric Statistical Test: A computational
       approach*, CRC Press, p. 150--152, 2012.

    .. [3] Statext LLC, "Page's L Trend Test - Easy Statistics", *Statext -
       Statistics Study*, https://www.statext.com/practice/PageTrendTest03.php,
       Accessed July 12, 2020.

    .. [4] "Page's Trend Test", *Wikipedia*, WikimediaFoundation,
       https://en.wikipedia.org/wiki/Page%27s_trend_test,
       Accessed July 12, 2020.

    .. [5] Robert E. Odeh, "The exact distribution of Page's L-statistic in
       the two-way layout", *Communications in Statistics - Simulation and
       Computation*,  6(1), p. 49--61, 1977.

    Examples
    --------
    We use the example from [3]_: 10 students are asked to rate three
    teaching methods - tutorial, lecture, and seminar - on a scale of 1-5,
    with 1 being the lowest and 5 being the highest. We have decided that
    a confidence level of 99% is required to reject the null hypothsis in favor
    of our alternative: that the seminar will have the highest ratings and the
    tutorial will have the lowest. Initially, the data have been tabulated with
    each row representing an individual student's ratings of the three methods
    in the following order: tutorial, lecture, seminar.

    >>> table = [[3, 4, 3],
    ...          [2, 2, 4],
    ...          [3, 3, 5],
    ...          [1, 3, 2],
    ...          [2, 3, 2],
    ...          [2, 4, 5],
    ...          [1, 2, 4],
    ...          [3, 4, 4],
    ...          [2, 4, 5],
    ...          [1, 3, 4]]

    We rank the data such that high scores correspond with high ranks, settling
    ties with an average rank:

    >>> from scipy.stats import rankdata
    >>> ranks = rankdata(table, axis=1)
    >>> ranks
    array([[1.5, 3. , 1.5],
           [1.5, 1.5, 3. ],
           [1.5, 1.5, 3. ],
           [1. , 3. , 2. ],
           [1.5, 3. , 1.5],
           [1. , 2. , 3. ],
           [1. , 2. , 3. ],
           [1. , 2.5, 2.5],
           [1. , 2. , 3. ],
           [1. , 2. , 3. ]])

    Because the seminar is hypothesized to have the highest
    ratings and thus the highest ranks, the column corresponding with
    seminar rankings should be first. In fact, according to our alternative
    hypothesis, the order of the columns should be reversed. We pass the
    rearranged table of ranks into `pagel`.

    >>> from scipy.stats import pagel
    >>> arranged_ranks = ranks[:, ::-1]
    >>> res = pagel(arranged_ranks, method="asymptotic")
    >>> res
    Page_L_Result(statistic=133.5, pvalue=0.0012693433690751756, method='asymptotic')

    The value of the :math:`L` statistic, 133.5, is as expected:

    >>> import numpy as np
    >>> m, n = arranged_ranks.shape
    >>> predicted_ranks = np.arange(n, 0, -1)
    >>> L = (predicted_ranks * np.sum(arranged_ranks, axis=0)).sum()
    >>> res.statistic == L
    True

    The *p*-value is the survival function of the normal distribution evaluated
    at the standardized test statistic:

    >>> from scipy.stats import norm
    >>> E0 = (m*n*(n+1)**2)/4
    >>> V0 = (m*n**2*(n+1)*(n**2-1))/144
    >>> Lambda = (L-E0)/np.sqrt(V0)
    >>> p = norm.sf(Lambda)
    >>> res.pvalue == p
    True

    This value indicates that there is a 0.1269% chance that
    the :math:`L` statistic would reach such an extreme value under the null
    hypothesis. Because 0.1269% is less than 1%, we have evidence to reject
    the null hypothesis in favor of our alternative at a 99% confidence level.

    Note that the we can also pass in the data as originally tabulated if we
    also provide ``ranked`` and ``predicted_ranks`` arguments:

    >>> res = pagel(table,                      # data as originally tabulated
    ...             ranked=False,               # originally, data was not ranked
    ...             predicted_ranks=[1, 2, 3],  # originally, data was in order of increasing rank
    ...             method="asymptotic"
    ...             )
    >>> res
    Page_L_Result(statistic=133.5, pvalue=0.0012693433690751756, method='asymptotic')

    As presented in [3]_, the *p*-value was calculated based on the asymptotic
    distribution of the :math:`L` statistic. However, the asymptotic
    distribution is not very accurate, nor conservative, for :math:`m \leq 12`
    and :math:`n \leq 8`. Rather than passing in ``method="asymptotic"``,
    we should leave ``pagel`` to select the method for us.

    >>> res = pagel(table, ranked=False, predicted_ranks=[1, 2, 3])
    >>> res
    Page_L_Result(statistic=133.5, pvalue=0.0018191161948127822, method='exact')

    Note that ``pagel`` chose to use ``method='exact'`` based on the dimensions
    of the table. This exact *p*-value is greater than that given by the
    asymptotic approach, but we can still reject the null hypothesis in favor
    of our alternative at the 99% confidence level.
    """

    # Possible values of the method parameter and the corresponding function
    # used to evaluate the p value
    method = method.lower()
    methods = {"asymptotic": _l_p_asymptotic,
               "mc": lambda L, m, n: _l_p_mc(L, m, n, n_s),
               "exact": _l_p_exact,
               "auto": None}
    if method not in methods:
        raise Exception(f"`method` must be in {set(methods)}")

    # ensure NumPy array and rank the data if it's not already ranked
    if ranked:
        # We're going to trust the user on this. Checking that the data is
        # properly ranked could take as much time as ranking it.
        ranks = np.array(data, copy=False)
    else:
        if np.any(np.isnan(data)):
            raise Exception("`data` contains NaNs, which cannot be ranked meaningfully")
        ranks = scipy.stats.rankdata(data, axis=-1)

    if ranks.ndim != 2:  # TODO: relax this to accept 3d arrays?
        raise Exception(f"`data` must be a 2d array.")
    m, n = ranks.shape
    if m < 2 or n < 3:
        raise Exception("Page's L is only appropriate for data with two "
                        "or more rows and three or more columns.")

    # generate predicted ranks if not provided, ensure valid NumPy array
    if predicted_ranks is None:
        predicted_ranks = np.arange(n, 0, -1)
    else:
        predicted_ranks = np.array(predicted_ranks, copy=False)
        if (set(predicted_ranks) != set(range(1, n+1)) or
                len(predicted_ranks) != n):
            raise Exception(f"`predicted_ranks` must include each integer "
                            f"from 1 to {n} (the number of columns in data) "
                            f"exactly once.")

    if type(ranked) is not bool:
        raise Exception("`ranked` must be boolean.")

    if int(n_s) != n_s:
        raise Exception("`n_s` must be an integer.")

    # Calculate the L statistic
    L = _l_vectorized(ranks, predicted_ranks)

    # Calculate the p-value
    if method == "auto":
        method = _choose_method(ranks)
    p_fun = methods[method]  # get the function corresponding with the method
    p = p_fun(L, m, n)

    return Page_L_Result(L, p, method)


def _choose_method(ranks):
    '''Choose method for computing p-value automatically'''
    m, n = ranks.shape
    if n > 8 or (m > 12 and n > 3) or m > 20:  # as in [1], [4]
        method = "asymptotic"
    else:
        method = "exact"
    return method


def _l_vectorized(ranks, predicted_ranks):
    '''Calculate's Page's L statistic for each page of a 3d array'''
    colsums = ranks.sum(axis=-2, keepdims=True)
    products = predicted_ranks * colsums
    Ls = products.sum(axis=-1)
    Ls = Ls[0] if Ls.size == 1 else Ls.ravel()
    return Ls


def _l_p_asymptotic(L, m, n):
    '''Calculate the p-value of Page's L from the asymptotic distribution'''
    # Using [1] as a reference, the asymptotic p-value would be calculated as:
    # chi_L = (12*L - 3*m*n*(n+1)**2)**2/(m*n**2*(n**2-1)*(n+1))
    # p = chi2.sf(chi_L, df=1, loc=0, scale=1)/2
    # but this is insentive to the direction of the hypothesized ranking

    # See [2] page 151
    E0 = (m*n*(n+1)**2)/4
    V0 = (m*n**2*(n+1)*(n**2-1))/144
    Lambda = (L-E0)/np.sqrt(V0)
    p = norm.sf(Lambda)
    return p


def _l_p_mc(L, m, n, n_s):
    '''Calculate the p-value of Page's L from a Monte Carlo distribution'''
    mc_data = np.random.rand(n_s, m, n)
    mc_ranks = scipy.stats.rankdata(mc_data, axis=-1)
    predicted_ranks = np.arange(n, 0, -1)
    Ls = _l_vectorized(mc_ranks, predicted_ranks)
    p = 1 - scipy.stats.percentileofscore(Ls, L)/100
    return p


def _l_p_exact(L, m, n):
    L, n, k = int(L), int(m), int(n)  # different papers use different symbols
    a, b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6

#    dir_path = os.path.dirname(os.path.realpath(__file__))
    dir_path = os.path.abspath(os.path.dirname(__file__))
    datafile = os.path.join(dir_path, "data_folder", "pagel_exact.npy")
    all_pmfs = np.load(datafile, allow_pickle=True).item()
    try:
        pmf = all_pmfs[n][k]
        return np.sum(pmf[L-a*n:])
    except KeyError:
        return _sf(L, k, n)


def _counts_l_k_1(k):
    '''Count frequency of each L value over all possible single rows'''
    # See [5] Equation (6)
    ranks = range(1, k+1)
    rank_perms = np.array(list(permutations(ranks)))  # generate all possible rows
    Ls = (ranks*rank_perms).sum(axis=1)  # compute Page's L for all rows
    a, b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6  # min and max possible L
    counts = np.histogram(Ls, np.arange(a-0.5, b+1.5))[0]  # count occurences of each L value
    return counts  # not saving the corresponding L values [a, b], just counts


def _generate_row_counts(k_max):
    '''Generate PMF data for hard-coding'''
    # The results of this function are copy-pasted into the code
    all_counts = [list(_counts_l_k_1(k)) for k in range(1, k_max+1)]
    print(repr(all_counts))

# _generate_row_counts(9)


def _get_p_l_k_1(k):
    '''Generates function to evaluate p(t, k, 1); see [5] Equation 6'''
    # Data in the following row was generated using _generate_row_counts
    all_counts = [[1], [1, 1], [1, 2, 0, 2, 1], [1, 3, 1, 4, 2, 2, 2, 4, 1, 3, 1], [1, 4, 3, 6, 7, 6, 4, 10, 6, 10, 6, 10, 6, 10, 4, 6, 7, 6, 3, 4, 1], [1, 5, 6, 9, 16, 12, 14, 24, 20, 21, 23, 28, 24, 34, 20, 32, 42, 29, 29, 42, 32, 20, 34, 24, 28, 23, 21, 20, 24, 14, 12, 16, 9, 6, 5, 1], [1, 6, 10, 14, 29, 26, 35, 46, 55, 54, 74, 70, 84, 90, 78, 90, 129, 106, 123, 134, 147, 98, 168, 130, 175, 144, 168, 144, 184, 144, 168, 144, 175, 130, 168, 98, 147, 134, 123, 106, 129, 90, 78, 90, 84, 70, 74, 54, 55, 46, 35, 26, 29, 14, 10, 6, 1], [1, 7, 15, 22, 47, 54, 70, 94, 129, 124, 178, 183, 237, 238, 276, 264, 379, 349, 380, 400, 517, 394, 542, 492, 640, 557, 666, 595, 776, 684, 786, 718, 922, 745, 917, 781, 982, 826, 950, 844, 1066, 845, 936, 845, 1066, 844, 950, 826, 982, 781, 917, 745, 922, 718, 786, 684, 776, 595, 666, 557, 640, 492, 542, 394, 517, 400, 380, 349, 379, 264, 276, 238, 237, 183, 178, 124, 129, 94, 70, 54, 47, 22, 15, 7, 1], [1, 8, 21, 34, 72, 102, 130, 190, 260, 284, 398, 454, 555, 616, 756, 744, 1022, 1042, 1159, 1282, 1555, 1392, 1719, 1758, 2009, 2032, 2282, 2214, 2676, 2590, 2878, 2928, 3397, 3138, 3647, 3568, 3921, 3866, 4311, 4050, 4852, 4492, 4816, 4784, 5505, 4954, 5638, 5304, 5890, 5486, 6188, 5502, 6436, 5822, 6233, 6024, 6697, 5720, 6672, 6020, 6688, 6020, 6672, 5720, 6697, 6024, 6233, 5822, 6436, 5502, 6188, 5486, 5890, 5304, 5638, 4954, 5505, 4784, 4816, 4492, 4852, 4050, 4311, 3866, 3921, 3568, 3647, 3138, 3397, 2928, 2878, 2590, 2676, 2214, 2282, 2032, 2009, 1758, 1719, 1392, 1555, 1282, 1159, 1042, 1022, 744, 756, 616, 555, 454, 398, 284, 260, 190, 130, 102, 72, 34, 21, 8, 1]]
    try:
        counts = all_counts[k-1]
    except IndexError:
        counts = _counts_l_k_1(k)
    ps = np.array(counts)/np.math.factorial(k)  # See [5] Equation (6)
    a, b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6  # min L, max L for a row

    def p_l_k_1(l):
        if l < a or l > b:  # 0 rows have L < a or L > b
            return 0
        return ps[l-a]

    return p_l_k_1


@lru_cache(maxsize=None)
def _pmf_recursive(l, k, n, a, b, p_l_k_1):
    '''Recursive function to evaluate p(l, k, n); see [5] Equation 1'''
    # If we care about users generating exact p-values on the fly,
    # this should look up whether p(*, k, n-1) is already tabulated
    # in pagel_exact.npy
    if n == 1:
        p = p_l_k_1(l)
        return(p)
    else:
        p = 0
        low = max(l-(n-1)*b, a)
        high = min(l-(n-1)*a, b)
        for t in range(low, high+1):
            p1 = _pmf_recursive(l-t, k, n-1, a, b, p_l_k_1)
            p2 = _pmf_recursive(t, k, 1, a, b, p_l_k_1)
            p += p1*p2
        return p


def _pmf(l, k, n):
    '''Probability mass function of Page's L statistic'''
    a, b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6
    p_l_k_1 = _get_p_l_k_1(k)
    return _pmf_recursive(l, k, n, a, b, p_l_k_1)


def _sf(l, k, n):
    '''Survival function of Page's L statistic'''
    a, b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6
    ps = [_pmf(l, k, n) for l in range(l, n*b + 1)]
    return np.sum(ps)


def _whole_pmf(k, n):
    '''Return list of all values of Page's L PMF'''
    a, b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6
    ps = [_pmf(l, k, n) for l in range(a*n, n*b + 1)]
    return ps


def _all_pmf():
    '''Generate all PMFs tabulated in [1]'''
    data = {}
    for n in range(2, 13):
        data[n] = {}
        for k in range(3, 9):
            print(n, k)
            data[n][k] = _whole_pmf(k, n)

    k = 3
    for n in range(13, 21):
        data[n] = {}
        data[n][k] = _whole_pmf(k, n)

    np.save('pagel_exact', data)
