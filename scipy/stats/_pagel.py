from collections import namedtuple
from itertools import permutations
from functools import lru_cache
import os
import numpy as np
from ._continuous_distns import norm
import scipy.stats


def pagel(data, ranked=True, predicted_ranks=None, method='auto'):
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

    method : {'auto', 'asymptotic', 'exact'}, optional
        Selects the method used to calculate the *p*-value. The following
        options are available.

            * 'auto': selects between 'exact' and 'asymptotic' to
              achieve reasonably accurate results in reasonable time
            * 'asymptotic': compares the standardized test statistic against
              the chi-square distribution with one degree of freedom
            * 'exact': computes the exact *p*-value by comparing the observed
              :math:`L` statistic against those realized by all possible
              permutations of ranks (under the null hypothesis that each
              permutation is equally likely)

    Returns
    -------
    res : StatsTestResult
        A :class:`scipy.stats.StatsTestResult` consisting of the fields:

            statistic : float
                Page's :math:`L` test statistic.
            pvalue : float
                The associated *p*-value
            method : {'asymptotic', 'exact'}
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
    normal under the null hypothesis".

    The *p*-value for ``method='exact'`` is generated by comparing the observed
    value of :math:`L` against the :math:`L` values generated for all
    :math:`(n!)^m` possible permutations of ranks. The calculation is performed
    using the recursive method of [5].

    The *p*-values are not adjusted for the possibility of ties. When
    ties are present, the reported  ``'exact'`` *p*-values may be somewhat
    larger (i.e. more conservative) than the true *p*-value [2]_. The
    ``'asymptotic'``` *p*-values, however, tend to be smaller (i.e. less
    conservative) than the ``'exact'`` *p*-values.


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
    >>> res = pagel(arranged_ranks, method='asymptotic')
    >>> res
        method: 'asymptotic'
        pvalue: 0.0012693433690751756
     statistic: 133.5

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
        method: 'asymptotic'
        pvalue: 0.0012693433690751756
     statistic: 133.5

    As presented in [3]_, the *p*-value was calculated based on the asymptotic
    distribution of the :math:`L` statistic. However, the asymptotic
    distribution is not very accurate, nor conservative, for :math:`m \leq 12`
    and :math:`n \leq 8`. Rather than passing in ``method="asymptotic"``,
    we can allow ``pagel`` to select the appropriate method.

    >>> res = pagel(table, ranked=False, predicted_ranks=[1, 2, 3])
    >>> res
        method: 'exact'
        pvalue: 0.0018191161948127822
     statistic: 133.5

    Note that ``pagel`` chose to use ``method='exact'`` based on the dimensions
    of the table and the recommendations in Page's original paper [1]_.
    This exact *p*-value is greater than that given by the
    asymptotic approach, but we can still reject the null hypothesis in favor
    of our alternative at the 99% confidence level.
    """

    # Possible values of the method parameter and the corresponding function
    # used to evaluate the p value
    method = method.lower()
    methods = {"asymptotic": _l_p_asymptotic,
               "exact": _l_p_exact,
               "auto": None}
    if method not in methods:
        raise ValueError(f"`method` must be in {set(methods)}")

    # ensure NumPy array and rank the data if it's not already ranked
    if ranked:
        # We're going to trust the user on this. Checking that the data is
        # properly ranked could take as much time as ranking it.
        ranks = np.array(data, copy=False)
    else:
        if np.any(np.isnan(data)):
            raise ValueError("`data` contains NaNs, which cannot be ranked meaningfully")
        ranks = scipy.stats.rankdata(data, axis=-1)

    if ranks.ndim != 2:  # TODO: relax this to accept 3d arrays?
        raise ValueError(f"`data` must be a 2d array.")
    m, n = ranks.shape
    if m < 2 or n < 3:
        raise ValueError("Page's L is only appropriate for data with two "
                        "or more rows and three or more columns.")

    # generate predicted ranks if not provided, ensure valid NumPy array
    if predicted_ranks is None:
        predicted_ranks = np.arange(n, 0, -1)
    else:
        predicted_ranks = np.array(predicted_ranks, copy=False)
        if (predicted_ranks.ndim < 1 or
                (set(predicted_ranks) != set(range(1, n+1)) or
                len(predicted_ranks) != n)):
            raise ValueError(f"`predicted_ranks` must include each integer "
                            f"from 1 to {n} (the number of columns in data) "
                            f"exactly once.")

    if type(ranked) is not bool:
        raise TypeError("`ranked` must be boolean.")

    # Calculate the L statistic
    L = _l_vectorized(ranks, predicted_ranks)

    # Calculate the p-value
    if method == "auto":
        method = _choose_method(ranks)
    p_fun = methods[method]  # get the function corresponding with the method
    p = p_fun(L, m, n)

    pagel_result = scipy.stats.StatsTestResult()
    pagel_result['statistic'] = L
    pagel_result['pvalue'] = p
    pagel_result['method'] = method
    return pagel_result


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


def _l_p_exact(L, m, n):
    L, n, k = int(L), int(m), int(n)  # different papers use different symbols

    _pagel_state.set_k(k)
    try:
        pmf = _pagel_state.all_pmfs[n][k]
        return np.sum(pmf[L-_pagel_state.a*n:])
    except KeyError:
        print('calculating')
        return _pagel_state.sf(L, n)


class _PageL:
    '''Class to maintain state of Page's L between `pagel` executions'''
    def __init__(self):
        self.all_pmfs = None
        self.all_counts = None
        self.k = None
        self.a = None
        self.b = None
        self.p_l_k_1 = None
        dir_path = os.path.abspath(os.path.dirname(__file__))
        self.all_pmfs_file = os.path.join(dir_path, "data_folder",
                                          "pagel_pmfs.npy")
        self.all_counts_file = os.path.join(dir_path, "data_folder",
                                            "pagel_counts.npy")

    def load_pagel_data(self):
        '''Loads PageL data from file'''
        # functions for creating these files are at the end of the class
        self.all_counts = np.load(self.all_counts_file, allow_pickle=True)
        self.all_pmfs = np.load(self.all_pmfs_file, allow_pickle=True).item()

    def set_k(self, k):
        '''Generates function to evaluate p(t, k, 1); see [5] Equation 6'''

        if self.k == k:
            return

        if self.all_counts is None:
            self.load_pagel_data()

        self.k = k
        # min L, max L of a row
        self.a, self.b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6

        try:
            counts = self.all_counts[k-1]
        except IndexError:
            counts = self.counts_l_k_1(k)

        # See [5] Equation (6)
        ps = np.array(counts)/np.math.factorial(k)

        def p_l_k_1(l):
            if l < self.a or l > self.b:  # 0 rows have L < a or L > b
                return 0
            return ps[l-self.a]

        self.p_l_k_1 = p_l_k_1

    def pmf(self, l, n):
        '''Probability mass function of Page's L statistic'''
        return _pmf_recursive(l, self.k, n, self.a, self.b, self.p_l_k_1)

    def sf(self, l, n):
        '''Survival function of Page's L statistic'''
        ps = [self.pmf(l, n) for l in range(l, n*self.b + 1)]
        return np.sum(ps)

    def whole_pmf(self, k, n):
        '''Return list of all values of Page's L PMF'''
        self.set_k(k)
        ps = [self.pmf(l, n) for l in range(self.a*n, self.b*n + 1)]
        return ps

    def all_pmf(self):
        '''Generate all PMFs tabulated in [1]'''
        data = {}
        for n in range(2, 13):
            data[n] = {}
            for k in range(3, 9):
                print(n, k)
                data[n][k] = self.whole_pmf(k, n)

        k = 3
        for n in range(13, 21):
            data[n] = {}
            data[n][k] = self.whole_pmf(k, n)

        np.save(self.all_pmfs_file, data)

    def generate_row_counts(self, k_max):
        '''Generate PMF data for hard-coding'''
        all_counts = [self.counts_l_k_1(k) for k in range(1, k_max+1)]
        np.save(self.all_counts_file, all_counts)

    def counts_l_k_1(self, k):
        '''Count frequency of each L value over all possible single rows'''
        # See [5] Equation (6)
        ranks = range(1, k+1)
        # generate all possible rows of length k
        rank_perms = np.array(list(permutations(ranks)))
        # compute Page's L for all possible rows
        Ls = (ranks*rank_perms).sum(axis=1)
        # min and max possible L for row of length k
        a, b = (k*(k+1)*(k+2))//6, (k*(k+1)*(2*k+1))//6
        # count occurences of each L value
        counts = np.histogram(Ls, np.arange(a-0.5, b+1.5))[0]
        # not saving the corresponding L values [a, b], just counts
        return counts


# left out of class to ensure that object state does not interfere with caching
@lru_cache(maxsize=None)
def _pmf_recursive(l, k, n, a, b, p_l_k_1):
    '''Recursive function to evaluate p(l, k, n); see [5] Equation 1'''
    # If we care about users generating exact p-values on the fly,
    # this should look up whether p(*, k, n-1) is already tabulated
    # in pagel_exact.npy and set p1 = all_pmfs[n-1][k][l-t-a*(n-1)]
    # Put all these functions in a PageLExact class, too, to avoid
    # passing around all these parameters.
    if n == 1:
        p = p_l_k_1(l)
        return(p)

    p = 0
    low = max(l-(n-1)*b, a)
    high = min(l-(n-1)*a, b)

    for t in range(low, high+1):
        p1 = _pmf_recursive(l-t, k, n-1, a, b, p_l_k_1)
        p2 = _pmf_recursive(t, k, 1, a, b, p_l_k_1)
        p += p1*p2
    return p


# Fast to instantiate, and only loads data from file once
# (on first use of `pagel` with `method='exact')
_pagel_state = _PageL()
