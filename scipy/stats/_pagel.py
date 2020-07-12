from collections import namedtuple
import numpy as np
from ._continuous_distns import chi2


Page_L_Result = namedtuple('Page_L_Result',
                           ('statistic', 'pvalue', 'method'))


def pagel(data, ranked=True, hypothesis=None, method='auto', samples=2000,
          ties=None):
    r"""
    Compute Page's L, a measure of trend in observations between treatments.

    Page's L is useful when:

        * there are :math:`n \geq 3` treatments,
        * :math:`k \geq 2` subjects are observed for each treatment, and
        * the observations are hypothesized to have a particular order.

    Specifically, the test considers the null hypothesis that

    .. math::

        m_1 = m_2 = m_3 \cdots = m_n,

    where :math:`m_j` is the mean of the observed quantity under treatment
    :math:`j`, against the alternative hypothesis that

    .. math::

        m_1 > m_2 > m_3 > \cdots > m_n.

    Parameters
    ----------
    data : array-like
        A :math:`k \times n` array; the element in row :math:`i` and
        column :math:`j` is the observation corresponding with subject
        :math:`i` and treatment :math:`j`. By default, the data
        is assumed to be ordinal or ranked, and the columns are assumed to
        be arranged in order of decreasing hypothesized mean.

    ranked : boolean, default = ``True``
        By default, the data is assumed to be ordinal or ranked. If ``False``,
        the ``data`` will be ranked with `scipy.stats.rankdata` along
        ``axis=1``.

    hypothesis : array-like, optional
        The hypothesized ordering of the column means; if specified the columns
        of ``data`` will be rearranged like ``data[:, hypothesis]``.

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
    rankdata, friedmanchisquare

    Notes
    -----
    As noted in [1]_, "the :math:`n` 'treatments' could just as well represent
    :math:`n` objects or events or performances or persons or trials ranked."
    Similarly, the :math:`k` 'subjects' could equally stand for :math:`k`
    "groupings by ability or some other control variable, or judges doing
    the ranking, or random replications of some other sort."

    The procedure for performing the asymptotic :math:`L` test, adapted from
    [1]_, is:

        1. Predetermine with careful logic the appropriate hypotheses
           concerning the predicted ording of the experimental results.
           If no reasonable basis for ordering any treatments is known, the
           :math:`L` test is not appropriate.
        2. As in other experiments, determine at what level of confidence
           you will reject the null hypothesis that there is no agreement of
           experimental results with the monotonic hypothesis.
        3. Cast the experimental material into a two-way table of :math:`n`
           columns (treatments, objects ranked, conditions) and :math:`m`
           rows (subjects, replication groups, levels of control variables).
        4. When experimental observations are recorded, rank them across each
           row, e.g. ``data = scipy.stats.rankdata(data, axis=1)``.
        5. Add the ranks in each column, e.g. ``a = np.sum(data, axis=0)``.
        6. Multiply each sum of ranks by the predicted rank for that same
           column, e.g. ``b = (hypothesis+1) * a``.
        7. Sum all such products, e.g. ``c = b.sum()``.
        8. For an asymptotic test, use the statistic

        .. math::

           \chi_L^2 = \frac{\left[12L-3kn(n+1)^2\right]^2}{kn^2(n^2-1)(n+1)}

           which is distributed approximately as chi-square with 1 degree of
           freedom. The ordinary use of :math:`\chi^2` tables would be
           equivalent to a two-sided test of agreement. If a one-sided test
           is desired, *as will almost always be the case*, the probability
           discovered in the chi-square table should be *halved*.

    Here, the *p*-value always corresponds with a one-sided test. The
    asymptotic calculation is not adjusted for the possibility of ties,
    which, according to [2]_, does not violate the significance level of the
    result.

    The *p*-value for ``method='exact'`` is generated by comparing the observed
    value of :math:`L` against the :math:`L` values generated for all
    :math:`(k!)^n` possible permutations of ranks. For ``method='mc'``,
    the observed value of :math:`L`  is compared against ``n_s``
    randomly selected rank permutations. If the calculation is to be adjusted
    for the possibility of ties, then the observed :math:`L` is compared
    against the :math:`L` values generated for all  :math:`(k^k)^n` possible
    tables of ranks (`exact`) or a randomly generated subset of them (`mc`).


    References
    ----------
    .. [1] Ellis Batten Page, "Ordered hypotheses for multiple treatments:
       a significant test for linear ranks", *Journal of the American
       Statistical Association* 58(301), p. 216--230, 1963.

    .. [2] Markus Neuhauser, *Nonparametric Statistical Test: A computational
       approach*, CRC Press, p. 150--152, 2012.

    """

    return Page_L_Result(None, None, None)
