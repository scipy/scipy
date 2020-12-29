import numpy as np
from dataclasses import make_dataclass
from collections import namedtuple
from scipy import special
from scipy import stats


def _broadcast_concatenate(x, y, axis):
    '''Broadcast then concatenate arrays, leaving concatenation axis last'''
    x = x.swapaxes(axis, -1)
    y = y.swapaxes(axis, -1)
    z = np.broadcast(x[..., 0], y[..., 0])
    x = np.broadcast_to(x, z.shape + (x.shape[-1],))
    y = np.broadcast_to(y, z.shape + (y.shape[-1],))
    z = np.concatenate((x, y), axis = -1)
    return x, y, z


class _MWU:
    '''Distribution of MWU statistic under the null hypothesis'''
    # Possible improvement: if m and n are small enough, use integer arithmetic

    def __init__(self):
        '''Minimal initializer'''
        self._fmnks = -np.ones((1, 1, 1))

    def pmf(self, k, m, n):
        '''Probability mass function'''
        self._resize_fmnks(m, n, np.max(k))
        # could loop over just the unique elements, but probably not worth
        # the time to find them
        for i in np.ravel(k):
            self._f(m, n, i)
        return self._fmnks[m, n, k] / special.binom(m + n, m)

    def cdf(self, k, m, n):
        '''Cumulative distribution function'''
        # We could use the fact that the distribution is symmetric to avoid
        # summing more than m*n/2 terms, but it might not be worth the
        # overhead. Let's leave that to an improvement.
        pmfs = self.pmf(np.arange(0, np.max(k) + 1), m, n)
        cdfs = np.cumsum(pmfs)
        return cdfs[k]

    def sf(self, k, m, n):
        '''Survival function'''
        # Use the fact that the distribution is symmetric; i.e.
        # _f(m, n, m*n-k) = _f(m, n, k), and sum from the left
        k = m*n - k
        # Note that both CDF and SF include the PMF at k. The p-value is
        # calculated from the SF and should include the mass at k, so this
        # is desirable
        return self.cdf(k, m, n)

    def _resize_fmnks(self, m, n, k):
        '''If necessary, expand the array that remembers PMF values'''
        # could probably use `np.pad` but I'm not sure it would save code
        shape_old = np.array(self._fmnks.shape)
        shape_new = np.array((m+1, n+1, k+1))
        if np.any(shape_new > shape_old):
            shape = np.maximum(shape_old, shape_new)
            fmnks = -np.ones(shape)             # create the new array
            m0, n0, k0 = shape_old
            fmnks[:m0, :n0, :k0] = self._fmnks  # copy remembered values
            self._fmnks = fmnks

    def _f(self, m, n, k):
        '''Recursive implementation of function of [3] Theorem 2.5'''

        fmnks = self._fmnks  # for convenience

        # [3] Theorem 2.5 Line 1
        if k < 0 or m < 0 or n < 0 or k > m*n:
            return 0

        # if already calculated, return the value
        if fmnks[m, n, k] >= 0:
            return fmnks[m, n, k]

        if k == 0 and m >= 0 and n >= 0:  # [3] Theorem 2.5 Line 2
            fmnk = 1
        else:   # [3] Theorem 2.5 Line 3 / Equation 3
            fmnk = self._f(m-1, n, k-n)  +  self._f(m, n-1, k)

        fmnks[m, n, k] = fmnk  # remember result

        return fmnk


# Maintain state for faster repeat calls to mannwhitneyu2 w/ method='exact'
_mwu_state = _MWU()


def _get_mwu_z(U, n1, n2, ranks, axis=0, continuity=True):
    '''Standardized MWU statistic'''
    # Follows mannwhitneyu2 [2]
    m_u = n1 * n2 / 2
    n = n1 + n2

    # Tie correction. scipy.stats.tiecorrect is not vectorized; this is.
    # element i of t is the number of subjects sharing rank i
    _, t = np.unique(ranks, return_counts=True, axis=-1)
    s = np.sqrt(n1*n2/12 * ((n + 1) - (t**3 - t).sum(axis=-1)/(n*(n-1))))

    c = -0.5 * continuity
    z = (U - m_u + c) / s
    return z


def _mwu_input_validation(x, y, use_continuity, alternative, axis, method):
    ''' Input validation and standardization for mannwhitneyu2 '''
    # Would use np.asarray_chkfinite, but infs are OK
    x, y = np.asarray(x), np.asarray(y)
    if np.isnan(x).any() or np.isnan(y).any():
        raise ValueError('`x` and `y` must not contain NaNs.')
    if np.size(x) == 0 or np.size(y) == 0:
        raise ValueError('`x` and `y` must be of nonzero size.')

    bools = {True, False}
    if use_continuity not in bools:
        raise ValueError(f'`use_continuity` must be one of {bools}.')

    alternatives = {"two-sided", "less", "greater"}
    alternative = alternative.lower()
    if alternative not in alternatives:
        raise ValueError(f'`alternative` must be one of {alternatives}.')

    axis_int = int(axis)
    if axis != axis_int or axis_int < 0:
        raise ValueError('`index` must be a positive integer.')

    methods = {"asymptotic", "exact"}
    method = method.lower()
    if method not in methods:
        raise ValueError(f'`method` must be one of {methods}.')

    return x, y, use_continuity, alternative, axis_int, method

mwu_result = make_dataclass("MannWhitneyUResult", ("statistic", "pvalue"))
# Using `nametuple` for now to pass existing mannwhitneyu tests without
# modification
MannwhitneyuResult = namedtuple('MannwhitneyuResult', ('statistic', 'pvalue'))


def mannwhitneyu2(x, y, use_continuity=True, alternative="two-sided",
                  axis=0, method="asymptotic"):
    '''Perform the Mann-Whitney U rank test on two independent samples.

    The Mann-Whitney U test is a nonparametric test of the null hypothesis
    that the distribution underlying sample `x` is neither stochastically
    greater nor stochastically less than the distribution underlying sample
    `y`.

    Parameters
    ----------
    x, y : array-like
        N-d arrays of samples. The arrays must be broadcastable except along
        the dimension given by `axis`.
    use_continuity : bool, optional
            Whether a continuity correction (1/2) should be applied.
            Default is True when `method` is `asymptotic`; has no effect
            otherwise.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available:

          * 'two-sided' (default)
          * 'less': one-sided
          * 'greater': one-sided

    axis : int, optional
        Axis along which to perform the test. Default is 0.
    method : {'asymptotic', 'exact'}, optional
        Selects the method used to calculate the *p*-value. The following
        options are available.
        * 'asymptotic': (default) compares the standardized test statistic
          against the normal distribution, correcting for ties.
        * 'exact': computes the exact *p*-value by comparing the observed
          :math:`L` statistic against those realized by all possible
          permutations of ranks (under the null hypothesis that each
          permutation is equally likely)

    Returns
    -------
    res : MannwhitneyuResult
        An object containing attributes:

        statistic : float
            The Mann-Whitney U statistic corresponding with sample `x`.
        pvalue : float
            The associated *p*-value for the chosen `alternative`.

    Notes
    -----
    `method` 'exact' is recommended when there are no ties and either sample
    size is less than 8 _[1]. The implementation follows the original
    recurrence relation proposed in _[1] as described in _[3].

    If ``U1`` is the statistic corresponding with sample `x`, then the
    statistic corresponding with sample `y` is
    ``x.shape[axis] * y.shape[axis] - U1``.

    `mannwhitneyu2` is for independent samples. For related samples,
    consider `wilcoxon`.

    See Also
    --------
    `wilcoxon`

    References
    ----------
    .. [1] H.B. Mann and D.R. Whitney, "On a test of whether one of two random
           variables is stochastically larger than the other", The Annals of
           Mathematical Statistics, Vol. 18, pp. 50-60, 1947.
    .. [2] Mann-Whitney U Test, Wikipedia,
           http://en.wikipedia.org/wiki/Mann-Whitney_U_test
    .. [3] A. Di Bucchianico, "Combinatorics, computer algebra, and the
           Wilcoxon-Mann-Whitney test", Journal of Statistical Planning and
           Inference, Vol. 79, pp. 349-364, 1999.
    '''

    x, y, use_continuity, alternative, axis_int, method = (
        _mwu_input_validation(x, y, use_continuity, alternative, axis, method))

    x, y, xy = _broadcast_concatenate(x, y, axis)

    # Follows [2]
    n1, n2 = x.shape[-1], y.shape[-1]
    ranks = stats.rankdata(xy, axis=-1)  # method 2, step 1
    R1 = ranks[..., :n1].sum(axis=-1)    # method 2, step 2
    U1 = R1 - n1*(n1+1)/2                # method 2, step 3
    U2 = n1 * n2 - U1                    # as U1 + U2 = n1 * n2

    if alternative == "greater":
        U, f = U1, 1  # U is the statistic to use for p-value, f is a factor
    elif alternative == "less":
        U, f = U2, 1  # Due to symmetry, use SF of U2 rather than CDF of U1
    else:
        U, f = np.maximum(U1, U2), 2  # multiply SF by two for two-sided test

    if method == "exact":
        p = _mwu_state.sf(U.astype(int), n1, n2)
    elif method == "asymptotic":
        z = _get_mwu_z(U, n1, n2, ranks, continuity=use_continuity)
        p = stats.norm.sf(z)
    p *= f

    # return mwu_result(U, p)
    return MannwhitneyuResult(U1, p)  # temporary to integrate with tests
