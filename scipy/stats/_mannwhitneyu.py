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

    def __init__(self):
        '''Lightweight initializer to remember previously calculated f/F'''
        self.fmnks = -np.ones((1, 1, 1))

    def pmf(self, k, m, n):
        '''Probability mass function'''
        self._resize_fmnks(m, n, k)
        return self._f(m, n, k) / special.binom(m + n, m)

    def cdf(self, k, m, n):
        '''Cumulative distribution function'''
        self._resize_fmnks(m, n, k)
        return (np.sum([self._f(m, n, i) for i in range(0, k + 1)])
                / special.binom(m + n, m))

    def sf(self, k, m, n):
        '''Survival function'''
        self._resize_fmnks(m, n, m*n)
        return (np.sum([self._f(m, n, i) for i in range(m*n, k, -1)])
                / special.binom(m + n, m))

    def _resize_fmnks(self, m, n, k):
        '''If necessary, expand the array that remembers PMF values'''
        # could probably use `np.pad` but I'm not sure it would save code
        shape_old = np.array(self.fmnks.shape)
        shape_new = np.array((m+1, n+1, k+1))
        if np.any(shape_new > shape_old):
            shape = np.maximum(shape_old, shape_new)
            fmnks = -np.ones(shape)            # create the new array
            m0, n0, k0 = shape_old
            fmnks[:m0, :n0, :k0] = self.fmnks  # copy remembered values
            self.fmnks = fmnks

    def _f(self, m, n, k):
        '''Recursive implementation of function of [3] Theorem 2.5'''

        fmnks = self.fmnks  # for convenience

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
    n = n1 + n2
    c = -0.5 * continuity
    m_u = n1 * n2 / 2

    # Tie correction. scipy.stats.tiecorrect is not vectorized; this is.
    _, t = np.unique(ranks, return_counts=True, axis=-1)
    s = np.sqrt(n1*n2/12 * (n1 + n2 + 1 - (t**3 - t).sum(axis=-1)/(n*(n-1))))

    z = (U + c - m_u) / s
    return z


mwu_result = make_dataclass("MannWhitneyUResult", ("statistic", "pvalue"))
# Using `nametuple` for now to pass existing mannwhitneyu tests without
# modificatino
MannwhitneyuResult = namedtuple('MannwhitneyuResult', ('statistic', 'pvalue'))


def mannwhitneyu2(x, y, continuity=True, alternative=None, axis=0,
                  exact=False):
    '''Mann-Whitney U Test

    Yet another implementation. Currently allows n-d x and y for asymptotic
    test only, but exact test can be vectorized, too. Calculation of
    exact distribution is memoized.

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
    x, y = np.asarray(x), np.asarray(y)
    x, y, xy = _broadcast_concatenate(x, y, axis)

    # Follows [2]
    n1, n2 = x.shape[-1], y.shape[-1]
    ranks = stats.rankdata(xy, axis=-1)
    R1 = ranks[..., :n1].sum(axis=-1)
    R2 = ranks[..., n1:].sum(axis=-1)
    U1 = R1 - n1*(n1+1)/2
    U2 = R2 - n2*(n2+1)/2

    if alternative == "greater":
        U, f = U1, 1  # U is the statistic to use for p-value, f is a factor
    elif alternative == "less":
        U, f = U2, 1
    else:
        U, f = np.maximum(U1, U2), 2  # multiply by two for two-sided test

    if exact:
        p = (_mwu_state.sf(int(U), n1, n2)  # because distribution is discrete,
             + _mwu_state.pmf(int(U), n1, n2))  # we need to add the PMF at U
    else:
        z = _get_mwu_z(U, n1, n2, ranks, continuity=continuity)
        p = stats.norm.sf(z)
    p *= f

    # return mwu_result(U, p)
    return MannwhitneyuResult(U1, p)  # temporary to integrate with tests