import warnings
import numpy as np

from scipy import stats
from ._stats_py import _normtest_finish, _rankdata
from . import _morestats
from ._axis_nan_policy import _broadcast_arrays
from ._hypotests import _get_wilcoxon_distr
from scipy._lib._util import _lazywhere


class WilcoxonDistribution:

    def __init__(self):
        self._dists = {}

    def _get_dist(self, n):
        if n not in self._dists:
            self._dists[n] = _get_wilcoxon_distr(n)
        return self._dists[n]

    def _cdf1(self, k, n):
        pmfs = self._get_dist(n)
        return pmfs[:k + 1].sum()

    def _cdf(self, k, n):
        return np.vectorize(self._cdf1, otypes=[float])(k, n)

    def _sf1(self, k, n):
        pmfs = self._get_dist(n)
        return pmfs[k:].sum()

    def _sf(self, k, n):
        return np.vectorize(self._sf1, otypes=[float])(k, n)

    def mean(self, n):
        return n * (n + 1) // 4

    def _prep(self, k, n):
        k = np.asarray(k).astype(int, copy=False)
        n = np.asarray(n).astype(int, copy=False)
        mn = self.mean(n)
        out = np.empty(k.shape, dtype=np.float64)
        return k, n, mn, out

    def cdf(self, k, n):
        k, n, mn, out = self._prep(k, n)
        return _lazywhere(k <= mn, (k, n), self._cdf,
                          f2=lambda k, n: 1 - self._sf(k+1, n))[()]

    def sf(self, k, n):
        k, n, mn, out = self._prep(k, n)
        return _lazywhere(k <= mn, (k, n), self._sf,
                          f2=lambda k, n: 1 - self._cdf(k-1, n))[()]


def _wilcoxon_iv(x, y, zero_method, correction, alternative, method, axis):

    axis = np.asarray(axis)[()]
    message = "`axis` must be an integer."
    if not np.issubdtype(axis.dtype, np.integer) or axis.ndim != 0:
        raise ValueError(message)

    message = '`axis` must be compatible with the shape(s) of `x` (and `y`)'
    try:
        if y is None:
            x = np.asarray(x)
            d = x
        else:
            x, y = _broadcast_arrays((x, y), axis=axis)
            d = x - y
        d = np.moveaxis(d, axis, -1)
    except np.AxisError as e:
        raise ValueError(message) from e

    message = "`x` and `y` must have the same length along `axis`."
    if y is not None and x.shape[axis] != y.shape[axis]:
        raise ValueError(message)

    message = "`x` (and `y`, if provided) must be an array of real numbers."
    if np.issubdtype(d.dtype, np.integer):
        d = d.astype(np.float64)
    if not np.issubdtype(d.dtype, np.floating):
        raise ValueError(message)

    zero_method = str(zero_method).lower()
    zero_methods = {"wilcox", "pratt", "zsplit"}
    message = f"`zero_method` must be one of {zero_methods}."
    if zero_method not in zero_methods:
        raise ValueError(message)

    corrections = {True, False}
    message = f"`correction` must be one of {corrections}."
    if correction not in corrections:
        raise ValueError(message)

    alternative = str(alternative).lower()
    alternatives = {"two-sided", "less", "greater"}
    message = f"`alternative` must be one of {alternatives}."
    if alternative not in alternatives:
        raise ValueError(message)

    if not isinstance(method, stats.PermutationMethod):
        methods = {"auto", "approx", "exact"}
        message = (f"`method` must be one of {methods} or "
                   "an instance of `stats.PermutationMethod`.")
        if method not in methods:
            raise ValueError(message)

    message = ("Zeros are present, but `method='exact'` does not compute "
               "accurate p-values in the presence of zeros.")
    n_zero = np.sum(d == 0, axis=-1)
    has_zeros = np.any(n_zero > 0)
    if has_zeros and method == "exact":
        warnings.warn(message, stacklevel=2)

    message = ("The sample size is small (fewer than 10 nonzero elements); "
               "so `method='approx'` may not produce accurate p-values.")
    count = d.shape[-1] - n_zero
    too_small = np.any(count < 10)
    if too_small and method == "approx":
        warnings.warn(message, stacklevel=2)

    if method == "auto":
        if d.shape[-1] <= 50 and not has_zeros:
            method = "exact"
        elif not too_small:
            method = "approx"
        else:
            method = stats.PermutationMethod()

    return d, zero_method, correction, alternative, method, axis


def _wilcoxon_statistic(d, zero_method='wilcox'):

    i_zeros = (d == 0)

    if zero_method == 'wilcox':
        # Wilcoxon's method for treating zeros was to remove them from
        # the calculation. We do this by replacing 0s with NaNs, which
        # are ignored anyway.
        d[i_zeros] = np.nan

    i_nan = np.isnan(d)
    n_nan = np.sum(i_nan, axis=-1)
    count = d.shape[-1] - n_nan

    r, t = _rankdata(abs(d), 'average', return_ties=True)
    r_plus = np.sum((d > 0) * r, axis=-1)
    r_minus = np.sum((d < 0) * r, axis=-1)

    if zero_method == "zsplit":
        # The "zero-split" method for treating zeros is to add half their contribution
        # to r_plus and half to r_minus.
        # See gh-2263 for the origin of this method.
        r_zero_2 = np.sum(i_zeros * r, axis=-1) / 2
        r_plus += r_zero_2
        r_minus += r_zero_2

    mn = count * (count + 1.) * 0.25
    se = count * (count + 1.) * (2. * count + 1.)

    if zero_method == "pratt":
        # Pratt's method for treating zeros was just to modify the z-statistic.

        # normal approximation needs to be adjusted, see Cureton (1967)
        n_zero = i_zeros.sum(axis=-1)
        mn -= n_zero * (n_zero + 1.) * 0.25
        se -= n_zero * (n_zero + 1.) * (2. * n_zero + 1.)

        # zeros are not to be included in tie-correction.
        # any tie counts corresponding with zeros are in the 0th column
        t[i_zeros.any(axis=-1), 0] = 0

    tie_correct = (t**3 - t).sum(axis=-1)
    se -= tie_correct/2
    se = np.sqrt(se / 24)

    # need to add continuity correction
    z = (r_plus - mn) / se

    return r_plus, r_minus, se, z, count


def _correction_sign(z, alternative):
    if alternative == 'greater':
        return 1
    elif alternative == 'less':
        return -1
    else:
        return np.sign(z)


def _wilcoxon_nd(x, y=None, zero_method='wilcox', correction=True, alternative='two-sided',
                 method='auto', axis=0):
    # test wilcox/pratt with all zeros
    # change: respects specified method
    temp = _wilcoxon_iv(x, y, zero_method, correction, alternative, method, axis)
    d, zero_method, correction, alternative, method, axis = temp

    r_plus, r_minus, se, z, count = _wilcoxon_statistic(d, zero_method)

    if method == 'approx':
        if correction:
            sign = _correction_sign(z, alternative)
            z -= sign * 0.5 / se
        z, p = _normtest_finish(z, alternative)
    elif method == 'exact':
        dist = WilcoxonDistribution()
        if alternative == 'less':
            p = dist.cdf(r_plus, count)
        elif alternative == 'greater':
            p = dist.sf(r_plus, count)
        else:
            p = 2 * np.minimum(dist.sf(r_plus, count),
                               dist.cdf(r_plus, count))
            p = np.clip(p, 0, 1)
    else:
        p = stats.permutation_test(
            (d,), lambda d: _wilcoxon_statistic(d, zero_method)[0],
            permutation_type='samples', alternative=alternative, axis=-1).pvalue

    res = _morestats.WilcoxonResult(statistic=r_plus[()], pvalue=p[()])
    if method == 'approx':
        res.zstatistic = z[()]
    return res
