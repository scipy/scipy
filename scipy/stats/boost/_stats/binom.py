'''Boost drop-in replacement for scipy.stats.binom.'''

import numpy as np
from scipy.stats import rv_discrete

from scipy.stats.boost.binom_ufunc import (
    _binom_pdf, _binom_cdf, _binom_icdf, _binom_quantile, _binom_iquantile,
    _binom_mean, _binom_variance, _binom_skewness, _binom_kurtosis_excess,
)


class binom_gen(rv_discrete):
    def _rvs(self, n, p, size=None, random_state=None):
        return random_state.binomial(n, p, size)

    def _argcheck(self, n, p):
        return (n >= 0) & (p >= 0) & (p <= 1)

    def _get_support(self, n, p):
        return self.a, n

    def _pmf(self, x, n, p):
        return _binom_pdf(x, n, p)

    def _cdf(self, x, n, p):
        return _binom_cdf(x, n, p)

    def _sf(self, x, n, p):
        return _binom_icdf(x, n, p)

    def _isf(self, x, n, p):
        return _binom_iquantile(x, n, p) + 1

    def _ppf(self, q, n, p):
        return _binom_quantile(q, n, p)

    def _stats(self, n, p):
        return(
            _binom_mean(n, p),
            _binom_variance(n, p),
            _binom_skewness(n, p),
            _binom_kurtosis_excess(n, p),
        )

binom = binom_gen(name='binom')
