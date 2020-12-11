'''Boost drop-in replacement class for scipy.stats.nbinom.'''

import numpy as np
from scipy.stats import rv_discrete

from scipy.stats.boost.nbinom_ufunc import (
    _nbinom_pdf, _nbinom_cdf, _nbinom_icdf, _nbinom_quantile, _nbinom_iquantile,
    _nbinom_mean, _nbinom_variance, _nbinom_skewness, _nbinom_kurtosis_excess,
)


class nbinom_gen(rv_discrete):
    def _rvs(self, n, p, size=None, random_state=None):
        return random_state.negative_binomial(n, p, size)

    def _argcheck(self, n, p):
        return (n > 0) & (p >= 0) & (p <= 1)

    def _pmf(self, x, n, p):
        return _nbinom_pdf(x, n, p)

    def _cdf(self, x, n, p):
        return _nbinom_cdf(x, n, p)

    def _sf(self, x, n, p):
        return _nbinom_icdf(x, n, p)

    def _isf(self, x, n, p):
        return _nbinom_iquantile(x, n, p) + 1

    def _ppf(self, q, n, p):
        return _nbinom_quantile(q, n, p)

    def _stats(self, n, p):
        return(
            _nbinom_mean(n, p),
            _nbinom_variance(n, p),
            _nbinom_skewness(n, p),
            _nbinom_kurtosis_excess(n, p),
        )

nbinom = nbinom_gen(name='nbinom')
