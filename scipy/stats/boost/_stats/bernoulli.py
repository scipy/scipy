'''Boost drop-in replacement for scipy.stats.bernoulli.'''

import numpy as np
from scipy.special import entr

from scipy.stats.boost._stats.binom import binom_gen

from scipy.stats.boost.bernoulli_ufunc import (
    _bernoulli_pdf, _bernoulli_cdf, _bernoulli_icdf, _bernoulli_quantile, _bernoulli_iquantile,
    _bernoulli_mean, _bernoulli_variance, _bernoulli_skewness, _bernoulli_kurtosis_excess,
)


class bernoulli_gen(binom_gen):
    def _rvs(self, p, size=None, random_state=None):
        return binom_gen._rvs(self, 1, p, size=size, random_state=random_state)

    def _argcheck(self, p):
        return (p >= 0) & (p <= 1)

    def _get_support(self, p):
        return self.a, self.b

    def _pmf(self, x, p):
        return _bernoulli_pdf(x, p)

    def _cdf(self, x, p):
        return _bernoulli_cdf(x, p)

    def _isf(self, x, p):
        return _bernoulli_iquantile(x, p)

    def _sf(self, x, p):
        return _bernoulli_icdf(x, p)

    def _ppf(self, q, p):
        return _bernoulli_quantile(q, p)

    def _stats(self, p):
        return(
            _bernoulli_mean(p),
            _bernoulli_variance(p),
            _bernoulli_skewness(p),
            _bernoulli_kurtosis_excess(p),
        )

    def median(self, p):
        return _bernoulli_quantile(0.5, p)

    def _entropy(self, p):
        return entr(p) + entr(1-p)


bernoulli = bernoulli_gen(b=1, name='bernoulli')
