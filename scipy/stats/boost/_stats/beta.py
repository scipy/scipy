from scipy.stats import rv_continuous

from scipy.stats.boost.beta_ufunc import (
    _beta_pdf, _beta_cdf, _beta_icdf, _beta_quantile, _beta_iquantile,
    _beta_mean, _beta_variance, _beta_skewness, _beta_kurtosis_excess,
)

class beta_gen(rv_continuous):
    def _rvs(self, a, b, size=None, random_state=None):
        return random_state.beta(a, b, size)

    def _pdf(self, x, a, b):
        return _beta_pdf(x, a, b)

    def _cdf(self, x, a, b):
        return _beta_cdf(x, a, b)

    def _sf(self, x, a, b):
        return _beta_icdf(x, a, b)

    def _isf(self, x, a, b):
        return _beta_iquantile(x, a, b)

    def _ppf(self, q, a, b):
        return _beta_quantile(q, a, b)

    def _stats(self, a, b):
        return(
            _beta_mean(a, b),
            _beta_variance(a, b),
            _beta_skewness(a, b),
            _beta_kurtosis_excess(a, b))


beta = beta_gen(a=0.0, b=1.0, name='beta')
