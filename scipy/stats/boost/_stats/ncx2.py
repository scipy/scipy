'''Boost drop-in replacement for  scipy.stats.ncx2.'''

import numpy as np
from scipy.stats import rv_continuous, chi2
from scipy._lib._util import _lazywhere

from scipy.stats.boost.ncx2_ufunc import (
    _ncx2_pdf, _ncx2_cdf, _ncx2_icdf, _ncx2_quantile, _ncx2_iquantile,
    _ncx2_mean, _ncx2_variance, _ncx2_skewness, _ncx2_kurtosis_excess,
)


class ncx2_gen(rv_continuous):
    '''Only requires substitution of _distn_infrastructure._ncx2_* -> boost._ncx2_*.'''
    def _argcheck(self, df, nc):
        return (df > 0) & (nc >= 0)

    def _rvs(self, df, nc, size=None, random_state=None):
        return random_state.noncentral_chisquare(df, nc, size)

    def _pdf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (x, df, nc), f=_ncx2_pdf, f2=chi2.pdf)
        # return _ncx2_pdf(x, df, nc)

    def _cdf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (x, df, nc), f=_ncx2_cdf, f2=chi2.cdf)
        # return _ncx2_cdf(x, df, nc)

    def _sf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (x, df, nc), f=_ncx2_icdf, f2=chi2.sf)
        # return _ncx2_icdf(x, df, nc)

    def _isf(self, x, df, nc):
        cond = np.ones_like(x, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (x, df, nc), f=_ncx2_iquantile, f2=chi2.isf)
        # return _ncx2_iquantile(x, df, nc)

    def _ppf(self, q, df, nc):
        cond = np.ones_like(q, dtype=bool) & (nc != 0)
        return _lazywhere(cond, (q, df, nc), f=_ncx2_quantile, f2=chi2.ppf)
        # return _ncx2_quantile(q, df, nc)

    def _stats(self, df, lamda):
        return(
            _ncx2_mean(df, lamda),
            _ncx2_variance(df, lamda),
            _ncx2_skewness(df, lamda),
            _ncx2_kurtosis_excess(df, lamda))

ncx2 = ncx2_gen(a=0.0, name='ncx2')
