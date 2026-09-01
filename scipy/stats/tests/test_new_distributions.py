# file for distribution-specific tests with new infrastructure (UnivariateDistribution)
import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy import stats

class TestBinomial:
    @pytest.mark.parametrize('fun', ['cdf', 'logcdf', 'ccdf', 'logccdf'])
    @pytest.mark.parametrize('method', ['quadrature', 'log/exp',
                                        'formula', 'complement'])
    def test_gh26072_non_integer_cdf_and_ccdf(self, fun, method):
        # gh-26072 found that cdf-like methods of discrete distributions
        # did not produce the expected step behavior
        n, p = 10, 0.3
        x = np.arange(n+1)
        x = np.concat((x, np.nextafter(x, np.inf), np.nextafter(x, -np.inf)))
        X = stats.Binomial(n=n, p=p)
        Y = stats.binom(n=n, p=p)
        X_fun = getattr(X, fun)
        Y_fun = getattr(Y, fun.replace('ccdf', 'sf'))
        assert_allclose(X_fun(x, method=method), Y_fun(x))
        assert_allclose(X_fun(x, method=method), X_fun(np.floor(x)))

    def test_gh23708_binomial_logcdf_method_complement(self):
        # gh-23708 found that `logcdf` method='complement' was inaccurate in the tails
        x = np.asarray([0., 18.])
        X = stats.Binomial(n=np.asarray([18.]), p=np.asarray(0.71022842))
        assert_allclose(X.logcdf(x, method='complement'), X.logcdf(x), rtol=1e-15)
        assert_allclose(X.logccdf(x, method='complement'), X.logccdf(x), rtol=1e-15)

        # going even deeper into the tails
        X = stats.Binomial(n=100, p=0.5)
        assert_allclose(X.logcdf(0, method='complement'), X.logpmf(0), rtol=1e-15)
        assert_allclose(X.logccdf(99, method='complement'), X.logpmf(100), rtol=1e-15)
