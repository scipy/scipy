# file for distribution-specific tests with new infrastructure (UnivariateDistribution)
import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy import stats
from scipy.stats.tests.test_continuous import DistributionsTest
from scipy.stats._new_distributions import StandardNormal


class TestBinomial(DistributionsTest):
    seed = 706381675
    family = stats.Binomial

    def is_degenerate(self, dist):
        return np.any((dist.p == 0) | (dist.p == 1) | np.isnan(dist.p))

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case):
        if self.is_degenerate(case.dist):
            with np.errstate(invalid='ignore'):
                return super().test_moment(case)
        super().test_moment(case)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_skewness(self, case):
        if self.is_degenerate(case.dist):
            with np.errstate(invalid='ignore'):
                return super().test_skewness(case)
        super().test_moment(case)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_kurtosis(self, case):
        if self.is_degenerate(case.dist):
            with np.errstate(invalid='ignore'):
                return super().test_kurtosis(case)
        super().test_moment(case)

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


class TestLogistic(DistributionsTest):
    seed = 389513556
    family = stats.Logistic

    @pytest.mark.filterwarnings("ignore:divide:RuntimeWarning")
    def test_cdf2(self, case):
        return super().test_cdf2(case)


class TestNormal(DistributionsTest):
    seed = 353965734
    family = stats.Normal

    @pytest.mark.filterwarnings("ignore:divide:RuntimeWarning")
    def test_logpdf(self, case):
        return super().test_logpdf(case)

    def test_lmoment(self, case):
        return super().test_lmoment(case, tol_override={'atol': 1e-8})


class TestStandardNormal(DistributionsTest):
    seed = 726527242
    family = StandardNormal

    @pytest.mark.filterwarnings("ignore:divide:RuntimeWarning")
    def test_cdf2(self, case):
        return super().test_cdf2(case)


class TestUniform(DistributionsTest):
    seed = 893709074
    family = stats.Uniform

    def test_mode(self, case):
        assert_allclose(case.dist.mode(), case.dist.a + case.dist.ab/2)

    @pytest.mark.thread_unsafe(reason="looks like an _rng_spawn issue?")
    @pytest.mark.fail_slow(10)
    def test_quasi_random_sample(self, case):
        return super().test_quasi_random_sample(case)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case):
        return super().test_moment(case, tol_override={'atol': 1e-9})
