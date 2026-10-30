# file for distribution-specific tests with new infrastructure (UnivariateDistribution)
import pytest
import numpy as np
from scipy._external import array_api_extra as xpx
from scipy._lib._array_api import make_xp_test_case, is_cupy
from scipy._lib._array_api_no_0d import xp_assert_close
from scipy import stats, special
from scipy.stats.tests.test_continuous import DistributionsTest
from scipy.stats._new_distributions import StandardNormal


@make_xp_test_case(stats.Binomial)
class TestBinomial(DistributionsTest):
    seed = 706381677
    family = stats.Binomial

    def is_degenerate(self, dist, xp):
        return xp.any((dist.p == 0) | (dist.p == 1) | xp.isnan(dist.p))

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        xp_assert_close(dist.pmf(x),
                        stats.binom(p=dist.p, n=dist.n).pmf(x))

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case, *, xp):
        if self.is_degenerate(case.dist, xp=xp):
            with np.errstate(invalid='ignore'):
                return super().test_moment(case, xp=xp)
        super().test_moment(case, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_skewness(self, case, *, xp):
        if self.is_degenerate(case.dist, xp=xp):
            with np.errstate(invalid='ignore'):
                return super().test_skewness(case, xp=xp)
        super().test_skewness(case, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_kurtosis(self, case, *, xp):
        if self.is_degenerate(case.dist, xp=xp):
            with np.errstate(invalid='ignore'):
                return super().test_kurtosis(case, xp=xp)
        super().test_kurtosis(case, xp=xp)

    @pytest.mark.parametrize('fun', ['cdf', 'logcdf', 'ccdf', 'logccdf'])
    @pytest.mark.parametrize('method', ['quadrature', 'log/exp',
                                        'formula', 'complement'])
    def test_gh26072_non_integer_cdf_and_ccdf(self, fun, method, xp):
        # gh-26072 found that cdf-like methods of discrete distributions
        # did not produce the expected step behavior
        n, p = 10, 0.3
        x = np.arange(n+1)
        x = np.concat((x, np.nextafter(x, np.inf), np.nextafter(x, -np.inf)))
        x_xp = xp.asarray(x)
        X = stats.Binomial(n=n, p=p, xp=xp)
        Y = stats.binom(n=n, p=p)
        X_fun = getattr(X, fun)
        Y_fun = getattr(Y, fun.replace('ccdf', 'sf'))
        xp_assert_close(X_fun(x_xp, method=method), xp.asarray(Y_fun(x)))
        xp_assert_close(X_fun(x_xp, method=method), X_fun(np.floor(x)))

    def test_gh23708_binomial_logcdf_method_complement(self, xp):
        # gh-23708 found that `logcdf` method='complement' was inaccurate in the tails
        x = xp.asarray([0., 18.])
        X = stats.Binomial(n=xp.asarray([18.]), p=xp.asarray(0.71022842))
        xp_assert_close(X.logcdf(x, method='complement'), X.logcdf(x), rtol=1e-15)
        xp_assert_close(X.logccdf(x, method='complement'), X.logccdf(x), rtol=1e-15)

        # going even deeper into the tails
        X = stats.Binomial(n=100, p=0.5)
        xp_assert_close(X.logcdf(0, method='complement'), X.logpmf(0), rtol=1e-15)
        xp_assert_close(X.logccdf(99, method='complement'), X.logpmf(100), rtol=1e-15)


@make_xp_test_case(stats.Logistic)
class TestLogistic(DistributionsTest):
    seed = 389513556
    family = stats.Logistic

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        xp_assert_close(dist.pdf(x),
                        xp.exp(-x) / (1+xp.exp(-x))**2)

    @pytest.mark.filterwarnings("ignore:divide:RuntimeWarning")
    def test_cdf2(self, case, xp):
        return super().test_cdf2(case, xp=xp)


@make_xp_test_case(stats.Normal)
class TestNormal(DistributionsTest):
    seed = 353965734
    family = stats.Normal

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        y = (x - dist.mu) / dist.sigma
        xp_assert_close(dist.pdf(x), xp.exp(-y**2/2)/xp.sqrt(dist.sigma**2*2*xp.pi))

    @pytest.mark.filterwarnings("ignore:divide:RuntimeWarning")
    def test_logpdf(self, case, xp):
        return super().test_logpdf(case, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case, xp):
        atol = 1e-8 if case.dist._dtype == xp.float64 else 5e-3
        return super().test_moment(case, tol_override={'atol': atol}, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_lmoment(self, case, xp):
        atol = 1e-8 if case.dist._dtype == xp.float64 else 5e-3
        return super().test_lmoment(case, tol_override={'atol': atol}, xp=xp)


@make_xp_test_case(StandardNormal)
class TestStandardNormal(DistributionsTest):
    seed = 726527245
    family = StandardNormal

    def tol(self, case, xp):
        return dict(atol=1e-6) if case.dist._dtype == xp.float32 else dict()

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        xp_assert_close(dist.pdf(x), xp.exp(-x**2/2)/(2*xp.pi)**0.5)

    def test_cdf(self, case, xp):
        return super().test_cdf(case, xp=xp, tol_override=self.tol(case, xp))

    def test_cdf2(self, case, xp):
        return super().test_cdf2(case, xp=xp, tol_override=self.tol(case, xp))

    def test_ccdf(self, case, xp):
        return super().test_ccdf(case, xp=xp, tol_override=self.tol(case, xp))

    def test_ccdf2(self, case, xp):
        return super().test_ccdf2(case, xp=xp, tol_override=self.tol(case, xp))

    @pytest.mark.filterwarnings("ignore:divide:RuntimeWarning")
    def test_logpdf(self, case, xp):
        return super().test_logpdf(case, xp=xp)

    def test_logcdf(self, case, xp):
        return super().test_logcdf(case, xp=xp, tol_override=self.tol(case, xp))

    def test_logccdf(self, case, xp):
        return super().test_logccdf(case, xp=xp, tol_override=self.tol(case, xp))

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case, xp):
        tol_override = dict(atol=5e-3) if case.dist._dtype == xp.float32 else None
        return super().test_moment(case, tol_override=tol_override, xp=xp)


@make_xp_test_case(stats.Uniform)
class TestUniform(DistributionsTest):
    seed = 893709074
    family = stats.Uniform

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        _, c = xp.broadcast_arrays(x, dist.b - dist.a)
        xp_assert_close(dist.pdf(x), 1 / c)

    def test_mode(self, case, xp):
        xp_assert_close(case.dist.mode(), case.dist.a + case.dist.ab/2)

    @pytest.mark.thread_unsafe(reason="looks like an _rng_spawn issue?")
    @pytest.mark.fail_slow(10)
    def test_quasi_random_sample(self, case, xp):
        return super().test_quasi_random_sample(case, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case, xp):
        if xpx.default_dtype(xp) == xp.float32:
            pytest.xfail("Significant accuracy issues; deserves investigation")
        return super().test_moment(case, tol_override={'atol': 1e-9}, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_lmoment(self, case, xp):
        if xpx.default_dtype(xp) == xp.float32:
            pytest.xfail("Significant accuracy issues; deserves investigation")
        return super().test_lmoment(case, xp=xp)


@make_xp_test_case(stats.VonMises)
class TestVonMises(DistributionsTest):
    seed = 6954568351
    family = stats.VonMises

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        ref = (xp.exp(dist.kappa * xp.cos(x - dist.mu))
               / (2 * xp.pi * special.i0(dist.kappa)))
        xp_assert_close(dist.pdf(x), ref)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case, xp):
        # Formula not available on CuPy. Correctly falls back to quadrature, but
        # test complains that formula is not implemented. Not sure the best way
        # to handle xp-dependent adjustments to tests.
        if is_cupy(xp):
            pytest.skip("moment formula not available for CuPy")
        return super().test_moment(case, xp=xp)

    def test_median(self, case, xp):
        # can only expect about half precision with optimization
        atol = 1e-6 if case.dist._dtype == xp.float64 else 1e-3
        return super().test_median(case, tol_override={'atol': atol}, xp=xp)

    def test_mode(self, case, xp):
        # can only expect about half precision with optimization
        atol = 1e-6 if case.dist._dtype == xp.float64 else 1e-3
        return super().test_mode(case, tol_override={'atol': atol}, xp=xp)
