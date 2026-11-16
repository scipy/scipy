import math
import os

import numpy as np
from numpy import inf
import pytest
from numpy.testing import assert_allclose, assert_equal

from scipy import special
from scipy import stats
from scipy.stats._distr_params import distcont, distdiscrete
from scipy.stats._univariate_distribution_util import (
    _RealInterval, _RealParameter)
from scipy.stats._univariate_distribution import ContinuousDistribution
from scipy.stats._univariate_distribution_extensions import Mixture
from scipy.stats._new_distributions import _Gamma, _StandardCauchy
from scipy.stats import Normal, Uniform
from scipy.stats.tests.test_univariate_distribution import draw_distribution_from_family
from scipy._lib._array_api import make_xp_test_case, is_numpy, is_cupy, xp_copy
from scipy._lib._array_api_no_0d import xp_assert_close


@make_xp_test_case(stats.make_distribution)
class TestMakeDistribution:
    @pytest.mark.parametrize('i, distdata', list(enumerate(distcont + distdiscrete)))
    def test_rv_generic(self, i, distdata):
        distname = distdata[0]

        slow = {'argus', 'exponpow', 'exponweib', 'genexpon', 'gompertz',
                'johnsonsb', 'kappa4', 'ksone', 'kstwo', 'kstwobign', 'norminvgauss',
                'powerlognorm', 'powernorm', 'recipinvgauss', 'studentized_range',
                'vonmises_line'}  # continuous
        if not int(os.environ.get('SCIPY_XSLOW', '0')) and distname in slow:
            pytest.skip('Skipping as XSLOW')

        if distname in {              # skip these distributions
            'levy_stable',            # private methods seem to require >= 1d args
            'vonmises',               # circular distribution; shouldn't work
            'poisson_binom',          # vector shape parameter
            'hypergeom',              # distribution functions need interpolation
            'nchypergeom_fisher',     # distribution functions need interpolation
            'nchypergeom_wallenius',  # distribution functions need interpolation
        }:
            return

        # skip single test, mostly due to slight disagreement
        custom_tolerances = {'ksone': 1e-5, 'kstwo': 1e-5}  # discontinuous PDF
        skip_entropy = {'kstwobign', 'pearson3'}  # tolerance issue
        skip_skewness = {'exponpow', 'ksone', 'nchypergeom_wallenius'}  # tolerance
        skip_kurtosis = {'chi', 'exponpow', 'invgamma',  # tolerance
                         'johnsonsb', 'ksone', 'kstwo',  # tolerance
                         'nchypergeom_wallenius'}  # tolerance
        skip_logccdf = {'arcsine', 'skewcauchy', 'trapezoid', 'triang',  # tolerance
                        'dlaplace'}  # slight tolerance issue, but only for array shape
        skip_raw = {2: {'alpha', 'foldcauchy', 'halfcauchy', 'levy', 'levy_l'},
                    3: {'pareto'},  # stats.pareto is just wrong
                    4: {'invgamma'}}  # tolerance issue
        skip_standardized = {'exponpow', 'ksone'}  # tolerances

        dist = getattr(stats, distname)
        params = dict(zip(dist.shapes.split(', '), distdata[1])) if dist.shapes else {}
        rng = np.random.default_rng(7548723590230982)
        CustomDistribution = stats.make_distribution(dist)
        params = {key: np.asarray([val, val]) for key, val in params.items()}
        X = CustomDistribution(**params)
        Y = dist(**params)
        x = X.sample(shape=10, rng=rng)
        p = X.cdf(x)
        rtol = custom_tolerances.get(distname, 1e-7)
        atol = 1e-12

        with np.errstate(divide='ignore', invalid='ignore'):
            m, v, s, k = Y.stats('mvsk')
            assert_allclose(X.support(), Y.support())
            if distname not in skip_entropy:
                assert_allclose(X.entropy(), Y.entropy(), rtol=rtol)
            if isinstance(Y, stats.rv_discrete):
                # some continuous distributions have trouble with `logentropy` because
                # it uses complex numbers
                assert_allclose(np.exp(X.logentropy()), Y.entropy(), rtol=rtol)
            assert_allclose(X.median(), Y.median(), rtol=rtol)
            assert_allclose(X.mean(), m, rtol=rtol, atol=atol)
            assert_allclose(X.variance(), v, rtol=rtol, atol=atol)
            if distname not in skip_skewness:
                assert_allclose(X.skewness(), s, rtol=rtol, atol=atol)
            if distname not in skip_kurtosis:
                assert_allclose(X.kurtosis(convention='excess'), k,
                                rtol=rtol, atol=atol)
            if isinstance(dist, stats.rv_continuous):
                assert_allclose(X.logpdf(x), Y.logpdf(x), rtol=rtol)
                assert_allclose(X.pdf(x), Y.pdf(x), rtol=rtol)
            else:
                assert_allclose(X.logpmf(x), Y.logpmf(x), rtol=rtol)
                assert_allclose(X.pmf(x), Y.pmf(x), rtol=rtol)
            assert_allclose(X.logcdf(x), Y.logcdf(x), rtol=rtol)
            assert_allclose(X.cdf(x), Y.cdf(x), rtol=rtol)
            if distname not in skip_logccdf:
                assert_allclose(X.logccdf(x), Y.logsf(x), rtol=rtol)
            assert_allclose(X.ccdf(x), Y.sf(x), rtol=rtol)

            # old infrastructure convention for ppf(p=0) and isf(p=1) is different than
            # new infrastructure. Adjust reference values accordingly.
            a, _ = Y.support()
            a, p = np.broadcast_arrays(a, p)
            ref_ppf = Y.ppf(p)
            ref_ppf[p == 0] = a[p == 0]
            ref_isf = Y.isf(p)
            ref_isf[p == 1] = a[p == 1]

            assert_allclose(X.icdf(p), ref_ppf, rtol=rtol)
            assert_allclose(X.iccdf(p), ref_isf, rtol=rtol)

            for order in range(5):
                if distname not in skip_raw.get(order, {}):
                    assert_allclose(X.moment(order, kind='raw'),
                                    Y.moment(order), rtol=rtol, atol=atol)
            for order in range(3, 4):
                if distname not in skip_standardized:
                    assert_allclose(X.moment(order, kind='standardized'),
                                    Y.stats('mvsk'[order-1]), rtol=rtol, atol=atol)
            if isinstance(dist, stats.rv_continuous):
                # For discrete distributions, these won't agree at the far left end
                # of the support, and the new infrastructure is slow there (for now).
                seed = 845298245687345
                assert_allclose(X.sample(shape=10, rng=seed),
                                Y.rvs(size=p.shape,
                                      random_state=np.random.default_rng(seed)),
                                rtol=rtol)

    def test_zipfian_gh23678(self):
        # Regression test for gh-23678: calling the cdf method at the end
        # point of the Zipfian distribution would generate a warning.
        Zipfian = stats.make_distribution(stats.zipfian)
        zdist = Zipfian(a=0.75, n=15)
        # This should not generate any warnings.
        assert_equal(zdist.cdf(15), 1.0)

    def test_custom(self, xp):
        rng = np.random.default_rng(7548723590230982)

        class MyLogUniform:
            @property
            def __make_distribution_version__(self):
                return "1.16.0"

            @property
            def parameters(self):
                return {'a': {'endpoints': (0, xp.inf), 'inclusive': (False, False)},
                        'b': {'endpoints': ('a', xp.inf), 'inclusive': (False, False)}}

            @property
            def support(self):
                return {'endpoints': ('a', 'b')}

            def pdf(self, x, a, b):
                return 1 / (x * (xp.log(b) - xp.log(a)))

            def sample(self, shape, *, a, b, rng=None):
                # hmm... self._xp does not exist. Makes sense, but how
                # can user-defined classes know which backend is being used?
                dtype = (np.float32 if 'float32' in str(xp.result_type(a, b))
                         else np.float64)
                p = xp.asarray(rng.random(size=shape, dtype=dtype))
                return xp.exp(xp.log(a) + p * (xp.log(b) - xp.log(a)))

            def moment(self, order, kind='raw', *, a, b):
                if order == 1 and kind == 'raw':
                    # quadrature is perfectly accurate here; add 1e-10 error so we
                    # can tell the difference between the two
                    return (b - a) / xp.log(b/a) + 1e-10

            def lmoment(self, order, *, a, b):
                s = xp.log(b/a)
                l1 = (b - a) / s
                l2 = ((s - 2) * b + (s + 2) * a) / s**2
                l3 = ((s**2 - 6*s + 12) * b - (s**2 + 6*s + 12) * a) / s ** 3
                l4 = ((s**3 - 12*s**2 + 60*s - 120) * b
                      + (s**3 + 12*s**2 + 60*s + 120) * a) / s ** 4
                ls = {1: l1, 2: l2, 3: l3, 4: l4}
                return ls.get(int(order), None)

        LogUniform = stats.make_distribution(MyLogUniform())

        X = LogUniform(a=1., b=xp.e, xp=xp)
        Y = stats.exp(Uniform(a=0., b=1., xp=xp))

        # assert str(X) == f"MyLogUniform(a=1.0, b={xp.e})"
        # assert repr(X) == f"MyLogUniform(a=np.float64(1.0), b=np.float64({np.e}))"

        x = xp.asarray(X.sample(shape=10, rng=rng))
        p = X.cdf(x)

        xp_assert_close(X.support()[0], Y.support()[0])
        xp_assert_close(X.support()[1], Y.support()[1])
        xp_assert_close(X.entropy(), Y.entropy())
        xp_assert_close(X.median(), Y.median())
        xp_assert_close(X.logpdf(x), Y.logpdf(x))
        xp_assert_close(X.pdf(x), Y.pdf(x))
        xp_assert_close(X.logcdf(x), Y.logcdf(x))
        xp_assert_close(X.cdf(x), Y.cdf(x))
        xp_assert_close(X.logccdf(x), Y.logccdf(x))
        xp_assert_close(X.ccdf(x), Y.ccdf(x))
        xp_assert_close(X.icdf(p), Y.icdf(p))
        xp_assert_close(X.iccdf(p), Y.iccdf(p))
        for kind in ['raw', 'central', 'standardized']:
            for order in range(5):
                xp_assert_close(X.moment(order, kind=kind),
                                Y.moment(order, kind=kind))
        for standardize in [False, True]:
            for order in range(1, 5):
                if standardize and order < 3:
                    continue
                xp_assert_close(X.lmoment(order, standardize=standardize),
                                Y.lmoment(order, standardize=standardize),
                                atol=5e-13 if X._dtype == xp.float64 else 5e-3)

        # Confirm that the `sample` and `moment` methods are overridden as expected
        sample_formula = X.sample(shape=10, rng=0, method='formula')
        sample_inverse = X.sample(shape=10, rng=0, method='inverse_transform')
        xp_assert_close(sample_formula, sample_inverse)
        assert not xp.all(sample_formula == sample_inverse)

        xp_assert_close(X.mean(method='formula'), X.mean(method='quadrature'))
        assert not X.mean(method='formula') == X.mean(method='quadrature')

        method = 'quadrature_icdf' if is_numpy(xp) else 'order_statistics'
        xp_assert_close(X.lmoment(method='formula'),
                        X.lmoment(method=method))
        assert not X.lmoment(method='formula') == X.lmoment(method=method)

    # pdf and cdf formulas below can warn on boundary of support in some cases.
    # See https://github.com/scipy/scipy/pull/22560#discussion_r1962763840.
    @pytest.mark.slow
    @pytest.mark.filterwarnings("ignore::RuntimeWarning")
    @pytest.mark.parametrize("c", [-1, 0, 1, np.asarray([-2.1, -1., 0., 1., 2.1])])
    def test_custom_variable_support(self, c, xp):
        rng = np.random.default_rng(7548723590230982)

        class MyGenExtreme:
            @property
            def __make_distribution_version__(self):
                return "1.16.0"

            @property
            def parameters(self):
                return {
                    'c': {'endpoints': (-xp.inf, xp.inf), 'inclusive': (False, False)},
                    'mu': {'endpoints': (-xp.inf, xp.inf), 'inclusive': (False, False)},
                    'sigma': {'endpoints': (0, xp.inf), 'inclusive': (False, False)}
                }

            @property
            def support(self):
                def left(*, c, mu, sigma):
                    c, mu, sigma = xp.broadcast_arrays(c, mu, sigma)
                    result = xp.empty_like(c)
                    result[c >= 0] = -xp.inf
                    result[c < 0] = mu[c < 0] + sigma[c < 0] / c[c < 0]
                    return result[()]

                def right(*, c, mu, sigma):
                    c, mu, sigma = xp.broadcast_arrays(c, mu, sigma)
                    result = xp.empty_like(c)
                    result[c <= 0] = xp.inf
                    result[c > 0] = mu[c > 0] + sigma[c > 0] / c[c > 0]
                    return result[()]

                return {"endpoints": (left, right), "inclusive": (False, False)}

            def pdf(self, x, *, c, mu, sigma):
                x, c, mu, sigma = xp.broadcast_arrays(x, c, mu, sigma)
                t = xp.empty_like(x)
                mask = (c == 0)
                t[mask] = xp.exp(-(x[mask] - mu[mask])/sigma[mask])
                t[~mask] = (
                    1  - c[~mask]*(x[~mask] - mu[~mask])/sigma[~mask]
                )**(1/c[~mask])
                result = 1/sigma * t**(1 - c)*xp.exp(-t)
                return result[()]

            def cdf(self, x, *, c, mu, sigma):
                x, c, mu, sigma = xp.broadcast_arrays(x, c, mu, sigma)
                t = xp.empty_like(x)
                mask = (c == 0)
                t[mask] = xp.exp(-(x[mask] - mu[mask])/sigma[mask])
                t[~mask] = (
                    1  - c[~mask]*(x[~mask] - mu[~mask])/sigma[~mask]
                )**(1/c[~mask])
                return xp.exp(-t)[()]

        GenExtreme1 = stats.make_distribution(MyGenExtreme())
        GenExtreme2 = stats.make_distribution(stats.genextreme)

        X1 = GenExtreme1(c=xp.asarray(c), mu=0, sigma=1)
        X2 = GenExtreme2(c=c)

        x = X2.sample(shape=10, rng=rng)
        p = X2.cdf(x)
        x_xp, p_xp = xp.asarray(x), xp.asarray(p)

        xp_assert_close(X1.support()[0], xp.asarray(X2.support()[0]))
        xp_assert_close(X1.support()[1], xp.asarray(X2.support()[1]))
        xp_assert_close(X1.entropy(), xp.asarray(X2.entropy()), rtol=5e-6)
        xp_assert_close(X1.median(), xp.asarray(X2.median()))
        xp_assert_close(X1.logpdf(x_xp), xp.asarray(X2.logpdf(x)))
        xp_assert_close(X1.pdf(x_xp), xp.asarray(X2.pdf(x)))
        xp_assert_close(X1.logcdf(x_xp), xp.asarray(X2.logcdf(x)))
        xp_assert_close(X1.cdf(x_xp), xp.asarray(X2.cdf(x)))
        xp_assert_close(X1.logccdf(x_xp), xp.asarray(X2.logccdf(x)), rtol=1e-7)
        xp_assert_close(X1.ccdf(x_xp), xp.asarray(X2.ccdf(x)))
        xp_assert_close(X1.icdf(p_xp), xp.asarray(X2.icdf(p)))
        xp_assert_close(X1.iccdf(p_xp), xp.asarray(X2.iccdf(p)))

    @pytest.mark.slow
    @pytest.mark.parametrize("a", [0.5, np.asarray([0.5, 1.0, 2.0, 4.0, 8.0])])
    @pytest.mark.parametrize("b", [0.5, np.asarray([0.5, 1.0, 2.0, 4.0, 8.0])])
    def test_custom_multiple_parameterizations(self, a, b, xp):
        rng = np.random.default_rng(7548723590230982)
        class MyBeta:
            @property
            def __make_distribution_version__(self):
                return "1.16.0"

            @property
            def parameters(self):
                return (
                    {"a": (0, xp.inf), "b": (0, xp.inf)},
                    {"mu": (0, 1), "nu": (0, xp.inf)},
                )

            def process_parameters(self, a=None, b=None, mu=None, nu=None):
                if a is not None and b is not None and mu is None and nu is None:
                    nu = a + b
                    mu = a / nu
                else:
                    a = mu * nu
                    b = nu - a
                return {"a": a, "b": b, "mu": mu, "nu": nu}

            @property
            def support(self):
                return {'endpoints': (0, 1)}

            def pdf(self, x, a, b, mu, nu):
                return x**(a-1)*(1-x)**(b-1) / special.beta(a, b)

            def cdf(self, x, a, b, mu, nu):
                return special.betainc(a, b, x)

        Beta = stats.make_distribution(stats.beta)
        MyBeta = stats.make_distribution(MyBeta())

        mu = a / (a + b)
        nu = a + b

        Z = Beta(a=a, b=b)

        dtype = xp.float64  # stats.beta doesn't do float32
        X = MyBeta(a=xp.asarray(a, dtype=dtype), b=xp.asarray(b, dtype=dtype))
        Y = MyBeta(mu=xp.asarray(mu, dtype=dtype), nu=xp.asarray(nu, dtype=dtype))

        x = Z.sample(shape=10, rng=rng)
        p = Z.cdf(x)
        x_xp, p_xp = xp.asarray(x), xp.asarray(p)

        xp_assert_close(X.support()[0], xp.asarray(xp_copy(Z.support()[0])))
        xp_assert_close(X.support()[1], xp.asarray(xp_copy(Z.support()[1])))
        xp_assert_close(X.median(), xp.asarray(Z.median()))
        xp_assert_close(X.pdf(x_xp), xp.asarray(Z.pdf(x)))
        xp_assert_close(X.cdf(x_xp), xp.asarray(Z.cdf(x)))
        xp_assert_close(X.ccdf(x_xp), xp.asarray(Z.ccdf(x)))
        xp_assert_close(X.icdf(p_xp), xp.asarray(Z.icdf(p)))
        xp_assert_close(X.iccdf(p_xp), xp.asarray(Z.iccdf(p)))

        xp_assert_close(Y.support()[0], xp.asarray(xp_copy(Z.support()[0])))
        xp_assert_close(Y.support()[1], xp.asarray(xp_copy(Z.support()[1])))
        xp_assert_close(Y.median(), xp.asarray(Z.median()))
        xp_assert_close(Y.pdf(x_xp), xp.asarray(Z.pdf(x)))
        xp_assert_close(Y.cdf(x_xp), xp.asarray(Z.cdf(x)))
        xp_assert_close(Y.ccdf(x_xp), xp.asarray(Z.ccdf(x)))
        xp_assert_close(Y.icdf(p_xp), xp.asarray(Z.icdf(p)))
        xp_assert_close(Y.iccdf(p_xp), xp.asarray(Z.iccdf(p)))

    @pytest.mark.parametrize("n", [10, [10., 11., 12.]])
    @pytest.mark.parametrize("p", [0.4, [0.4, 0.5, 0.6]])
    def test_custom_discrete(self, n, p, xp):
        rng = np.random.default_rng(5875373614)

        class MyBinomial:
            @property
            def __make_distribution_version__(self):
                return "1.16.0"

            @property
            def parameters(self):
                return {'n': {'domain_type': 'discrete',
                              'endpoints': (0., np.inf),
                              'inclusive': (True, False)},
                        'p': {'domain_type': 'continuous',
                              'endpoints': (0., 1.),
                              'inclusive': (True, True)}}

            @property
            def support(self):
                return {'endpoints': (0, 'n'), 'inclusive': (True, True)}

            def pmf(self, x, n, p):
                return special.binom(n, x) * p**x * (1 - p)**(n - x)

        Y = stats.make_distribution(stats.binom)(n=n, p=p)
        Binomial = stats.make_distribution(MyBinomial())
        dtype = xp.float64  # binom doesn't do float32
        n, p = xp.asarray(n, dtype=dtype), xp.asarray(p, dtype=dtype)
        X = Binomial(n=n, p=p)

        x = Y.sample(shape=10, rng=rng)
        p = rng.random(x.shape)
        x_np, p_np = xp.asarray(x), xp.asarray(p)

        xp_assert_close(X.support()[0], xp.asarray(xp_copy(Y.support()[0])))
        xp_assert_close(X.support()[1], xp.asarray(xp_copy(Y.support()[1])))
        xp_assert_close(X.entropy(), xp.asarray(Y.entropy()))
        xp_assert_close(X.median(), xp.asarray(Y.median()))
        xp_assert_close(X.logpdf(x_np), xp.asarray(Y.logpdf(x)))
        xp_assert_close(X.pdf(x_np), xp.asarray(Y.pdf(x)))
        xp_assert_close(X.logpmf(x_np), xp.asarray(Y.logpmf(x)))
        xp_assert_close(X.pmf(x_np), xp.asarray(Y.pmf(x)))
        xp_assert_close(X.logcdf(x_np), xp.asarray(Y.logcdf(x)))
        xp_assert_close(X.cdf(x_np), xp.asarray(Y.cdf(x)))
        xp_assert_close(X.logccdf(x_np), xp.asarray(Y.logccdf(x)))
        xp_assert_close(X.ccdf(x_np), xp.asarray(Y.ccdf(x)))
        xp_assert_close(X.icdf(p_np), xp.asarray(Y.icdf(p)))
        xp_assert_close(X.iccdf(p_np), xp.asarray(Y.iccdf(p)))

    def test_input_validation(self):
        message = '`levy_stable` is not supported.'
        with pytest.raises(NotImplementedError, match=message):
            stats.make_distribution(stats.levy_stable)

        message = '`vonmises` is not supported.'
        with pytest.raises(NotImplementedError, match=message):
            stats.make_distribution(stats.vonmises)

        message = "The argument must be an instance of..."
        with pytest.raises(ValueError, match=message):
            stats.make_distribution(object())

        message = "If specified, the `domain_type` value..."
        class MyTestDistribution:
            __make_distribution_version__ = "1.16.0"
            parameters = {'n': {'domain_type': 'other', 'endpoints': (0, np.inf)}}
            support = {'endpoints': (0, 'n'), 'inclusive': (True, True)}
            def pmf(self, x, n):
                return np.full_like(x, 1/n)
        with pytest.raises(ValueError, match=message):
            stats.make_distribution(MyTestDistribution())(n=10)

        message = "The argument of `make_distribution` must implement at least..."
        class MyTestDistribution:
            __make_distribution_version__ = "1.16.0"
            parameters = {'n': {'endpoints': (0, np.inf)}}
            support = {'endpoints': (0, 'n'), 'inclusive': (True, True)}
        with pytest.raises(ValueError, match=message):
            stats.make_distribution(MyTestDistribution())(n=10)

    def test_repr_str_docs(self):
        from scipy.stats._univariate_distribution_extensions import _distribution_names
        for dist in _distribution_names.keys():
            assert hasattr(stats, dist)

        dist = stats.make_distribution(stats.gamma)
        assert str(dist(a=2)) == "Gamma(a=2)"
        assert repr(dist(a=2)) == "Gamma(a=2)"
        assert 'Gamma' in dist.__doc__

        dist = stats.make_distribution(stats.halfgennorm)
        assert str(dist(beta=2)) == "HalfGeneralizedNormal(beta=2)"
        assert repr(dist(beta=2)) == "HalfGeneralizedNormal(beta=2)"
        assert 'HalfGeneralizedNormal' in dist.__doc__

    @pytest.mark.slow  # just in case
    @pytest.mark.parametrize('seed', 77760325 + np.arange(20))
    def test_draw_distribution(self, seed, xp):
        # `draw_distribution_from_family` is a private function right now, but we may
        # want that functionality to be public someday. It was broken for custom
        # distributions because the `typical` parameter of the support was ignored.
        # Check that this is resolved.
        rng = np.random.default_rng(seed)
        u_typical = tuple(np.sort(rng.standard_normal(2)))
        s_typical = tuple(np.sort(rng.random(2)*2))
        x_typical = tuple(np.sort(rng.standard_normal(2)))

        class MyNormal:
            __make_distribution_version__ = "1.16.0"
            parameters = {'u': {'endpoints': (-np.inf, np.inf), 'typical': u_typical},
                          's': {'endpoints': (0, np.inf), 'typical': s_typical}}
            support = {'endpoints': (-np.inf, np.inf), 'typical': x_typical}

            def pdf(self, x, u, s):
                return 1 / xp.sqrt(2*xp.pi) / s * xp.exp(-((x-u)/s)**2/2)

        family = stats.make_distribution(MyNormal())
        proportions = (1.0, 0., 0., 0.)
        tmp = draw_distribution_from_family(family, rng, {'min_side': 1},
                                            proportions, xp=xp)
        dist, x = tmp.dist, tmp.x
        assert u_typical[0] < xp.min(dist.u) and xp.max(dist.u) < u_typical[1]
        assert s_typical[0] < xp.min(dist.s) and xp.max(dist.s) < s_typical[1]
        assert x_typical[0] < xp.min(x) and xp.max(x) < x_typical[1]

    @pytest.mark.parametrize('p', [None, ()])
    def test_no_parameters(self, p, xp):
        # To create a distribution without parameters, it is natural to try an empty
        # dictionary (since a dictionary with entries is valid if there is only
        # one parameterization), but that didn't work; only an empty tuple (no
        # parameterizations) worked originally. Check that this is resolved.
        class MyStandardNormal:
            __make_distribution_version__ = "1.16.0"
            parameters = {} if p is None else p
            support = {'endpoints': (-xp.inf, xp.inf)}

            def pdf(self, x):
                return 1 / xp.sqrt(2*xp.pi) * xp.exp(-x**2/2)

        StandardNormal = stats.make_distribution(MyStandardNormal())
        X = StandardNormal(xp=xp)
        assert X.support() == (-xp.inf, xp.inf)


@make_xp_test_case(stats.truncate)
class TestTruncate:
    def test_truncate(self, xp):
        rng = np.random.default_rng(81345982345826)
        lb = rng.random((3, 1))
        ub = rng.random((3, 1))
        lb, ub = np.minimum(lb, ub), np.maximum(lb, ub)
        Y0 = stats.truncnorm(lb, ub)

        lb, ub = xp.asarray(lb), xp.asarray(ub)
        Y = stats.truncate(Normal(xp=xp), lb=lb, ub=ub)

        y = Y0.rvs((3, 10), random_state=rng)
        p = Y0.cdf(y)
        y_xp, p_xp = xp.asarray(y), xp.asarray(p)

        xp_assert_close(Y.logentropy(), xp.log(xp.asarray(Y0.entropy() + 0j)))
        xp_assert_close(Y.entropy(), xp.asarray(Y0.entropy()))
        xp_assert_close(Y.median(), xp.asarray(Y0.ppf(0.5)))
        xp_assert_close(Y.mean(), xp.asarray(Y0.mean()))
        xp_assert_close(Y.variance(), xp.asarray(Y0.var()))
        xp_assert_close(Y.standard_deviation(), xp.asarray(Y0.var()**0.5))
        xp_assert_close(Y.skewness(), xp.asarray(Y0.stats('s')))
        xp_assert_close(Y.kurtosis(), xp.asarray(Y0.stats('k') + 3))
        xp_assert_close(Y.support()[0], xp.asarray(Y0.support()[0]))
        xp_assert_close(Y.support()[1], xp.asarray(Y0.support()[1]))
        xp_assert_close(Y.pdf(y_xp), xp.asarray(Y0.pdf(y)))
        xp_assert_close(Y.cdf(y_xp), xp.asarray(Y0.cdf(y)))
        xp_assert_close(Y.ccdf(y_xp), xp.asarray(Y0.sf(y)))
        xp_assert_close(Y.icdf(p_xp), xp.asarray(Y0.ppf(p)))
        xp_assert_close(Y.iccdf(p_xp), xp.asarray(Y0.isf(p)))
        xp_assert_close(Y.logpdf(y_xp), xp.asarray(Y0.logpdf(y)))
        xp_assert_close(Y.logcdf(y_xp), xp.asarray(Y0.logcdf(y)))
        xp_assert_close(Y.logccdf(y_xp), xp.asarray(Y0.logsf(y)))
        xp_assert_close(Y.ilogcdf(xp.log(p_xp)), xp.asarray(Y0.ppf(p)))
        xp_assert_close(Y.ilogccdf(xp.log(p_xp)), xp.asarray(Y0.isf(p)))
        sample = Y.sample(10)
        assert xp.all((sample > lb) & (sample < ub))


@make_xp_test_case(stats.log)
class TestLog:
    @pytest.mark.fail_slow(5)
    def test_log(self, xp):
        rng = np.random.default_rng(81345982345826)
        a = rng.random((3, 1))

        Y0 = stats.loggamma(a)
        X = _Gamma(a=xp.asarray(a))
        Y = stats.log(X)
        y = Y0.rvs((3, 10), random_state=rng)
        p = Y0.cdf(y)
        y_xp, p_xp = xp.asarray(y), xp.asarray(p)

        xp_assert_close(Y.logentropy(), xp.asarray(np.log(Y0.entropy()+0j)))
        xp_assert_close(Y.entropy(), xp.asarray(Y0.entropy()))
        xp_assert_close(Y.median(), xp.asarray(Y0.ppf(0.5)))
        xp_assert_close(Y.mean(), xp.asarray(Y0.mean()))
        xp_assert_close(Y.variance(), xp.asarray(Y0.var()))
        xp_assert_close(Y.standard_deviation(), xp.asarray(np.sqrt(Y0.var())))
        xp_assert_close(Y.skewness(), xp.asarray(Y0.stats('s')))
        xp_assert_close(Y.kurtosis(), xp.asarray(Y0.stats('k') + 3))
        xp_assert_close(Y.support()[0], xp.asarray(Y0.support()[0]))
        xp_assert_close(Y.support()[1], xp.asarray(Y0.support()[1]))
        xp_assert_close(Y.pdf(y_xp), xp.asarray(Y0.pdf(y)))
        xp_assert_close(Y.cdf(y_xp), xp.asarray(Y0.cdf(y)))
        xp_assert_close(Y.ccdf(y_xp), xp.asarray(Y0.sf(y)))
        xp_assert_close(Y.icdf(p_xp), xp.asarray(Y0.ppf(p)))
        xp_assert_close(Y.iccdf(p_xp), xp.asarray(Y0.isf(p)))
        xp_assert_close(Y.logpdf(y_xp), xp.asarray(Y0.logpdf(y)))
        xp_assert_close(Y.logcdf(y_xp), xp.asarray(Y0.logcdf(y)))
        xp_assert_close(Y.logccdf(y_xp), xp.asarray(Y0.logsf(y)))
        with np.errstate(invalid='ignore'):
            xp_assert_close(Y.ilogcdf(xp.log(p_xp)), xp.asarray(Y0.ppf(p)))
            xp_assert_close(Y.ilogccdf(xp.log(p_xp)), xp.asarray(Y0.isf(p)))
        seed = 3984593485
        xp_assert_close(Y.sample(rng=seed), xp.asarray(np.log(X.sample(rng=seed))))


@make_xp_test_case(stats.exp)
class TestExpPow:
    @pytest.mark.fail_slow(5)
    @pytest.mark.parametrize('exp_pow', ['exp', 'pow'])
    def test_exp_pow(self, exp_pow, xp):
        rng = np.random.default_rng(81345982345826)
        mu = rng.random((3, 1))
        sigma = rng.random((3, 1))

        X = Normal(xp=xp)*xp.asarray(sigma) + xp.asarray(mu)
        if exp_pow == 'exp':
            Y = stats.exp(X)
        else:
            Y = xp.e ** X
        Y0 = stats.lognorm(sigma, scale=np.exp(mu))

        y = Y0.rvs((3, 10), random_state=rng)
        p = Y0.cdf(y)
        y_xp, p_xp = xp.asarray(y), xp.asarray(p)

        xp_assert_close(Y.logentropy(), xp.asarray(np.log(Y0.entropy()) + 0j))
        xp_assert_close(Y.entropy(), xp.asarray(Y0.entropy()))
        xp_assert_close(Y.median(), xp.asarray(Y0.ppf(0.5)))
        xp_assert_close(Y.mean(), xp.asarray(Y0.mean()))
        xp_assert_close(Y.variance(), xp.asarray(Y0.var()))
        xp_assert_close(Y.standard_deviation(), xp.asarray(Y0.var()**0.5))
        xp_assert_close(Y.skewness(), xp.asarray(Y0.stats('s')))
        xp_assert_close(Y.kurtosis(), xp.asarray(Y0.stats('k') + 3))
        xp_assert_close(Y.support()[0], xp.asarray(Y0.support()[0]))
        xp_assert_close(Y.support()[1], xp.asarray(Y0.support()[1]))
        xp_assert_close(Y.pdf(y_xp), xp.asarray(Y0.pdf(y)))
        xp_assert_close(Y.cdf(y_xp), xp.asarray(Y0.cdf(y)))
        xp_assert_close(Y.ccdf(y_xp), xp.asarray(Y0.sf(y)))
        xp_assert_close(Y.icdf(p_xp), xp.asarray(Y0.ppf(p)))
        xp_assert_close(Y.iccdf(p_xp), xp.asarray(Y0.isf(p)))
        xp_assert_close(Y.logpdf(y_xp), xp.asarray(Y0.logpdf(y)))
        xp_assert_close(Y.logcdf(y_xp), xp.asarray(Y0.logcdf(y)))
        xp_assert_close(Y.logccdf(y_xp), xp.asarray(Y0.logsf(y)))
        xp_assert_close(Y.ilogcdf(xp.log(p_xp)), xp.asarray(Y0.ppf(p)))
        xp_assert_close(Y.ilogccdf(xp.log(p_xp)), xp.asarray(Y0.isf(p)))
        seed = 3984593485
        xp_assert_close(Y.sample(rng=seed), xp.asarray(np.exp(X.sample(rng=seed))))


class MixedDist(ContinuousDistribution):
    _variable = _RealParameter('x', domain=_RealInterval(endpoints=(-inf, inf)))
    def _pdf_formula(self, x, *args, **kwargs):
        xp = self._xp
        return (0.4 * 1/(1.1 * math.sqrt(2*math.pi)) * xp.exp(-0.5*((x+0.25)/1.1)**2)
                + 0.6 * 1/(0.9 * math.sqrt(2*math.pi)) * xp.exp(-0.5*((x-0.5)/0.9)**2))


@make_xp_test_case(stats.Mixture)
class TestMixture:
    def test_input_validation(self, xp):
        message = "`components` must contain at least one random variable."
        with pytest.raises(ValueError, match=message):
            Mixture([])

        message = "Each element of `components` must be an instance..."
        with pytest.raises(ValueError, match=message):
            Mixture((1, 2, 3))

        message = "All elements of `components` must have scalar shapes."
        with pytest.raises(ValueError, match=message):
            Mixture([Normal(mu=[1, 2], xp=xp), Normal(xp=xp)])

        message = "`components` and `weights` must have the same length."
        with pytest.raises(ValueError, match=message):
            Mixture([Normal(xp=xp)], weights=[0.5, 0.5])

        message = "`weights` must have floating point dtype."
        with pytest.raises(ValueError, match=message):
            Mixture([Normal(xp=xp)], weights=[1])

        message = "`weights` must have floating point dtype."
        with pytest.raises(ValueError, match=message):
            Mixture([Normal(xp=xp)], weights=[1])

        message = "`weights` must sum to 1.0."
        with pytest.raises(ValueError, match=message):
            Mixture([Normal(xp=xp), Normal(xp=xp)], weights=[0.5, 1.0])

        message = "All `weights` must be non-negative."
        with pytest.raises(ValueError, match=message):
            Mixture([Normal(xp=xp), Normal(xp=xp)], weights=[1.5, -0.5])

    @pytest.mark.parametrize('shape', [(), (10,)])
    def test_basic(self, shape, xp):
        rng = np.random.default_rng(582348972387243524)
        X = Mixture((Normal(mu=-0.25, sigma=1.1, xp=xp),
                     Normal(mu=0.5, sigma=0.9, xp=xp)),
                    weights=(0.4, 0.6))
        Y = MixedDist(xp=xp)
        x = xp.asarray(rng.random(shape), dtype=X._dtype)

        xp_assert_close(X.logentropy(), Y.logentropy())
        xp_assert_close(X.entropy(), Y.entropy())
        xp_assert_close(X.mode(), Y.mode())
        xp_assert_close(X.median(), Y.median())
        xp_assert_close(X.mean(), Y.mean())
        xp_assert_close(X.variance(), Y.variance())
        xp_assert_close(X.standard_deviation(), Y.standard_deviation())
        xp_assert_close(X.skewness(), Y.skewness())
        xp_assert_close(X.kurtosis(), Y.kurtosis())
        xp_assert_close(X.logpdf(x), Y.logpdf(x))
        xp_assert_close(X.pdf(x), Y.pdf(x))
        xp_assert_close(X.logcdf(x), Y.logcdf(x))
        xp_assert_close(X.cdf(x), Y.cdf(x))
        xp_assert_close(X.logccdf(x), Y.logccdf(x))
        xp_assert_close(X.ccdf(x), Y.ccdf(x))
        xp_assert_close(X.ilogcdf(x), Y.ilogcdf(x))
        xp_assert_close(X.icdf(x), Y.icdf(x))
        xp_assert_close(X.ilogccdf(x), Y.ilogccdf(x))
        xp_assert_close(X.iccdf(x), Y.iccdf(x))
        for kind in ['raw', 'central', 'standardized']:
            for order in range(5):
                xp_assert_close(X.moment(order, kind=kind),
                                Y.moment(order, kind=kind),
                                atol=1e-15 if X._dtype == xp.float64 else 1e-7)
        message = "L-moments are not currently available..."
        with pytest.raises(NotImplementedError, match=message):
            X.lmoment(1)

        # weak test of `sample`
        if is_numpy(xp):
            shape = (10, 20, 5)
            y = X.sample(shape, rng=rng)
            assert y.shape == shape
            assert stats.ks_1samp(y.ravel(), X.cdf).pvalue > 0.05

    def test_default_weights(self, xp):
        a = 1.1
        X = _Gamma(a=a, xp=xp)
        Y = stats.Mixture((X, -X))
        x = np.linspace(-4, 4, 300)
        xp_assert_close(Y.pdf(xp.asarray(x, dtype=X._dtype)),
                        xp.asarray(stats.dgamma(a=a).pdf(x), dtype=X._dtype))

    def test_properties(self, xp):
        components = [Normal(mu=-0.25, sigma=1.1, xp=xp),
                      Normal(mu=0.5, sigma=0.9, xp=xp)]
        weights = (0.4, 0.6)
        X = Mixture(components, weights=weights)

        # Replacing properties doesn't work
        # Different version of Python have different messages
        with pytest.raises(AttributeError):
            X.components = 10
        with pytest.raises(AttributeError):
            X.weights = 10

        # Mutation doesn't work
        X.components[0] = components[1]
        assert X.components[0] == components[0]
        X.weights[0] = weights[1]
        assert X.weights[0] == weights[0]

    def test_inverse(self, xp):
        # Originally, inverse relied on the mean to start the bracket search.
        # This didn't work for distributions with non-finite mean. Check that
        # this is resolved.
        rng = np.random.default_rng(24358934657854237863456)
        X0 = _StandardCauchy(xp=xp)
        X = stats.Mixture([X0, X0])
        p = xp.asarray(rng.random(size=10).tolist())
        logp = xp.log(p)
        xp_assert_close(X.icdf(p), X0.icdf(p))
        xp_assert_close(X.iccdf(p), X0.iccdf(p))
        xp_assert_close(X.ilogcdf(logp), X0.ilogcdf(logp))
        xp_assert_close(X.ilogccdf(logp), X0.ilogccdf(logp))

    def test_mixed_continuous_discrete(self, xp):
        if is_cupy(xp):
            pytest.skip("`Binomial` incompatible with CuPy.")

        X = stats.Binomial(n=0, p=0, xp=xp)
        Y = stats.Normal(xp=xp)
        Z = stats.Mixture((X, Y))

        # These are just spot-checks; reference values are supposed to be
        # obviously correct. Not the strongest tests, but they should be
        # sensitive enough to detect major regressions.
        xp_assert_close(Z.support()[0], xp.asarray(-xp.inf))
        xp_assert_close(Z.support()[1], xp.asarray(xp.inf))
        xp_assert_close(Z.mean(), xp.asarray(0.))
        xp_assert_close(Z.variance(), xp.asarray(0.5))
        xp_assert_close(Z.standard_deviation(), xp.sqrt(Z.variance()))
        xp_assert_close(Z.skewness(), xp.asarray(0.0))
        xp_assert_close(Z.kurtosis(), xp.asarray(6.0))
        xp_assert_close(Z.logpmf(0), xp.asarray(math.log(0.5)))
        xp_assert_close(Z.logpmf(1), xp.asarray(-math.inf))
        xp_assert_close(Z.pmf(0), xp.asarray(0.5))
        xp_assert_close(Z.pmf(1), xp.asarray(0.))
        xp_assert_close(Z.logpdf(0), xp.asarray(math.inf))
        xp_assert_close(Z.logpdf(1), xp.log(Y.pdf(1)/2))
        xp_assert_close(Z.pdf(0), xp.asarray(math.inf))
        xp_assert_close(Z.pdf(1), xp.asarray(Y.pdf(1)/2))
        xp_assert_close(Z.logcdf(0), xp.asarray(math.log(0.75)))
        xp_assert_close(Z.logcdf(-1e-10), xp.asarray(math.log(0.25)))
        xp_assert_close(Z.logcdf(1e-10), xp.asarray(math.log(0.75)))
        xp_assert_close(Z.cdf(0), xp.asarray(0.75))
        xp_assert_close(Z.cdf(-1e-10), xp.asarray(0.25))
        xp_assert_close(Z.cdf(1e-10), xp.asarray(0.75))
        xp_assert_close(Z.logccdf(0), xp.asarray(math.log(0.25)))
        xp_assert_close(Z.logccdf(-1e-10), xp.asarray(math.log(0.75)))
        xp_assert_close(Z.logccdf(1e-10), xp.asarray(math.log(0.25)))
        xp_assert_close(Z.ccdf(0), xp.asarray(0.25))
        xp_assert_close(Z.ccdf(-1e-10), xp.asarray(0.75))
        xp_assert_close(Z.ccdf(1e-10), xp.asarray(0.25))

        rng = np.random.default_rng(847823487109293)
        n = 100000
        sample = Z.sample(shape=n, rng=rng)
        atol = n/500
        xp_assert_close(xp.count_nonzero(sample == 0)[()], xp.asarray(n//2), atol=atol)
        xp_assert_close(xp.count_nonzero(sample > 0)[()], xp.asarray(n//4), atol=atol)
        xp_assert_close(xp.count_nonzero(sample < 0)[()], xp.asarray(n//4), atol=atol)

    @pytest.mark.parametrize('method', ['median', 'mode', 'entropy', 'logentropy',
                                        'icdf', 'iccdf', 'ilogcdf', 'ilogccdf'])
    def test_mixed_not_implemented(self, method, xp):
        X = stats.Binomial(n=0, p=0, xp=xp)
        Y = stats.Normal(xp=xp)
        Z = stats.Mixture((X, Y))
        args = (0.,) if method.startswith('i') else ()
        f = getattr(Z, method)
        message = f"`{method}` is implemented only..."
        with pytest.raises(NotImplementedError, match=message):
            f(*args)


@make_xp_test_case(stats.order_statistic)
class TestOrderStatistic:
    @pytest.mark.fail_slow(20)  # Moments require integration
    def test_order_statistic(self, xp):
        rng = np.random.default_rng(7546349802439582)
        dtype = xp.float64
        X = Uniform(a=0, b=1, dtype=dtype, xp=xp)
        n = 5
        r = np.asarray([[1], [3], [5]])
        Y = stats.order_statistic(X, n=n, r=xp.asarray(r))
        Y0 = stats.beta(r, n + 1 - r)

        y = Y0.rvs((3, 10), random_state=rng)
        p = Y0.cdf(y)
        y_xp, p_xp = xp.asarray(y), xp.asarray(p)

        # log methods need some attention before merge
        xp_assert_close(xp.exp(Y.logentropy()), xp.asarray(Y0.entropy() + 0j))
        xp_assert_close(Y.entropy(), xp.asarray(Y0.entropy()))
        xp_assert_close(Y.mean(), xp.asarray(Y0.mean()))
        xp_assert_close(Y.variance(), xp.asarray(Y0.var()))
        xp_assert_close(Y.skewness(), xp.asarray(Y0.stats('s')), atol=1e-15)
        xp_assert_close(Y.kurtosis(), xp.asarray(Y0.stats('k') + 3), atol=1e-15)
        xp_assert_close(Y.median(), xp.asarray(Y0.ppf(0.5)))
        xp_assert_close(Y.support()[0], xp.asarray(Y0.support()[0]))
        xp_assert_close(Y.support()[1], xp.asarray(Y0.support()[1]))
        xp_assert_close(Y.pdf(y_xp), xp.asarray(Y0.pdf(y)))
        xp_assert_close(Y.cdf(y_xp, method='formula'),
                        xp.asarray(Y.cdf(y_xp, method='quadrature')))
        xp_assert_close(Y.ccdf(y_xp, method='formula'),
                        xp.asarray(Y.ccdf(y_xp, method='quadrature')))
        xp_assert_close(Y.icdf(p_xp, method='formula'),
                        xp.asarray(Y.icdf(p_xp, method='inversion')))
        if is_numpy(xp):  # CuPy needs betainccinv
            xp_assert_close(Y.iccdf(p_xp, method='formula'),
                            xp.asarray(Y.iccdf(p, method='inversion')))
        xp_assert_close(Y.logpdf(y_xp), xp.asarray(Y0.logpdf(y)))
        xp_assert_close(Y.logcdf(y_xp), xp.asarray(Y0.logcdf(y)))
        xp_assert_close(Y.logccdf(y_xp), xp.asarray(Y0.logsf(y)))
        with np.errstate(invalid='ignore', divide='ignore'):
            xp_assert_close(Y.ilogcdf(xp.log(p_xp),), xp.asarray(Y0.ppf(p)))
            xp_assert_close(Y.ilogccdf(xp.log(p_xp)), xp.asarray(Y0.isf(p)))

        message = "`r` and `n` must contain only positive integers."
        with pytest.raises(ValueError, match=message):
            stats.order_statistic(X, n=n, r=-1)
        with pytest.raises(ValueError, match=message):
            stats.order_statistic(X, n=-1, r=r)
        with pytest.raises(ValueError, match=message):
            stats.order_statistic(X, n=n, r=1.5)
        with pytest.raises(ValueError, match=message):
            stats.order_statistic(X, n=1.5, r=r)

    def test_support_gh22037(self, xp):
        # During review of gh-22037, it was noted that the `support` of
        # an `OrderStatisticDistribution` returned incorrect results;
        # this was resolved by overriding `_support`.
        if is_numpy(xp):
            Uniform = stats.make_distribution(stats.uniform)
            X = Uniform()
        else:
            X = stats.Uniform(a=0, b=1, xp=xp)
        Y = X*5 + 2
        Z = stats.order_statistic(Y, r=3, n=5)
        xp_assert_close(Z.support()[0], Y.support()[0])
        xp_assert_close(Z.support()[1], Y.support()[1])

    def test_composition_gh22037(self, xp):
        # During review of gh-22037, it was noted that an error was
        # raised when creating an `OrderStatisticDistribution` from
        # a `TruncatedDistribution`. This was resolved by overriding
        # `_update_parameters`.
        if is_numpy(xp):
            Normal = stats.make_distribution(stats.norm)
        else:
            Normal = stats.Normal
        TruncatedNormal = stats.make_distribution(stats.truncnorm)
        a, b = [-2, -1], 1
        r, n = 3, [[4], [5]]
        x = [[[-0.3]], [[0.1]]]
        X1 = Normal(xp=xp)
        Y1 = stats.truncate(X1, a, b)
        Z1 = stats.order_statistic(Y1, r=r, n=n)
        X2 = TruncatedNormal(a=a, b=b)
        Z2 = stats.order_statistic(X2, r=r, n=n)
        xp_assert_close(Z1.cdf(xp.asarray(x)), xp.asarray(Z2.cdf(x)).tolist())

