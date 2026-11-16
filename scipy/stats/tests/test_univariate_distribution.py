import operator
import pickle
from copy import deepcopy

import numpy as np
from numpy import inf
import pytest
from numpy.testing import assert_equal

from scipy import special
from scipy import stats
from scipy.stats._fit import _kolmogorov_smirnov
from scipy.stats._ksstats import kolmogn
from scipy.stats import qmc
from scipy.stats._univariate_distribution_util import (
    _RealInterval, _Parameterization, _RealParameter,
    _generate_domain_support)
from scipy.stats._univariate_distribution import (
    ContinuousDistribution, DiscreteDistribution, CircularDistribution,
    ShiftedScaledDistribution, TransformedDistribution,
    MonotonicTransformedDistribution)
from scipy.stats._new_distributions import (
    StandardNormal, _LogUniform, _Gamma, _TestCircular)
from scipy.stats import Normal, Uniform
from scipy._external import array_api_extra as xpx
from scipy._lib._testutils import mutually_broadcastable_shapes
from scipy._lib._util import _RichResult
from scipy._lib._array_api import (
    make_xp_test_case, is_numpy, is_cupy, is_torch, xp_compat_namespace)
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal, xp_assert_less


def draw_distribution_from_family(family, rng, shape_options, proportions=None, *, xp):
    # If the distribution has parameters, choose a parameterization and
    # draw broadcastable shapes for the parameter arrays.
    proportions = (0.7, 0.1, 0.1, 0.1) if proportions is None else proportions
    n_parameterizations = family._num_parameterizations()
    if n_parameterizations > 0:
        i = rng.integers(0, n_parameterizations)
        n_parameters = family._num_parameters(i)
        shapes = mutually_broadcastable_shapes(num_shapes=n_parameters,
                                               rng=rng, **shape_options)
        result_shape = xp.broadcast_shapes(*shapes)
        dist = family._draw(shapes, rng=rng, proportions=proportions,
                            i_parameterization=i, xp=xp)
    else:
        dist = family._draw(rng=rng, xp=xp)
        result_shape = tuple()

    # Draw a broadcastable shape for the arguments, and draw values for the
    # arguments.
    x_shape, = mutually_broadcastable_shapes(1, base_shape=result_shape,
                                             rng=rng, **shape_options)
    x = dist._variable.draw(x_shape, parameter_values=dist._parameters,
                            proportions=proportions, rng=rng, region='typical', xp=xp)
    x_result_shape = xp.broadcast_shapes(x_shape, result_shape)
    y_shape, = mutually_broadcastable_shapes(1, base_shape=x_result_shape,
                                             rng=rng, **shape_options)
    y = dist._variable.draw(y_shape, parameter_values=dist._parameters,
                            proportions=proportions, rng=rng, region='typical', xp=xp)
    xy_result_shape = xp.broadcast_shapes(y_shape, x_result_shape)
    p_domain = _RealInterval((0, 1), (True, True))
    p_var = _RealParameter('p', domain=p_domain)
    p = p_var.draw(x_shape, proportions=proportions, rng=rng, xp=xp)
    with np.errstate(divide='ignore', invalid='ignore'):
        logp = xp.log(p)

    return _RichResult(dist=dist, x=x, y=y, p=p, logp=logp, result_shape=result_shape,
                       x_result_shape=x_result_shape, xy_result_shape=xy_result_shape)


class DistributionsTest:
    _options = [
        (dict(min_dims=0, max_dims=0), [1, 0, 0, 0]), # all valid scalar
        (dict(min_dims=0, max_dims=0), [0, 0, 0, 1]), # all nan scalar
        (dict(min_dims=1, max_dims=1,
            min_side=5, max_side=6), [1, 0, 0, 0]), # all valid array
        (dict(min_dims=1, max_dims=1,
            min_side=5, max_side=6), [0, 0, 0, 1]), # all nan array
        (dict(min_dims=1, max_dims=1,  # mixed valid, invalid, and edge cases
            min_side=20, max_side=25), [0.25, 0.25, 0.25, 0.25]),
    ] + [({}, None)]*15  # other, randomly-generated shape options

    @pytest.fixture(params=list(enumerate(_options)), ids=range(len(_options)))
    def options(self, request):
        return request.param

    @pytest.fixture
    def case(self, options, xp):
        i, _options = options
        shape_options, proportions = _options
        rng = np.random.default_rng([i, self.seed])
        tmp = draw_distribution_from_family(self.family, rng,
                                            shape_options, proportions, xp=xp)
        return _RichResult(family=self.family, rng=rng, **tmp)

    @pytest.fixture
    def valid_dist_x(self, xp):
        shape_options = dict(min_dims=1, max_dims=1, min_side=20, max_side=21)
        proportions = (1, 0, 0, 0)  # all valid
        rng = np.random.default_rng(self.seed)
        tmp = draw_distribution_from_family(self.family, rng,
                                            shape_options, proportions, xp=xp)
        return tmp.dist, tmp.x

    def test_purported_distribution(self, valid_dist_x):
        message = ("This method must be overridden to demonstrate that the "
                   "distribution is the purported one.")
        raise NotImplementedError(message)

    def test_support(self, case, *, xp):
        check_support(case.dist, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case, *, tol_override=None, xp):
        check_moment_funcs(case.dist, case.result_shape,
                           tol_override=tol_override, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_lmoment(self, case, *, tol_override=None, xp):
        if isinstance(case.dist, CircularDistribution):
            pytest.skip("`CircularDistributions` don't support `lmoment`.")
        check_lmoment_funcs(case.dist, case.result_shape,
                            tol_override=tol_override, xp=xp)

    def test_random_sample(self, case, *, xp):
        sample_shape, = mutually_broadcastable_shapes(1, max_side=20, rng=case.rng)
        check_sample_shape_NaNs(case.dist, 'sample', sample_shape,
                                case.result_shape, case.rng, xp=xp)

    @pytest.mark.thread_unsafe(reason="looks like an _rng_spawn issue?")
    def test_quasi_random_sample(self, case, *, xp):
        sample_shape, = mutually_broadcastable_shapes(1, max_side=20, rng=case.rng)
        qrng = qmc.Halton(d=1, seed=case.rng)
        check_sample_shape_NaNs(case.dist, 'sample', sample_shape,
                                case.result_shape, qrng, xp=xp)

    def test_entropy(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'entropy', None, case.result_shape,
                        {'log/exp', 'quadrature'},
                        tol_override=tol_override, xp=xp)

    def test_logentropy(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'logentropy', None, case.result_shape,
                        {'log/exp', 'quadrature'},
                        tol_override=tol_override, xp=xp)

    def test_median(self, case, *, tol_override=None, xp):
        methods = ({'optimization'} if isinstance(case.dist, CircularDistribution)
                   else {'icdf'})
        check_dist_func(case.dist, 'median', None, case.result_shape, methods,
                        tol_override=tol_override, xp=xp)

    def test_mode(self, case, *, tol_override=None, xp):
        default_atol = 1e-6 if case.dist._dtype == xp.float64 else 1e-4
        tol_override = {'atol': default_atol} if tol_override is None else tol_override
        check_dist_func(case.dist, 'mode', None, case.result_shape,
                        {'optimization'}, tol_override=tol_override, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_mean(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'mean', None, case.result_shape, {'cache'},
                        tol_override=tol_override, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_variance(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'variance', None, case.result_shape, {'cache'},
                        tol_override=tol_override, xp=xp)

    def test_standard_deviation(self, case, *, tol_override=None, xp):
        tol_override = {} if tol_override is None else tol_override
        std = case.dist.standard_deviation()
        res = (1 - xp.exp(-0.5*std**2) if isinstance(case.dist, CircularDistribution)
               else std**2)
        xp_assert_close(res, case.dist.variance(), **tol_override)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_skewness(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'skewness', None, case.result_shape,
                        {'cache'}, tol_override=tol_override, xp=xp)

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_kurtosis(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'kurtosis', None, case.result_shape,
                        {'cache'}, tol_override=tol_override, xp=xp)

    def test_pdf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'pdf', case.x, case.x_result_shape,
                        {'log/exp'}, tol_override=tol_override, xp=xp)

    def test_logpdf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'logpdf', case.x, case.x_result_shape,
                        {'log/exp'}, tol_override=tol_override, xp=xp)

    def test_pmf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'pmf', case.x, case.x_result_shape,
                        {'log/exp'}, tol_override=tol_override, xp=xp)

    def test_logpmf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'logpmf', case.x, case.x_result_shape,
                        {'log/exp'}, tol_override=tol_override, xp=xp)

    def test_logcdf(self, case, *, tol_override=None, xp):
        if isinstance(case.dist, CircularDistribution):
            pytest.skip("`CircularDistributions` don't support `logcdf`.")
        check_dist_func(case.dist, 'logcdf', case.x, case.x_result_shape,
                        {'log/exp', 'complement', 'quadrature'},
                        tol_override=tol_override, xp=xp)

    def test_cdf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'cdf', case.x, case.x_result_shape,
                        {'log/exp', 'complement', 'quadrature'},
                        tol_override=tol_override, xp=xp)

    def test_cdf2(self, case, *, tol_override=None, xp):
        if isinstance(case.dist, CircularDistribution):
            pytest.skip("`CircularDistributions` don't support two-arg `cdf`.")
        check_cdf2(case.dist, False, case.x, case.y,
                    case.xy_result_shape, {'quadrature'},
                    tol_override=tol_override, xp=xp)
        check_cdf2(case.dist, True, case.x, case.y,
                    case.xy_result_shape, {'quadrature'},
                    tol_override=tol_override, xp=xp)

    def test_logccdf(self, case, *, tol_override=None, xp):
        if isinstance(case.dist, CircularDistribution):
            pytest.skip("`CircularDistributions` don't support `logccdf`.")
        check_dist_func(case.dist, 'logccdf', case.x, case.x_result_shape,
                        {'log/exp', 'complement', 'quadrature'},
                        tol_override=tol_override, xp=xp)

    def test_ccdf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'ccdf', case.x, case.x_result_shape,
                        {'log/exp', 'complement', 'quadrature'},
                        tol_override=tol_override, xp=xp)

    def test_ccdf2(self, case, *, tol_override=None, xp):
        if isinstance(case.dist, CircularDistribution):
            pytest.skip("`CircularDistributions` don't support two-arg `ccdf`.")
        check_ccdf2(case.dist, False, case.x, case.y,
                    case.xy_result_shape, {'addition'},
                    tol_override=tol_override, xp=xp)
        check_ccdf2(case.dist, True, case.x, case.y,
                    case.xy_result_shape, {'addition'},
                    tol_override=tol_override, xp=xp)

    def test_ilogcdf(self, case, *, tol_override=None, xp):
        if isinstance(case.dist, CircularDistribution):
            pytest.skip("`CircularDistributions` don't support `ilogcdf`.")
        with np.errstate(divide='ignore', over='ignore'):
            check_dist_func(case.dist, 'ilogcdf', case.logp, case.x_result_shape,
                {'complement', 'inversion'}, tol_override=tol_override, xp=xp)

    def test_icdf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'icdf', case.p, case.x_result_shape,
            {'complement', 'inversion'}, tol_override=tol_override, xp=xp)

    def test_ilogccdf(self, case, *, tol_override=None, xp):
        if isinstance(case.dist, CircularDistribution):
            pytest.skip("`CircularDistributions` don't support `ilogccdf`.")
        with np.errstate(divide='ignore', over='ignore'):
            check_dist_func(case.dist, 'ilogccdf', case.logp, case.x_result_shape,
                {'complement', 'inversion'}, tol_override=tol_override, xp=xp)

    def test_iccdf(self, case, *, tol_override=None, xp):
        check_dist_func(case.dist, 'iccdf', case.p, case.x_result_shape,
            {'complement', 'inversion'}, tol_override=tol_override, xp=xp)


# Since this distributions is not public, the tests is really included to test the
# distribution infrastructure more than the distribution itself. It also avoids
# `DistributionTest` being detected as a misnamed test class by `check_test_name.py`
@make_xp_test_case(_LogUniform)
class Test_LogUniform(DistributionsTest):
    seed = 260607439
    family = _LogUniform

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        xp_assert_close(dist.pdf(x),
                        1 / (x * (xp.log(dist.b) - xp.log(dist.a))))

    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_moment(self, case, xp):
        atol = 1e-8 if case.dist._dtype == xp.float64 else 5e-3
        return super().test_moment(case, tol_override={'atol': atol}, xp=xp)

    @pytest.mark.fail_slow(10)
    @pytest.mark.thread_unsafe(reason="tests cache of shared `case.dist`")
    def test_lmoment(self, case, xp):
        return super().test_lmoment(case, xp=xp)

    def test_gh_25180(self, xp):
        X = _LogUniform(log_a=xp.asarray(-0.28915434544814245, dtype=xp.float64),
                        log_b=xp.asarray(0.3085224670197856, dtype=xp.float64))
        actual = X.mode(method='optimization')
        assert xp.isfinite(actual)


# Since this distributions is not public, the tests is really included to test the
# distribution infrastructure more than the distribution itself. It also avoids
# `DistributionTest` being detected as a misnamed test class by `check_test_name.py`
@make_xp_test_case(_TestCircular)
class Test_TestCircular(DistributionsTest):
    seed = 586247944
    family = _TestCircular

    def test_purported_distribution(self, valid_dist_x, xp):
        dist, x = valid_dist_x
        xp_assert_close(dist.pdf(x), 0.75 * (1 - x**2))

    def test_median(self, case, xp):
        # can only expect about half precision with optimization
        atol = 1e-6 if case.dist._dtype == xp.float64 else 1e-3
        return super().test_median(case, tol_override={'atol': atol}, xp=xp)


@make_xp_test_case(stats.Normal)
class TestOtherMethods:
    @make_xp_test_case((stats.Uniform, 'plot'))
    def test_plot(self, xp):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return

        X = Uniform(a=0., b=1., xp=xp)
        ax = X.plot()
        assert ax == plt.gca()

    @pytest.mark.parametrize('method_name', ['cdf', 'ccdf'])
    def test_complement_safe(self, method_name, xp):
        X = stats.Normal(mu=1, sigma=2, xp=xp)
        X.tol = xp.finfo(X._dtype).eps * 3e3
        p = xp.asarray([1e-4, 1e-3])
        func = getattr(X, method_name)
        ifunc = getattr(X, 'i'+method_name)
        x = ifunc(p, method='formula')
        p1 = func(x, method='complement_safe')
        p2 = func(x, method='complement')
        xp_assert_equal(p1[1], p2[1])
        assert p1[0] != p2[0]
        xp_assert_close(p1[0], p[0], rtol=X.tol)

    @pytest.mark.parametrize('method_name', ['cdf', 'ccdf'])
    def test_icomplement_safe(self, method_name, xp):
        X = stats.Normal(mu=1, sigma=2, xp=xp)
        X.tol = xp.finfo(X._dtype).eps * 3e3
        p = xp.asarray([1e-4, 1e-3])
        func = getattr(X, method_name)
        ifunc = getattr(X, 'i'+method_name)
        x1 = ifunc(p, method='complement_safe')
        x2 = ifunc(p, method='complement')
        xp_assert_equal(x1[1], x2[1])
        assert x1[0] != x2[0]
        xp_assert_close(func(x1[0]), p[0], rtol=X.tol)

    def test_subtraction_safe(self, xp):
        X = stats.Normal(xp=xp)
        X.tol = 1e-12 if xpx.default_dtype(xp) == xp.float64 else 1e-6

        # Regular subtraction is fine in either tail (and of course, across tails)
        x = xp.asarray([-11, -10, 10, 11])
        y = xp.asarray([-10, -11, 11, 10])
        p0 = X.cdf(x, y, method='quadrature')
        p1 = X.cdf(x, y, method='subtraction_safe')
        p2 = X.cdf(x, y, method='subtraction')
        xp_assert_close(p1, p0, rtol=X.tol)

        if is_torch(xp):
            pass  # https://github.com/pytorch/pytorch/issues/192295
        else:
            xp_assert_equal(p2, p1)

        # Safe subtraction is needed in special cases
        x = xp.asarray([-1e-20, -1e-21, 1e-20, 1e-21, -1e-20])
        y = xp.asarray([-1e-21, -1e-20, 1e-21, 1e-20, 1e-20])

        p0 = X.pdf(0)*(y-x)
        p1 = X.cdf(x, y, method='subtraction_safe')
        p2 = X.cdf(x, y, method='subtraction')
        xp_assert_equal(p2, xp.zeros_like(x))
        xp_assert_close(p1, p0, rtol=X.tol)

    def test_logentropy_safe(self, xp):
        # simulate an `entropy` calculation over/underflowing with extreme parameters
        class _Normal(stats.Normal):
            def _entropy_formula(self, **params):
                out = xp.asarray(super()._entropy_formula(**params))
                out[0] = 0
                out[-1] = xp.inf
                return out

        X = _Normal(sigma=xp.asarray([1, 2, 3]))
        with np.errstate(divide='ignore'):
            res1 = X.logentropy(method='logexp_safe')
            res2 = X.logentropy(method='logexp')
        ref = X.logentropy(method='quadrature')
        i_fl = xp.asarray([0, -1])  # first and last
        assert xp.all(xp.isinf(res2[i_fl]))
        assert res1[1] == res2[1]
        # quadrature happens to be perfectly accurate on some platforms
        # assert res1[1] != ref[1]
        xp_assert_equal(res1[i_fl], ref[i_fl])

    def test_logcdf2_safe(self, xp):
        # test what happens when 2-arg `cdf` underflows
        X = stats.Normal(sigma=xp.asarray([1, 2, 3]))
        x = xp.asarray([-301, 1, 300])
        y = xp.asarray([-300, 2, 301])
        with np.errstate(divide='ignore'):
            res1 = X.logcdf(x, y, method='logexp_safe')
            res2 = X.logcdf(x, y, method='logexp')
        ref = X.logcdf(x, y, method='quadrature')
        i_fl = xp.asarray([0, -1])  # first and last
        assert xp.all(xp.isinf(res2[i_fl]))
        assert res1[1] == res2[1]
        # quadrature happens to be perfectly accurate on some platforms
        # assert res1[1] != ref[1]
        xp_assert_equal(res1[i_fl], ref[i_fl])

    @pytest.mark.parametrize('method_name', ['logcdf', 'logccdf'])
    def test_logexp_safe(self, method_name, xp):
        # test what happens when `cdf`/`ccdf` underflows
        X = stats.Normal(sigma=2, xp=xp)
        x = xp.asarray([-301, 1] if method_name == 'logcdf' else [301, 1])
        func = getattr(X, method_name)
        with np.errstate(divide='ignore'):
            res1 = func(x, method='logexp_safe')
            res2 = func(x, method='logexp')
        ref = func(x, method='quadrature')
        assert res1[0] == ref[0]
        assert res1[0] != res2[0]
        assert res1[1] == res2[1]
        assert res1[1] != ref[1]


def check_sample_shape_NaNs(dist, fname, sample_shape, result_shape, rng, *, xp):
    full_shape = sample_shape + result_shape
    if fname == 'sample':
        sample_method = dist.sample

    methods = {'inverse_transform'}
    if (is_numpy(xp)
            and dist._overrides(f'_{fname}_formula')
            and not isinstance(rng, qmc.QMCEngine)):
        methods.add('formula')

    for method in methods:
        res = sample_method(sample_shape, method=method, rng=rng)
        valid_parameters = xp.broadcast_to(get_valid_parameters(dist), res.shape)
        assert res.shape == full_shape
        assert res.dtype == dist._dtype

        if is_numpy(xp) and full_shape == ():
            # NumPy random makes a distinction between a 0d array and a scalar.
            # In stats, we consistently turn 0d arrays into scalars, so
            # maintain that behavior here. (With Array API arrays, this will
            # change.)
            assert np.isscalar(res)

        strict = dict(check_dtype=False, check_shape=False, check_0d=False)
        assert xp.all(xp.isfinite(res[valid_parameters]))
        xp_assert_equal(res[~valid_parameters], np.nan, **strict)

        sample1 = sample_method(sample_shape, method=method, rng=42)
        sample2 = sample_method(sample_shape, method=method, rng=42)
        if not isinstance(dist, DiscreteDistribution):
            # The idea is that it's very unlikely that the random sample
            # for a randomly chosen seed will match that for seed 42,
            # but it is not so unlikely if `dist` is a discrete distribution.
            assert not xp.any(xp.equal(res, sample1))
        xp_assert_equal(sample1, sample2)


def check_support(dist, *, xp):
    a, b = dist.support()
    check_nans_and_edges(dist, 'support', None, a, xp=xp)
    check_nans_and_edges(dist, 'support', None, b, xp=xp)
    assert a.shape == dist._shape
    assert b.shape == dist._shape
    assert a.dtype == dist._dtype
    assert b.dtype == dist._dtype


def check_dist_func(dist, fname, arg, result_shape, methods, tol_override=None, *, xp):
    # Check that all computation methods of all distribution functions agree
    # with one another, effectively testing the correctness of the generic
    # computation methods and confirming the consistency of specific
    # distributions with their pdf/logpdf.
    args = tuple() if arg is None else (arg,)
    methods = methods.copy()

    if "cache" in methods:
        # If "cache" is specified before the value has been evaluated, it
        # raises an error. After the value is evaluated, it will succeed.
        with pytest.raises(NotImplementedError):
            getattr(dist, fname)(*args, method="cache")

    ref = getattr(dist, fname)(*args)
    check_nans_and_edges(dist, fname, arg, ref, xp=xp)

    if dist._overrides(f'_{fname}_formula'):
        methods.add('formula')

    assert ref.shape == result_shape
    # Until we convert to array API, let's do the familiar thing:
    # 0d things are scalars, not arrays
    if is_numpy(xp) and result_shape == tuple():
        assert np.isscalar(ref)

    for method in methods:
        atol = 1e-6 if (method in {'complement', 'log/exp'}
                        and ref.dtype == xp.float32) else 0
        _tol_override = {'atol': atol} if tol_override is None else tol_override

        res = getattr(dist, fname)(*args, method=method)
        if 'log' in fname:
            xp_assert_close(xp.exp(res), xp.exp(ref), **_tol_override)
        else:
            if fname in {'median', 'mode'} and isinstance(dist, CircularDistribution):
                # wrap to center around 0
                a, b = dist.support()
                period = b - a
                res = xp.where(res > period/2, res - period, res)[()]
                ref = xp.where(ref > period/2, ref - period, ref)[()]
            xp_assert_close(res, ref, **_tol_override)

        # for now, make sure dtypes are consistent; later, we can check whether
        # they are correct.
        assert res.dtype == ref.dtype
        assert res.shape == result_shape
        if is_numpy(xp) and result_shape == tuple():
            assert np.isscalar(res)


def check_cdf2(dist, log, x, y, result_shape, methods, *, tol_override=None, xp):
    # Specialized test for 2-arg cdf since the interface is a bit different
    # from the other methods. Here, we'll use 1-arg cdf as a reference, and
    # since we have already checked 1-arg cdf in `check_nans_and_edges`, this
    # checks the equivalent of both `check_dist_func` and
    # `check_nans_and_edges`.
    tol_override = {'atol' : 0} if tol_override is None else tol_override
    methods = methods.copy()

    if log:
        if dist._overrides('_logcdf2_formula'):
            methods.add('formula')
        if dist._overrides('_logcdf_formula') or dist._overrides('_logccdf_formula'):
            methods.add('subtraction')
        if (dist._overrides('_cdf_formula')
                or dist._overrides('_ccdf_formula')):
            methods.add('log/exp')
    else:
        if dist._overrides('_cdf2_formula'):
            methods.add('formula')
        if dist._overrides('_cdf_formula') or dist._overrides('_ccdf_formula'):
            methods.add('subtraction')
        if (dist._overrides('_logcdf_formula')
                or dist._overrides('_logccdf_formula')):
            methods.add('log/exp')

    ref = dist.cdf(y) - dist.cdf(x)
    ref = ref + 0j if log else ref
    assert ref.shape == result_shape

    if is_numpy(xp) and result_shape == tuple():
        assert np.isscalar(ref)

    for method in methods:
        if isinstance(dist, DiscreteDistribution):
            message = ("Two argument cdf functions are currently only supported for "
                       "continuous distributions.")
            with pytest.raises(NotImplementedError, match=message):
                res = (xp.exp(dist.logcdf(x, y, method=method)) if log
                       else dist.cdf(x, y, method=method))
            continue
        with np.errstate(invalid='ignore'):
            # np.exp(np.nan) raises on some platforms?
            res = (xp.exp(dist.logcdf(x, y, method=method)) if log
                else dist.cdf(x, y, method=method))
        xp_assert_close(res, ref, **tol_override)
        assert res.dtype == ref.dtype
        assert res.shape == result_shape
        if is_numpy(xp) and result_shape == tuple():
            assert np.isscalar(res)


def check_ccdf2(dist, log, x, y, result_shape, methods, *, tol_override=None, xp):
    # Specialized test for 2-arg ccdf since the interface is a bit different
    # from the other methods. Could be combined with check_cdf2 above, but
    # writing it separately is simpler.
    tol_override = {'atol' : 0} if tol_override is None else tol_override
    methods = methods.copy()

    if dist._overrides(f'_{"log" if log else ""}ccdf2_formula'):
        methods.add('formula')

    ref = dist.cdf(x) + dist.ccdf(y)
    assert ref.shape == result_shape

    if is_numpy(xp) and result_shape == tuple():
        assert np.isscalar(ref)

    for method in methods:
        message = ("Two argument cdf functions are currently only supported for "
                   "continuous distributions.")
        if isinstance(dist, DiscreteDistribution):
            with pytest.raises(NotImplementedError, match=message):
                res = (xp.exp(dist.logccdf(x, y, method=method)) if log
                       else dist.ccdf(x, y, method=method))
            continue
        res = (xp.exp(dist.logccdf(x, y, method=method)) if log
               else dist.ccdf(x, y, method=method))
        xp_assert_close(res, ref, **tol_override)
        assert res.dtype == ref.dtype
        assert res.shape == result_shape
        if is_numpy(xp) and result_shape == tuple():
            assert np.isscalar(res)


def check_nans_and_edges(dist, fname, arg, res, *, xp):

    valid_parameters = get_valid_parameters(dist)
    if fname in {'icdf', 'iccdf'}:
        arg_domain = _RealInterval(endpoints=(0, 1), inclusive=(True, True))
    elif fname in {'ilogcdf', 'ilogccdf'}:
        arg_domain = _RealInterval(endpoints=(-inf, 0), inclusive=(True, True))
    else:
        arg_domain = dist._variable.domain

    classified_args = classify_arg(dist, arg, arg_domain, xp=xp)
    valid_parameters, *classified_args = xp.broadcast_arrays(valid_parameters,
                                                             *classified_args)
    valid_arg, endpoint_arg, outside_arg, nan_arg = classified_args
    all_valid = valid_arg & valid_parameters

    strict = dict(check_dtype=False, check_shape=False, check_0d=False)
    # Check NaN pattern and edge cases
    xp_assert_equal(res[~valid_parameters], xp.nan, **strict)
    xp_assert_equal(res[nan_arg], xp.nan, **strict)

    a, b = dist.support()
    a = xp.broadcast_to(a, res.shape)
    b = xp.broadcast_to(b, res.shape)

    if isinstance(dist, CircularDistribution):
        return check_periodicity(dist, fname, arg, res, a, b, valid_parameters,
                                 classified_args, xp=xp)

    outside_arg_minus = (outside_arg == -1) & valid_parameters
    outside_arg_plus = (outside_arg == 1) & valid_parameters
    endpoint_arg_minus = (endpoint_arg == -1) & valid_parameters
    endpoint_arg_plus = (endpoint_arg == 1) & valid_parameters

    is_discrete = isinstance(dist, DiscreteDistribution)
    # Writing this independently of how the are set in the distribution
    # infrastructure. That is very compact; this is very verbose.
    if fname in {'logpdf'}:
        xp_assert_equal(res[outside_arg_minus], -np.inf, **strict)
        xp_assert_equal(res[outside_arg_plus], -np.inf, **strict)
        ref = -np.inf if not is_discrete else np.inf
        xp_assert_equal(res[endpoint_arg_minus & ~valid_arg], ref, **strict)
        xp_assert_equal(res[endpoint_arg_plus & ~valid_arg], ref, **strict)
    elif fname in {'pdf'}:
        xp_assert_equal(res[outside_arg_minus], 0, **strict)
        xp_assert_equal(res[outside_arg_plus], 0, **strict)
        ref = 0 if not is_discrete else np.inf
        xp_assert_equal(res[endpoint_arg_minus & ~valid_arg], ref, **strict)
        xp_assert_equal(res[endpoint_arg_plus & ~valid_arg], ref, **strict)
    elif fname in {'logcdf'}:
        xp_assert_equal(res[outside_arg_minus], -inf, **strict)
        xp_assert_equal(res[outside_arg_plus], 0, **strict)
        if not is_discrete:
            xp_assert_equal(res[endpoint_arg_minus], -inf, **strict)
        xp_assert_equal(res[endpoint_arg_plus], 0, **strict)
    elif fname in {'cdf'}:
        xp_assert_equal(res[outside_arg_minus], 0, **strict)
        xp_assert_equal(res[outside_arg_plus], 1, **strict)
        if not is_discrete:
            xp_assert_equal(res[endpoint_arg_minus], 0, **strict)
        xp_assert_equal(res[endpoint_arg_plus], 1, **strict)
    elif fname in {'logccdf'}:
        xp_assert_equal(res[outside_arg_minus], 0, **strict)
        xp_assert_equal(res[outside_arg_plus], -inf, **strict)
        if not is_discrete:
            xp_assert_equal(res[endpoint_arg_minus], 0, **strict)
        xp_assert_equal(res[endpoint_arg_plus], -inf, **strict)
    elif fname in {'ccdf'}:
        xp_assert_equal(res[outside_arg_minus], 1, **strict)
        xp_assert_equal(res[outside_arg_plus], 0, **strict)
        if not is_discrete:
            xp_assert_equal(res[endpoint_arg_minus], 1, **strict)
        xp_assert_equal(res[endpoint_arg_plus], 0, **strict)
    elif fname in {'ilogcdf', 'icdf'}:
        xp_assert_equal(res[outside_arg == -1], np.nan, **strict)
        xp_assert_equal(res[outside_arg == 1], np.nan, **strict)
        xp_assert_equal(res[endpoint_arg == -1], a[endpoint_arg == -1], **strict)
        xp_assert_equal(res[endpoint_arg == 1], b[endpoint_arg == 1], **strict)
    elif fname in {'ilogccdf', 'iccdf'}:
        xp_assert_equal(res[outside_arg == -1], np.nan, **strict)
        xp_assert_equal(res[outside_arg == 1], np.nan, **strict)
        xp_assert_equal(res[endpoint_arg == -1], b[endpoint_arg == -1], **strict)
        xp_assert_equal(res[endpoint_arg == 1], a[endpoint_arg == 1], **strict)

    exclude = {'mean', 'skewness', 'support'}
    if isinstance(dist, ContinuousDistribution):
        # Continuous distributions zero PMF everywhere (-inf in log-space)
        exclude.update({'logpmf'})
    if isinstance(dist, DiscreteDistribution):
        # Discrete distributions have infinite pdf at supported points
        exclude.update({'pdf', 'logpdf'})

    if fname not in exclude:
        dist.cache_policy = 'no_cache'
        variance = dist.variance()
        dist.cache_policy = None
        mask_finite = all_valid & (endpoint_arg == 0) & (variance > 0)
        assert xp.all(xp.isfinite(res[mask_finite]))


def check_periodicity(dist, fname, arg, res, a, b,
                      valid_parameters, classified_args, *, xp):
    inside_arg, endpoint_arg, outside_arg, nan_arg = classified_args
    # we've already checked NaN patterns; don't include known NaNs in these checks
    inside_arg = xp.where(valid_parameters, inside_arg, False)
    endpoint_arg = xp.where(valid_parameters, endpoint_arg, 0)
    method = getattr(dist, fname)
    period = (b - a)
    if fname in {'pdf', 'logpdf'}:
        xp_assert_close(method(arg + period), res)
        xp_assert_close(method(arg - period), res)
    elif fname in {'cdf'}:
        xp_assert_close(method(arg + period), res + 1)
        xp_assert_close(method(arg - period), res - 1)
        assert xp.all(res[inside_arg] <= 1.)
        assert xp.all(res[inside_arg] >= 0.)
        assert xp.all(res[endpoint_arg == 1] == 1.)
        assert xp.all(res[endpoint_arg == -1] == 0.)
    elif fname in {'icdf'}:
        xp_assert_close(method(arg + 1), res + period)
        xp_assert_close(method(arg - 1), res - period)
        assert xp.all(res[inside_arg] <= b[inside_arg])
        assert xp.all(res[inside_arg] >= a[inside_arg])
        assert xp.all(res[endpoint_arg == 1] == b[endpoint_arg == 1])
        assert xp.all(res[endpoint_arg == -1] == a[endpoint_arg == -1])
    elif fname in {'ccdf'}:
        xp_assert_close(method(arg + period), res - 1)
        xp_assert_close(method(arg - period), res + 1)
        assert xp.all(res[inside_arg] <= 1.)
        assert xp.all(res[inside_arg] >= 0.)
        assert xp.all(res[endpoint_arg == 1] == 0.)
        assert xp.all(res[endpoint_arg == -1] == 1.)
    elif fname in {'iccdf'}:
        xp_assert_close(method(arg + 1), res - period)
        xp_assert_close(method(arg - 1), res + period)
        assert xp.all(res[inside_arg] <= b[inside_arg])
        assert xp.all(res[inside_arg] >= a[inside_arg])
        assert xp.all(res[endpoint_arg == 1] == a[endpoint_arg == 1])
        assert xp.all(res[endpoint_arg == -1] == b[endpoint_arg == -1])


def check_moment_funcs(dist, result_shape, *, tol_override=None, xp):
    # Check that all computation methods of all distribution functions agree
    # with one another, effectively testing the correctness of the generic
    # computation methods and confirming the consistency of specific
    # distributions with their pdf/logpdf.
    default_atol = 1e-13 if dist._dtype == xp.float64 else 1e-7
    default_rtol = 1e-7 if dist._dtype == xp.float64 else 1e-5
    _tol_override = {'rtol': default_rtol, 'atol': default_atol}
    if tol_override is not None:
        _tol_override.update(tol_override)
    tol_override = _tol_override

    def check(order, kind, method=None, ref=None, success=True):
        if success:
            res = dist.moment(order, kind, method=method)
            rtol = tol_override['rtol']
            atol = tol_override['atol']
            xp_assert_close(res, ref, rtol=rtol, atol=atol*10**order)
            assert res.shape == ref.shape
        else:
            with pytest.raises(NotImplementedError):
                dist.moment(order, kind, method=method)

    def has_formula(order, kind):
        formula_name = f'_moment_{kind}_formula'
        overrides = dist._overrides(formula_name)
        if not overrides:
            return False
        formula = getattr(dist, formula_name)
        orders = getattr(formula, 'orders', set(range(6)))
        return order in orders

    dist.reset_cache()

    ### Check Raw Moments ###
    for i in range(6):
        check(i, 'raw', 'cache', success=False)  # not cached yet
        ref = dist.moment(i, 'raw', method='quadrature')
        check_nans_and_edges(dist, 'moment', None, ref, xp=xp)
        assert ref.shape == result_shape
        check(i, 'raw','cache', ref, success=True)  # cached now
        check(i, 'raw', 'formula', ref, success=has_formula(i, 'raw'))
        check(i, 'raw', 'general', ref, success=(i == 0))
        if dist.__class__ == stats.Normal:
            check(i, 'raw', 'quadrature_icdf', ref, success=True)


    # Clearing caches to better check their behavior
    dist.reset_cache()

    # If we have central or standard moment formulas, or if there are
    # values in their cache, we can use method='transform'
    dist.moment(0, 'central')  # build up the cache
    dist.moment(1, 'central')
    for i in range(2, 6):
        ref = dist.moment(i, 'raw', method='quadrature')
        check(i, 'raw', 'transform', ref,
              success=has_formula(i, 'central') or has_formula(i, 'standardized'))
        dist.moment(i, 'central')  # build up the cache
        check(i, 'raw', 'transform', ref)

    dist.reset_cache()

    ### Check Central Moments ###

    for i in range(6):
        check(i, 'central', 'cache', success=False)
        ref = dist.moment(i, 'central', method='quadrature')
        if i == 0:  # taking central moment of any order with quadrature requires mean
            dist._moment_central_cache.pop(1, None)
        assert ref.shape == result_shape
        check(i, 'central', 'cache', ref, success=True)
        check(i, 'central', 'formula', ref, success=has_formula(i, 'central'))
        check(i, 'central', 'general', ref,
              success=(i == 0) if isinstance(dist, CircularDistribution) else (i <= 1))
        if dist.__class__ == stats.Normal:
            check(i, 'central', 'quadrature_icdf', ref, success=True)
        if not (dist.__class__ == stats.Uniform and i == 5):
            # Quadrature is not super accurate for 5th central moment when the
            # support is really big. Skip this one failing test. We need to come
            # up with a better system of skipping individual failures w/ hypothesis.
            check(i, 'central', 'transform', ref,
                  success=has_formula(i, 'raw') or (i <= 1))
        if not has_formula(i, 'raw'):
            dist.moment(i, 'raw')
            check(i, 'central', 'transform', ref)

    # Circular distributions don't have kind='standardized'
    if isinstance(dist, CircularDistribution):
        return

    variance = dist.variance()
    dist.reset_cache()

    # If we have standard moment formulas, or if there are
    # values in their cache, we can use method='normalize'.
    if not xp.any(variance == 0):  # ...unless the variance is zero.
        # (In that case the standardized moment is NaN; we can't recover the central
        #  moment from that.)
        dist.moment(0, 'standardized')  # build up the cache
        dist.moment(1, 'standardized')
        dist.moment(2, 'standardized')
        for i in range(3, 6):
            ref = dist.moment(i, 'central', method='quadrature')
            check(i, 'central', 'normalize', ref,
                  success=has_formula(i, 'standardized'))
            dist.moment(i, 'standardized')  # build up the cache
            check(i, 'central', 'normalize', ref)

    ### Check Standardized Moments ###

    var = dist.moment(2, 'central', method='quadrature')
    dist.reset_cache()

    for i in range(6):
        check(i, 'standardized', 'cache', success=False)
        ref = dist.moment(i, 'central', method='quadrature') / var ** (i / 2)
        assert ref.shape == result_shape
        check(i, 'standardized', 'formula', ref,
              success=has_formula(i, 'standardized'))
        if not xp.any(variance == 0):
            # `method='general'` results just aren't correct for degenerate
            # distributions or when moments are not defined. That's OK.
            check(i, 'standardized', 'general', ref, success=i <= 2)
        check(i, 'standardized', 'normalize', ref)

    dist.reset_cache()


def check_lmoment_funcs(dist, result_shape, *, tol_override=None, xp):
    # Perform consistency check for L-moments similar to check_moment_funcs above
    default_atol = 5e-13 if dist._dtype == xp.float64 else 5e-7
    tol_override = {'atol' : default_atol} if tol_override is None else tol_override

    if not isinstance(dist, ContinuousDistribution):
        message = "L-moments are currently available only for continuous..."
        with pytest.raises(NotImplementedError, match=message):
            dist.lmoment(1)
        return

    def check(order, standardize=False, method=None, ref=None, success=True):
        if success:
            res = dist.lmoment(order, standardize=standardize, method=method)
            xp_assert_close(res, ref, **tol_override)
            assert res.shape == ref.shape
        else:
            with pytest.raises(NotImplementedError):
                dist.lmoment(order, standardize=standardize, method=method)

    ### Check L-Moments ###

    standardize = False
    for i in range(1, 6):
        check(i, standardize, 'cache', success=standardize)  # not cached yet
        ref = dist.lmoment(i, standardize=standardize, method='order_statistics')
        check_nans_and_edges(dist, 'lmoment', None, ref, xp=xp)
        assert ref.shape == result_shape
        check(i, standardize, 'cache', ref, success=True)  # cached now
        check(i, standardize, 'formula', ref,
              success=dist._overrides('_lmoment_formula')
                      and (i < 5 or dist.__class__.__name__ == "Uniform"))
        check(i, standardize, 'general', ref, success=(i == 1))
        if is_numpy(xp) and dist._overrides('_icdf_formula'):
            check(i, standardize, 'quadrature_icdf', ref, success=True)

    standardize=True
    for i in range(3, 6):
        ref = dist.lmoment(i, standardize=standardize, method='order_statistics')
        assert ref.shape == result_shape
        check(i, standardize, 'formula', ref,
              success=dist._overrides('_lmoment_formula')
                      and (i < 5 or dist.__class__.__name__ == "Uniform"))
        check(i, standardize, 'general', ref, success=False)
        if is_numpy(xp) and dist._overrides('_icdf_formula'):
            check(i, standardize, 'quadrature_icdf', ref, success=True)


# NumPy-only due to use of KS-test
@pytest.mark.parametrize('family', (Normal,))
@pytest.mark.parametrize('x_shape', [tuple(), (2, 3)])
@pytest.mark.parametrize('dist_shape', [tuple(), (4, 1)])
@pytest.mark.parametrize('fname', ['sample'])
@pytest.mark.parametrize('rng_type', [np.random.Generator, qmc.Halton, qmc.Sobol])
def test_sample_against_cdf(family, dist_shape, x_shape, fname, rng_type):
    rng = np.random.default_rng(842582438235635)
    num_parameters = family._num_parameters()

    if dist_shape and num_parameters == 0:
        pytest.skip("Distribution can't have a shape without parameters.")

    dist = family._draw(dist_shape, rng, xp=xp_compat_namespace(np))

    n = 1024
    sample_size = (n,) + x_shape
    sample_array_shape = sample_size + dist_shape

    if fname == 'sample':
        sample_method = dist.sample

    if rng_type != np.random.Generator:
        rng = rng_type(d=1, seed=rng)
    x = sample_method(sample_size, rng=rng)
    assert x.shape == sample_array_shape

    # probably should give `axis` argument to ks_1samp, review that separately
    statistic = _kolmogorov_smirnov(dist, x, axis=0)
    pvalue = kolmogn(x.shape[0], statistic, cdf=False)
    p_threshold = 0.01
    num_pvalues = pvalue.size
    num_small_pvalues = np.sum(pvalue < p_threshold)
    assert num_small_pvalues < p_threshold * num_pvalues


def get_valid_parameters(dist):
    # Given a distribution, return a logical array that is true where all
    # distribution parameters are within their respective domains. The code
    # here is probably quite similar to that used to form the `_invalid`
    # attribute of the distribution, but this was written about a week later
    # without referring to that code, so it is a somewhat independent check.

    # Get all parameter values and `_Parameter` objects
    xp = dist._xp
    parameter_values = dist._parameters
    parameters = {}
    for parameterization in dist._parameterizations:
        parameters.update(parameterization.parameters)

    all_valid = xp.ones(dist._shape, dtype=xp.bool)
    for name, value in parameter_values.items():
        if name not in parameters:  # cached value not part of parameterization
            continue
        parameter = parameters[name]

        # Check that the numerical endpoints and inclusivity attribute
        # agree with the `contains` method about which parameter values are
        # within the domain.
        a, b = parameter.domain.get_numerical_endpoints(
            parameter_values=parameter_values, xp=xp)
        a_included, b_included = parameter.domain.inclusive
        valid = (a <= value) if a_included else a < value
        valid &= (value <= b) if b_included else value < b
        xp_assert_equal(valid, parameter.domain.contains(
            value, parameter_values=parameter_values, xp=xp))

        # Form `all_valid` mask that is True where *all* parameters are valid
        all_valid &= valid

    # Check that the `all_valid` mask formed here is the complement of the
    # `dist._invalid` mask stored by the infrastructure
    xp_assert_equal(~all_valid, dist._invalid)

    return all_valid


def classify_arg(dist, arg, arg_domain, *, xp):
    if arg is None:
        valid_args = xp.ones(dist._shape, dtype=bool)
        endpoint_args = xp.zeros(dist._shape, dtype=bool)
        outside_args = xp.zeros(dist._shape, dtype=bool)
        nan_args = xp.zeros(dist._shape, dtype=bool)
        return valid_args, endpoint_args, outside_args, nan_args

    a, b = arg_domain.get_numerical_endpoints(
        parameter_values=dist._parameters, xp=xp)

    a, b, arg = xp.broadcast_arrays(a, b, arg)
    a_included, b_included = arg_domain.inclusive

    inside = (a <= arg) if a_included else a < arg
    inside &= (arg <= b) if b_included else arg < b
    # TODO: add `supported` method and check here
    on = xp.zeros(a.shape, dtype=xp.int64)
    on[a == arg] = -1
    on[b == arg] = 1
    outside = xp.zeros(a.shape, dtype=xp.int64)
    outside[(arg < a) if a_included else arg <= a] = -1
    outside[(b < arg) if b_included else b <= arg] = 1
    nan = xp.isnan(arg)

    return inside, on, outside, nan


@make_xp_test_case(stats.Normal)
def test_input_validation(xp):
    class Test(ContinuousDistribution):
        _variable = _RealParameter('x', domain=_RealInterval())
        def __init__(self, *args, **kwargs):
            return super().__init__(*args, xp=xp, **kwargs)

    message = ("The `Test` distribution family does not accept parameters, "
               "but parameters `{'a'}` were provided.")
    with pytest.raises(ValueError, match=message):
        Test(a=1, )

    message = "Attribute `tol` of `Test` must be a positive float, if specified."
    with pytest.raises(ValueError, match=message):
        Test(tol=xp.asarray([]))
    with pytest.raises(ValueError, match=message):
        Test(tol=xp.asarray([1, 2, 3]))
    with pytest.raises(ValueError, match=message):
        Test(tol=xp.nan)
    with pytest.raises(ValueError, match=message):
        Test(tol=-1)

    message = ("Argument `order` of `Test.moment` must be a "
               "finite integer greater than or equal to 0.")
    with pytest.raises(ValueError, match=message):
        Test().moment(-1)
    with pytest.raises(ValueError, match=message):
        Test().moment(xp.inf)

    message = ("Argument `order` of `Test.lmoment` must be a "
               "finite integer greater than or equal to 1.")
    with pytest.raises(ValueError, match=message):
        Test().lmoment(0)

    message = "Argument `kind` of `Test.moment` must be one of..."
    with pytest.raises(ValueError, match=message):
        Test().moment(2, kind='coconut')

    class Test2(ContinuousDistribution):
        _p1 = _RealParameter('c', domain=_RealInterval())
        _p2 = _RealParameter('d', domain=_RealInterval())
        _parameterizations = [_Parameterization(_p1, _p2)]
        _variable = _RealParameter('x', domain=_RealInterval())
        def __init__(self, *args, **kwargs):
            return super().__init__(*args, xp=xp, **kwargs)

    message = ("The provided parameters `{a}` do not match a supported "
               "parameterization of the `Test2` distribution family.")
    with pytest.raises(ValueError, match=message):
        Test2(a=1)

    message = ("The `Test2` distribution family requires parameters, but none "
               "were provided.")
    with pytest.raises(ValueError, match=message):
        Test2()

    message = ("The parameters `{c, d}` provided to the `Test2` "
               "distribution family cannot be broadcast to the same shape.")
    with pytest.raises(ValueError, match=message):
        Test2(c=[1, 2], d=[1, 2, 3])

    message = ("The argument provided to `Test2.pdf` cannot be broadcast to "
              "the same shape as the distribution parameters.")
    with pytest.raises(ValueError, match=message):
        dist = Test2(c=[1, 2, 3], d=[1, 2, 3])
        dist.pdf([1, 2])

    message = "Parameter `c` must be of real dtype."
    with pytest.raises(TypeError, match=message):
        Test2(c=[1, 1+0j], d=[1, 2])

    message = "`dtype` must be..."
    with pytest.raises(ValueError, match=message):
        Test2(c=1, d=1, dtype=xp.int64)

    message = "Parameter `convention` of `Test2.kurtosis` must be one of..."
    with pytest.raises(ValueError, match=message):
        dist = Test2(c=[1, 2, 3], d=[1, 2, 3])
        dist.kurtosis(convention='coconut')


@make_xp_test_case(stats.Uniform)
def test_rng_deepcopy_pickle(xp):
    # test behavior of `rng` attribute and copy behavior
    kwargs = dict(a=xp.asarray([-1, 2]), b=xp.asarray(10))
    dist1 = Uniform(**kwargs)
    dist2 = deepcopy(dist1)
    dist3 = pickle.loads(pickle.dumps(dist1))

    res1, res2, res3 = dist1.sample(), dist2.sample(), dist3.sample()
    assert xp.all(res2 != res1)
    assert xp.all(res3 != res1)

    res1, res2, res3 = dist1.sample(rng=42), dist2.sample(rng=42), dist3.sample(rng=42)
    xp_assert_equal(res2, res1)
    xp_assert_equal(res3, res1)


@make_xp_test_case(stats.Normal)
class TestAttributes:
    def test_cache_policy(self, xp):
        dist = StandardNormal(cache_policy="no_cache", xp=xp)
        # make error message more appropriate
        message = "`StandardNormal` does not provide an accurate implementation of the "
        with pytest.raises(NotImplementedError, match=message):
            dist.mean(method='cache')
        mean = dist.mean()
        with pytest.raises(NotImplementedError, match=message):
            dist.mean(method='cache')

        # add to enum
        dist.cache_policy = None
        with pytest.raises(NotImplementedError, match=message):
            dist.mean(method='cache')
        mean = dist.mean()  # method is 'formula' by default
        cached_mean = dist.mean(method='cache')
        xp_assert_equal(cached_mean, mean)

        # cache is not cleared when used in a TransformedDistribution
        dist2 = stats.exp(dist)
        xp_assert_equal(dist.mean(method='cache'), cached_mean)
        dist2.reset_cache()
        xp_assert_equal(dist.mean(method='cache'), cached_mean)

        # cache is overridden by latest evaluation
        quadrature_mean = dist.mean(method='quadrature')
        cached_mean = dist.mean(method='cache')
        xp_assert_equal(cached_mean, quadrature_mean)
        assert not xp.all(mean == quadrature_mean)

        # We can turn the cache off, and it won't change, but the old cache is
        # still available
        dist.cache_policy = "no_cache"
        mean = dist.mean(method='formula')
        cached_mean = dist.mean(method='cache')
        xp_assert_equal(cached_mean, quadrature_mean)
        assert not xp.all(mean == quadrature_mean)

        dist.reset_cache()
        with pytest.raises(NotImplementedError, match=message):
            dist.mean(method='cache')

        message = "Attribute `cache_policy` of `StandardNormal`..."
        with pytest.raises(ValueError, match=message):
            dist.cache_policy = "invalid"

    def test_tol(self, xp):
        x = 3.
        X = stats.Normal(xp=xp)

        message = "Attribute `tol` of `StandardNormal` must..."
        with pytest.raises(ValueError, match=message):
            X.tol = -1.
        with pytest.raises(ValueError, match=message):
            X.tol = (0.1,)
        with pytest.raises(ValueError, match=message):
            X.tol = xp.nan

        fine_tol = 1e-12 if xpx.default_dtype(xp) == xp.float64 else 1e-6
        X1 = stats.Normal(tol=1e-1, xp=xp)
        X2 = stats.Normal(tol=fine_tol, xp=xp)
        ref = X.cdf(x)
        res1 = X1.cdf(x, method='quadrature')
        res2 = X2.cdf(x, method='quadrature')
        xp_assert_close(res1, ref, rtol=X1.tol)
        xp_assert_close(res2, ref, rtol=X2.tol)
        xp_assert_less(abs(res2 - ref), abs(res1 - ref))

        p = 0.99
        X1.tol, X2.tol = X2.tol, X1.tol
        ref = X.icdf(p)
        res1 = X1.icdf(p, method='inversion')
        res2 = X2.icdf(p, method='inversion')
        xp_assert_close(res1, ref, rtol=X1.tol)
        xp_assert_close(res2, ref, rtol=X2.tol)
        xp_assert_less(abs(res1 - ref), abs(res2 - ref))

    @pytest.mark.skip_xp_backends("torch", reason="can't make string array")
    def test_iv_policy(self, xp):
        X = Uniform(a=0, b=1, xp=xp)
        assert X.pdf(2) == 0

        X.validation_policy = 'skip_all'
        assert X.pdf(xp.asarray(2.)) == 1

        # Tests _set_invalid_nan
        a, b = xp.asarray(1.), xp.asarray(0.)  # invalid parameters
        X = Uniform(a=a, b=b, validation_policy='skip_all', xp=xp)
        assert X.pdf(xp.asarray(2.)) == -1

        # Tests _set_invalid_nan_property
        class MyUniform(Uniform):
            def _entropy_formula(self, *args, **kwargs):
                return 'incorrect'

            def _moment_raw_formula(self, order, **params):
                return 'incorrect'

        X = MyUniform(a=a, b=b, validation_policy='skip_all', xp=xp)
        assert X.entropy() == 'incorrect'

        # Tests _validate_order_kind
        assert X.moment(kind='raw', order=-1) == 'incorrect'

        # Test input validation
        message = "Attribute `validation_policy` of `MyUniform`..."
        with pytest.raises(ValueError, match=message):
            X.validation_policy = "invalid"

    def test_shapes(self, xp):
        X = stats.Normal(mu=1, sigma=2, xp=xp)
        Y = stats.Normal(mu=[2], sigma=3, xp=xp)

        # Check that attributes are available as expected
        assert X.mu == 1
        assert X.sigma == 2
        assert Y.mu[0] == 2
        assert Y.sigma[0] == 3

        # Trying to set an attribute raises
        # message depends on Python version
        with pytest.raises(AttributeError):
            X.mu = 2

        # Trying to mutate an attribute really mutates a copy
        Y.mu[0] = 10
        assert Y.mu[0] == 2


@make_xp_test_case(TransformedDistribution)
class TestTransforms:
    def test_ContinuousDistribution_only(self, xp):
        X = stats.Binomial(n=10, p=0.5, xp=xp)
        # This is applied at the top level TransformedDistribution,
        # so testing one subclass is enough
        message = "Transformations are currently supported only for..."
        with pytest.raises(NotImplementedError, match=message):
            stats.exp(X)


@make_xp_test_case(ShiftedScaledDistribution)
class TestLocScale:
    @pytest.mark.fail_slow(10)
    @pytest.mark.parametrize('seed', 911679639 + np.arange(20))
    def test_loc_scale(self, seed, xp):
        # Need tests with negative scale
        rng = np.random.default_rng(seed)

        class TransformedNormal(ShiftedScaledDistribution):
            def __init__(self, *args, **kwargs):
                super().__init__(StandardNormal(xp=xp), *args, **kwargs)

        tmp = draw_distribution_from_family(
            TransformedNormal, rng, {'min_side': 1}, proportions=(1, 0, 0, 0), xp=xp)
        dist, x, y, p, logp = tmp.dist, tmp.x, tmp.y, tmp.p, tmp.logp

        loc = dist.loc
        # negative scale tested in test_abs_finite_support, test_reciprocal, etc.
        scale = dist.scale
        dist0 = StandardNormal(xp=xp)
        dist_ref = Normal(mu=loc, sigma=scale, xp=xp)
        z = xp.zeros_like(scale)  # compact broadcasting

        x0 = (x - loc) / scale
        y0 = (y - loc) / scale

        a, b = dist.support()
        a0, b0 = dist0.support()
        xp_assert_close(a, a0 + loc)
        xp_assert_close(b, b0 + loc)

        with np.errstate(invalid='ignore', divide='ignore'):
            xp_assert_close(xp.exp(dist.logentropy()), dist.entropy() + 0j)
            xp_assert_close(dist.entropy(), dist_ref.entropy())
            xp_assert_close(dist.median(), dist0.median() + loc)
            xp_assert_close(dist.mode(), dist0.mode() + loc)
            xp_assert_close(dist.mean(), dist0.mean() + loc)
            xp_assert_close(dist.variance(), dist0.variance() * scale**2)
            xp_assert_close(dist.standard_deviation(), dist.variance()**0.5)
            xp_assert_close(dist.skewness(), dist0.skewness() * xp.sign(scale))
            xp_assert_close(dist.kurtosis(), dist0.kurtosis() + z)
            xp_assert_close(dist.logpdf(x), dist0.logpdf(x0) - xp.log(scale))
            xp_assert_close(dist.pdf(x), dist0.pdf(x0) / scale)
            xp_assert_close(dist.logcdf(x), dist0.logcdf(x0))
            xp_assert_close(dist.cdf(x), dist0.cdf(x0))
            xp_assert_close(dist.logccdf(x), dist0.logccdf(x0))
            xp_assert_close(dist.ccdf(x), dist0.ccdf(x0))
            xp_assert_close(dist.logcdf(x, y), dist0.logcdf(x0, y0))
            xp_assert_close(dist.cdf(x, y), dist0.cdf(x0, y0))
            xp_assert_close(dist.logccdf(x, y), dist0.logccdf(x0, y0))
            xp_assert_close(dist.ccdf(x, y), dist0.ccdf(x0, y0))
            xp_assert_close(dist.ilogcdf(logp), dist0.ilogcdf(logp)*scale + loc)
            xp_assert_close(dist.icdf(p), dist0.icdf(p)*scale + loc)
            xp_assert_close(dist.ilogccdf(logp), dist0.ilogccdf(logp)*scale + loc)
            xp_assert_close(dist.iccdf(p), dist0.iccdf(p)*scale + loc)
            for i in range(1, 5):
                xp_assert_close(dist.moment(i, 'raw'), dist_ref.moment(i))
                xp_assert_close(dist.moment(i, 'central'),
                                dist0.moment(i, 'central') * scale**i)
                xp_assert_close(dist.moment(i, 'standardized'),
                                dist0.moment(i, 'standardized') * np.sign(scale)**i)
            for i in range(1, 5):
                xp_assert_close(dist.lmoment(i), dist_ref.lmoment(i), atol=1e-8)
                if i >= 3:
                    xp_assert_close(dist.lmoment(i, standardize=True),
                                    dist_ref.lmoment(i, standardize=True),
                                    atol=1e-8)

        # Transform back to the original distribution using all arithmetic
        # operations; check that it behaves as expected.
        dist = (dist - 2*loc) + loc
        dist = dist/scale**2 * scale

        a, b = dist.support()
        a0, b0 = dist0.support()
        xp_assert_close(a, a0 + z)
        xp_assert_close(b, b0 + z)

        with np.errstate(invalid='ignore', divide='ignore'):
            xp_assert_close(dist.logentropy(), dist0.logentropy() + z)
            xp_assert_close(dist.entropy(), dist0.entropy() + z)
            xp_assert_close(dist.median(), dist0.median() + z)
            xp_assert_close(dist.mode(), dist0.mode() + z)
            xp_assert_close(dist.mean(), dist0.mean() + z)
            xp_assert_close(dist.variance(), dist0.variance() + z)
            xp_assert_close(dist.standard_deviation(), dist0.standard_deviation() + z)
            xp_assert_close(dist.skewness(), dist0.skewness() + z)
            xp_assert_close(dist.kurtosis(), dist0.kurtosis() + z)
            xp_assert_close(dist.logpdf(x), dist0.logpdf(x)+z)
            xp_assert_close(dist.pdf(x), dist0.pdf(x) + z)
            xp_assert_close(dist.logcdf(x), dist0.logcdf(x) + z)
            xp_assert_close(dist.cdf(x), dist0.cdf(x) + z)
            xp_assert_close(dist.logccdf(x), dist0.logccdf(x) + z)
            xp_assert_close(dist.ccdf(x), dist0.ccdf(x) + z)
            xp_assert_close(dist.ilogcdf(logp), dist0.ilogcdf(logp) + z)
            xp_assert_close(dist.icdf(p), dist0.icdf(p) + z)
            xp_assert_close(dist.ilogccdf(logp), dist0.ilogccdf(logp) + z)
            xp_assert_close(dist.iccdf(p), dist0.iccdf(p) + z)
            for i in range(1, 5):
                xp_assert_close(dist.moment(i, 'raw'),
                                dist0.moment(i, 'raw') + z)
                xp_assert_close(dist.moment(i, 'central'),
                                dist0.moment(i, 'central') + z)
                xp_assert_close(dist.moment(i, 'standardized'),
                                dist0.moment(i, 'standardized') + z)
            for i in range(1, 5):
                xp_assert_close(dist.lmoment(i), dist0.lmoment(i) + z)
                if i >= 3:
                    xp_assert_close(dist.lmoment(i, standardize=True),
                                    dist0.lmoment(i, standardize=True) + z)

            # These are tough to compare because of the way the shape works
            # rng = np.random.default_rng(seed)
            # rng0 = np.random.default_rng(seed)
            # assert_allclose(dist.sample(x_result_shape, rng=rng),
            #                 dist0.sample(x_result_shape, rng=rng0) * scale + loc)
            # Should also try to test fit, plot?

@make_xp_test_case(TransformedDistribution)
class TestReciprocal:
    @pytest.mark.fail_slow(10)
    @pytest.mark.parametrize('scale', [1, 2, -1])
    @pytest.mark.xfail_on_32bit("`scale=-1` fails on 32-bit; needs investigation")
    def test_reciprocal(self, scale, xp):
        if scale == -1 and (is_cupy(xp) or is_torch(xp)):
            pytest.skip("`scale=-1' fails w/ CuPy; needs investigation.")
        rng = np.random.default_rng(81345982345826)
        a = rng.random((3, 1))

        # Separate sign from scale. It's easy to scale the resulting
        # RV with negative scale; we want to test the ability to divide
        # by a RV with negative support
        sign, scale = np.sign(scale), float(abs(scale))

        # Reference distribution
        InvGamma = stats.make_distribution(stats.invgamma)
        Y0 = sign * scale * InvGamma(a=a)

        # Test distribution
        a = xp.asarray(a)
        X = _Gamma(a=a) if sign > 0 else -_Gamma(a=a)
        # shouldn't matter whether scale is a scalar or array...
        Y = (scale if is_cupy(xp) else xp.asarray(scale)) / X

        y = Y0.sample(shape=(3, 10), rng=rng)
        p = Y0.cdf(y)
        logp = np.log(p)
        y_xp, p_xp, logp_xp = xp.asarray(y), xp.asarray(p), xp.asarray(logp)

        xp_assert_close(Y.logentropy(), xp.asarray(np.log(Y0.entropy()+0j)))
        xp_assert_close(Y.entropy(), xp.asarray(Y0.entropy()))
        xp_assert_close(Y.median(), xp.asarray(Y0.median()))
        # moments are not finite
        xp_assert_close(Y.support()[0], xp.asarray(Y0.support()[0]))
        xp_assert_close(Y.support()[1], xp.asarray(Y0.support()[1]))
        xp_assert_close(Y.pdf(y_xp), xp.asarray(Y0.pdf(y)))
        xp_assert_close(Y.cdf(y_xp), xp.asarray(Y0.cdf(y)))
        xp_assert_close(Y.ccdf(y_xp), xp.asarray(Y0.ccdf(y)))
        xp_assert_close(Y.icdf(p_xp), xp.asarray(Y0.icdf(p)), rtol=1e-7)
        xp_assert_close(Y.iccdf(p_xp), xp.asarray(Y0.iccdf(p)))
        xp_assert_close(Y.logpdf(y_xp), xp.asarray(Y0.logpdf(y)))
        xp_assert_close(Y.logcdf(y_xp), xp.asarray(Y0.logcdf(y)))
        xp_assert_close(Y.logccdf(y_xp), xp.asarray(Y0.logccdf(y)))
        with np.errstate(divide='ignore', invalid='ignore'):
            xp_assert_close(Y.ilogcdf(logp_xp), xp.asarray(Y0.ilogcdf(logp)), rtol=1e-7)
            xp_assert_close(Y.ilogccdf(logp_xp), xp.asarray(Y0.ilogccdf(logp)))
        seed = 3984593485
        xp_assert_close(Y.sample(rng=seed), scale/(X.sample(rng=seed)))


@make_xp_test_case(MonotonicTransformedDistribution)
class TestMonotonicTransform:
    def test_monotonic_transforms_str(self):
        X = Uniform(a=1, b=2)
        X_str = "Uniform(a=1, b=2)"

        assert str(stats.log(X)) == f"log({X_str})"
        assert str(1 / X) == f"1/({X_str})"
        assert str(stats.exp(X)) == f"exp({X_str})"

    def test_monotonic_transforms_input_validation(self, xp):
        X = Uniform(a=-1, b=2, xp=xp)
        message = "Division by a random variable is only implemented when the..."
        with pytest.raises(NotImplementedError, match=message):
            1 / X
        message = "The logarithm of a random variable is only implemented when the..."
        with pytest.raises(NotImplementedError, match=message):
            stats.log(X)
        message = "Raising an argument to the power of a random variable is only..."
        with pytest.raises(NotImplementedError, match=message):
            (-2) ** X
        with pytest.raises(NotImplementedError, match=message):
            1 ** X
        with pytest.raises(NotImplementedError, match=message):
            [0.5, 1.5] ** X

        message = "Raising a random variable to the power of an argument is only"
        with pytest.raises(NotImplementedError, match=message):
            X ** (-2)
        with pytest.raises(NotImplementedError, match=message):
            X ** 0
        with pytest.raises(NotImplementedError, match=message):
            X ** [0.5, 1.5]


@make_xp_test_case(TransformedDistribution)
class TestArithmetic:
    def test_arithmetic_operators(self, xp):
        rng = np.random.default_rng(2348923495832349834)

        a, b, loc, scale = 0.294, 1.34, 0.57, 1.16

        Y = _LogUniform(a=a, b=b, xp=xp)
        x = xp.asarray(rng.uniform(-3, 3, 100), dtype=Y._dtype)

        X = scale*Y + loc
        xp_assert_close(X.cdf(x), Y.cdf((x - loc) / scale))
        X = loc + Y*scale
        xp_assert_close(X.cdf(x), Y.cdf((x - loc) / scale))

        X = Y/scale - loc
        xp_assert_close(X.cdf(x), Y.cdf((x + loc) * scale))
        X = loc - Y/scale
        xp_assert_close(X.cdf(x), Y.ccdf((-x + loc)*scale))


@make_xp_test_case(stats.abs)
class TestAbs:
    def test_abs(self, xp):
        rng = np.random.default_rng(81345982345826)
        loc = rng.random((3, 1))

        Y = stats.abs(Normal(xp=xp) + xp.asarray(loc))
        Y0 = stats.foldnorm(loc)

        y = Y0.rvs((3, 10), random_state=rng)
        p = Y0.cdf(y)
        y_xp, p_xp = xp.asarray(y), xp.asarray(p)

        xp_assert_close(Y.logentropy(), xp.asarray(np.log(Y0.entropy() + 0j)))
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
        xp_assert_close(Y.ilogcdf(xp.log(p_xp)), xp.asarray(Y0.ppf(p)))
        xp_assert_close(Y.ilogccdf(xp.log(p_xp)), xp.asarray(Y0.isf(p)))
        sample = Y.sample(10)
        assert xp.all(sample > 0)

    def test_abs_finite_support(self):
        # The original implementation of `FoldedDistribution` might evaluate
        # the private distribution methods outside the support. Check that this
        # is resolved.
        Weibull = stats.make_distribution(stats.weibull_min)
        X = Weibull(c=2)
        Y = abs(-X)
        assert_equal(X.logpdf(1), Y.logpdf(1))
        assert_equal(X.pdf(1), Y.pdf(1))
        assert_equal(X.logcdf(1), Y.logcdf(1))
        assert_equal(X.cdf(1), Y.cdf(1))
        assert_equal(X.logccdf(1), Y.logccdf(1))
        assert_equal(X.ccdf(1), Y.ccdf(1))


@make_xp_test_case(TransformedDistribution)
class TestPow:
    def test_pow(self, xp):
        rng = np.random.default_rng(81345982345826)

        dtype = xp.float64  # stats.chi2 doesn't do float32
        Y = Normal(xp=xp, dtype=dtype)**2
        Y0 = stats.chi2(df=1)

        y = Y0.rvs(10, random_state=rng)
        p = Y0.cdf(y)
        y_xp, p_xp = xp.asarray(y), xp.asarray(p)

        xp_assert_close(Y.logentropy(),
                        xp.asarray(np.log(Y0.entropy() + 0j)),
                        rtol=1e-6)
        xp_assert_close(Y.entropy(), xp.asarray(Y0.entropy()), rtol=1e-6)
        xp_assert_close(Y.median(), xp.asarray(Y0.median()))
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
        xp_assert_close(Y.ilogcdf(xp.log(p_xp)), xp.asarray(Y0.ppf(p)))
        xp_assert_close(Y.ilogccdf(xp.log(p_xp)), xp.asarray(Y0.isf(p)))
        sample = Y.sample(10)
        assert xp.all(sample > 0)


@make_xp_test_case(stats.make_distribution)
class TestQuantileDefinedDistribution:
    @pytest.mark.slow
    @pytest.mark.parametrize('fun', ['cdf', 'ccdf', 'icdf', 'iccdf'])
    def test_tukey_lambda(self, fun, xp):
        if fun in {'cdf', 'ccdf'} and is_cupy(xp):
            pytest.skip("special.tklmbda doesn't work with CuPy arrays")

        class MyTukeyLambda:
            __make_distribution_version__ = "1.16.0"

            @property
            def parameters(self):
                return {'lam': {'endpoints': (-xp.inf, xp.inf),
                                'inclusive': (False, False)}}

            @property
            def support(self):
                return {'endpoints': (lambda *, lam: xp.where(lam > 0, -1/lam, -xp.inf),
                                      lambda *, lam: xp.where(lam > 0, 1/lam, xp.inf)),
                        'inclusive': (False, False)}

        def cdf(self, x, lam):
            return special.tklmbda(x, lam)

        def ccdf(self, x, lam):
            return 1 - special.tklmbda(x, lam)

        def icdf(self, x, lam):
            with np.errstate(divide='ignore'):
                return 1/lam * (x**lam - (1 - x)**lam)

        def iccdf(self, x, lam):
            with np.errstate(divide='ignore'):
                return 1/lam * ((1-x)**lam - x**lam)

        funs = {'cdf': cdf, 'ccdf': ccdf, 'icdf': icdf, 'iccdf': iccdf}
        setattr(MyTukeyLambda, fun, funs[fun])

        TukeyLambda1 = stats.make_distribution(MyTukeyLambda())
        TukeyLambda2 = stats.make_distribution(stats.tukeylambda)

        xp_dtype = xp.float64  # stats.tukeylambda only does float64
        lam = 5
        X = TukeyLambda1(lam=lam, dtype=xp_dtype, xp=xp)
        Y = TukeyLambda2(lam=lam)

        x = np.linspace(-1/lam, 1/lam, 10)
        p = Y.cdf(x)
        x_xp, p_xp = xp.asarray(x, dtype=X._dtype), xp.asarray(p, dtype=X._dtype)

        with np.errstate(divide='ignore'):  # stats.tukeylambda is noisy
            xp_assert_close(X.support()[0], xp.asarray(Y.support()[0]))
            xp_assert_close(X.support()[1], xp.asarray(Y.support()[1]))
            xp_assert_close(X.logentropy(), xp.asarray(Y.logentropy()+0j))
            xp_assert_close(X.entropy(), xp.asarray(Y.entropy()))
            xp_assert_close(X.mode(), xp.asarray(Y.mode()), atol=5e-7)
            xp_assert_close(X.median(), xp.asarray(Y.median()), atol=1e-10)
            xp_assert_close(X.mean(), xp.asarray(Y.mean()), atol=1e-10)
            xp_assert_close(X.variance(), xp.asarray(Y.variance()))
            xp_assert_close(X.standard_deviation(), xp.asarray(Y.standard_deviation()))
            xp_assert_close(X.skewness(), xp.asarray(Y.skewness()), atol=1e-10)
            xp_assert_close(X.kurtosis(), xp.asarray(Y.kurtosis()))
            xp_assert_close(X.logpdf(x_xp), xp.asarray(Y.logpdf(x)))
            xp_assert_close(X.pdf(x_xp), xp.asarray(Y.pdf(x)))
            xp_assert_close(X.logcdf(x_xp), xp.asarray(Y.logcdf(x)))
            xp_assert_close(X.cdf(x_xp), xp.asarray(Y.cdf(x)))
            xp_assert_close(X.logccdf(x_xp), xp.asarray(Y.logccdf(x)))
            xp_assert_close(X.ccdf(x_xp), xp.asarray(Y.ccdf(x)))
            xp_assert_close(X.ilogcdf(p_xp), xp.asarray(Y.ilogcdf(p)))
            xp_assert_close(X.icdf(p_xp), xp.asarray(Y.icdf(p)))
            xp_assert_close(X.ilogccdf(p_xp), xp.asarray(Y.ilogccdf(p)))
            xp_assert_close(X.iccdf(p_xp), xp.asarray(Y.iccdf(p)))
            for kind in ['raw', 'central', 'standardized']:
                for order in range(5):
                    xp_assert_close(X.moment(order, kind=kind),
                                    xp.asarray(Y.moment(order, kind=kind)),
                                    atol=1e-8)


# Adds tests just to get to 100% test coverage; this way it's more obvious
# if new lines are untested.
class TestFullCoverage:
    def test_ContinuousDistribution_set_invalid_nan(self):
        # Exercise code paths when formula returns wrong shape and dtype
        # We could consider making this raise an error to force authors
        # to return the right shape and dytpe, but this would need to be
        # configurable.
        class TestDist(ContinuousDistribution):
            _variable = _RealParameter('x', domain=_RealInterval(endpoints=(0., 1.)))
            def _logpdf_formula(self, x, *args, **kwargs):
                return 0

        X = TestDist()
        dtype = np.float32
        X._dtype = dtype
        x = np.asarray([0.5], dtype=dtype)
        assert X.logpdf(x).dtype == dtype

    def test_generate_domain_support(self):
        msg = _generate_domain_support(StandardNormal)
        assert "accepts no distribution parameters" in msg

        msg = _generate_domain_support(Normal)
        assert "accepts one parameterization" in msg

        msg = _generate_domain_support(_LogUniform)
        assert "accepts two parameterizations" in msg

    def test_ContinuousDistribution__repr__(self):
        X = Uniform(a=0, b=1)
        assert repr(X) == "Uniform(a=0, b=1)"
        assert repr(X*3 + 2) == "np.float64(3.0)*Uniform(a=0, b=1) + np.float64(2.0)"

        X = Uniform(a=np.zeros(4), b=1)
        assert repr(X) == "Uniform(a=array([0., 0., 0., 0.]), b=1)"

        X = Uniform(a=np.zeros(4, dtype=np.float32), b=np.ones(4, dtype=np.float32))
        assert repr(X) == (
            "Uniform(a=array([0., 0., 0., 0.], dtype=float32),"
            " b=array([1., 1., 1., 1.], dtype=float32))"
        )


class TestReprs:
    U = Uniform(a=0, b=1)
    V = Uniform(a=np.float32(0.0), b=np.float32(1.0))
    X = Normal(mu=-1, sigma=1)
    Y = Normal(mu=1, sigma=1)
    Z = Normal(mu=np.zeros(1000), sigma=1)

    @pytest.mark.parametrize(
        "dist",
        [
            U,
            U - np.array([1.0, 2.0]),
            V,
            np.ones(2, dtype=np.float32)*V + np.zeros(2, dtype=np.float64),
            3*U + 2,
            U**4,
            (3*U + 2)**4,
            (3*U + 2)**3,
            2**U,
            2**(3*U + 1),
            1 / (1 + U),
            stats.order_statistic(U, r=3, n=5),
            stats.truncate(U, 0.2, 0.8),
            stats.Mixture([X, Y], weights=[0.3, 0.7]),
            abs(U),
            stats.exp(U),
            stats.log(1 + U),
            np.array([1.0, 2.0])*U + np.array([2.0, 3.0]),
        ]
    )
    def test_executable(self, dist):
        # Test that reprs actually evaluate to proper distribution
        # provided relevant imports are made.
        from numpy import array  # noqa: F401
        from numpy import float32  # noqa: F401
        from scipy.stats import abs, exp, log, order_statistic, truncate # noqa: F401
        from scipy.stats import Mixture, Normal # noqa: F401
        from scipy.stats._new_distributions import Uniform # noqa: F401
        new_dist = eval(repr(dist))
        # A basic check that the distributions are the same
        sample1 = dist.sample(shape=10, rng=1234)
        sample2 = new_dist.sample(shape=10, rng=1234)
        assert_equal(sample1, sample2)
        assert sample1.dtype is sample2.dtype

    @pytest.mark.parametrize(
        "dist",
        [
            Z,
            np.full(1000, 2.0) * X + 1.0,
            2.0 * X + np.full(1000, 1.0),
            np.full(1000, 2.0) * X + 1.0,
            stats.truncate(Z, -1, 1),
            stats.truncate(Z, -np.ones(1000), np.ones(1000)),
            stats.order_statistic(X, r=np.arange(1, 1000), n=1000),
            Z**2,
            1.0 / (1 + stats.exp(Z)),
            2**Z,
        ]
    )
    def test_not_too_long(self, dist):
        # Tests that array summarization is working to ensure reprs aren't too long.
        # None of the reprs above will be executable.
        assert len(repr(dist)) < 250


@make_xp_test_case(stats.VonMises)
class TestCircular:
    def test_input_validation(self, xp):
        X = stats.VonMises(mu=0, kappa=1, xp=xp)

        message = "Circular distributions do not support "
        with pytest.raises(NotImplementedError, match=message + '`lmoment`.'):
            X.lmoment()
        with pytest.raises(NotImplementedError, match=message + '`logcdf`.'):
            X.logcdf(1)
        with pytest.raises(NotImplementedError, match=message + '`logcdf`.'):
            X.logcdf(1, 2)
        with pytest.raises(NotImplementedError, match=message + '`logccdf`.'):
            X.logccdf(1)
        with pytest.raises(NotImplementedError, match=message + '`logccdf`.'):
            X.logccdf(1, 2)
        with pytest.raises(NotImplementedError, match=message + '`ilogcdf`.'):
            X.ilogcdf(1)
        with pytest.raises(NotImplementedError, match=message + '`ilogccdf`.'):
            X.ilogccdf(1)

        message = "Circular distributions do not support two-argument "
        with pytest.raises(NotImplementedError, match=message + '`cdf`.'):
            X.cdf(1, 2)
        with pytest.raises(NotImplementedError, match=message + '`ccdf`.'):
            X.ccdf(1, 2)

        message = "`VonMises.kurtosis` supports only the default value of `convention`."
        with pytest.raises(ValueError, match=message):
            X.kurtosis(convention='excess')

        message = "Argument `kind` of `VonMises.moment` must be one of..."
        with pytest.raises(ValueError, match=message):
            X.moment(3, kind='standardized')

        message = "`VonMises` does not provide an accurate implementation of..."
        with pytest.raises(NotImplementedError, match=message):
            X.moment(3, kind='central', method='normalize')

    def test_mixture(self, xp):
        X = stats.VonMises(mu=0, kappa=1, xp=xp)
        message = "`Mixture` does not currently support circular components."
        with pytest.raises(NotImplementedError, match=message):
            stats.Mixture([X, X])

    @pytest.mark.parametrize("op, args", [
        (stats.abs, ()), (stats.exp, ()), (stats.log, ()), (operator.neg, ()),
        (operator.add, (1,)), (operator.sub, (1,)), (operator.pow, (1,)),
        (operator.mul, (1,)), (operator.truediv, (1,)), (stats.truncate, (1, 2))])
    def test_transform(self, op, args, xp):
        X = stats.VonMises(mu=0, kappa=1, xp=xp)
        message = "Transformations are currently supported only for continuous..."
        with pytest.raises(NotImplementedError, match=message):
            op(X, *args)
