import itertools as it
import math
import operator
import os
import pickle
from copy import deepcopy

import numpy as np
from numpy import inf
import pytest
from numpy.testing import assert_allclose, assert_equal
from hypothesis import strategies, given
import hypothesis.extra.numpy as npst

from scipy import special
from scipy import stats
from scipy.stats._fit import _kolmogorov_smirnov
from scipy.stats._ksstats import kolmogn
from scipy.stats import qmc
from scipy.stats._distr_params import distcont, distdiscrete
from scipy.stats._distribution_infrastructure import (
    _Domain, _RealInterval, _Parameter, _Parameterization, _RealParameter,
    ContinuousDistribution, DiscreteDistribution, CircularDistribution, Mixture,
    ShiftedScaledDistribution, _generate_domain_support, _logexpxmexpy)
from scipy.stats._new_distributions import (
    StandardNormal, _LogUniform, _Gamma, _StandardCauchy, _TestCircular)
from scipy.stats import Normal, Uniform
from scipy._lib._testutils import mutually_broadcastable_shapes
from scipy._lib._util import _RichResult
from scipy._lib._array_api import (
    make_xp_test_case, is_numpy, is_cupy, is_torch, xp_copy)
from scipy._lib._array_api_no_0d import xp_assert_close, xp_assert_equal, xp_assert_less


import torch
torch.set_default_dtype(torch.float32)

class Test_RealInterval:
    rng = np.random.default_rng(349849812549824)

    def test_iv(self):
        domain = _RealInterval(endpoints=('a', 'b'))
        message = "The endpoints of the distribution are defined..."
        with pytest.raises(TypeError, match=message):
            domain.get_numerical_endpoints(dict, xp=np)

    @pytest.mark.parametrize('x', [rng.uniform(10, 10, size=(2, 3, 4)),
                                   -np.inf, np.pi])
    def test_contains_simple(self, x):
        # Test `contains` when endpoints are defined by constants
        a, b = -np.inf, np.pi
        domain = _RealInterval(endpoints=(a, b), inclusive=(False, True))
        assert_equal(domain.contains(x, xp=np), (a < x) & (x <= b))

    @pytest.mark.slow
    @given(shapes=npst.mutually_broadcastable_shapes(num_shapes=3, min_side=0),
           inclusive_a=strategies.booleans(),
           inclusive_b=strategies.booleans(),
           data=strategies.data())
    def test_contains(self, shapes, inclusive_a, inclusive_b, data):
        # Test `contains` when endpoints are defined by parameters
        shapes = [tuple(shape) for shape in shapes]
        input_shapes, result_shape = shapes
        shape_a, shape_b, shape_x = input_shapes

        # Without defining min and max values, I spent forever trying to set
        # up a valid test without overflows or similar just drawing arrays.
        a_elements = dict(allow_nan=False, allow_infinity=False,
                          min_value=-1e3, max_value=1)
        b_elements = dict(allow_nan=False, allow_infinity=False,
                          min_value=2, max_value=1e3)
        a = data.draw(npst.arrays(npst.floating_dtypes(),
                                  shape_a, elements=a_elements))
        b = data.draw(npst.arrays(npst.floating_dtypes(),
                                  shape_b, elements=b_elements))
        # ensure some points are to the left, some to the right, and some
        # are exactly on the boundary
        d = b - a
        x = np.concatenate([np.linspace(a-d, a, 10),
                            np.linspace(a, b, 10),
                            np.linspace(b, b+d, 10)])
        # Domain is defined by two parameters, 'a' and 'b'
        domain = _RealInterval(endpoints=('a', 'b'),
                             inclusive=(inclusive_a, inclusive_b))
        domain.define_parameters(_RealParameter('a', domain=_RealInterval()),
                                 _RealParameter('b', domain=_RealInterval()))
        # Check that domain and string evaluation give the same result
        res = domain.contains(x, dict(a=a, b=b), xp=np)

        # Apparently, `np.float16([2]) < np.float32(2.0009766)` is False
        # but `np.float16([2]) < np.float32([2.0009766])` is True
        # dtype = np.result_type(a.dtype, b.dtype, x.dtype)
        # a, b, x = a.astype(dtype), b.astype(dtype), x.astype(dtype)
        # unclear whether we should be careful about this, since it will be
        # fixed with NEP50. Just do what makes the test pass.
        left_comparison = '<=' if inclusive_a else '<'
        right_comparison = '<=' if inclusive_b else '<'
        ref = eval(f'(a {left_comparison} x) & (x {right_comparison} b)')
        assert_equal(res, ref)

    @pytest.mark.parametrize("inclusive", list(it.product([True, False], repeat=2)))
    @pytest.mark.parametrize("a,b", [(0, 1), (3, 1)])
    def test_contains_function_endpoints(self, inclusive, a, b):
        # Test `contains` when endpoints are defined by functions.
        endpoints = (lambda a, b: (a - b) / 2, lambda a, b: (a + b) / 2)
        domain = _RealInterval(endpoints=endpoints, inclusive=inclusive)
        x = np.asarray([(a - 2*b)/2, (a - b)/2, a/2, (a + b)/2, (a + 2*b)/2])
        params = dict(a=a, b=b)
        res = domain.contains(x, params, xp=np)

        numerical_endpoints = ((a - b) / 2, (a + b) / 2)
        assert numerical_endpoints == domain.get_numerical_endpoints(params, xp=np)
        alpha, beta = numerical_endpoints

        above_left = alpha <= x if inclusive[0] else alpha < x
        below_right = x <= beta if inclusive[1] else x < beta
        ref = above_left & below_right
        assert_equal(res, ref)


    @pytest.mark.parametrize('case', [
        (-np.inf, np.pi, False, True, r"(-\infty, \pi]"),
        ('a', 5, True, False, "[a, 5)")
    ])
    def test_str(self, case):
        domain = _RealInterval(endpoints=case[:2], inclusive=case[2:4])
        assert str(domain) == case[4]

    @pytest.mark.slow
    @given(a=strategies.one_of(
        strategies.decimals(allow_nan=False),
        strategies.characters(whitelist_categories="L"),
        strategies.sampled_from(list(_Domain.symbols))),
           b=strategies.one_of(
        strategies.decimals(allow_nan=False),
        strategies.characters(whitelist_categories="L"),
        strategies.sampled_from(list(_Domain.symbols))),
           inclusive_a=strategies.booleans(),
           inclusive_b=strategies.booleans(),
           )
    def test_str2(self, a, b, inclusive_a, inclusive_b):
        # I wrote this independently from the implementation of __str__, but
        # I imagine it looks pretty similar to __str__.
        a = _Domain.symbols.get(a, a)
        b = _Domain.symbols.get(b, b)
        left_bracket = '[' if inclusive_a else '('
        right_bracket = ']' if inclusive_b else ')'
        domain = _RealInterval(endpoints=(a, b),
                             inclusive=(inclusive_a, inclusive_b))
        ref = f"{left_bracket}{a}, {b}{right_bracket}"
        assert str(domain) == ref

    def test_symbols_gh22137(self):
        # `symbols` was accidentally shared between instances originally
        # Check that this is no longer the case
        domain1 = _RealInterval(endpoints=(0, 1))
        domain2 = _RealInterval(endpoints=(0, 1))
        assert domain1.symbols is not domain2.symbols


def draw_distribution_from_family(family, rng, shape_options, proportions=None, xp=np):
    # If the distribution has parameters, choose a parameterization and
    # draw broadcastable shapes for the parameter arrays.
    proportions = (0.7, 0.1, 0.1, 0.1) if proportions is None else proportions
    n_parameterizations = family._num_parameterizations()
    if n_parameterizations > 0:
        i = rng.integers(0, n_parameterizations)
        n_parameters = family._num_parameters(i)
        shapes = mutually_broadcastable_shapes(num_shapes=n_parameters,
                                               rng=rng, **shape_options)
        result_shape = tuple(xp.broadcast_shapes(*shapes))
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
        return super().test_median(case, tol_override={'atol': 1e-6}, xp=xp)


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
        X.tol = 1e-12
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
        X.tol = 1e-12
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
        X.tol = 1e-12

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
    tol_override = {'atol': 0} if tol_override is None else tol_override

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
        res = getattr(dist, fname)(*args, method=method)
        if 'log' in fname:
            xp_assert_close(xp.exp(res), xp.exp(ref), **tol_override)
        else:
            if fname in {'median', 'mode'} and isinstance(dist, CircularDistribution):
                # wrap to center around 0
                a, b = dist.support()
                period = b - a
                res = xp.where(res > period/2, res - period, res)[()]
                ref = xp.where(ref > period/2, ref - period, ref)[()]
            xp_assert_close(res, ref, **tol_override)

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

    dist = family._draw(dist_shape, rng, xp=np)

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
        Test2(c=[1, object()], d=[1, 2])

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

        X1 = stats.Normal(tol=1e-1, xp=xp)
        X2 = stats.Normal(tol=1e-12, xp=xp)
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
                p = xp.asarray(rng.uniform(size=shape))
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
                                Y.lmoment(order, standardize=standardize))

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

        X = MyBeta(a=xp.asarray(a), b=xp.asarray(b))
        Y = MyBeta(mu=xp.asarray(mu), nu=xp.asarray(nu))

        x = Z.sample(shape=10, rng=rng)
        p = Z.cdf(x)
        x_xp, p_xp = xp.asarray(x), xp.asarray(p)

        xp_assert_close(X.support()[0], xp.asarray(Z.support()[0]))
        xp_assert_close(X.support()[1], xp.asarray(Z.support()[1]))
        xp_assert_close(X.median(), xp.asarray(Z.median()))
        xp_assert_close(X.pdf(x_xp), xp.asarray(Z.pdf(x)))
        xp_assert_close(X.cdf(x_xp), xp.asarray(Z.cdf(x)))
        xp_assert_close(X.ccdf(x_xp), xp.asarray(Z.ccdf(x)))
        xp_assert_close(X.icdf(p_xp), xp.asarray(Z.icdf(p)))
        xp_assert_close(X.iccdf(p_xp), xp.asarray(Z.iccdf(p)))

        xp_assert_close(Y.support()[0], xp.asarray(Z.support()[0]))
        xp_assert_close(Y.support()[1], xp.asarray(Z.support()[1]))
        xp_assert_close(Y.median(), xp.asarray(Z.median()))
        xp_assert_close(Y.pdf(x_xp), xp.asarray(Z.pdf(x)))
        xp_assert_close(Y.cdf(x_xp), xp.asarray(Z.cdf(x)))
        xp_assert_close(Y.ccdf(x_xp), xp.asarray(Z.ccdf(x)))
        xp_assert_close(Y.icdf(p_xp), xp.asarray(Z.icdf(p)))
        xp_assert_close(Y.iccdf(p_xp), xp.asarray(Z.iccdf(p)))

    @pytest.mark.parametrize("n", [10, np.asarray([10., 11., 12.])])
    @pytest.mark.parametrize("p", [0.4, np.asarray([0.4, 0.5, 0.6])])
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
        X = Binomial(n=xp.asarray(n), p=xp.asarray(p))

        x = Y.sample(shape=10, rng=rng)
        p = rng.random(x.shape)
        x_np, p_np = xp.asarray(x), xp.asarray(p)

        xp_assert_close(X.support()[0], xp.asarray(xp_copy(Y.support()[0])))
        xp_assert_close(X.support()[1], xp.asarray(xp_copy(Y.support()[1])))
        xp_assert_close(X.entropy(), xp.asarray(Y.entropy()))
        xp_assert_close(X.median(), xp.asarray(Y.median()))
        xp_assert_close(X.logpdf(x_np), xp.asarray(Y.logpdf(x)))
        xp_assert_close(X.pdf(x_np), xp.asarray(Y.pdf(x)))
        xp_assert_close(X.logpmf(x), xp.asarray(Y.logpmf(x)))
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
        from scipy.stats._distribution_infrastructure import _distribution_names
        for dist in _distribution_names.keys():
            assert hasattr(stats, dist)

        dist = stats.make_distribution(stats.gamma)
        assert str(dist(a=2)) == "Gamma(a=2.0)"
        assert repr(dist(a=2)) == "Gamma(a=np.float64(2.0))"
        assert 'Gamma' in dist.__doc__

        dist = stats.make_distribution(stats.halfgennorm)
        assert str(dist(beta=2)) == "HalfGeneralizedNormal(beta=2.0)"
        assert repr(dist(beta=2)) == "HalfGeneralizedNormal(beta=np.float64(2.0))"
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


@make_xp_test_case(stats._distribution_infrastructure.TransformedDistribution)
class TestTransforms:
    def test_ContinuousDistribution_only(self, xp):
        X = stats.Binomial(n=10, p=0.5, xp=xp)
        # This is applied at the top level TransformedDistribution,
        # so testing one subclass is enough
        message = "Transformations are currently supported only for..."
        with pytest.raises(NotImplementedError, match=message):
            stats.exp(X)


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


@make_xp_test_case(stats._distribution_infrastructure.ShiftedScaledDistribution)
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


@make_xp_test_case(stats._distribution_infrastructure.TransformedDistribution)
class TestReciprocal:
    @pytest.mark.fail_slow(10)
    @pytest.mark.parametrize('scale', [1, 2, -1])
    @pytest.mark.xfail_on_32bit("`scale=-1` fails on 32-bit; needs investigation")
    def test_reciprocal(self, scale, xp):
        if scale == -1 and is_cupy(xp):
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
        Y = scale / X

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


@make_xp_test_case(stats._distribution_infrastructure.MonotonicTransformedDistribution)
class TestMonotonicTransform:
    def test_monotonic_transforms_str(self):
        X = Uniform(a=1, b=2)
        X_str = "Uniform(a=1.0, b=2.0)"

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


@make_xp_test_case(stats._distribution_infrastructure.TransformedDistribution)
class TestArithmetic:
    def test_arithmetic_operators(self, xp):
        rng = np.random.default_rng(2348923495832349834)

        a, b, loc, scale = 0.294, 1.34, 0.57, 1.16

        x = rng.uniform(-3, 3, 100)
        Y = _LogUniform(a=a, b=b, xp=xp)

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


@make_xp_test_case(stats._distribution_infrastructure.TransformedDistribution)
class TestPow:
    def test_pow(self, xp):
        rng = np.random.default_rng(81345982345826)

        Y = Normal(xp=xp)**2
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


@make_xp_test_case(stats.order_statistic)
class TestOrderStatistic:
    @pytest.mark.fail_slow(20)  # Moments require integration
    def test_order_statistic(self, xp):
        rng = np.random.default_rng(7546349802439582)
        X = Uniform(a=0, b=1, xp=xp)
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
                        xp.asarray(Y.cdf(y, method='quadrature')))
        xp_assert_close(Y.ccdf(y_xp, method='formula'),
                        xp.asarray(Y.ccdf(y, method='quadrature')))
        xp_assert_close(Y.icdf(p_xp, method='formula'),
                        xp.asarray(Y.icdf(p, method='inversion')))
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
        xp_assert_close(Z1.cdf(xp.asarray(x)), xp.asarray(Z2.cdf(x)))


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

        lam = 5  # test `_derivative`` `initial_step` logic
        X = TukeyLambda1(lam=lam, xp=xp)
        Y = TukeyLambda2(lam=lam)

        x = np.linspace(-1/lam, 1/lam, 10)
        p = Y.cdf(x)
        x_xp, p_xp = xp.asarray(x), xp.asarray(p)

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


class TestFullCoverage:
    # Adds tests just to get to 100% test coverage; this way it's more obvious
    # if new lines are untested.
    def test_Domain(self):
        with pytest.raises(NotImplementedError):
            _Domain.contains(None, 1.)
        with pytest.raises(NotImplementedError):
            _Domain.get_numerical_endpoints(None, 1.)
        with pytest.raises(NotImplementedError):
            _Domain.__str__(None)

    def test_Parameter(self):
        with pytest.raises(NotImplementedError):
            _Parameter.validate(None, 1.)

    @pytest.mark.parametrize(("dtype_in", "dtype_out"),
                              [(np.float16, np.float16),
                               (np.int16, np.float64)])
    def test_RealParameter_uncommon_dtypes(self, dtype_in, dtype_out):
        domain = _RealInterval((-1, 1))
        parameter = _RealParameter('x', domain=domain)

        x = np.asarray([0.5, 2.5], dtype=dtype_in)
        arr, dtype, valid = parameter.validate(x, parameter_values={}, xp=np)
        assert_equal(arr, x)
        assert dtype == dtype_out
        assert_equal(valid, [True, False])

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
        assert repr(X) == "Uniform(a=np.float64(0.0), b=np.float64(1.0))"
        assert repr(X*3 + 2) == (
            "np.float64(3.0)*Uniform(a=np.float64(0.0), b=np.float64(1.0))"
            " + np.float64(2.0)"
        )

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
        x = xp.asarray(rng.random(shape))

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
                                atol=1e-15)
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
        xp_assert_close(Y.pdf(xp.asarray(x)), xp.asarray(stats.dgamma(a=a).pdf(x)))

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
        p = xp.asarray(rng.random(size=10))
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
        xp_assert_close(Z.logpmf(0), xp.asarray(np.log(0.5)))
        xp_assert_close(Z.logpmf(1), xp.asarray(-np.inf))
        xp_assert_close(Z.pmf(0), xp.asarray(0.5))
        xp_assert_close(Z.pmf(1), xp.asarray(0.))
        xp_assert_close(Z.logpdf(0), xp.asarray(np.inf))
        xp_assert_close(Z.logpdf(1), xp.asarray(np.log(Y.pdf(1)/2)))
        xp_assert_close(Z.pdf(0), xp.asarray(np.inf))
        xp_assert_close(Z.pdf(1), xp.asarray(Y.pdf(1)/2))
        xp_assert_close(Z.logcdf(0), xp.asarray(np.log(0.75)))
        xp_assert_close(Z.logcdf(-1e-10), xp.asarray(np.log(0.25)))
        xp_assert_close(Z.logcdf(1e-10), xp.asarray(np.log(0.75)))
        xp_assert_close(Z.cdf(0), xp.asarray(0.75))
        xp_assert_close(Z.cdf(-1e-10), xp.asarray(0.25))
        xp_assert_close(Z.cdf(1e-10), xp.asarray(0.75))
        xp_assert_close(Z.logccdf(0), xp.asarray(np.log(0.25)))
        xp_assert_close(Z.logccdf(-1e-10), xp.asarray(np.log(0.75)))
        xp_assert_close(Z.logccdf(1e-10), xp.asarray(np.log(0.25)))
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


@make_xp_test_case(_logexpxmexpy)
class Test_logexpxmexpy:
    # Regression tests for some cases that the simplest version of `_logexpmexpy`
    # did not originally handle correctly.
    # This is a private function. Internally, we want NumPy arrays, not NumPy scalars,
    # hence `check_0d=True`.

    def test_x_equals_y(self, xp):
        # Test x - x == 0 in log-space
        x = xp.asarray(2.)
        xp_assert_equal(_logexpxmexpy(x, x, xp=xp), xp.asarray(-xp.inf+0j),
                        check_0d=True)

    def test_y_neg_inf(self, xp):
        # Test x - 0 == x in log-space
        x, y = xp.asarray(2.), xp.asarray(-inf)
        xp_assert_equal(_logexpxmexpy(x, y, xp=xp), xp.asarray(x+0j),
                        check_0d=True)

    def test_nan(self, xp):
        # operations involving NaNs should not produce warnings
        x = xp.asarray(xp.nan)
        xp_assert_equal(_logexpxmexpy(x, x, xp=xp), xp.asarray(x+0j),
                        check_0d=True)
