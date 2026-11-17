import itertools as it

import numpy as np
from numpy import inf
import pytest
from numpy.testing import assert_equal
from hypothesis import strategies, given
import hypothesis.extra.numpy as npst

from scipy.stats._univariate_distribution._univariate_distribution_util import (
    _Domain, _RealInterval, _Parameter, _RealParameter, _logexpxmexpy)
from scipy._lib._array_api import make_xp_test_case
from scipy._lib._array_api_no_0d import xp_assert_equal

class Test_RealInterval:
    rng = np.random.default_rng(349849812549824)

    def test_iv(self):
        domain = _RealInterval(endpoints=('a', 'b'))
        message = "The endpoints of the distribution are defined..."
        with pytest.raises(TypeError, match=message):
            domain.get_numerical_endpoints(dict, xp=np)

    @pytest.mark.parametrize('x', [rng.uniform(10, 10, size=(2, 3, 4)),
                                   -np.asarray(np.inf), np.asarray(np.pi)])
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


class TestDomain:
    def test_Domain(self):
        with pytest.raises(NotImplementedError):
            _Domain.contains(None, 1.)
        with pytest.raises(NotImplementedError):
            _Domain.get_numerical_endpoints(None, 1.)
        with pytest.raises(NotImplementedError):
            _Domain.__str__(None)


class TestParameter:
    def test_Parameter(self):
        with pytest.raises(NotImplementedError):
            _Parameter.validate(None, 1.)


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
