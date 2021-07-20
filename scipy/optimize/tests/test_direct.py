"""
Unit test for DIRECT optimization algorithm.
"""
from numpy.testing import (assert_, assert_array_almost_equal,
                           assert_allclose, assert_equal)
from pytest import raises as assert_raises
import pytest
import numpy as np

from scipy.optimize import minimize

class TestDIRECT:

    def get_test_func(self, name):
        if name == 'goldstein_price_function':
            def gp_func(x):
                x1, x2 = x[0], x[1]
                fact1a = (x1 + x2 + 1)**2
                fact1b = 19 - 14*x1 + 3*x1**2 - 14*x2 + 6*x1*x2 + 3*x2**2
                fact1 = 1 + fact1a*fact1b

                fact2a = (2*x1 - 3*x2)**2
                fact2b = 18 - 32*x1 + 12*x1**2 + 48*x2 - 36*x1*x2 + 27*x2**2
                fact2 = 30 + fact2a*fact2b

                return fact1*fact2
            return gp_func
    
    @pytest.mark.parametrize(
        ("function_name, bounds, arg_min, min, "
         "arg_min_rtol, min_rtol, arg_min_atol, min_atol"), [
        ('goldstein_price_function', [[-2.0, 2.0], [-2.0, 2.0]], 
            np.array([0.0, -1.0]), 3.0, 1e-5, 1e-7, 1e-3, 0)
        ])
    def test_bounds(self, function_name, bounds, 
                    arg_min, min, arg_min_rtol, min_rtol,
                    arg_min_atol, min_atol):
        func = self.get_test_func(function_name)
        res = minimize(func, None, bounds=bounds, method='direct')
        assert_allclose(res.x, arg_min,
                        rtol=arg_min_rtol, atol=arg_min_atol)
        assert_allclose(res.fun, min,
                        rtol=min_rtol, atol=min_atol)
        assert_(res['success'] is False)
        assert_equal(res['status'], 1)
        assert_equal(res['nfev'], 20005)
        assert_equal(res['nit'], 460)
