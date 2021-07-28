"""
Unit test for DIRECT optimization algorithm.
"""
from numpy.testing import (assert_, assert_array_almost_equal,
                           assert_equal)
from numpy.testing._private.utils import assert_almost_equal
import pytest
import numpy as np
from scipy.optimize import minimize


class TestDIRECT:

    MAXFEVAL = 20000
    MAXITER = 6000

    def gp_func(x):
        x1, x2 = x[0], x[1]
        fact1a = (x1 + x2 + 1)**2
        fact1b = 19 - 14*x1 + 3*x1**2 - 14*x2 + 6*x1*x2 + 3*x2**2
        fact1 = 1 + fact1a*fact1b

        fact2a = (2*x1 - 3*x2)**2
        fact2b = 18 - 32*x1 + 12*x1**2 + 48*x2 - 36*x1*x2 + 27*x2**2
        fact2 = 30 + fact2a*fact2b

        return fact1*fact2

    def dot_func(x):
        if np.sum(np.abs(x)) > 20:
            return np.nan
        x -= np.array([-1., 2., -4., 3.])
        return np.dot(x, x)

    def neg_inv_func(x):
        if np.sum(x) == 0:
            return -np.inf
        return -1/np.sum(x)

    @pytest.mark.parametrize(
        ("func, bounds, result"), [
         (gp_func, [[-2.0, 2.0], [-2.0, 2.0]],
          {'arg_min': np.array([0.0, -1.0]), 'min': 3.0,
           'arg_decimal': 4, 'decimal': 7,
           'status': 1, 'success': False}),
         (neg_inv_func, 4*[(-10, 10)],
          {'arg_min': np.array([0., 0., 0., 0.]), 'min': -np.inf,
           'arg_decimal': 7, 'decimal': 7,
           'status': 1, 'success': False}),
         (dot_func, 4*[(-10, 10)],
          {'arg_min': np.array([-1., 2., -4., 3.]), 'min': 0.0,
           'arg_decimal': 7, 'decimal': 7,
           'status': 1, 'success': False}),
        ])
    def test_algorithm(self, func, bounds, result):
        res = minimize(func, None, bounds=bounds, method='direct')
        assert_array_almost_equal(res.x, result['arg_min'],
                                  decimal=result['arg_decimal'])
        _bounds = np.asarray(bounds)
        assert_(np.all(res.x >= _bounds[:, 0]))
        assert_(np.all(res.x <= _bounds[:, 1]))
        assert_almost_equal(res.fun, result['min'],
                            decimal=result['decimal'])
        assert_equal(res['success'], result['success'])
        assert_equal(res['status'], result['status'])
        assert_(res['nfev'] >= self.MAXFEVAL)
        assert_(res['nit'] <= self.MAXITER)

        def callback(x):
            print("DIRECT minimization algorithm callback test")

        minimize(func, None, bounds=bounds, method='direct',
                 callback=callback, options={'maxiter': 2})

    def inv(self, x):
        if np.sum(x) == 0:
            raise ZeroDivisionError()
        return 1/np.sum(x)

    def test_exception(self):
        bounds = 4*[(-10, 10)]
        with pytest.raises(ZeroDivisionError):
            minimize(self.inv, None, bounds=bounds,
                     method='direct')
