"""
Unit test for DIRECT optimization algorithm.
"""
from numpy.testing import (assert_, assert_array_almost_equal,
                           assert_equal)
from numpy.testing._private.utils import assert_almost_equal
import pytest
import numpy as np
from scipy.optimize import direct


class TestDIRECT:

    MAXFEVAL = 20000
    MAXITER = 6000

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
         (neg_inv_func, 4*[(-10, 10)],
          {'arg_min': np.array([0., 0., 0., 0.]), 'min': -np.inf,
           'arg_decimal': 7, 'decimal': 7,
           'status': 4, 'success': True}),
         (dot_func, 4*[(-10, 10)],
          {'arg_min': np.array([-1., 2., -4., 3.]), 'min': 0.0,
           'arg_decimal': 3, 'decimal': 6,
           'status': 4, 'success': True}),
        ])
    def test_algorithm(self, func, bounds, result):
        res = direct(func, bounds=bounds)
        assert_array_almost_equal(res.x, result['arg_min'],
                                  decimal=result['arg_decimal'])
        _bounds = np.asarray(bounds)
        assert_(np.all(res.x >= _bounds[:, 0]))
        assert_(np.all(res.x <= _bounds[:, 1]))
        assert_almost_equal(res.fun, result['min'],
                            decimal=result['decimal'])
        assert_equal(res['success'], result['success'])
        assert_equal(res['status'], result['status'])
        assert_(res['nfev'] <= self.MAXFEVAL)
        assert_(res['nit'] <= self.MAXITER)

        def callback(x):
            print("DIRECT minimization algorithm callback test")

        direct(func, bounds=bounds, callback=callback,
               maxiter=2, disp=True, maxfun=1)

    def inv(self, x):
        if np.sum(x) == 0:
            raise ZeroDivisionError()
        return 1/np.sum(x)

    def test_exception(self):
        bounds = 4*[(-10, 10)]
        with pytest.raises(ZeroDivisionError):
            direct(self.inv, bounds=bounds)
