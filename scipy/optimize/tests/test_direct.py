"""
Unit test for DIRECT optimization algorithm.
"""
from numpy.testing import (assert_allclose,
                           assert_array_less)
import pytest
import numpy as np
from scipy.optimize import direct, Bounds


class TestDIRECT:

    # actual maxeval is always a bit higher than 20000
    MAXFEVAL = 20100
    MAXITER = 6000

    # test functions
    def sphere(x):
        return np.square(x).sum()

    def neg_inv_func(x):
        if np.sum(x) == 0:
            return -np.inf
        return -1/np.sum(x)

    def sphere_2(self, x):
        return np.square(x).sum()

    def styblinski_tang(self, pos):
        x, y = pos
        return 0.5 * (x**4 - 16 * x**2 + 5 * x + y**4 - 16 * y**2 + 5 * y)

    @pytest.mark.parametrize(
        ("func, bounds, result"), [
         (neg_inv_func, 4*[(-10, 10)],
          {'arg_min': np.array([0., 0., 0., 0.]), 'min': -np.inf,
           'status': 4, 'success': True}),
         (sphere, 4*[(-10, 10)],
          {'arg_min': np.zeros((4, )), 'min': 0.0,
           'status': 4, 'success': True}),
        ])
    def test_direct_locally_biased(self, func, bounds, result):
        res = direct(func, bounds=bounds)
        assert_allclose(res.x, result['arg_min'])
        _bounds = np.asarray(bounds)

        # test that result lies within bounds
        assert_array_less(_bounds[:, 0], res.x)
        assert_array_less(res.x, _bounds[:, 1])

        # test accuracy
        assert_allclose(res.fun, result['min'])
        assert res.success == result['success']
        assert res.status == result['status']
        assert res.nfev <= self.MAXFEVAL
        assert res.nit <= self.MAXITER

        def callback(x):
            print("DIRECT minimization algorithm callback test")

        direct(func, bounds=bounds, callback=callback,
               maxiter=2, disp=True, maxfun=1)

    @pytest.mark.parametrize(
        ("func, bounds, result"), [
         (neg_inv_func, 4*[(-10, 10)],
          {'arg_min': np.array([0., 0., 0., 0.]), 'min': -np.inf,
           'status': 4, 'success': True}),
         (sphere, 4*[(-10, 10)],
          {'arg_min': np.zeros((4, )), 'min': 0.0,
           'status': 4, 'success': True}),
        ])
    def test_direct_original(self, func, bounds, result):
        res = direct(func, bounds=bounds, locally_biased=False)
        assert_allclose(res.x, result['arg_min'])
        _bounds = np.asarray(bounds)

        # test that result lies within bounds
        assert_array_less(_bounds[:, 0], res.x)
        assert_array_less(res.x, _bounds[:, 1])

        # test accuracy
        assert_allclose(res.fun, result['min'])
        assert res.success == result['success']
        assert res.status == result['status']
        assert res.nfev <= self.MAXFEVAL
        assert res.nit <= self.MAXITER

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

    def test_len_tol(self):
        bounds = 4*[(-10., 10.)]
        res = direct(self.sphere_2, bounds=bounds, len_tol=1e-3)
        assert res.status == 5
        assert res.success
        assert_allclose(res.x, np.zeros((4, )))

    @pytest.mark.parametrize("f_min_rtol", [1e-3, 1e-5, 1e-7])
    def test_f_min(self, f_min_rtol):
        # test that desired function value is reached within
        # relative tolerance of f_min_rtol
        f_min = 1.
        bounds = 4*[(-2., 10.)]
        res = direct(self.sphere_2, bounds=bounds, f_min=f_min,
                     f_min_rtol=f_min_rtol)
        assert res.status == 3
        assert res.success
        assert res.fun < f_min + f_min_rtol * f_min

    def circle_with_args(self, x, a, b):
        return np.square(x[0] - a) + np.square(x[1] - b).sum()

    def test_f_circle_with_args(self):
        bounds = 2*[(-2.0, 2.0)]

        res = direct(self.circle_with_args, bounds, args=(1, 1), maxfun=1250)
        assert res.nfev < 1250
        assert res.status == 4
        assert res.success
        assert_allclose(res.x, np.array([1., 1.]))

    def test_failure_maxfun(self):
        # test that if optimization runs for the maximal number of
        # evaluations, success = False is returned

        maxfun = 100
        bounds = Bounds([-4., -4.], [4., 4.])
        result = direct(self.styblinski_tang, bounds, maxfun=maxfun)
        assert result.success is False
        assert result.status == 1
        assert result.nfev >= maxfun

    def test_failure_maxiter(self):
        # test that if optimization runs for the maximal number of
        # iterations, success = False is returned

        maxiter = 50
        bounds = Bounds([-4., -4.], [4., 4.])
        result = direct(self.styblinski_tang, bounds, maxiter=maxiter)
        assert result.success is False
        assert result.status == 2
        assert result.nit >= maxiter

    def test_bounds_variants(self):
        # test that new and old bounds yield same result

        lb = [-6., 1., -5.]
        ub = [-1., 3., 5.]
        x_opt = np.array([-1., 1., 0.])
        bounds_old = list(zip(lb, ub))
        bounds_new = Bounds(lb, ub)

        res_old_bounds = direct(self.sphere_2, bounds_old)
        res_new_bounds = direct(self.sphere_2, bounds_new)

        assert res_new_bounds.nfev == res_old_bounds.nfev
        assert res_new_bounds.message == res_old_bounds.message
        assert res_new_bounds.success == res_old_bounds.success
        assert res_new_bounds.nit == res_old_bounds.nit
        assert_allclose(res_new_bounds.x, res_old_bounds.x)
        assert_allclose(res_new_bounds.x, x_opt, rtol=1e-4)

    @pytest.mark.parametrize("eps", [1e-5, 1e-4, 1e-3])
    def test_epsilon(self, eps):
        bounds = Bounds([-4., -4.], [4., 4.])
        result = direct(self.styblinski_tang, bounds, eps=eps, vol_tol=1e-6)
        assert result.status == 4
        assert result.success

    def test_segmentation_fault(self):
        # test that an excessive number of function evaluations
        # does not result in segmentation fault

        bounds = [(-5., 20.)] * 100
        result = direct(self.sphere_2, bounds, maxfun=10000000,
                        maxiter=1000000)
        assert result is not None
