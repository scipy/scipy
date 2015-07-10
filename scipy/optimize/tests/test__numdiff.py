from __future__ import division

import math

import numpy as np
from numpy.testing import (assert_raises, assert_allclose, assert_equal,
                           assert_, run_module_suite)

from scipy.optimize._numdiff import (
    _adjust_scheme_to_bounds, approx_derivative, check_derivative)


class TestAdjustSchemeToBounds(object):
    def test_no_bounds(self):
        x0 = np.zeros(3)
        h = np.ones(3) * 1e-2
        inf_lower = np.empty_like(x0)
        inf_upper = np.empty_like(x0)
        inf_lower.fill(-np.inf)
        inf_upper.fill(np.inf)

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '1-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '1-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '2-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(~one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '2-sided', inf_lower, inf_upper)
        assert_allclose(h_adjusted, h)
        assert_(np.all(~one_sided))

    def test_with_bound(self):
        x0 = np.array([0.0, 0.85, -0.85])
        lb = -np.ones(3)
        ub = np.ones(3)
        h = np.array([1, 1, -1]) * 1e-1

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 1, '1-sided', lb, ub)
        assert_allclose(h_adjusted, h)

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 2, '1-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([1, -1, 1]) * 1e-1)

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.abs(h))
        assert_(np.all(~one_sided))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([1, -1, 1]) * 1e-1)
        assert_equal(one_sided, np.array([False, True, True]))

    def test_tight_bounds(self):
        lb = np.array([-0.03, -0.03])
        ub = np.array([0.05, 0.05])
        x0 = np.array([0.0, 0.03])
        h = np.array([-0.1, -0.1])

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 1, '1-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.05, -0.06]))

        h_adjusted, _ = _adjust_scheme_to_bounds(x0, h, 2, '1-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.025, -0.03]))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 1, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.03, -0.03]))
        assert_equal(one_sided, np.array([False, True]))

        h_adjusted, one_sided = _adjust_scheme_to_bounds(
            x0, h, 2, '2-sided', lb, ub)
        assert_allclose(h_adjusted, np.array([0.015, -0.015]))
        assert_equal(one_sided, np.array([False, True]))


class TestApproxDerivatives(object):
    def fun_scalar_scalar(self, x):
        return np.sinh(x)

    def jac_scalar_scalar(self, x):
        return np.cosh(x)

    def fun_scalar_vector(self, x):
        return np.array([x[0]**2, np.tan(x[0]), np.exp(x[0])])

    def jac_scalar_vector(self, x):
        return np.array(
            [2 * x[0], np.cos(x[0]) ** -2, np.exp(x[0])]).reshape(-1, 1)

    def fun_vector_scalar(self, x):
        return np.sin(x[0] * x[1]) * np.log(x[0])

    def wrong_dimensions_fun(self, x):
        return np.array([x**2, np.tan(x), np.exp(x)])

    def jac_vector_scalar(self, x):
        return np.array([
            x[1] * np.cos(x[0] * x[1]) * np.log(x[0]) +
            np.sin(x[0] * x[1]) / x[0],
            x[0] * np.cos(x[0] * x[1]) * np.log(x[0])
        ])

    def fun_vector_vector(self, x):
        return np.array([
            x[0] * np.sin(x[1]),
            x[1] * np.cos(x[0]),
            x[0] ** 3 * x[1] ** -0.5
        ])

    def jac_vector_vector(self, x):
        return np.array([
            [np.sin(x[1]), x[0] * np.cos(x[1])],
            [-x[1] * np.sin(x[0]), np.cos(x[0])],
            [3 * x[0] ** 2 * x[1] ** -0.5, -0.5 * x[0] ** 3 * x[1] ** -1.5]
        ])

    def fun_parametrized(self, x, c0, c1=1.0):
        return np.array([np.exp(c0 * x[0]), np.exp(c1 * x[1])])

    def jac_parametrized(self, x, c0, c1=0.1):
        return np.array([
            [c0 * np.exp(c0 * x[0]), 0],
            [0, c1 * np.exp(c1 * x[1])]
        ])

    def fun_with_nan(self, x):
        return x if np.abs(x) <= 1e-8 else np.nan

    def jac_with_nan(self, x):
        return 1.0 if np.abs(x) <= 1e-8 else np.nan

    def fun_zero_jacobian(self, x):
        return np.array([x[0] * x[1], np.cos(x[0] * x[1])])

    def jac_zero_jacobian(self, x):
        return np.array([
            [x[1], x[0]],
            [-x[1] * np.sin(x[0] * x[1]), -x[0] * np.sin(x[0] * x[1])]
        ])

    def fun_non_numpy(self, x):
        return math.exp(x)

    def jac_non_numpy(self, x):
        return math.exp(x)

    def test_scalar_scalar(self):
        x0 = 1.0
        jac_diff_2 = approx_derivative(self.fun_scalar_scalar, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.fun_scalar_scalar, x0)
        jac_true = self.jac_scalar_scalar(x0)
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

    def test_scalar_vector(self):
        x0 = 0.5
        jac_diff_2 = approx_derivative(self.fun_scalar_vector, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.fun_scalar_vector, x0)
        jac_true = self.jac_scalar_vector(np.atleast_1d(x0))
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

    def test_vector_scalar(self):
        x0 = np.array([100.0, -0.5])
        jac_diff_2 = approx_derivative(self.fun_vector_scalar, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.fun_vector_scalar, x0)
        jac_true = self.jac_vector_scalar(x0)
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-7)

    def test_vector_vector(self):
        x0 = np.array([-100.0, 0.2])
        jac_diff_2 = approx_derivative(self.fun_vector_vector, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.fun_vector_vector, x0)
        jac_true = self.jac_vector_vector(x0)
        assert_allclose(jac_diff_2, jac_true, rtol=1e-5)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-6)

    def test_wrong_dimensions(self):
        x0 = 1.0
        assert_raises(RuntimeError, approx_derivative,
                      self.wrong_dimensions_fun, x0)
        f0 = self.wrong_dimensions_fun(np.atleast_1d(x0))
        assert_raises(ValueError, approx_derivative,
                      self.wrong_dimensions_fun, x0, f0=f0)

    def test_custom_rel_step(self):
        x0 = np.array([-0.1, 0.1])
        jac_diff_2 = approx_derivative(self.fun_vector_vector, x0,
                                       method='2-point', rel_step=1e-4)
        jac_diff_3 = approx_derivative(self.fun_vector_vector, x0,
                                       rel_step=1e-4)
        jac_true = self.jac_vector_vector(x0)
        assert_allclose(jac_diff_2, jac_true, rtol=1e-2)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-4)

    def test_options(self):
        x0 = np.array([1.0, 1.0])
        c0 = -1.0
        c1 = 1.0
        lb = 0.0
        ub = 2.0
        f0 = self.fun_parametrized(x0, c0, c1=c1)
        rel_step = np.array([-1e-6, 1e-7])
        jac_true = self.jac_parametrized(x0, c0, c1)
        jac_diff_2 = approx_derivative(
            self.fun_parametrized, x0, method='2-point', rel_step=rel_step,
            f0=f0, args=(c0,), kwargs=dict(c1=c1), bounds=(lb, ub))
        jac_diff_3 = approx_derivative(
            self.fun_parametrized, x0, rel_step=rel_step,
            f0=f0, args=(c0,), kwargs=dict(c1=c1), bounds=(lb, ub))
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

    def test_with_bounds_2_point(self):
        lb = -np.ones(2)
        ub = np.ones(2)

        x0 = np.array([-2.0, 0.2])
        assert_raises(ValueError, approx_derivative,
                      self.fun_vector_vector, x0, bounds=(lb, ub))

        x0 = np.array([-1.0, 1.0])
        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     method='2-point', bounds=(lb, ub))
        jac_true = self.jac_vector_vector(x0)
        assert_allclose(jac_diff, jac_true, rtol=1e-6)

    def test_with_bounds_3_point(self):
        lb = np.array([1.0, 1.0])
        ub = np.array([2.0, 2.0])

        x0 = np.array([1.0, 2.0])
        jac_true = self.jac_vector_vector(x0)

        jac_diff = approx_derivative(self.fun_vector_vector, x0)
        assert_allclose(jac_diff, jac_true, rtol=1e-9)

        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     bounds=(lb, np.inf))
        assert_allclose(jac_diff, jac_true, rtol=1e-9)

        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     bounds=(-np.inf, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-9)

        jac_diff = approx_derivative(self.fun_vector_vector, x0,
                                     bounds=(lb, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-9)

    def test_tight_bounds(self):
        x0 = np.array([10.0, 10.0])
        lb = x0 - 3e-9
        ub = x0 + 2e-9
        jac_true = self.jac_vector_vector(x0)
        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, method='2-point', bounds=(lb, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-6)
        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, method='2-point',
            rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-6)

        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, bounds=(lb, ub))
        assert_allclose(jac_diff, jac_true, rtol=1e-6)
        jac_diff = approx_derivative(
            self.fun_vector_vector, x0, rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_true, jac_diff, rtol=1e-6)

    def test_bound_switches(self):
        lb = -1e-8
        ub = 1e-8
        x0 = 0.0
        jac_true = self.jac_with_nan(x0)
        jac_diff_2 = approx_derivative(
            self.fun_with_nan, x0, method='2-point', rel_step=1e-6,
            bounds=(lb, ub))
        jac_diff_3 = approx_derivative(
            self.fun_with_nan, x0, rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

        x0 = 1e-8
        jac_true = self.jac_with_nan(x0)
        jac_diff_2 = approx_derivative(
            self.fun_with_nan, x0, method='2-point', rel_step=1e-6,
            bounds=(lb, ub))
        jac_diff_3 = approx_derivative(
            self.fun_with_nan, x0, rel_step=1e-6, bounds=(lb, ub))
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-9)

    def test_non_numpy(self):
        x0 = 1.0
        jac_true = self.jac_non_numpy(x0)
        jac_diff_2 = approx_derivative(self.jac_non_numpy, x0,
                                       method='2-point')
        jac_diff_3 = approx_derivative(self.jac_non_numpy, x0)
        assert_allclose(jac_diff_2, jac_true, rtol=1e-6)
        assert_allclose(jac_diff_3, jac_true, rtol=1e-8)

    def test_check_derivative(self):
        x0 = np.array([-10.0, 10])
        accuracy = check_derivative(self.fun_vector_vector,
                                    self.jac_vector_vector, x0)
        assert_(accuracy < 1e-9)
        accuracy = check_derivative(self.fun_vector_vector,
                                    self.jac_vector_vector, x0)
        assert_(accuracy < 1e-6)

        x0 = np.array([0.0, 0.0])
        accuracy = check_derivative(self.fun_zero_jacobian,
                                    self.jac_zero_jacobian, x0)
        assert_(accuracy == 0)
        accuracy = check_derivative(self.fun_zero_jacobian,
                                    self.jac_zero_jacobian, x0)
        assert_(accuracy == 0)


if __name__ == '__main__':
    run_module_suite()
