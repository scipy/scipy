from numpy.testing import (run_module_suite, assert_, assert_allclose,
                           assert_raises, assert_equal)
import numpy as np

from scipy.optimize._lsq.common import (
    step_size_to_bound, find_active_constraints, make_strictly_feasible,
    scaling_vector, intersect_trust_region, build_quadratic_1d,
    minimize_quadratic_1d, evaluate_quadratic)


class TestBounds(object):
    def test_step_size_to_bounds(self):
        lb = np.array([-1.0, 2.5, 10.0])
        ub = np.array([1.0, 5.0, 100.0])
        x = np.array([0.0, 2.5, 12.0])

        s = np.array([0.1, 0.0, 0.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, 10)
        assert_equal(hits, [1, 0, 0])

        s = np.array([0.01, 0.05, -1.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, 2)
        assert_equal(hits, [0, 0, -1])

        s = np.array([10.0, -0.0001, 100.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, np.array(-0))
        assert_equal(hits, [0, -1, 0])

        s = np.array([1.0, 0.5, -2.0])
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, 1.0)
        assert_equal(hits, [1, 0, -1])

        s = np.zeros(3)
        step, hits = step_size_to_bound(x, s, lb, ub)
        assert_equal(step, np.inf)
        assert_equal(hits, [0, 0, 0])

    def test_find_active_constraints(self):
        lb = np.array([0.0, -10.0, 1.0])
        ub = np.array([1.0, 0.0, 100.0])

        x = np.array([0.5, -5.0, 2.0])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [0, 0, 0])

        x = np.array([0.0, 0.0, 10.0])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [-1, 1, 0])

        active = find_active_constraints(x, lb, ub, rtol=0)
        assert_equal(active, [-1, 1, 0])

        x = np.array([1e-9, -1e-8, 100 - 1e-9])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [0, 0, 1])

        active = find_active_constraints(x, lb, ub, rtol=1.5e-9)
        assert_equal(active, [-1, 0, 1])

        lb = np.array([1.0, -np.inf, -np.inf])
        ub = np.array([np.inf, 10.0, np.inf])

        x = np.ones(3)
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [-1, 0, 0])

        # Handles out-of-bound cases.
        x = np.array([0.0, 11.0, 0.0])
        active = find_active_constraints(x, lb, ub)
        assert_equal(active, [-1, 1, 0])

        active = find_active_constraints(x, lb, ub, rtol=0)
        assert_equal(active, [-1, 1, 0])

    def test_make_strictly_feasible(self):
        lb = np.array([-0.5, -0.8, 2.0])
        ub = np.array([0.8, 1.0, 3.0])

        x = np.array([-0.5, 0.0, 2 + 1e-10])

        x_new = make_strictly_feasible(x, lb, ub, rstep=0)
        assert_(x_new[0] > -0.5)
        assert_equal(x_new[1:], x[1:])

        x_new = make_strictly_feasible(x, lb, ub, rstep=1e-4)
        assert_equal(x_new, [-0.5 + 1e-4, 0.0, 2 * (1 + 1e-4)])

        x = np.array([-0.5, -1, 3.1])
        x_new = make_strictly_feasible(x, lb, ub)
        assert_(np.all((x_new >= lb) & (x_new <= ub)))

        x_new = make_strictly_feasible(x, lb, ub, rstep=0)
        assert_(np.all((x_new >= lb) & (x_new <= ub)))

        lb = np.array([-1, 100.0])
        ub = np.array([1, 100.0 + 1e-10])
        x = np.array([0, 100.0])
        x_new = make_strictly_feasible(x, lb, ub, rstep=1e-8)
        assert_equal(x_new, [0, 100.0 + 0.5e-10])

    def test_scaling_vector(self):
        lb = np.array([-np.inf, -5.0, 1.0, -np.inf])
        ub = np.array([1.0, np.inf, 10.0, np.inf])
        x = np.array([0.5, 2.0, 5.0, 0.0])
        g = np.array([1.0, 0.1, -10.0, 0.0])
        v, dv = scaling_vector(x, g, lb, ub)
        assert_equal(v, [1.0, 7.0, 5.0, 1.0])
        assert_equal(dv, [0.0, 1.0, -1.0, 0.0])


class TestQuadraticFunction(object):
    def __init__(self):
        self.J = np.array([
            [0.1, 0.2],
            [-1.0, 1.0],
            [0.5, 0.2]])
        self.g = np.array([0.8, -2.0])
        self.diag = np.array([1.0, 2.0])

    def test_build_quadratic_1d(self):
        s = np.zeros(2)
        a, b = build_quadratic_1d(self.J, self.g, s)
        assert_equal(a, 0)
        assert_equal(b, 0)

        a, b = build_quadratic_1d(self.J, self.g, s, diag=self.diag)
        assert_equal(a, 0)
        assert_equal(b, 0)

        s = np.array([1.0, -1.0])
        a, b = build_quadratic_1d(self.J, self.g, s)
        assert_equal(a, 2.05)
        assert_equal(b, 2.8)

        a, b = build_quadratic_1d(self.J, self.g, s, diag=self.diag)
        assert_equal(a, 3.55)
        assert_equal(b, 2.8)

        s0 = np.array([0.5, 0.5])
        a, b = build_quadratic_1d(self.J, self.g, s, diag=self.diag, s0=s0)
        assert_equal(a, 3.55)
        assert_allclose(b, 2.39)

    def test_minimize_quadratic_1d(self):
        a = 5
        b = -1

        t, y = minimize_quadratic_1d(a, b, 1, 2)
        assert_equal(t, 1)
        assert_equal(y, a * t**2 + b * t)

        t, y = minimize_quadratic_1d(a, b, -2, -1)
        assert_equal(t, -1)
        assert_equal(y, a * t**2 + b * t)

        t, y = minimize_quadratic_1d(a, b, -1, 1)
        assert_equal(t, 0.1)
        assert_equal(y, a * t**2 + b * t)

    def test_evaluate_quadratic(self):
        s = np.array([1.0, -1.0])

        value = evaluate_quadratic(self.J, self.g, s)
        assert_equal(value, 4.85)

        value = evaluate_quadratic(self.J, self.g, s, diag=self.diag)
        assert_equal(value, 6.35)

        s = np.array([[1.0, -1.0],
                     [1.0, 1.0],
                     [0.0, 0.0]])

        values = evaluate_quadratic(self.J, self.g, s)
        assert_allclose(values, [4.85, -0.91, 0.0])

        values = evaluate_quadratic(self.J, self.g, s, diag=self.diag)
        assert_allclose(values, [6.35, 0.59, 0.0])


class TestTrustRegion(object):
    def test_intersect(self):
        Delta = 1.0

        x = np.zeros(3)
        s = np.array([1.0, 0.0, 0.0])
        t_neg, t_pos = intersect_trust_region(x, s, Delta)
        assert_equal(t_neg, -1)
        assert_equal(t_pos, 1)

        s = np.array([-1.0, 1.0, -1.0])
        t_neg, t_pos = intersect_trust_region(x, s, Delta)
        assert_allclose(t_neg, -3**-0.5)
        assert_allclose(t_pos, 3**-0.5)

        x = np.array([0.5, -0.5, 0])
        s = np.array([0, 0, 1.0])
        t_neg, t_pos = intersect_trust_region(x, s, Delta)
        assert_allclose(t_neg, -2**-0.5)
        assert_allclose(t_pos, 2**-0.5)

        x = np.ones(3)
        assert_raises(ValueError, intersect_trust_region, x, s, Delta)

        x = np.zeros(3)
        s = np.zeros(3)
        assert_raises(ValueError, intersect_trust_region, x, s, Delta)


if __name__ == '__main__':
    run_module_suite()
