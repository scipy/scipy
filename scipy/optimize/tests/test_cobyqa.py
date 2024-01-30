import math
import numpy as np

from numpy.testing import assert_allclose, assert_, assert_array_equal

from scipy.optimize import minimize, Bounds, NonlinearConstraint


class TestCOBYQA:

    def setup_method(self):
        self.x0 = [4.95, 0.66]
        self.solution = [math.sqrt(25 - (2.0/3)**2), 2.0/3]
        self.opts = {
            'disp': False,
            'radius_init': 1,
            'tol': 1e-5,
            'maxfev': 100,
        }

    def fun(self, x):
        return x[0]**2 + abs(x[1])**3

    def con(self, x):
        return x[0]**2 + x[1]**2 - 25

    def test_minimize_simple(self):
        class Callback:
            def __init__(self):
                self.n_calls = 0

            def __call__(self, x):
                self.n_calls += 1

        callback = Callback()

        # Minimize with method='COBYQA'
        cons = NonlinearConstraint(self.con, 0, 0)
        sol = minimize(self.fun, self.x0, method='cobyqa', constraints=cons,
                       callback=callback, options=self.opts)
        assert_allclose(sol.x, self.solution, atol=1e-4)
        assert_(sol.success, sol.message)
        assert_(sol.maxcv < 1e-5, sol)
        assert_(sol.nfev < 70, sol)
        assert_(sol.fun < self.fun(self.solution) + 1e-3, sol)
        assert_(sol.nfev == callback.n_calls,
                "Callback is not called exactly once for every function eval.")

    def test_minimize_constraint_violation(self):
        np.random.seed(1234)
        pb = np.random.rand(10, 10)
        spread = np.random.rand(10)

        def p(w):
            return pb.dot(w)

        def f(w):
            return -(w * spread).sum()

        def c(w):
            return abs(p(w)).sum()

        cons = (
            NonlinearConstraint(c, -np.inf, 5),
            NonlinearConstraint(c, 6, np.inf),
        )
        w0 = np.zeros((10,))
        sol = minimize(f, w0, method='cobyqa', constraints=cons,
                       options={'feasibility_tol': 1e-6})
        assert_(sol.maxcv > 1e-6)
        assert_(not sol.success)


def test_vector_constraints():
    # test that minimize can take a combination of constraints, some returning a
    # number and others an array
    def fun(x):
        return (x[0] - 1)**2 + (x[1] - 2.5)**2

    def fmin(x):
        return fun(x) - 1

    def cons1(x):
        a = np.array([[1, -2, 2], [-1, -2, 6], [-1, 2, 2]])
        return np.array([a[i, 0] * x[0] + a[i, 1] * x[1] +
                         a[i, 2] for i in range(len(a))])

    def cons2(x):
        return x     # identity, acts as bounds x > 0

    x0 = np.array([2, 0])
    cons_list = [fun, cons1, cons2]

    xsol = [1.4, 1.7]
    fsol = 0.8

    # testing minimize
    constraints = [{'type': 'ineq', 'fun': cons} for cons in cons_list]
    sol = minimize(fun, x0, constraints=constraints, tol=1e-5)
    assert_allclose(sol.x, xsol, atol=1e-4)
    assert_(sol.success, sol.message)
    assert_allclose(sol.fun, fsol, atol=1e-4)

    constraints = {'type': 'ineq', 'fun': fmin}
    sol = minimize(fun, x0, constraints=constraints, tol=1e-5)
    assert_allclose(sol.fun, 1, atol=1e-4)


class TestBounds:
    # Test cobyqa support for bounds
    # Invalid bounds is tested in
    # test_optimize.TestOptimizeSimple.test_minimize_invalid_bounds

    def test_basic(self):
        def f(x):
            return np.sum(x**2)

        lb = [-1, None, 1, None, -0.5]
        ub = [-0.5, -0.5, None, None, -0.5]
        bounds = [(a, b) for a, b in zip(lb, ub)]
        # these are converted to Bounds internally

        res = minimize(f, x0=[1, 2, 3, 4, 5], method='cobyqa', bounds=bounds)
        ref = [-0.5, -0.5, 1, 0, -0.5]
        assert res.success
        assert_allclose(res.x, ref, atol=1e-3)

    def test_unbounded(self):
        def f(x):
            return np.sum(x**2)

        bounds = Bounds([-np.inf, -np.inf], [np.inf, np.inf])
        res = minimize(f, x0=[1, 2], method='cobyqa', bounds=bounds)
        assert res.success
        assert_allclose(res.x, 0, atol=1e-3)

        bounds = Bounds([1, -np.inf], [np.inf, np.inf])
        res = minimize(f, x0=[1, 2], method='cobyqa', bounds=bounds)
        assert res.success
        assert_allclose(res.x, [1, 0], atol=1e-3)
