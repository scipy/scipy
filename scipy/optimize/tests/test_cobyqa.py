import numpy as np
from numpy.testing import assert_allclose, assert_

from scipy.optimize import (
    Bounds,
    LinearConstraint,
    NonlinearConstraint,
    minimize,
)


class TestCOBYQA:

    def setup_method(self):
        self.x0 = [4.95, 0.66]
        self.options = {'maxfev': 100}

    @staticmethod
    def fun(x, c=1.0):
        return x[0]**2 + c * abs(x[1])**3

    @staticmethod
    def con(x):
        return x[0]**2 + x[1]**2 - 25.0

    def test_minimize_simple(self):
        class Callback:
            def __init__(self):
                self.n_calls = 0

            def __call__(self, x):
                self.n_calls += 1

        callback = Callback()

        # Minimize with method='cobyqa'.
        constraints = NonlinearConstraint(self.con, 0.0, 0.0)
        sol = minimize(
            self.fun,
            self.x0,
            method='cobyqa',
            constraints=constraints,
            callback=callback,
            options=self.options,
        )
        solution = [np.sqrt(25.0 - 4.0 / 9.0), 2.0 / 3.0]
        assert_allclose(sol.x, solution, atol=1e-4)
        assert_(sol.success, sol.message)
        assert_(sol.maxcv < 1e-8, sol)
        assert_(sol.nfev <= 100, sol)
        assert_(sol.fun < self.fun(solution) + 1e-3, sol)
        assert_(
            sol.nfev == callback.n_calls,
            "Callback is not called exactly once for every function eval.",
        )

    def test_minimize_bounds(self):
        # Case where the bounds are not active at the solution.
        bounds = Bounds([4.5, 0.6], [5.0, 0.7])
        constraints = NonlinearConstraint(self.con, 0.0, 0.0)
        sol = minimize(
            self.fun,
            self.x0,
            method='cobyqa',
            bounds=bounds,
            constraints=constraints,
            options=self.options,
        )
        solution = [np.sqrt(25.0 - 4.0 / 9.0), 2.0 / 3.0]
        assert_allclose(sol.x, solution, atol=1e-4)
        assert_(sol.success, sol.message)
        assert_(sol.maxcv < 1e-8, sol)
        assert_(sol.nfev <= 100, sol)
        assert_(sol.fun < self.fun(solution) + 1e-3, sol)

        # Case where the bounds are active at the solution.
        bounds = Bounds([5.0, 0.6], [5.5, 0.65])
        sol = minimize(
            self.fun,
            self.x0,
            method='cobyqa',
            bounds=bounds,
            constraints=constraints,
            options=self.options,
        )
        assert_(not sol.success, sol.message)
        assert_(sol.maxcv > 0.35, sol)
        assert_(sol.nfev <= 100, sol)

    def test_minimize_linear_constraints(self):
        constraints = LinearConstraint([1.0, 1.0], 1.0, 1.0)
        sol = minimize(
            self.fun,
            self.x0,
            method='cobyqa',
            constraints=constraints,
            options=self.options,
        )
        solution = [(4 - np.sqrt(7)) / 3, (np.sqrt(7) - 1) / 3]
        assert_allclose(sol.x, solution, atol=1e-4)
        assert_(sol.success, sol.message)
        assert_(sol.maxcv < 1e-8, sol)
        assert_(sol.nfev <= 100, sol)
        assert_(sol.fun < self.fun(solution) + 1e-3, sol)

    def test_minimize_args(self):
        constraints = NonlinearConstraint(self.con, 0.0, 0.0)
        sol = minimize(
            self.fun,
            self.x0,
            (2.0,),
            'cobyqa',
            constraints=constraints,
            options=self.options,
        )
        solution = [np.sqrt(25.0 - 4.0 / 36.0), 2.0 / 6.0]
        assert_allclose(sol.x, solution, atol=1e-4)
        assert_(sol.success, sol.message)
        assert_(sol.maxcv < 1e-8, sol)
        assert_(sol.nfev <= 100, sol)
        assert_(sol.fun < self.fun(solution, 2.0) + 1e-3, sol)
