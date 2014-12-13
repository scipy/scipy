from __future__ import division, print_function, absolute_import

import math
import numpy as np

from numpy.testing import assert_allclose, TestCase, run_module_suite, \
     assert_

from scipy.optimize import fmin_cobyla, minimize


class TestCobyla(TestCase):
    def setUp(self):
        self.x0 = [4.95,0.66]
        self.solution = [math.sqrt(25 - (2.0/3)**2), 2.0/3]
        self.opts = {'disp': False, 'rhobeg': 1, 'tol': 1e-5,
                     'maxiter': 100}

    def fun(self, x):
        return x[0]**2 + abs(x[1])**3

    def con1(self, x):
        return x[0]**2 + x[1]**2 - 25

    def con2(self, x):
        return -self.con1(x)

    def test_simple(self):
        x = fmin_cobyla(self.fun, self.x0, [self.con1, self.con2], rhobeg=1,
                        rhoend=1e-5, iprint=0, maxfun=100)
        assert_allclose(x, self.solution, atol=1e-4)

    def test_minimize_simple(self):
        # Minimize with method='COBYLA'
        cons = ({'type': 'ineq', 'fun': self.con1},
                {'type': 'ineq', 'fun': self.con2})
        sol = minimize(self.fun, self.x0, method='cobyla', constraints=cons,
                       options=self.opts)
        assert_allclose(sol.x, self.solution, atol=1e-4)
        assert_(sol.success, sol.message)
        assert_(sol.maxcv < 1e-5, sol)
        assert_(sol.nfev < 70, sol)
        assert_(sol.fun < self.fun(self.solution) + 1e-3, sol)

    def test_minimize_constraint_violation(self):
        np.random.seed(1234)
        pb = np.random.rand(10, 10)
        spread = np.random.rand(10)

        def p(w):
            return pb.dot(w)

        def f(w):
            return -(w * spread).sum()

        def c1(w):
            return 500 - abs(p(w)).sum()

        def c2(w):
            return 5 - abs(p(w).sum())

        def c3(w):
            return 5 - abs(p(w)).max()

        cons = ({'type': 'ineq', 'fun': c1},
                {'type': 'ineq', 'fun': c2},
                {'type': 'ineq', 'fun': c3})
        w0 = np.zeros((10, 1))
        sol = minimize(f, w0, method='cobyla', constraints=cons,
                       options={'catol': 1e-6})
        assert_(sol.maxcv > 1e-6)
        assert_(not sol.success)


if __name__ == "__main__":
    run_module_suite()
