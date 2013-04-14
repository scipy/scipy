from __future__ import division, print_function, absolute_import

import math

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
        """ fmin_cobyla """
        x = fmin_cobyla(self.fun, self.x0, [self.con1, self.con2], rhobeg=1,
                        rhoend=1e-5, iprint=0, maxfun=100)
        assert_allclose(x, self.solution, atol=1e-4)

    def test_minimize_simple(self):
        """ Minimize with method='COBYLA' """
        cons = ({'type': 'ineq', 'fun': self.con1},
                {'type': 'ineq', 'fun': self.con2})
        sol = minimize(self.fun, self.x0, method='cobyla', constraints=cons,
                       options=self.opts)
        assert_allclose(sol.x, self.solution, atol=1e-4)
        assert_(sol.success, sol.message)
        assert_(sol.maxcv < 1e-5, sol)
        assert_(sol.nfev < 70, sol)
        assert_(sol.fun < self.fun(self.solution) + 1e-3, sol)

if __name__ == "__main__":
    run_module_suite()
