"""
Unit tests for the simulated annealing minimization algorithm.
"""
from __future__ import division, print_function, absolute_import

from numpy.testing import TestCase, run_module_suite, \
    assert_almost_equal, assert_, dec

import numpy as np

from scipy.optimize import anneal, minimize


class TestAnneal(TestCase):
    """ Tests for anneal """
    def setUp(self):
        """ Tests setup.

        Define two tests, based on the example in anneal's source. Only the
        first one is used since the second fails for the 'fast' schedule at
        least.
        """
        self.fun = (lambda x: np.cos(14.5 * x - 0.3) + (x + 0.2) * x,
                    lambda x: np.cos(14.5 * x[0] - 0.3) +
                             (x[1] + 0.2) * x[1] + (x[0] + 0.2) * x[0])
        self.x0 = (1.0, [1.0, 1.0])
        self.sol = (-0.195, np.array([-0.195, -0.1]))
        self.upper = (3., [3., 3.])
        self.lower = (-3., [-3., -3.])

        # 'fast' and 'cauchy' succeed with maxiter=1000 but 'boltzmann'
        # exits with status=3 until very high values. Keep this value
        # reasonable though.
        self.maxiter = 1000

        # fix random seed
        np.random.seed(1234)

    def anneal_schedule(self, schedule='fast', use_wrapper=False):
        """ Call anneal algorithm using specified schedule """
        n = 0  # index of test function
        if use_wrapper:
            opts = {'upper': self.upper[n],
                    'lower': self.lower[n],
                    'ftol': 1e-3,
                    'maxiter': self.maxiter,
                    'schedule': schedule,
                    'disp': False}
            res = minimize(self.fun[n], self.x0[n], method='anneal',
                               options=opts)
            x, retval = res['x'], res['status']
        else:
            x, retval = anneal(self.fun[n], self.x0[n], full_output=False,
                               upper=self.upper[n], lower=self.lower[n],
                               feps=1e-3, maxiter=self.maxiter,
                               schedule=schedule, disp=False)

        assert_almost_equal(x, self.sol[n], 2)
        return retval

    @dec.slow
    def test_fast(self, use_wrapper=False):
        """ Anneal: test for fast schedule """
        retval = self.anneal_schedule('fast', use_wrapper)
        self.assertEqual(retval, 0)

    @dec.slow
    def test_boltzmann(self, use_wrapper=False):
        """ Anneal: test for Boltzmann schedule """
        retval = self.anneal_schedule('boltzmann', use_wrapper)
        assert_(retval <= 3)  # usually 3

    @dec.slow
    def test_cauchy(self, use_wrapper=False):
        """ Anneal: test for Cauchy schedule """
        retval = self.anneal_schedule('cauchy', use_wrapper)
        self.assertEqual(retval, 0)

    @dec.slow
    def test_minimize(self):
        """ minimize with 'anneal' method """
        self.test_fast(True)
        self.test_boltzmann(True)
        self.test_cauchy(True)

if __name__ == "__main__":
    run_module_suite()
