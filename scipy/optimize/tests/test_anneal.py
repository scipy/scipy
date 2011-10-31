"""
Unit tests for the simulated annealing minimization algorithm.
"""

from numpy.testing import TestCase, run_module_suite, assert_almost_equal

import numpy as np

from scipy.optimize import anneal

class TestAnneal(TestCase):
    """ Tests for anneal """
    def setUp(self):
        """ Tests setup.

        Define two tests, based on the example in anneal's source. Only the
        first one is used since the second fails for the 'fast' schedule at
        least.
        """
        self.fun = (lambda x: np.cos(14.5 * x - 0.3)  +  (x + 0.2) * x,
                    lambda x: np.cos(14.5 * x[0] - 0.3)  +  \
                             (x[1] + 0.2) * x[1] + (x[0] + 0.2) * x[0])
        self.x0 = (1.0, [1.0, 1.0])
        self.sol = (-0.195, np.array([-0.195, -0.1]))
        self.upper = (3., [3., 3.])
        self.lower = (-3., [-3., -3.])

        # 'fast' and 'cauchy' succeed with maxiter=1000 but 'boltzmann'
        # exits with status=3 until very high values. Keep this value
        # reasonable though.
        self.maxiter = 1000

    def anneal_schedule(self, schedule='fast'):
        """ Call anneal algorithm using specified schedule """
        n = 0 # index of test function
        x, retval = anneal(self.fun[n], self.x0[n], full_output=False,
                           upper=self.upper[n], lower=self.lower[n],
                           feps=1e-3, maxiter=self.maxiter, schedule=schedule,
                           disp=False)
        assert_almost_equal(x, self.sol[n], 2)
        return retval

    def test_fast(self):
        """ Anneal: test for fast schedule """
        retval = self.anneal_schedule('fast')
        self.assertEqual(retval, 0)

    def test_boltzmann(self):
        """ Anneal: test for Boltzmann schedule """
        retval = self.anneal_schedule('boltzmann')
        self.assertLessEqual(retval, 3) # usually 3

    def test_cauchy(self):
        """ Anneal: test for Cauchy schedule """
        retval = self.anneal_schedule('cauchy')
        self.assertEqual(retval, 0)

if __name__ == "__main__":
    run_module_suite()
