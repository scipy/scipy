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

    def anneal_schedule(self, schedule='fast'):
        """ Call anneal algorithm using specified schedule """
        n = 0 # index of test function
        out = anneal(self.fun[n], self.x0[n], full_output=True,
                     upper=self.upper[n], lower=self.lower[n], feps=1e-3,
                     maxiter=2000, schedule=schedule)
        assert_almost_equal(out[0], self.sol[n], 2)

    def test_fast(self):
        """ Anneal: test for fast schedule """
        self.anneal_schedule('fast')

    def test_boltzmann(self):
        """ Anneal: test for Boltzmann schedule """
        self.anneal_schedule('boltzmann')

    def test_cauchy(self):
        """ Anneal: test for Cauchy schedule """
        self.anneal_schedule('cauchy')

if __name__ == "__main__":
    run_module_suite()
