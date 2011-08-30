"""Regression tests for optimize.

"""

from numpy.testing import TestCase, run_module_suite, assert_almost_equal
import scipy.optimize

class TestRegression(TestCase):

    def test_newton_x0_is_0(self):
        """Ticket #1074"""

        tgt = 1
        res = scipy.optimize.newton(lambda x: x - 1, 0)
        assert_almost_equal(res, tgt)

    def test_newton_integers(self):
        """Ticket #1214"""
        root = scipy.optimize.newton(lambda x: x**2 - 1, x0=2,
                                    fprime=lambda x: 2*x)
        assert_almost_equal(root, 1.0)

if __name__ == "__main__":
    run_module_suite()
