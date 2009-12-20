"""Regression tests for optimize.

"""
from numpy.testing import *
import numpy as np


class TestRegression(TestCase):
    def test_newton_x0_is_0(self):
        """Ticket #1074"""
        import scipy.optimize
        tgt = 1
        res = scipy.optimize.newton(lambda x: x - 1, 0)
        assert_almost_equal(res, tgt)
