"""
Test functions for models.utils
"""

import numpy as N
import numpy.random as R
from numpy.testing import *

from scipy.stats.models import utils

class TestUtils(TestCase):

    def test_recipr(self):
        X = N.array([[2,1],[-1,0]])
        Y = utils.recipr(X)
        assert_almost_equal(Y, N.array([[0.5,1],[0,0]]))

    def test_recipr0(self):
        X = N.array([[2,1],[-4,0]])
        Y = utils.recipr0(X)
        assert_almost_equal(Y, N.array([[0.5,1],[-0.25,0]]))

    def test_rank(self):
        X = R.standard_normal((40,10))
        self.assertEquals(utils.rank(X), 10)

        X[:,0] = X[:,1] + X[:,2]
        self.assertEquals(utils.rank(X), 9)

    def test_fullrank(self):
        X = R.standard_normal((40,10))
        X[:,0] = X[:,1] + X[:,2]

        Y = utils.fullrank(X)
        self.assertEquals(Y.shape, (40,9))
        self.assertEquals(utils.rank(Y), 9)

        X[:,5] = X[:,3] + X[:,4]
        Y = utils.fullrank(X)
        self.assertEquals(Y.shape, (40,8))
        self.assertEquals(utils.rank(Y), 8)

    def test_StepFunction(self):
        x = N.arange(20)
        y = N.arange(20)
        f = utils.StepFunction(x, y)
        assert_almost_equal(f( N.array([[3.2,4.5],[24,-3.1]]) ), [[ 3, 4], [19, 0]])

    def test_StepFunctionBadShape(self):
        x = N.arange(20)
        y = N.arange(21)
        self.assertRaises(ValueError, utils.StepFunction, x, y)
        x = N.zeros((2, 2))
        y = N.zeros((2, 2))
        self.assertRaises(ValueError, utils.StepFunction, x, y)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
