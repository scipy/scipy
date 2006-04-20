import unittest, csv, os
import numpy as N
import numpy.random as R
import scipy

from scipy_stats_models import utils

class UtilsTest(unittest.TestCase):

    def test_recipr(self):
        X = N.array([[2,1],[-1,0]])
        Y = utils.recipr(X)
        scipy.testing.assert_almost_equal(Y, N.array([[0.5,1],[0,0]]))

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

def suite():
    suite = unittest.makeSuite(UtilsTest)
    return suite

if __name__ == '__main__':
    unittest.main()
