import unittest

import numpy.random as R

import scipy.sandbox.models as S

W = R.standard_normal

class RegressionTest(unittest.TestCase):

    def testRobust(self):
        X = W((40,10))
        Y = W((40,))
        model = S.rlm(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testRobustdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = S.rlm(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)


if __name__ == '__main__':
    unittest.main()
