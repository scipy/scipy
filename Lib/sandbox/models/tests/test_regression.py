#import unittest

from numpy.random import standard_normal
from numpy.testing import NumpyTest, NumpyTestCase

from scipy.sandbox.models.regression import ols_model, ar_model

W = standard_normal

class test_Regression(NumpyTestCase):

    def testOLS(self):
        X = W((40,10))
        Y = W((40,))
        model = ols_model(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testAR(self):
        X = W((40,10))
        Y = W((40,))
        model = ar_model(design=X, rho=0.4)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testOLSdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = ols_model(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)

    def testARdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = ar_model(design=X, rho=0.9)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)

if __name__ == "__main__":
    NumpyTest().run()
#if __name__ == '__main__':
#    unittest.main()
