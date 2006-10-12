import unittest
from numpy.random import standard_normal
from scipy.sandbox.models.regression import OLSModel, ARModel
from numpy.testing import *

W = standard_normal

class test_Regression(ScipyTestCase):

    def testOLS(self):
        X = W((40,10))
        Y = W((40,))
        model = OLSModel(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testAR(self):
        X = W((40,10))
        Y = W((40,))
        model = ARModel(design=X, rho=0.4)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testOLSdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = OLSModel(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)

    def testARdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = ARModel(design=X, rho=0.9)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)

def suite():
    suite = unittest.makeSuite(RegressionTest)
    return suite
        

if __name__ == '__main__':
    ScipyTest.run()
