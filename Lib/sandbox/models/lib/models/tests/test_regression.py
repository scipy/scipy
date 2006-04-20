import models as S
import unittest
import numpy.random as R
import numpy as N

W = R.standard_normal

class RegressionTest(unittest.TestCase):

    def testOLS(self):
        X = W((40,10))
        Y = W((40,))
        model = S.regression.OLSModel(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testcovbeta1(self):
        X = W((40,10))
        Y = W((40,))
        model = S.regression.OLSModel(design=X)
        results = model.fit(Y)
        v = N.asarray(results.cov_beta(column=2))
        print v
        self.assertEquals(v.shape, ())

    def testcovbeta2(self):
        X = W((40,10))
        Y = W((40,))
        model = S.regression.OLSModel(design=X)
        results = model.fit(Y)
        v = results.cov_beta(column=[2,3])
        self.assertEquals(v.shape, (2,2))

    def testcovbeta3(self):
        X = W((40,10))
        Y = W((40,))
        model = S.regression.OLSModel(design=X)
        results = model.fit(Y)
        v = results.cov_beta()
        self.assertEquals(v.shape, (10,10))

    def testAR(self):
        X = W((40,10))
        Y = W((40,))
        model = S.regression.ARModel(design=X, rho=0.4)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testOLSdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = S.regression.OLSModel(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)

    def testcovbeta4(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = S.regression.OLSModel(design=X)
        results = model.fit(Y)
        v = results.cov_beta()
        self.assertEquals(S.utils.rank(v), 9)
        self.assertEquals(v.shape, (10,10))

    def testARdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = S.regression.ARModel(design=X, rho=0.9)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)

def suite():
    suite = unittest.makeSuite(RegressionTest)
    return suite
        

if __name__ == '__main__':
    unittest.main()
