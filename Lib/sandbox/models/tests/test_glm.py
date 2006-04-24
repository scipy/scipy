import models as S
import unittest
import numpy.random as R
import numpy as N

W = R.standard_normal

class RegressionTest(unittest.TestCase):

    def testLogistic(self):
        X = W((40,10))
        Y = N.greater(W((40,)), 0)
        family = S.family.Binomial()
        model = S.glm.GeneralizedLinearModel(design=X, family=S.family.Binomial())
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def testLogisticdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = N.greater(W((40,)), 0)
        family = S.family.Binomial()
        model = S.glm.GeneralizedLinearModel(design=X, family=S.family.Binomial())
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)


def suite():
    suite = unittest.makeSuite(RegressionTest)
    return suite
        

if __name__ == '__main__':
    unittest.main()
