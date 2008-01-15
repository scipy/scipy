"""
Test functions for models.GLM
"""

import numpy as N
import numpy.random as R
from scipy.testing import *

import scipy.stats.models as S
import scipy.stats.models.glm as models

W = R.standard_normal

class TestRegression(TestCase):

    def test_Logistic(self):
        X = W((40,10))
        Y = N.greater(W((40,)), 0)
        family = S.family.Binomial()
        cmodel = models(design=X, family=S.family.Binomial())
        results = cmodel.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def test_Logisticdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = N.greater(W((40,)), 0)
        family = S.family.Binomial()
        cmodel = models(design=X, family=S.family.Binomial())
        results = cmodel.fit(Y)
        self.assertEquals(results.df_resid, 31)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
