"""
Test functions for models.rlm
"""

import numpy.random as R
from numpy.testing import NumpyTest, NumpyTestCase

import scipy.sandbox.models.rlm as models

W = R.standard_normal

class test_Regression(NumpyTestCase):

    def test_Robust(self):
        X = W((40,10))
        Y = W((40,))
        model = models(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 30)

    def test_Robustdegenerate(self):
        X = W((40,10))
        X[:,0] = X[:,1] + X[:,2]
        Y = W((40,))
        model = models(design=X)
        results = model.fit(Y)
        self.assertEquals(results.df_resid, 31)

if __name__ == "__main__":
    NumpyTest().run()
