"""
Test functions for models.robust.scale
"""

import numpy.random as R
from scipy.testing import *

import scipy.stats.models.robust.scale as scale

W = R.standard_normal

class TestScale(TestCase):

    def test_MAD(self):
        X = W((40,10))
        m = scale.MAD(X)
        self.assertEquals(m.shape, (10,))

    def test_MADaxes(self):
        X = W((40,10,30))
        m = scale.MAD(X, axis=0)
        self.assertEquals(m.shape, (10,30))

        m = scale.MAD(X, axis=1)
        self.assertEquals(m.shape, (40,30))

        m = scale.MAD(X, axis=2)
        self.assertEquals(m.shape, (40,10))

        m = scale.MAD(X, axis=-1)
        self.assertEquals(m.shape, (40,10))

    def test_huber(self):
        X = W((40,10))
        m = scale.huber(X)
        self.assertEquals(m.shape, (10,))

    def test_huberaxes(self):
        X = W((40,10,30))
        m = scale.huber(X, axis=0)
        self.assertEquals(m.shape, (10,30))

        m = scale.huber(X, axis=1)
        self.assertEquals(m.shape, (40,30))

        m = scale.huber(X, axis=2)
        self.assertEquals(m.shape, (40,10))

        m = scale.huber(X, axis=-1)
        self.assertEquals(m.shape, (40,10))

if __name__ == "__main__":
    unittest.main()
