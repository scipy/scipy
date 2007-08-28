"""
Test functions for models.glm
"""

import numpy as N
from numpy.testing import NumpyTest, NumpyTestCase

import scipy.sandbox.models as S
import scipy.sandbox.models.bspline as B


class test_BSpline(NumpyTestCase):

    def test1(self):
        b = B.BSpline(N.linspace(0,10,11), x=N.linspace(0,10,101))
        old = b._basisx.shape
        b.x = N.linspace(0,10,51)
        new = b._basisx.shape
        self.assertEqual((old[0], 51), new)


if __name__ == "__main__":
    NumpyTest().run()
