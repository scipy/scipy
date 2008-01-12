"""
Test functions for models.bspline
"""

import numpy as N
from scipy.testing import *

import scipy.stats.models as S
import scipy.stats.models.bspline as B


class TestBSpline(TestCase):

    def test1(self):
        b = B.BSpline(N.linspace(0,10,11), x=N.linspace(0,10,101))
        old = b._basisx.shape
        b.x = N.linspace(0,10,51)
        new = b._basisx.shape
        self.assertEqual((old[0], 51), new)


if __name__ == "__main__":
    unittest.main()
