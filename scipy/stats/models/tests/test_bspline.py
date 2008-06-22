"""
Test functions for models.bspline
"""

import numpy as N
from scipy.testing import *

import scipy.stats.models as S
try:
    import scipy.stats.models.bspline as B
except ImportError:
    B = None


class TestBSpline(TestCase):

    def test1(self):
        if B:
            b = B.BSpline(N.linspace(0,10,11), x=N.linspace(0,10,101))
            old = b._basisx.shape
            b.x = N.linspace(0,10,51)
            new = b._basisx.shape
            self.assertEqual((old[0], 51), new)


if __name__ == "__main__":
    nose.run(argv=['', __file__])
