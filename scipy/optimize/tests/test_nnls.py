""" Unit tests for nonnegative least squares
Author: Uwe Schmitt
Sep 2008
"""

from numpy.testing import *

from scipy.optimize import nnls
from numpy import arange, dot
from numpy.linalg import norm


class TestNNLS(TestCase):

    def test_nnls(self):
        a=arange(25.0).reshape(-1,5)
        x=arange(5.0)
        y=dot(a,x)
        x, res= nnls(a,y)
        assert res<1e-7
        assert norm(dot(a,x)-y)<1e-7

if __name__ == "__main__":
    run_module_suite()
