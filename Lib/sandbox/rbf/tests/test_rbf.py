#!/usr/bin/env python
# Created by John Travers, February 2007
""" Test functions for rbf module """

from numpy.testing import *
import numpy as n

set_package_path()
from rbf.rbf import Rbf
restore_path()

class test_Rbf1D(NumpyTestCase):
    def check_multiquadrics(self):
        x = n.linspace(0,10,9)
        y = n.sin(x) 
        rbf = Rbf(x, y)
        yi = rbf(x)
        assert_array_almost_equal(y.flatten(), yi)

class test_Rbf2D(NumpyTestCase):
    def check_multiquadrics(self):
        x = n.random.rand(50,1)*4-2
        y = n.random.rand(50,1)*4-2
        z = x*n.exp(-x**2-y**2)
        rbf = Rbf(n.c_[x.flatten(),y.flatten()].T,z.T,constant=2)
        zi = rbf(n.c_[x.flatten(), y.flatten()].T)
        zi.shape = x.shape
        assert_array_almost_equal(z, zi)

if __name__ == "__main__":
    NumpyTest().run()