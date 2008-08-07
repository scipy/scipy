""" test fitpack_wrapper2d.py
"""

# hack to test on Field's computer
import sys
sys.path.append('c:/home/python/Interpolate1d')

# testing
import unittest
import time
from numpy import arange, allclose, ones, meshgrid, ravel, array
import numpy as np
from fitpack_wrapper import Spline2d

class Test(unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(x, y))
        
    def test_k_1(self):
        """ make sure : linear interpolation (spline with order = 1, s = 0)works
        """
        N = 10.
        X, Y = meshgrid( arange(N), arange(N) )
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z])
        
        newx = arange(N-1) +.5
        newy = newx
        
        interp_func = Spline2d(x, y, z, kx=1, ky=1)
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx+newy)
        
    def test_k_2(self):
        """ make sure : quadratic interpolation (spline with order = 2, s = 0)works
        """
        N = 10.
        X, Y = meshgrid( arange(N), arange(N) )
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z])
        
        newx = arange(N-1)+0.5
        newy = newx
        
        interp_func = Spline2d(x, y, z, kx=2, ky=2)
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx+newy)
        
        
    def test_instantiate_without_init(self):
        """ make sure : it's possible to instantiate Spline2d without setting x and y
        """
        N = 10.
        X, Y = meshgrid( arange(N), arange(N) )
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z])
        
        newx = np.arange(N-1)+0.5
        newy = newx
        
        interp_func = Spline2d(kx=1, ky=3)
        interp_func.init_xyz(x, y, z)
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx+newy)
        
    def test_extrap(self):
        """ make sure : linear extrapolation works
        """
        N = 10.
        X, Y = meshgrid( arange(N), arange(N) )
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z])
        
        interp_func = Spline2d(x, y, z, kx=1, ky=1)
        
        # upper-right region of R2
        self.assertAllclose(interp_func(np.array([N+1.]),np.array([N+1.])) , 2*N-2)
        
        # directly above interpolation region; only extrapolating in one variable
        self.assertAllclose(interp_func(np.array([N])/2.,2.*np.array([N])) , N/2. + (N-1.))
    
    def runTest(self):
        test_list = [name for name in dir(self) if name.find('test_')==0]
        for test_name in test_list:
            exec("self.%s()" % test_name)
        
        
if __name__ == '__main__':
    unittest.main()
    
    