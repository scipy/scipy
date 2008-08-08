""" Test module for interpolate1d.py
"""

# hack to test on Field's computer
import sys
sys.path.append('c:/home/python/Interpolate1d')

import numpy as np
from numpy import arange, meshgrid, ravel

from interpolate2d import interp2d, Interpolate2d, atleast_1d_and_contiguous
from fitpack_wrapper import Spline2d

# unit testing
import unittest, time
class Test(unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(atleast_1d_and_contiguous(x), atleast_1d_and_contiguous(y)))
        
    def test_callable_class_interpolation(self):
        """ make sure : instance of callable class with xyz not initiated works
        """
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        newx = np.arange(N-1)+0.5
        newy = newx
        
        interp_func = Interpolate2d(x, y, z, kind=Spline2d(kx=1, ky=1))
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx+newy)
        
    def test_no_out_of_range_args(self):
        """ make sure : having no out-of-range elements in new_x is fine
        """
        # There was a bug with this earlier.        
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        newx = arange(1,N-1)+.2
        newy = arange(1,N-1)+.3
        
        interp_func = Interpolate2d(x, y, z, kind='linear', out = 'linear')
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx + newy)
        
    def test_remove_bad_data(self):
        """make sure : works with bad data
        """
        N = 5.
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        x[2] = 55.0
        z[4] = np.NaN
        
        newx = arange(1, N-1)-0.5
        newy = newx
        
        newz = interp2d(x, y, z, newx, newy, kind='linear', bad_data = [55.0])
        self.assertAllclose(newz, newx+newy)
        
    def test_interp2d(self):
        """ make sure : interp2d works, at least in the linear case
        """
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        newx = arange(1, N)-0.5
        newy = newx
        
        newz = interp2d(x, y, z, newx, newy, kind='linear', out='linear')    
        
        self.assertAllclose(newz, newx+newy)
        
    def test_scalar_input(self):
        """ make sure : scalar input or a 0-degree array works
        """
        
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
                
        interp_func = Interpolate2d(x, y, z, kind='linear', out='linear')
        
        # scalar input
        newx1 = 0.5
        newy1 = 0.7
        
        newz1 = interp_func(newx1, newy1)
        self.assert_( np.isscalar(newz1) )
        
        # zero-degree array
        newx2 = 0.5
        newy2 = 0.7
        
        newz2 = interp_func(newx2, newy2)
        self.assert_( np.isscalar(newz2) )
        
    def test_extrapolation(self):
        """ make sure : default extrapolation works
        """
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        newx = np.arange(0, N+1)-0.5
        newy = newx
        
        interp_func = Interpolate2d(x, y, z, kind=Spline2d(kx=1, ky=1) )
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz[1:5], newx[1:5]+newy[1:5] )
        self.assert_( np.isnan(newz[0]) )
        self.assert_( np.isnan(newz[-1]) )
        
    def test_string_linear(self):
        """ make sure : string 'linear' works
        """
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        newx = np.arange(1, N)-0.5
        newy = newx
        
        interp_func = Interpolate2d(x, y, z, kind='linear', out='linear')
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx+newy)
        
    def test_string_quadratic(self):
        """ make sure : string 'quadratic' works
        """
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        newx = np.arange(1, N)-0.5
        newy = newx
        
        interp_func = Interpolate2d(x, y, z, kind='quadratic', out='quad')
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx+newy)
        
    def test_string_cubic(self):
        """make sure : string "cubic" works
        """
        N = 7
        X, Y = meshgrid(arange(N), arange(N))
        Z = X + Y
        x, y, z = map(ravel, [X, Y, Z] )
        
        newx = np.arange(1, N)-0.5
        newy = newx
        
        interp_func = Interpolate2d(x, y, z, kind='cubic', out='cubic')
        newz = interp_func(newx, newy)
        
        self.assertAllclose(newz, newx+newy)
        
if __name__ == '__main__':
    unittest.main()                 