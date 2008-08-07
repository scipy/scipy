""" module for testing ndimage_wrapper
"""

# hack to test on Field's computer
import sys
sys.path.append('c:/home/python/Interpolate1d')

import unittest
import time
from numpy import arange, allclose, ones, array
import numpy as np
import ndimage_wrapper as nd

class Test (unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(x, y))
    
    def test_linear(self):
        """ Make sure : basic linear works
        """
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data, kind = 'linear')
        self.assertAllclose( interp(np.array([[2.3], [1.0], [3.9]])) , 1.0 )
        
    def test_linear_not_1(self):
        """ Make sure : linear interpolation works on a general dataset
        """
        X, Y = np.meshgrid(arange(10.), arange(10.))
        interesting_data = X+Y
        interp = nd.InterpolateNd(interesting_data, kind = 'linear')
        self.assertAllclose( interp(np.array([[2.3], [1.0]])) , 3.3 )
        
    def test_data_is_list(self):
        """ Make sure : data can be entered as a list
        """
        boring_data = [ [1.0, 1.0, 1.0],
                              [1.0, 1.0, 1.0],
                              [1.0, 1.0, 1.0]]
        interp = nd.InterpolateNd(boring_data)
        self.assertAllclose( interp(np.array([[1.3], [1.0]])) , 1.0 )
        
    def test_coords_is_1d(self):
        """ Make sure : coordinates for a single point can be entered as a 1D array
        """
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data)
        self.assertAllclose( interp(np.array([2.3, 1.0, 3.9])) , 1.0 )
        
    def test_coords_is_list(self):
        """ Make sure : coordinates for a single point can be entered as a list
        """
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data)
        self.assertAllclose( interp([2.3, 1.0, 3.9]) , 1.0 )
            
    def test_order2(self):
        """ Make sure : quadratic interpolation works
        """
        X, Y = np.meshgrid(arange(10.), arange(10.))
        interesting_data = X+Y
        interp = nd.InterpolateNd(interesting_data, kind = 2)
        self.assertAllclose( interp(np.array([[2.3], [1.0]])) , 3.3 )
        
    def test_order0(self):
        """ Make sure : block interpolation works
        """
        X, Y = np.meshgrid(arange(10.), arange(10.))
        interesting_data = X+Y
        interp = nd.InterpolateNd(interesting_data, kind = 0)
        self.assertAllclose( interp(np.array([[2.3], [1.1]])) , 3.0 )
        
    def test_order3(self):
        """ Make sure : cubic interpolation works
        """
        X, Y = np.meshgrid(arange(10.), arange(10.))
        interesting_data = X+Y
        interp = nd.InterpolateNd(interesting_data, kind = 3)
        self.assertAllclose( interp(np.array([[4.3], [4.1]])) , 8.4 )
        
    def test_out(self):
        """ Make sure : out-of-bounds returns NaN
        """
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data, kind = 'linear')
        self.assert_( np.isnan(interp(  np.array([[7.3], [1.0], [3.9]])  )))
        
    def test_starting_coords(self):
        """ Make sure : non-zero starting coordinates work correctly
        """
        X, Y = np.meshgrid(arange(10.), arange(10.))
        interesting_data = X+Y
        interp = nd.InterpolateNd(interesting_data, starting_coords = array([2, 1]))
        self.assertAllclose( interp(np.array([[2.3], [1.0]])) , 0.3 )
    
    def test_spacings(self):
        """ Make sure : spacings other than 1 work correctly
        """
        X, Y = np.meshgrid(arange(10.), arange(10.))
        interesting_data = X+Y
        interp = nd.InterpolateNd(interesting_data, spacings = array([2, 1]))
        self.assertAllclose( interp(np.array([[2.4], [1.0]])) , 2.2 )
        
    def runTest(self):
        """ run all tests
        """
        test_list = [method_name for method_name in dir(self) if method_name.find('test')==0]
        for test_name in test_list:
            exec("self.%s()" % test_name)
        
        
if __name__ == '__main__':
    unittest.main()