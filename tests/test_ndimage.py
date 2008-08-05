""" module for testing ndimage_wrapper
"""

# hack to test on Field's computer
import sys
sys.path.append('c:/home/python/Interpolate1d')

import unittest
import time
from numpy import arange, allclose, ones
import numpy as np
import ndimage_wrapper as nd

class Test (unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(x, y))
    
    def test_linear(self):
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data, order = 1)
        self.assertAllclose( interp(np.array([[2.3], [1.0], [3.9]])) , 1.0 )
        
    def test_data_is_list(self):
        boring_data = [ [1.0, 1.0, 1.0],
                              [1.0, 1.0, 1.0],
                              [1.0, 1.0, 1.0]]
        interp = nd.InterpolateNd(boring_data, order = 1)
        self.assertAllclose( interp(np.array([[1.3], [1.0]])) , 1.0 )
        
    def test_coords_is_1d(self):
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data, order = 1)
        self.assertAllclose( interp(np.array([2.3, 1.0, 3.9])) , 1.0 )
        
    def test_coords_is_list(self):
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data, order = 1)
        self.assertAllclose( interp([2.3, 1.0, 3.9]) , 1.0 )
        
    def runTest(self):
        test_list = [method_name for method_name in dir(self) if method_name.find('test')==0]
        for test_name in test_list:
            exec("self.%s()" % test_name)
            
    def test_order2(self):
        boring_data = np.ones((5,5,5))
        interp = nd.InterpolateNd(boring_data, order = 2)
        self.assertAllclose( interp(np.array([[2.3], [1.0], [3.9]])) , 1.0 )
        
    def test_out(self):
        pass
        
        
if __name__ == '__main__':
    unittest.main()