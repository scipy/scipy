""" module to test fitpack_wrapper.py
"""

# hack to test on Field's computer
import sys
sys.path.append('c:/home/python/Interpolate1d')

# testing
import unittest
import time
from numpy import arange, allclose, ones
import numpy as np
from fitpack_wrapper import Spline

class Test(unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(x, y))
        
    def test_k_1(self):
        """ make sure : linear interpolation (spline with order = 1, s = 0)works
        """
        N = 3000.
        x = np.arange(N)
        y = np.arange(N)
        #T1 = time.clock()
        interp_func = Spline(x, y, k=1)
        #T2 = time.clock()
        #print "time to create order 1 spline interpolation function with N = %i:" % N, T2 - T1
        new_x = np.arange(N)+0.5
        #t1 = time.clock()
        new_y = interp_func(new_x)
        #t2 = time.clock()
        #print "time for order 1 spline interpolation with N = %i:" % N, t2 - t1
        self.assertAllclose(new_y[:5], [0.5, 1.5, 2.5, 3.5, 4.5])
        
    def test_k_2(self):
        """ make sure : quadratic interpolation (spline with order = 2, s = 0)works
        """
        N = 3000.
        x = np.arange(N)
        y = x**2
        interp_func = Spline(x, y, k=2)
        #print "time to create order 1 spline interpolation function with N = %i:" % N, T2 - T1
        new_x = np.arange(N)+0.5
        #t1 = time.clock()
        new_y = interp_func(x)
        #t2 = time.clock()
        #print "time for order 1 spline interpolation with N = %i:" % N, t2 - t1
        self.assertAllclose(new_y, y)
        
    def test_extrap(self):
        """ make sure 1D extrapolation works
        """
        N = 3000.
        x = np.arange(N)
        y = np.arange(N)
        
        interp_func = Spline(x, y, k=1)
        newx = np.arange(-2, N)+0.5
        newy = interp_func(newx)
        
        self.assertAllclose(newx, newy)
        
    def test_inputFormat(self):
        """ make sure : it's possible to instantiate Spline without x and y
        """
        #print "testing input format"
        N = 3000.
        x = np.arange(N)
        y = np.arange(N)
        interp_func = Spline(k=1)
        interp_func.init_xy(x, y)
        new_x = np.arange(N)+0.5
        new_y = interp_func(new_x)
        self.assertAllclose(new_y[:5], [0.5, 1.5, 2.5, 3.5, 4.5])
    
    def runTest(self):
        test_list = [name for name in dir(self) if name.find('test_')==0]
        for test_name in test_list:
            exec("self.%s()" % test_name)
           
        
                             
if __name__ == '__main__':
    unittest.main()
    
    