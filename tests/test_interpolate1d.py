""" Test module for interpolate1d.py
"""

# hack to test on Field's computer
import sys
sys.path.append('c:/home/python/Interpolate1d')

import numpy as np
from numpy import arange

from interpolate1d import interp1d, Interpolate1d, atleast_1d_and_contiguous
from fitpack_wrapper import Spline

# unit testing
import unittest, time
class Test(unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(atleast_1d_and_contiguous(x), atleast_1d_and_contiguous(y)))
        
    def test_interpolate_wrapper(self):
        """ run unit test contained in interpolate_wrapper.py
        """
        #print "\n\nTESTING _interpolate_wrapper MODULE"
        from test_interpolate_wrapper import Test
        T = Test()
        T.runTest()
        
    def test_fitpack_wrapper(self):
        """ run unit test contained in fitpack_wrapper.py
        """
        #print "\n\nTESTING _fitpack_wrapper MODULE"
        from test_fitpack_wrapper import Test
        T = Test()
        T.runTest()
        
    def test_callable_class_interpolation(self):
        """ make sure : an instance of a callable class in which
            x and y haven't been initiated works
        """
        N = 7 #must be > 5
        x = np.arange(N)
        y = np.arange(N)
        interp_func = Interpolate1d(x, y, kind=Spline(k=2), low=Spline, high=Spline)
        new_x = np.arange(N+1)-0.5
        new_y = interp_func(new_x)
        self.assertAllclose(new_x, new_y)
        
    def test_no_out_of_range_args(self):
        """ make sure : having no out-of-range elements in new_x is fine
        """
        # There was a bug with this earlier.        
        N = 5
        x = arange(N)
        y = arange(N)
        new_x = arange(1,N-1)+.2
        interp_func = Interpolate1d(x, y, kind='linear', low='linear', high=np.NaN)
        new_y = interp_func(new_x)
        self.assertAllclose(new_x, new_y)
        
    def test_remove_bad_data(self):
        """make sure : interp1d works with bad data
        """
        N = 7.0 # must be >=5
        x = arange(N); x[2] = np.NaN
        y = arange(N); y[4] = 52.3; y[0]=np.NaN
        new_x = arange(N+1)-0.5
        new_y = interp1d(x, y, new_x, kind='linear', low='linear', 
                                high='linear', bad_data = [52.3])
        self.assertAllclose(new_x, new_y)
        
    def test_interp1d(self):
        """ make sure : interp1d works, at least in the linear case
        """
        N = 7
        x = arange(N)
        y = arange(N)
        new_x = arange(N+1)-0.5
        new_y = interp1d(x, y, new_x, kind='linear', low='linear', high='linear')        
        self.assertAllclose(new_x, new_y)
        
        
    def test_spline1_default_extrapolation(self):
        """ make sure : spline order 1 (linear) interpolation works correctly
            make sure : default extrapolation works
        """
        #print "\n\nTESTING LINEAR (1st ORDER) SPLINE"
        N = 7 # must be > 5
        x = np.arange(N)
        y = np.arange(N)
        interp_func = Interpolate1d(x, y, kind=Spline(k=1), low=None, high=599.73)
        new_x = np.arange(N+1)-0.5
        new_y = interp_func(new_x)
        
        self.assertAllclose(new_y[1:5], [0.5, 1.5, 2.5, 3.5])
        self.assert_(new_y[0] == None)
        self.assert_(new_y[-1] == 599.73)
        
    def test_spline2(self):
        """ make sure : order-2 splines work on linear data
            make sure : order-2 splines work on non-linear data
            make sure : 'cubic' and 'quad' as arguments yield
                                the desired spline
        """
        #print "\n\nTESTING 2nd ORDER SPLINE"
        N = 7 #must be > 5
        x = np.arange(N)
        y = np.arange(N)
        T1 = time.clock()
        interp_func = Interpolate1d(x, y, kind=Spline(k=2), low='spline', high='spline')
        T2 = time.clock()
        #print "time to create 2nd order spline interp function with N = %i: " % N, T2 - T1
        new_x = np.arange(N+1)-0.5
        t1 = time.clock()
        new_y = interp_func(new_x)
        t2 = time.clock()
        #print "time to evaluate 2nd order spline interp function with N = %i: " % N, t2 - t1
        self.assertAllclose(new_x, new_y)
        
        # make sure for non-linear data
        N = 7
        x = np.arange(N)
        y = x**2
        interp_func = Interpolate1d(x, y, kind=Spline(k=2), low='quad', high='cubic')
        new_x = np.arange(N+1)-0.5
        new_y = interp_func(new_x)
        self.assertAllclose(new_x**2, new_y)
        
        
    def test_linear(self):
        """ make sure : linear interpolation works 
            make sure : linear extrapolation works
        """
        #print "\n\nTESTING LINEAR INTERPOLATION"
        N = 7
        x = arange(N)
        y = arange(N)
        new_x = arange(N+1)-0.5
        T1 = time.clock()
        interp_func = Interpolate1d(x, y, kind='linear', low='linear', high='linear')
        T2 = time.clock()
        #print "time to create linear interp function with N = %i: " % N, T2 - T1
        t1 = time.clock()
        new_y = interp_func(new_x)
        t2 = time.clock()
        #print "time to create linear interp function with N = %i: " % N, t2 - t1
        
        self.assertAllclose(new_x, new_y)
        
    def test_scalar_input(self):
        """ make sure : newx being a scalar or a 0-degree array
                            makes the output a scalar
        """
        N = 7
        x = arange(N)
        y = arange(N)
        interp_func = Interpolate1d(x, y, kind='linear', low='linear', high='linear')
        
        # scalar input
        newx1 = 0.5
        newy1 = interp_func(newx1)
        self.assert_( np.isscalar(newy1) )
        
        # zero-degree array
        newx2 = np.array(0.5)
        newy2 = interp_func(newx2)
        self.assert_( np.isscalar(newy2) )
        
if __name__ == '__main__':
    unittest.main()                 