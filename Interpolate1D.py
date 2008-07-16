""" A module for intepolation
"""

# fixme: information strings giving mathematical descriptions of the actions
#     of the functions.

from interpolate_wrapper import linear, logarithmic, block, block_average_above
from fitpack_wrapper import Spline
import numpy as np
from numpy import array, arange, empty, float64, NaN

# fixme: use this to ensure proper type of all inputs and outputs in Interpolate1D
def make_array_safe(ary, typecode=np.float64):
    ary = np.atleast_1d(np.asarray(ary, typecode))
    if not ary.flags['CONTIGUOUS']:
        ary = ary.copy()
    return ary
    

class Interpolate1D(object):
    # see enthought.interpolate
    
    # fixme: Handle other data types.
    
    def __init__(self, x, y, k=1, kind='linear', low=None, high=None):

        # fixme: Handle checking if they are the correct size.
        self._x = make_array_safe(x).copy()
        self._y = make_array_safe(y).copy()
        self._xdtype = type(self._x[0])
        self._ydtype = type(self._y[0])

        assert( len(x) == len(y) , "x and y must be of the same length" )
        assert( x.ndim == 1 , "x must be one-dimensional" )
        assert( y.ndim == 1 , "y must be one-dimensional" )
        # fixme: let y be 2-dimensional.  Involves reworking of Interpolate1D.__call__
        # because Spline enumerates y along the last, rather then first, axis,
        # while concatenate works along first axis
        
        self.kind = self._init_interp_method(self._x, self._y, k, kind)
        self.low = self._init_interp_method(self._x, self._y, k, low)
        self.high = self._init_interp_method(self._x, self._y, k, high)

    def _init_interp_method(self, x, y, k, interp_arg):
        from inspect import isclass, isfunction
        
        if interp_arg in ['linear', 'logarithmic', 'block', 'block_average_above']:
            func = {'linear':linear, 'logarithmic':logarithmic, 'block':block, \
                        'block_average_above':block_average_above}[interp_arg]
            result = lambda new_x : func(self._x, self._y, new_x)
        elif interp_arg in ['Spline', Spline, 'spline']:
            result = Spline(self._x, self._y, k=k)
        elif isfunction(interp_arg):
            result = interp_arg
        elif isclass(interp_arg):
            result = interp_arg(x, y)
        else:
            print "warning: defaulting on extrapolation"
            result = np.vectorize(lambda new_x : interp_arg)
        return result

    def __call__(self, x):
        
        x = make_array_safe(x)
        low_mask = x<self._x[0]
        high_mask = x>self._x[-1]
        interp_mask = (~low_mask) & (~high_mask)

        # hack, since getting an error when self.low or self.high gets 0-length array
        # and they return None or NaN
        if len(x[low_mask]) == 0: new_low=np.array([])
        else: new_low = self.low(x[low_mask])
        if len(x[interp_mask])==0: new_interp=np.array([])
        else: new_interp = self.kind(x[interp_mask])
        if len(x[high_mask]) == 0: new_high = np.array([])
        else: new_high = self.high(x[high_mask])
        
        result = np.concatenate((new_low, new_interp, new_high))
        
        return result
        
# unit testing
import unittest, time
class Test(unittest.TestCase):
    
    def assertAllclose(self, x, y):
        self.assert_(np.allclose(make_array_safe(x), make_array_safe(y)))
        
    # fixme: run the test contained in the wrapper modules
        
    def test_Interp_linearSpl(self):
        #return
        print "\n\nTESTING LINEAR (1st ORDER) SPLINE"
        N = 7
        x = np.arange(N)
        y = np.arange(N)
        T1 = time.clock()
        interp_func = Interpolate1D(x, y, k=1, kind='Spline', low=None, high=None)
        T2 = time.clock()
        print 'time to create linear interp function: ', T2 - T1
        new_x = np.arange(N)-0.5
        t1 = time.clock()
        new_y = interp_func(new_x)
        t2 = time.clock()
        print '1d interp (sec):', t2 - t1
        
        print "new_y: ", new_y
        self.assertAllclose(new_y[1:5], [0.5, 1.5, 2.5, 3.5])
        self.assert_(new_y[0] == None) 
        
    def test_linear(self):
        print "\n\nTESTING LINEAR INTERPOLATION"
        N = 7
        x = arange(N)
        y = arange(N)
        new_x = arange(N+1)-0.5
        T1 = time.clock()
        interp_func = Interpolate1D(x, y, kind='linear', low=None, high=None)
        T2 = time.clock()
        print 'time to create linear interp function: ', T2 - T1
        t1 = time.clock()
        new_y = interp_func(new_x)
        t2 = time.clock()
        print '1d interp (sec):', t2 - t1
        
        self.assertAllclose(new_y[1:5], [0.5, 1.5, 2.5, 3.5])
        self.assert_(new_y[0] == None)
        
if __name__ == '__main__':
    unittest.main()                 