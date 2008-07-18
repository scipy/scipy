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
    # fixme: could pick correct typecode
    ary = np.atleast_1d(np.asarray(ary, typecode))
    if not ary.flags['CONTIGUOUS']:
        ary = ary.copy()
    return ary
    

class Interpolate1D(object):
    
    
    def __init__(self, x, y, kind='linear', low=np.NaN, high=np.NaN, kindkw={}, lowkw={}, highkw={}, missing_data=[None, np.NaN]):
        """
        Object for interpolation of 1D data.
        
        REQUIRED ARGUMENTS:
        
        x -- list or NumPy array
            x includes the x-values for the data set to
            interpolate from.  It must be sorted in
            ascending order
            
        y -- list or NumPy array
            y includes the y-values for the data set  to
            interpolate from.
            
        OPTIONAL ARGUMENTS:
        
        kind -- Usu. function or string.  But can be any type.
            Specifies the type of extrapolation to use for values within
            the range of x.  If a string is passed, it will look for an object
            or function with that name and call it when evaluating.  If 
            a function or object is passed, it will be called when interpolating.
            If nothing else, assumes the argument is intended as a value
            to be returned for all arguments.  Defaults to linear interpolation.
            
        kindkw -- dictionary
            If kind is a class, function or string, additional keyword arguments
            may be needed (example: if you want a 2nd order spline, kind = 'spline'
            and kindkw = {'k' : 2}.
            
        low (high) -- same as for kind
            Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
            a number') for all values outside the range of x.
            
        
    
        """
        # fixme: Handle checking if they are the correct size.
        self._x = make_array_safe(x).copy()
        self._y = make_array_safe(y).copy()
        
        assert len(x) == len(y) , "x and y must be of the same length"
        assert x.ndim == 1 , "x must be one-dimensional"
        assert y.ndim == 1 , "y must be one-dimensional"
        # fixme: let y be 2-dimensional.  Involves reworking of Interpolate1D.__call__
        # because Spline enumerates y along the last, rather then first, axis,
        # while concatenate works along first axis
        
        self.kind = self._init_interp_method(self._x, self._y, kind, kindkw)
        self.low = self._init_interp_method(self._x, self._y, low, lowkw)
        self.high = self._init_interp_method(self._x, self._y, high, highkw)

    def _format_array(x, y, missing_data=[None, np.NaN]):
        # fixme: don't allow copying multiple times.
                        
        assert len(x) > 0 and len(y) > 0 , "interpolation does not support\
                                        array of length 0"
        assert len(x) == len(y) , "x and y must be of the same length"
        mask = [((xi not in missing_data) and (y[i] not in missing_data)) \
                    for i, xi in enumerate(x) ]
        if isinstance(x, list): 
            x = [x[i] for (i, good_data) in enumerate(mask) if good_data]
        else: 
            x = x[mask]
        if isinstance(y, list): 
            y = [y[i] for (i, good_data) in enumerate(mask) if good_data]
        else: 
            y = y[mask]
        self._xdtype = type(x[0])
        self._x = make_array_safe(x, _xdtype).copy()
        self._ydtype = type(y[0])
        self._y = make_array_safe(y, _ydtype).copy()
            
        assert self._x.ndim == 1 , "x must be one-dimensional"
        assert self._y.ndim == 1 , "y must be one-dimensional"    
    
    def _init_interp_method(self, x, y, interp_arg, kw):
        from inspect import isclass, isfunction
        
        if interp_arg in ['linear', 'logarithmic', 'block', 'block_average_above']:
            func = {'linear':linear, 'logarithmic':logarithmic, 'block':block, \
                        'block_average_above':block_average_above}[interp_arg]
            result = lambda new_x : func(self._x, self._y, new_x, **kw)
        elif interp_arg in ['Spline', Spline, 'spline']:
            result = Spline(self._x, self._y, **kw)
        elif isfunction(interp_arg):
            result = lambda new_x : interp_arg(new_x, **kw)
        elif isclass(interp_arg):
            result = interp_arg(x, y, **kw)
        else:
            result = np.vectorize(lambda new_x : interp_arg)
        return result

    def __call__(self, x):
        
        x = make_array_safe(x)
        low_mask = x<self._x[0]
        high_mask = x>self._x[-1]
        interp_mask = (~low_mask) & (~high_mask)
        
        if len(x[low_mask]) == 0: new_low=np.array([]) # hack, since vectorize is failing
                                                                            # work on lists/arrays of length 0
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
        
    def test__interpolate_wrapper(self):
        print "\n\nTESTING _interpolate_wrapper MODULE"
        from interpolate_wrapper import Test
        T = Test()
        T.runTest()
        
    def test__fitpack_wrapper(self):
        print "\n\nTESTING _fitpack_wrapper MODULE"
        from fitpack_wrapper import Test
        T = Test()
        T.runTest()
        
    def test_spline1_defaultExt(self):
        # make sure : spline order 1 (linear) interpolation works correctly
        # make sure : default extrapolation works
        print "\n\nTESTING LINEAR (1st ORDER) SPLINE"
        N = 7 # must be > 5
        x = np.arange(N)
        y = np.arange(N)
        interp_func = Interpolate1D(x, y, kind='Spline', kindkw={'k':1}, low=None, high=599.73)
        new_x = np.arange(N+1)-0.5
        new_y = interp_func(new_x)
        
        self.assertAllclose(new_y[1:5], [0.5, 1.5, 2.5, 3.5])
        self.assert_(new_y[0] == None)
        self.assert_(new_y[-1] == 599.73)
        
    def test_spline2(self):
        print "\n\nTESTING 2nd ORDER SPLINE"
        # make sure : order-2 splines work on linear data
        N = 7 #must be > 5
        x = np.arange(N)
        y = np.arange(N)
        T1 = time.clock()
        interp_func = Interpolate1D(x, y, kind='Spline', kindkw={'k':2}, low='spline', high='spline')
        T2 = time.clock()
        print "time to create 2nd order spline interp function with N = %i: " % N, T2 - T1
        new_x = np.arange(N+1)-0.5
        t1 = time.clock()
        new_y = interp_func(new_x)
        t2 = time.clock()
        print "time to evaluate 2nd order spline interp function with N = %i: " % N, t2 - t1
        self.assertAllclose(new_x, new_y)
        
        # make sure for non-linear data
        N = 7
        x = np.arange(N)
        y = x**2
        interp_func = Interpolate1D(x, y, kind='Spline', kindkw={'k':2}, low='spline', high='spline')
        new_x = np.arange(N+1)-0.5
        new_y = interp_func(new_x)
        self.assertAllclose(new_x**2, new_y)
        
        
    def test_linear(self):
        # make sure : linear interpolation works 
        # make sure : linear extrapolation works
        print "\n\nTESTING LINEAR INTERPOLATION"
        N = 7
        x = arange(N)
        y = arange(N)
        new_x = arange(N+1)-0.5
        T1 = time.clock()
        interp_func = Interpolate1D(x, y, kind='linear', low='linear', high='linear')
        T2 = time.clock()
        print "time to create linear interp function with N = %i: " % N, T2 - T1
        t1 = time.clock()
        new_y = interp_func(new_x)
        t2 = time.clock()
        print "time to create linear interp function with N = %i: " % N, t2 - t1
        
        self.assertAllclose(new_x, new_y)
        
    def test_noLow(self):
        # make sure : having the out-of-range elements in new_x is fine
        # there was a bug with this
        N = 5
        x = arange(N)
        y = arange(N)
        new_x = arange(1,N-1)+.2
        interp_func = Interpolate1D(x, y, kind='linear', low='linear', high=np.NaN)
        new_y = interp_func(new_x)
        self.assertAllclose(new_x, new_y)
        
if __name__ == '__main__':
    unittest.main()                 