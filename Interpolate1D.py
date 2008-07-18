"""
    Interpolation of 1D data

    This module provides several functions and classes for interpolation
    and extrapolation of 1D data (1D in both input and output).  The
    primary function provided is:

        interp1d(x, y, new_x) : from data points x and y, interpolates
                                        values for points in new_x and
                                        returns them as an array.

    Classes provided include:

        Interpolate1D  :   an object for interpolation of
                                various kinds.  interp1d is a wrapper
                                around this class.
                                
        Spline : an object for spline interpolation
        
    Functions provided include:

        linear : linear interpolation
        logarithmic :  logarithmic interpolation
        block : block interpolation
        block_average_above : block average above interpolation

"""

# FIXME: information strings giving mathematical descriptions of the actions
#     of the functions.

from interpolate_wrapper import linear, logarithmic, block, block_average_above
from fitpack_wrapper import Spline
import numpy as np
from numpy import array, arange, empty, float64, NaN

# FIXME: use this to ensure proper type of all inputs and outputs in Interpolate1D
def make_array_safe(ary, typecode=np.float64):
    # FIXME: could pick correct typecode
    ary = np.atleast_1d(np.asarray(ary, typecode))
    if not ary.flags['CONTIGUOUS']:
        ary = ary.copy()
    return ary
    
def interp1d(x, y, new_x, kind='linear', low=np.NaN, high=np.NaN, kindkw={}, lowkw={}, highkw={}, \
                        remove_bad_data = False, bad_data=[]):
    """ A function for interpolation of 1D data.
        
        REQUIRED ARGUMENTS:
            
        x -- list or NumPy array
            x includes the x-values for the data set to
            interpolate from.  It must be sorted in
            ascending order
                
        y -- list or NumPy array
            y includes the y-values for the data set  to
            interpolate from.  Note that y must be
            one-dimensional.
            
        new_x -- list or NumPy array
            points whose value is to be interpolated from x and y.
            new_x must be in sorted order, lowest to highest.
                
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
        
        remove_bad_data -- bool
            indicates whether to remove bad data.
            
        bad_data -- list
            List of values (in x or y) which indicate unacceptable data. All points
            that have x or y value in missing_data will be removed before
            any interpolatin is performed if remove_bad_data is true.
            
            numpy.NaN is always considered bad data.
            
        SAMPLE ACCEPTABLE ARGUMENTS:
            "linear" -- linear interpolation : default
            "logarithmic" -- logarithmic interpolation : linear in log space?
            "block" --
            "block_average_above' -- block average above
            "Spline" -- spline interpolation.  keyword k (defaults to 3) 
                indicates order of spline
            numpy.NaN -- return numpy.NaN
        
        EXAMPLES:
            >>> import numpy
            >>> from Interpolate1D import interp1d
            >>> x = range(5)        # note list is permitted
            >>> y = numpy.arange(5.)
            >>> new_x = [.2, 2.3, 5.6]
            >>> interp1d(x, y, new_x)
            array([.2, 2.3, 5.6, NaN])
    """
    return Interpolate1D(x, y, kind=kind, low=low, high=high, kindkw=kindkw, lowkw=lowkw, highkw=highkw, \
                        remove_bad_data = remove_bad_data, bad_data=bad_data)(new_x)

class Interpolate1D(object):
    """ An object for interpolation of 1D data.
        
    REQUIRED ARGUMENTS:
        
    x -- list or NumPy array
        x includes the x-values for the data set to
        interpolate from.  It must be sorted in
        ascending order
            
    y -- list or NumPy array
        y includes the y-values for the data set  to
        interpolate from.  Note that y must be
        one-dimensional.
            
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
    
    remove_bad_data -- bool
        indicates whether to remove bad data.
        
    bad_data -- list
        List of values (in x or y) which indicate unacceptable data. All points
        that have x or y value in missing_data will be removed before
        any interpolatin is performed if remove_bad_data is true.
        
        numpy.NaN is always considered bad data.
        
    SAMPLE ACCEPTABLE ARGUMENTS:
        "linear" -- linear interpolation : default
        "logarithmic" -- logarithmic interpolation : linear in log space?
        "block" --
        "block_average_above' -- block average above
        "Spline" -- spline interpolation.  keyword k (defaults to 3) 
            indicates order of spline
        numpy.NaN -- return numpy.NaN
    
    EXAMPLES:
        >>> import numpy
        >>> from Interpolate1D import Interpolate1D
        >>> x = range(5)        # note list is permitted
        >>> y = numpy.arange(5.)
        >>> interp = Interpolate1D(x, y)
        >>> new_x = [.2, 2.3, 5.6]
        >>> interp(new_x)
        array([.2, 2.3, 5.6, NaN])
        
    """
    # FIXME: more informative descriptions of sample arguments
    # FIXME: examples in doc string
    
    def __init__(self, x, y, kind='linear', low=np.NaN, high=np.NaN, kindkw={}, lowkw={}, highkw={}, \
                        remove_bad_data = False, bad_data=[]):
        
        #print "size of bad_data: ", len(bad_data)
        
        self._format_array(x, y, remove_bad_data = remove_bad_data, bad_data = bad_data)
        # FIXME: Handle checking if they are the correct size.
        #self._x = make_array_safe(x).copy()
        #self._y = make_array_safe(y).copy()
        
        #assert len(x) == len(y) , "x and y must be of the same length"
        #assert x.ndim == 1 , "x must be one-dimensional"
        #assert y.ndim == 1 , "y must be one-dimensional"
        # FIXME: let y be 2-dimensional.  Involves reworking of Interpolate1D.__call__
        #   because Spline enumerates y along the last, rather then first, axis,
        #   while concatenate works along first axis
        
        self.kind = self._init_interp_method(self._x, self._y, kind, kindkw)
        self.low = self._init_interp_method(self._x, self._y, low, lowkw)
        self.high = self._init_interp_method(self._x, self._y, high, highkw)

    def _format_array(self, x, y, remove_bad_data = False, bad_data = []):#=[None, np.NaN]):#=[None, np.NaN]):
        """
        Assigns properly formatted versions of x and y to self._x and self._y.
        Also records data types.
        
        Formatting includes removal of all points whose x or y coordinate
        is in missing_data.  This is the primary difference from
        make_array_safe.
        
        """
        # FIXME: don't allow copying multiple times.
         
        # check acceptable lengths for x and y
        assert len(x) > 0 and len(y) > 0 , "Interpolate1D does not support\
                                        arrays of length 0"
        assert len(x) == len(y) , "x and y must be of the same length"
        x = np.array(x)
        y = np.array(y)
        
        # remove bad data
        if remove_bad_data:
            mask = np.array([  (xi not in bad_data) and (not np.isnan(xi)) and (y[i] not in bad_data) and (not np.isnan(y[i])) \
                    for i, xi in enumerate(x) ])
            print 'mask equals: ', mask, type(mask)
            print 'x equals: ', x
            print 'x[mask] is: ', x[mask]
            x = x[mask]
            y = y[mask]
            
        # collect dataypes and make arrays
        self._xdtype = {np.float32 : np.float32}.setdefault(type(x[0]), np.float64) # unless data is float32,  cast to float64
        self._ydtype = {np.float32 : np.float32}.setdefault(type(y[0]), np.float64)
        self._x = make_array_safe(x, self._xdtype).copy()
        self._y = make_array_safe(y, self._ydtype).copy()
            
        # check dimensionality
        assert self._x.ndim == 1 , "x must be one-dimensional"
        assert self._y.ndim == 1 , "y must be one-dimensional"    
    
    def _init_interp_method(self, x, y, interp_arg, kw):
        """
        User provides interp_arg and dictionary kw.  _init_interp_method
        returns the interpolating function from x and y specified by interp_arg,
        possibly with extra keyword arguments given in kw.
        
        """
        # FIXME : error checking specific to interpolation method.  x and y long
        #   enough for order-3 spline if that's indicated, etc.  Functions should throw
        #   errors themselves when Interpolate1D is called, but errors at instantiation
        #   would be nice.
        
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
        """
            
        """
        
        x = make_array_safe(x)
        low_mask = x<self._x[0]
        high_mask = x>self._x[-1]
        interp_mask = (~low_mask) & (~high_mask)
        
        if len(x[low_mask]) == 0: new_low=np.array([]) # hack, since vectorize is failing
                                                                            # work on lists/arrays of length 0
                                                                            # on computer where this is being
                                                                            # developed
        else: new_low = self.low(x[low_mask])
        if len(x[interp_mask])==0: new_interp=np.array([])
        else: new_interp = self.kind(x[interp_mask])
        if len(x[high_mask]) == 0: new_high = np.array([])
        else: new_high = self.high(x[high_mask])
        
        result = np.concatenate((new_low, new_interp, new_high)) # FIXME : deal with mixed datatypes
        
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
        # make sure : having no out-of-range elements in new_x is fine
        # There was a bug with this earlier.
        N = 5
        x = arange(N)
        y = arange(N)
        new_x = arange(1,N-1)+.2
        interp_func = Interpolate1D(x, y, kind='linear', low='linear', high=np.NaN)
        new_y = interp_func(new_x)
        self.assertAllclose(new_x, new_y)
        
    def test_intper1d(self):
        # make sure : interp1d works, at least in the linear case
        N = 7
        x = arange(N)
        y = arange(N)
        new_x = arange(N+1)-0.5
        new_y = interp1d(x, y, new_x, kind='linear', low='linear', high='linear')        
        self.assertAllclose(new_x, new_y)
        
if __name__ == '__main__':
    unittest.main()                 