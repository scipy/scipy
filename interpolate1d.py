# FIXME: information strings giving mathematical descriptions of the actions
#     of the functions.

from interpolate_wrapper import linear, logarithmic, block, block_average_above, atleast_1d_and_contiguous
from fitpack_wrapper import Spline
import numpy as np
from numpy import array, arange, empty, float64, NaN
    
def interp1d(x, y, new_x, interp = 'linear', low = NaN, high = NaN,
                    interpkw = {}, lowkw={}, highkw={},
                    bad_data = None):
    """ A function for interpolation of 1D data.
        
        Parameters
        -----------
            
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
                
        Optional Arguments
        -------------------
        
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
            
        Acceptable Input Strings
        ------------------------
        
            "linear" -- linear interpolation : default
            "logarithmic" -- logarithmic interpolation : linear in log space?
            "block" --
            "block_average_above' -- block average above
            "Spline" -- spline interpolation.  keyword k (defaults to 3) 
                indicates order of spline
            numpy.NaN -- return numpy.NaN
        
        Examples
        ---------
        
            >>> import numpy
            >>> from Interpolate1D import interp1d
            >>> x = range(5)        # note list is permitted
            >>> y = numpy.arange(5.)
            >>> new_x = [.2, 2.3, 5.6]
            >>> interp1d(x, y, new_x)
            array([.2, 2.3, 5.6, NaN])
    """
    return Interpolate1d(x, y, 
                                interp = interp,
                                low = low,
                                high = high,
                                interpkw = interpkw,
                                lowkw = lowkw,
                                highkw = highkw,
                                bad_data = bad_data
                                )(new_x)

class Interpolate1d(object):
    """ A callable class for interpolation of 1D, real-valued data.
        
        Parameters
        -----------
            
        x -- list or 1D NumPy array
            x includes the x-values for the data set to
            interpolate from.  It must be sorted in
            ascending order.
                
        y -- list or 1D NumPy array
            y includes the y-values for the data set  to
            interpolate from.  Note that 2-dimensional
            y is not supported.
                
        Optional Arguments
        -------------------
        
        kind -- Usu. string or function.  But can be any type.
            Specifies the type of interpolation to use for values within
            the range of x.
            
            If a string is passed, it will look for an object
            or function with that name and call it when evaluating.
            This is the primary mode of operation.  See below for list
            of acceptable strings.
            
            By default, linear interpolation is used.
            
            Other options are also available:
            
                If a callable class is passed, it is assumed to have format
                    instance = Class(x, y, **kw).
                It is instantiated and used for interpolation when the instance
                of Interpolate1d is called.
                
                If a callable object with method "init_xy" or "set_xy" is
                passed, that method will be used to set x and y, and the
                object will be called during interpolation.
                
                If a function is passed, it will be called when interpolating.
                It is assumed to have the form 
                    newy = kind(x, y, newx), 
                where x, y, newx, and newy are all numpy arrays.
                
                A primitive type which is not a string signifies a function
                which is identically that value (e.g. val and 
                lambda x, y, newx : val are equivalent).
            
        low  -- same as for kind
            How to extrapolate values for inputs below the range of x.
            Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
            a number') for all values below the range of x.
            
        high  -- same as for kind
            How to extrapolate values for inputs above the range of x.
            Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
            a number') for all values above the range of x.
            
        kindkw -- dictionary
            If kind is a class, function or string, additional keyword arguments
            may be needed (example: if you want a 2nd order spline, you could
            set kind = 'spline' and kindkw = {'k' : 2}.)
            
        lowkw -- like kindkw, but for low extrapolation
            
        highkw -- like kindkw, except for high extrapolation
        
        remove_bad_data -- bool
            indicates whether to remove bad data points from x and y.
            
        bad_data -- list
            List of values (in x or y) which indicate unacceptable data. All points
            that have x or y value in missing_data will be removed before
            any interpolatin is performed if remove_bad_data is true.
            
            numpy.NaN is always considered bad data.
            
        Some Acceptable Input Strings
        ------------------------
        
            "linear" -- linear interpolation : default
            "logarithmic" -- logarithmic interpolation : linear in log space?
            "block" --
            "block_average_above' -- block average above
            "Spline" -- spline interpolation.  keyword k (defaults to 3) 
                indicates order of spline
            "quad", "quadratic" -- spline interpolation order 2
            "cubic" -- spline interpolation order 3
            "quartic" -- spline interpolation order 4
            "quintic" -- spline interpolation order 5
        
        Examples
        ---------
        
            >>> import numpy
            >>> from interpolate1d import Interpolate1d
            >>> x = range(5)        # note list is permitted
            >>> y = numpy.arange(5.)
            >>> new_x = [.2, 2.3, 5.6]
            >>> interp_func = Interpolate1d(x, y)
            >>> interp_fuc(new_x)
            array([.2, 2.3, 5.6, NaN])
    """
    # FIXME: more informative descriptions of sample arguments
    # FIXME: examples in doc string
    # FIXME : Allow copying or not of arrays.  non-copy + remove_bad_data should flash 
    #           a warning (esp if we interpolate missing values), but work anyway.
    
    def __init__(self, x, y, interp = 'linear', low = NaN, high = NaN,
                        interpkw={}, lowkw={}, highkw={}, bad_data = None):
        # FIXME: don't allow copying multiple times.
        # FIXME : allow no copying, in case user has huge dataset
        
        # remove bad data, is there is any
        if bad_data is not None:
            x, y = self._remove_bad_data(x, y, bad_data)
        
        # check acceptable size and dimensions
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        assert len(x) > 0 and len(y) > 0 , "Arrays cannot be of zero length"
        assert x.ndim == 1 , "x must be one-dimensional"
        assert y.ndim == 1 , "y must be one-dimensional" 
        assert len(x) == len(y) , "x and y must be of the same length"
        
        # store data, and remove bad data points is applicable
        # FIXME : may be good to let x and y be initialized later, or changed after-the-fact
        self._init_xy(x, y)
        
        # store interpolation functions for each range
        self.interp = self._init_interp_method(interp, interpkw)
        self.low = self._init_interp_method(low, lowkw)
        self.high = self._init_interp_method(high, highkw)

    def _init_xy(self, x, y):
        
        # select proper dataypes and make arrays
        self._xdtype = {np.float32 : np.float32}.setdefault(type(x[0]), np.float64) # unless data is float32,  cast to float64
        self._ydtype = {np.float32 : np.float32}.setdefault(type(y[0]), np.float64)
        self._x = atleast_1d_and_contiguous(x, self._xdtype).copy()
        self._y = atleast_1d_and_contiguous(y, self._ydtype).copy()

    def _remove_bad_data(self, x, y, bad_data = [None, NaN]):
        """ removes data points whose x or y coordinate is
            either in bad_data or is a NaN.
        """
        # FIXME : In the future, it may be good to just replace the bad points with good guesses.
        #       Especially in generalizing the higher dimensions
        # FIXME : This step is very inefficient because it iterates over the array
        mask = np.array([  (xi not in bad_data) and (not np.isnan(xi)) and \
                                    (y[i] not in bad_data) and (not np.isnan(y[i])) \
                                for i, xi in enumerate(x) ])
        x = x[mask]
        y = y[mask]
        return x, y
        
    def _init_interp_method(self, interp_arg, kw):
        """
            User provides interp_arg and dictionary kw.  _init_interp_method
            returns the interpolating function specified by interp_arg,
            possibly with extra keyword arguments given in kw.
        
        """
        # FIXME : error checking specific to interpolation method.  x and y long
        #   enough for order-3 spline if that's indicated, etc.  Functions should throw
        #   errors themselves, but errors at instantiation would be nice.
        
        from inspect import isclass, isfunction
        
        # primary usage : user passes a string indicating a known function
        if interp_arg in ['linear', 'logarithmic', 'block', 'block_average_above']:
            # string used to indicate interpolation method,  Select appropriate function
            func = {'linear':linear, 'logarithmic':logarithmic, 'block':block, \
                        'block_average_above':block_average_above}[interp_arg]
            result = lambda new_x : func(self._x, self._y, new_x, **kw)
        elif interp_arg in ['Spline', 'spline']:
            # use the Spline class from fitpack_wrapper
            # k = 3 unless otherwise specified
            result = Spline(self._x, self._y, **kw)
        elif interp_arg in ['Quadratic', 'quadratic', 'Quad', 'quad', \
                                'Cubic', 'cubic', \
                                'Quartic', 'quartic', 'Quar', 'quar',\
                                'Quintic', 'quintic', 'Quin', 'quin']:
            # specify specific kinds of splines
            if interp_arg in ['Quadratic', 'quadratic', 'Quad', 'quad']:
                result = Spline(self._x, self._y, k=2)
            elif interp_arg in ['Cubic', 'cubic']:
                result = Spline(self._x, self._y, k=3)
            elif interp_arg in ['Quartic', 'quartic', 'Quar', 'quar']:
                result = Spline(self._x, self._y, k=4)
            elif interp_arg in ['Quintic', 'quintic', 'Quin', 'quin']:
                result = Spline(self._x, self._y, k=5)
        elif isinstance(interp_arg, basestring):
            raise TypeError, "input string %s not valid" % interp_arg
        
        # secondary usage : user passes a callable class
        elif isclass(interp_arg) and hasattr(interp_arg, '__call__'):
            if hasattr(interp_arg, 'init_xy'):
                result = interp_arg(**kw)
                result.init_xy(self._x, self._y)
            elif hasattr(interp_arg, 'set_xy'):
                result = interp_arg(**kw)
                result.set_xy(self._x, self._y)
            else:
                result = interp_arg(x, y, **kw)
                
        # user passes an instance of a callable class which has yet
        # to have its x and y initialized.
        elif hasattr(interp_arg, 'init_xy') and hasattr(interp_arg, '__call__'):
            result = interp_arg
            result.init_xy(self._x, self._y)
        elif hasattr(interp_arg, 'set_xy') and hasattr(interp_arg, '__call__'):
            result = interp_arg
            result.set_xy(self._x, self._y)
                
        # user passes a function to be called
        # Assume function has form of f(x, y, newx, **kw)
        elif isfunction(interp_arg):
            result = lambda new_x : interp_arg(self._x, self._y, new_x, **kw)
        
        # default : user has passed a default value to always be returned
        else:
            result = np.vectorize(lambda new_x : interp_arg)
            
        return result

    def __call__(self, newx):
        """
            Input x must be a list or NumPy array in sorted order.
            
            Breaks x into pieces in-range, below-range, and above range.
            Performs appropriate operation on each and concatenates results.
        """
        # FIXME : atleast_1d_and_contiguous may also be called within the interpolation technique.
        #   waste of time, but ok for the time being.
        
        # if input is scalar or 0-dimemsional array, output will be scalar
        input_is_scalar = np.isscalar(newx) or (isinstance(newx, type(np.array([1.0]))) and np.shape(newx) == ())
        
        newx_array = atleast_1d_and_contiguous(newx)
        
        # masks indicate which elements fall into which interpolation region
        low_mask = newx_array<self._x[0]
        high_mask = newx_array>self._x[-1]
        interp_mask = (~low_mask) & (~high_mask)
                
        # use correct function for x values in each region
        if len(newx_array[low_mask]) == 0: new_low=np.array([])  # FIXME : remove need for if/else.
                                                                            # if/else is a hack, since vectorize is failing
                                                                            # to work on lists/arrays of length 0
                                                                            # on the computer where this is being
                                                                            # developed
        else: new_low = self.low(newx_array[low_mask])
        if len(newx_array[interp_mask])==0: new_interp=np.array([])
        else: new_interp = self.interp(newx_array[interp_mask])
        if len(newx_array[high_mask]) == 0: new_high = np.array([])
        else: new_high = self.high(newx_array[high_mask])
        
        result_array = np.concatenate((new_low, new_interp, new_high)) # FIXME : deal with mixed datatypes
                                                                                          # Would be nice to say result = zeros(dtype=?)
                                                                                          # and fill in
        
        if input_is_scalar:
            result = float(result_array)
        else:
            result = result_array
        
        return result
  