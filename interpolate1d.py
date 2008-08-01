# FIXME: information strings giving mathematical descriptions of the actions
#     of the functions.

from interpolate_wrapper import atleast_1d_and_contiguous, \
                linear, logarithmic, block, block_average_above, nearest
from fitpack_wrapper import Spline
import numpy as np
from numpy import array, arange, empty, float64, NaN

# dictionary of interpolation functions/classes/objects
method_register = \
                { # functions
                    'linear' : linear,  'Linear' : linear, 
                    'logarithmic' : logarithmic, 'Logarithmic' : logarithmic, 
                    'block' : block, 'Block' : block, 
                    'block_average_above' : block_average_above, 
                    'Block_average_above' : block_average_above, 
                    'nearest' : nearest, 'Nearest' : nearest,
                    
                    # Splines
                    'Spline' : Spline, 'spline' : Spline,
                    'Quadratic' : Spline(k=2), 'quadratic' : Spline(k=2),
                    'Quad' : Spline(k=2), 'quad' : Spline(k=2),
                    'Cubic' : Spline(k=3), 'cubic' : Spline(k=3),
                    'Quartic' : Spline(k=4), 'quartic' : Spline(k=4),
                    'Quar' : Spline(k=4), 'quar' : Spline(k=4),
                    'Quintic' : Spline(k=5), 'quintic' : Spline(k=5),
                    'Quin' : Spline(k=5), 'quin' : Spline(k=5)
                }
                
# dictionary of types for casting.  key = possible datatype, value = datatype it is cast to
# BEWARE : if you cast things to integers, you will lose interpolation ability
dtype_register = {np.float32 : np.float32, 
                            np.float64 : np.float64
                            }
dtype_default = np.float64

def interp1d(x, y, new_x, 
                    kind = 'linear', low = NaN, high = NaN,
                    bad_data = None):
    # FIXME : all y to be multi-dimensional
    # NOTE : This docstring is considered suboordinate to that for Interpolate1d.
    #       That is, update Interpolate1d and copy-and-paste
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
                y is not currently supported.
                
        Optional Arguments
        -------------------
        
            kind -- Usually a string.  But can be any type.
                Specifies the type of interpolation to use for values within
                the range of x.
                
                By default, linear interpolation is used.
                
                See below for details on other options.
                
            low  -- same as for kind
                How to extrapolate values for inputs below the range of x.
                Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
                a number') for all values below the range of x.
                
            high  -- same as for kind
                How to extrapolate values for inputs above the range of x.
                Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
                a number') for all values above the range of x.
                
            bad_data -- list of numbers
                List of numerical values (in x or y) which indicate unacceptable data. 
                
                If bad_data is not None (its default), all points whose x or y coordinate is in
                bad_data, OR ones of whose coordinates is NaN, will be removed.  Note that
                bad_data != None means NaNs will be removed even if they are not in
                bad_data.
                
        Some Acceptable Input Strings
        ------------------------
        
            "linear" -- linear interpolation : default
            "logarithmic" -- logarithmic interpolation : linear in log space?
            "block" --
            "block_average_above' -- block average above
            "Spline" -- spline interpolation of default order
            "quad", "quadratic" -- spline interpolation order 2
            "cubic" -- spline interpolation order 3
            "quartic" -- spline interpolation order 4
            "quintic" -- spline interpolation order 5
            
        Other options for kind, low, and high
        ---------------------------------------------------
        
            If you choose to use a non-string argument, you must
            be careful to use correct formatting.
            
            If a function is passed, it will be called when interpolating.
            It is assumed to have the form 
                newy = interp(x, y, newx), 
            where x, y, newx, and newy are all numpy arrays.
            
            If a callable class is passed, it is assumed to have format
                instance = Class(x, y).
            which can then be called by
                new_y = instance(new_x)
            
            If a callable object with method "init_xy" or "set_xy" is
            passed, that method will be used to set x and y as follows
                instance.set_xy(x, y)
            and the object will be called during interpolation.
                new_y = instance(new_x)
            If the "init_xy" and "set_xy" are not present, it will be called as
                new_y = argument(x, y, new_x)
                
            A primitive type which is not a string signifies a function
            which is identically that value (e.g. val and 
            lambda x, y, newx : val are equivalent).
            
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
                                kind = kind,
                                low = low,
                                high = high,
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
        
            kind -- Usually a string.  But can be any type.
                Specifies the type of interpolation to use for values within
                the range of x.
                
                By default, linear interpolation is used.
                
                See below for details on other options.
                
            low  -- same as for kind
                How to extrapolate values for inputs below the range of x.
                Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
                a number') for all values below the range of x.
                
            high  -- same as for kind
                How to extrapolate values for inputs above the range of x.
                Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
                a number') for all values above the range of x.
                
            bad_data -- list of numbers
                List of numerical values (in x or y) which indicate unacceptable data. 
                
                If bad_data is not None (its default), all points whose x or y coordinate is in
                bad_data, OR ones of whose coordinates is NaN, will be removed.  Note that
                bad_data != None means NaNs will be removed even if they are not in
                bad_data.
                
        Some Acceptable Input Strings
        ------------------------
        
            "linear" -- linear interpolation : default
            "logarithmic" -- logarithmic interpolation : linear in log space?
            "block" --
            "block_average_above' -- block average above
            "Spline" -- spline interpolation of default order
            "quad", "quadratic" -- spline interpolation order 2
            "cubic" -- spline interpolation order 3
            "quartic" -- spline interpolation order 4
            "quintic" -- spline interpolation order 5
            
        Other options for kind, low, and high
        ---------------------------------------------------
        
            If you choose to use a non-string argument, you must
            be careful to use correct formatting.
            
            If a function is passed, it will be called when interpolating.
            It is assumed to have the form 
                newy = interp(x, y, newx), 
            where x, y, newx, and newy are all numpy arrays.
            
            If a callable class is passed, it is assumed to have format
                instance = Class(x, y).
            which can then be called by
                new_y = instance(new_x)
            
            If a callable object with method "init_xy" or "set_xy" is
            passed, that method will be used to set x and y as follows
                instance.set_xy(x, y)
            and the object will be called during interpolation.
                new_y = instance(new_x)
            If the "init_xy" and "set_xy" are not present, it will be called as
                new_y = argument(x, y, new_x)
                
            A primitive type which is not a string signifies a function
            which is identically that value (e.g. val and 
            lambda x, y, newx : val are equivalent).
            
        Example
        ---------
        
            >>> import numpy
            >>> from interpolate import Interpolate1d
            >>> x = range(5)        # note list is permitted
            >>> y = numpy.arange(5.)
            >>> new_x = [.2, 2.3, 5.6, 7.0]
            >>> interp_func = Interpolate1d(x, y)
            >>> interp_fuc(new_x)
            array([.2, 2.3, 5.6, NaN])
            
    """
    
    def __init__(self, x, y, 
                        kind = 'linear', 
                        low = NaN, 
                        high = NaN,
                        bad_data = None):
        
        # put data into nice format and store it
        self._init_xy(x, y, bad_data)
        
        # store interpolation functions for each range
        self.interp = self._init_interp_method(kind)
        self.extrap_low = self._init_interp_method(low)
        self.extrap_high = self._init_interp_method(high)

    def _init_xy(self, x, y, bad_data):
        # FIXME : no-copying option, in case user has huge dataset.  non-copy + remove bad data should
        #               flash a warning (esp if, in the future, bad values are interpolated).
        
        # remove bad data if applicable
        if bad_data is not None:
            try: # check that bad_data contains only numerical values
                sum_of_bad_data = sum(bad_data)
            except:
                raise TypeError, "bad_data must be either None \
                        or a list of numbers"            
            x, y = self._remove_bad_data(x, y, bad_data)
        
        # check acceptable size and dimensions
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        assert len(x) > 0 and len(y) > 0 , "Arrays cannot be of zero length"
        assert x.ndim == 1 , "x must be one-dimensional"
        assert y.ndim == 1 , "y must be one-dimensional" 
        assert len(x) == len(y) , "x and y must be of the same length"
        
        # select proper dataypes and make arrays
        self._xdtype = dtype_register.setdefault(type(x[0]), dtype_default)
        self._ydtype = dtype_register.setdefault(type(y[0]), dtype_default)
        self._x = atleast_1d_and_contiguous(x, self._xdtype).copy()
        self._y = atleast_1d_and_contiguous(y, self._ydtype).copy()

    def _remove_bad_data(self, x, y, bad_data = []):
        """ removes data points whose x or y coordinate is
            either in bad_data or is a NaN.
        """
        
        bad_data_mask = np.isnan(x) | np.isnan(y)
        for bad_num in bad_data:
              bad_data_mask =  bad_data_mask | (x==bad_num) | (y==bad_num)
              
        x = x[~bad_data_mask]
        y = y[~bad_data_mask]
        return x, y
        
    def _init_interp_method(self, interp_arg):
        """ returns the interpolating function specified by interp_arg.
        """        
        from inspect import isclass, isfunction
        
        # primary usage : user passes a string indicating a known function
        # pick interpolator accordingly
        if isinstance(interp_arg, basestring):
            interpolator = method_register.setdefault(interp_arg, None )
            if interpolator is None: 
                raise TypeError, "input string %s not valid" % interp_arg
        else:
            interpolator = interp_arg
        
        # interpolator is a callable : function, class, or instance of class
        if hasattr(interpolator, '__call__'):
            # function
            if isfunction(interpolator):
                result = lambda newx : interpolator(self._x, self._y, newx)
                
            # callable class 
            elif isclass(interpolator):
                if hasattr(interpolator, 'set_xy'):
                    result = interpolator()
                    result.set_xy(self._x, self._y)
                if hasattr(interpolator, 'init_xy'):
                    result = interpolator()
                    result.init_xy(self._x, self._y)
                else:
                    result = interpolator(self._x, self._y)
                
            # instance of callable class
            else:
                if hasattr(interpolator, 'init_xy'):
                    result = interpolator
                    result.init_xy(self._x, self._y)
                elif hasattr(interpolator, 'set_xy'):
                    result = interpolator
                    result.set_xy(self._x, self._y)
                else:
                    result = lambda new_x : interpolator(self._x, self._y, new_x)
            
        # non-callable : user has passed a default value to always be returned
        else:
            result = np.vectorize(lambda new_x : interp_arg)
        
        return result

    def __call__(self, newx):
        """ Input x must be a list or NumPy array in sorted order.
            
            Breaks x into pieces in-range, below-range, and above range.
            Performs appropriate operation on each and concatenates results.
        """
        
        # record if input is scalar or 0-dimemsional array, in which case output will be scalar
        input_is_scalar = np.isscalar(newx) or \
                                    (
                                        isinstance(  newx , np.ndarray  ) and 
                                        np.shape(newx) == ()
                                    )
        
        # make input into a nice 1d, contiguous array
        newx_array = atleast_1d_and_contiguous(newx, dtype=self._xdtype)
        assert newx_array.ndim == 1, "new_x can be at most 1-dimensional"
        
        # masks indicate which elements fall into which interpolation region
        low_mask = newx_array<self._x[0]
        high_mask = newx_array>self._x[-1]
        interp_mask = (~low_mask) & (~high_mask)
                
        # use correct function for x values in each region and create output as an array
        if len(newx_array[low_mask]) == 0: new_low=np.array([])  # FIXME : remove need for if/else hack.
                                                                            # it's there since vectorize is failing on arrays of zero length
        else: new_low = self.extrap_low(newx_array[low_mask])
        if len(newx_array[interp_mask])==0: new_interp=np.array([])
        else: new_interp = self.interp(newx_array[interp_mask])
        if len(newx_array[high_mask]) == 0: new_high = np.array([])
        else: new_high = self.extrap_high(newx_array[high_mask])
        result_array = np.concatenate((new_low, new_interp, new_high))
        
        # convert to scalar if scalar was passed in
        if input_is_scalar:
            result = float(result_array)
        else:
            result = result_array
        
        return result
  