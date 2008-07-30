# FIXME: information strings giving mathematical descriptions of the actions
#     of the functions.

from interpolate_wrapper import linear, logarithmic, block, block_average_above, atleast_1d_and_contiguous
from fitpack_wrapper import Spline
import numpy as np
from numpy import array, arange, empty, float64, NaN

# dictionary of tuples.  First element is a callable (class, instance of a class, or function
# second argument is dictionary of additional keywords, if any
dict_of_interp_types = \
                { 'linear' : (linear, {}), 
                    'logarithmic' : (logarithmic, {}), 
                    'block' : (block, {}),
                    'block_average_above' : (block_average_above, {}),
                    'Spline' : (Spline, {}), 'spline' : (Spline, {}),
                    'Quadratic' : (Spline, {'k':2}), 'quadratic' : (Spline, {'k':2}),
                    'Quad' : (Spline, {'k':2}), 'quad' : (Spline, {'k':2}),
                    'Cubic' : (Spline, {'k':3}), 'cubic' : (Spline, {'k':3}),
                    'Quartic' : (Spline, {'k':4}), 'quartic' : (Spline, {'k':4}),
                    'Quar' : (Spline, {'k':4}), 'quar' : (Spline, {'k':4}),
                    'Quintic' : (Spline, {'k':5}), 'quintic' : (Spline, {'k':5}),
                    'Quin' : (Spline, {'k':5}), 'quin' : (Spline, {'k':5})
                }

def interp1d(x, y, new_x, 
                    interp = 'linear', extrap_low = NaN, extrap_high = NaN,
                    interpkw = {}, lowkw = {}, highkw ={},
                    bad_data = None):
    # FIXME : all y to be multi-dimensional
    # FIXME : update the doc string to match that of Interpolate1d
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
        
        interp -- Usu. function or string.  But can be any type.
            Specifies the type of extrapolation to use for values within
            the range of x.  If a string is passed, it will look for an object
            or function with that name and call it when evaluating.  If 
            a function or object is passed, it will be called when interpolating.
            If nothing else, assumes the argument is intended as a value
            to be returned for all arguments.  Defaults to linear interpolation.
            
        low (high) -- same as for interp
            Same options as for 'interp'.  Defaults to returning numpy.NaN ('not 
            a number') for all values outside the range of x.
        
        interpkw -- dictionary
            If 
            
        bad_data -- list
            List of values (in x or y) which indicate unacceptable data. All points
            that have x or y value in missing_data will be removed before
            any interpolatin is performed if bad_data is not None.
            
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
                                extrap_low = extrap_low,
                                extrap_high = extrap_high,
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
        
            interp -- Usually a string.  But can be any type.
                Specifies the type of interpolation to use for values within
                the range of x.
                
                By default, linear interpolation is used.
                
                See below for details on other options.
                
            extrap_low  -- same as for kind
                How to extrapolate values for inputs below the range of x.
                Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
                a number') for all values below the range of x.
                
            extrap_high  -- same as for kind
                How to extrapolate values for inputs above the range of x.
                Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
                a number') for all values above the range of x.
                
            bad_data -- list
                List of numerical values (in x or y) which indicate unacceptable data. 
                
                If bad_data is not None (its default), all points whose x or y coordinate is in
                bad_data, OR ones of whose coordinates is NaN, will be removed.
                
            interpkw -- dictionary
                If interp is set to a function, class or callable object, this contains
                additional keywords.
                
            lowkw (highkw) -- dictionary
                like interpkw, but for extrap_low and extrap_high
                
            
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
            
        Other options for interp, extrap_low, and extrap_high
        ---------------------------------------------------
        
            If you choose to use a non-string argument, you must
            be careful to use correct formatting.
            
            If a function is passed, it will be called when interpolating.
            It is assumed to have the form 
                newy = interp(x, y, newx, **kw), 
            where x, y, newx, and newy are all numpy arrays.
            
            If a callable class is passed, it is assumed to have format
                instance = Class(x, y, **kw).
            which can then be called by
                new_y = instance(new_x)
            
            If a callable object with method "init_xy" or "set_xy" is
            passed, that method will be used to set x and y as follows
                instance.set_xy(x, y, **kw)
            and the object will be called during interpolation.
                new_y = instance(new_x)
            If the "init_xy" and "set_xy" are not present, it will be called as
                new_y = argument(new_x)
                
            A primitive type which is not a string signifies a function
            which is identically that value (e.g. val and 
            lambda x, y, newx : val are equivalent).
            
        Example
        ---------
        
            >>> import numpy
            >>> from interpolate1d import Interpolate1d
            >>> x = range(5)        # note list is permitted
            >>> y = numpy.arange(5.)
            >>> new_x = [.2, 2.3, 5.6, 7.0]
            >>> interp_func = Interpolate1d(x, y)
            >>> interp_fuc(new_x)
            array([.2, 2.3, 5.6, NaN])
            
    """
    # FIXME: more informative descriptions of sample arguments
    # FIXME: examples in doc string
    # FIXME : Allow copying or not of arrays.  non-copy + remove_bad_data should flash 
    #           a warning (esp if we interpolate missing values), but work anyway.
    
    def __init__(self, x, y, 
                        interp = 'linear', 
                        extrap_low = NaN, 
                        extrap_high = NaN,
                        interpkw = {},
                        lowkw = {},
                        highkw = {},
                        bad_data = None):
        # FIXME: don't allow copying multiple times.
        # FIXME : allow no copying, in case user has huge dataset
        
        # remove bad data, is there is any
        if bad_data is not None:
            try:
                sum_of_bad_data = sum(bad_data)
            except:
                raise TypeError, "bad_data must be either None \
                        or a list of numerical types"
            
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
        self.extrap_low = self._init_interp_method(extrap_low, lowkw)
        self.extrap_high = self._init_interp_method(extrap_high, highkw)

    def _init_xy(self, x, y):
        
        # select proper dataypes and make arrays
        self._xdtype = {np.float32 : np.float32}.setdefault(type(x[0]), np.float64) # unless data is float32,  cast to float64
        self._ydtype = {np.float32 : np.float32}.setdefault(type(y[0]), np.float64)
        self._x = atleast_1d_and_contiguous(x, self._xdtype).copy()
        self._y = atleast_1d_and_contiguous(y, self._ydtype).copy()

    def _remove_bad_data(self, x, y, bad_data = []):
        """ removes data points whose x or y coordinate is
            either in bad_data or is a NaN.
        """
        # FIXME : In the future, it may be good to just replace the bad points with good guesses.
        #       Especially in generalizing the higher dimensions
        # FIXME : This step is very inefficient because it iterates over the array
        
        bad_data_mask = np.isnan(x) | np.isnan(y)
        for bad_num in bad_data:
              bad_data_mask =  bad_data_mask | (x==bad_num) | (y==bad_num)
              
        x = x[~bad_data_mask]
        y = y[~bad_data_mask]
        return x, y
        
    def _init_interp_method(self, interp_arg, kw):
        """
            returns the interpolating function specified by interp_arg.
        """
        # FIXME : error checking specific to interpolation method.  x and y long
        #   enough for order-3 spline if that's indicated, etc.  Functions should throw
        #   errors themselves, but errors at instantiation would be nice.
        
        from inspect import isclass, isfunction
        
        # primary usage : user passes a string indicating a known function
        if isinstance(interp_arg, basestring):
            interpolator, kw = dict_of_interp_types.setdefault(interp_arg, (None, {}) )
            
            if interpolator is None: 
                raise TypeError, "input string %s not valid" % interp_arg
        else:
            interpolator = interp_arg
        
        # interpolator is a callable : function, class, or instance of class
        if hasattr(interpolator, '__call__'):
            # function
            if isfunction(interpolator):
                result = lambda newx : interpolator(self._x, self._y, newx, **kw)
                
            # callable class 
            elif isclass(interpolator):
                if hasattr(interpolator, 'set_xy'):
                    result = interpolator(**kw)
                    result.set_xy(self._x, self._y)
                if hasattr(interpolator, 'init_xy'):
                    result = interpolator(**kw)
                    result.init_xy(self._x, self._y)
                else:
                    result = interpolator(self._x, self._y, **kw)
                
            # instance of callable class
            else:
                if hasattr(interpolator, 'init_xy'):
                    result = interpolator
                    result.init_xy(self._x, self._y, **kw)
                elif hasattr(interpolator, 'set_xy'):
                    result = interpolator
                    result.set_xy(self._x, self._y, **kw)
                else:
                    result = interpolator
            
        # non-callable : user has passed a default value to always be returned
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
        input_is_scalar = np.isscalar(newx) or \
                                    (
                                        isinstance(  newx , np.ndarray  ) and 
                                        np.shape(newx) == ()
                                    )
        
        # make 
        newx_array = atleast_1d_and_contiguous(newx)
        
        # masks indicate which elements fall into which interpolation region
        low_mask = newx_array<self._x[0]
        high_mask = newx_array>self._x[-1]
        interp_mask = (~low_mask) & (~high_mask)
                
        type(newx_array[low_mask])
                
                
        # use correct function for x values in each region
        if len(newx_array[low_mask]) == 0: new_low=np.array([])  # FIXME : remove need for if/else.
                                                                            # if/else is a hack, since vectorize is failing
                                                                            # to work on lists/arrays of length 0
                                                                            # on the computer where this is being
                                                                            # developed
        else: new_low = self.extrap_low(newx_array[low_mask])
        if len(newx_array[interp_mask])==0: new_interp=np.array([])
        else: new_interp = self.interp(newx_array[interp_mask])
        if len(newx_array[high_mask]) == 0: new_high = np.array([])
        else: new_high = self.extrap_high(newx_array[high_mask])
        
        result_array = np.concatenate((new_low, new_interp, new_high)) # FIXME : deal with mixed datatypes
                                                                                          # Would be nice to say result = zeros(dtype=?)
                                                                                          # and fill in
        
        # convert to scalar if scalar was passed in
        if input_is_scalar:
            result = float(result_array)
        else:
            result = result_array
        
        return result
  