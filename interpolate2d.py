
from numpy import NaN, array
import numpy as np
from fitpack_wrapper import Spline2d
from algorithm526_wrapper import algorithm526

def atleast_1d_and_contiguous(ary, dtype = np.float64):
    # FIXME : don't have in 2 places
    return np.atleast_1d( np.ascontiguousarray(ary, dtype) )

# dictionary of interpolation functions/classes/objects
# keys are possible values of keyword "kind"
method_register = \
                { 
                    'linear' : Spline2d(kx=1, ky=1),
                    'spline' : Spline2d(),
                    'quadratic' : Spline2d(kx=2, ky=2),
                    'quad' : Spline2d(kx=2, ky=2),
                    'cubic' : Spline2d(kx=3, ky=3),
                    'natural' : Spline2d(kx=3, ky=3),
                    'quartic' : Spline2d(kx=4, ky=4),
                    'quar' : Spline2d(kx=4, ky=4),
                    'quintic' : Spline2d(kx=5, ky=5),
                    'quin' : Spline2d(kx=5, ky=5),
                    '526' : algorithm526, 'algorithm526':algorithm526,
                }
                
# dictionary of types for casting.  key = possible datatype, value = datatype it is cast to
# BEWARE : if you cast things to integers, you will lose interpolation ability
dtype_register = {   
                        np.float32 : np.float32, 
                        np.float64 : np.float64
                        }
# input will be cast to this type if it's not a key in dtype_register
dtype_default = np.float64

# functional interface: creates and calls an instance of objective interface
def interp2d(x, y, z, newx, newy, kind='linear', out=NaN, bad_data=None):
    return Interpolate2d(x, y, z, kind=kind, out=out, bad_data=bad_data)(newx, newy)

# objective interface
class Interpolate2d:
    """ A callable class for interpolation of 1D, real-valued data.
        
        Parameters
        -----------
            
            x -- list or 1D NumPy array
                x includes the x-values for the data set to
                interpolate from.
                    
            y -- list or 1D NumPy array
                y includes the y-values for the data set  to
                interpolate from.
                
            z -- list or 1D NumPy array
                z includes the z-values for the data set to
                interpolate from.
                
        Optional Arguments
        -------------------
        
            kind -- Usually a string.  But can be any type.
                Specifies the type of interpolation to use for values within
                the range of x.
                
                By default, linear interpolation is used.
                
                See below for details on other options.
                
            out  -- same as for kind
                How to extrapolate values for outside the rectangle defined by
                    min(x) <= newx[i] <= max(x)  ,  min(y) <= newy[i] <= max(y)
                Same options as for 'kind'.  Defaults to returning numpy.NaN ('not 
                a number') for all values below the region.
                
            bad_data -- list of numbers
                List of numerical values (in x, y or z) which indicate unacceptable data. 
                
                If bad_data is not None (its default), all points whose x, y or z coordinate is in
                bad_data, OR ones of whose coordinates is NaN, will be removed.  Note that
                bad_data != None means NaNs will be removed even if they are not in
                bad_data.
                
        Some Acceptable Input Strings
        ------------------------
        
            "linear" -- linear interpolation : default
            "spline" -- spline interpolation of default order
            "quad", "quadratic" -- spline interpolation order 2
            "cubic" -- spline interpolation order 3
            "quartic" -- spline interpolation order 4
            "quintic" -- spline interpolation order 5
            
        Other options for kind and out
        ---------------------------------------------------
        
            If you choose to use a non-string argument, you must
            be careful to use correct formatting.
            
            If a function is passed, it will be called when interpolating.
            It is assumed to have the form 
                newz = interp(x, y, z, newx, newy), 
            where x, y, newx, and newy are all numpy arrays.
            
            If a callable class is passed, it is assumed to have format
                instance = Class(x, y, z).
            which can then be called by
                newz = instance(newx, newy)
            
            If a callable object with method "init_xyz" or "set_xyz" is
            passed, that method will be used to set x and y as follows
                instance.set_xy(x, y)
            and the object will be called during interpolation.
                newz = instance(newx, newy)
            If the "init_xyz" and "set_xyz" are not present, it will be called as
                newz = argument(x, y, z, newx, newy)
                
            A primitive type which is not a string signifies a function
            which is identically that value (e.g. val and 
            lambda x, y, newx : val are equivalent).
            
        Example
        ---------
        
            >>> import numpy
            >>> from interpolate import Interpolate2d
            >>> x = range(5)        # note list is permitted
            >>> y = numpy.arange(5.)
            >>> z = x+y
            >>> newx = [.2, 2.3, 2.6, 7.0]
            >>> newy = [1, 1, 1, 1]
            >>> interp_func = Interpolate2d(x, y, z)
            >>> interp_fuc(newx, newy)
            array([1.2, 3.3, 3.6, NaN])
            
    """
    def __init__(self, x, y, z, kind='linear', out=NaN, bad_data=None):
        
        self._init_xyz(x, y, z, bad_data)
        
        self.kind = self._init_interp_method(kind)
        
        self.out = self._init_interp_method(out)
        
    def _init_xyz(self, x, y, z, bad_data):

        # FIXME : perhaps allow 2D input if it is inthe form of meshgrid
         
        # check acceptable sizes and dimensions
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        z = np.atleast_1d(z)
        assert len(x) > 0 and len(y) > 0 and len(z)>0, "Arrays cannot be of zero length"
        assert x.ndim == 1 , "x must be one-dimensional"
        assert y.ndim == 1 , "y must be one-dimensional"
        assert z.ndim == 1 , "z must be one-dimensional" 
        assert len(x) == len(y) , "x and y must be of the same length"
        assert len(x) == len(z) , "x and z must be of the same length"
        
        # remove bad data if applicable
        if bad_data is not None:
            try: # check that bad_data contains only numerical values
                sum_of_bad_data = sum(bad_data)
            except:
                raise TypeError, "bad_data must be either None \
                        or a list of numbers.  Sorry."  
            x, y, z = self._remove_bad_data(x, y, z, bad_data)
            
        # select proper dataypes and make arrays
        self._xdtype = dtype_register.setdefault(type(x[0]), dtype_default)
        self._ydtype = dtype_register.setdefault(type(y[0]), dtype_default)
        self._zdtype = dtype_register.setdefault(type(z[0]), dtype_default)
        self._x = atleast_1d_and_contiguous(x, self._xdtype).copy()
        self._y = atleast_1d_and_contiguous(y, self._ydtype).copy()
        self._z = atleast_1d_and_contiguous(z, self._zdtype).copy()
                
    def _init_interp_method(self, method):
        """ returns the interpolating function specified by interp_arg.
        """        
        from inspect import isclass, isfunction
        
        # primary usage : user passes a string indicating a known function
        # pick interpolator accordingly
        if isinstance(method, basestring):
            interpolator = method_register.setdefault(method.lower(), None )
            if interpolator is None: 
                raise TypeError, "input string %s not valid" % method
        else:
            interpolator = method
        
        # interpolator is a callable : function, class, or instance of class
        if hasattr(interpolator, '__call__'):
            # function
            if isfunction(interpolator):
                result = lambda newx, newy : interpolator(self._x, self._y, self._z, newx, newy)
                
            # callable class 
            elif isclass(interpolator):
                if hasattr(interpolator, 'set_xyz'):
                    result = interpolator()
                    result.set_xyz(self._x, self._y, self._z)
                if hasattr(interpolator, 'init_xyz'):
                    result = interpolator()
                    result.init_xyz(self._x, self._y, self._z)
                else:
                    result = interpolator(self._x, self._y, self._z)
                
            # instance of callable class
            else:
                if hasattr(interpolator, 'init_xyz'):
                    result = interpolator
                    result.init_xyz(self._x, self._y, self._z)
                elif hasattr(interpolator, 'set_xyz'):
                    result = interpolator
                    result.set_xyz(self._x, self._y, self._z)
                else:
                    result = lambda newx, newy : interpolator(self._x, self._y, self._z, newx, newy)
            
        # non-callable : user has passed a default value to always be returned
        else:
            result = np.vectorize(lambda newx, newy : interpolator)
        
        return result
        
    def _remove_bad_data(self, x, y, z, bad_data):
        """ removes data points whose x or y coordinate is
            either in bad_data or is a NaN.
        """

        bad_data_mask = np.isnan(x) | np.isnan(y) | np.isnan(z)
        for bad_num in bad_data:
            bad_data_mask =  \
                    bad_data_mask | (x==bad_num) | (y==bad_num) | (z==bad_num)
        
        x = x[~bad_data_mask]
        y = y[~bad_data_mask]
        z = z[~bad_data_mask]
        
        return x, y, z
        
    def __call__(self, newx, newy):
        
        # record if input is scalar or 0-dimemsional array, in which case output will be scalar
        input_is_scalar = np.isscalar(newx) or  np.isscalar(newy) or \
                                isinstance(  newx , np.ndarray  ) and np.shape(newx) == () or \
                                isinstance(  newy , np.ndarray  ) and np.shape(newy) == ()
        
        # make input into a nice 1d, contiguous array
        newx = atleast_1d_and_contiguous(newx, dtype=self._xdtype)
        newy = atleast_1d_and_contiguous(newy, dtype=self._ydtype)
        assert newx.ndim == 1, "newx can be at most 1-dimensional"
        assert newy.ndim == 1, "newy can be at most 1-dimensional"
        assert len(newx) == len(newy), "newx and newy must be the same length"
        
        in_range_mask = (min(self._x) <= newx)  & (newx <= max(self._x)) & \
                                (min(self._y) <= newy) & (newy <= max(self._y))        
        
        # filling array of interpolated z-values
        result = np.zeros(np.shape(newx), dtype = self._zdtype)
        if sum(in_range_mask) > 0:  # if there are in-range values.  hack to deal
                                               # with behavior of np.vectorize on arrays of length 0
            result[in_range_mask] = self.kind(newx[in_range_mask], newy[in_range_mask])        
        if sum(~in_range_mask) > 0:
            result[~in_range_mask] = self.out(newx[~in_range_mask], newy[~in_range_mask])
        
        # revert to scalar if applicable
        if input_is_scalar:
            result = result[0]
        
        return result