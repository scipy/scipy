""" Class for interpolating values

    !! Need to find argument for keeping initialize.  If it isn't
    !! found, get rid of it!
"""

__all__ = ['interp1d']

from scipy_base import *
from scipy_base.fastumath import *

import fitpack

# The following are cluges to fix brain-deadness of take and
# sometrue when dealing with 0 dimensional arrays.
# Shouldn't they go to scipy_base??

_take = take
def take(a,indices,axis=0):    
    x = asarray(a); y = asarray(indices)
    if shape(x) == (): x = x.flat
    if shape(y) == (): y = y.flat
    return _take(x,y,axis)

_sometrue = sometrue
def sometrue(a,axis=0):    
    x = asarray(a)
    if shape(x) == (): x = x.flat
    return _sometrue(x)

def reduce_sometrue(a):
    all = a
    while len(shape(all)) > 1:    
        all = sometrue(all)
    return all

class interp2d:
    def __init__(self,x,y,z,kind='linear',
                 copy=1,bounds_error=0,fill_value=None):
        """
        Input:
          x,y  - 1-d arrays defining 2-d grid (or 2-d meshgrid arrays)
          z    - 2-d array of grid values
          kind - interpolation type ('nearest', 'linear', 'cubic', 'spline')
          copy - if true then data is copied into class, otherwise only a
                   reference is held.
          bounds_error - if true, then when out_of_bounds occurs, an error is
                          raised otherwise, the output is filled with
                          fill_value.
          fill_value - if None, then NaN, otherwise the value to fill in
                        outside defined region.
        """
        self.x = atleast_1d(x).copy()
        self.y = atleast_1d(y).copy()
        if rank(self.x) > 2 or rank(self.y) > 2:
            raise ValueError, "One of the input arrays is not 1-d or 2-d."
        if rank(self.x) == 2:
            self.x = self.x[:,0]
        if rank(self.y) == 2:
            self.y = self.y[0]
        self.z = array(z,copy=1)
        if rank(z) != 2:
            raise ValueError, "Grid values is not a 2-d array."


        

    def __call__(self,x,y,dx=0,dy=0):
        """
        Input:
          x,y   - 1-d arrays defining points to interpolate.
          dx,dy - order of partial derivatives in x and y, respectively.
                  0<=dx<kx, 0<=dy<ky
        Output:
          z     - 2-d array of interpolated values
        """
        x = atleast_1d(x)
        y = atleast_1d(y)
        z,ier=fitpack._fitpack._bispev(*(self.tck+[x,y,dx,dy]))
        if ier==10: raise ValueError,"Invalid input data"
        if ier: raise TypeError,"An error occurred"
        z.shape=len(x),len(y)
        z = transpose(z)
        if len(z)==1: z = z[0]
        return array(z)

class interp1d:
    interp_axis = -1 # used to set which is default interpolation
                     # axis.  DO NOT CHANGE OR CODE WILL BREAK.
                     
    def __init__(self,x,y,kind='linear',axis = -1,
                 copy = 1,bounds_error=1, fill_value=None):
        """Initialize a 1d linear interpolation class

        Description:
          x and y are arrays of values used to approximate some function f:
            y = f(x)
          This class returns a function whose call method uses linear
          interpolation to find the value of new points.

        Inputs:
            x -- a 1d array of monotonically increasing real values.
                 x cannot include duplicate values. (otherwise f is
                 overspecified)
            y -- an nd array of real values.  y's length along the
                 interpolation axis must be equal to the length
                 of x.
            kind -- specify the kind of interpolation: 'nearest', 'linear',
                    'cubic', or 'spline'
            axis -- specifies the axis of y along which to 
                    interpolate. Interpolation defaults to the last
                    axis of y.  (default: -1)
            copy -- If 1, the class makes internal copies of x and y.
                    If 0, references to x and y are used. The default 
                    is to copy. (default: 1)
            bounds_error -- If 1, an error is thrown any time interpolation
                            is attempted on a value outside of the range
                            of x (where extrapolation is necessary).
                            If 0, out of bounds values are assigned the
                            NaN (#INF) value.  By default, an error is
                            raised, although this is prone to change.
                            (default: 1)
        """      
        self.axis = axis
        self.copy = copy
        self.bounds_error = bounds_error
        if fill_value is None:
            self.fill_value = array(0.0) / array(0.0)
        else:
            self.fill_value = fill_value

        if kind != 'linear':
            raise NotImplementedError, "Only linear supported for now. Use fitpack routines for other types."
        
        # Check that both x and y are at least 1 dimensional.    
        if len(shape(x)) == 0 or len(shape(y)) == 0:
            raise ValueError, "x and y arrays must have at least one dimension."  
        # make a "view" of the y array that is rotated to the
        # interpolation axis.  
        oriented_x = x
        oriented_y = swapaxes(y,self.interp_axis,axis)            
        interp_axis = self.interp_axis        
        len_x,len_y = shape(oriented_x)[interp_axis], shape(oriented_y)[interp_axis]
        if len_x != len_y:
            raise ValueError, "x and y arrays must be equal in length along "\
                              "interpolation axis."
        if len_x < 2 or len_y < 2:
            raise ValueError, "x and y arrays must have more than 1 entry"            
        self.x = array(oriented_x,copy=self.copy)
        self.y = array(oriented_y,copy=self.copy)       
        
    def __call__(self,x_new):
        """Find linearly interpolated y_new = <name>(x_new).

        Inputs:        
          x_new -- New independent variables.

        Outputs:
          y_new -- Linearly interpolated values corresponding to x_new.
        """
        # 1. Handle values in x_new that are outside of x.  Throw error,
        #    or return a list of mask array indicating the outofbounds values.
        #    The behavior is set by the bounds_error variable.
        x_new = atleast_1d(x_new)
        out_of_bounds = self._check_bounds(x_new)
        # 2. Find where in the orignal data, the values to interpolate
        #    would be inserted.  
        #    Note: If x_new[n] = x[m], then m is returned by searchsorted.
        x_new_indices = searchsorted(self.x,x_new)
        # 3. Clip x_new_indices so that they are within the range of 
        #    self.x indices and at least 1.  Removes mis-interpolation
        #    of x_new[n] = x[0] 
        x_new_indices = clip(x_new_indices,1,len(self.x)-1).astype(Int)
        # 4. Calculate the slope of regions that each x_new value falls in.
        lo = x_new_indices - 1; hi = x_new_indices        
        
        # !! take() should default to the last axis (IMHO) and remove
        # !! the extra argument.
        x_lo = take(self.x,lo,axis=self.interp_axis)
        x_hi = take(self.x,hi,axis=self.interp_axis);
        y_lo = take(self.y,lo,axis=self.interp_axis)
        y_hi = take(self.y,hi,axis=self.interp_axis);
        slope = (y_hi-y_lo)/(x_hi-x_lo)
        # 5. Calculate the actual value for each entry in x_new.
        y_new = slope*(x_new-x_lo) + y_lo 
        # 6. Fill any values that were out of bounds with NaN
        # !! Need to think about how to do this efficiently for 
        # !! mutli-dimensional Cases.
        yshape = y_new.shape
        y_new = y_new.flat
        new_shape = list(yshape)
        new_shape[self.interp_axis] = 1
        sec_shape = [1]*len(new_shape)
        sec_shape[self.interp_axis] = len(out_of_bounds)
        out_of_bounds.shape = sec_shape
        new_out = ones(new_shape)*out_of_bounds
        putmask(y_new, new_out.flat, self.fill_value)
        y_new.shape = yshape      
        # Rotate the values of y_new back so that they coorespond to the
        # correct x_new values.
        result = swapaxes(y_new,self.interp_axis,self.axis)
        try:
            len(x_new)
            return result
        except TypeError:
            return result[0]
        return result
    
    def _check_bounds(self,x_new):
        # If self.bounds_error = 1, we raise an error if any x_new values
        # fall outside the range of x.  Otherwise, we return an array indicating
        # which values are outside the boundary region.  
        # !! Needs some work for multi-dimensional x !!
        below_bounds = less(x_new,self.x[0])
        above_bounds = greater(x_new,self.x[-1])
        #  Note: sometrue has been redefined to handle length 0 arrays
        # !! Could provide more information about which values are out of bounds
        if self.bounds_error and sometrue(below_bounds):
            raise ValueError, " A value in x_new is below the"\
                              " interpolation range."
        if self.bounds_error and sometrue(above_bounds):
            raise ValueError, " A value in x_new is above the"\
                              " interpolation range."
        # !! Should we emit a warning if some values are out of bounds.
        # !! matlab does not.
        out_of_bounds = logical_or(below_bounds,above_bounds)
        return out_of_bounds
       
    def model_error(self,x_new,y_new):
        # How well do x_new,yy points fit the model?
        # Return an array of error values.
        pass

    
#assumes module test_xxx is in python path
#def test():
#    test_module = 'test_' + __name__ # __name__ is name of this module
#    test_string = 'import %s;reload(%s);%s.test()' % ((test_module,)*3)
#    exec(test_string)

#if __name__ == '__main__':
#    test()
