""" Classes for interpolating values.
"""

# !! Need to find argument for keeping initialize.  If it isn't
# !! found, get rid of it!

__all__ = ['interp1d', 'interp2d']

from numpy import shape, sometrue, rank, array, transpose, \
     swapaxes, searchsorted, clip, take, ones, putmask, less, greater, \
     logical_or, atleast_1d, atleast_2d, meshgrid
import numpy as np

import fitpack

def reduce_sometrue(a):
    all = a
    while len(shape(all)) > 1:
        all = sometrue(all,axis=0)
    return all

class interp2d(object):
    """ Interpolate over a 2D grid.

    See Also
    --------
    bisplrep, bisplev - spline interpolation based on FITPACK
    BivariateSpline - a more recent wrapper of the FITPACK routines
    """

    def __init__(self, x, y, z, kind='linear', copy=True, bounds_error=False,
        fill_value=np.nan):
        """ Initialize a 2D interpolator.

        Parameters
        ----------
        x : 1D array or 2D meshgrid array
        y : 1D array or 2D meshgrid array
            Arrays defining the coordinates of a 2D grid.
        z : 2D array
            The values of the interpolated function on the grid points.
        kind : 'linear', 'cubic', 'quintic'
            The kind of interpolation to use.
        copy : bool
            If True, then data is copied, otherwise only a reference is held.
        bounds_error : bool
            If True, when interoplated values are requested outside of the
            domain of the input data, an error is raised.
            If False, then fill_value is used.
        fill_value : number
            If provided, the value to use for points outside of the
            interpolation domain. Defaults to NaN.

        Raises
        ------
        ValueError when inputs are invalid.

        """        
        self.x = atleast_1d(x).copy()
        self.y = atleast_1d(y).copy()
        if rank(self.x) == 1 and rank(self.y) == 1:
            # Prepare meshgrid for bisplrep
            x, y = meshgrid(x,y)        
        if rank(self.x) > 2 or rank(self.y) > 2:
            raise ValueError("One of the input arrays is not 1-d or 2-d.")
        if rank(self.x) == 2:
            self.x = self.x[:,0]
        if rank(self.y) == 2:
            self.y = self.y[0]
        self.z = array(z,copy=True)
        if rank(z) != 2:
            raise ValueError("Grid values is not a 2-d array.")
        try:
            kx = ky = {'linear' : 1,
                       'cubic' : 3,
                       'quintic' : 5}[kind]
        except KeyError:
            raise ValueError("Unsupported interpolation type.")
        self.tck = fitpack.bisplrep(x, y, z, kx=kx, ky=ky, s=0.)

    def __call__(self,x,y,dx=0,dy=0):
        """ Interpolate the function.

        Parameters
        ----------
        x : 1D array
        y : 1D array
            The points to interpolate.
        dx : int >= 0, < kx
        dy : int >= 0, < ky
            The order of partial derivatives in x and y, respectively.

        Returns
        -------
        z : 2D array with shape (len(y), len(x))
            The interpolated values.
        """

        x = atleast_1d(x)
        y = atleast_1d(y)
        z = fitpack.bisplev(x, y, self.tck, dx, dy)
        z = atleast_2d(z)
        z = transpose(z)
        if len(z)==1:
            z = z[0]
        return array(z)


class interp1d(object):
    """ Interpolate a 1D function.
    
    See Also
    --------
    splrep, splev - spline interpolation based on FITPACK
    UnivariateSpline - a more recent wrapper of the FITPACK routines
    """

    _interp_axis = -1 # used to set which is default interpolation
                      # axis.  DO NOT CHANGE OR CODE WILL BREAK.

    def __init__(self, x, y, kind='linear', axis=-1,
                 copy=True, bounds_error=True, fill_value=np.nan):
        """ Initialize a 1D linear interpolation class.

        Description
        -----------
        x and y are arrays of values used to approximate some function f:
            y = f(x)
        This class returns a function whose call method uses linear
        interpolation to find the value of new points.

        Parameters
        ----------
        x : array
            A 1D array of monotonically increasing real values.  x cannot
            include duplicate values (otherwise f is overspecified)
        y : array
            An N-D array of real values.  y's length along the interpolation
            axis must be equal to the length of x.
        kind : str
            Specifies the kind of interpolation. At the moment, only 'linear' is
            implemented.
        axis : int
            Specifies the axis of y along which to interpolate. Interpolation
            defaults to the last axis of y.
        copy : bool
            If True, the class makes internal copies of x and y.  
            If False, references to x and y are used.
            The default is to copy.
        bounds_error : bool
            If True, an error is thrown any time interpolation is attempted on
            a value outside of the range of x (where extrapolation is
            necessary).
            If False, out of bounds values are assigned fill_value.
            By default, an error is raised.
        fill_value : float
            If provided, then this value will be used to fill in for requested
            points outside of the data range.
            If not provided, then the default is NaN.
        """

        self.copy = copy
        self.bounds_error = bounds_error
        self.fill_value = fill_value

        if kind != 'linear':
            raise NotImplementedError("Only linear supported for now. Use "
                "fitpack routines for other types.")

        x = array(x, copy=self.copy)
        y = array(y, copy=self.copy)

        if len(x.shape) != 1:
            raise ValueError("the x array must have exactly one dimension.")
        if len(y.shape) == 0:
            raise ValueError("the y array must have at least one dimension.")

        # Normalize the axis to ensure that it is positive.
        self.axis = axis % len(y.shape)

        # Make a "view" of the y array that is rotated to the interpolation
        # axis.
        oriented_y = y.swapaxes(self._interp_axis, axis)
        len_x = len(x)
        len_y = oriented_y.shape[self._interp_axis]
        if len_x != len_y:
            raise ValueError("x and y arrays must be equal in length along"
                "interpolation axis.")
        if len_x < 2 or len_y < 2:
            raise ValueError("x and y arrays must have more than 1 entry")
        self.x = x
        self.y = oriented_y

    def __call__(self, x_new):
        """ Find linearly interpolated y_new = f(x_new).

        Parameters
        ----------
        x_new : number or array
            New independent variable(s).

        Returns
        -------
        y_new : number or array
            Linearly interpolated value(s) corresponding to x_new.
        """

        # 1. Handle values in x_new that are outside of x.  Throw error,
        #    or return a list of mask array indicating the outofbounds values.
        #    The behavior is set by the bounds_error variable.
        x_new = atleast_1d(x_new)
        out_of_bounds = self._check_bounds(x_new)

        # 2. Find where in the orignal data, the values to interpolate
        #    would be inserted.
        #    Note: If x_new[n] == x[m], then m is returned by searchsorted.
        x_new_indices = searchsorted(self.x, x_new)

        # 3. Clip x_new_indices so that they are within the range of
        #    self.x indices and at least 1.  Removes mis-interpolation
        #    of x_new[n] = x[0]
        x_new_indices = x_new_indices.clip(1, len(self.x)-1).astype(int)

        # 4. Calculate the slope of regions that each x_new value falls in.
        lo = x_new_indices - 1
        hi = x_new_indices

        x_lo = self.x[lo]
        x_hi = self.x[hi]
        y_lo = self.y[..., lo]
        y_hi = self.y[..., hi]

        # Note that the following two expressions rely on the specifics of the
        # broadcasting semantics.
        slope = (y_hi-y_lo) / (x_hi-x_lo)

        # 5. Calculate the actual value for each entry in x_new.
        y_new = slope*(x_new-x_lo) + y_lo

        # 6. Fill any values that were out of bounds with fill_value.
        y_new[..., out_of_bounds] = self.fill_value

        # Rotate the values of y_new back so that they correspond to the
        # correct x_new values. For N-D x_new, take the last N axes from y_new
        # and insert them where self.axis was in the list of axes.
        nx = len(x_new.shape)
        ny = len(y_new.shape)
        axes = range(ny - nx)
        axes[self.axis:self.axis] = range(ny - nx, ny)
        result = y_new.transpose(axes)

        return result

    def _check_bounds(self, x_new):
        """ Check the inputs for being in the bounds of the interpolated data.

        Parameters
        ----------
        x_new : array

        Returns
        -------
        out_of_bounds : bool array
            The mask on x_new of values that are out of the bounds.
        """

        # If self.bounds_error is True, we raise an error if any x_new values
        # fall outside the range of x.  Otherwise, we return an array indicating
        # which values are outside the boundary region.
        below_bounds = x_new < self.x[0]
        above_bounds = x_new > self.x[-1]

        # !! Could provide more information about which values are out of bounds
        if self.bounds_error and below_bounds.any():
            raise ValueError("A value in x_new is below the interpolation "
                "range.")
        if self.bounds_error and above_bounds.any():
            raise ValueError("A value in x_new is above the interpolation "
                "range.")

        # !! Should we emit a warning if some values are out of bounds?
        # !! matlab does not.
        out_of_bounds = logical_or(below_bounds, above_bounds)
        return out_of_bounds

