""" Classes for interpolating values.
"""

# !! Need to find argument for keeping initialize.  If it isn't
# !! found, get rid of it!

__all__ = ['interp1d', 'interp2d', 'cspline', 'cspeval', 'csprep', 'csp2pp',
           'ppval']

from numpy import shape, sometrue, rank, array, transpose, \
     swapaxes, searchsorted, clip, take, ones, putmask, less, greater, \
     logical_or, atleast_1d, atleast_2d, meshgrid, ravel
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
        x : 1D array
        y : 1D array
            Arrays defining the coordinates of a 2D grid.  If the
            points lie on a regular grid, x and y can simply specify
            the rows and colums, i.e.

            x = [0,1,2]  y = [0,1,2]
            
            otherwise x and y must specify the full coordinates, i.e.
            
            x = [0,1,2,0,1.5,2,0,1,2]  y = [0,1,2,0,1,2,0,1,2]
            
            If x and y are 2-dimensional, they are flattened (allowing
            the use of meshgrid, for example).
            
        z : 1D array
            The values of the interpolated function on the grid
            points. If z is a 2-dimensional array, it is flattened.

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
        self.x, self.y, self.z = map(ravel, map(array, [x, y, z]))
        if not map(rank, [self.x, self.y, self.z]) == [1,1,1]:
            raise ValueError("One of the input arrays is not 1-d.")
        if len(self.x) != len(self.y):
            raise ValueError("x and y must have equal lengths")
        if len(self.z) == len(self.x) * len(self.y):
            self.x, self.y = meshgrid(x,y)
            self.x, self.y = map(ravel, [self.x, self.y])
        if len(self.z) != len(self.x):
            raise ValueError("Invalid length for input z")
        
        try:
            kx = ky = {'linear' : 1,
                       'cubic' : 3,
                       'quintic' : 5}[kind]
        except KeyError:
            raise ValueError("Unsupported interpolation type.")
        
        self.tck = fitpack.bisplrep(self.x, self.y, self.z, kx=kx, ky=ky, s=0.)

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


def _get_cspline_Bb(xk, yk, kind, conds):
    # internal function to compute different tri-diagonal system
    # depending on the kind of spline requested.
    # conds is only used for 'second' and 'first' 
    Np1 = len(xk)
    if kind in ['natural', 'second']:
        if kind == 'natural':
            m0, mN = 0.0, 0.0
        else:
            m0, mN = conds

        # the matrix to invert is (N-1,N-1)
        beta = 2*(xk[2:]-xk[:-2])
        alpha = xk[1:]-xk[:-1]
        B = np.diag(alpha[1:-1],k=-1) + np.diag(beta) + np.diag(alpha[2:],k=1)
        dyk = yk[1:]-yk[:-1]
        b = (dyk[1:]/alpha[1:] - dyk[:-1]/alpha[:-1])
        b *= 6
        b[0] -= m0
        b[-1] -= mN

        # put m0 and mN into the correct shape for
        #  concatenation
        m0 = array(m0,copy=0,ndmin=yk.ndim)
        mN = array(mN,copy=0,ndmin=yk.ndim)
        if m0.shape[1:] != yk.shape[1:]:
            m0 = m0*(ones(yk.shape[1:])[newaxis,...])
        if mN.shape[1:] != yk.shape[1:]:
            mN = mN*(ones(yk.shape[1:])[newaxis,...])

        return B, b, m0, mN           


    elif kind in ['clamped', 'endslope', 'first', 'not-a-knot', 'runout',
                  'parabolic']:
        if kind == 'endslope':
            # match slope of lagrange interpolating polynomial of
            # order 3 at end-points. 
            x0,x1,x2,x3 = xk[:4]
            sl_0 = (1./(x0-x1)+1./(x0-x2)+1./(x0-x3))*yk[0]
            sl_0 += (x0-x2)*(x0-x3)/((x1-x0)*(x1-x2)*(x1-x3))*yk[1]
            sl_0 += (x0-x1)*(x0-x3)/((x2-x0)*(x2-x1)*(x3-x2))*yk[2]
            sl_0 += (x0-x1)*(x0-x2)/((x3-x0)*(x3-x1)*(x3-x2))*yk[3]

            xN3,xN2,xN1,xN0 = xk[-4:]
            sl_N = (1./(xN0-xN1)+1./(xN0-xN2)+1./(xN0-xN3))*yk[-1]
            sl_N += (xN0-xN2)*(xN0-xN3)/((xN1-xN0)*(xN1-xN2)*(xN1-xN3))*yk[-2]
            sl_N += (xN0-xN1)*(xN0-xN3)/((xN2-xN0)*(xN2-xN1)*(xN3-xN2))*yk[-3]
            sl_N += (xN0-xN1)*(xN0-xN2)/((xN3-xN0)*(xN3-xN1)*(xN3-xN2))*yk[-4]
        elif kind == 'clamped':
            sl_0, sl_N = 0.0, 0.0
        elif kind == 'first':
            sl_0, sl_N = conds

        # Now set up the (N+1)x(N+1) system of equations
        beta = np.r_[0,2*(xk[2:]-xk[:-2]),0]
        alpha = xk[1:]-xk[:-1]
        gamma = np.r_[0,alpha[1:]]
        B = np.diag(alpha,k=-1) + np.diag(beta) + np.diag(gamma,k=1)
        d1 = alpha[0]
        dN = alpha[-1]
        if kind == 'not-a-knot':
            d2 = alpha[1]
            dN1 = alpha[-2]
            B[0,:3] = [d2,-d1-d2,d1]
            B[-1,-3:] = [dN,-dN1-dN,dN1]
        elif kind == 'runout':
            B[0,:3] = [1,-2,1]
            b[-1,-3:] = [1,-2,1]
        elif kind == 'parabolic':
            B[0,:2] = [1,-1]
            B[-1,-2:] = [-1,1]
        elif kind == 'periodic':
            raise NotImplementedError
        elif kind == 'symmetric':
            raise NotImplementedError
        else:
            B[0,:2] = [2*d1,d1]
            B[-1,-2:] = [dN,2*dN]

        # Set up RHS (b)
        b = np.empty((Np1,)+yk.shape[1:])
        dyk = (yk[1:]-yk[:-1])*1.0
        if kind in ['not-a-knot', 'runout', 'parabolic']:
            b[0] = b[-1] = 0.0
        elif kind == 'periodic':
            raise NotImplementedError
        elif kind == 'symmetric':
            raise NotImplementedError
        else:
            b[0] = (dyk[0]/d1 - sl_0)
            b[-1] = -(dyk[-1]/dN - sl_N)
        b[1:-1,...] = (dyk[1:]/alpha[1:]-dyk[:-1]/alpha[:-1])
        b *= 6.0
        return B, b, None, None
    else:
        raise ValueError, "%s not supported" % kind
        

def cspeval((mk,xk,yk),xnew):
    """Evaluate a cubic-spline representation of the points (xk,yk)
    at the new values xnew.  The mk values are the second derivatives at xk
    The xk vector must be sorted.
    
    More than one curve can be represented using 2-d arrays for mk and yk.
    However, the last dimension must have the same shape as the 1-d array xk.
    The first-dimension will be considered the interpolating dimension. 
    """
    indxs = np.searchsorted(xk, xnew)
    indxs[indxs==0] = 1
    indxsm1 = indxs-1
    xkm1 = xk[indxsm1]
    xkvals = xk[indxs]
    dm1 = xnew - xkm1
    d = xkvals - xnew
    mk0 = mk[indxs]
    mkm1 = mk[indxsm1]
    dk = xkvals-xkm1
    val = (mk0*dm1**3. + mkm1*d**3.)/(6*dk)
    val += (yk[indxsm1]/dk - mkm1*dk/6.)*d
    val += (yk[indxs]/dk - mk0*dk/6.)*dm1
    return val

def csp2pp(mk,xk,yk):
    """Return an N-d array providing the piece-wise polynomial form.

    mk - second derivative at the knots
    xk - knot-points
    yk - values of the curve at the knots
    
    The first 2 dimensions are the polynomial for a particular
    curve.  The first dimension provides the coefficients for the
    polynomial and the second dimension provides the different pieces
    """
    dk = xk[1:] - xk[:-1]
    temp1 = mk[1:] - mk[:-1]
    temp2 = mk[1:]*xk[:-1]-mk[:-1]*xk[1:]
    c3 = temp1/(6*dk)
    c2 = -temp2/(2*dk)
    c1 = (mk[1:]*xk[:-1]**2 - mk[:-1]*xk[1:]**2)/(2*dk)
    c1 -= temp1*dk/6.
    c1 += (yk[1:]-yk[:-1])/dk
    c0 = (mk[:-1]*xk[1:]**3 - mk[1:]*xk[:-1]**3)/(6*dk)
    c0 += temp2*dk/6.
    c0 += (yk[:-1]*xk[1:] - yk[1:]*xk[:-1])/dk
    return np.array([c3,c2,c1,c0])    

def ppval(pp, xk, xnew):
    """Compute a piece-wise polynomial defined by the array of
    coefficents pp and the break-points xk on the grid xnew
    """
    indxs = numpy.searchsorted(xk, xnew)-1
    indxs[indxs<0]=0
    return array([numpy.polyval(pp[:,k],xnew[i]) for i,k in enumerate(indxs)])

def csprep(xk,yk,kind='not-a-knot',conds=None):
    """Return a (Spp,xk,yk) representation of a cubic spline given
    data-points

    yk can be an N-d array to represent more than one curve, through
    the same xk points. The first dimension is assumed to be the
    interpolating dimenion.

    kind can be 'natural', 'second', 'clamped', 'endslope', 'periodic',
                'symmetric', 'parabolic', 'not-a-knot', 'runout'

    for 'second', and 'clamped' conditions can be given which should
    be the desired second and first derivatives at
    the end-points, respectively. 
    """
    yk = np.asanyarray(yk)
    N = yk.shape[0]-1
    B,b,first,last = _get_cspline_Bb(xk, yk, kind, conds)
    mk = np.dual.solve(B,b)
    if first is not None:
        mk = np.concatenate((first, mk), axis=0)
    if last is not None:
        mk = np.concatenate((mk, last), axis=0)
    return mk, xk, yk

def cspline(xk,yk,xnew,kind='not-a-knot',conds=None):
    return cspeval(csprep(xk,yk,kind=kind,conds=conds),xnew)
