

""" Classes for interpolating values.
"""

__all__ = ['interp1d', 'interp2d', 'spline', 'spleval', 'splmake', 'spltopp',
           'ppform', 'lagrange']

from numpy import shape, sometrue, array, transpose, searchsorted, \
                  ones, logical_or, atleast_1d, atleast_2d, meshgrid, ravel, \
                  dot, poly1d, asarray, intp
import numpy as np
import scipy.special as spec
import math

import fitpack
import _fitpack

def reduce_sometrue(a):
    all = a
    while len(shape(all)) > 1:
        all = sometrue(all,axis=0)
    return all

def lagrange(x, w):
    """
    Return a Lagrange interpolating polynomial.

    Given two 1-D arrays `x` and `w,` returns the Lagrange interpolating
    polynomial through the points ``(x, w)``.

    Warning: This implementation is numerically unstable. Do not expect to
    be able to use more than about 20 points even if they are chosen optimally.

    Parameters
    ----------
    x : array_like
        `x` represents the x-coordinates of a set of datapoints.
    w : array_like
        `w` represents the y-coordinates of a set of datapoints, i.e. f(`x`).

    """
    M = len(x)
    p = poly1d(0.0)
    for j in xrange(M):
        pt = poly1d(w[j])
        for k in xrange(M):
            if k == j: continue
            fac = x[j]-x[k]
            pt *= poly1d([1.0,-x[k]])/fac
        p += pt
    return p


# !! Need to find argument for keeping initialize.  If it isn't
# !! found, get rid of it!

class interp2d(object):
    """
    interp2d(x, y, z, kind='linear', copy=True, bounds_error=False,
             fill_value=nan)

    Interpolate over a 2-D grid.

    `x`, `y` and `z` are arrays of values used to approximate some function
    f: ``z = f(x, y)``. This class returns a function whose call method uses
    spline interpolation to find the value of new points.

    Methods
    -------
    __call__

    Parameters
    ----------
    x, y : 1-D ndarrays
        Arrays defining the data point coordinates.

        If the points lie on a regular grid, `x` can specify the column
        coordinates and `y` the row coordinates, for example::

          >>> x = [0,1,2];  y = [0,3]; z = [[1,2,3], [4,5,6]]

        Otherwise, x and y must specify the full coordinates for each point,
        for example::

          >>> x = [0,1,2,0,1,2];  y = [0,0,0,3,3,3]; z = [1,2,3,4,5,6]

        If `x` and `y` are multi-dimensional, they are flattened before use.

    z : 1-D ndarray
        The values of the function to interpolate at the data points. If
        `z` is a multi-dimensional array, it is flattened before use.
    kind : {'linear', 'cubic', 'quintic'}, optional
        The kind of spline interpolation to use. Default is 'linear'.
    copy : bool, optional
        If True, then data is copied, otherwise only a reference is held.
    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, an error is raised.
        If False, then `fill_value` is used.
    fill_value : number, optional
        If provided, the value to use for points outside of the
        interpolation domain. Defaults to NaN.

    See Also
    --------
    bisplrep, bisplev
        Spline interpolation based on FITPACK
    BivariateSpline : a more recent wrapper of the FITPACK routines
    interp1d

    Notes
    -----
    The minimum number of data points required along the interpolation
    axis is ``(k+1)**2``, with k=1 for linear, k=3 for cubic and k=5 for
    quintic interpolation.

    The interpolator is constructed by `bisplrep`, with a smoothing factor
    of 0. If more control over smoothing is needed, `bisplrep` should be
    used directly.

    Examples
    --------
    Construct a 2-D grid and interpolate on it:

    >>> x = np.arange(-5.01, 5.01, 0.25)
    >>> y = np.arange(-5.01, 5.01, 0.25)
    >>> xx, yy = np.meshgrid(x, y)
    >>> z = np.sin(xx**2+yy**2)
    >>> f = sp.interpolate.interp2d(x, y, z, kind='cubic')

    Now use the obtained interpolation function and plot the result:

    >>> xnew = np.arange(-5.01, 5.01, 1e-2)
    >>> ynew = np.arange(-5.01, 5.01, 1e-2)
    >>> znew = f(xnew, ynew)
    >>> plt.plot(x, z[:, 0], 'ro-', xnew, znew[:, 0], 'b-')

    """

    def __init__(self, x, y, z, kind='linear', copy=True, bounds_error=False,
                 fill_value=np.nan):
        self.x, self.y, self.z = map(ravel, map(asarray, [x, y, z]))

        if len(self.z) == len(self.x) * len(self.y):
            self.x, self.y = meshgrid(x,y)
            self.x, self.y = map(ravel, [self.x, self.y])
        if len(self.x) != len(self.y):
            raise ValueError("x and y must have equal lengths")
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
        """Interpolate the function.

        Parameters
        ----------
        x : 1D array
            x-coordinates of the mesh on which to interpolate.
        y : 1D array
            y-coordinates of the mesh on which to interpolate.
        dx : int >= 0, < kx
            Order of partial derivatives in x.
        dy : int >= 0, < ky
            Order of partial derivatives in y.

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
    """
    interp1d(x, y, kind='linear', axis=-1, copy=True, bounds_error=True,
             fill_value=np.nan)

    Interpolate a 1-D function.

    `x` and `y` are arrays of values used to approximate some function f:
    ``y = f(x)``.  This class returns a function whose call method uses
    interpolation to find the value of new points.

    Parameters
    ----------
    x : array_like
        A 1-D array of monotonically increasing real values.
    y : array_like
        A N-D array of real values. The length of `y` along the interpolation
        axis must be equal to the length of `x`.
    kind : str or int, optional
        Specifies the kind of interpolation as a string
        ('linear','nearest', 'zero', 'slinear', 'quadratic, 'cubic')
        or as an integer specifying the order of the spline interpolator
        to use. Default is 'linear'.
    axis : int, optional
        Specifies the axis of `y` along which to interpolate.
        Interpolation defaults to the last axis of `y`.
    copy : bool, optional
        If True, the class makes internal copies of x and y.
        If False, references to `x` and `y` are used. The default is to copy.
    bounds_error : bool, optional
        If True, an error is thrown any time interpolation is attempted on
        a value outside of the range of x (where extrapolation is
        necessary). If False, out of bounds values are assigned `fill_value`.
        By default, an error is raised.
    fill_value : float, optional
        If provided, then this value will be used to fill in for requested
        points outside of the data range. If not provided, then the default
        is NaN.

    See Also
    --------
    UnivariateSpline : A more recent wrapper of the FITPACK routines.
    splrep, splev
        Spline interpolation based on FITPACK.
    interp2d

    Examples
    --------
    >>> import scipy.interpolate
    >>> x = np.arange(0, 10)
    >>> y = np.exp(-x/3.0)
    >>> f = sp.interpolate.interp1d(x, y)

    >>> xnew = np.arange(0,9, 0.1)
    >>> ynew = f(xnew)   # use interpolation function returned by `interp1d`
    >>> plt.plot(x, y, 'o', xnew, ynew, '-')

    """

    def __init__(self, x, y, kind='linear', axis=-1,
                 copy=True, bounds_error=True, fill_value=np.nan):
        """ Initialize a 1D linear interpolation class."""

        self.copy = copy
        self.bounds_error = bounds_error
        self.fill_value = fill_value

        if kind in ['zero', 'slinear', 'quadratic', 'cubic']:
            order = {'nearest':0, 'zero':0,'slinear':1,
                     'quadratic':2, 'cubic':3}[kind]
            kind = 'spline'
        elif isinstance(kind, int):
            order = kind
            kind = 'spline'
        elif kind not in ('linear', 'nearest'):
            raise NotImplementedError("%s is unsupported: Use fitpack "\
                                      "routines for other types." % kind)
        x = array(x, copy=self.copy)
        y = array(y, copy=self.copy)

        if x.ndim != 1:
            raise ValueError("the x array must have exactly one dimension.")
        if y.ndim == 0:
            raise ValueError("the y array must have at least one dimension.")

        # Force-cast y to a floating-point type, if it's not yet one
        if not issubclass(y.dtype.type, np.inexact):
            y = y.astype(np.float_)

        # Normalize the axis to ensure that it is positive.
        self.axis = axis % len(y.shape)
        self._kind = kind

        if kind in ('linear', 'nearest'):
            # Make a "view" of the y array that is rotated to the interpolation
            # axis.
            axes = range(y.ndim)
            del axes[self.axis]
            axes.append(self.axis)
            oriented_y = y.transpose(axes)
            minval = 2
            len_y = oriented_y.shape[-1]
            if kind == 'linear':
                self._call = self._call_linear
            elif kind == 'nearest':
                self.x_bds = (x[1:] + x[:-1]) / 2.0
                self._call = self._call_nearest
        else:
            axes = range(y.ndim)
            del axes[self.axis]
            axes.insert(0, self.axis)
            oriented_y = y.transpose(axes)
            minval = order + 1
            len_y = oriented_y.shape[0]
            self._call = self._call_spline
            self._spline = splmake(x,oriented_y,order=order)

        len_x = len(x)
        if len_x != len_y:
            raise ValueError("x and y arrays must be equal in length along "
                             "interpolation axis.")
        if len_x < minval:
            raise ValueError("x and y arrays must have at "
                             "least %d entries" % minval)
        self.x = x
        self.y = oriented_y

    def _call_linear(self, x_new):

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

        return y_new

    def _call_nearest(self, x_new):
        """ Find nearest neighbour interpolated y_new = f(x_new)."""

        # 2. Find where in the averaged data the values to interpolate
        #    would be inserted.
        #    Note: use side='left' (right) to searchsorted() to define the
        #    halfway point to be nearest to the left (right) neighbour
        x_new_indices = searchsorted(self.x_bds, x_new, side='left')

        # 3. Clip x_new_indices so that they are within the range of x indices.
        x_new_indices = x_new_indices.clip(0,  len(self.x)-1).astype(intp)

        # 4. Calculate the actual value for each entry in x_new.
        y_new = self.y[..., x_new_indices]

        return y_new

    def _call_spline(self, x_new):
        x_new =np.asarray(x_new)
        result = spleval(self._spline,x_new.ravel())
        return result.reshape(x_new.shape+result.shape[1:])

    def __call__(self, x_new):
        """Find interpolated y_new = f(x_new).

        Parameters
        ----------
        x_new : number or array
            New independent variable(s).

        Returns
        -------
        y_new : ndarray
            Interpolated value(s) corresponding to x_new.

        """

        # 1. Handle values in x_new that are outside of x.  Throw error,
        #    or return a list of mask array indicating the outofbounds values.
        #    The behavior is set by the bounds_error variable.
        x_new = asarray(x_new)
        out_of_bounds = self._check_bounds(x_new)

        y_new = self._call(x_new)

        # Rotate the values of y_new back so that they correspond to the
        # correct x_new values. For N-D x_new, take the last (for linear)
        # or first (for other splines) N axes
        # from y_new and insert them where self.axis was in the list of axes.
        nx = x_new.ndim
        ny = y_new.ndim

        # 6. Fill any values that were out of bounds with fill_value.
        # and
        # 7. Rotate the values back to their proper place.

        if nx == 0:
            # special case: x is a scalar
            if out_of_bounds:
                if ny == 0:
                    return asarray(self.fill_value)
                else:
                    y_new[...] = self.fill_value
            return asarray(y_new)
        elif self._kind in ('linear', 'nearest'):
            y_new[..., out_of_bounds] = self.fill_value
            axes = range(ny - nx)
            axes[self.axis:self.axis] = range(ny - nx, ny)
            return y_new.transpose(axes)
        else:
            y_new[out_of_bounds] = self.fill_value
            axes = range(nx, ny)
            axes[self.axis:self.axis] = range(nx)
            return y_new.transpose(axes)

    def _check_bounds(self, x_new):
        """Check the inputs for being in the bounds of the interpolated data.

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

class ppform(object):
    """The ppform of the piecewise polynomials is given in terms of coefficients
    and breaks.  The polynomial in the ith interval is
    x_{i} <= x < x_{i+1}

    S_i = sum(coefs[m,i]*(x-breaks[i])^(k-m), m=0..k)
    where k is the degree of the polynomial.
    """
    def __init__(self, coeffs, breaks, fill=0.0, sort=False):
        self.coeffs = np.asarray(coeffs)
        if sort:
            self.breaks = np.sort(breaks)
        else:
            self.breaks = np.asarray(breaks)
        self.K = self.coeffs.shape[0]
        self.fill = fill
        self.a = self.breaks[0]
        self.b = self.breaks[-1]

    def __call__(self, xnew):
        saveshape = np.shape(xnew)
        xnew = np.ravel(xnew)
        res = np.empty_like(xnew)
        mask = (xnew >= self.a) & (xnew <= self.b)
        res[~mask] = self.fill
        xx = xnew.compress(mask)
        indxs = np.searchsorted(self.breaks, xx)-1
        indxs = indxs.clip(0,len(self.breaks))
        pp = self.coeffs
        diff = xx - self.breaks.take(indxs)
        V = np.vander(diff,N=self.K)
        # values = np.diag(dot(V,pp[:,indxs]))
        values = array([dot(V[k,:],pp[:,indxs[k]]) for k in xrange(len(xx))])
        res[mask] = values
        res.shape = saveshape
        return res

    def fromspline(cls, xk, cvals, order, fill=0.0):
        N = len(xk)-1
        sivals = np.empty((order+1,N), dtype=float)
        for m in xrange(order,-1,-1):
            fact = spec.gamma(m+1)
            res = _fitpack._bspleval(xk[:-1], xk, cvals, order, m)
            res /= fact
            sivals[order-m,:] = res
        return cls(sivals, xk, fill=fill)
    fromspline = classmethod(fromspline)


def _dot0(a, b):
    """Similar to numpy.dot, but sum over last axis of a and 1st axis of b"""
    if b.ndim <= 2:
        return dot(a, b)
    else:
        axes = range(b.ndim)
        axes.insert(-1, 0)
        axes.pop(0)
        return dot(a, b.transpose(axes))

def _find_smoothest(xk, yk, order, conds=None, B=None):
    # construct Bmatrix, and Jmatrix
    # e = J*c
    # minimize norm(e,2) given B*c=yk
    # if desired B can be given
    # conds is ignored
    N = len(xk)-1
    K = order
    if B is None:
        B = _fitpack._bsplmat(order, xk)
    J = _fitpack._bspldismat(order, xk)
    u,s,vh = np.dual.svd(B)
    ind = K-1
    V2 = vh[-ind:,:].T
    V1 = vh[:-ind,:].T
    A = dot(J.T,J)
    tmp = dot(V2.T,A)
    Q = dot(tmp,V2)
    p = np.dual.solve(Q,tmp)
    tmp = dot(V2,p)
    tmp = np.eye(N+K) - tmp
    tmp = dot(tmp,V1)
    tmp = dot(tmp,np.diag(1.0/s))
    tmp = dot(tmp,u.T)
    return _dot0(tmp, yk)

def _setdiag(a, k, v):
    if not a.ndim == 2:
        raise ValueError("Input array should be 2-D.")
    M,N = a.shape
    if k > 0:
        start = k
        num = N-k
    else:
        num = M+k
        start = abs(k)*N
    end = start + num*(N+1)-1
    a.flat[start:end:(N+1)] = v

# Return the spline that minimizes the dis-continuity of the
# "order-th" derivative; for order >= 2.

def _find_smoothest2(xk, yk):
    N = len(xk)-1
    Np1 = N+1
    # find pseudo-inverse of B directly.
    Bd = np.empty((Np1,N))
    for k in range(-N,N):
        if (k<0):
            l = np.arange(-k,Np1)
            v = (l+k+1)
            if ((k+1) % 2):
                v = -v
        else:
            l = np.arange(k,N)
            v = N-l
            if ((k % 2)):
                v = -v
        _setdiag(Bd,k,v)
    Bd /= (Np1)
    V2 = np.ones((Np1,))
    V2[1::2] = -1
    V2 /= math.sqrt(Np1)
    dk = np.diff(xk)
    b = 2*np.diff(yk, axis=0)/dk
    J = np.zeros((N-1,N+1))
    idk = 1.0/dk
    _setdiag(J,0,idk[:-1])
    _setdiag(J,1,-idk[1:]-idk[:-1])
    _setdiag(J,2,idk[1:])
    A = dot(J.T,J)
    val = dot(V2,dot(A,V2))
    res1 = dot(np.outer(V2,V2)/val,A)
    mk = dot(np.eye(Np1)-res1, _dot0(Bd,b))
    return mk

def _get_spline2_Bb(xk, yk, kind, conds):
    Np1 = len(xk)
    dk = xk[1:]-xk[:-1]
    if kind == 'not-a-knot':
        # use banded-solver
        nlu = (1,1)
        B = ones((3,Np1))
        alpha = 2*(yk[1:]-yk[:-1])/dk
        zrs = np.zeros((1,)+yk.shape[1:])
        row = (Np1-1)//2
        b = np.concatenate((alpha[:row],zrs,alpha[row:]),axis=0)
        B[0,row+2:] = 0
        B[2,:(row-1)] = 0
        B[0,row+1] = dk[row-1]
        B[1,row] = -dk[row]-dk[row-1]
        B[2,row-1] = dk[row]
        return B, b, None, nlu
    else:
        raise NotImplementedError("quadratic %s is not available" % kind)

def _get_spline3_Bb(xk, yk, kind, conds):
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
        # use banded solver
        beta = 2*(xk[2:]-xk[:-2])
        alpha = xk[1:]-xk[:-1]
        nlu = (1,1)
        B = np.empty((3,Np1-2))
        B[0,1:] = alpha[2:]
        B[1,:] = beta
        B[2,:-1] = alpha[1:-1]
        dyk = yk[1:]-yk[:-1]
        b = (dyk[1:]/alpha[1:] - dyk[:-1]/alpha[:-1])
        b *= 6
        b[0] -= m0
        b[-1] -= mN

        def append_func(mk):
            # put m0 and mN into the correct shape for
            #  concatenation
            ma = array(m0,copy=0,ndmin=yk.ndim)
            mb = array(mN,copy=0,ndmin=yk.ndim)
            if ma.shape[1:] != yk.shape[1:]:
                ma = ma*(ones(yk.shape[1:])[np.newaxis,...])
            if mb.shape[1:] != yk.shape[1:]:
                mb = mb*(ones(yk.shape[1:])[np.newaxis,...])
            mk = np.concatenate((ma,mk),axis=0)
            mk = np.concatenate((mk,mb),axis=0)
            return mk

        return B, b, append_func, nlu


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
            B[-1,-3:] = [1,-2,1]
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
        raise ValueError("%s not supported" % kind)

# conds is a tuple of an array and a vector
#  giving the left-hand and the right-hand side
#  of the additional equations to add to B
def _find_user(xk, yk, order, conds, B):
    lh = conds[0]
    rh = conds[1]
    B = np.concatenate((B,lh),axis=0)
    w = np.concatenate((yk,rh),axis=0)
    M,N = B.shape
    if (M>N):
        raise ValueError("over-specification of conditions")
    elif (M<N):
        return _find_smoothest(xk, yk, order, None, B)
    else:
        return np.dual.solve(B, w)

# If conds is None, then use the not_a_knot condition
#  at K-1 farthest separated points in the interval
def _find_not_a_knot(xk, yk, order, conds, B):
    raise NotImplementedError
    return _find_user(xk, yk, order, conds, B)

# If conds is None, then ensure zero-valued second
#  derivative at K-1 farthest separated points
def _find_natural(xk, yk, order, conds, B):
    raise NotImplementedError
    return _find_user(xk, yk, order, conds, B)

# If conds is None, then ensure zero-valued first
#  derivative at K-1 farthest separated points
def _find_clamped(xk, yk, order, conds, B):
    raise NotImplementedError
    return _find_user(xk, yk, order, conds, B)

def _find_fixed(xk, yk, order, conds, B):
    raise NotImplementedError
    return _find_user(xk, yk, order, conds, B)

# If conds is None, then use coefficient periodicity
# If conds is 'function' then use function periodicity
def _find_periodic(xk, yk, order, conds, B):
    raise NotImplementedError
    return _find_user(xk, yk, order, conds, B)

# Doesn't use conds
def _find_symmetric(xk, yk, order, conds, B):
    raise NotImplementedError
    return _find_user(xk, yk, order, conds, B)

# conds is a dictionary with multiple values
def _find_mixed(xk, yk, order, conds, B):
    raise NotImplementedError
    return _find_user(xk, yk, order, conds, B)


def splmake(xk,yk,order=3,kind='smoothest',conds=None):
    """Return a (xk, cvals, k) representation of a spline given
    data-points where the (internal) knots are at the data-points.

    yk can be an N-d array to represent more than one curve, through
    the same xk points. The first dimension is assumed to be the
    interpolating dimension.

    kind can be 'smoothest', 'not_a_knot', 'fixed',
                'clamped', 'natural', 'periodic', 'symmetric',
                'user', 'mixed'

                it is ignored if order < 2
    """
    yk = np.asanyarray(yk)
    N = yk.shape[0]-1

    order = int(order)
    if order < 0:
        raise ValueError("order must not be negative")
    if order == 0:
        return xk, yk[:-1], order
    elif order == 1:
        return xk, yk, order

    try:
        func = eval('_find_%s' % kind)
    except:
        raise NotImplementedError

    # the constraint matrix
    B = _fitpack._bsplmat(order, xk)
    coefs = func(xk, yk, order, conds, B)
    return xk, coefs, order

def spleval((xj,cvals,k),xnew,deriv=0):
    """Evaluate a fixed spline represented by the given tuple at the new
    x-values. The xj values are the interior knot points.  The approximation
    region is xj[0] to xj[-1].  If N+1 is the length of xj, then cvals should
    have length N+k where k is the order of the spline.

    Internally, an additional k-1 knot points are added on either side of
    the spline.

    If cvals represents more than one curve (cvals.ndim > 1) and/or xnew is
    N-d, then the result is xnew.shape + cvals.shape[1:] providing the
    interpolation of multiple curves.
    """
    oldshape = np.shape(xnew)
    xx = np.ravel(xnew)
    sh = cvals.shape[1:]
    res = np.empty(xx.shape + sh, dtype=cvals.dtype)
    for index in np.ndindex(*sh):
        sl = (slice(None),)+index
        if issubclass(cvals.dtype.type, np.complexfloating):
            res[sl].real = _fitpack._bspleval(xx,xj,cvals.real[sl],k,deriv)
            res[sl].imag = _fitpack._bspleval(xx,xj,cvals.imag[sl],k,deriv)
        else:
            res[sl] = _fitpack._bspleval(xx,xj,cvals[sl],k,deriv)
    res.shape = oldshape + sh
    return res

def spltopp(xk,cvals,k):
    """Return a piece-wise polynomial object from a fixed-spline tuple.
    """
    return ppform.fromspline(xk, cvals, k)

def spline(xk,yk,xnew,order=3,kind='smoothest',conds=None):
    """Interpolate a curve (xk,yk) at points xnew using a spline fit.
    """
    return spleval(splmake(xk,yk,order=order,kind=kind,conds=conds),xnew)
