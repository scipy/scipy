""" Classes for interpolating values.
"""
from __future__ import division, print_function, absolute_import

__all__ = ['interp1d', 'interp2d', 'spline', 'spleval', 'splmake', 'spltopp',
           'ppform', 'lagrange', 'PPoly', 'BPoly', 'RegularGridInterpolator',
           'interpn']

import itertools

from numpy import (shape, sometrue, array, transpose, searchsorted,
                   ones, logical_or, atleast_1d, atleast_2d, ravel,
                   dot, poly1d, asarray, intp)
import numpy as np
import scipy.linalg
import scipy.special as spec
from scipy.misc import comb
import math
import warnings
import functools
import operator

from scipy.lib.six import xrange, integer_types

from . import fitpack
from . import dfitpack
from . import _fitpack
from .polyint import _Interpolator1D
from . import _ppoly
from .fitpack2 import RectBivariateSpline
from .interpnd import _ndim_coords_from_arrays


def reduce_sometrue(a):
    all = a
    while len(shape(all)) > 1:
        all = sometrue(all, axis=0)
    return all


def prod(x):
    """Product of a list of numbers; ~40x faster vs np.prod for Python tuples"""
    if len(x) == 0:
        return 1
    return functools.reduce(operator.mul, x)


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

    Returns
    -------
    lagrange : numpy.poly1d instance
        The Lagrange interpolating polynomial.

    """
    M = len(x)
    p = poly1d(0.0)
    for j in xrange(M):
        pt = poly1d(w[j])
        for k in xrange(M):
            if k == j:
                continue
            fac = x[j]-x[k]
            pt *= poly1d([1.0, -x[k]])/fac
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

    If `x` and `y` represent a regular grid, consider using
    RectBivariateSpline.

    Methods
    -------
    __call__

    Parameters
    ----------
    x, y : array_like
        Arrays defining the data point coordinates.

        If the points lie on a regular grid, `x` can specify the column
        coordinates and `y` the row coordinates, for example::

          >>> x = [0,1,2];  y = [0,3]; z = [[1,2,3], [4,5,6]]

        Otherwise, `x` and `y` must specify the full coordinates for each
        point, for example::

          >>> x = [0,1,2,0,1,2];  y = [0,0,0,3,3,3]; z = [1,2,3,4,5,6]

        If `x` and `y` are multi-dimensional, they are flattened before use.
    z : array_like
        The values of the function to interpolate at the data points. If
        `z` is a multi-dimensional array, it is flattened before use.  The
        length of a flattened `z` array is either
        len(`x`)*len(`y`) if `x` and `y` specify the column and row coordinates
        or ``len(z) == len(x) == len(y)`` if `x` and `y` specify coordinates
        for each point.
    kind : {'linear', 'cubic', 'quintic'}, optional
        The kind of spline interpolation to use. Default is 'linear'.
    copy : bool, optional
        If True, the class makes internal copies of x, y and z.
        If False, references may be used. The default is to copy.
    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data (x,y), a ValueError is raised.
        If False, then `fill_value` is used.
    fill_value : number, optional
        If provided, the value to use for points outside of the
        interpolation domain. If omitted (None), values outside
        the domain are extrapolated.

    Returns
    -------
    values_x : ndarray, shape xi.shape[:-1] + values.shape[ndim:]
        Interpolated values at input coordinates.

    See Also
    --------
    RectBivariateSpline :
        Much faster 2D interpolation if your input data is on a grid
    bisplrep, bisplev :
        Spline interpolation based on FITPACK
    BivariateSpline : a more recent wrapper of the FITPACK routines
    interp1d : one dimension version of this function

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

    >>> from scipy import interpolate
    >>> x = np.arange(-5.01, 5.01, 0.25)
    >>> y = np.arange(-5.01, 5.01, 0.25)
    >>> xx, yy = np.meshgrid(x, y)
    >>> z = np.sin(xx**2+yy**2)
    >>> f = interpolate.interp2d(x, y, z, kind='cubic')

    Now use the obtained interpolation function and plot the result:

    >>> xnew = np.arange(-5.01, 5.01, 1e-2)
    >>> ynew = np.arange(-5.01, 5.01, 1e-2)
    >>> znew = f(xnew, ynew)
    >>> plt.plot(x, z[0, :], 'ro-', xnew, znew[0, :], 'b-')
    >>> plt.show()

    """

    def __init__(self, x, y, z, kind='linear', copy=True, bounds_error=False,
                 fill_value=None):
        x = ravel(x)
        y = ravel(y)
        z = asarray(z)

        rectangular_grid = (z.size == len(x) * len(y))
        if rectangular_grid:
            if z.ndim == 2:
                if z.shape != (len(y), len(x)):
                    raise ValueError("When on a regular grid with x.size = m "
                                     "and y.size = n, if z.ndim == 2, then z "
                                     "must have shape (n, m)")
            if not np.all(x[1:] >= x[:-1]):
                j = np.argsort(x)
                x = x[j]
                z = z[:, j]
            if not np.all(y[1:] >= y[:-1]):
                j = np.argsort(y)
                y = y[j]
                z = z[j, :]
            z = ravel(z.T)
        else:
            z = ravel(z)
            if len(x) != len(y):
                raise ValueError(
                    "x and y must have equal lengths for non rectangular grid")
            if len(z) != len(x):
                raise ValueError(
                    "Invalid length for input z for non rectangular grid")

        try:
            kx = ky = {'linear': 1,
                       'cubic': 3,
                       'quintic': 5}[kind]
        except KeyError:
            raise ValueError("Unsupported interpolation type.")

        if not rectangular_grid:
            # TODO: surfit is really not meant for interpolation!
            self.tck = fitpack.bisplrep(x, y, z, kx=kx, ky=ky, s=0.0)
        else:
            nx, tx, ny, ty, c, fp, ier = dfitpack.regrid_smth(
                x, y, z, None, None, None, None,
                kx=kx, ky=ky, s=0.0)
            self.tck = (tx[:nx], ty[:ny], c[:(nx - kx - 1) * (ny - ky - 1)],
                        kx, ky)

        self.bounds_error = bounds_error
        self.fill_value = fill_value
        self.x, self.y, self.z = [array(a, copy=copy) for a in (x, y, z)]

        self.x_min, self.x_max = np.amin(x), np.amax(x)
        self.y_min, self.y_max = np.amin(y), np.amax(y)

    def __call__(self, x, y, dx=0, dy=0, assume_sorted=False):
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
        assume_sorted : bool, optional
            If False, values of `x` and `y` can be in any order and they are
            sorted first.
            If True, `x` and `y` have to be arrays of monotonically
            increasing values.

        Returns
        -------
        z : 2D array with shape (len(y), len(x))
            The interpolated values.

        """

        x = atleast_1d(x)
        y = atleast_1d(y)

        if x.ndim != 1 or y.ndim != 1:
            raise ValueError("x and y should both be 1-D arrays")

        if not assume_sorted:
            x = np.sort(x)
            y = np.sort(y)

        if self.bounds_error or self.fill_value is not None:
            out_of_bounds_x = (x < self.x_min) | (x > self.x_max)
            out_of_bounds_y = (y < self.y_min) | (y > self.y_max)

            any_out_of_bounds_x = np.any(out_of_bounds_x)
            any_out_of_bounds_y = np.any(out_of_bounds_y)

        if self.bounds_error and (any_out_of_bounds_x or any_out_of_bounds_y):
            raise ValueError("Values out of range; x must be in %r, y in %r"
                             % ((self.x_min, self.x_max),
                                (self.y_min, self.y_max)))

        z = fitpack.bisplev(x, y, self.tck, dx, dy)
        z = atleast_2d(z)
        z = transpose(z)

        if self.fill_value is not None:
            if any_out_of_bounds_x:
                z[:, out_of_bounds_x] = self.fill_value
            if any_out_of_bounds_y:
                z[out_of_bounds_y, :] = self.fill_value

        if len(z) == 1:
            z = z[0]
        return array(z)


class interp1d(_Interpolator1D):
    """
    Interpolate a 1-D function.

    `x` and `y` are arrays of values used to approximate some function f:
    ``y = f(x)``.  This class returns a function whose call method uses
    interpolation to find the value of new points.

    Parameters
    ----------
    x : (N,) array_like
        A 1-D array of real values.
    y : (...,N,...) array_like
        A N-D array of real values. The length of `y` along the interpolation
        axis must be equal to the length of `x`.
    kind : str or int, optional
        Specifies the kind of interpolation as a string
        ('linear', 'nearest', 'zero', 'slinear', 'quadratic, 'cubic'
        where 'slinear', 'quadratic' and 'cubic' refer to a spline
        interpolation of first, second or third order) or as an integer
        specifying the order of the spline interpolator to use.
        Default is 'linear'.
    axis : int, optional
        Specifies the axis of `y` along which to interpolate.
        Interpolation defaults to the last axis of `y`.
    copy : bool, optional
        If True, the class makes internal copies of x and y.
        If False, references to `x` and `y` are used. The default is to copy.
    bounds_error : bool, optional
        If True, a ValueError is raised any time interpolation is attempted on
        a value outside of the range of x (where extrapolation is
        necessary). If False, out of bounds values are assigned `fill_value`.
        By default, an error is raised.
    fill_value : float, optional
        If provided, then this value will be used to fill in for requested
        points outside of the data range. If not provided, then the default
        is NaN.
    assume_sorted : bool, optional
        If False, values of `x` can be in any order and they are sorted first.
        If True, `x` has to be an array of monotonically increasing values.

    Methods
    -------
    __call__

    See Also
    --------
    splrep, splev
        Spline interpolation/smoothing based on FITPACK.
    UnivariateSpline : An object-oriented wrapper of the FITPACK routines.
    interp2d : 2-D interpolation

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy import interpolate
    >>> x = np.arange(0, 10)
    >>> y = np.exp(-x/3.0)
    >>> f = interpolate.interp1d(x, y)

    >>> xnew = np.arange(0, 9, 0.1)
    >>> ynew = f(xnew)   # use interpolation function returned by `interp1d`
    >>> plt.plot(x, y, 'o', xnew, ynew, '-')
    >>> plt.show()

    """

    def __init__(self, x, y, kind='linear', axis=-1,
                 copy=True, bounds_error=True, fill_value=np.nan,
                 assume_sorted=False):
        """ Initialize a 1D linear interpolation class."""
        _Interpolator1D.__init__(self, x, y, axis=axis)

        self.copy = copy
        self.bounds_error = bounds_error
        self.fill_value = fill_value

        if kind in ['zero', 'slinear', 'quadratic', 'cubic']:
            order = {'nearest': 0, 'zero': 0,'slinear': 1,
                     'quadratic': 2, 'cubic': 3}[kind]
            kind = 'spline'
        elif isinstance(kind, int):
            order = kind
            kind = 'spline'
        elif kind not in ('linear', 'nearest'):
            raise NotImplementedError("%s is unsupported: Use fitpack "
                                      "routines for other types." % kind)
        x = array(x, copy=self.copy)
        y = array(y, copy=self.copy)

        if not assume_sorted:
            ind = np.argsort(x)
            x = x[ind]
            y = np.take(y, ind, axis=axis)

        if x.ndim != 1:
            raise ValueError("the x array must have exactly one dimension.")
        if y.ndim == 0:
            raise ValueError("the y array must have at least one dimension.")

        # Force-cast y to a floating-point type, if it's not yet one
        if not issubclass(y.dtype.type, np.inexact):
            y = y.astype(np.float_)

        # Backward compatibility
        self.axis = axis % y.ndim

        # Interpolation goes internally along the first axis
        self.y = y
        y = self._reshape_yi(y)

        # Adjust to interpolation kind; store reference to *unbound*
        # interpolation methods, in order to avoid circular references to self
        # stored in the bound instance methods, and therefore delayed garbage
        # collection.  See: http://docs.python.org/2/reference/datamodel.html
        if kind in ('linear', 'nearest'):
            # Make a "view" of the y array that is rotated to the interpolation
            # axis.
            minval = 2
            if kind == 'nearest':
                self.x_bds = (x[1:] + x[:-1]) / 2.0
                self._call = self.__class__._call_nearest
            else:
                self._call = self.__class__._call_linear
        else:
            minval = order + 1
            self._spline = splmake(x, y, order=order)
            self._call = self.__class__._call_spline

        if len(x) < minval:
            raise ValueError("x and y arrays must have at "
                             "least %d entries" % minval)

        self._kind = kind
        self.x = x
        self._y = y

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
        y_lo = self._y[lo]
        y_hi = self._y[hi]

        # Note that the following two expressions rely on the specifics of the
        # broadcasting semantics.
        slope = (y_hi - y_lo) / (x_hi - x_lo)[:, None]

        # 5. Calculate the actual value for each entry in x_new.
        y_new = slope*(x_new - x_lo)[:, None] + y_lo

        return y_new

    def _call_nearest(self, x_new):
        """ Find nearest neighbour interpolated y_new = f(x_new)."""

        # 2. Find where in the averaged data the values to interpolate
        #    would be inserted.
        #    Note: use side='left' (right) to searchsorted() to define the
        #    halfway point to be nearest to the left (right) neighbour
        x_new_indices = searchsorted(self.x_bds, x_new, side='left')

        # 3. Clip x_new_indices so that they are within the range of x indices.
        x_new_indices = x_new_indices.clip(0, len(self.x)-1).astype(intp)

        # 4. Calculate the actual value for each entry in x_new.
        y_new = self._y[x_new_indices]

        return y_new

    def _call_spline(self, x_new):
        return spleval(self._spline, x_new)

    def _evaluate(self, x_new):
        # 1. Handle values in x_new that are outside of x.  Throw error,
        #    or return a list of mask array indicating the outofbounds values.
        #    The behavior is set by the bounds_error variable.
        x_new = asarray(x_new)
        out_of_bounds = self._check_bounds(x_new)
        y_new = self._call(self, x_new)
        if len(y_new) > 0:
            y_new[out_of_bounds] = self.fill_value
        return y_new

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


class _PPolyBase(object):
    """
    Base class for piecewise polynomials.
    """
    __slots__ = ('c', 'x', 'extrapolate')

    def __init__(self, c, x, extrapolate=None):
        self.c = np.asarray(c)
        self.x = np.ascontiguousarray(x, dtype=np.float64)
        if extrapolate is None:
            extrapolate = True
        self.extrapolate = bool(extrapolate)

        if self.x.ndim != 1:
            raise ValueError("x must be 1-dimensional")
        if self.x.size < 2:
            raise ValueError("at least 2 breakpoints are needed")
        if self.c.ndim < 2:
            raise ValueError("c must have at least 2 dimensions")
        if self.c.shape[0] == 0:
            raise ValueError("polynomial must be at least of order 0")
        if self.c.shape[1] != self.x.size-1:
            raise ValueError("number of coefficients != len(x)-1")
        if np.any(self.x[1:] - self.x[:-1] < 0):
            raise ValueError("x-coordinates are not in increasing order")

        dtype = self._get_dtype(self.c.dtype)
        self.c = np.ascontiguousarray(self.c, dtype=dtype)

    def _get_dtype(self, dtype):
        if np.issubdtype(dtype, np.complexfloating) \
               or np.issubdtype(self.c.dtype, np.complexfloating):
            return np.complex_
        else:
            return np.float_

    @classmethod
    def construct_fast(cls, c, x, extrapolate=None):
        """
        Construct the piecewise polynomial without making checks.

        Takes the same parameters as the constructor. Input arguments
        `c` and `x` must be arrays of the correct shape and type.  The
        `c` array can only be of dtypes float and complex, and `x`
        array must have dtype float.

        """
        self = object.__new__(cls)
        self.c = c
        self.x = x
        if extrapolate is None:
            extrapolate = True
        self.extrapolate = extrapolate
        return self

    def _ensure_c_contiguous(self):
        """
        c and x may be modified by the user. The Cython code expects
        that they are C contiguous.
        """
        if not self.x.flags.c_contiguous:
            self.x = self.x.copy()
        if not self.c.flags.c_contiguous:
            self.c = self.c.copy()

    def extend(self, c, x, right=True):
        """
        Add additional breakpoints and coefficients to the polynomial.

        Parameters
        ----------
        c : ndarray, size (k, m, ...)
            Additional coefficients for polynomials in intervals
            ``self.x[-1] <= x < x_right[0]``, ``x_right[0] <= x < x_right[1]``,
            ..., ``x_right[m-2] <= x < x_right[m-1]``
        x : ndarray, size (m,)
            Additional breakpoints. Must be sorted and either to
            the right or to the left of the current breakpoints.
        right : bool, optional
            Whether the new intervals are to the right or to the left
            of the current intervals.

        """
        c = np.asarray(c)
        x = np.asarray(x)

        if c.ndim < 2:
            raise ValueError("invalid dimensions for c")
        if x.ndim != 1:
            raise ValueError("invalid dimensions for x")
        if x.shape[0] != c.shape[1]:
            raise ValueError("x and c have incompatible sizes")
        if c.shape[2:] != self.c.shape[2:] or c.ndim != self.c.ndim:
            raise ValueError("c and self.c have incompatible shapes")
        if right:
            if x[0] < self.x[-1]:
                raise ValueError("new x are not to the right of current ones")
        else:
            if x[-1] > self.x[0]:
                raise ValueError("new x are not to the left of current ones")

        if c.size == 0:
            return

        dtype = self._get_dtype(c.dtype)

        k2 = max(c.shape[0], self.c.shape[0])
        c2 = np.zeros((k2, self.c.shape[1] + c.shape[1]) + self.c.shape[2:],
                      dtype=dtype)

        if right:
            c2[k2-self.c.shape[0]:, :self.c.shape[1]] = self.c
            c2[k2-c.shape[0]:, self.c.shape[1]:] = c
            self.x = np.r_[self.x, x]
        else:
            c2[k2-self.c.shape[0]:, :c.shape[1]] = c
            c2[k2-c.shape[0]:, c.shape[1]:] = self.c
            self.x = np.r_[x, self.x]

        self.c = c2

    def __call__(self, x, nu=0, extrapolate=None):
        """
        Evaluate the piecewise polynomial or its derivative

        Parameters
        ----------
        x : array-like
            Points to evaluate the interpolant at.
        nu : int, optional
            Order of derivative to evaluate. Must be non-negative.
        extrapolate : bool, optional
            Whether to extrapolate to ouf-of-bounds points based on first
            and last intervals, or to return NaNs.

        Returns
        -------
        y : array-like
            Interpolated values. Shape is determined by replacing
            the interpolation axis in the original array with the shape of x.

        Notes
        -----
        Derivatives are evaluated piecewise for each polynomial
        segment, even if the polynomial is not differentiable at the
        breakpoints. The polynomial intervals are considered half-open,
        ``[a, b)``, except for the last interval which is closed
        ``[a, b]``.

        """
        if extrapolate is None:
            extrapolate = self.extrapolate
        x = np.asarray(x)
        x_shape = x.shape
        x = np.ascontiguousarray(x.ravel(), dtype=np.float_)
        out = np.empty((len(x), prod(self.c.shape[2:])), dtype=self.c.dtype)
        self._ensure_c_contiguous()
        self._evaluate(x, nu, extrapolate, out)
        return out.reshape(x_shape + self.c.shape[2:])


class PPoly(_PPolyBase):
    """
    Piecewise polynomial in terms of coefficients and breakpoints

    The polynomial in the ith interval is ``x[i] <= xp < x[i+1]``::

        S = sum(c[m, i] * (xp - x[i])**(k-m) for m in range(k+1))

    where ``k`` is the degree of the polynomial. This representation
    is the local power basis.

    Parameters
    ----------
    c : ndarray, shape (k, m, ...)
        Polynomial coefficients, order `k` and `m` intervals
    x : ndarray, shape (m+1,)
        Polynomial breakpoints. These must be sorted in
        increasing order.
    extrapolate : bool, optional
        Whether to extrapolate to ouf-of-bounds points based on first
        and last intervals, or to return NaNs. Default: True.

    Attributes
    ----------
    x : ndarray
        Breakpoints.
    c : ndarray
        Coefficients of the polynomials. They are reshaped
        to a 3-dimensional array with the last dimension representing
        the trailing dimensions of the original coefficient array.

    Methods
    -------
    __call__
    derivative
    antiderivative
    integrate
    roots
    extend
    from_spline
    from_bernstein_basis
    construct_fast

    See also
    --------
    BPoly : piecewise polynomials in the Bernstein basis

    Notes
    -----
    High-order polynomials in the power basis can be numerically
    unstable.  Precision problems can start to appear for orders
    larger than 20-30.

    """

    def _evaluate(self, x, nu, extrapolate, out):
        _ppoly.evaluate(self.c.reshape(self.c.shape[0], self.c.shape[1], -1),
                        self.x, x, nu, bool(extrapolate), out)

    def derivative(self, nu=1):
        """
        Construct a new piecewise polynomial representing the derivative.

        Parameters
        ----------
        n : int, optional
            Order of derivative to evaluate. (Default: 1)
            If negative, the antiderivative is returned.

        Returns
        -------
        pp : PPoly
            Piecewise polynomial of order k2 = k - n representing the derivative
            of this polynomial.

        Notes
        -----
        Derivatives are evaluated piecewise for each polynomial
        segment, even if the polynomial is not differentiable at the
        breakpoints. The polynomial intervals are considered half-open,
        ``[a, b)``, except for the last interval which is closed
        ``[a, b]``.

        """
        if nu < 0:
            return self.antiderivative(-nu)

        # reduce order
        if nu == 0:
            c2 = self.c.copy()
        else:
            c2 = self.c[:-nu,:].copy()

        if c2.shape[0] == 0:
            # derivative of order 0 is zero
            c2 = np.zeros((1,) + c2.shape[1:], dtype=c2.dtype)

        # multiply by the correct rising factorials
        factor = spec.poch(np.arange(c2.shape[0], 0, -1), nu)
        c2 *= factor[(slice(None),) + (None,)*(c2.ndim-1)]

        # construct a compatible polynomial
        return self.construct_fast(c2, self.x, self.extrapolate)

    def antiderivative(self, nu=1):
        """
        Construct a new piecewise polynomial representing the antiderivative.

        Antiderivativative is also the indefinite integral of the function,
        and derivative is its inverse operation.

        Parameters
        ----------
        n : int, optional
            Order of antiderivative to evaluate. (Default: 1)
            If negative, the derivative is returned.

        Returns
        -------
        pp : PPoly
            Piecewise polynomial of order k2 = k + n representing
            the antiderivative of this polynomial.

        Notes
        -----
        The antiderivative returned by this function is continuous and
        continuously differentiable to order n-1, up to floating point
        rounding error.

        """
        if nu <= 0:
            return self.derivative(-nu)

        c = np.zeros((self.c.shape[0] + nu, self.c.shape[1]) + self.c.shape[2:],
                     dtype=self.c.dtype)
        c[:-nu] = self.c

        # divide by the correct rising factorials
        factor = spec.poch(np.arange(self.c.shape[0], 0, -1), nu)
        c[:-nu] /= factor[(slice(None),) + (None,)*(c.ndim-1)]

        # fix continuity of added degrees of freedom
        self._ensure_c_contiguous()
        _ppoly.fix_continuity(c.reshape(c.shape[0], c.shape[1], -1),
                              self.x, nu)

        # construct a compatible polynomial
        return self.construct_fast(c, self.x, self.extrapolate)

    def integrate(self, a, b, extrapolate=None):
        """
        Compute a definite integral over a piecewise polynomial.

        Parameters
        ----------
        a : float
            Lower integration bound
        b : float
            Upper integration bound
        extrapolate : bool, optional
            Whether to extrapolate to ouf-of-bounds points based on first
            and last intervals, or to return NaNs.

        Returns
        -------
        ig : array_like
            Definite integral of the piecewise polynomial over [a, b]

        """
        if extrapolate is None:
            extrapolate = self.extrapolate

        # Swap integration bounds if needed
        sign = 1
        if b < a:
            a, b = b, a
            sign = -1

        # Compute the integral
        range_int = np.empty((prod(self.c.shape[2:]),), dtype=self.c.dtype)
        self._ensure_c_contiguous()
        _ppoly.integrate(self.c.reshape(self.c.shape[0], self.c.shape[1], -1),
                         self.x, a, b, bool(extrapolate),
                         out=range_int)

        # Return
        range_int *= sign
        return range_int.reshape(self.c.shape[2:])

    def roots(self, discontinuity=True, extrapolate=None):
        """
        Find real roots of the piecewise polynomial.

        Parameters
        ----------
        discontinuity : bool, optional
            Whether to report sign changes across discontinuities at
            breakpoints as roots.
        extrapolate : bool, optional
            Whether to return roots from the polynomial extrapolated
            based on first and last intervals.

        Returns
        -------
        roots : ndarray
            Roots of the polynomial(s).

            If the PPoly object describes multiple polynomials, the
            return value is an object array whose each element is an
            ndarray containing the roots.

        Notes
        -----
        This routine works only on real-valued polynomials.

        If the piecewise polynomial contains sections that are
        identically zero, the root list will contain the start point
        of the corresponding interval, followed by a ``nan`` value.

        If the polynomial is discontinuous across a breakpoint, and
        there is a sign change across the breakpoint, this is reported
        if the `discont` parameter is True.

        Examples
        --------

        Finding roots of ``[x**2 - 1, (x - 1)**2]`` defined on intervals
        ``[-2, 1], [1, 2]``:

        >>> from scipy.interpolate import PPoly
        >>> pp = PPoly(np.array([[1, 0, -1], [1, 0, 0]]).T, [-2, 1, 2])
        >>> pp.roots()
        array([-1.,  1.])

        """
        if extrapolate is None:
            extrapolate = self.extrapolate

        self._ensure_c_contiguous()

        if np.issubdtype(self.c.dtype, np.complexfloating):
            raise ValueError("Root finding is only for "
                             "real-valued polynomials")

        r = _ppoly.real_roots(self.c.reshape(self.c.shape[0], self.c.shape[1], -1),
                              self.x, bool(discontinuity),
                              bool(extrapolate))
        if self.c.ndim == 2:
            return r[0]
        else:
            r2 = np.empty(prod(self.c.shape[2:]), dtype=object)
            # this for-loop is equivalent to ``r2[...] = r``, but that's broken
            # in numpy 1.6.0
            for ii, root in enumerate(r):
                r2[ii] = root

            return r2.reshape(self.c.shape[2:])

    @classmethod
    def from_spline(cls, tck, extrapolate=None):
        """
        Construct a piecewise polynomial from a spline

        Parameters
        ----------
        tck
            A spline, as returned by `splrep`
        extrapolate : bool, optional
            Whether to extrapolate to ouf-of-bounds points based on first
            and last intervals, or to return NaNs. Default: True.

        """
        t, c, k = tck

        cvals = np.empty((k + 1, len(t)-1), dtype=c.dtype)
        for m in xrange(k, -1, -1):
            y = fitpack.splev(t[:-1], tck, der=m)
            cvals[k - m, :] = y/spec.gamma(m+1)

        return cls.construct_fast(cvals, t, extrapolate)

    @classmethod
    def from_bernstein_basis(cls, bp, extrapolate=None):
        """
        Construct a piecewise polynomial in the power basis
        from a polynomial in Bernstein basis.

        Parameters
        ----------
        bp : BPoly
            A Bernstein basis polynomial, as created by BPoly
        extrapolate : bool, optional
            Whether to extrapolate to ouf-of-bounds points based on first
            and last intervals, or to return NaNs. Default: True.

        """
        dx = np.diff(bp.x)
        k = bp.c.shape[0] - 1  # polynomial order

        rest = (None,)*(bp.c.ndim-2)

        c = np.zeros_like(bp.c)
        for a in range(k+1):
            factor = (-1)**(a) * comb(k, a) * bp.c[a]
            for s in range(a, k+1):
                val = comb(k-a, s-a) * (-1)**s
                c[k-s] += factor * val / dx[(slice(None),)+rest]**s

        if extrapolate is None:
            extrapolate = bp.extrapolate

        return cls.construct_fast(c, bp.x, extrapolate)


class BPoly(_PPolyBase):
    """
    Piecewise polynomial in terms of coefficients and breakpoints

    The polynomial in the ``i``-th interval ``x[i] <= xp < x[i+1]``
    is written in the Bernstein polynomial basis::

        S = sum(c[a, i] * b(a, k; x) for a in range(k+1))

    where ``k`` is the degree of the polynomial, and::

        b(a, k; x) = comb(k, a) * t**k * (1 - t)**(k - a)

    with ``t = (x - x[i]) / (x[i+1] - x[i])``.

    Parameters
    ----------
    c : ndarray, shape (k, m, ...)
        Polynomial coefficients, order `k` and `m` intervals
    x : ndarray, shape (m+1,)
        Polynomial breakpoints. These must be sorted in
        increasing order.
    extrapolate : bool, optional
        Whether to extrapolate to ouf-of-bounds points based on first
        and last intervals, or to return NaNs. Default: True.

    Attributes
    ----------
    x : ndarray
        Breakpoints.
    c : ndarray
        Coefficients of the polynomials. They are reshaped
        to a 3-dimensional array with the last dimension representing
        the trailing dimensions of the original coefficient array.

    Methods
    -------
    __call__
    extend
    derivative
    construct_fast
    from_power_basis
    from_derivatives

    See also
    --------
    PPoly : piecewise polynomials in the power basis

    Notes
    -----
    Properties of Bernstein polynomials are well documented in the literature.
    Here's a non-exhaustive list:

    .. [1] http://en.wikipedia.org/wiki/Bernstein_polynomial

    .. [2] Kenneth I. Joy, Bernstein polynomials,
      http://www.idav.ucdavis.edu/education/CAGDNotes/Bernstein-Polynomials.pdf

    .. [3] E. H. Doha, A. H. Bhrawy, and M. A. Saker, Boundary Value Problems,
         vol 2011, article ID 829546, doi:10.1155/2011/829543

    Examples
    --------

    >>> x = [0, 1]
    >>> c = [[1], [2], [3]]
    >>> bp = BPoly(c, x)

    This creates a 2nd order polynomial

    .. math::

        B(x) = 1 \\times b_{0, 2}(x) + 2 \\times b_{1, 2}(x) + 3 \\times b_{2, 2}(x) \\\\
             = 1 \\times (1-x)^2 + 2 \\times 2 x (1 - x) + 3 \\times x^2

    """

    def _evaluate(self, x, nu, extrapolate, out):
        _ppoly.evaluate_bernstein(
            self.c.reshape(self.c.shape[0], self.c.shape[1], -1),
            self.x, x, nu, bool(extrapolate), out, self.c.dtype)

    def derivative(self, nu=1):
        """
        Construct a new piecewise polynomial representing the derivative.

        Parameters
        ----------
        nu : int, optional
            Order of derivative to evaluate. (Default: 1)
            If negative, the antiderivative is returned.

        Returns
        -------
        bp : BPoly
            Piecewise polynomial of order k2 = k - nu representing the derivative
            of this polynomial.

        """
        if nu < 0:
            raise NotImplementedError('Antiderivative not implemented.')

        if nu > 1:
            bp = self
            for k in range(nu):
                bp = bp.derivative()
            return bp

        # reduce order
        if nu == 0:
            c2 = self.c.copy()
        else:
            # For a polynomial
            #    B(x) = \sum_{a=0}^{k} c_a b_{a, k}(x),
            # we use the fact that
            #   b'_{a, k} = k ( b_{a-1, k-1} - b_{a, k-1} ),
            # which leads to
            #   B'(x) = \sum_{a=0}^{k-1} (c_{a+1} - c_a) b_{a, k-1}
            #
            # finally, for an interval [y, y + dy] with dy != 1,
            # we need to correct for an extra power of dy

            rest = (None,)*(self.c.ndim-2)

            k = self.c.shape[0] - 1
            dx = np.diff(self.x)[(None, slice(None))+rest]
            c2 = k * np.diff(self.c, axis=0) / dx

        if c2.shape[0] == 0:
            # derivative of order 0 is zero
            c2 = np.zeros((1,) + c2.shape[1:], dtype=c2.dtype)

        # construct a compatible polynomial
        return self.construct_fast(c2, self.x, self.extrapolate)

    def extend(self, c, x, right=True):
        k = max(self.c.shape[0], c.shape[0])
        self.c = self._raise_degree(self.c, k - self.c.shape[0])
        c = self._raise_degree(c, k - c.shape[0])
        return _PPolyBase.extend(self, c, x, right)
    extend.__doc__ = _PPolyBase.extend.__doc__

    @classmethod
    def from_power_basis(cls, pp, extrapolate=None):
        """
        Construct a piecewise polynomial in Bernstein basis
        from a power basis polynomial.

        Parameters
        ----------
        pp : PPoly
            A piecewise polynomial in the power basis
        extrapolate : bool, optional
            Whether to extrapolate to ouf-of-bounds points based on first
            and last intervals, or to return NaNs. Default: True.

        """
        dx = np.diff(pp.x)
        k = pp.c.shape[0] - 1   # polynomial order

        rest = (None,)*(pp.c.ndim-2)

        c = np.zeros_like(pp.c)
        for a in range(k+1):
            factor = pp.c[a] / comb(k, k-a) * dx[(slice(None),)+rest]**(k-a)
            for j in range(k-a, k+1):
                c[j] += factor * comb(j, k-a)

        if extrapolate is None:
            extrapolate = pp.extrapolate

        return cls.construct_fast(c, pp.x, extrapolate)

    @classmethod
    def from_derivatives(cls, xi, yi, orders=None, extrapolate=None):
        """Construct a piecewise polynomial in the Bernstein basis,
        compatible with the specified values and derivatives at breakpoints.

        Parameters
        ----------
        xi : array_like
            sorted 1D array of x-coordinates
        yi : array_like or list of array-likes
            ``yi[i][j]`` is the ``j``-th derivative known at ``xi[i]``
        orders : None or int or array_like of ints. Default: None.
            Specifies the degree of local polynomials. If not None, some
            derivatives are ignored.
        extrapolate : bool, optional
            Whether to extrapolate to ouf-of-bounds points based on first
            and last intervals, or to return NaNs. Default: True.

        Notes
        -----
        If ``k`` derivatives are specified at a breakpoint ``x``, the
        constructed polynomial is exactly ``k`` times continuously
        differentiable at ``x``, unless the ``order`` is provided explicitly.
        In the latter case, the smoothness of the polynomial at
        the breakpoint is controlled by the ``order``.

        Deduces the number of derivatives to match at each end
        from ``order`` and the number of derivatives available. If
        possible it uses the same number of derivatives from
        each end; if the number is odd it tries to take the
        extra one from y2. In any case if not enough derivatives
        are available at one end or another it draws enough to
        make up the total from the other end.

        If the order is too high and not enough derivatives are available,
        an exception is raised.

        Examples
        --------

        >>> BPoly.from_derivatives([0, 1], [[1, 2], [3, 4]])

        Creates a polynomial `f(x)` of degree 3, defined on `[0, 1]`
        such that `f(0) = 1, df/dx(0) = 2, f(1) = 3, df/dx(1) = 4`

        >>> BPoly.from_derivatives([0, 1, 2], [[0, 1], [0], [2]])

        Creates a piecewise polynomial `f(x)`, such that
        `f(0) = f(1) = 0`, `f(2) = 2`, and `df/dx(0) = 1`.
        Based on the number of derivatives provided, the order of the
        local polynomials is 2 on `[0, 1]` and 1 on `[1, 2]`.
        Notice that no restriction is imposed on the derivatives at
        `x = 1` and `x = 2`.

        Indeed, the explicit form of the polynomial is::

            f(x) = | x * (1 - x),  0 <= x < 1
                   | 2 * (x - 1),  1 <= x <= 2

        So that f'(1-0) = -1 and f'(1+0) = 2

        """
        xi = np.asarray(xi)
        if len(xi) != len(yi):
            raise ValueError("xi and yi need to have the same length")
        if np.any(xi[1:] - xi[:1] <= 0):
            raise ValueError("x coordinates are not in increasing order")

        # number of intervals
        m = len(xi) - 1

        # global poly order is k-1, local orders are <=k and can vary
        try:
            k = max(len(yi[i]) + len(yi[i+1]) for i in range(m))
        except TypeError:
            raise ValueError("Using a 1D array for y? Please .reshape(-1, 1).")

        if orders is None:
            orders = [None] * m
        else:
            if isinstance(orders, integer_types):
                orders = [orders] * m
            k = max(k, max(orders))

            if any(o <= 0 for o in orders):
                raise ValueError("Orders must be positive.")

        c = []
        for i in range(m):
            y1, y2 = yi[i], yi[i+1]
            if orders[i] is None:
                n1, n2 = len(y1), len(y2)
            else:
                n = orders[i]+1
                n1 = min(n//2, len(y1))
                n2 = min(n - n1, len(y2))
                n1 = min(n - n2, len(y2))
                if n1+n2 != n:
                    raise ValueError("Point %g has %d derivatives, point %g"
                            " has %d derivatives, but order %d requested" %
                            (xi[i], len(y1), xi[i+1], len(y2), orders[i]))
                if not (n1 <= len(y1) and n2 <= len(y2)):
                    raise ValueError("`order` input incompatible with"
                            " length y1 or y2.")

            b = BPoly._construct_from_derivatives(xi[i], xi[i+1], y1[:n1], y2[:n2])
            if len(b) < k:
                b = BPoly._raise_degree(b, k - len(b))
            c.append(b)

        c = np.asarray(c)
        return cls(c.swapaxes(0, 1), xi, extrapolate)

    @staticmethod
    def _construct_from_derivatives(xa, xb, ya, yb):
        """Compute the coefficients of a polynomial in the Bernstein basis
        given the values and derivatives at the edges.

        Return the coefficients of a polynomial in the Bernstein basis
        defined on `[xa, xb]` and having the values and derivatives at the
        endpoints ``xa`` and ``xb`` as specified by ``ya`` and ``yb``.
        The polynomial constructed is of the minimal possible degree, i.e.,
        if the lengths of ``ya`` and ``yb`` are ``na`` and ``nb``, the degree
        of the polynomial is ``na + nb - 1``.

        Parameters
        ----------
        xa : float
            Left-hand end point of the interval
        xb : float
            Right-hand end point of the interval
        ya : array_like
            Derivatives at ``xa``. ``ya[0]`` is the value of the function, and
            ``ya[i]`` for ``i > 0`` is the value of the ``i``-th derivative.
        yb : array_like
            Derivatives at ``xb``.

        Returns
        -------
        array
            coefficient array of a polynomial having specified derivatives

        Notes
        -----
        This uses several facts from life of Bernstein basis functions.
        First of all,

            .. math:: b'_{a, n} = n (b_{a-1, n-1} - b_{a, n-1})

        If B(x) is a linear combination of the form

            .. math:: B(x) = \sum_{a=0}^{n} c_a b_{a, n},

        then :math: B'(x) = n \sum_{a=0}^{n-1} (c_{a+1} - c_{a}) b_{a, n-1}.
        Iterating the latter one, one finds for the q-th derivative

            .. math:: B^{q}(x) = n!/(n-q)! \sum_{a=0}^{n-q} Q_a b_{a, n-q},

        with

          .. math:: Q_a = \sum_{j=0}^{q} (-)^{j+q} comb(q, j) c_{j+a}

        This way, only `a=0` contributes to :math: `B^{q}(x = xa)`, and
        `c_q` are found one by one by iterating `q = 0, ..., na`.

        At `x = xb` it's the same with `a = n - q`.

        """
        ya, yb = np.asarray(ya), np.asarray(yb)
        if ya.shape[1:] != yb.shape[1:]:
            raise ValueError('ya and yb have incompatible dimensions.')

        dta, dtb = ya.dtype, yb.dtype
        if (np.issubdtype(dta, np.complexfloating)
               or np.issubdtype(dtb, np.complexfloating)):
            dt = np.complex_
        else:
            dt = np.float_

        na, nb = len(ya), len(yb)
        n = na + nb

        c = np.empty((na+nb,) + ya.shape[1:], dtype=dt)

        # compute coefficients of a polynomial degree na+nb-1
        # walk left-to-right
        for q in range(0, na):
            c[q] = ya[q] / spec.poch(n - q, q) * (xb - xa)**q
            for j in range(0, q):
                c[q] -= (-1)**(j+q) * comb(q, j) * c[j]

        # now walk right-to-left
        for q in range(0, nb):
            c[-q-1] = yb[q] / spec.poch(n - q, q) * (-1)**q * (xb - xa)**q
            for j in range(0, q):
                c[-q-1] -= (-1)**(j+1) * comb(q, j+1) * c[-q+j]

        return c

    @staticmethod
    def _raise_degree(c, d):
        """Raise a degree of a polynomial in the Bernstein basis.

        Given the coefficients of a polynomial degree `k`, return (the
        coefficients of) the equivalent polynomial of degree `k+d`.

        Parameters
        ----------
        c : array_like
            coefficient array, 1D
        d : integer

        Returns
        -------
        array
            coefficient array, 1D array of length `c.shape[0] + d`

        Notes
        -----
        This uses the fact that a Bernstein polynomial `b_{a, k}` can be
        identically represented as a linear combination of polynomials of
        a higher degree `k+d`:

            .. math:: b_{a, k} = comb(k, a) \sum_{j=0}^{d} b_{a+j, k+d} \
                                comb(d, j) / comb(k+d, a+j)

        """
        if d == 0:
            return c

        k = c.shape[0] - 1
        out = np.zeros((c.shape[0] + d,) + c.shape[1:], dtype=c.dtype)

        for a in range(c.shape[0]):
            f = c[a] * comb(k, a)
            for j in range(d+1):
                out[a+j] += f * comb(d, j) / comb(k+d, a+j)
        return out


class RegularGridInterpolator(object):
    """
    Interpolation on a regular grid in arbitrary dimensions

    The data must be defined on a regular grid; the grid spacing however may be
    uneven.  Linear and nearest-neighbour interpolation are supported. After
    setting up the interpolator object, the interpolation method (*linear* or
    *nearest*) may be chosen at each evaluation.

    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, )
        The points defining the regular grid in n dimensions.

    values : array_like, shape (m1, ..., mn, ...)
        The data on the regular grid in n dimensions.

    method : str
        The method of interpolation to perform. Supported are "linear" and
        "nearest". This parameter will become the default for the object's
        ``__call__`` method.

    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then `fill_value` is used.

    fill_value : number, optional
        If provided, the value to use for points outside of the
        interpolation domain. If None, values outside
        the domain are extrapolated.

    Methods
    -------
    __call__

    Notes
    -----
    Contrary to LinearNDInterpolator and NearestNDInterpolator, this class
    avoids expensive triangulation of the input data by taking advantage of the
    regular grid structure.

    .. versionadded:: 0.14

    See also
    --------
    NearestNDInterpolator : Nearest neighbour interpolation on unstructured
                            data in N dimensions

    LinearNDInterpolator : Piecewise linear interpolant on unstructured data
                           in N dimensions

    References
    ----------
    .. [1] Python package *regulargrid* by Johannes Buchner, see
           https://pypi.python.org/pypi/regulargrid/
    .. [2] Trilinear interpolation. (2013, January 17). In Wikipedia, The Free
           Encyclopedia. Retrieved 27 Feb 2013 01:28.
           http://en.wikipedia.org/w/index.php?title=Trilinear_interpolation&oldid=533448871
    .. [3] Weiser, Alan, and Sergio E. Zarantonello. "A note on piecewise linear
           and multilinear table interpolation in many dimensions." MATH.
           COMPUT. 50.181 (1988): 189-196.
           http://www.ams.org/journals/mcom/1988-50-181/S0025-5718-1988-0917826-0/S0025-5718-1988-0917826-0.pdf

    """
    # this class is based on code originally programmed by Johannes Buchner,
    # see https://github.com/JohannesBuchner/regulargrid

    def __init__(self, points, values, method="linear", bounds_error=True,
                 fill_value=np.nan):
        if method not in ["linear", "nearest"]:
            raise ValueError("Method '%s' is not defined" % method)
        self.method = method
        self.bounds_error = bounds_error

        if not hasattr(values, 'ndim'):
            # allow reasonable duck-typed values
            values = np.asarray(values)

        if len(points) > values.ndim:
            raise ValueError("There are %d point arrays, but values has %d "
                             "dimensions" % (len(points), values.ndim))

        if hasattr(values, 'dtype') and hasattr(values, 'astype'):
            if not np.issubdtype(values.dtype, np.inexact):
                values = values.astype(float)

        self.fill_value = fill_value
        if fill_value is not None:
            fill_value_dtype = np.asarray(fill_value).dtype
            if (hasattr(values, 'dtype')
                    and not np.can_cast(fill_value_dtype, values.dtype,
                                        casting='same_kind')):
                raise ValueError("fill_value must be either 'None' or "
                                 "of a type compatible with values")

        for i, p in enumerate(points):
            if not np.all(np.diff(p) > 0.):
                raise ValueError("The points in dimension %d must be strictly "
                                 "ascending" % i)
            if not np.asarray(p).ndim == 1:
                raise ValueError("The points in dimension %d must be "
                                 "1-dimensional" % i)
            if not values.shape[i] == len(p):
                raise ValueError("There are %d points and %d values in "
                                 "dimension %d" % (len(p), values.shape[i], i))
        self.grid = tuple([np.asarray(p) for p in points])
        self.values = values

    def __call__(self, xi, method=None):
        """
        Interpolation at coordinates

        Parameters
        ----------
        xi : ndarray of shape (..., ndim)
            The coordinates to sample the gridded data at

        method : str
            The method of interpolation to perform. Supported are "linear" and
            "nearest".

        """
        method = self.method if method is None else method
        if method not in ["linear", "nearest"]:
            raise ValueError("Method '%s' is not defined" % method)

        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError("The requested sample points xi have dimension "
                             "%d, but this RegularGridInterpolator has "
                             "dimension %d" % (xi.shape[1], ndim))

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        if self.bounds_error:
            for i, p in enumerate(xi.T):
                if not np.logical_and(np.all(self.grid[i][0] <= p),
                                      np.all(p <= self.grid[i][-1])):
                    raise ValueError("One of the requested xi is out of bounds "
                                     "in dimension %d" % i)

        indices, norm_distances, out_of_bounds = self._find_indices(xi.T)
        if method == "linear":
            result = self._evaluate_linear(indices, norm_distances, out_of_bounds)
        elif method == "nearest":
            result = self._evaluate_nearest(indices, norm_distances, out_of_bounds)
        if not self.bounds_error and self.fill_value is not None:
            result[out_of_bounds] = self.fill_value

        return result.reshape(xi_shape[:-1] + self.values.shape[ndim:])

    def _evaluate_linear(self, indices, norm_distances, out_of_bounds):
        # slice for broadcasting over trailing dimensions in self.values
        vslice = (slice(None),) + (None,)*(self.values.ndim - len(indices))

        # find relevant values
        # each i and i+1 represents a edge
        edges = itertools.product(*[[i, i + 1] for i in indices])
        values = 0.
        for edge_indices in edges:
            weight = 1.
            for ei, i, yi in zip(edge_indices, indices, norm_distances):
                weight *= np.where(ei == i, 1 - yi, yi)
            values += np.asarray(self.values[edge_indices]) * weight[vslice]
        return values

    def _evaluate_nearest(self, indices, norm_distances, out_of_bounds):
        idx_res = []
        for i, yi in zip(indices, norm_distances):
            idx_res.append(np.where(yi <= .5, i, i + 1))
        return self.values[idx_res]

    def _find_indices(self, xi):
        # find relevant edges between which xi are situated
        indices = []
        # compute distance to lower edge in unity units
        norm_distances = []
        # check for out of bounds xi
        out_of_bounds = np.zeros((xi.shape[1]), dtype=bool)
        # iterate through dimensions
        for x, grid in zip(xi, self.grid):
            i = np.searchsorted(grid, x) - 1
            i[i < 0] = 0
            i[i > grid.size - 2] = grid.size - 2
            indices.append(i)
            norm_distances.append((x - grid[i]) /
                                  (grid[i + 1] - grid[i]))
            if not self.bounds_error:
                out_of_bounds += x < grid[0]
                out_of_bounds += x > grid[-1]
        return indices, norm_distances, out_of_bounds


def interpn(points, values, xi, method="linear", bounds_error=True,
            fill_value=np.nan):
    """
    Multidimensional interpolation on regular grids.

    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, )
        The points defining the regular grid in n dimensions.

    values : array_like, shape (m1, ..., mn, ...)
        The data on the regular grid in n dimensions.

    xi : ndarray of shape (..., ndim)
        The coordinates to sample the gridded data at

    method : str
        The method of interpolation to perform. Supported are "linear" and
        "nearest", and "splinef2d". "splinef2d" is only supported for
        2-dimensional data.

    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then `fill_value` is used.

    fill_value : number, optional
        If provided, the value to use for points outside of the
        interpolation domain. If None, values outside
        the domain are extrapolated.  Extrapolation is not supported by method
        "splinef2d".

    Returns
    -------
    values_x : ndarray, shape xi.shape[:-1] + values.shape[ndim:]
        Interpolated values at input coordinates.

    Notes
    -----

    .. versionadded:: 0.14

    See also
    --------
    NearestNDInterpolator : Nearest neighbour interpolation on unstructured
                            data in N dimensions

    LinearNDInterpolator : Piecewise linear interpolant on unstructured data
                           in N dimensions

    RegularGridInterpolator : Linear and nearest-neighbor Interpolation on a
                              regular grid in arbitrary dimensions

    RectBivariateSpline : Bivariate spline approximation over a rectangular mesh

    """
    # sanity check 'method' kwarg
    if method not in ["linear", "nearest", "splinef2d"]:
        raise ValueError("interpn only understands the methods 'linear', "
                         "'nearest', and 'splinef2d'. You provided %s." %
                         method)

    if not hasattr(values, 'ndim'):
        values = np.asarray(values)

    ndim = values.ndim
    if ndim > 2 and method == "splinef2d":
        raise ValueError("The method spline2fd can only be used for "
                         "2-dimensional input data")
    if not bounds_error and fill_value is None and method == "splinef2d":
        raise ValueError("The method spline2fd does not support extrapolation.")

    # sanity check consistency of input dimensions
    if len(points) > ndim:
        raise ValueError("There are %d point arrays, but values has %d "
                         "dimensions" % (len(points), ndim))
    if len(points) != ndim and method == 'splinef2d':
        raise ValueError("The method spline2fd can only be used for "
                         "scalar data with one point per coordinate")

    # sanity check input grid
    for i, p in enumerate(points):
        if not np.all(np.diff(p) > 0.):
            raise ValueError("The points in dimension %d must be strictly "
                             "ascending" % i)
        if not np.asarray(p).ndim == 1:
            raise ValueError("The points in dimension %d must be "
                             "1-dimensional" % i)
        if not values.shape[i] == len(p):
            raise ValueError("There are %d points and %d values in "
                             "dimension %d" % (len(p), values.shape[i], i))
    grid = tuple([np.asarray(p) for p in points])

    # sanity check requested xi
    xi = _ndim_coords_from_arrays(xi, ndim=len(grid))
    if xi.shape[-1] != len(grid):
        raise ValueError("The requested sample points xi have dimension "
                         "%d, but this RegularGridInterpolator has "
                         "dimension %d" % (xi.shape[1], len(grid)))

    for i, p in enumerate(xi.T):
        if bounds_error and not np.logical_and(np.all(grid[i][0] <= p),
                                               np.all(p <= grid[i][-1])):
            raise ValueError("One of the requested xi is out of bounds "
                             "in dimension %d" % i)

    # perform interpolation
    if method == "linear":
        interp = RegularGridInterpolator(points, values, method="linear",
                                         bounds_error=bounds_error,
                                         fill_value=fill_value)
        return interp(xi)
    elif method == "nearest":
        interp = RegularGridInterpolator(points, values, method="nearest",
                                         bounds_error=bounds_error,
                                         fill_value=fill_value)
        return interp(xi)
    elif method == "splinef2d":
        xi_shape = xi.shape
        xi = xi.reshape(-1, xi.shape[-1])

        # RectBivariateSpline doesn't support fill_value; we need to wrap here
        idx_valid = np.all((grid[0][0] <= xi[:, 0], xi[:, 0] <= grid[0][-1],
                            grid[1][0] <= xi[:, 1], xi[:, 1] <= grid[1][-1]),
                           axis=0)
        result = np.empty_like(xi[:, 0])

        # make a copy of values for RectBivariateSpline
        interp = RectBivariateSpline(points[0], points[1], values[:])
        result[idx_valid] = interp.ev(xi[idx_valid, 0], xi[idx_valid, 1])
        result[np.logical_not(idx_valid)] = fill_value

        return result.reshape(xi_shape[:-1])


# backward compatibility wrapper
class ppform(PPoly):
    """
    Deprecated piecewise polynomial class.

    New code should use the `PPoly` class instead.

    """

    def __init__(self, coeffs, breaks, fill=0.0, sort=False):
        warnings.warn("ppform is deprecated -- use PPoly instead",
                      category=DeprecationWarning)

        if sort:
            breaks = np.sort(breaks)
        else:
            breaks = np.asarray(breaks)

        PPoly.__init__(self, coeffs, breaks)

        self.coeffs = self.c
        self.breaks = self.x
        self.K = self.coeffs.shape[0]
        self.fill = fill
        self.a = self.breaks[0]
        self.b = self.breaks[-1]

    def __call__(self, x):
        return PPoly.__call__(self, x, 0, False)

    def _evaluate(self, x, nu, extrapolate, out):
        PPoly._evaluate(self, x, nu, extrapolate, out)
        out[~((x >= self.a) & (x <= self.b))] = self.fill
        return out

    @classmethod
    def fromspline(cls, xk, cvals, order, fill=0.0):
        # Note: this spline representation is incompatible with FITPACK
        N = len(xk)-1
        sivals = np.empty((order+1, N), dtype=float)
        for m in xrange(order, -1, -1):
            fact = spec.gamma(m+1)
            res = _fitpack._bspleval(xk[:-1], xk, cvals, order, m)
            res /= fact
            sivals[order-m, :] = res
        return cls(sivals, xk, fill=fill)


def _dot0(a, b):
    """Similar to numpy.dot, but sum over last axis of a and 1st axis of b"""
    if b.ndim <= 2:
        return dot(a, b)
    else:
        axes = list(range(b.ndim))
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
    u, s, vh = scipy.linalg.svd(B)
    ind = K-1
    V2 = vh[-ind:,:].T
    V1 = vh[:-ind,:].T
    A = dot(J.T,J)
    tmp = dot(V2.T,A)
    Q = dot(tmp,V2)
    p = scipy.linalg.solve(Q, tmp)
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
        num = N - k
    else:
        num = M + k
        start = abs(k)*N
    end = start + num*(N+1)-1
    a.flat[start:end:(N+1)] = v

# Return the spline that minimizes the dis-continuity of the
# "order-th" derivative; for order >= 2.


def _find_smoothest2(xk, yk):
    N = len(xk) - 1
    Np1 = N + 1
    # find pseudo-inverse of B directly.
    Bd = np.empty((Np1, N))
    for k in range(-N,N):
        if (k < 0):
            l = np.arange(-k, Np1)
            v = (l+k+1)
            if ((k+1) % 2):
                v = -v
        else:
            l = np.arange(k,N)
            v = N - l
            if ((k % 2)):
                v = -v
        _setdiag(Bd, k, v)
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
    B = np.concatenate((B, lh), axis=0)
    w = np.concatenate((yk, rh), axis=0)
    M, N = B.shape
    if (M > N):
        raise ValueError("over-specification of conditions")
    elif (M < N):
        return _find_smoothest(xk, yk, order, None, B)
    else:
        return scipy.linalg.solve(B, w)

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


def splmake(xk, yk, order=3, kind='smoothest', conds=None):
    """
    Return a representation of a spline given data-points at internal knots

    Parameters
    ----------
    xk : array_like
        The input array of x values of rank 1
    yk : array_like
        The input array of y values of rank N. `yk` can be an N-d array to
        represent more than one curve, through the same `xk` points. The first
        dimension is assumed to be the interpolating dimension and is the same
        length of `xk`.
    order : int, optional
        Order of the spline
    kind : str, optional
        Can be 'smoothest', 'not_a_knot', 'fixed', 'clamped', 'natural',
        'periodic', 'symmetric', 'user', 'mixed' and it is ignored if order < 2
    conds : optional
        Conds

    Returns
    -------
    splmake : tuple
        Return a (`xk`, `cvals`, `k`) representation of a spline given
        data-points where the (internal) knots are at the data-points.

    """
    yk = np.asanyarray(yk)

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


def spleval(xck, xnew, deriv=0):
    """
    Evaluate a fixed spline represented by the given tuple at the new x-values

    The `xj` values are the interior knot points.  The approximation
    region is `xj[0]` to `xj[-1]`.  If N+1 is the length of `xj`, then `cvals`
    should have length N+k where `k` is the order of the spline.

    Parameters
    ----------
    (xj, cvals, k) : tuple
        Parameters that define the fixed spline
    xj : array_like
        Interior knot points
    cvals : array_like
        Curvature
    k : int
        Order of the spline
    xnew : array_like
        Locations to calculate spline
    deriv : int
        Deriv

    Returns
    -------
    spleval : ndarray
        If `cvals` represents more than one curve (`cvals.ndim` > 1) and/or
        `xnew` is N-d, then the result is `xnew.shape` + `cvals.shape[1:]`
        providing the interpolation of multiple curves.

    Notes
    -----
    Internally, an additional `k`-1 knot points are added on either side of
    the spline.

    """
    (xj,cvals,k) = xck
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


def spltopp(xk, cvals, k):
    """Return a piece-wise polynomial object from a fixed-spline tuple.
    """
    return ppform.fromspline(xk, cvals, k)


def spline(xk, yk, xnew, order=3, kind='smoothest', conds=None):
    """
    Interpolate a curve at new points using a spline fit

    Parameters
    ----------
    xk, yk : array_like
        The x and y values that define the curve.
    xnew : array_like
        The x values where spline should estimate the y values.
    order : int
        Default is 3.
    kind : string
        One of {'smoothest'}
    conds : Don't know
        Don't know

    Returns
    -------
    spline : ndarray
        An array of y values; the spline evaluated at the positions `xnew`.

    """
    return spleval(splmake(xk,yk,order=order,kind=kind,conds=conds),xnew)
