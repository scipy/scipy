from __future__ import division, print_function, absolute_import

import numpy as np

from . import BPoly, PPoly
from .polyint import _isscalar


__all__ = ["PchipInterpolator", "pchip_interpolate", "pchip",
           "Akima1DInterpolator"]


class PchipInterpolator(object):
    """PCHIP 1-d monotonic cubic interpolation

    x and y are arrays of values used to approximate some function f,
    with ``y = f(x)``.  The interpolant uses monotonic cubic splines
    to find the value of new points. (PCHIP stands for Piecewise Cubic
    Hermite Interpolating Polynomial).

    Parameters
    ----------
    x : ndarray
        A 1-D array of monotonically increasing real values.  `x` cannot
        include duplicate values (otherwise f is overspecified)
    y : ndarray
        A 1-D array of real values.  `y`'s length along the interpolation
        axis must be equal to the length of `x`. If N-D array, use axis 
        parameter to select correct axis.
    axis : int, optional
        Axis in the y array corresponding to the x-coordinate values.
    extrapolate : bool, optional
        Whether to extrapolate to ouf-of-bounds points based on first
        and last intervals, or to return NaNs.

    Methods
    -------
    __call__
    derivative

    See Also
    --------
    Akima1DInterpolator

    Notes
    -----
    The first derivatives are guaranteed to be continuous, but the second
    derivatives may jump at x_k. 

    Preserves monotonicity in the interpolation data and does not overshoot
    if the data is not smooth.

    Determines the derivatives at the points x_k, d_k, by using PCHIP algorithm:
      
    Let m_k be the slope of the kth segment (between k and k+1)
    If m_k=0 or m_{k-1}=0 or sgn(m_k) != sgn(m_{k-1}) then d_k == 0
    else use weighted harmonic mean:

       w_1 = 2h_k + h_{k-1}, w_2 = h_k + 2h_{k-1}
       1/d_k = 1/(w_1 + w_2)*(w_1 / m_k + w_2 / m_{k-1})

    where h_k is the spacing between x_k and x_{k+1}.

    """
    def __init__(self, x, y, axis=0, extrapolate=None):
        x = np.asarray(x)
        if not np.issubdtype(x.dtype, np.inexact):
            x = x.astype(float)

        y = np.asarray(y)
        if not np.issubdtype(y.dtype, np.inexact):
            y = y.astype(float)

        axis = axis % y.ndim
        
        xp = x.reshape((x.shape[0],) + (1,)*(y.ndim-1))
        yp = np.rollaxis(y, axis)

        dk = self._find_derivatives(xp, yp)
        data = np.hstack((yp[:, None, ...], dk[:, None, ...]))

        self._bpoly = BPoly.from_derivatives(x, data, orders=None,
                extrapolate=extrapolate)
        self.axis = axis

    def __call__(self, x, der=0, extrapolate=None):
        """
        Evaluate the PCHIP interpolant or its derivative.

        Parameters
        ----------
        x : array-like
            Points to evaluate the interpolant at.
        der : int, optional
            Order of derivative to evaluate. Must be non-negative.
        extrapolate : bool, optional
            Whether to extrapolate to ouf-of-bounds points based on first
            and last intervals, or to return NaNs.

        Returns
        -------
        y : array-like
            Interpolated values. Shape is determined by replacing
            the interpolation axis in the original array with the shape of x.

        """
        out = self._bpoly(x, der, extrapolate)
        return self._reshaper(x, out)

    def derivative(self, der=1):
        """
        Construct a piecewise polynomial representing the derivative.

        Parameters
        ----------
        der : int, optional
            Order of derivative to evaluate. (Default: 1)
            If negative, the antiderivative is returned.

        Returns
        ------- 
        Piecewise polynomial of order k2 = k - der representing the derivative
        of this polynomial.

        """
        t = object.__new__(self.__class__)
        t.axis = self.axis
        t._bpoly = self._bpoly.derivative(der)
        return t

    def roots(self):
        """
        Return the roots of the interpolated function.
        """
        return (PPoly.from_bernstein_basis(self._bpoly)).roots()

    def _reshaper(self, x, out):
        x = np.asarray(x)
        l = x.ndim
        transp = (tuple(range(l, l+self.axis)) + tuple(range(l)) +
                tuple(range(l+self.axis, out.ndim)))
        return out.transpose(transp)

    @staticmethod
    def _edge_case(m0, d1, out):
        m0 = np.atleast_1d(m0)
        d1 = np.atleast_1d(d1)
        mask = (d1 != 0) & (m0 != 0)
        out[mask] = 1.0/(1.0/m0[mask]+1.0/d1[mask])

    @staticmethod
    def _find_derivatives(x, y):
        # Determine the derivatives at the points y_k, d_k, by using
        #  PCHIP algorithm is:
        # We choose the derivatives at the point x_k by
        # Let m_k be the slope of the kth segment (between k and k+1)
        # If m_k=0 or m_{k-1}=0 or sgn(m_k) != sgn(m_{k-1}) then d_k == 0
        # else use weighted harmonic mean:
        #   w_1 = 2h_k + h_{k-1}, w_2 = h_k + 2h_{k-1}
        #   1/d_k = 1/(w_1 + w_2)*(w_1 / m_k + w_2 / m_{k-1})
        #   where h_k is the spacing between x_k and x_{k+1}
        y_shape = y.shape
        if y.ndim == 1:
            # So that _edge_case doesn't end up assigning to scalars
            x = x[:,None]
            y = y[:,None]

        hk = x[1:] - x[:-1]
        mk = (y[1:] - y[:-1]) / hk
        smk = np.sign(mk)
        condition = ((smk[1:] != smk[:-1]) | (mk[1:] == 0) | (mk[:-1] == 0))

        w1 = 2*hk[1:] + hk[:-1]
        w2 = hk[1:] + 2*hk[:-1]
        # values where division by zero occurs will be excluded
        # by 'condition' afterwards
        with np.errstate(divide='ignore'):
            whmean = 1.0/(w1+w2)*(w1/mk[1:] + w2/mk[:-1])
        dk = np.zeros_like(y)
        dk[1:-1][condition] = 0.0
        dk[1:-1][~condition] = 1.0/whmean[~condition]

        # For end-points choose d_0 so that 1/d_0 = 1/m_0 + 1/d_1 unless
        #  one of d_1 or m_0 is 0, then choose d_0 = 0
        PchipInterpolator._edge_case(mk[0],dk[1], dk[0])
        PchipInterpolator._edge_case(mk[-1],dk[-2], dk[-1])

        return dk.reshape(y_shape)


def pchip_interpolate(xi, yi, x, der=0, axis=0):
    """
    Convenience function for pchip interpolation. 
    xi and yi are arrays of values used to approximate some function f,
    with ``yi = f(xi)``.  The interpolant uses monotonic cubic splines
    to find the value of new points x and the derivatives there.

    See `PchipInterpolator` for details.

    Parameters
    ----------
    xi : array_like
        A sorted list of x-coordinates, of length N.
    yi :  array_like
        A 1-D array of real values.  `yi`'s length along the interpolation
        axis must be equal to the length of `xi`. If N-D array, use axis 
        parameter to select correct axis.
    x : scalar or array_like
        Of length M.
    der : integer or list
        How many derivatives to extract; None for all potentially
        nonzero derivatives (that is a number equal to the number
        of points), or a list of derivatives to extract. This number
        includes the function value as 0th derivative.
    axis : int, optional
        Axis in the yi array corresponding to the x-coordinate values.

    See Also
    --------
    PchipInterpolator

    Returns
    -------
    y : scalar or array_like
        The result, of length R or length M or M by R,

    """
    P = PchipInterpolator(xi, yi, axis=axis)
    if der == 0:
        return P(x)
    elif _isscalar(der):
        return P(x, der=der)
    else:
        return [P(x, nu) for nu in der]


# Backwards compatibility
pchip = PchipInterpolator


class Akima1DInterpolator(PPoly):
    """
    Akima interpolator

    Fit piecewise cubic polynomials, given vectors x and y. The interpolation
    method by Akima uses a continuously differentiable sub-spline built from
    piecewise cubic polynomials. The resultant curve passes through the given
    data points and will appear smooth and natural.

    Parameters
    ----------
    x : ndarray, shape (m, )
        1-D array of monotonically increasing real values.
    y : ndarray, shape (m, ...)
        N-D array of real values. The length of *y* along the first axis must
        be equal to the length of *x*.
    axis : int, optional
        Specifies the axis of *y* along which to interpolate. Interpolation
        defaults to the last axis of *y*.

    Methods
    -------
    __call__

    See Also
    --------
    PchipInterpolator

    Notes
    -----
    .. versionadded:: 0.14

    Use only for precise data, as the fitted curve passes through the given
    points exactly. This routine is useful for plotting a pleasingly smooth
    curve through a few given points for purposes of plotting.

    References
    ----------
    [1] A new method of interpolation and smooth curve fitting based
        on local procedures. Hiroshi Akima, J. ACM, October 1970, 17(4),
        589-602.

    """

    def __init__(self, x, y):
        # Original implementation in MATLAB by N. Shamsundar (BSD licensed), see
        # http://www.mathworks.de/matlabcentral/fileexchange/1814-akima-interpolation
        if np.any(np.diff(x) < 0.):
            raise ValueError("x must be strictly ascending")
        if x.ndim != 1:
            raise ValueError("x must be 1-dimensional")
        if x.size < 2:
            raise ValueError("at least 2 breakpoints are needed")
        if x.size != y.shape[0]:
            raise ValueError("x.shape must equal y.shape[0]")

        # determine slopes between breakpoints
        m = np.empty((x.size + 3, ) + y.shape[1:])
        dx = np.diff(x)
        dx = dx[(slice(None), ) + (None, ) * (y.ndim - 1)]
        m[2:-2] = np.diff(y, axis=0) / dx

        # add two additional points on the left ...
        m[1] = 2. * m[2] - m[3]
        m[0] = 2. * m[1] - m[2]
        # ... and on the right
        m[-2] = 2. * m[-3] - m[-4]
        m[-1] = 2. * m[-2] - m[-3]

        # if m1 == m2 != m3 == m4, the slope at the breakpoint is not defined.
        # This is the fill value:
        t = .5 * (m[3:] + m[:-3])
        # get the denominator of the slope t
        dm = np.abs(np.diff(m, axis=0))
        f1 = dm[2:]
        f2 = dm[:-2]
        f12 = f1 + f2
        # These are the indices where the the slope at breakpoint is defined:
        id_ = np.nonzero(f12 > 1e-9 * np.max(f12))[0]
        # set the slope at breakpoint
        t[id_] = (f1[id_] * m[id_ + 1] + f2[id_] * m[id_ + 2]) / f12[id_]

        # calculate the higher order coefficients
        c = (3. * m[2:-2] - 2. * t[:-1] - t[1:]) / dx
        d = (t[:-1] + t[1:] - 2. * m[2:-2]) / dx ** 2

        coeff = np.zeros((4, x.size - 1) + y.shape[1:])
        coeff[3] = y[:-1]
        coeff[2] = t[:-1]
        coeff[1] = c
        coeff[0] = d

        super(Akima1DInterpolator, self).__init__(coeff, x, extrapolate=False)

    def extend(self):
        raise NotImplementedError("Extending a 1D Akima interpolator is not "
                "yet implemented")

    # These are inherited from PPoly, but they do not produce an Akima
    # interpolor. Hence stub them out.
    @classmethod    
    def from_spline(cls, tck, extrapolate=None):
        raise NotImplementedError("This method does not make sense for "
                "an Akima interpolator.")

    @classmethod
    def from_bernstein_basis(cls, bp, extrapolate=None):
        raise NotImplementedError("This method does not make sense for "
                "an Akima interpolator.")

    
