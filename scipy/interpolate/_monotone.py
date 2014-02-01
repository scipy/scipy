from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.misc import factorial

from . import polyint

__all__ = ["PchipInterpolator", "pchip_interpolate", "pchip"]

class PchipInterpolator(polyint.PiecewisePolynomial):
    """PCHIP 1-d monotonic cubic interpolation

    x and y are arrays of values used to approximate some function f,
    with ``y = f(x)``.  The interpolant uses monotonic cubic splines
    to find the value of new points.

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

    Notes
    -----
    Assumes x is sorted in monotonic order (e.g. ``x[1] > x[0]``).

    """
    def __init__(self, x, y, axis=0):
        x = np.asarray(x)
        y = np.asarray(y)

        axis = axis % y.ndim

        xp = x.reshape((x.shape[0],) + (1,)*(y.ndim-1))
        yp = np.rollaxis(y, axis)

        data = np.empty((yp.shape[0], 2) + yp.shape[1:], y.dtype)
        data[:,0] = yp
        data[:,1] = PchipInterpolator._find_derivatives(xp, yp)

        s = list(range(2, y.ndim + 1))
        s.insert(axis, 1)
        s.insert(axis, 0)
        data = data.transpose(s)

        PiecewisePolynomial.__init__(self, x, data, orders=3, direction=None,
                                     axis=axis)

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
        return P.derivative(x,der=der)
    else:
        return P.derivatives(x,der=np.amax(der)+1)[der]

# Backwards compatibility
pchip = PchipInterpolator
