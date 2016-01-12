from __future__ import division, print_function, absolute_import

import numpy as np

from . import BPoly, PPoly
from .polyint import _isscalar
from scipy._lib._util import _asarray_validated


__all__ = ["PchipInterpolator", "pchip_interpolate", "pchip",
           "Akima1DInterpolator"]


class PchipInterpolator(BPoly):
    r"""PCHIP 1-d monotonic cubic interpolation.

    `x` and `y` are arrays of values used to approximate some function f,
    with ``y = f(x)``. The interpolant uses monotonic cubic splines
    to find the value of new points. (PCHIP stands for Piecewise Cubic
    Hermite Interpolating Polynomial).

    Parameters
    ----------
    x : ndarray
        A 1-D array of monotonically increasing real values.  `x` cannot
        include duplicate values (otherwise f is overspecified)
    y : ndarray
        A 1-D array of real values. `y`'s length along the interpolation
        axis must be equal to the length of `x`. If N-D array, use `axis`
        parameter to select correct axis.
    axis : int, optional
        Axis in the y array corresponding to the x-coordinate values.
    extrapolate : bool, optional
        Whether to extrapolate to out-of-bounds points based on first
        and last intervals, or to return NaNs.

    Methods
    -------
    __call__
    derivative
    antiderivative
    roots

    See Also
    --------
    Akima1DInterpolator

    Notes
    -----
    The interpolator preserves monotonicity in the interpolation data and does
    not overshoot if the data is not smooth.

    The first derivatives are guaranteed to be continuous, but the second
    derivatives may jump at :math:`x_k`.

    Determines the derivatives at the points :math:`x_k`, :math:`f'_k`,
    by using PCHIP algorithm [1]_.

    Let :math:`h_k = x_{k+1} - x_k`, and  :math:`d_k = (y_{k+1} - y_k) / h_k`
    are the slopes at internal points :math:`x_k`.
    If the signs of :math:`d_k` and :math:`d_{k-1}` are different or either of
    them equals zero, then :math:`f'_k = 0`. Otherwise, it is given by the
    weighted harmonic mean

    .. math::

        \frac{w_1 + w_2}{f'_k} = \frac{w_1}{d_{k-1}} + \frac{w_2}{d_k}

    where :math:`w_1 = 2 h_k + h_{k-1}` and :math:`w_2 = h_k + 2 h_{k-1}`.

    The end slopes are set using a one-sided scheme [2]_.


    References
    ----------
    .. [1] F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic Interpolation,
           SIAM J. Numer. Anal., 17(2), 238 (1980).
           DOI:10.1137/0717021
    .. [2] see, e.g., C. Moler, Numerical Computing with Matlab, 2004.
           DOI: http://dx.doi.org/10.1137/1.9780898717952


    """
    def __init__(self, x, y, axis=0, extrapolate=None):
        x = _asarray_validated(x, check_finite=False, as_inexact=True)
        y = _asarray_validated(y, check_finite=False, as_inexact=True)

        axis = axis % y.ndim

        xp = x.reshape((x.shape[0],) + (1,)*(y.ndim-1))
        yp = np.rollaxis(y, axis)

        dk = self._find_derivatives(xp, yp)
        data = np.hstack((yp[:, None, ...], dk[:, None, ...]))

        _b = BPoly.from_derivatives(x, data, orders=None)
        super(PchipInterpolator, self).__init__(_b.c, _b.x,
                                                extrapolate=extrapolate)
        self.axis = axis

    def roots(self):
        """
        Return the roots of the interpolated function.
        """
        return (PPoly.from_bernstein_basis(self._bpoly)).roots()

    @staticmethod
    def _edge_case(h0, h1, m0, m1):
        # one-sided three-point estimate for the derivative
        d = ((2*h0 + h1)*m0 - h0*m1) / (h0 + h1)

        # try to preserve shape
        mask = np.sign(d) != np.sign(m0)
        mask2 = (np.sign(m0) != np.sign(m1)) & (np.abs(d) > 3.*np.abs(m0))
        mmm = (~mask) & mask2

        d[mask] = 0.
        d[mmm] = 3.*m0[mmm]

        return d

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
            x = x[:, None]
            y = y[:, None]

        hk = x[1:] - x[:-1]
        mk = (y[1:] - y[:-1]) / hk
        smk = np.sign(mk)
        condition = (smk[1:] != smk[:-1]) | (mk[1:] == 0) | (mk[:-1] == 0)

        w1 = 2*hk[1:] + hk[:-1]
        w2 = hk[1:] + 2*hk[:-1]

        # values where division by zero occurs will be excluded
        # by 'condition' afterwards
        with np.errstate(divide='ignore'):
            whmean = (w1/mk[:-1] + w2/mk[1:]) / (w1 + w2)

        dk = np.zeros_like(y)
        dk[1:-1][condition] = 0.0
        dk[1:-1][~condition] = 1.0 / whmean[~condition]

        # special case endpoints, as suggested in 
        # Cleve Moler, Numerical Computing with MATLAB, Chap 3.4
        dk[0] = PchipInterpolator._edge_case(hk[0], hk[1], mk[0], mk[1])
        dk[-1] = PchipInterpolator._edge_case(hk[-1], hk[-2], mk[-1], mk[-2])

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
    der : int or list, optional
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
        N-D array of real values. The length of `y` along the first axis must
        be equal to the length of `x`.
    axis : int, optional
        Specifies the axis of `y` along which to interpolate. Interpolation
        defaults to the first axis of `y`.

    Methods
    -------
    __call__
    derivative
    antiderivative
    roots

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

    def __init__(self, x, y, axis=0):
        # Original implementation in MATLAB by N. Shamsundar (BSD licensed), see
        # http://www.mathworks.de/matlabcentral/fileexchange/1814-akima-interpolation
        x, y = map(np.asarray, (x, y))
        axis = axis % y.ndim

        if np.any(np.diff(x) < 0.):
            raise ValueError("x must be strictly ascending")
        if x.ndim != 1:
            raise ValueError("x must be 1-dimensional")
        if x.size < 2:
            raise ValueError("at least 2 breakpoints are needed")
        if x.size != y.shape[axis]:
            raise ValueError("x.shape must equal y.shape[%s]" % axis)

        # move interpolation axis to front
        y = np.rollaxis(y, axis)

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
        self.axis = axis

    def extend(self):
        raise NotImplementedError("Extending a 1D Akima interpolator is not "
                                  "yet implemented")

    # These are inherited from PPoly, but they do not produce an Akima
    # interpolator. Hence stub them out.
    @classmethod
    def from_spline(cls, tck, extrapolate=None):
        raise NotImplementedError("This method does not make sense for "
                                  "an Akima interpolator.")

    @classmethod
    def from_bernstein_basis(cls, bp, extrapolate=None):
        raise NotImplementedError("This method does not make sense for "
                                  "an Akima interpolator.")
