
import numpy as np
from scipy.linalg import lstsq
from math import factorial
from scipy.ndimage import convolve1d


def savgol_coeffs(window_length, polyorder, deriv=0, delta=1.0):
    """Compute the coefficients for a 1-d Savitzky-Golay FIR filter.

    A symmetric window is assumed.

    Generally an odd integer is used for the window length, but the
    function allows an even value.  The phase delay of the causal filter
    is delta * 0.5 * (window_length - 1).  If an even window length is
    used, the filtered values are shifted by an extra half step, which
    might not be desirable.

    Parameters
    ----------
    window_length : int
        The length of the filter window (i.e. the number of coefficients).
        `window_length` must be a positive integer.
    polyorder : int
        The order of the polynomial used to fit the samples.
        `polyorder` must be less than `window_length`.
    deriv : int, optional
        The order of the derivative to compute.  This must be a
        nonnegative integer.  The default is 0, which means to filter
        the data without differentiating.
    delta : float, optional
        The spacing of the samples to which the filter will be applied.
        This is only used if deriv > 0.

    Returns
    -------
    c : 1-d ndarray
        The filter coefficients.

    Notes
    -----
    Use any convolution function  to apply the filter to an array
    (e.g. numpy.convolve, scipy.ndimage.convolve1d).

    References
    ----------
    A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of Data by
    Simplified Least Squares Procedures. Analytical Chemistry, 1964, 36 (8),
    pp 1627-1639.

    See Also
    --------
    savgol_filter

    Examples
    --------
    >>> savitzky_golay_fir(5, 2)
    array([-0.08571429,  0.34285714,  0.48571429,  0.34285714, -0.08571429])
    >>> savitzky_golay_fir(5, 2, deriv=1)
    array([  2.00000000e-01,   1.00000000e-01,   2.00607895e-16,
            -1.00000000e-01,  -2.00000000e-01])
    """

    if polyorder >= window_length:
        raise ValueError("polyorder must be less than window_length.")

    k = window_length // 2

    if window_length % 2 and polyorder == window_length - 1:
        # Trivial case; don't bother computing it.
        coeffs = np.zeros(window_length)
        coeffs[k] = 1.0
        return coeffs

    # Form the "design matrix" A.  The columns of A are powers
    # of the integers -k to k.  The powers (i.e. rows) range from
    # 0 to polyorder.
    x = np.linspace(-k, k, window_length)
    order = np.arange(polyorder + 1)
    A = x ** order[:, np.newaxis]

    # y determines which order derivative is returned.
    y = np.zeros(polyorder + 1)
    y[deriv] = (-1.0) ** deriv

    # Find the least-squares solution of A*c = y
    c, _, _, _ = lstsq(A, y)

    # Make the result 1-d, and adjust the coefficients for derivatives,
    # if needed.
    coeffs = factorial(deriv) * c.squeeze() / (delta ** deriv)

    return coeffs


def savgol_filter(x, window_length, polyorder, deriv=0, delta=1.0,
                          axis=-1, mode='constant', cval=0.0):
    """ Apply a Savitzky-Golay filter to an array.

    This is a 1-d filter.  If `x`  has dimension greater than 1, `axis`
    determines the axis along which the filter is applied.

    This function is a convenience wrapper around `savgol_coeffs()` and
    `ndimage.convolve1d()`.  The parameters `window_length`, `polyorder`
    and `deriv` are passed to `savitzky_golay_coeffs()`, and the parameters
    `mode` and `cval` are passed to `ndimage.convolve1d`.

    Parameters
    ----------
    x : array_like
        The data to be filtered.
    window_length : int
        The length of the filter window (i.e. the number of coefficients).
        `window_length` must be a positive integer.
    polyorder : int
        The order of the polynomial used to fit the samples.
        `polyorder` must be less than `window_length`.
    deriv : int, optional
        The order of the derivative to compute.  This must be a
        nonnegative integer.  The default is 0, which means to filter
        the data without differentiating.
    delta : float, optional
        The spacing of the samples to which the filter will be applied.
        This is only used if deriv > 0.  Default is 1.0.
    axis : int, optional
        The axis of the array `x` along which the filter is to be applied.
        Default is -1.
    padtype : str, optional
        Must be 'mirror', 'constant', or 'nearest'.  This
        determines the type of extension to use for the padded signal to
        which the filter is applied.  When `mode` is 'constant' (the
        default), the padding value is given by `cval`.  The default
        for `cval` is 0.0.
    cval : scalar, optional
        Value to fill past edges of input if `padtype` is 'constant'.
        Default is 0.0.

    Returns
    -------
    y : ndarray, same shape as `x`
        The filtered data.

    See Also
    --------
    savgol_coeffs

    """

    x = np.asarray(x)

    coeffs = savgol_coeffs(window_length, polyorder,
                                   deriv=deriv, delta=delta)

    # Use ndimage.convolve1d to apply the FIR filter.
    y = convolve1d(x, coeffs, axis=axis, mode=mode, cval=cval)

    return y
