# Author: Travis Oliphant
# 1999 -- 2002

from __future__ import division, print_function, absolute_import

import warnings

from . import sigtools
from scipy.lib.six import callable
from scipy import linalg
from scipy.fftpack import fft, ifft, ifftshift, fft2, ifft2, fftn, \
        ifftn, fftfreq
from numpy.fft import rfftn, irfftn
from numpy import polyadd, polymul, polydiv, polysub, roots, \
        poly, polyval, polyder, cast, asarray, isscalar, atleast_1d, \
        ones, real_if_close, zeros, array, arange, where, rank, \
        newaxis, product, ravel, sum, r_, iscomplexobj, take, \
        argsort, allclose, expand_dims, unique, prod, sort, reshape, \
        transpose, dot, mean, ndarray, atleast_2d
import numpy as np
from scipy.misc import factorial
from .windows import get_window
from ._arraytools import axis_slice, axis_reverse, odd_ext, even_ext, const_ext

__all__ = ['correlate', 'fftconvolve', 'convolve', 'convolve2d', 'correlate2d',
           'order_filter', 'medfilt', 'medfilt2d', 'wiener', 'lfilter',
           'lfiltic', 'deconvolve', 'hilbert', 'hilbert2', 'cmplx_sort',
           'unique_roots', 'invres', 'invresz', 'residue', 'residuez',
           'resample', 'detrend', 'lfilter_zi', 'filtfilt', 'decimate']


_modedict = {'valid': 0, 'same': 1, 'full': 2}

_boundarydict = {'fill': 0, 'pad': 0, 'wrap': 2, 'circular': 2, 'symm': 1,
                 'symmetric': 1, 'reflect': 4}


def _valfrommode(mode):
    try:
        val = _modedict[mode]
    except KeyError:
        if mode not in [0, 1, 2]:
            raise ValueError("Acceptable mode flags are 'valid' (0),"
                    " 'same' (1), or 'full' (2).")
        val = mode
    return val


def _bvalfromboundary(boundary):
    try:
        val = _boundarydict[boundary] << 2
    except KeyError:
        if val not in [0, 1, 2]:
            raise ValueError("Acceptable boundary flags are 'fill', 'wrap'"
                    " (or 'circular'), \n  and 'symm' (or 'symmetric').")
        val = boundary << 2
    return val


def _check_valid_mode_shapes(shape1, shape2):
    for d1, d2 in zip(shape1, shape2):
        if not d1 >= d2:
            raise ValueError(
                "in1 should have at least as many items as in2 in "
                "every dimension for 'valid' mode.")


def correlate(in1, in2, mode='full'):
    """
    Cross-correlate two N-dimensional arrays.

    Cross-correlate `in1` and `in2`, with the output size determined by the
    `mode` argument.

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`;
        if sizes of `in1` and `in2` are not equal then `in1` has to be the
        larger array.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear cross-correlation
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.

    Returns
    -------
    correlate : array
        An N-dimensional array containing a subset of the discrete linear
        cross-correlation of `in1` with `in2`.

    Notes
    -----
    The correlation z of two arrays x and y of rank d is defined as:

      z[...,k,...] = sum[..., i_l, ...]
                         x[..., i_l,...] * conj(y[..., i_l + k,...])

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    val = _valfrommode(mode)

    if rank(in1) == rank(in2) == 0:
        return in1 * in2
    elif not in1.ndim == in2.ndim:
        raise ValueError("in1 and in2 should have the same rank")

    if mode == 'valid':
        _check_valid_mode_shapes(in1.shape, in2.shape)
        ps = [i - j + 1 for i, j in zip(in1.shape, in2.shape)]
        out = np.empty(ps, in1.dtype)

        z = sigtools._correlateND(in1, in2, out, val)
    else:
        ps = [i + j - 1 for i, j in zip(in1.shape, in2.shape)]
        # zero pad input
        in1zpadded = np.zeros(ps, in1.dtype)
        sc = [slice(0, i) for i in in1.shape]
        in1zpadded[sc] = in1.copy()

        if mode == 'full':
            out = np.empty(ps, in1.dtype)
        elif mode == 'same':
            out = np.empty(in1.shape, in1.dtype)

        z = sigtools._correlateND(in1zpadded, in2, out, val)

    return z


def _centered(arr, newsize):
    # Return the center newsize portion of the array.
    newsize = asarray(newsize)
    currsize = array(arr.shape)
    startind = (currsize - newsize) // 2
    endind = startind + newsize
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


def fftconvolve(in1, in2, mode="full"):
    """Convolve two N-dimensional arrays using FFT.

    Convolve `in1` and `in2` using the fast Fourier transform method, with
    the output size determined by the `mode` argument.

    This is generally much faster than `convolve` for large arrays (n > ~500),
    but can be slower when only a few output values are needed, and can only
    output float arrays (int or object array inputs will be cast to float).

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`;
        if sizes of `in1` and `in2` are not equal then `in1` has to be the
        larger array.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.

    Returns
    -------
    out : array
        An N-dimensional array containing a subset of the discrete linear
        convolution of `in1` with `in2`.

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    if rank(in1) == rank(in2) == 0:  # scalar inputs
        return in1 * in2
    elif not in1.ndim == in2.ndim:
        raise ValueError("in1 and in2 should have the same rank")
    elif in1.size == 0 or in2.size == 0:  # empty arrays
        return array([])

    s1 = array(in1.shape)
    s2 = array(in2.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex) or
                      np.issubdtype(in2.dtype, np.complex))
    size = s1 + s2 - 1

    if mode == "valid":
        _check_valid_mode_shapes(s1, s2)

    # Always use 2**n-sized FFT
    fsize = 2 ** np.ceil(np.log2(size)).astype(int)
    fslice = tuple([slice(0, int(sz)) for sz in size])
    if not complex_result:
        ret = irfftn(rfftn(in1, fsize) *
                     rfftn(in2, fsize), fsize)[fslice].copy()
        ret = ret.real
    else:
        ret = ifftn(fftn(in1, fsize) * fftn(in2, fsize))[fslice].copy()

    if mode == "full":
        return ret
    elif mode == "same":
        return _centered(ret, s1)
    elif mode == "valid":
        return _centered(ret, s1 - s2 + 1)


def convolve(in1, in2, mode='full'):
    """
    Convolve two N-dimensional arrays.

    Convolve `in1` and `in2`, with the output size determined by the
    `mode` argument.

    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`;
        if sizes of `in1` and `in2` are not equal then `in1` has to be the
        larger array.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.

    Returns
    -------
    convolve : array
        An N-dimensional array containing a subset of the discrete linear
        convolution of `in1` with `in2`.

    """
    volume = asarray(in1)
    kernel = asarray(in2)

    if rank(volume) == rank(kernel) == 0:
        return volume * kernel

    slice_obj = [slice(None, None, -1)] * len(kernel.shape)

    if np.iscomplexobj(kernel):
        return correlate(volume, kernel[slice_obj].conj(), mode)
    else:
        return correlate(volume, kernel[slice_obj], mode)


def order_filter(a, domain, rank):
    """
    Perform an order filter on an N-dimensional array.

    Perform an order filter on the array in.  The domain argument acts as a
    mask centered over each pixel.  The non-zero elements of domain are
    used to select elements surrounding each input pixel which are placed
    in a list.   The list is sorted, and the output for that pixel is the
    element corresponding to rank in the sorted list.

    Parameters
    ----------
    a : ndarray
        The N-dimensional input array.
    domain : array_like
        A mask array with the same number of dimensions as `in`.
        Each dimension should have an odd number of elements.
    rank : int
        A non-negative integer which selects the element from the
        sorted list (0 corresponds to the smallest element, 1 is the
        next smallest element, etc.).

    Returns
    -------
    out : ndarray
        The results of the order filter in an array with the same
        shape as `in`.

    Examples
    --------
    >>> from scipy import signal
    >>> x = np.arange(25).reshape(5, 5)
    >>> domain = np.identity(3)
    >>> x
    array([[ 0,  1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14],
           [15, 16, 17, 18, 19],
           [20, 21, 22, 23, 24]])
    >>> signal.order_filter(x, domain, 0)
    array([[  0.,   0.,   0.,   0.,   0.],
           [  0.,   0.,   1.,   2.,   0.],
           [  0.,   5.,   6.,   7.,   0.],
           [  0.,  10.,  11.,  12.,   0.],
           [  0.,   0.,   0.,   0.,   0.]])
    >>> signal.order_filter(x, domain, 2)
    array([[  6.,   7.,   8.,   9.,   4.],
           [ 11.,  12.,  13.,  14.,   9.],
           [ 16.,  17.,  18.,  19.,  14.],
           [ 21.,  22.,  23.,  24.,  19.],
           [ 20.,  21.,  22.,  23.,  24.]])

    """
    domain = asarray(domain)
    size = domain.shape
    for k in range(len(size)):
        if (size[k] % 2) != 1:
            raise ValueError("Each dimension of domain argument "
                    " should have an odd number of elements.")
    return sigtools._order_filterND(a, domain, rank)


def medfilt(volume, kernel_size=None):
    """
    Perform a median filter on an N-dimensional array.

    Apply a median filter to the input array using a local window-size
    given by `kernel_size`.

    Parameters
    ----------
    volume : array_like
        An N-dimensional input array.
    kernel_size : array_like, optional
        A scalar or an N-length list giving the size of the median filter
        window in each dimension.  Elements of `kernel_size` should be odd.
        If `kernel_size` is a scalar, then this scalar is used as the size in
        each dimension. Default size is 3 for each dimension.

    Returns
    -------
    out : ndarray
        An array the same size as input containing the median filtered
        result.

    """
    volume = atleast_1d(volume)
    if kernel_size is None:
        kernel_size = [3] * len(volume.shape)
    kernel_size = asarray(kernel_size)
    if kernel_size.shape == ():
        kernel_size = np.repeat(kernel_size.item(), volume.ndim)

    for k in range(len(volume.shape)):
        if (kernel_size[k] % 2) != 1:
            raise ValueError("Each element of kernel_size should be odd.")

    domain = ones(kernel_size)

    numels = product(kernel_size, axis=0)
    order = numels // 2
    return sigtools._order_filterND(volume, domain, order)


def wiener(im, mysize=None, noise=None):
    """
    Perform a Wiener filter on an N-dimensional array.

    Apply a Wiener filter to the N-dimensional array `im`.

    Parameters
    ----------
    im : ndarray
        An N-dimensional array.
    mysize : int or arraylike, optional
        A scalar or an N-length list giving the size of the Wiener filter
        window in each dimension.  Elements of mysize should be odd.
        If mysize is a scalar, then this scalar is used as the size
        in each dimension.
    noise : float, optional
        The noise-power to use. If None, then noise is estimated as the
        average of the local variance of the input.

    Returns
    -------
    out : ndarray
        Wiener filtered result with the same shape as `im`.

    """
    im = asarray(im)
    if mysize is None:
        mysize = [3] * len(im.shape)
    mysize = asarray(mysize)
    if mysize.shape == ():
        mysize = np.repeat(mysize.item(), im.ndim)

    # Estimate the local mean
    lMean = correlate(im, ones(mysize), 'same') / product(mysize, axis=0)

    # Estimate the local variance
    lVar = (correlate(im ** 2, ones(mysize), 'same') / product(mysize, axis=0)
            - lMean ** 2)

    # Estimate the noise power if needed.
    if noise is None:
        noise = mean(ravel(lVar), axis=0)

    res = (im - lMean)
    res *= (1 - noise / lVar)
    res += lMean
    out = where(lVar < noise, lMean, res)

    return out


def convolve2d(in1, in2, mode='full', boundary='fill', fillvalue=0):
    """
    Convolve two 2-dimensional arrays.

    Convolve `in1` and `in2` with output size determined by `mode`, and
    boundary conditions determined by `boundary` and `fillvalue`.

    Parameters
    ----------
    in1, in2 : array_like
        Two-dimensional input arrays to be convolved.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.

    boundary : str {'fill', 'wrap', 'symm'}, optional
        A flag indicating how to handle boundaries:

        ``fill``
           pad input arrays with fillvalue. (default)
        ``wrap``
           circular boundary conditions.
        ``symm``
           symmetrical boundary conditions.

    fillvalue : scalar, optional
        Value to fill pad input arrays with. Default is 0.

    Returns
    -------
    out : ndarray
        A 2-dimensional array containing a subset of the discrete linear
        convolution of `in1` with `in2`.

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    if mode == 'valid':
        _check_valid_mode_shapes(in1.shape, in2.shape)

    val = _valfrommode(mode)
    bval = _bvalfromboundary(boundary)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.ComplexWarning)
        # FIXME: some cast generates a warning here
        out = sigtools._convolve2d(in1, in2, 1, val, bval, fillvalue)

    return out


def correlate2d(in1, in2, mode='full', boundary='fill', fillvalue=0):
    """
    Cross-correlate two 2-dimensional arrays.

    Cross correlate `in1` and `in2` with output size determined by `mode`, and
    boundary conditions determined by `boundary` and `fillvalue`.

    Parameters
    ----------
    in1, in2 : array_like
        Two-dimensional input arrays to be convolved.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:

        ``full``
           The output is the full discrete linear cross-correlation
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.

    boundary : str {'fill', 'wrap', 'symm'}, optional
        A flag indicating how to handle boundaries:

        ``fill``
           pad input arrays with fillvalue. (default)
        ``wrap``
           circular boundary conditions.
        ``symm``
           symmetrical boundary conditions.

    fillvalue : scalar, optional
        Value to fill pad input arrays with. Default is 0.

    Returns
    -------
    correlate2d : ndarray
        A 2-dimensional array containing a subset of the discrete linear
        cross-correlation of `in1` with `in2`.

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    if mode == 'valid':
        _check_valid_mode_shapes(in1.shape, in2.shape)

    val = _valfrommode(mode)
    bval = _bvalfromboundary(boundary)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.ComplexWarning)
        # FIXME: some cast generates a warning here
        out = sigtools._convolve2d(in1, in2, 0, val, bval, fillvalue)

    return out


def medfilt2d(input, kernel_size=3):
    """
    Median filter a 2-dimensional array.

    Apply a median filter to the `input` array using a local window-size
    given by `kernel_size` (must be odd).

    Parameters
    ----------
    input : array_like
        A 2-dimensional input array.
    kernel_size : array_like, optional
        A scalar or a list of length 2, giving the size of the
        median filter window in each dimension.  Elements of
        `kernel_size` should be odd.  If `kernel_size` is a scalar,
        then this scalar is used as the size in each dimension.
        Default is a kernel of size (3, 3).

    Returns
    -------
    out : ndarray
        An array the same size as input containing the median filtered
        result.

    """
    image = asarray(input)
    if kernel_size is None:
        kernel_size = [3] * 2
    kernel_size = asarray(kernel_size)
    if kernel_size.shape == ():
        kernel_size = np.repeat(kernel_size.item(), 2)

    for size in kernel_size:
        if (size % 2) != 1:
            raise ValueError("Each element of kernel_size should be odd.")

    return sigtools._medfilt2d(image, kernel_size)


def lfilter(b, a, x, axis=-1, zi=None):
    """
    Filter data along one-dimension with an IIR or FIR filter.

    Filter a data sequence, `x`, using a digital filter.  This works for many
    fundamental data types (including Object type).  The filter is a direct
    form II transposed implementation of the standard difference equation
    (see Notes).

    Parameters
    ----------
    b : array_like
        The numerator coefficient vector in a 1-D sequence.
    a : array_like
        The denominator coefficient vector in a 1-D sequence.  If ``a[0]``
        is not 1, then both `a` and `b` are normalized by ``a[0]``.
    x : array_like
        An N-dimensional input array.
    axis : int
        The axis of the input data array along which to apply the
        linear filter. The filter is applied to each subarray along
        this axis.  Default is -1.
    zi : array_like, optional
        Initial conditions for the filter delays.  It is a vector
        (or array of vectors for an N-dimensional input) of length
        ``max(len(a),len(b))-1``.  If `zi` is None or is not given then
        initial rest is assumed.  See `lfiltic` for more information.

    Returns
    -------
    y : array
        The output of the digital filter.
    zf : array, optional
        If `zi` is None, this is not returned, otherwise, `zf` holds the
        final filter delay values.

    Notes
    -----
    The filter function is implemented as a direct II transposed structure.
    This means that the filter implements::

       a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[nb]*x[n-nb]
                               - a[1]*y[n-1] - ... - a[na]*y[n-na]

    using the following difference equations::

         y[m] = b[0]*x[m] + z[0,m-1]
         z[0,m] = b[1]*x[m] + z[1,m-1] - a[1]*y[m]
         ...
         z[n-3,m] = b[n-2]*x[m] + z[n-2,m-1] - a[n-2]*y[m]
         z[n-2,m] = b[n-1]*x[m] - a[n-1]*y[m]

    where m is the output sample number and n=max(len(a),len(b)) is the
    model order.

    The rational transfer function describing this filter in the
    z-transform domain is::

                             -1               -nb
                 b[0] + b[1]z  + ... + b[nb] z
         Y(z) = ---------------------------------- X(z)
                             -1               -na
                 a[0] + a[1]z  + ... + a[na] z

    """
    if isscalar(a):
        a = [a]
    if zi is None:
        return sigtools._linear_filter(b, a, x, axis)
    else:
        return sigtools._linear_filter(b, a, x, axis, zi)


def lfiltic(b, a, y, x=None):
    """
    Construct initial conditions for lfilter.

    Given a linear filter (b, a) and initial conditions on the output `y`
    and the input `x`, return the inital conditions on the state vector zi
    which is used by `lfilter` to generate the output given the input.

    Parameters
    ----------
    b : array_like
        Linear filter term.
    a : array_like
        Linear filter term.
    y : array_like
        Initial conditions.

        If ``N=len(a) - 1``, then ``y = {y[-1], y[-2], ..., y[-N]}``.

        If `y` is too short, it is padded with zeros.
    x : array_like, optional
        Initial conditions.

        If ``M=len(b) - 1``, then ``x = {x[-1], x[-2], ..., x[-M]}``.

        If `x` is not given, its initial conditions are assumed zero.

        If `x` is too short, it is padded with zeros.

    Returns
    -------
    zi : ndarray
        The state vector ``zi``.
        ``zi = {z_0[-1], z_1[-1], ..., z_K-1[-1]}``, where ``K = max(M,N)``.

    See Also
    --------
    lfilter

    """
    N = np.size(a) - 1
    M = np.size(b) - 1
    K = max(M, N)
    y = asarray(y)
    zi = zeros(K, y.dtype.char)
    if x is None:
        x = zeros(M, y.dtype.char)
    else:
        x = asarray(x)
        L = np.size(x)
        if L < M:
            x = r_[x, zeros(M - L)]
    L = np.size(y)
    if L < N:
        y = r_[y, zeros(N - L)]

    for m in range(M):
        zi[m] = sum(b[m + 1:] * x[:M - m], axis=0)

    for m in range(N):
        zi[m] -= sum(a[m + 1:] * y[:N - m], axis=0)

    return zi


def deconvolve(signal, divisor):
    """Deconvolves `divisor` out of `signal`.

    Parameters
    ----------
    signal : array
        Signal input
    divisor : array
        Divisor input

    Returns
    -------
    q : array
        Quotient of the division
    r : array
        Remainder

    Examples
    --------
    >>> from scipy import signal
    >>> sig = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1,])
    >>> filter = np.array([1,1,0])
    >>> res = signal.convolve(sig, filter)
    >>> signal.deconvolve(res, filter)
    (array([ 0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.]),
     array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]))

    """
    num = atleast_1d(signal)
    den = atleast_1d(divisor)
    N = len(num)
    D = len(den)
    if D > N:
        quot = []
        rem = num
    else:
        input = ones(N - D + 1, float)
        input[1:] = 0
        quot = lfilter(num, den, input)
        rem = num - convolve(den, quot, mode='full')
    return quot, rem


def hilbert(x, N=None, axis=-1):
    """
    Compute the analytic signal, using the Hilbert transform.

    The transformation is done along the last axis by default.

    Parameters
    ----------
    x : array_like
        Signal data.  Must be real.
    N : int, optional
        Number of Fourier components.  Default: ``x.shape[axis]``
    axis : int, optional
        Axis along which to do the transformation.  Default: -1.

    Returns
    -------
    xa : ndarray
        Analytic signal of `x`, of each 1-D array along `axis`

    Notes
    -----
    The analytic signal ``x_a(t)`` of signal ``x(t)`` is:

    .. math:: x_a = F^{-1}(F(x) 2U) = x + i y

    where `F` is the Fourier transform, `U` the unit step function,
    and `y` the Hilbert transform of `x`. [1]_

    In other words, the negative half of the frequency spectrum is zeroed
    out, turning the real-valued signal into a complex signal.  The Hilbert
    transformed signal can be obtained from ``np.imag(hilbert(x))``, and the
    original signal from ``np.real(hilbert(x))``.

    References
    ----------
    .. [1] Wikipedia, "Analytic signal".
           http://en.wikipedia.org/wiki/Analytic_signal

    """
    x = asarray(x)
    if iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    Xf = fft(x, N, axis=axis)
    h = zeros(N)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2

    if len(x.shape) > 1:
        ind = [newaxis] * x.ndim
        ind[axis] = slice(None)
        h = h[ind]
    x = ifft(Xf * h, axis=axis)
    return x


def hilbert2(x, N=None):
    """
    Compute the '2-D' analytic signal of `x`

    Parameters
    ----------
    x : array_like
        2-D signal data.
    N : int or tuple of two ints, optional
        Number of Fourier components. Default is ``x.shape``

    Returns
    -------
    xa : ndarray
        Analytic signal of `x` taken along axes (0,1).

    References
    ----------
    .. [1] Wikipedia, "Analytic signal",
        http://en.wikipedia.org/wiki/Analytic_signal

    """
    x = atleast_2d(x)
    if len(x.shape) > 2:
        raise ValueError("x must be rank 2.")
    if iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape
    elif isinstance(N, int):
        if N <= 0:
            raise ValueError("N must be positive.")
        N = (N, N)
    elif len(N) != 2 or np.any(np.asarray(N) <= 0):
        raise ValueError("When given as a tuple, N must hold exactly "
                         "two positive integers")

    Xf = fft2(x, N, axes=(0, 1))
    h1 = zeros(N[0], 'd')
    h2 = zeros(N[1], 'd')
    for p in range(2):
        h = eval("h%d" % (p + 1))
        N1 = N[p]
        if N1 % 2 == 0:
            h[0] = h[N1 // 2] = 1
            h[1:N1 // 2] = 2
        else:
            h[0] = 1
            h[1:(N1 + 1) // 2] = 2
        exec("h%d = h" % (p + 1), globals(), locals())

    h = h1[:, newaxis] * h2[newaxis, :]
    k = len(x.shape)
    while k > 2:
        h = h[:, newaxis]
        k -= 1
    x = ifft2(Xf * h, axes=(0, 1))
    return x


def cmplx_sort(p):
    "sort roots based on magnitude."
    p = asarray(p)
    if iscomplexobj(p):
        indx = argsort(abs(p))
    else:
        indx = argsort(p)
    return take(p, indx, 0), indx


def unique_roots(p, tol=1e-3, rtype='min'):
    """
    Determine unique roots and their multiplicities from a list of roots.

    Parameters
    ----------
    p : array_like
        The list of roots.
    tol : float, optional
        The tolerance for two roots to be considered equal. Default is 1e-3.
    rtype : {'max', 'min, 'avg'}, optional
        How to determine the returned root if multiple roots are within
        `tol` of each other.

          - 'max': pick the maximum of those roots.
          - 'min': pick the minimum of those roots.
          - 'avg': take the average of those roots.

    Returns
    -------
    pout : ndarray
        The list of unique roots, sorted from low to high.
    mult : ndarray
        The multiplicity of each root.

    Notes
    -----
    This utility function is not specific to roots but can be used for any
    sequence of values for which uniqueness and multiplicity has to be
    determined. For a more general routine, see `numpy.unique`.

    Examples
    --------
    >>> from scipy import signal
    >>> vals = [0, 1.3, 1.31, 2.8, 1.25, 2.2, 10.3]
    >>> uniq, mult = signal.unique_roots(vals, tol=2e-2, rtype='avg')

    Check which roots have multiplicity larger than 1:

    >>> uniq[mult > 1]
    array([ 1.305])

    """
    if rtype in ['max', 'maximum']:
        comproot = np.max
    elif rtype in ['min', 'minimum']:
        comproot = np.min
    elif rtype in ['avg', 'mean']:
        comproot = np.mean
    else:
        raise ValueError("`rtype` must be one of "
                         "{'max', 'maximum', 'min', 'minimum', 'avg', 'mean'}")
    p = asarray(p) * 1.0
    tol = abs(tol)
    p, indx = cmplx_sort(p)
    pout = []
    mult = []
    indx = -1
    curp = p[0] + 5 * tol
    sameroots = []
    for k in range(len(p)):
        tr = p[k]
        if abs(tr - curp) < tol:
            sameroots.append(tr)
            curp = comproot(sameroots)
            pout[indx] = curp
            mult[indx] += 1
        else:
            pout.append(tr)
            curp = tr
            sameroots = [tr]
            indx += 1
            mult.append(1)
    return array(pout), array(mult)


def invres(r, p, k, tol=1e-3, rtype='avg'):
    """
    Compute b(s) and a(s) from partial fraction expansion: r,p,k

    If ``M = len(b)`` and ``N = len(a)``::

                b(s)     b[0] x**(M-1) + b[1] x**(M-2) + ... + b[M-1]
        H(s) = ------ = ----------------------------------------------
                a(s)     a[0] x**(N-1) + a[1] x**(N-2) + ... + a[N-1]

                 r[0]       r[1]             r[-1]
             = -------- + -------- + ... + --------- + k(s)
               (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like::

          r[i]      r[i+1]              r[i+n-1]
        -------- + ----------- + ... + -----------
        (s-p[i])  (s-p[i])**2          (s-p[i])**n

    Parameters
    ----------
    r : ndarray
        Residues.
    p : ndarray
        Poles.
    k : ndarray
        Coefficients of the direct polynomial term.
    tol : float, optional
        The tolerance for two roots to be considered equal. Default is 1e-3.
    rtype : {'max', 'min, 'avg'}, optional
        How to determine the returned root if multiple roots are within
        `tol` of each other.

          'max': pick the maximum of those roots.

          'min': pick the minimum of those roots.

          'avg': take the average of those roots.

    See Also
    --------
    residue, unique_roots

    """
    extra = k
    p, indx = cmplx_sort(p)
    r = take(r, indx, 0)
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for k in range(len(pout)):
        p.extend([pout[k]] * mult[k])
    a = atleast_1d(poly(p))
    if len(extra) > 0:
        b = polymul(extra, a)
    else:
        b = [0]
    indx = 0
    for k in range(len(pout)):
        temp = []
        for l in range(len(pout)):
            if l != k:
                temp.extend([pout[l]] * mult[l])
        for m in range(mult[k]):
            t2 = temp[:]
            t2.extend([pout[k]] * (mult[k] - m - 1))
            b = polyadd(b, r[indx] * poly(t2))
            indx += 1
    b = real_if_close(b)
    while allclose(b[0], 0, rtol=1e-14) and (b.shape[-1] > 1):
        b = b[1:]
    return b, a


def residue(b, a, tol=1e-3, rtype='avg'):
    """
    Compute partial-fraction expansion of b(s) / a(s).

    If ``M = len(b)`` and ``N = len(a)``, then the partial-fraction
    expansion H(s) is defined as::

              b(s)     b[0] s**(M-1) + b[1] s**(M-2) + ... + b[M-1]
      H(s) = ------ = ----------------------------------------------
              a(s)     a[0] s**(N-1) + a[1] s**(N-2) + ... + a[N-1]

               r[0]       r[1]             r[-1]
           = -------- + -------- + ... + --------- + k(s)
             (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer together than `tol`), then H(s)
    has terms like::

            r[i]      r[i+1]              r[i+n-1]
          -------- + ----------- + ... + -----------
          (s-p[i])  (s-p[i])**2          (s-p[i])**n

    Returns
    -------
    r : ndarray
        Residues.
    p : ndarray
        Poles.
    k : ndarray
        Coefficients of the direct polynomial term.

    See Also
    --------
    invres, numpy.poly, unique_roots

    """

    b, a = map(asarray, (b, a))
    rscale = a[0]
    k, b = polydiv(b, a)
    p = roots(a)
    r = p * 0.0
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for n in range(len(pout)):
        p.extend([pout[n]] * mult[n])
    p = asarray(p)
    # Compute the residue from the general formula
    indx = 0
    for n in range(len(pout)):
        bn = b.copy()
        pn = []
        for l in range(len(pout)):
            if l != n:
                pn.extend([pout[l]] * mult[l])
        an = atleast_1d(poly(pn))
        # bn(s) / an(s) is (s-po[n])**Nn * b(s) / a(s) where Nn is
        # multiplicity of pole at po[n]
        sig = mult[n]
        for m in range(sig, 0, -1):
            if sig > m:
                # compute next derivative of bn(s) / an(s)
                term1 = polymul(polyder(bn, 1), an)
                term2 = polymul(bn, polyder(an, 1))
                bn = polysub(term1, term2)
                an = polymul(an, an)
            r[indx + m - 1] = polyval(bn, pout[n]) / polyval(an, pout[n]) \
                          / factorial(sig - m)
        indx += sig
    return r / rscale, p, k


def residuez(b, a, tol=1e-3, rtype='avg'):
    """
    Compute partial-fraction expansion of b(z) / a(z).

    If ``M = len(b)`` and ``N = len(a)``::

                b(z)     b[0] + b[1] z**(-1) + ... + b[M-1] z**(-M+1)
        H(z) = ------ = ----------------------------------------------
                a(z)     a[0] + a[1] z**(-1) + ... + a[N-1] z**(-N+1)

                 r[0]                   r[-1]
         = --------------- + ... + ---------------- + k[0] + k[1]z**(-1) ...
           (1-p[0]z**(-1))         (1-p[-1]z**(-1))

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like::

             r[i]              r[i+1]                    r[i+n-1]
        -------------- + ------------------ + ... + ------------------
        (1-p[i]z**(-1))  (1-p[i]z**(-1))**2         (1-p[i]z**(-1))**n

    See also
    --------
    invresz, unique_roots

    """
    b, a = map(asarray, (b, a))
    gain = a[0]
    brev, arev = b[::-1], a[::-1]
    krev, brev = polydiv(brev, arev)
    if krev == []:
        k = []
    else:
        k = krev[::-1]
    b = brev[::-1]
    p = roots(a)
    r = p * 0.0
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for n in range(len(pout)):
        p.extend([pout[n]] * mult[n])
    p = asarray(p)
    # Compute the residue from the general formula (for discrete-time)
    #  the polynomial is in z**(-1) and the multiplication is by terms
    #  like this (1-p[i] z**(-1))**mult[i].  After differentiation,
    #  we must divide by (-p[i])**(m-k) as well as (m-k)!
    indx = 0
    for n in range(len(pout)):
        bn = brev.copy()
        pn = []
        for l in range(len(pout)):
            if l != n:
                pn.extend([pout[l]] * mult[l])
        an = atleast_1d(poly(pn))[::-1]
        # bn(z) / an(z) is (1-po[n] z**(-1))**Nn * b(z) / a(z) where Nn is
        # multiplicity of pole at po[n] and b(z) and a(z) are polynomials.
        sig = mult[n]
        for m in range(sig, 0, -1):
            if sig > m:
                # compute next derivative of bn(s) / an(s)
                term1 = polymul(polyder(bn, 1), an)
                term2 = polymul(bn, polyder(an, 1))
                bn = polysub(term1, term2)
                an = polymul(an, an)
            r[indx + m - 1] = (polyval(bn, 1.0 / pout[n]) /
                               polyval(an, 1.0 / pout[n]) /
                               factorial(sig - m) / (-pout[n]) ** (sig - m))
        indx += sig
    return r / gain, p, k


def invresz(r, p, k, tol=1e-3, rtype='avg'):
    """
    Compute b(z) and a(z) from partial fraction expansion: r,p,k

    If ``M = len(b)`` and ``N = len(a)``::

                b(z)     b[0] + b[1] z**(-1) + ... + b[M-1] z**(-M+1)
        H(z) = ------ = ----------------------------------------------
                a(z)     a[0] + a[1] z**(-1) + ... + a[N-1] z**(-N+1)

                     r[0]                   r[-1]
             = --------------- + ... + ---------------- + k[0] + k[1]z**(-1) ...
               (1-p[0]z**(-1))         (1-p[-1]z**(-1))

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like::

             r[i]              r[i+1]                    r[i+n-1]
        -------------- + ------------------ + ... + ------------------
        (1-p[i]z**(-1))  (1-p[i]z**(-1))**2         (1-p[i]z**(-1))**n

    See Also
    --------
    residuez, unique_roots

    """
    extra = asarray(k)
    p, indx = cmplx_sort(p)
    r = take(r, indx, 0)
    pout, mult = unique_roots(p, tol=tol, rtype=rtype)
    p = []
    for k in range(len(pout)):
        p.extend([pout[k]] * mult[k])
    a = atleast_1d(poly(p))
    if len(extra) > 0:
        b = polymul(extra, a)
    else:
        b = [0]
    indx = 0
    brev = asarray(b)[::-1]
    for k in range(len(pout)):
        temp = []
        # Construct polynomial which does not include any of this root
        for l in range(len(pout)):
            if l != k:
                temp.extend([pout[l]] * mult[l])
        for m in range(mult[k]):
            t2 = temp[:]
            t2.extend([pout[k]] * (mult[k] - m - 1))
            brev = polyadd(brev, (r[indx] * poly(t2))[::-1])
            indx += 1
    b = real_if_close(brev[::-1])
    return b, a


def resample(x, num, t=None, axis=0, window=None):
    """
    Resample `x` to `num` samples using Fourier method along the given axis.

    The resampled signal starts at the same value as `x` but is sampled
    with a spacing of ``len(x) / num * (spacing of x)``.  Because a
    Fourier method is used, the signal is assumed to be periodic.

    Parameters
    ----------
    x : array_like
        The data to be resampled.
    num : int
        The number of samples in the resampled signal.
    t : array_like, optional
        If `t` is given, it is assumed to be the sample positions
        associated with the signal data in `x`.
    axis : int, optional
        The axis of `x` that is resampled.  Default is 0.
    window : array_like, callable, string, float, or tuple, optional
        Specifies the window applied to the signal in the Fourier
        domain.  See below for details.

    Returns
    -------
    resampled_x or (resampled_x, resampled_t)
        Either the resampled array, or, if `t` was given, a tuple
        containing the resampled array and the corresponding resampled
        positions.

    Notes
    -----
    The argument `window` controls a Fourier-domain window that tapers
    the Fourier spectrum before zero-padding to alleviate ringing in
    the resampled values for sampled signals you didn't intend to be
    interpreted as band-limited.

    If `window` is a function, then it is called with a vector of inputs
    indicating the frequency bins (i.e. fftfreq(x.shape[axis]) ).

    If `window` is an array of the same length as `x.shape[axis]` it is
    assumed to be the window to be applied directly in the Fourier
    domain (with dc and low-frequency first).

    For any other type of `window`, the function `scipy.signal.get_window`
    is called to generate the window.

    The first sample of the returned vector is the same as the first
    sample of the input vector.  The spacing between samples is changed
    from dx to:

        dx * len(x) / num

    If `t` is not None, then it represents the old sample positions,
    and the new sample positions will be returned as well as the new
    samples.

    """
    x = asarray(x)
    X = fft(x, axis=axis)
    Nx = x.shape[axis]
    if window is not None:
        if callable(window):
            W = window(fftfreq(Nx))
        elif isinstance(window, ndarray) and window.shape == (Nx,):
            W = window
        else:
            W = ifftshift(get_window(window, Nx))
        newshape = ones(len(x.shape))
        newshape[axis] = len(W)
        W.shape = newshape
        X = X * W
    sl = [slice(None)] * len(x.shape)
    newshape = list(x.shape)
    newshape[axis] = num
    N = int(np.minimum(num, Nx))
    Y = zeros(newshape, 'D')
    sl[axis] = slice(0, (N + 1) // 2)
    Y[sl] = X[sl]
    sl[axis] = slice(-(N - 1) // 2, None)
    Y[sl] = X[sl]
    y = ifft(Y, axis=axis) * (float(num) / float(Nx))

    if x.dtype.char not in ['F', 'D']:
        y = y.real

    if t is None:
        return y
    else:
        new_t = arange(0, num) * (t[1] - t[0]) * Nx / float(num) + t[0]
        return y, new_t


def detrend(data, axis=-1, type='linear', bp=0):
    """
    Remove linear trend along axis from data.

    Parameters
    ----------
    data : array_like
        The input data.
    axis : int, optional
        The axis along which to detrend the data. By default this is the
        last axis (-1).
    type : {'linear', 'constant'}, optional
        The type of detrending. If ``type == 'linear'`` (default),
        the result of a linear least-squares fit to `data` is subtracted
        from `data`.
        If ``type == 'constant'``, only the mean of `data` is subtracted.
    bp : array_like of ints, optional
        A sequence of break points. If given, an individual linear fit is
        performed for each part of `data` between two break points.
        Break points are specified as indices into `data`.

    Returns
    -------
    ret : ndarray
        The detrended input data.

    Examples
    --------
    >>> from scipy import signal
    >>> randgen = np.random.RandomState(9)
    >>> npoints = 1e3
    >>> noise = randgen.randn(npoints)
    >>> x = 3 + 2*np.linspace(0, 1, npoints) + noise
    >>> (signal.detrend(x) - noise).max() < 0.01
    True

    """
    if type not in ['linear', 'l', 'constant', 'c']:
        raise ValueError("Trend type must be 'linear' or 'constant'.")
    data = asarray(data)
    dtype = data.dtype.char
    if dtype not in 'dfDF':
        dtype = 'd'
    if type in ['constant', 'c']:
        ret = data - expand_dims(mean(data, axis), axis)
        return ret
    else:
        dshape = data.shape
        N = dshape[axis]
        bp = sort(unique(r_[0, bp, N]))
        if np.any(bp > N):
            raise ValueError("Breakpoints must be less than length "
                    "of data along given axis.")
        Nreg = len(bp) - 1
        # Restructure data so that axis is along first dimension and
        #  all other dimensions are collapsed into second dimension
        rnk = len(dshape)
        if axis < 0:
            axis = axis + rnk
        newdims = r_[axis, 0:axis, axis + 1:rnk]
        newdata = reshape(transpose(data, tuple(newdims)),
                          (N, prod(dshape, axis=0) // N))
        newdata = newdata.copy()  # make sure we have a copy
        if newdata.dtype.char not in 'dfDF':
            newdata = newdata.astype(dtype)
        # Find leastsq fit and remove it for each piece
        for m in range(Nreg):
            Npts = bp[m + 1] - bp[m]
            A = ones((Npts, 2), dtype)
            A[:, 0] = cast[dtype](arange(1, Npts + 1) * 1.0 / Npts)
            sl = slice(bp[m], bp[m + 1])
            coef, resids, rank, s = linalg.lstsq(A, newdata[sl])
            newdata[sl] = newdata[sl] - dot(A, coef)
        # Put data back in original shape.
        tdshape = take(dshape, newdims, 0)
        ret = reshape(newdata, tuple(tdshape))
        vals = list(range(1, rnk))
        olddims = vals[:axis] + [0] + vals[axis:]
        ret = transpose(ret, tuple(olddims))
        return ret


def lfilter_zi(b, a):
    """
    Compute an initial state `zi` for the lfilter function that corresponds
    to the steady state of the step response.

    A typical use of this function is to set the initial state so that the
    output of the filter starts at the same value as the first element of
    the signal to be filtered.

    Parameters
    ----------
    b, a : array_like (1-D)
        The IIR filter coefficients. See `lfilter` for more
        information.

    Returns
    -------
    zi : 1-D ndarray
        The initial state for the filter.

    Notes
    -----
    A linear filter with order m has a state space representation (A, B, C, D),
    for which the output y of the filter can be expressed as::

        z(n+1) = A*z(n) + B*x(n)
        y(n)   = C*z(n) + D*x(n)

    where z(n) is a vector of length m, A has shape (m, m), B has shape
    (m, 1), C has shape (1, m) and D has shape (1, 1) (assuming x(n) is
    a scalar).  lfilter_zi solves::

        zi = A*zi + B

    In other words, it finds the initial condition for which the response
    to an input of all ones is a constant.

    Given the filter coefficients `a` and `b`, the state space matrices
    for the transposed direct form II implementation of the linear filter,
    which is the implementation used by scipy.signal.lfilter, are::

        A = scipy.linalg.companion(a).T
        B = b[1:] - a[1:]*b[0]

    assuming `a[0]` is 1.0; if `a[0]` is not 1, `a` and `b` are first
    divided by a[0].

    Examples
    --------
    The following code creates a lowpass Butterworth filter. Then it
    applies that filter to an array whose values are all 1.0; the
    output is also all 1.0, as expected for a lowpass filter.  If the
    `zi` argument of `lfilter` had not been given, the output would have
    shown the transient signal.

    >>> from numpy import array, ones
    >>> from scipy.signal import lfilter, lfilter_zi, butter
    >>> b, a = butter(5, 0.25)
    >>> zi = lfilter_zi(b, a)
    >>> y, zo = lfilter(b, a, ones(10), zi=zi)
    >>> y
    array([1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.])

    Another example:

    >>> x = array([0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0])
    >>> y, zf = lfilter(b, a, x, zi=zi*x[0])
    >>> y
    array([ 0.5       ,  0.5       ,  0.5       ,  0.49836039,  0.48610528,
        0.44399389,  0.35505241])

    Note that the `zi` argument to `lfilter` was computed using
    `lfilter_zi` and scaled by `x[0]`.  Then the output `y` has no
    transient until the input drops from 0.5 to 0.0.

    """

    # FIXME: Can this function be replaced with an appropriate
    # use of lfiltic?  For example, when b,a = butter(N,Wn),
    #    lfiltic(b, a, y=numpy.ones_like(a), x=numpy.ones_like(b)).
    #

    # We could use scipy.signal.normalize, but it uses warnings in
    # cases where a ValueError is more appropriate, and it allows
    # b to be 2D.
    b = np.atleast_1d(b)
    if b.ndim != 1:
        raise ValueError("Numerator b must be rank 1.")
    a = np.atleast_1d(a)
    if a.ndim != 1:
        raise ValueError("Denominator a must be rank 1.")

    while len(a) > 1 and a[0] == 0.0:
        a = a[1:]
    if a.size < 1:
        raise ValueError("There must be at least one nonzero `a` coefficient.")

    if a[0] != 1.0:
        # Normalize the coefficients so a[0] == 1.
        a = a / a[0]
        b = b / a[0]

    n = max(len(a), len(b))

    # Pad a or b with zeros so they are the same length.
    if len(a) < n:
        a = np.r_[a, np.zeros(n - len(a))]
    elif len(b) < n:
        b = np.r_[b, np.zeros(n - len(b))]

    IminusA = np.eye(n - 1) - linalg.companion(a).T
    B = b[1:] - a[1:] * b[0]
    # Solve zi = A*zi + B
    zi = np.linalg.solve(IminusA, B)

    # For future reference: we could also use the following
    # explicit formulas to solve the linear system:
    #
    # zi = np.zeros(n - 1)
    # zi[0] = B.sum() / IminusA[:,0].sum()
    # asum = 1.0
    # csum = 0.0
    # for k in range(1,n-1):
    #     asum += a[k]
    #     csum += b[k] - a[k]*b[0]
    #     zi[k] = asum*zi[0] - csum

    return zi


def filtfilt(b, a, x, axis=-1, padtype='odd', padlen=None):
    """
    A forward-backward filter.

    This function applies a linear filter twice, once forward
    and once backwards.  The combined filter has linear phase.

    Before applying the filter, the function can pad the data along the
    given axis in one of three ways: odd, even or constant.  The odd
    and even extensions have the corresponding symmetry about the end point
    of the data.  The constant extension extends the data with the values
    at end points.  On both the forward and backwards passes, the
    initial condition of the filter is found by using `lfilter_zi` and
    scaling it by the end point of the extended data.

    Parameters
    ----------
    b : (N,) array_like
        The numerator coefficient vector of the filter.
    a : (N,) array_like
        The denominator coefficient vector of the filter.  If a[0]
        is not 1, then both a and b are normalized by a[0].
    x : array_like
        The array of data to be filtered.
    axis : int, optional
        The axis of `x` to which the filter is applied.
        Default is -1.
    padtype : str or None, optional
        Must be 'odd', 'even', 'constant', or None.  This determines the
        type of extension to use for the padded signal to which the filter
        is applied.  If `padtype` is None, no padding is used.  The default
        is 'odd'.
    padlen : int or None, optional
        The number of elements by which to extend `x` at both ends of
        `axis` before applying the filter. This value must be less than
        `x.shape[axis]-1`.  `padlen=0` implies no padding.
        The default value is 3*max(len(a),len(b)).

    Returns
    -------
    y : ndarray
        The filtered output, an array of type numpy.float64 with the same
        shape as `x`.

    See Also
    --------
    lfilter_zi, lfilter

    Examples
    --------
    First we create a one second signal that is the sum of two pure sine
    waves, with frequencies 5 Hz and 250 Hz, sampled at 2000 Hz.

    >>> t = np.linspace(0, 1.0, 2001)
    >>> xlow = np.sin(2 * np.pi * 5 * t)
    >>> xhigh = np.sin(2 * np.pi * 250 * t)
    >>> x = xlow + xhigh

    Now create a lowpass Butterworth filter with a cutoff of 0.125 times
    the Nyquist rate, or 125 Hz, and apply it to x with filtfilt.  The
    result should be approximately xlow, with no phase shift.

    >>> from scipy import signal
    >>> b, a = signal.butter(8, 0.125)
    >>> y = signal.filtfilt(b, a, x, padlen=150)
    >>> np.abs(y - xlow).max()
    9.1086182074789912e-06

    We get a fairly clean result for this artificial example because
    the odd extension is exact, and with the moderately long padding,
    the filter's transients have dissipated by the time the actual data
    is reached.  In general, transient effects at the edges are
    unavoidable.

    """

    if padtype not in ['even', 'odd', 'constant', None]:
        raise ValueError(("Unknown value '%s' given to padtype.  padtype must "
                         "be 'even', 'odd', 'constant', or None.") %
                            padtype)

    b = np.asarray(b)
    a = np.asarray(a)
    x = np.asarray(x)

    ntaps = max(len(a), len(b))

    if padtype is None:
        padlen = 0

    if padlen is None:
        # Original padding; preserved for backwards compatibility.
        edge = ntaps * 3
    else:
        edge = padlen

    # x's 'axis' dimension must be bigger than edge.
    if x.shape[axis] <= edge:
        raise ValueError("The length of the input vector x must be at least "
                         "padlen, which is %d." % edge)

    if padtype is not None and edge > 0:
        # Make an extension of length `edge` at each
        # end of the input array.
        if padtype == 'even':
            ext = even_ext(x, edge, axis=axis)
        elif padtype == 'odd':
            ext = odd_ext(x, edge, axis=axis)
        else:
            ext = const_ext(x, edge, axis=axis)
    else:
        ext = x

    # Get the steady state of the filter's step response.
    zi = lfilter_zi(b, a)

    # Reshape zi and create x0 so that zi*x0 broadcasts
    # to the correct value for the 'zi' keyword argument
    # to lfilter.
    zi_shape = [1] * x.ndim
    zi_shape[axis] = zi.size
    zi = np.reshape(zi, zi_shape)
    x0 = axis_slice(ext, stop=1, axis=axis)

    # Forward filter.
    (y, zf) = lfilter(b, a, ext, axis=axis, zi=zi * x0)

    # Backward filter.
    # Create y0 so zi*y0 broadcasts appropriately.
    y0 = axis_slice(y, start=-1, axis=axis)
    (y, zf) = lfilter(b, a, axis_reverse(y, axis=axis), axis=axis, zi=zi * y0)

    # Reverse y.
    y = axis_reverse(y, axis=axis)

    if edge > 0:
        # Slice the actual signal from the extended signal.
        y = axis_slice(y, start=edge, stop=-edge, axis=axis)

    return y


from scipy.signal.filter_design import cheby1
from scipy.signal.fir_filter_design import firwin


def decimate(x, q, n=None, ftype='iir', axis=-1):
    """
    Downsample the signal by using a filter.

    By default, an order 8 Chebyshev type I filter is used.  A 30 point FIR
    filter with hamming window is used if `ftype` is 'fir'.

    Parameters
    ----------
    x : ndarray
        The signal to be downsampled, as an N-dimensional array.
    q : int
        The downsampling factor.
    n : int, optional
        The order of the filter (1 less than the length for 'fir').
    ftype : str {'iir', 'fir'}, optional
        The type of the lowpass filter.
    axis : int, optional
        The axis along which to decimate.

    Returns
    -------
    y : ndarray
        The down-sampled signal.

    See also
    --------
    resample

    """

    if not isinstance(q, int):
        raise TypeError("q must be an integer")

    if n is None:
        if ftype == 'fir':
            n = 30
        else:
            n = 8

    if ftype == 'fir':
        b = firwin(n + 1, 1. / q, window='hamming')
        a = 1.
    else:
        b, a = cheby1(n, 0.05, 0.8 / q)

    y = lfilter(b, a, x, axis=axis)

    sl = [slice(None)] * y.ndim
    sl[axis] = slice(None, None, q)
    return y[sl]
