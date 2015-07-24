# Author: Travis Oliphant
# 1999 -- 2002

from __future__ import division, print_function, absolute_import

import warnings
import threading

from . import sigtools
from scipy._lib.six import callable
from scipy._lib._version import NumpyVersion
from scipy import linalg
from scipy.fftpack import (fft, ifft, ifftshift, fft2, ifft2, fftn,
                           ifftn, fftfreq)
from numpy.fft import rfftn, irfftn
from numpy import (allclose, angle, arange, argsort, array, asarray,
                   atleast_1d, atleast_2d, cast, dot, exp, expand_dims,
                   iscomplexobj, mean, ndarray, newaxis, ones, pi,
                   poly, polyadd, polyder, polydiv, polymul, polysub, polyval,
                   prod, product, r_, ravel, real_if_close, reshape,
                   roots, sort, sum, take, transpose, unique, where, zeros,
                   zeros_like)
import numpy as np
from scipy.special import factorial
from .windows import get_window
from ._arraytools import axis_slice, axis_reverse, odd_ext, even_ext, const_ext


__all__ = ['correlate', 'fftconvolve', 'convolve', 'convolve2d', 'correlate2d',
           'order_filter', 'medfilt', 'medfilt2d', 'wiener', 'lfilter',
           'lfiltic', 'sosfilt', 'deconvolve', 'hilbert', 'hilbert2',
           'cmplx_sort', 'unique_roots', 'invres', 'invresz', 'residue',
           'residuez', 'resample', 'detrend', 'lfilter_zi', 'sosfilt_zi',
           'filtfilt', 'decimate', 'vectorstrength']


_modedict = {'valid': 0, 'same': 1, 'full': 2}

_boundarydict = {'fill': 0, 'pad': 0, 'wrap': 2, 'circular': 2, 'symm': 1,
                 'symmetric': 1, 'reflect': 4}


_rfft_mt_safe = (NumpyVersion(np.__version__) >= '1.9.0.dev-e24486e')

_rfft_lock = threading.Lock()


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
                             " (or 'circular'), \n  and 'symm'"
                             " (or 'symmetric').")
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
    The correlation z of two d-dimensional arrays x and y is defined as:

      z[...,k,...] = sum[..., i_l, ...]
                         x[..., i_l,...] * conj(y[..., i_l + k,...])

    Examples
    --------
    Implement a matched filter using cross-correlation, to recover a signal
    that has passed through a noisy channel.

    >>> from scipy import signal
    >>> sig = np.repeat([0., 1., 1., 0., 1., 0., 0., 1.], 128)
    >>> sig_noise = sig + np.random.randn(len(sig))
    >>> corr = signal.correlate(sig_noise, np.ones(128), mode='same') / 128

    >>> import matplotlib.pyplot as plt
    >>> clock = np.arange(64, len(sig), 128)
    >>> fig, (ax_orig, ax_noise, ax_corr) = plt.subplots(3, 1, sharex=True)
    >>> ax_orig.plot(sig)
    >>> ax_orig.plot(clock, sig[clock], 'ro')
    >>> ax_orig.set_title('Original signal')
    >>> ax_noise.plot(sig_noise)
    >>> ax_noise.set_title('Signal with noise')
    >>> ax_corr.plot(corr)
    >>> ax_corr.plot(clock, corr[clock], 'ro')
    >>> ax_corr.axhline(0.5, ls=':')
    >>> ax_corr.set_title('Cross-correlated with rectangular pulse')
    >>> ax_orig.margins(0, 0.1)
    >>> fig.tight_layout()
    >>> fig.show()

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    # Don't use _valfrommode, since correlate should not accept numeric modes
    try:
        val = _modedict[mode]
    except KeyError:
        raise ValueError("Acceptable mode flags are 'valid',"
                         " 'same', or 'full'.")

    if in1.ndim == in2.ndim == 0:
        return in1 * in2
    elif not in1.ndim == in2.ndim:
        raise ValueError("in1 and in2 should have the same dimensionality")

    if mode == 'valid':
        _check_valid_mode_shapes(in1.shape, in2.shape)
        ps = [i - j + 1 for i, j in zip(in1.shape, in2.shape)]
        out = np.empty(ps, in1.dtype)

        z = sigtools._correlateND(in1, in2, out, val)
    else:
        # _correlateND is far slower when in2.size > in1.size, so swap them
        # and then undo the effect afterward
        swapped_inputs = (mode == 'full') and (in2.size > in1.size)
        if swapped_inputs:
            in1, in2 = in2, in1

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

        # Reverse and conjugate to undo the effect of swapping inputs
        if swapped_inputs:
            slice_obj = [slice(None, None, -1)] * len(z.shape)
            z = z[slice_obj].conj()

    return z


def _centered(arr, newsize):
    # Return the center newsize portion of the array.
    newsize = asarray(newsize)
    currsize = array(arr.shape)
    startind = (currsize - newsize) // 2
    endind = startind + newsize
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]


def _next_regular(target):
    """
    Find the next regular number greater than or equal to target.
    Regular numbers are composites of the prime factors 2, 3, and 5.
    Also known as 5-smooth numbers or Hamming numbers, these are the optimal
    size for inputs to FFTPACK.

    Target must be a positive integer.
    """
    if target <= 6:
        return target

    # Quickly check if it's already a power of 2
    if not (target & (target-1)):
        return target

    match = float('inf')  # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)

            # Quickly find next power of 2 >= quotient
            try:
                p2 = 2**((quotient - 1).bit_length())
            except AttributeError:
                # Fallback for Python <2.7
                p2 = 2**(len(bin(quotient - 1)) - 2)

            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match


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

    Examples
    --------
    Autocorrelation of white noise is an impulse.  (This is at least 100 times
    as fast as `convolve`.)

    >>> from scipy import signal
    >>> sig = np.random.randn(1000)
    >>> autocorr = signal.fftconvolve(sig, sig[::-1], mode='full')

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_mag) = plt.subplots(2, 1)
    >>> ax_orig.plot(sig)
    >>> ax_orig.set_title('White noise')
    >>> ax_mag.plot(np.arange(-len(sig)+1,len(sig)), autocorr)
    >>> ax_mag.set_title('Autocorrelation')
    >>> fig.tight_layout()
    >>> fig.show()

    Gaussian blur implemented using FFT convolution.  Notice the dark borders
    around the image, due to the zero-padding beyond its boundaries.
    The `convolve2d` function allows for other types of image boundaries,
    but is far slower.

    >>> from scipy import misc
    >>> lena = misc.lena()
    >>> kernel = np.outer(signal.gaussian(70, 8), signal.gaussian(70, 8))
    >>> blurred = signal.fftconvolve(lena, kernel, mode='same')

    >>> fig, (ax_orig, ax_kernel, ax_blurred) = plt.subplots(1, 3)
    >>> ax_orig.imshow(lena, cmap='gray')
    >>> ax_orig.set_title('Original')
    >>> ax_orig.set_axis_off()
    >>> ax_kernel.imshow(kernel, cmap='gray')
    >>> ax_kernel.set_title('Gaussian kernel')
    >>> ax_kernel.set_axis_off()
    >>> ax_blurred.imshow(blurred, cmap='gray')
    >>> ax_blurred.set_title('Blurred')
    >>> ax_blurred.set_axis_off()
    >>> fig.show()

    """
    in1 = asarray(in1)
    in2 = asarray(in2)

    if in1.ndim == in2.ndim == 0:  # scalar inputs
        return in1 * in2
    elif not in1.ndim == in2.ndim:
        raise ValueError("in1 and in2 should have the same dimensionality")
    elif in1.size == 0 or in2.size == 0:  # empty arrays
        return array([])

    s1 = array(in1.shape)
    s2 = array(in2.shape)
    complex_result = (np.issubdtype(in1.dtype, complex) or
                      np.issubdtype(in2.dtype, complex))
    shape = s1 + s2 - 1

    if mode == "valid":
        _check_valid_mode_shapes(s1, s2)

    # Speed up FFT by padding to optimal size for FFTPACK
    fshape = [_next_regular(int(d)) for d in shape]
    fslice = tuple([slice(0, int(sz)) for sz in shape])
    # Pre-1.9 NumPy FFT routines are not threadsafe.  For older NumPys, make
    # sure we only call rfftn/irfftn from one thread at a time.
    if not complex_result and (_rfft_mt_safe or _rfft_lock.acquire(False)):
        try:
            ret = irfftn(rfftn(in1, fshape) *
                         rfftn(in2, fshape), fshape)[fslice].copy()
        finally:
            if not _rfft_mt_safe:
                _rfft_lock.release()
    else:
        # If we're here, it's either because we need a complex result, or we
        # failed to acquire _rfft_lock (meaning rfftn isn't threadsafe and
        # is already in use by another thread).  In either case, use the
        # (threadsafe but slower) SciPy complex-FFT routines instead.
        ret = ifftn(fftn(in1, fshape) * fftn(in2, fshape))[fslice].copy()
        if not complex_result:
            ret = ret.real

    if mode == "full":
        return ret
    elif mode == "same":
        return _centered(ret, s1)
    elif mode == "valid":
        return _centered(ret, s1 - s2 + 1)
    else:
        raise ValueError("Acceptable mode flags are 'valid',"
                         " 'same', or 'full'.")


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

    See also
    --------
    numpy.polymul : performs polynomial multiplication (same operation, but
                    also accepts poly1d objects)

    Examples
    --------
    Smooth a square pulse using a Hann window:

    >>> from scipy import signal
    >>> sig = np.repeat([0., 1., 0.], 100)
    >>> win = signal.hann(50)
    >>> filtered = signal.convolve(sig, win, mode='same') / sum(win)

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_win, ax_filt) = plt.subplots(3, 1, sharex=True)
    >>> ax_orig.plot(sig)
    >>> ax_orig.set_title('Original pulse')
    >>> ax_orig.margins(0, 0.1)
    >>> ax_win.plot(win)
    >>> ax_win.set_title('Filter impulse response')
    >>> ax_win.margins(0, 0.1)
    >>> ax_filt.plot(filtered)
    >>> ax_filt.set_title('Filtered signal')
    >>> ax_filt.margins(0, 0.1)
    >>> fig.tight_layout()
    >>> fig.show()

    """
    volume = asarray(in1)
    kernel = asarray(in2)

    if volume.ndim == kernel.ndim == 0:
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

    Examples
    --------
    Compute the gradient of an image by 2D convolution with a complex Scharr
    operator.  (Horizontal operator is real, vertical is imaginary.)  Use
    symmetric boundary condition to avoid creating edges at the image
    boundaries.

    >>> from scipy import signal
    >>> from scipy import misc
    >>> lena = misc.lena()
    >>> scharr = np.array([[ -3-3j, 0-10j,  +3 -3j],
    ...                    [-10+0j, 0+ 0j, +10 +0j],
    ...                    [ -3+3j, 0+10j,  +3 +3j]]) # Gx + j*Gy
    >>> grad = signal.convolve2d(lena, scharr, boundary='symm', mode='same')

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_mag, ax_ang) = plt.subplots(1, 3)
    >>> ax_orig.imshow(lena, cmap='gray')
    >>> ax_orig.set_title('Original')
    >>> ax_orig.set_axis_off()
    >>> ax_mag.imshow(np.absolute(grad), cmap='gray')
    >>> ax_mag.set_title('Gradient magnitude')
    >>> ax_mag.set_axis_off()
    >>> ax_ang.imshow(np.angle(grad), cmap='hsv') # hsv is cyclic, like angles
    >>> ax_ang.set_title('Gradient orientation')
    >>> ax_ang.set_axis_off()
    >>> fig.show()

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

    Examples
    --------
    Use 2D cross-correlation to find the location of a template in a noisy
    image:

    >>> from scipy import signal
    >>> from scipy import misc
    >>> lena = misc.lena() - misc.lena().mean()
    >>> template = np.copy(lena[235:295, 310:370]) # right eye
    >>> template -= template.mean()
    >>> lena = lena + np.random.randn(*lena.shape) * 50 # add noise
    >>> corr = signal.correlate2d(lena, template, boundary='symm', mode='same')
    >>> y, x = np.unravel_index(np.argmax(corr), corr.shape) # find the match

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_orig, ax_template, ax_corr) = plt.subplots(1, 3)
    >>> ax_orig.imshow(lena, cmap='gray')
    >>> ax_orig.set_title('Original')
    >>> ax_orig.set_axis_off()
    >>> ax_template.imshow(template, cmap='gray')
    >>> ax_template.set_title('Template')
    >>> ax_template.set_axis_off()
    >>> ax_corr.imshow(corr, cmap='gray')
    >>> ax_corr.set_title('Cross-correlation')
    >>> ax_corr.set_axis_off()
    >>> ax_orig.plot(x, y, 'ro')
    >>> fig.show()

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
    axis : int, optional
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
    a = np.atleast_1d(a)
    if len(a) == 1:
        # This path only supports types fdgFDGO to mirror _linear_filter below.
        # Any of b, a, x, or zi can set the dtype, but there is no default 
        # casting of other types; instead a NotImplementedError is raised.
        b = np.asarray(b)
        a = np.asarray(a)
        if b.ndim != 1 and a.ndim != 1:
            raise ValueError('object of too small depth for desired array')
        x = np.asarray(x)
        inputs = [b, a, x]
        if zi is not None:
            # _linear_filter does not broadcast zi, but does do expansion of singleton dims.
            zi = np.asarray(zi)
            if zi.ndim != x.ndim:
                raise ValueError('object of too small depth for desired array')
            expected_shape = list(x.shape)
            expected_shape[axis] = b.shape[0] - 1
            expected_shape = tuple(expected_shape)
            # check the trivial case where zi is the right shape first
            if zi.shape != expected_shape:
                strides = zi.ndim * [None]
                if axis < 0:
                    axis += zi.ndim
                for k in range(zi.ndim):
                    if k == axis and zi.shape[k] == expected_shape[k]:
                        strides[k] = zi.strides[k]
                    elif k != axis and zi.shape[k] == expected_shape[k]:
                        strides[k] = zi.strides[k]
                    elif k != axis and zi.shape[k] == 1:
                        strides[k] = 0
                    else:
                        raise ValueError('Unexpected shape for zi: expected '
                                         '%s, found %s.' %
                                         (expected_shape, zi.shape))
                zi = np.lib.stride_tricks.as_strided(zi, expected_shape, strides)
            inputs.append(zi)
        dtype = np.result_type(*inputs)

        if dtype.char not in 'fdgFDGO':
            raise NotImplementedError("input type '%s' not supported" % dtype)

        b = np.array(b, dtype=dtype, copy=False)
        a = np.array(a, dtype=dtype, copy=False)
        b /= a[0]
        x = np.array(x, dtype=dtype, copy=False)

        out_full = np.apply_along_axis(lambda y: np.convolve(b, y), axis, x)
        ind = out_full.ndim * [slice(None)]
        ind[axis] = slice(out_full.shape[axis] - len(b) + 1)
        out = out_full[ind]
        if zi is None:
            return out
        else:
            ind[axis] = slice(out_full.shape[axis] - len(b) + 1, None)
            zf = out_full[ind]
            ind[axis] = slice(zi.shape[axis])
            out[ind] += zi
            return out, zf
    else:
        if zi is None:
            return sigtools._linear_filter(b, a, x, axis)
        else:
            return sigtools._linear_filter(b, a, x, axis, zi)


def lfiltic(b, a, y, x=None):
    """
    Construct initial conditions for lfilter.

    Given a linear filter (b, a) and initial conditions on the output `y`
    and the input `x`, return the initial conditions on the state vector zi
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
    if y.dtype.kind in 'bui':
        # ensure calculations are floating point
        y = y.astype(np.float64)
    zi = zeros(K, y.dtype)
    if x is None:
        x = zeros(M, y.dtype)
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
    """Deconvolves ``divisor`` out of ``signal``.

    Returns the quotient and remainder such that
    ``signal = convolve(divisor, quotient) + remainder``

    Parameters
    ----------
    signal : array_like
        Signal data, typically a recorded signal
    divisor : array_like
        Divisor data, typically an impulse response or filter that was
        applied to the original signal

    Returns
    -------
    quotient : ndarray
        Quotient, typically the recovered original signal
    remainder : ndarray
        Remainder

    Examples
    --------
    Deconvolve a signal that's been filtered:

    >>> from scipy import signal
    >>> original = [0, 1, 0, 0, 1, 1, 0, 0]
    >>> impulse_response = [2, 1]
    >>> recorded = signal.convolve(impulse_response, original)
    >>> recorded
    array([0, 2, 1, 0, 2, 3, 1, 0, 0])
    >>> recovered, remainder = signal.deconvolve(recorded, impulse_response)
    >>> recovered
    array([ 0.,  1.,  0.,  0.,  1.,  1.,  0.,  0.])

    See also
    --------
    numpy.polydiv : performs polynomial division (same operation, but
                    also accepts poly1d objects)

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

    Examples
    ---------
    In this example we use the Hilbert transform to determine the amplitude
    envelope and instantaneous frequency of an amplitude-modulated signal.
        
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import hilbert, chirp

    >>> duration = 1.0
    >>> fs = 400.0
    >>> samples = int(fs*duration)
    >>> t = np.arange(samples) / fs

    We create a chirp of which the frequency increases from 20 Hz to 100 Hz and 
    apply an amplitude modulation.
    
    >>> signal = chirp(t, 20.0, t[-1], 100.0)    
    >>> signal *= (1.0 + 0.5 * np.sin(2.0*np.pi*3.0*t) )

    The amplitude envelope is given by magnitude of the analytic signal. The 
    instantaneous frequency can be obtained by differentiating the instantaneous 
    phase in respect to time. The instantaneous phase corresponds to the phase 
    angle of the analytic signal.

    >>> analytic_signal = hilbert(signal)
    >>> amplitude_envelope = np.abs(analytic_signal)
    >>> instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    >>> instantaneous_frequency = np.diff(instantaneous_phase) / (2.0*np.pi) * fs

    >>> fig = plt.figure()
    >>> ax0 = fig.add_subplot(211)
    >>> ax0.plot(t, signal, label='signal')
    >>> ax0.plot(t, amplitude_envelope, label='envelope')
    >>> ax0.set_xlabel("time in seconds")
    >>> ax0.legend()
    >>> ax1 = fig.add_subplot(212)
    >>> ax1.plot(t[1:], instantaneous_frequency)
    >>> ax1.set_xlabel("time in seconds")
    >>> ax1.set_ylim(0.0, 120.0)

    References
    ----------
    .. [1] Wikipedia, "Analytic signal".
           http://en.wikipedia.org/wiki/Analytic_signal
    .. [2] Leon Cohen, "Time-Frequency Analysis", 1995. Chapter 2.
    .. [3] Alan V. Oppenheim, Ronald W. Schafer. Discrete-Time Signal Processing, 
           Third Edition, 2009. Chapter 12. ISBN 13: 978-1292-02572-8

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
        raise ValueError("x must be 2-D.")
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
    """Sort roots based on magnitude.

    Parameters
    ----------
    p : array_like
        The roots to sort, as a 1-D array.

    Returns
    -------
    p_sorted : ndarray
        Sorted roots.
    indx : ndarray
        Array of indices needed to sort the input `p`.

    """
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
    Compute b(s) and a(s) from partial fraction expansion.

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
            b = polyadd(b, r[indx] * atleast_1d(poly(t2)))
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
            r[indx + m - 1] = (polyval(bn, pout[n]) / polyval(an, pout[n])
                               / factorial(sig - m))
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
    Compute b(z) and a(z) from partial fraction expansion.

    If ``M = len(b)`` and ``N = len(a)``::

                b(z)     b[0] + b[1] z**(-1) + ... + b[M-1] z**(-M+1)
        H(z) = ------ = ----------------------------------------------
                a(z)     a[0] + a[1] z**(-1) + ... + a[N-1] z**(-N+1)

                     r[0]                   r[-1]
             = --------------- + ... + ---------------- + k[0] + k[1]z**(-1)...
               (1-p[0]z**(-1))         (1-p[-1]z**(-1))

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like::

             r[i]              r[i+1]                    r[i+n-1]
        -------------- + ------------------ + ... + ------------------
        (1-p[i]z**(-1))  (1-p[i]z**(-1))**2         (1-p[i]z**(-1))**n

    See Also
    --------
    residuez, unique_roots, invres

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
            brev = polyadd(brev, (r[indx] * atleast_1d(poly(t2)))[::-1])
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
    from ``dx`` to ``dx * len(x) / num``.

    If `t` is not None, then it represents the old sample positions,
    and the new sample positions will be returned as well as the new
    samples.

    As noted, `resample` uses FFT transformations, which can be very
    slow if the number of input or output samples is large and prime;
    see `scipy.fftpack.fft`.

    Examples
    --------
    Note that the end of the resampled data rises to meet the first
    sample of the next cycle:

    >>> from scipy import signal
    
    >>> x = np.linspace(0, 10, 20, endpoint=False)
    >>> y = np.cos(-x**2/6.0)
    >>> f = signal.resample(y, 100)
    >>> xnew = np.linspace(0, 10, 100, endpoint=False)
    
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, y, 'go-', xnew, f, '.-', 10, y[0], 'ro')
    >>> plt.legend(['data', 'resampled'], loc='best')
    >>> plt.show()
    """
    x = asarray(x)
    X = fft(x, axis=axis)
    Nx = x.shape[axis]
    if window is not None:
        if callable(window):
            W = window(fftfreq(Nx))
        elif isinstance(window, ndarray):
            if window.shape != (Nx,):
                raise ValueError('window must have the same length as data')
            W = window
        else:
            W = ifftshift(get_window(window, Nx))
        newshape = [1] * x.ndim
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


def vectorstrength(events, period):
    '''
    Determine the vector strength of the events corresponding to the given
    period.

    The vector strength is a measure of phase synchrony, how well the
    timing of the events is synchronized to a single period of a periodic
    signal.

    If multiple periods are used, calculate the vector strength of each.
    This is called the "resonating vector strength".

    Parameters
    ----------
    events : 1D array_like
        An array of time points containing the timing of the events.
    period : float or array_like
        The period of the signal that the events should synchronize to.
        The period is in the same units as `events`.  It can also be an array
        of periods, in which case the outputs are arrays of the same length.

    Returns
    -------
    strength : float or 1D array
        The strength of the synchronization.  1.0 is perfect synchronization
        and 0.0 is no synchronization.  If `period` is an array, this is also
        an array with each element containing the vector strength at the
        corresponding period.
    phase : float or array
        The phase that the events are most strongly synchronized to in radians.
        If `period` is an array, this is also an array with each element
        containing the phase for the corresponding period.

    References
    ----------
    van Hemmen, JL, Longtin, A, and Vollmayr, AN. Testing resonating vector
        strength: Auditory system, electric fish, and noise.
        Chaos 21, 047508 (2011);
        doi: 10.1063/1.3670512
    van Hemmen, JL.  Vector strength after Goldberg, Brown, and von Mises:
        biological and mathematical perspectives.  Biol Cybern.
        2013 Aug;107(4):385-96. doi: 10.1007/s00422-013-0561-7.
    van Hemmen, JL and Vollmayr, AN.  Resonating vector strength: what happens
        when we vary the "probing" frequency while keeping the spike times
        fixed.  Biol Cybern. 2013 Aug;107(4):491-94.
        doi: 10.1007/s00422-013-0560-8
    '''
    events = asarray(events)
    period = asarray(period)
    if events.ndim > 1:
        raise ValueError('events cannot have dimensions more than 1')
    if period.ndim > 1:
        raise ValueError('period cannot have dimensions more than 1')

    # we need to know later if period was originally a scalar
    scalarperiod = not period.ndim

    events = atleast_2d(events)
    period = atleast_2d(period)
    if (period <= 0).any():
        raise ValueError('periods must be positive')

    # this converts the times to vectors
    vectors = exp(dot(2j*pi/period.T, events))

    # the vector strength is just the magnitude of the mean of the vectors
    # the vector phase is the angle of the mean of the vectors
    vectormean = mean(vectors, axis=1)
    strength = abs(vectormean)
    phase = angle(vectormean)

    # if the original period was a scalar, return scalars
    if scalarperiod:
        strength = strength[0]
        phase = phase[0]
    return strength, phase


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
        raise ValueError("Numerator b must be 1-D.")
    a = np.atleast_1d(a)
    if a.ndim != 1:
        raise ValueError("Denominator a must be 1-D.")

    while len(a) > 1 and a[0] == 0.0:
        a = a[1:]
    if a.size < 1:
        raise ValueError("There must be at least one nonzero `a` coefficient.")

    if a[0] != 1.0:
        # Normalize the coefficients so a[0] == 1.
        b = b / a[0]
        a = a / a[0]

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


def sosfilt_zi(sos):
    """
    Compute an initial state `zi` for the sosfilt function that corresponds
    to the steady state of the step response.

    A typical use of this function is to set the initial state so that the
    output of the filter starts at the same value as the first element of
    the signal to be filtered.

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    Returns
    -------
    zi : ndarray
        Initial conditions suitable for use with ``sosfilt``, shape
        ``(n_sections, 2)``.

    See Also
    --------
    sosfilt, zpk2sos

    Notes
    -----
    .. versionadded:: 0.16.0

    Examples
    --------
    Filter a rectangular pulse that begins at time 0, with and without
    the use of the `zi` argument of `scipy.signal.sosfilt`.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> sos = signal.butter(9, 0.125, output='sos')
    >>> zi = signal.sosfilt_zi(sos)
    >>> x = (np.arange(250) < 100).astype(int)
    >>> f1 = signal.sosfilt(sos, x)
    >>> f2, zo = signal.sosfilt(sos, x, zi=zi)

    >>> plt.plot(x, 'k--', label='x')
    >>> plt.plot(f1, 'b', alpha=0.5, linewidth=2, label='filtered')
    >>> plt.plot(f2, 'g', alpha=0.25, linewidth=4, label='filtered with zi')
    >>> plt.legend(loc='best')
    >>> plt.show()

    """
    sos = np.asarray(sos)
    if sos.ndim != 2 or sos.shape[1] != 6:
        raise ValueError('sos must be shape (n_sections, 6)')

    n_sections = sos.shape[0]
    zi = np.empty((n_sections, 2))
    scale = 1.0
    for section in range(n_sections):
        b = sos[section, :3]
        a = sos[section, 3:]
        zi[section] = scale * lfilter_zi(b, a)
        # If H(z) = B(z)/A(z) is this section's transfer function, then
        # b.sum()/a.sum() is H(1), the gain at omega=0.  That's the steady
        # state value of this section's step response.
        scale *= b.sum() / a.sum()

    return zi


def _filtfilt_gust(b, a, x, axis=-1, irlen=None):
    """Forward-backward IIR filter that uses Gustafsson's method.

    Apply the IIR filter defined by `(b,a)` to `x` twice, first forward
    then backward, using Gustafsson's initial conditions [1]_.

    Let ``y_fb`` be the result of filtering first forward and then backward,
    and let ``y_bf`` be the result of filtering first backward then forward.
    Gustafsson's method is to compute initial conditions for the forward
    pass and the backward pass such that ``y_fb == y_bf``.

    Parameters
    ----------
    b : scalar or 1-D ndarray
        Numerator coefficients of the filter.
    a : scalar or 1-D ndarray
        Denominator coefficients of the filter.
    x : ndarray
        Data to be filtered.
    axis : int, optional
        Axis of `x` to be filtered.  Default is -1.
    irlen : int or None, optional
        The length of the nonnegligible part of the impulse response.
        If `irlen` is None, or if the length of the signal is less than
        ``2 * irlen``, then no part of the impulse response is ignored.

    Returns
    -------
    y : ndarray
        The filtered data.
    x0 : ndarray
        Initial condition for the forward filter.
    x1 : ndarray
        Initial condition for the backward filter.

    Notes
    -----
    Typically the return values `x0` and `x1` are not needed by the
    caller.  The intended use of these return values is in unit tests.

    References
    ----------
    .. [1] F. Gustaffson. Determining the initial states in forward-backward
           filtering. Transactions on Signal Processing, 46(4):988-992, 1996.

    """
    # In the comments, "Gustafsson's paper" and [1] refer to the
    # paper referenced in the docstring.

    b = np.atleast_1d(b)
    a = np.atleast_1d(a)

    order = max(len(b), len(a)) - 1
    if order == 0:
        # The filter is just scalar multiplication, with no state.
        scale = (b[0] / a[0])**2
        y = scale * x
        return y, np.array([]), np.array([])

    if axis != -1 or axis != x.ndim - 1:
        # Move the axis containing the data to the end.
        x = np.swapaxes(x, axis, x.ndim - 1)

    # n is the number of samples in the data to be filtered.
    n = x.shape[-1]

    if irlen is None or n <= 2*irlen:
        m = n
    else:
        m = irlen

    # Create Obs, the observability matrix (called O in the paper).
    # This matrix can be interpreted as the operator that propagates
    # an arbitrary initial state to the output, assuming the input is
    # zero.
    # In Gustafsson's paper, the forward and backward filters are not
    # necessarily the same, so he has both O_f and O_b.  We use the same
    # filter in both directions, so we only need O. The same comment
    # applies to S below.
    Obs = np.zeros((m, order))
    zi = np.zeros(order)
    zi[0] = 1
    Obs[:, 0] = lfilter(b, a, np.zeros(m), zi=zi)[0]
    for k in range(1, order):
        Obs[k:, k] = Obs[:-k, 0]

    # Obsr is O^R (Gustafsson's notation for row-reversed O)
    Obsr = Obs[::-1]

    # Create S.  S is the matrix that applies the filter to the reversed
    # propagated initial conditions.  That is,
    #     out = S.dot(zi)
    # is the same as
    #     tmp, _ = lfilter(b, a, zeros(), zi=zi)  # Propagate ICs.
    #     out = lfilter(b, a, tmp[::-1])          # Reverse and filter.

    # Equations (5) & (6) of [1]
    S = lfilter(b, a, Obs[::-1], axis=0)

    # Sr is S^R (row-reversed S)
    Sr = S[::-1]

    # M is [(S^R - O), (O^R - S)]
    if m == n:
        M = np.hstack((Sr - Obs, Obsr - S))
    else:
        # Matrix described in section IV of [1].
        M = np.zeros((2*m, 2*order))
        M[:m, :order] = Sr - Obs
        M[m:, order:] = Obsr - S

    # Naive forward-backward and backward-forward filters.
    # These have large transients because the filters use zero initial
    # conditions.
    y_f = lfilter(b, a, x)
    y_fb = lfilter(b, a, y_f[..., ::-1])[..., ::-1]

    y_b = lfilter(b, a, x[..., ::-1])[..., ::-1]
    y_bf = lfilter(b, a, y_b)

    delta_y_bf_fb = y_bf - y_fb
    if m == n:
        delta = delta_y_bf_fb
    else:
        start_m = delta_y_bf_fb[..., :m]
        end_m = delta_y_bf_fb[..., -m:]
        delta = np.concatenate((start_m, end_m), axis=-1)

    # ic_opt holds the "optimal" initial conditions.
    # The following code computes the result shown in the formula
    # of the paper between equations (6) and (7).
    if delta.ndim == 1:
        ic_opt = linalg.lstsq(M, delta)[0]
    else:
        # Reshape delta so it can be used as an array of multiple
        # right-hand-sides in linalg.lstsq.
        delta2d = delta.reshape(-1, delta.shape[-1]).T
        ic_opt0 = linalg.lstsq(M, delta2d)[0].T
        ic_opt = ic_opt0.reshape(delta.shape[:-1] + (M.shape[-1],))

    # Now compute the filtered signal using equation (7) of [1].
    # First, form [S^R, O^R] and call it W.
    if m == n:
        W = np.hstack((Sr, Obsr))
    else:
        W = np.zeros((2*m, 2*order))
        W[:m, :order] = Sr
        W[m:, order:] = Obsr

    # Equation (7) of [1] says
    #     Y_fb^opt = Y_fb^0 + W * [x_0^opt; x_{N-1}^opt]
    # `wic` is (almost) the product on the right.
    # W has shape (m, 2*order), and ic_opt has shape (..., 2*order),
    # so we can't use W.dot(ic_opt).  Instead, we dot ic_opt with W.T,
    # so wic has shape (..., m).
    wic = ic_opt.dot(W.T)

    # `wic` is "almost" the product of W and the optimal ICs in equation
    # (7)--if we're using a truncated impulse response (m < n), `wic`
    # contains only the adjustments required for the ends of the signal.
    # Here we form y_opt, taking this into account if necessary.
    y_opt = y_fb
    if m == n:
        y_opt += wic
    else:
        y_opt[..., :m] += wic[..., :m]
        y_opt[..., -m:] += wic[..., -m:]

    x0 = ic_opt[..., :order]
    x1 = ic_opt[..., -order:]
    if axis != -1 or axis != x.ndim - 1:
        # Restore the data axis to its original position.
        x0 = np.swapaxes(x0, axis, x.ndim - 1)
        x1 = np.swapaxes(x1, axis, x.ndim - 1)
        y_opt = np.swapaxes(y_opt, axis, x.ndim - 1)

    return y_opt, x0, x1


def filtfilt(b, a, x, axis=-1, padtype='odd', padlen=None, method='pad',
             irlen=None):
    """
    A forward-backward filter.

    This function applies a linear filter twice, once forward and once
    backwards.  The combined filter has linear phase.

    The function provides options for handling the edges of the signal.

    When `method` is "pad", the function pads the data along the given axis
    in one of three ways: odd, even or constant.  The odd and even extensions
    have the corresponding symmetry about the end point of the data.  The
    constant extension extends the data with the values at the end points. On
    both the forward and backward passes, the initial condition of the
    filter is found by using `lfilter_zi` and scaling it by the end point of
    the extended data.

    When `method` is "gust", Gustafsson's method [1]_ is used.  Initial
    conditions are chosen for the forward and backward passes so that the
    forward-backward filter gives the same result as the backward-forward
    filter.

    Parameters
    ----------
    b : (N,) array_like
        The numerator coefficient vector of the filter.
    a : (N,) array_like
        The denominator coefficient vector of the filter.  If ``a[0]``
        is not 1, then both `a` and `b` are normalized by ``a[0]``.
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
        `axis` before applying the filter.  This value must be less than
        ``x.shape[axis] - 1``.  ``padlen=0`` implies no padding.
        The default value is ``3 * max(len(a), len(b))``.
    method : str, optional
        Determines the method for handling the edges of the signal, either
        "pad" or "gust".  When `method` is "pad", the signal is padded; the
        type of padding is determined by `padtype` and `padlen`, and `irlen`
        is ignored.  When `method` is "gust", Gustafsson's method is used,
        and `padtype` and `padlen` are ignored.
    irlen : int or None, optional
        When `method` is "gust", `irlen` specifies the length of the
        impulse response of the filter.  If `irlen` is None, no part
        of the impulse response is ignored.  For a long signal, specifying
        `irlen` can significantly improve the performance of the filter.

    Returns
    -------
    y : ndarray
        The filtered output, an array of type numpy.float64 with the same
        shape as `x`.

    See Also
    --------
    lfilter_zi, lfilter

    Notes
    -----
    The option to use Gustaffson's method was added in scipy version 0.16.0.

    References
    ----------
    .. [1] F. Gustaffson, "Determining the initial states in forward-backward
           filtering", Transactions on Signal Processing, Vol. 46, pp. 988-992,
           1996.

    Examples
    --------
    The examples will use several functions from `scipy.signal`.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    First we create a one second signal that is the sum of two pure sine
    waves, with frequencies 5 Hz and 250 Hz, sampled at 2000 Hz.

    >>> t = np.linspace(0, 1.0, 2001)
    >>> xlow = np.sin(2 * np.pi * 5 * t)
    >>> xhigh = np.sin(2 * np.pi * 250 * t)
    >>> x = xlow + xhigh

    Now create a lowpass Butterworth filter with a cutoff of 0.125 times
    the Nyquist rate, or 125 Hz, and apply it to ``x`` with `filtfilt`.
    The result should be approximately ``xlow``, with no phase shift.

    >>> b, a = signal.butter(8, 0.125)
    >>> y = signal.filtfilt(b, a, x, padlen=150)
    >>> np.abs(y - xlow).max()
    9.1086182074789912e-06

    We get a fairly clean result for this artificial example because
    the odd extension is exact, and with the moderately long padding,
    the filter's transients have dissipated by the time the actual data
    is reached.  In general, transient effects at the edges are
    unavoidable.

    The following example demonstrates the option ``method="gust"``.

    First, create a filter.

    >>> b, a = signal.ellip(4, 0.01, 120, 0.125)  # Filter to be applied.
    >>> np.random.seed(123456)

    `sig` is a random input signal to be filtered.

    >>> n = 60
    >>> sig = np.random.randn(n)**3 + 3*np.random.randn(n).cumsum()

    Apply `filtfilt` to `sig`, once using the Gustafsson method, and
    once using padding, and plot the results for comparison.

    >>> fgust = signal.filtfilt(b, a, sig, method="gust")
    >>> fpad = signal.filtfilt(b, a, sig, padlen=50)
    >>> plt.plot(sig, 'k-', label='input')
    >>> plt.plot(fgust, 'b-', linewidth=4, label='gust')
    >>> plt.plot(fpad, 'c-', linewidth=1.5, label='pad')
    >>> plt.legend(loc='best')
    >>> plt.show()

    The `irlen` argument can be used to improve the performance
    of Gustafsson's method.

    Estimate the impulse response length of the filter.

    >>> z, p, k = signal.tf2zpk(b, a)
    >>> eps = 1e-9
    >>> r = np.max(np.abs(p))
    >>> approx_impulse_len = int(np.ceil(np.log(eps) / np.log(r)))
    >>> approx_impulse_len
    137

    Apply the filter to a longer signal, with and without the `irlen`
    argument.  The difference between `y1` and `y2` is small.  For long
    signals, using `irlen` gives a significant performance improvement.

    >>> x = np.random.randn(5000)
    >>> y1 = signal.filtfilt(b, a, x, method='gust')
    >>> y2 = signal.filtfilt(b, a, x, method='gust', irlen=approx_impulse_len)
    >>> print(np.max(np.abs(y1 - y2)))
    1.80056858312e-10

    """
    b = np.atleast_1d(b)
    a = np.atleast_1d(a)
    x = np.asarray(x)

    if method not in ["pad", "gust"]:
        raise ValueError("method must be 'pad' or 'gust'.")

    if method == "gust":
        y, z1, z2 = _filtfilt_gust(b, a, x, axis=axis, irlen=irlen)
        return y

    # `method` is "pad"...

    ntaps = max(len(a), len(b))

    if padtype not in ['even', 'odd', 'constant', None]:
        raise ValueError(("Unknown value '%s' given to padtype.  padtype "
                          "must be 'even', 'odd', 'constant', or None.") %
                         padtype)

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


def sosfilt(sos, x, axis=-1, zi=None):
    """
    Filter data along one dimension using cascaded second-order sections

    Filter a data sequence, `x`, using a digital IIR filter defined by
    `sos`. This is implemented by performing `lfilter` for each
    second-order section.  See `lfilter` for details.

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. Each row corresponds to a second-order
        section, with the first three columns providing the numerator
        coefficients and the last three providing the denominator
        coefficients.
    x : array_like
        An N-dimensional input array.
    axis : int, optional
        The axis of the input data array along which to apply the
        linear filter. The filter is applied to each subarray along
        this axis.  Default is -1.
    zi : array_like, optional
        Initial conditions for the cascaded filter delays.  It is a (at
        least 2D) vector of shape ``(n_sections, ..., 2, ...)``, where
        ``..., 2, ...`` denotes the shape of `x`, but with ``x.shape[axis]``
        replaced by 2.  If `zi` is None or is not given then initial rest
        (i.e. all zeros) is assumed.
        Note that these initial conditions are *not* the same as the initial
        conditions given by `lfiltic` or `lfilter_zi`.

    Returns
    -------
    y : ndarray
        The output of the digital filter.
    zf : ndarray, optional
        If `zi` is None, this is not returned, otherwise, `zf` holds the
        final filter delay values.

    See Also
    --------
    zpk2sos, sos2zpk, sosfilt_zi

    Notes
    -----
    The filter function is implemented as a series of second-order filters
    with direct-form II transposed structure. It is designed to minimize
    numerical precision errors for high-order filters.

    .. versionadded:: 0.16.0

    Examples
    --------
    Plot a 13th-order filter's impulse response using both `lfilter` and
    `sosfilt`, showing the instability that results from trying to do a
    13th-order filter in a single stage (the numerical error pushes some poles
    outside of the unit circle):

    >>> import matplotlib.pyplot as plt
    >>> from scipy import signal
    >>> b, a = signal.ellip(13, 0.009, 80, 0.05, output='ba')
    >>> sos = signal.ellip(13, 0.009, 80, 0.05, output='sos')
    >>> x = np.zeros(700)
    >>> x[0] = 1.
    >>> y_tf = signal.lfilter(b, a, x)
    >>> y_sos = signal.sosfilt(sos, x)
    >>> plt.plot(y_tf, 'r', label='TF')
    >>> plt.plot(y_sos, 'k', label='SOS')
    >>> plt.legend(loc='best')
    >>> plt.show()

    """
    x = np.asarray(x)

    sos = atleast_2d(sos)
    if sos.ndim != 2:
        raise ValueError('sos array must be 2D')

    n_sections, m = sos.shape
    if m != 6:
        raise ValueError('sos array must be shape (n_sections, 6)')

    use_zi = zi is not None
    if use_zi:
        zi = np.asarray(zi)
        x_zi_shape = list(x.shape)
        x_zi_shape[axis] = 2
        x_zi_shape = tuple([n_sections] + x_zi_shape)
        if zi.shape != x_zi_shape:
            raise ValueError('Invalid zi shape.  With axis=%r, an input with '
                             'shape %r, and an sos array with %d sections, zi '
                             'must have shape %r.' %
                             (axis, x.shape, n_sections, x_zi_shape))
        zf = zeros_like(zi)

    for section in range(n_sections):
        if use_zi:
            x, zf[section] = lfilter(sos[section, :3], sos[section, 3:],
                                     x, axis, zi=zi[section])
        else:
            x = lfilter(sos[section, :3], sos[section, 3:], x, axis)
    out = (x, zf) if use_zi else x
    return out


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
