from . import _pocketfft
from ._pocketfft import helper as _helper
from bisect import bisect_left
import numpy as np


def next_fast_len(target, kind='C2C'):
    """Find the next fast size of input data to ``fft``, for zero-padding, etc.

    SciPy's FFT algorithms gain their speed by a recursive divide and conquer
    strategy. This relies on efficient functions for small prime factors of the
    input length. Thus, the transforms are fastest when using composites of the
    prime factors handled by the fft implementation. If there are efficient
    functions for all radices <= `n` then the result will be a number `x`
    >= ``target`` with only prime factors < `n`. (Also known as `n`-smooth
    numbers)

    Parameters
    ----------
    target : int
        Length to start searching from.  Must be a positive integer.
    kind : {'C2C', 'C2R', 'R2C'}, optional
        Kind of transform to be performed
        - 'C2C': Complex to complex (e.g. `fft`, `ifft`)
        - 'C2R': Complex to real (e.g. `irfft`, `hfft`)
        - 'R2C': Real to complex (e.g. `rfft`, `ihfft`)

    Returns
    -------
    out : int
        The smallest fast length greater than or equal to ``target``.

    Notes
    -----
    The result of this function may change in future as performance
    considerations change, for example if new prime factors are added.

    Calling `fft` or `ifft` with real input data performs an ``'R2C'``
    transform internally.

    Examples
    --------
    On a particular machine, an FFT of prime length takes 17 ms:

    >>> from scipy import fft
    >>> min_len = 93059  # prime length is worst case for speed
    >>> a = np.random.randn(min_len)
    >>> b = fft.fft(a)

    Zero-padding to the next regular length reduces computation time to
    1.3 ms, a speedup of 13 times:

    >>> fft.next_fast_len(min_len)
    93312
    >>> b = fft.fft(a, 93312)

    Rounding up to the next power of 2 is not optimal, taking 1.9 ms to
    compute; 1.3 times longer than the size given by ``next_fast_len``:

    >>> b = fft.fft(a, 131072)

    """
    return _helper.next_fast_len(target, kind)


def _init_nd_shape_and_axes(x, shape, axes):
    """Handle shape and axes arguments for n-dimensional transforms.

    Returns the shape and axes in a standard form, taking into account negative
    values and checking for various potential errors.

    Parameters
    ----------
    x : array_like
        The input array.
    shape : int or array_like of ints or None
        The shape of the result.  If both `shape` and `axes` (see below) are
        None, `shape` is ``x.shape``; if `shape` is None but `axes` is
        not None, then `shape` is ``scipy.take(x.shape, axes, axis=0)``.
        If `shape` is -1, the size of the corresponding dimension of `x` is
        used.
    axes : int or array_like of ints or None
        Axes along which the calculation is computed.
        The default is over all axes.
        Negative indices are automatically converted to their positive
        counterpart.

    Returns
    -------
    shape : array
        The shape of the result. It is a 1D integer array.
    axes : array
        The shape of the result. It is a 1D integer array.

    """
    return _helper._init_nd_shape_and_axes(x, shape, axes)
