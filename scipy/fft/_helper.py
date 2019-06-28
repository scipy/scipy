import scipy.fftpack as _fftpack


def next_fast_len(target, dtype=None):
    """Find the next fast size of input data to ``fft``, for zero-padding, etc.

    SciPy's FFT algorithms gain their speed by a recursive divide and conquer
    strategy. This relies on efficient functions for small prime factors of the
    input length. Thus, the transforms are fastest when using composites of the
    prime factors handled by the fft implementation. If there are efficient
    functions for all radices <= `n` then the result is the smallest number `x`
    >= ``target`` with only prime factors < `n`. (Also known as `n`-smooth
    numbers)

    Parameters
    ----------
    target : int
        Length to start searching from.  Must be a positive integer.

    Returns
    -------
    out : int
        The smallest fast length greater than or equal to ``target``.

    Notes
    -----
    The result of this function may change if new prime factors are added.

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

    # TODO: Complex transforms are also fast for 11-smooth lengths
    return _fftpack.next_fast_len(target)
