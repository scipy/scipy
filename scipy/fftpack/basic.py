"""
Discrete Fourier Transforms - basic.py
"""
# Created by Pearu Peterson, August,September 2002

__all__ = ['fft','ifft','fftn','ifftn','rfft','irfft',
           'fft2','ifft2', 'rfftfreq']

from numpy import zeros, swapaxes, integer, array
import numpy
import _fftpack

import atexit
atexit.register(_fftpack.destroy_zfft_cache)
atexit.register(_fftpack.destroy_zfftnd_cache)
atexit.register(_fftpack.destroy_drfft_cache)
atexit.register(_fftpack.destroy_cfft_cache)
atexit.register(_fftpack.destroy_cfftnd_cache)
atexit.register(_fftpack.destroy_rfft_cache)
del atexit

def istype(arr, typeclass):
    return issubclass(arr.dtype.type, typeclass)

# XXX: single precision FFTs partially disabled due to accuracy issues
#      for large prime-sized inputs.
#
#      See http://permalink.gmane.org/gmane.comp.python.scientific.devel/13834
#      ("fftpack test failures for 0.8.0b1", Ralf Gommers, 17 Jun 2010,
#       @ scipy-dev)
#
#      These should be re-enabled once the problems are resolved

def _is_safe_size(n):
    """
    Is the size of FFT such that FFTPACK can handle it in single precision
    with sufficient accuracy?

    Composite numbers of 2, 3, and 5 are accepted, as FFTPACK has those
    """
    n = int(n)
    for c in (2, 3, 5):
        while n % c == 0:
            n /= c
    return (n <= 1)

def _fake_crfft(x, n, *a, **kw):
    if _is_safe_size(n):
        return _fftpack.crfft(x, n, *a, **kw)
    else:
        return _fftpack.zrfft(x, n, *a, **kw).astype(numpy.complex64)

def _fake_cfft(x, n, *a, **kw):
    if _is_safe_size(n):
        return _fftpack.cfft(x, n, *a, **kw)
    else:
        return _fftpack.zfft(x, n, *a, **kw).astype(numpy.complex64)

def _fake_rfft(x, n, *a, **kw):
    if _is_safe_size(n):
        return _fftpack.rfft(x, n, *a, **kw)
    else:
        return _fftpack.drfft(x, n, *a, **kw).astype(numpy.float32)

def _fake_cfftnd(x, shape, *a, **kw):
    if numpy.all(map(_is_safe_size, shape)):
        return _fftpack.cfftnd(x, shape, *a, **kw)
    else:
        return _fftpack.zfftnd(x, shape, *a, **kw).astype(numpy.complex64)

_DTYPE_TO_FFT = {
#        numpy.dtype(numpy.float32): _fftpack.crfft,
        numpy.dtype(numpy.float32): _fake_crfft,
        numpy.dtype(numpy.float64): _fftpack.zrfft,
#        numpy.dtype(numpy.complex64): _fftpack.cfft,
        numpy.dtype(numpy.complex64): _fake_cfft,
        numpy.dtype(numpy.complex128): _fftpack.zfft,
}

_DTYPE_TO_RFFT = {
#        numpy.dtype(numpy.float32): _fftpack.rfft,
        numpy.dtype(numpy.float32): _fake_rfft,
        numpy.dtype(numpy.float64): _fftpack.drfft,
}

_DTYPE_TO_FFTN = {
#        numpy.dtype(numpy.complex64): _fftpack.cfftnd,
        numpy.dtype(numpy.complex64): _fake_cfftnd,
        numpy.dtype(numpy.complex128): _fftpack.zfftnd,
#        numpy.dtype(numpy.float32): _fftpack.cfftnd,
        numpy.dtype(numpy.float32): _fake_cfftnd,
        numpy.dtype(numpy.float64): _fftpack.zfftnd,
}

def _asfarray(x):
    """Like numpy asfarray, except that it does not modify x dtype if x is
    already an array with a float dtype, and do not cast complex types to
    real."""
    if hasattr(x, "dtype") and x.dtype.char in numpy.typecodes["AllFloat"]:
        return x
    else:
        # We cannot use asfarray directly because it converts sequences of
        # complex to sequence of real
        ret = numpy.asarray(x)
        if not ret.dtype.char in numpy.typecodes["AllFloat"]:
            return numpy.asfarray(x)
        return ret

def _fix_shape(x, n, axis):
    """ Internal auxiliary function for _raw_fft, _raw_fftnd."""
    s = list(x.shape)
    if s[axis] > n:
        index = [slice(None)]*len(s)
        index[axis] = slice(0,n)
        x = x[index]
    else:
        index = [slice(None)]*len(s)
        index[axis] = slice(0,s[axis])
        s[axis] = n
        z = zeros(s,x.dtype.char)
        z[index] = x
        x = z
    return x


def _raw_fft(x, n, axis, direction, overwrite_x, work_function):
    """ Internal auxiliary function for fft, ifft, rfft, irfft."""
    if n is None:
        n = x.shape[axis]
    elif n != x.shape[axis]:
        x = _fix_shape(x,n,axis)
        overwrite_x = 1
    if axis == -1 or axis == len(x.shape)-1:
        r = work_function(x,n,direction,overwrite_x=overwrite_x)
    else:
        x = swapaxes(x, axis, -1)
        r = work_function(x,n,direction,overwrite_x=overwrite_x)
        r = swapaxes(r, axis, -1)
    return r


def fft(x, n=None, axis=-1, overwrite_x=0):
    """
    Return discrete Fourier transform of arbitrary type sequence x.

    Parameters
    ----------
    x : array-like
        array to fourier transform.
    n : int, optional
        Length of the Fourier transform. If n<x.shape[axis],
        x is truncated. If n>x.shape[axis], x is zero-padded.
        (Default n=x.shape[axis]).
    axis : int, optional
        Axis along which the fft's are computed. (default=-1)
    overwrite_x : bool, optional
        If True the contents of x can be destroyed. (default=False)

    Returns
    -------
    z : complex ndarray
        with the elements:
            [y(0),y(1),..,y(n/2-1),y(-n/2),...,y(-1)]        if n is even
            [y(0),y(1),..,y((n-1)/2),y(-(n-1)/2),...,y(-1)]  if n is odd
        where
            y(j) = sum[k=0..n-1] x[k] * exp(-sqrt(-1)*j*k* 2*pi/n), j = 0..n-1
        Note that y(-j) = y(n-j).conjugate().

    See Also
    --------
    ifft : Inverse FFT
    rfft : FFT of a real sequence

    Notes
    -----
    The packing of the result is "standard": If A = fft(a, n), then A[0]
    contains the zero-frequency term, A[1:n/2+1] contains the
    positive-frequency terms, and A[n/2+1:] contains the negative-frequency
    terms, in order of decreasingly negative frequency. So for an 8-point
    transform, the frequencies of the result are [ 0, 1, 2, 3, 4, -3, -2, -1].

    This is most efficient for n a power of two.

    .. note:: In scipy 0.8.0 `fft` in single precision is available, but *only*
        for input array sizes which can be factorized into (combinations of) 2,
        3 and 5. For other sizes the computation will be done in double
        precision.

    Examples
    --------
    >>> x = np.arange(5)
    >>> np.all(np.abs(x-fft(ifft(x))<1.e-15) #within numerical accuracy.
    True

    """
    tmp = _asfarray(x)

    try:
        work_function = _DTYPE_TO_FFT[tmp.dtype]
    except KeyError:
        raise ValueError("type %s is not supported" % tmp.dtype)

    if istype(tmp, numpy.complex128):
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
    elif istype(tmp, numpy.complex64):
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
    else:
        overwrite_x = 1

    #return _raw_fft(tmp,n,axis,1,overwrite_x,work_function)
    if n is None:
        n = tmp.shape[axis]
    elif n != tmp.shape[axis]:
        tmp = _fix_shape(tmp,n,axis)
        overwrite_x = 1

    if axis == -1 or axis == len(tmp.shape) - 1:
        return work_function(tmp,n,1,0,overwrite_x)

    tmp = swapaxes(tmp, axis, -1)
    tmp = work_function(tmp,n,1,0,overwrite_x)
    return swapaxes(tmp, axis, -1)

def ifft(x, n=None, axis=-1, overwrite_x=0):
    """ ifft(x, n=None, axis=-1, overwrite_x=0) -> y

    Return inverse discrete Fourier transform of arbitrary type
    sequence x.

    The returned complex array contains
      [y(0),y(1),...,y(n-1)]
    where
      y(j) = 1/n sum[k=0..n-1] x[k] * exp(sqrt(-1)*j*k* 2*pi/n)

    Optional input: see fft.__doc__
    """
    tmp = _asfarray(x)

    try:
        work_function = _DTYPE_TO_FFT[tmp.dtype]
    except KeyError:
        raise ValueError("type %s is not supported" % tmp.dtype)

    if istype(tmp, numpy.complex128):
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
    elif istype(tmp, numpy.complex64):
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
    else:
        overwrite_x = 1

    #return _raw_fft(tmp,n,axis,-1,overwrite_x,work_function)
    if n is None:
        n = tmp.shape[axis]
    elif n != tmp.shape[axis]:
        tmp = _fix_shape(tmp,n,axis)
        overwrite_x = 1

    if axis == -1 or axis == len(tmp.shape) - 1:
        return work_function(tmp,n,-1,1,overwrite_x)

    tmp = swapaxes(tmp, axis, -1)
    tmp = work_function(tmp,n,-1,1,overwrite_x)
    return swapaxes(tmp, axis, -1)


def rfft(x, n=None, axis=-1, overwrite_x=0):
    """ rfft(x, n=None, axis=-1, overwrite_x=0) -> y

    Return discrete Fourier transform of real sequence x.

    The returned real arrays contains
      [y(0),Re(y(1)),Im(y(1)),...,Re(y(n/2))]              if n is even
      [y(0),Re(y(1)),Im(y(1)),...,Re(y(n/2)),Im(y(n/2))]   if n is odd
    where
      y(j) = sum[k=0..n-1] x[k] * exp(-sqrt(-1)*j*k* 2*pi/n)
      j = 0..n-1
    Note that y(-j) = y(n-j).conjugate().

    Optional input:
      n
        Defines the length of the Fourier transform. If n is not
        specified then n=x.shape[axis] is set. If n<x.shape[axis],
        x is truncated. If n>x.shape[axis], x is zero-padded.
      axis
        The transform is applied along the given axis of the input
        array (or the newly constructed array if n argument was used).
      overwrite_x
        If set to true, the contents of x can be destroyed.

    Notes:
      y == rfft(irfft(y)) within numerical accuracy.
    """
    tmp = _asfarray(x)

    if not numpy.isrealobj(tmp):
        raise TypeError,"1st argument must be real sequence"

    try:
        work_function = _DTYPE_TO_RFFT[tmp.dtype]
    except KeyError:
        raise ValueError("type %s is not supported" % tmp.dtype)

    return _raw_fft(tmp,n,axis,1,overwrite_x,work_function)


def rfftfreq(n,d=1.0):
    """ rfftfreq(n, d=1.0) -> f

    DFT sample frequencies (for usage with rfft,irfft).

    The returned float array contains the frequency bins in
    cycles/unit (with zero at the start) given a window length n and a
    sample spacing d:

      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2]/(d*n)   if n is even
      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2,n/2]/(d*n)   if n is odd
    """
    assert isinstance(n,int) or isinstance(n,integer)
    return array(range(1,n+1),dtype=int)/2/float(n*d)


def irfft(x, n=None, axis=-1, overwrite_x=0):
    """ irfft(x, n=None, axis=-1, overwrite_x=0) -> y

    Return inverse discrete Fourier transform of real sequence x.
    The contents of x is interpreted as the output of rfft(..)
    function.

    The returned real array contains
      [y(0),y(1),...,y(n-1)]
    where for n is even
      y(j) = 1/n (sum[k=1..n/2-1] (x[2*k-1]+sqrt(-1)*x[2*k])
                                   * exp(sqrt(-1)*j*k* 2*pi/n)
                  + c.c. + x[0] + (-1)**(j) x[n-1])
    and for n is odd
      y(j) = 1/n (sum[k=1..(n-1)/2] (x[2*k-1]+sqrt(-1)*x[2*k])
                                   * exp(sqrt(-1)*j*k* 2*pi/n)
                  + c.c. + x[0])
    c.c. denotes complex conjugate of preceeding expression.

    Optional input: see rfft.__doc__
    """
    tmp = _asfarray(x)
    if not numpy.isrealobj(tmp):
        raise TypeError,"1st argument must be real sequence"

    try:
        work_function = _DTYPE_TO_RFFT[tmp.dtype]
    except KeyError:
        raise ValueError("type %s is not supported" % tmp.dtype)

    return _raw_fft(tmp,n,axis,-1,overwrite_x,work_function)

def _raw_fftnd(x, s, axes, direction, overwrite_x, work_function):
    """ Internal auxiliary function for fftnd, ifftnd."""
    if s is None:
        if axes is None:
            s = x.shape
        else:
            s = numpy.take(x.shape, axes)

    s = tuple(s)
    if axes is None:
        noaxes = True
        axes = range(-x.ndim, 0)
    else:
        noaxes = False
    if len(axes) != len(s):
        raise ValueError("when given, axes and shape arguments "\
                         "have to be of the same length")

    # No need to swap axes, array is in C order
    if noaxes:
        for i in axes:
            x = _fix_shape(x, s[i], i)
        #print x.shape, s
        return work_function(x,s,direction,overwrite_x=overwrite_x)

    # We ordered axes, because the code below to push axes at the end of the
    # array assumes axes argument is in ascending order.
    id = numpy.argsort(axes)
    axes = [axes[i] for i in id]
    s = [s[i] for i in id]

    # Swap the request axes, last first (i.e. First swap the axis which ends up
    # at -1, then at -2, etc...), such as the request axes on which the
    # operation is carried become the last ones
    for i in range(1, len(axes)+1):
        x = numpy.swapaxes(x, axes[-i], -i)

    # We can now operate on the axes waxes, the p last axes (p = len(axes)), by
    # fixing the shape of the input array to 1 for any axis the fft is not
    # carried upon.
    waxes = range(x.ndim - len(axes), x.ndim)
    shape = numpy.ones(x.ndim)
    shape[waxes] = s

    for i in range(len(waxes)):
        x = _fix_shape(x, s[i], waxes[i])

    r = work_function(x, shape, direction, overwrite_x=overwrite_x)

    # reswap in the reverse order (first axis first, etc...) to get original
    # order
    for i in range(len(axes), 0, -1):
        r = numpy.swapaxes(r, -i, axes[-i])

    return r


def fftn(x, shape=None, axes=None, overwrite_x=0):
    """ fftn(x, shape=None, axes=None, overwrite_x=0) -> y

    Return multi-dimensional discrete Fourier transform of arbitrary
    type sequence x.

    The returned array contains

      y[j_1,..,j_d] = sum[k_1=0..n_1-1, ..., k_d=0..n_d-1]
         x[k_1,..,k_d] * prod[i=1..d] exp(-sqrt(-1)*2*pi/n_i * j_i * k_i)

    where d = len(x.shape) and n = x.shape.
    Note that y[..., -j_i, ...] = y[..., n_i-j_i, ...].conjugate().

    Optional input:
      shape
        Defines the shape of the Fourier transform. If shape is not
        specified then shape=take(x.shape,axes,axis=0).
        If shape[i]>x.shape[i] then the i-th dimension is padded with
        zeros. If shape[i]<x.shape[i], then the i-th dimension is
        truncated to desired length shape[i].
      axes
        The transform is applied along the given axes of the input
        array (or the newly constructed array if shape argument was
        used).
      overwrite_x
        If set to true, the contents of x can be destroyed.

    Notes:
      y == fftn(ifftn(y)) within numerical accuracy.
    """
    return _raw_fftn_dispatch(x, shape, axes, overwrite_x, 1)

def _raw_fftn_dispatch(x, shape, axes, overwrite_x, direction):
    tmp = _asfarray(x)

    try:
        work_function = _DTYPE_TO_FFTN[tmp.dtype]
    except KeyError:
        raise ValueError("type %s is not supported" % tmp.dtype)

    if istype(tmp, numpy.complex128):
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
    elif istype(tmp, numpy.complex64):
        pass
    else:
        overwrite_x = 1
    return _raw_fftnd(tmp,shape,axes,direction,overwrite_x,work_function)


def ifftn(x, shape=None, axes=None, overwrite_x=0):
    """
    Return inverse multi-dimensional discrete Fourier transform of
    arbitrary type sequence x.

    The returned array contains::

      y[j_1,..,j_d] = 1/p * sum[k_1=0..n_1-1, ..., k_d=0..n_d-1]
         x[k_1,..,k_d] * prod[i=1..d] exp(sqrt(-1)*2*pi/n_i * j_i * k_i)

    where ``d = len(x.shape)``, ``n = x.shape``, and ``p = prod[i=1..d] n_i``.

    For description of parameters see `fftn`.

    See Also
    --------
    fftn : for detailed information.

    """
    return _raw_fftn_dispatch(x, shape, axes, overwrite_x, -1)

def fft2(x, shape=None, axes=(-2,-1), overwrite_x=0):
    """
    2-D discrete Fourier transform.

    Return the two-dimensional discrete Fourier transform of the 2-D argument
    `x`.

    See Also
    --------
    fftn : for detailed information.

    """
    return fftn(x,shape,axes,overwrite_x)


def ifft2(x, shape=None, axes=(-2,-1), overwrite_x=0):
    """ ifft2(x, shape=None, axes=(-2,-1), overwrite_x=0) -> y

    Return inverse two-dimensional discrete Fourier transform of
    arbitrary type sequence x.

    See ifftn.__doc__ for more information.
    """
    return ifftn(x,shape,axes,overwrite_x)
