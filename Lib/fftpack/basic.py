"""
Discrete Fourier Transforms - basic.py 
"""
# Created by Pearu Peterson, August,September 2002

__all__ = ['fft','ifft','fftn','ifftn','rfft','irfft',
           'fft2','ifft2']

from scipy_base import asarray, zeros, swapaxes
import _fftpack as fftpack

import atexit
atexit.register(fftpack.destroy_zfft_cache)
atexit.register(fftpack.destroy_zfftnd_cache)
atexit.register(fftpack.destroy_drfft_cache)
del atexit

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
        z = zeros(s,x.typecode())
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
    """ fft(x, n=None, axis=-1, overwrite_x=0) -> y

    Return discrete Fourier transform of arbitrary type sequence x.

    The returned complex array contains
      [y(0),y(1),..,y(n/2-1),y(-n/2),...,y(-1)]        if n is even
      [y(0),y(1),..,y((n-1)/2),y(-(n-1)/2),...,y(-1)]  if n is odd
    where
      y(j) = sum[k=0..n-1] x[k] * exp(-sqrt(-1)*j*k* 2*pi/n)
      j = 0..n-1
    Note that y(-j) = y(n-j).

    Optional input:
      n
        Defines the length of the Fourier transform. If n is not
        specified then n=x.shape[axis] is set. If n<x.shape[axis],
        x is truncated. If n>x.shape[axis], x is zero-padded.
      axis
        The trasnform is applied along the given axis of the input
        array (or the newly constructed array if n argument was used).
      overwrite_x
        If set to true, the contents of x can be destroyed.

    Notes:
      y == fft(ifft(y)) within numerical accuracy.
    """
    tmp = asarray(x)
    t = tmp.typecode()
    if t=='D':
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
        work_function = fftpack.zfft
    elif t=='F':
        raise NotImplementedError
    else:
        overwrite_x = 1
        work_function = fftpack.zrfft

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
    tmp = asarray(x)
    t = tmp.typecode()
    if t=='D':
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
        work_function = fftpack.zfft
    elif t=='F':
        raise NotImplementedError
    else:
        overwrite_x = 1
        work_function = fftpack.zrfft

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
    Note that y(-j) = y(n-j).

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
    tmp = asarray(x)
    t = tmp.typecode()
    if t in 'DF':
        raise TypeError,"1st argument must be real sequence"
    work_function = fftpack.drfft
    return _raw_fft(tmp,n,axis,1,overwrite_x,work_function)


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
    tmp = asarray(x)
    t = tmp.typecode()
    if t in 'DF':
        raise TypeError,"1st argument must be real sequence"
    work_function = fftpack.drfft
    return _raw_fft(tmp,n,axis,-1,overwrite_x,work_function)


def _raw_fftnd(x, s, axes, direction, overwrite_x, work_function):
    """ Internal auxiliary function for fftnd, ifftnd."""
    if s is None:
        s = x.shape
    s = tuple(s)
    if s!=x.shape:
        assert len(s)<=len(x.shape)
        for i in range(-len(s),0):
            if x.shape[i]!=s[i]:
                x = _fix_shape(x,s[i],i)
    if axes is None:
        return work_function(x,s,direction,overwrite_x=overwrite_x)

    #XXX: should we allow/check for repeated indices in axes?
    # If allowed then using it has an effect of reducing the shape
    # implicitly. 
    s = [x.shape[i] for i in axes]
    s = [1]*(len(x.shape)-len(s)) + s
    swaps = []
    state = range(-len(s),0)
    for i in range(-1,-len(axes)-1,-1):
        j = state[axes[i]]
        if i==j: continue
        state[i],state[j]=state[j],state[i]
        swaps.append((j,i))
        x = swapaxes(x, j,i)
    r = work_function(x,s,direction,overwrite_x=overwrite_x)
    swaps.reverse()
    for i,j in swaps:
        r = swapaxes(r,i,j)
    return r


def fftn(x, shape=None, axes=None, overwrite_x=0):
    """ fftn(x, shape=None, axes=None, overwrite_x=0) -> y
    
    Return multi-dimensional discrete Fourier transform of arbitrary
    type sequence x.

    The returned array contains

      y[j_1,..,j_d] = sum[k_1=0..n_1-1, ..., k_d=0..n_d-1]
         x[k_1,..,k_d] * prod[i=1..d] exp(-sqrt(-1)*2*pi/n_i * j_i * k_i)

    where d = len(x.shape) and n = x.shape.
    Note that y[..., -j_i, ...] = y[..., n_i-j_i, ...].

    Optional input:
      shape
        Defines the shape of the Fourier transform. If shape is not
        specified then shape=take(x.shape,axes).
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
    tmp = asarray(x)
    t = tmp.typecode()
    if t=='D':
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
        work_function = fftpack.zfftnd
    elif t=='F':
        raise NotImplementedError
    else:
        overwrite_x = 1
        work_function = fftpack.zfftnd
    return _raw_fftnd(tmp,shape,axes,1,overwrite_x,work_function)


def ifftn(x, shape=None, axes=None, overwrite_x=0):
    """ ifftn(x, s=None, axes=None, overwrite_x=0) -> y
    
    Return inverse multi-dimensional discrete Fourier transform of
    arbitrary type sequence x.

    The returned array contains

      y[j_1,..,j_d] = 1/p * sum[k_1=0..n_1-1, ..., k_d=0..n_d-1]
         x[k_1,..,k_d] * prod[i=1..d] exp(sqrt(-1)*2*pi/n_i * j_i * k_i)

    where d = len(x.shape), n = x.shape, and p = prod[i=1..d] n_i.

    Optional input: see fftn.__doc__
    """
    tmp = asarray(x)
    t = tmp.typecode()
    if t=='D':
        overwrite_x = overwrite_x or (tmp is not x and not \
                                      hasattr(x,'__array__'))
        work_function = fftpack.zfftnd
    elif t=='F':
        raise NotImplementedError
    else:
        overwrite_x = 1
        work_function = fftpack.zfftnd
    return _raw_fftnd(tmp,shape,axes,-1,overwrite_x,work_function)


def fft2(x, shape=None, axes=(-2,-1), overwrite_x=0):
    """ fft2(x, shape=None, axes=(-2,-1), overwrite_x=0) -> y
    
    Return two-dimensional discrete Fourier transform of
    arbitrary type sequence x.

    See fftn.__doc__ for more information.
    """
    return fftn(x,shape,axes,overwrite_x)


def ifft2(x, shape=None, axes=(-2,-1), overwrite_x=0):
    """ ifft2(x, shape=None, axes=(-2,-1), overwrite_x=0) -> y
    
    Return inverse two-dimensional discrete Fourier transform of
    arbitrary type sequence x.

    See ifftn.__doc__ for more information.
    """
    return ifftn(x,shape,axes,overwrite_x)
