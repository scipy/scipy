"""
Discrete Fourier Transforms - helper.py
"""
# Created by Pearu Peterson, September 2002

__all__ = ['fftshift','ifftshift','fftfreq','rfftfreq']

import Numeric
import types

def fftshift(x,axes=None):
    """ fftshift(x, axes=None) -> y

    Shift zero-frequency component to center of spectrum.

    This function swaps half-spaces for all axes listed (defaults to all).

    Notes:
      If len(x) is even then the Nyquist component is y[0].
    """
    tmp = Numeric.asarray(x)
    ndim = len(tmp.shape)
    if axes is None:
        axes = range(ndim)
    y = tmp
    for k in axes:
        n = tmp.shape[k]
        p2 = (n+1)/2
        mylist = Numeric.concatenate((Numeric.arange(p2,n),Numeric.arange(p2)))
        y = Numeric.take(y,mylist,k)
    return y


def ifftshift(x,axes=None):
    """ ifftshift(x,axes=None) - > y

    Inverse of fftshift.
    """
    tmp = Numeric.asarray(x)
    ndim = len(tmp.shape)
    if axes is None:
        axes = range(ndim)
    y = tmp
    for k in axes:
        n = tmp.shape[k]
        p2 = n-(n+1)/2
        mylist = Numeric.concatenate((Numeric.arange(p2,n),Numeric.arange(p2)))
        y = Numeric.take(y,mylist,k)
    return y

def fftfreq(n,d=1.0):
    """ fftfreq(n, d=1.0) -> f

    DFT sample frequencies

    The returned float array contains the frequency bins in
    cycles/unit (with zero at the start) given a window length n and a
    sample spacing d:

      f = [0,1,...,n/2-1,-n/2,...,-1]/(d*n)         if n is even
      f = [0,1,...,(n-1)/2,-(n-1)/2,...,-1]/(d*n)   if n is odd
    """
    assert isinstance(n,types.IntType)
    k = range(0,(n-1)/2+1)+range(-(n/2),0)
    return Numeric.array(k,'d')/(n*d)

def rfftfreq(n,d=1.0):
    """ rfftfreq(n, d=1.0) -> f

    DFT sample frequencies (for usage with rfft,irfft).

    The returned float array contains the frequency bins in
    cycles/unit (with zero at the start) given a window length n and a
    sample spacing d:

      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2]/(d*n)   if n is even
      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2,n/2]/(d*n)   if n is odd
    """
    assert isinstance(n,types.IntType)
    return Numeric.array(range(1,n+1),'i')/2/float(n*d)

##########################################################################

def test(level=10):
    from scipy_test.testing import module_test
    module_test(__name__,__file__,level=level)

def test_suite(level=1):
    from scipy_test.testing import module_test_suite
    return module_test_suite(__name__,__file__,level=level)
