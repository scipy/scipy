"""
Differential and pseudo-differential operators.
"""
# Created by Pearu Peterson, September 2002

__all__ = ['diff',
           'tilbert','itilbert','hilbert','ihilbert',
           'cs_diff','cc_diff','sc_diff','ss_diff',
           'shift']

from scipy_base import pi, asarray, sin, cos, sinh, cosh, tanh
import convolve

import atexit
atexit.register(convolve.destroy_convolve_cache)
del atexit


_cache = {}
def diff(x,order=1,period=None,
            _cache = _cache):
    """ diff(x, order=1, period=2*pi) -> y

    Return k-th derivative (or integral) of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = pow(sqrt(-1)*j*2*pi/period, order) * x_j
      y_0 = 0 if order is not 0.

    Optional input:
      order
        The order of differentiation. Default order is 1. If order is
        negative, then integration is carried out under the assumption
        that x_0==0.
      period
        The assumed period of the sequence. Default is 2*pi.

    Notes:
      If sum(x)=0 then
          diff(diff(x,k),-k)==x (within numerical accuracy)
      For odd order and even len(x), the Nyquist mode is taken zero.
    """
    tmp = asarray(x)
    if order==0:
        return tmp
    if tmp.typecode() in ['F','D']:
        return diff(tmp.real,order,period)+1j*diff(tmp.imag,order,period)
    if period is not None:
        c = 2*pi/period
    else:
        c = 1.0
    n = len(x)
    omega = _cache.get((n,order,c))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k,order=order,c=c):
            if k:
                return pow(c*k,order)
            return 0
        omega = convolve.init_convolution_kernel(n,kernel,d=order,
                                                 zero_nyquist=1)
        _cache[(n,order,c)] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve(tmp,omega,swap_real_imag=order%2,
                             overwrite_x=overwrite_x)
del _cache


_cache = {}
def tilbert(x,h,period=None,
            _cache = _cache):
    """ tilbert(x, h, period=2*pi) -> y

    Return h-Tilbert transform of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = -sqrt(-1)*coth(j*h*2*pi/period) * x_j
      y_0 = 0

    Input:
      h
        Defines the parameter of the Tilbert transform.
      period
        The assumed period of the sequence. Default period is 2*pi.

    Notes:
      If sum(x)==0 and n=len(x) is odd then
        tilbert(itilbert(x)) == x
      If 2*pi*h/period is approximately 10 or larger then numerically
        tilbert == hilbert
      (theoretically oo-Tilbert == Hilbert).
      For even len(x), the Nyquist mode of x is taken zero.
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return tilbert(tmp.real,h,period)+\
               1j*tilbert(tmp.imag,h,period)
    if period is not None:
        h = h*2*pi/period
    n = len(x)
    omega = _cache.get((n,h))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k,h=h):
            if k: return -1/tanh(h*k)
            return 0
        omega = convolve.init_convolution_kernel(n,kernel,d=1)
        _cache[(n,h)] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve(tmp,omega,swap_real_imag=1,overwrite_x=overwrite_x)
del _cache


_cache = {}
def itilbert(x,h,period=None,
            _cache = _cache):
    """ itilbert(x, h, period=2*pi) -> y

    Return inverse h-Tilbert transform of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = sqrt(-1)*tanh(j*h*2*pi/period) * x_j
      y_0 = 0

    Optional input: see tilbert.__doc__
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return itilbert(tmp.real,h,period)+\
               1j*itilbert(tmp.imag,h,period)
    if period is not None:
        h = h*2*pi/period
    n = len(x)
    omega = _cache.get((n,h))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k,h=h):
            if k: return tanh(h*k)
            return 0
        omega = convolve.init_convolution_kernel(n,kernel,d=1)
        _cache[(n,h)] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve(tmp,omega,swap_real_imag=1,overwrite_x=overwrite_x)
del _cache


_cache = {}
def hilbert(x,
            _cache=_cache):
    """ hilbert(x) -> y

    Return Hilbert transform of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = -sqrt(-1)*sign(j) * x_j
      y_0 = 0

    Notes:
      If sum(x)==0 then
        hilbert(ihilbert(x)) == x
      For even len(x), the Nyquist mode of x is taken zero.
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return hilbert(tmp.real,tol)+1j*hilbert(tmp.imag)
    n = len(x)
    omega = _cache.get(n)
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k):
            if k>0: return -1
            elif k<0: return 1
            return 0
        omega = convolve.init_convolution_kernel(n,kernel,d=1)
        _cache[n] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')        
    return convolve.convolve(tmp,omega,swap_real_imag=1,overwrite_x=overwrite_x)
del _cache


def ihilbert(x):
    """ ihilbert(x) -> y

    Return inverse Hilbert transform of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = sqrt(-1)*sign(j) * x_j
      y_0 = 0
    """
    return -hilbert(x)


_cache = {}
def cs_diff(x, a, b, period=None,
            _cache = _cache):
    """ cs_diff(x, a, b, period=2*pi) -> y

    Return (a,b)-cosh/sinh pseudo-derivative of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = -sqrt(-1)*cosh(j*a*2*pi/period)/sinh(j*b*2*pi/period) * x_j
      y_0 = 0

    Input:
      a,b
        Defines the parameters of the cosh/sinh pseudo-differential
        operator.
      period
        The period of the sequence. Default period is 2*pi.

    Notes:
      For even len(x), the Nyquist mode of x is taken zero.
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return cs_diff(tmp.real,a,b,period)+\
               1j*cs_diff(tmp.imag,a,b,period)
    if period is not None:
        a = a*2*pi/period
        b = b*2*pi/period
    n = len(x)
    omega = _cache.get((n,a,b))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k,a=a,b=b):
            if k: return -cosh(a*k)/sinh(b*k)
            return 0
        omega = convolve.init_convolution_kernel(n,kernel,d=1)
        _cache[(n,a,b)] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve(tmp,omega,swap_real_imag=1,overwrite_x=overwrite_x)
del _cache


_cache = {}
def sc_diff(x, a, b, period=None,
            _cache = _cache):
    """ sc_diff(x, a, b, period=2*pi) -> y

    Return (a,b)-sinh/cosh pseudo-derivative of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = sqrt(-1)*sinh(j*a*2*pi/period)/cosh(j*b*2*pi/period) * x_j
      y_0 = 0

    Input:
      a,b
        Defines the parameters of the sinh/cosh pseudo-differential
        operator.
      period
        The period of the sequence x. Default is 2*pi.

    Notes:
      sc_diff(cs_diff(x,a,b),b,a) == x
      For even len(x), the Nyquist mode of x is taken zero.
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return sc_diff(tmp.real,a,b,period)+\
               1j*sc_diff(tmp.imag,a,b,period)
    if period is not None:
        a = a*2*pi/period
        b = b*2*pi/period
    n = len(x)
    omega = _cache.get((n,a,b))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k,a=a,b=b):
            if k: return sinh(a*k)/cosh(b*k)
            return 0
        omega = convolve.init_convolution_kernel(n,kernel,d=1)
        _cache[(n,a,b)] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve(tmp,omega,swap_real_imag=1,overwrite_x=overwrite_x)
del _cache


_cache = {}
def ss_diff(x, a, b, period=None,
            _cache = _cache):
    """ ss_diff(x, a, b, period=2*pi) -> y

    Return (a,b)-sinh/sinh pseudo-derivative of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = sinh(j*a*2*pi/period)/sinh(j*b*2*pi/period) * x_j
      y_0 = a/b * x_0

    Input:
      a,b
        Defines the parameters of the sinh/sinh pseudo-differential
        operator.
      period
        The period of the sequence x. Default is 2*pi.

    Notes:
      ss_diff(ss_diff(x,a,b),b,a) == x
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return ss_diff(tmp.real,a,b,period)+\
               1j*ss_diff(tmp.imag,a,b,period)
    if period is not None:
        a = a*2*pi/period
        b = b*2*pi/period
    n = len(x)
    omega = _cache.get((n,a,b))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k,a=a,b=b):
            if k: return sinh(a*k)/sinh(b*k)
            return float(a)/b
        omega = convolve.init_convolution_kernel(n,kernel)
        _cache[(n,a,b)] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve(tmp,omega,overwrite_x=overwrite_x)
del _cache


_cache = {}
def cc_diff(x, a, b, period=None,
            _cache = _cache):
    """ cc_diff(x, a, b, period=2*pi) -> y

    Return (a,b)-cosh/cosh pseudo-derivative of a periodic sequence x.

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

      y_j = cosh(j*a*2*pi/period)/cosh(j*b*2*pi/period) * x_j

    Input:
      a,b
        Defines the parameters of the sinh/sinh pseudo-differential
        operator.

    Optional input:
      period
        The period of the sequence x. Default is 2*pi.

    Notes:
      cc_diff(cc_diff(x,a,b),b,a) == x
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return cc_diff(tmp.real,a,b,period)+\
               1j*cc_diff(tmp.imag,a,b,period)
    if period is not None:
        a = a*2*pi/period
        b = b*2*pi/period
    n = len(x)
    omega = _cache.get((n,a,b))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel(k,a=a,b=b):
            return cosh(a*k)/cosh(b*k)
        omega = convolve.init_convolution_kernel(n,kernel)
        _cache[(n,a,b)] = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve(tmp,omega,overwrite_x=overwrite_x)
del _cache

_cache = {}
def shift(x, a, period=None,
          _cache = _cache):
    """ shift(x, a, period=2*pi) -> y

    Shift periodic sequence x by a: y(u) = x(u+a).

    If x_j and y_j are Fourier coefficients of periodic functions x
    and y, respectively, then

          y_j = exp(j*a*2*pi/period*sqrt(-1)) * x_f

    Optional input:
      period
        The period of the sequences x and y. Default period is 2*pi.
    """
    tmp = asarray(x)
    if tmp.typecode() in ['F','D']:
        return shift(tmp.real,a,period)+1j*shift(tmp.imag,a,period)
    if period is not None:
        a = a*2*pi/period
    n = len(x)
    omega = _cache.get((n,a))
    if omega is None:
        if len(_cache)>20:
            while _cache: _cache.popitem()
        def kernel_real(k,a=a): return cos(a*k)
        def kernel_imag(k,a=a): return sin(a*k)
        omega_real = convolve.init_convolution_kernel(n,kernel_real,d=0,
                                                      zero_nyquist=0)
        omega_imag = convolve.init_convolution_kernel(n,kernel_imag,d=1,
                                                      zero_nyquist=0)
        _cache[(n,a)] = omega_real,omega_imag
    else:
        omega_real,omega_imag = omega
    overwrite_x = tmp is not x and not hasattr(x,'__array__')
    return convolve.convolve_z(tmp,omega_real,omega_imag,
                               overwrite_x=overwrite_x)

del _cache

