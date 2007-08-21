__all__ = ['fftshift', 'ifftshift', 'fftfreq', 'rfftfreq']

from numpy import array
from numpy.fft.helper import fftshift, ifftshift, fftfreq

def rfftfreq(n,d=1.0):
    """ rfftfreq(n, d=1.0) -> f

    DFT sample frequencies (for usage with rfft,irfft).

    The returned float array contains the frequency bins in
    cycles/unit (with zero at the start) given a window length n and a
    sample spacing d:

      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2]/(d*n)   if n is even
      f = [0,1,1,2,2,...,n/2-1,n/2-1,n/2,n/2]/(d*n)   if n is odd
    """
    assert isinstance(n,int)
    return array(range(1,n+1),dtype=int)/2/float(n*d)
