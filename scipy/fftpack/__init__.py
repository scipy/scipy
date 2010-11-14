"""
Home of discrete Fourier transform algorithms

Modules
=======

.. autosummary::
   :toctree: generated/

   basic - Basic discrete Fourier transform operators
   convolve - Convolution functions
   helper - TODO
   pseudo_diffs - Differential and pseudo-differential operators
   realtransforms - Real spectrum tranforms (DCT, DST, MDCT)

Functions
=========

Fast Fourier Transforms (FFTs)
------------------------------

.. autosummary::
   :toctree: generated/

   fft - Fast (discrete) Fourier Transform (FFT)
   ifft - Inverse FFT
   fft2 - Two dimensional FFT
   ifft2 - Two dimensional inverse FFT
   fftn - n-dimensional FFT
   ifftn - n-dimensional inverse FFT
   rfft - FFT of strictly real-valued sequence
   irfft - Inverse of rfft
   rfftfreq - DFT sample frequencies (specific to rfft and irfft)
   dct - Discrete cosine transform
   idct - Inverse discrete cosine transform

Differential and pseudo-differential operators
----------------------------------------------

.. autosummary::
   :toctree: generated/

   diff - Differentiation and integration of periodic sequences
   tilbert - Tilbert transform:         cs_diff(x,h,h)
   itilbert - Inverse Tilbert transform: sc_diff(x,h,h)
   hilbert - Hilbert transform:         cs_diff(x,inf,inf)
   ihilbert - Inverse Hilbert transform: sc_diff(x,inf,inf)
   cs_diff - cosh/sinh pseudo-derivative of periodic sequences
   sc_diff - sinh/cosh pseudo-derivative of periodic sequences
   ss_diff - sinh/sinh pseudo-derivative of periodic sequences
   cc_diff - cosh/cosh pseudo-derivative of periodic sequences
   shift - Shift periodic sequences

"""
#
# fftpack - Discrete Fourier Transform algorithms.
#
# Created: Pearu Peterson, August,September 2002

from info import __all__,__doc__

from fftpack_version import fftpack_version as __version__

from basic import *
from pseudo_diffs import *
from helper import *

from numpy.dual import register_func
for k in ['fft', 'ifft', 'fftn', 'ifftn', 'fft2', 'ifft2']:
    register_func(k, eval(k))
del k, register_func

from realtransforms import *
__all__.extend(['dct', 'idct'])

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
