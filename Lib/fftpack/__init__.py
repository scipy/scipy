"""
fftpack - Discrete Fourier Transform package.

Fast Fourier Transforms:

  fft       --- FFT of arbitrary type periodic sequences
  ifft      --- Inverse of fft
  fftn      --- Multi-dimensional FFT
  ifftn     --- Inverse of fftn
  fft2      --- Two-dimensional FFT
  ifft2     --- Inverse of fft2
  rfft      --- FFT of real periodic sequences
  irfft     --- Inverse of rfft

Differential and pseudo-differential operators:

  diff      --- Differentiation and integration of periodic sequences
  tilbert   --- Tilbert transform:         cs_diff(x,h,h)
  itilbert  --- Inverse Tilbert transform: sc_diff(x,h,h)
  hilbert   --- Hilbert transform:         cs_diff(x,inf,inf)
  ihilbert  --- Inverse Hilbert transform: sc_diff(x,inf,inf)
  cs_diff   --- cosh/sinh pseudo-derivative of periodic sequences
  sc_diff   --- sinh/cosh pseudo-derivative of periodic sequences
  ss_diff   --- sinh/sinh pseudo-derivative of periodic sequences
  cc_diff   --- cosh/cosh pseudo-derivative of periodic sequences
  shift     --- Shift periodic sequences

Helper functions:

  fftshift  --- Shift zero-frequency component to center of spectrum
  ifftshift --- Inverse of freqshift
  dftfreq   --- DFT sample frequencies
  rfftfreq  --- DFT sample frequencies (specific to rfft,irfft)

Extension modules:

  fftpack   --- Provides functions zfft, drfft, zrfft, zfftnd,
                destroy_*_cache
  convolve  --- Provides functions convolve, convolve_z,
                init_convolution_kernel, destroy_convolve_cache
"""
# Created by Pearu Peterson, August,September 2002

__all__ = ['fft','ifft','fftn','ifftn','rfft','irfft',
           'fft2','ifft2',
           'diff',
           'tilbert','itilbert','hilbert','ihilbert',
           'sc_diff','cs_diff','cc_diff','ss_diff',
           'shift',
           'fftshift','ifftshift',
           'fftfreq','rfftfreq'
           ]

from basic import *
from pseudo_diffs import *
from helper import *

################## test functions #########################

def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    exec 'import %s as this_mod' % (__name__)
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
