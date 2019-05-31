"""
==================================================
Discrete Fourier transforms (:mod:`scipy.fft`)
==================================================

Fast Fourier Transforms (FFTs)
==============================

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
   rfft2 - Two dimensional FFT of real sequence
   irfft2 - Inverse of rfft2
   rfftn - n-dimensional FFT of real sequence
   irfftn - Inverse of rfftn
   dct - Discrete cosine transform
   idct - Inverse discrete cosine transform
   dctn - n-dimensional Discrete cosine transform
   idctn - n-dimensional Inverse discrete cosine transform
   dst - Discrete sine transform
   idst - Inverse discrete sine transform
   dstn - n-dimensional Discrete sine transform
   idstn - n-dimensional Inverse discrete sine transform

Helper functions
================

.. autosummary::
   :toctree: generated/

   fftshift - Shift the zero-frequency component to the center of the spectrum
   ifftshift - The inverse of `fftshift`
   fftfreq - Return the Discrete Fourier Transform sample frequencies
   rfftfreq - DFT sample frequencies (for usage with rfft, irfft)
   next_fast_len - Find the optimal length to zero-pad an FFT for speed

Note that ``fftshift``, ``ifftshift`` and ``fftfreq`` are numpy functions
exposed by ``fftpack``; importing them from ``numpy`` should be preferred.

"""

from __future__ import division, print_function, absolute_import

from scipy.fft._pocketfft import (
    fft, ifft, fft2,ifft2, fftn, ifftn,
    rfft, irfft, rfft2, irfft2, rfftn, irfftn)

from scipy.fft._fftpack import(
    shift,
    fftfreq, rfftfreq,
    fftshift, ifftshift,
    next_fast_len,
    dct, idct, dst, idst, dctn, idctn, dstn, idstn)

from numpy.dual import register_func
for k in ['fft', 'ifft', 'fftn', 'ifftn', 'fft2', 'ifft2']:
    register_func(k, eval(k))
del k, register_func

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
