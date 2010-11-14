# This file is executed by __init__.py.
"""
Discrete Fourier transform algorithms

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

__all__ = ['fft','ifft','fftn','ifftn','rfft','irfft',
           'fft2','ifft2',
           'diff',
           'tilbert','itilbert','hilbert','ihilbert',
           'sc_diff','cs_diff','cc_diff','ss_diff',
           'shift',
           'rfftfreq'
           ]

if __doc__:
    __doc_title__ = __doc__.lstrip().split('\n',1)[0]
else:
    __doc_title__ = None

postpone_import = 1

global_symbols = ['fft','fftn','fft2','ifft','ifft2','ifftn',
                  'fftshift','ifftshift','fftfreq']
