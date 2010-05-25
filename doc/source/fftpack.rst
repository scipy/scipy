Fourier transforms (:mod:`scipy.fftpack`)
=========================================

.. module:: scipy.fftpack

Fast Fourier transforms
-----------------------

.. autosummary::
   :toctree: generated/

   fft
   ifft
   fftn
   ifftn
   fft2
   ifft2
   rfft
   irfft

Differential and pseudo-differential operators
----------------------------------------------

.. autosummary::
   :toctree: generated/

   diff
   tilbert
   itilbert
   hilbert
   ihilbert
   cs_diff
   sc_diff
   ss_diff
   cc_diff
   shift

Helper functions
----------------

.. autosummary::
   :toctree: generated/

   fftshift
   ifftshift
   dftfreq
   rfftfreq

Convolutions (:mod:`scipy.fftpack.convolve`)
--------------------------------------------

.. module:: scipy.fftpack.convolve

.. autosummary::
   :toctree: generated/

   convolve
   convolve_z
   init_convolution_kernel
   destroy_convolve_cache


Other (:mod:`scipy.fftpack._fftpack`)
-------------------------------------

.. module:: scipy.fftpack._fftpack

.. autosummary::
   :toctree: generated/

   drfft
   zfft
   zrfft
   zfftnd
   destroy_drfft_cache
   destroy_zfft_cache
   destroy_zfftnd_cache
