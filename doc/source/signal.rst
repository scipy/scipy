=======================================
Signal processing (:mod:`scipy.signal`)
=======================================

.. module:: scipy.signal

Convolution
===========

.. autosummary::
   :toctree: generated/

   convolve
   correlate
   fftconvolve
   convolve2d
   correlate2d
   sepfir2d

B-splines
=========

.. autosummary::
   :toctree: generated/

   bspline
   gauss_spline
   cspline1d
   qspline1d
   cspline2d
   qspline2d
   spline_filter

Filtering
=========

.. autosummary::
   :toctree: generated/

   order_filter
   medfilt
   medfilt2d
   wiener

   symiirorder1
   symiirorder2
   lfilter
   lfiltic

   deconvolve

   hilbert
   get_window

   decimate
   detrend
   resample

Filter design
=============

.. autosummary::
   :toctree: generated/

   bilinear
   firwin
   freqs
   freqz
   iirdesign
   iirfilter
   kaiserord
   remez

   unique_roots
   residue
   residuez
   invres

Matlab-style IIR filter design
==============================

.. autosummary::
   :toctree: generated/

   butter
   buttord
   cheby1
   cheb1ord
   cheby2
   cheb2ord
   ellip
   ellipord
   bessel

Linear Systems
==============

.. autosummary::
   :toctree: generated/

   lti
   lsim
   lsim2
   impulse
   impulse2
   step
   step2

LTI Representations
===================

.. autosummary::
   :toctree: generated/

   tf2zpk
   zpk2tf
   tf2ss
   ss2tf
   zpk2ss
   ss2zpk

Waveforms
=========

.. autosummary::
   :toctree: generated/

   chirp
   gausspulse
   sawtooth
   square
   sweep_poly

Window functions
================

.. autosummary::
   :toctree: generated/

   get_window
   barthann
   bartlett
   blackman
   blackmanharris
   bohman
   boxcar
   chebwin
   flattop
   gaussian
   general_gaussian
   hamming
   hann
   kaiser
   nuttall
   parzen
   slepian
   triang

Wavelets
========

.. autosummary::
   :toctree: generated/

   cascade
   daub
   morlet
   qmf
