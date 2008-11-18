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
   medfilt2
   wiener

   symiirorder1
   symiirorder2
   lfilter

   deconvolve

   hilbert
   get_window

   detrend
   resample

Filter design
=============

.. autosummary::
   :toctree: generated/

   remez
   firwin
   iirdesign
   iirfilter
   freqs
   freqz

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
   impulse
   step

LTI Reresentations
==================

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

   sawtooth
   square
   gausspulse
   chirp

Window functions
================

.. autosummary::
   :toctree: generated/

   boxcar
   triang
   parzen
   bohman
   blackman
   blackmanharris
   nuttall
   flattop
   bartlett
   hann
   barthann
   hamming
   kaiser
   gaussian
   general_gaussian
   slepian

Wavelets
========

.. autosummary::
   :toctree: generated/

   daub
   qmf
   cascade
