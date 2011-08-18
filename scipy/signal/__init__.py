"""
=======================================
Signal processing (:mod:`scipy.signal`)
=======================================

.. module:: scipy.signal

Convolution
===========

.. autosummary::
   :toctree: generated/

   convolve    -- N-dimensional convolution.
   correlate   -- N-dimensional correlation.
   fftconvolve -- N-dimensional convolution using the FFT.
   convolve2d  -- 2-dimensional convolution (more options).
   correlate2d -- 2-dimensional correlation (more options).
   sepfir2d    -- Convolve with a 2-D separable FIR filter.

B-splines
=========

.. autosummary::
   :toctree: generated/

   bspline       -- B-spline basis function of order n.
   gauss_spline  -- Gaussian approximation to the B-spline basis function.
   cspline1d     -- Coefficients for 1-D cubic (3rd order) B-spline.
   qspline1d     -- Coefficients for 1-D quadratic (2nd order) B-spline.
   cspline2d     -- Coefficients for 2-D cubic (3rd order) B-spline.
   qspline2d     -- Coefficients for 2-D quadratic (2nd order) B-spline.
   spline_filter -- Smoothing spline (cubic) filtering of a rank-2 array.

Filtering
=========

.. autosummary::
   :toctree: generated/

   order_filter  -- N-dimensional order filter.
   medfilt       -- N-dimensional median filter.
   medfilt2d     -- 2-dimensional median filter (faster).
   wiener        -- N-dimensional wiener filter.

   symiirorder1  -- 2nd-order IIR filter (cascade of first-order systems).
   symiirorder2  -- 4th-order IIR filter (cascade of second-order systems).
   lfilter       -- 1-dimensional FIR and IIR digital linear filtering.
   lfiltic       -- Construct initial conditions for `lfilter`.
   lfilter_zi    -- Compute an initial state zi for the lfilter function that
                 -- corresponds to the steady state of the step response.
   filtfilt      -- A forward-backward filter.

   deconvolve    -- 1-d deconvolution using lfilter.

   hilbert       -- Compute the analytic signal of a 1-d signal.
   get_window    -- Create FIR window.

   decimate      -- Downsample a signal.
   detrend       -- Remove linear and/or constant trends from data.
   resample      -- Resample using Fourier method.

Filter design
=============

.. autosummary::
   :toctree: generated/

   bilinear      -- Digital filter from an analog filter using
                    -- the bilinear transform.
   firwin        -- Windowed FIR filter design, with frequency response
                    -- defined as pass and stop bands.
   firwin2       -- Windowed FIR filter design, with arbitrary frequency
                    -- response.
   freqs         -- Analog filter frequency response.
   freqz         -- Digital filter frequency response.
   iirdesign     -- IIR filter design given bands and gains.
   iirfilter     -- IIR filter design given order and critical frequencies.
   kaiser_atten  -- Compute the attenuation of a Kaiser FIR filter, given
                    -- the number of taps and the transition width at
                    -- discontinuities in the frequency response.
   kaiser_beta   -- Compute the Kaiser parameter beta, given the desired
                    -- FIR filter attenuation.
   kaiserord     -- Design a Kaiser window to limit ripple and width of
                    -- transition region.
   remez         -- Optimal FIR filter design.

   unique_roots  -- Unique roots and their multiplicities.
   residue       -- Partial fraction expansion of b(s) / a(s).
   residuez      -- Partial fraction expansion of b(z) / a(z).
   invres        -- Inverse partial fraction expansion.

Matlab-style IIR filter design
==============================

.. autosummary::
   :toctree: generated/

   butter -- Butterworth
   buttord
   cheby1 -- Chebyshev Type I
   cheb1ord
   cheby2 -- Chebyshev Type II
   cheb2ord
   ellip -- Elliptic (Cauer)
   ellipord
   bessel -- Bessel (no order selection available -- try butterod)

Continuous-Time Linear Systems
==============================

.. autosummary::
   :toctree: generated/

   lti      -- linear time invariant system object.
   lsim     -- continuous-time simulation of output to linear system.
   lsim2    -- like lsim, but `scipy.integrate.odeint` is used.
   impulse  -- impulse response of linear, time-invariant (LTI) system.
   impulse2 -- like impulse, but `scipy.integrate.odeint` is used.
   step     -- step response of continous-time LTI system.
   step2    -- like step, but `scipy.integrate.odeint` is used.

Discrete-Time Linear Systems
============================
   dlsim    -- simulation of output to a discrete-time linear system.
   dimpulse -- impulse response of a discrete-time LTI system.
   dstep    -- step response of a discrete-time LTI system.

LTI Representations
===================

.. autosummary::
   :toctree: generated/

   tf2zpk        -- transfer function to zero-pole-gain.
   zpk2tf        -- zero-pole-gain to transfer function.
   tf2ss         -- transfer function to state-space.
   ss2tf         -- state-pace to transfer function.
   zpk2ss        -- zero-pole-gain to state-space.
   ss2zpk        -- state-space to pole-zero-gain.
   cont2discrete -- continuous-time to discrete-time LTI conversion.

Waveforms
=========

.. autosummary::
   :toctree: generated/

   chirp       -- Frequency swept cosine signal, with several freq functions.
   gausspulse  -- Gaussian modulated sinusoid
   sawtooth    -- Periodic sawtooth
   square      -- Square wave
   sweep_poly  -- Frequency swept cosine signal; freq is arbitrary polynomial

Window functions
================

.. autosummary::
   :toctree: generated/

   get_window        -- Return a window of a given length and type.
   barthann          -- Bartlett-Hann window
   bartlett          -- Bartlett window
   blackman          -- Blackman window
   blackmanharris    -- Minimum 4-term Blackman-Harris window
   bohman            -- Bohman window
   boxcar            -- Boxcar window
   chebwin           -- Dolph-Chebyshev window
   flattop           -- Flat top window
   gaussian          -- Gaussian window
   general_gaussian  -- Generalized Gaussian window
   hamming           -- Hamming window
   hann              -- Hann window
   kaiser            -- Kaiser window
   nuttall           -- Nuttall's minimum 4-term Blackman-Harris window
   parzen            -- Parzen window
   slepian           -- Slepian window
   triang            -- Triangular window

Wavelets
========

.. autosummary::
   :toctree: generated/

   cascade  -- compute scaling function and wavelet from coefficients
   daub     -- return low-pass
   morlet   -- Complex Morlet wavelet.
   qmf      -- return quadrature mirror filter from low-pass

"""

import sigtools
from waveforms import *

# The spline module (a C extension) provides:
#     cspline2d, qspline2d, sepfir2d, symiirord1, symiirord2
from spline import *

from bsplines import *
from cont2discrete import *
from dltisys import *
from filter_design import *
from fir_filter_design import *
from ltisys import *
from windows import *
from signaltools import *
from spectral import *
from wavelets import *

__all__ = filter(lambda s: not s.startswith('_'), dir())
from numpy.testing import Tester
test = Tester().test
