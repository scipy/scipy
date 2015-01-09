"""
=======================================
Signal processing (:mod:`scipy.signal`)
=======================================

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

   bspline        -- B-spline basis function of order n.
   cubic          -- B-spline basis function of order 3.
   quadratic      -- B-spline basis function of order 2.
   gauss_spline   -- Gaussian approximation to the B-spline basis function.
   cspline1d      -- Coefficients for 1-D cubic (3rd order) B-spline.
   qspline1d      -- Coefficients for 1-D quadratic (2nd order) B-spline.
   cspline2d      -- Coefficients for 2-D cubic (3rd order) B-spline.
   qspline2d      -- Coefficients for 2-D quadratic (2nd order) B-spline.
   cspline1d_eval -- Evaluate a cubic spline at the given points.
   qspline1d_eval -- Evaluate a quadratic spline at the given points.
   spline_filter  -- Smoothing spline (cubic) filtering of a rank-2 array.

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
   savgol_filter -- Filter a signal using the Savitzky-Golay filter.

   deconvolve    -- 1-d deconvolution using lfilter.

   sosfilt       -- 1-dimensional IIR digital linear filtering using
                 -- a second-order-sections filter representation.
   sosfilt_zi    -- Compute an initial state zi for the sosfilt function that
                 -- corresponds to the steady state of the step response.
   hilbert       -- Compute 1-D analytic signal, using the Hilbert transform.
   hilbert2      -- Compute 2-D analytic signal, using the Hilbert transform.

   decimate      -- Downsample a signal.
   detrend       -- Remove linear and/or constant trends from data.
   resample      -- Resample using Fourier method.

Filter design
=============

.. autosummary::
   :toctree: generated/

   bilinear      -- Digital filter from an analog filter using
                    -- the bilinear transform.
   findfreqs     -- Find array of frequencies for computing filter response.
   firwin        -- Windowed FIR filter design, with frequency response
                    -- defined as pass and stop bands.
   firwin2       -- Windowed FIR filter design, with arbitrary frequency
                    -- response.
   freqs         -- Analog filter frequency response.
   freqz         -- Digital filter frequency response.
   group_delay   -- Digital filter group delay.
   iirdesign     -- IIR filter design given bands and gains.
   iirfilter     -- IIR filter design given order and critical frequencies.
   kaiser_atten  -- Compute the attenuation of a Kaiser FIR filter, given
                    -- the number of taps and the transition width at
                    -- discontinuities in the frequency response.
   kaiser_beta   -- Compute the Kaiser parameter beta, given the desired
                    -- FIR filter attenuation.
   kaiserord     -- Design a Kaiser window to limit ripple and width of
                    -- transition region.
   savgol_coeffs -- Compute the FIR filter coefficients for a Savitzky-Golay
                    -- filter.
   remez         -- Optimal FIR filter design.

   unique_roots  -- Unique roots and their multiplicities.
   residue       -- Partial fraction expansion of b(s) / a(s).
   residuez      -- Partial fraction expansion of b(z) / a(z).
   invres        -- Inverse partial fraction expansion for analog filter.
   invresz       -- Inverse partial fraction expansion for digital filter.

Lower-level filter design functions:

.. autosummary::
   :toctree: generated/

   abcd_normalize -- Check state-space matrices and ensure they are rank-2.
   band_stop_obj  -- Band Stop Objective Function for order minimization.
   besselap       -- Return (z,p,k) for analog prototype of Bessel filter.
   buttap         -- Return (z,p,k) for analog prototype of Butterworth filter.
   cheb1ap        -- Return (z,p,k) for type I Chebyshev filter.
   cheb2ap        -- Return (z,p,k) for type II Chebyshev filter.
   cmplx_sort     -- Sort roots based on magnitude.
   ellipap        -- Return (z,p,k) for analog prototype of elliptic filter.
   lp2bp          -- Transform a lowpass filter prototype to a bandpass filter.
   lp2bs          -- Transform a lowpass filter prototype to a bandstop filter.
   lp2hp          -- Transform a lowpass filter prototype to a highpass filter.
   lp2lp          -- Transform a lowpass filter prototype to a lowpass filter.
   normalize      -- Normalize polynomial representation of a transfer function.



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

   freqresp -- frequency response of a continuous-time LTI system.
   lti      -- linear time invariant system object.
   lsim     -- continuous-time simulation of output to linear system.
   lsim2    -- like lsim, but `scipy.integrate.odeint` is used.
   impulse  -- impulse response of linear, time-invariant (LTI) system.
   impulse2 -- like impulse, but `scipy.integrate.odeint` is used.
   step     -- step response of continous-time LTI system.
   step2    -- like step, but `scipy.integrate.odeint` is used.
   bode     -- Calculate Bode magnitude and phase data.

Discrete-Time Linear Systems
============================

.. autosummary::
   :toctree: generated/

   dlsim    -- simulation of output to a discrete-time linear system.
   dimpulse -- impulse response of a discrete-time LTI system.
   dstep    -- step response of a discrete-time LTI system.

LTI Representations
===================

.. autosummary::
   :toctree: generated/

   tf2zpk        -- transfer function to zero-pole-gain.
   tf2sos        -- transfer function to second-order sections.
   tf2ss         -- transfer function to state-space.
   zpk2tf        -- zero-pole-gain to transfer function.
   zpk2sos       -- zero-pole-gain to second-order sections.
   zpk2ss        -- zero-pole-gain to state-space.
   ss2tf         -- state-pace to transfer function.
   ss2zpk        -- state-space to pole-zero-gain.
   sos2zpk       -- second-order-sections to zero-pole-gain.
   sos2tf        -- second-order-sections to transfer function.
   cont2discrete -- continuous-time to discrete-time LTI conversion.
   place_poles   -- pole placement.
   
Waveforms
=========

.. autosummary::
   :toctree: generated/

   chirp       -- Frequency swept cosine signal, with several freq functions.
   gausspulse  -- Gaussian modulated sinusoid
   max_len_seq -- Maximum length sequence
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
   cosine            -- Cosine window
   exponential       -- Exponential window
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
   ricker   -- return ricker wavelet
   cwt      -- perform continuous wavelet transform

Peak finding
============

.. autosummary::
   :toctree: generated/

   find_peaks_cwt -- Attempt to find the peaks in the given 1-D array
   argrelmin      -- Calculate the relative minima of data
   argrelmax      -- Calculate the relative maxima of data
   argrelextrema  -- Calculate the relative extrema of data

Spectral Analysis
=================

.. autosummary::
   :toctree: generated/

   periodogram    -- Computes a (modified) periodogram
   welch          -- Compute a periodogram using Welch's method
   lombscargle    -- Computes the Lomb-Scargle periodogram
   vectorstrength -- Computes the vector strength

"""
from __future__ import division, print_function, absolute_import

from . import sigtools
from .waveforms import *
from ._max_len_seq import max_len_seq

# The spline module (a C extension) provides:
#     cspline2d, qspline2d, sepfir2d, symiirord1, symiirord2
from .spline import *

from .bsplines import *
from .cont2discrete import *
from .dltisys import *
from .filter_design import *
from .fir_filter_design import *
from .ltisys import *
from .windows import *
from .signaltools import *
from ._savitzky_golay import savgol_coeffs, savgol_filter
from .spectral import *
from .wavelets import *
from ._peak_finding import *

__all__ = [s for s in dir() if not s.startswith('_')]
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
