"""
Signal Processing Tools
=======================

 Convolution:
 
    convolve      --  N-dimensional convolution.
    correlate     --  N-dimensional correlation.
    convolve2d    --  2-dimensional convolution (more options).
    correlate2d   --  2-dimensional correlation (more options).
    sepfir2d      --  Convolve with a 2-D separable FIR filter.

 B-splines:
    
    bspline       --  B-spline basis function of order n.
    gauss_spline  --  Gaussian approximation to the B-spline basis function.
    cspline1d     --  Coefficients for 1-D cubic (3rd order) B-spline.
    qspline1d     --  Coefficients for 1-D quadratic (2nd order) B-spline.
    cspline2d     --  Coefficients for 2-D cubic (3rd order) B-spline.
    qspline2d     --  Coefficients for 2-D quadratic (2nd order) B-spline.
    spline_filter --  Smoothing spline (cubic) filtering of a rank-2 array. 

 Filtering:

    order_filter  --  N-dimensional order filter.
    medfilt       --  N-dimensional median filter.
    medfilt2      --  2-dimensional median filter (faster).
    wiener        --  N-dimensional wiener filter.

    symiirorder1  --  2nd-order IIR filter (cascade of first-order systems).
    symiirorder2  --  4th-order IIR filter (cascade of second-order systems).
    lfilter       --  1-dimensional FIR and IIR digital linear filtering.

    hilbert       --- Compute the analytic signal of a 1-d signal.
    get_window    --- Create FIR window.

    detrend       --- Remove linear and/or constant trends from data.

 Filter design:
 
    remez         --  Optimal FIR filter design.
    firwin        --- Windowed FIR filter design.    
    iirdesign     --- IIR filter design given bands and gains
    iirfilter     --- IIR filter design given order and critical frequencies
    freqs         --- Analog filter frequency response
    freqz         --- Digital filter frequency response


 matlab-style IIR filter design:
 
    butter (buttord)  -- Butterworth
    cheby1 (cheb1ord) -- Chebyshev Type I
    cheby2 (cheb2ord) -- Chebyshev Type II
    ellip (ellipord)  -- Elliptic (Cauer)
    bessel            -- Bessel (no order selection available -- try butterod)

 Linear Systems:

    lti -- linear time invariant system object.
    lsim -- continuous-time simulation of output to linear system.
    impulse -- impulse response of linear, time-invariant (LTI) system.
    step -- step response of continous-time LTI system.

 LTI Reresentations:
 
    tf2zpk -- transfer function to zero-pole-gain.
    zpk2tf -- zero-pole-gain to transfer function.
    tf2ss -- transfer function to state-space.
    ss2tf -- state-pace to transfer function.
    zpk2ss -- zero-pole-gain to state-space.
    ss2zpk -- state-space to pole-zero-gain.    
"""

postpone_import = 1
