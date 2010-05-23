"""
Signal Processing Tools
=======================

 Convolution:

    convolve:
        N-dimensional convolution.

    correlate:
        N-dimensional correlation.
    fftconvolve:
        N-dimensional convolution using the FFT.
    convolve2d:
        2-dimensional convolution (more options).
    correlate2d:
        2-dimensional correlation (more options).
    sepfir2d:
        Convolve with a 2-D separable FIR filter.

 B-splines:

    bspline:
        B-spline basis function of order n.
    gauss_spline:
        Gaussian approximation to the B-spline basis function.
    cspline1d:
        Coefficients for 1-D cubic (3rd order) B-spline.
    qspline1d:
        Coefficients for 1-D quadratic (2nd order) B-spline.
    cspline2d:
        Coefficients for 2-D cubic (3rd order) B-spline.
    qspline2d:
        Coefficients for 2-D quadratic (2nd order) B-spline.
    spline_filter:
        Smoothing spline (cubic) filtering of a rank-2 array.

 Filtering:

    order_filter:
        N-dimensional order filter.
    medfilt:
        N-dimensional median filter.
    medfilt2:
        2-dimensional median filter (faster).
    wiener:
        N-dimensional wiener filter.
    symiirorder1:
        2nd-order IIR filter (cascade of first-order systems).
    symiirorder2:
        4th-order IIR filter (cascade of second-order systems).
    lfilter:
        1-dimensional FIR and IIR digital linear filtering.
    lfiltic:
        Construct initial conditions for `lfilter`.
    deconvolve:
        1-d deconvolution using lfilter.
    hilbert:
        Compute the analytic signal of a 1-d signal.
    get_window:
        Create FIR window.
    decimate:
        Downsample a signal.
    detrend:
        Remove linear and/or constant trends from data.
    resample:
        Resample using Fourier method.

 Filter design:

    bilinear:
        Return a digital filter from an analog filter using the bilinear transform.
    firwin:
        Windowed FIR filter design.
    freqs:
        Analog filter frequency response.
    freqz:
        Digital filter frequency response.
    iirdesign:
        IIR filter design given bands and gains.
    iirfilter:
        IIR filter design given order and critical frequencies.
    invres:
        Inverse partial fraction expansion.
    kaiserord:
        Design a Kaiser window to limit ripple and width of transition region.
    remez:
        Optimal FIR filter design.
    residue:
        Partial fraction expansion of b(s) / a(s).
    residuez:
        Partial fraction expansion of b(z) / a(z).
    unique_roots:
        Unique roots and their multiplicities.

 Matlab-style IIR filter design:

    butter (buttord):
        Butterworth
    cheby1 (cheb1ord):
        Chebyshev Type I
    cheby2 (cheb2ord):
        Chebyshev Type II
    ellip (ellipord):
        Elliptic (Cauer)
    bessel:
        Bessel (no order selection available -- try butterod)

 Linear Systems:

    lti:
        linear time invariant system object.
    lsim:
        continuous-time simulation of output to linear system.
    lsim2:
        like lsim, but `scipy.integrate.odeint` is used.
    impulse:
        impulse response of linear, time-invariant (LTI) system.
    impulse2:
        like impulse, but `scipy.integrate.odeint` is used.
    step:
        step response of continous-time LTI system.
    step2:
        like step, but `scipy.integrate.odeint` is used.

 LTI Representations:

    tf2zpk:
        transfer function to zero-pole-gain.
    zpk2tf:
        zero-pole-gain to transfer function.
    tf2ss:
        transfer function to state-space.
    ss2tf:
        state-pace to transfer function.
    zpk2ss:
        zero-pole-gain to state-space.
    ss2zpk:
        state-space to pole-zero-gain.

 Waveforms:

    sawtooth:
        Periodic sawtooth
    square:
        Square wave
    gausspulse:
        Gaussian modulated sinusoid
    chirp:
        Frequency swept cosine signal, with several frequency functions.
    sweep_poly:
        Frequency swept cosine signal; frequency is arbitrary polynomial.

 Window functions:

    get_window:
        Return a window of a given length and type.
    barthann:
        Bartlett-Hann window
    bartlett:
        Bartlett window
    blackman:
        Blackman window
    blackmanharris:
        Minimum 4-term Blackman-Harris window
    bohman:
        Bohman window
    boxcar:
        Boxcar window
    chebwin:
        Dolph-Chebyshev window
    flattop:
        Flat top window
    gaussian:
        Gaussian window
    general_gaussian:
        Generalized Gaussian window
    hamming:
        Hamming window
    hann:
        Hann window
    kaiser:
        Kaiser window
    nuttall:
        Nuttall's minimum 4-term Blackman-Harris window
    parzen:
        Parzen window
    slepian:
        Slepian window
    triang:
        Triangular window

 Wavelets:

    daub:
        return low-pass
    qmf:
        return quadrature mirror filter from low-pass
    cascade:
        compute scaling function and wavelet from coefficients
    morlet:
        Complex Morlet wavelet.
"""

postpone_import = 1
