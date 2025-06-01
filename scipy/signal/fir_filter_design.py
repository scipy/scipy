# This file is not meant for public use and will be removed in SciPy v2.0.0.
# Use the `scipy.signal` namespace for importing the functions
# included below.

from scipy._lib.deprecation import _sub_module_deprecation

__all__ = [  # noqa: F822
    'kaiser_beta', 'kaiser_atten', 'kaiserord',
    'firwin', 'firwin2', 'remez', 'remezord', 'firls', 'minimum_phase',
    'firwin_2d',
]


def __dir__():
    return __all__


def __getattr__(name):
    return _sub_module_deprecation(sub_package="signal", module="fir_filter_design",
                                   private_modules=["_fir_filter_design"], all=__all__,
                                   attribute=name)
def kaiser_beta(a):
    """Compute the Kaiser parameter `beta`, given the attenuation `a`.

    Parameters
    ----------
    a : float
        The desired attenuation in the stopband and maximum ripple in
        the passband, in dB.  This should be a *positive* number.

    Returns
    -------
    beta : float
        The `beta` parameter to be used in the formula for a Kaiser window.

    References
    ----------
    Oppenheim, Schafer, "Discrete-Time Signal Processing", p.475-476.
    """
    if a > 50:
        beta = 0.1102 * (a - 8.7)
    elif a > 21:
        beta = 0.5842 * (a - 21) ** 0.4 + 0.07886 * (a - 21)
    else:
        beta = 0.0
    return beta


def kaiser_atten(numtaps, width):
    """Compute the attenuation of a Kaiser FIR filter.

    Given the number of taps `N` and the transition width `width`, compute the
    attenuation `a` in dB, given by Kaiser's formula:

        a = 2.285 * (N - 1) * pi * width + 7.95

    Parameters
    ----------
    numtaps : int
        The number of taps in the FIR filter.
    width : float
        The desired width of the transition region between passband and
        stopband (or, in general, at any discontinuity) for the filter.

    Returns
    -------
    a : float
        The attenuation of the ripple, in dB.

    See Also
    --------
    kaiserord, kaiser_beta
    """
    a = 2.285 * (numtaps - 1) * np.pi * width + 7.95
    return a


def kaiserord(ripple, width):
    """
    Design a Kaiser window to limit ripple and width of transition region.

    Parameters
    ----------
    ripple : float
        Positive number specifying maximum ripple in passband (dB) and minimum
        ripple in stopband.
    width : float
        Width of transition region (normalized so that 1 corresponds to pi
        radians / sample).

    Returns
    -------
    numtaps : int
        The length of the kaiser window.
    beta : float
        The beta parameter for the kaiser window.

    See Also
    --------
    kaiser_beta, kaiser_atten

    Notes
    -----
    There are several ways to obtain the Kaiser window:

    - ``signal.kaiser(numtaps, beta, sym=0)``
    - ``signal.get_window(beta, numtaps)``
    - ``signal.get_window(('kaiser', beta), numtaps)``

    The empirical equations discovered by Kaiser are used.

    References
    ----------
    Oppenheim, Schafer, "Discrete-Time Signal Processing", p.475-476.

    """
    A = abs(ripple)  # in case somebody is confused as to what's meant
    if A < 8:
        # Formula for N is not valid in this range.
        raise ValueError("Requested maximum ripple attentuation %f is too "
                         "small for the Kaiser formula." % A)
    beta = kaiser_beta(A)

    # Kaiser's formula (as given in Oppenheim and Schafer) is for the filter
    # order, so we have to add 1 to get the number of taps.
    numtaps = (A - 7.95) / 2.285 / (np.pi * width) + 1

    return int(ceil(numtaps)), beta


def firwin(numtaps, cutoff, width=None, window='hamming', pass_zero=True,
           scale=True, nyq=1.0):
    """
    FIR filter design using the window method.

    This function computes the coefficients of a finite impulse response
    filter.  The filter will have linear phase; it will be Type I if
    `numtaps` is odd and Type II if `numtaps` is even.

    Type II filters always have zero response at the Nyquist rate, so a
    ValueError exception is raised if firwin is called with `numtaps` even and
    having a passband whose right end is at the Nyquist rate.

    Parameters
    ----------
    numtaps : int
        Length of the filter (number of coefficients, i.e. the filter
        order + 1).  `numtaps` must be even if a passband includes the
        Nyquist frequency.
    cutoff : float or 1D array_like
        Cutoff frequency of filter (expressed in the same units as `nyq`)
        OR an array of cutoff frequencies (that is, band edges). In the
        latter case, the frequencies in `cutoff` should be positive and
        monotonically increasing between 0 and `nyq`.  The values 0 and
        `nyq` must not be included in `cutoff`.
    width : float or None, optional
        If `width` is not None, then assume it is the approximate width
        of the transition region (expressed in the same units as `nyq`)
        for use in Kaiser FIR filter design.  In this case, the `window`
        argument is ignored.
    window : string or tuple of string and parameter values, optional
        Desired window to use. See `scipy.signal.get_window` for a list
        of windows and required parameters.
    pass_zero : bool, optional
        If True, the gain at the frequency 0 (i.e. the "DC gain") is 1.
        Otherwise the DC gain is 0.
    scale : bool, optional
        Set to True to scale the coefficients so that the frequency
        response is exactly unity at a certain frequency.
        That frequency is either:

        - 0 (DC) if the first passband starts at 0 (i.e. pass_zero
          is True)
        - `nyq` (the Nyquist rate) if the first passband ends at
          `nyq` (i.e the filter is a single band highpass filter);
          center of first passband otherwise

    nyq : float, optional
        Nyquist frequency.  Each frequency in `cutoff` must be between 0
        and `nyq`.

    Returns
    -------
    h : (numtaps,) ndarray
        Coefficients of length `numtaps` FIR filter.

    Raises
    ------
    ValueError
        If any value in `cutoff` is less than or equal to 0 or greater
        than or equal to `nyq`, if the values in `cutoff` are not strictly
        monotonically increasing, or if `numtaps` is even but a passband
        includes the Nyquist frequency.

    See also
    --------
    scipy.signal.firwin2

    Examples
    --------
    Low-pass from 0 to f::

    >>> from scipy import signal
    >>> signal.firwin(numtaps, f)

    Use a specific window function::

    >>> signal.firwin(numtaps, f, window='nuttall')

    High-pass ('stop' from 0 to f)::

    >>> signal.firwin(numtaps, f, pass_zero=False)

    Band-pass::

    >>> signal.firwin(numtaps, [f1, f2], pass_zero=False)

    Band-stop::

    >>> signal.firwin(numtaps, [f1, f2])

    Multi-band (passbands are [0, f1], [f2, f3] and [f4, 1])::

    >>> signal.firwin(numtaps, [f1, f2, f3, f4])

    Multi-band (passbands are [f1, f2] and [f3,f4])::

    >>> signal.firwin(numtaps, [f1, f2, f3, f4], pass_zero=False)

    """

    # The major enhancements to this function added in November 2010 were
    # developed by Tom Krauss (see ticket #902).

    cutoff = np.atleast_1d(cutoff) / float(nyq)

    # Check for invalid input.
    if cutoff.ndim > 1:
        raise ValueError("The cutoff argument must be at most "
                         "one-dimensional.")
    if cutoff.size == 0:
        raise ValueError("At least one cutoff frequency must be given.")
    if cutoff.min() <= 0 or cutoff.max() >= 1:
        raise ValueError("Invalid cutoff frequency: frequencies must be "
                         "greater than 0 and less than nyq.")
    if np.any(np.diff(cutoff) <= 0):
        raise ValueError("Invalid cutoff frequencies: the frequencies "
                         "must be strictly increasing.")

    if width is not None:
        # A width was given.  Find the beta parameter of the Kaiser window
        # and set `window`.  This overrides the value of `window` passed in.
        atten = kaiser_atten(numtaps, float(width) / nyq)
        beta = kaiser_beta(atten)
        window = ('kaiser', beta)

    pass_nyquist = bool(cutoff.size & 1) ^ pass_zero
    if pass_nyquist and numtaps % 2 == 0:
        raise ValueError("A filter with an even number of coefficients must "
                         "have zero response at the Nyquist rate.")

    # Insert 0 and/or 1 at the ends of cutoff so that the length of cutoff
    # is even, and each pair in cutoff corresponds to passband.
    cutoff = np.hstack(([0.0] * pass_zero, cutoff, [1.0] * pass_nyquist))

    # `bands` is a 2D array; each row gives the left and right edges of
    # a passband.
    bands = cutoff.reshape(-1, 2)

    # Build up the coefficients.
    alpha = 0.5 * (numtaps - 1)
    m = np.arange(0, numtaps) - alpha
    h = 0
    for left, right in bands:
        h += right * sinc(right * m)
        h -= left * sinc(left * m)

    # Get and apply the window function.
    from .signaltools import get_window
    win = get_window(window, numtaps, fftbins=False)
    h *= win

    # Now handle scaling if desired.
    if scale:
        # Get the first passband.
        left, right = bands[0]
        if left == 0:
            scale_frequency = 0.0
        elif right == 1:
            scale_frequency = 1.0
        else:
            scale_frequency = 0.5 * (left + right)
        c = np.cos(np.pi * m * scale_frequency)
        s = np.sum(h * c)
        h /= s

    return h


# Original version of firwin2 from scipy ticket #457, submitted by "tash".
#
# Rewritten by Warren Weckesser, 2010.

def firwin2(numtaps, freq, gain, nfreqs=None, window='hamming', nyq=1.0,
            antisymmetric=False):
    """
    FIR filter design using the window method.

    From the given frequencies `freq` and corresponding gains `gain`,
    this function constructs an FIR filter with linear phase and
    (approximately) the given frequency response.

    Parameters
    ----------
    numtaps : int
        The number of taps in the FIR filter.  `numtaps` must be less than
        `nfreqs`.
    freq : array_like, 1D
        The frequency sampling points. Typically 0.0 to 1.0 with 1.0 being
        Nyquist.  The Nyquist frequency can be redefined with the argument
        `nyq`.
        The values in `freq` must be nondecreasing.  A value can be repeated
        once to implement a discontinuity.  The first value in `freq` must
        be 0, and the last value must be `nyq`.
    gain : array_like
        The filter gains at the frequency sampling points. Certain
        constraints to gain values, depending on the filter type, are applied,
        see Notes for details.
    nfreqs : int, optional
        The size of the interpolation mesh used to construct the filter.
        For most efficient behavior, this should be a power of 2 plus 1
        (e.g, 129, 257, etc).  The default is one more than the smallest
        power of 2 that is not less than `numtaps`.  `nfreqs` must be greater
        than `numtaps`.
    window : string or (string, float) or float, or None, optional
        Window function to use. Default is "hamming".  See
        `scipy.signal.get_window` for the complete list of possible values.
        If None, no window function is applied.
    nyq : float, optional
        Nyquist frequency.  Each frequency in `freq` must be between 0 and
        `nyq` (inclusive).
    antisymmetric : bool, optional
        Whether resulting impulse response is symmetric/antisymmetric.
        See Notes for more details.

    Returns
    -------
    taps : ndarray
        The filter coefficients of the FIR filter, as a 1-D array of length
        `numtaps`.

    See also
    --------
    scipy.signal.firwin

    Notes
    -----
    From the given set of frequencies and gains, the desired response is
    constructed in the frequency domain.  The inverse FFT is applied to the
    desired response to create the associated convolution kernel, and the
    first `numtaps` coefficients of this kernel, scaled by `window`, are
    returned.

    The FIR filter will have linear phase. The type of filter is determined by
    the value of 'numtaps` and `antisymmetric` flag.
    There are four possible combinations:

       - odd  `numtaps`, `antisymmetric` is False, type I filter is produced
       - even `numtaps`, `antisymmetric` is False, type II filter is produced
       - odd  `numtaps`, `antisymmetric` is True, type III filter is produced
       - even `numtaps`, `antisymmetric` is True, type IV filter is produced

    Magnitude response of all but type I filters are subjects to following
    constraints:

       - type II  -- zero at the Nyquist frequency
       - type III -- zero at zero and Nyquist frequencies
       - type IV  -- zero at zero frequency

    .. versionadded:: 0.9.0

    References
    ----------
    .. [1] Oppenheim, A. V. and Schafer, R. W., "Discrete-Time Signal
       Processing", Prentice-Hall, Englewood Cliffs, New Jersey (1989).
       (See, for example, Section 7.4.)

    .. [2] Smith, Steven W., "The Scientist and Engineer's Guide to Digital
       Signal Processing", Ch. 17. http://www.dspguide.com/ch17/1.htm

    Examples
    --------
    A lowpass FIR filter with a response that is 1 on [0.0, 0.5], and
    that decreases linearly on [0.5, 1.0] from 1 to 0:

    >>> from scipy import signal
    >>> taps = signal.firwin2(150, [0.0, 0.5, 1.0], [1.0, 1.0, 0.0])
    >>> print(taps[72:78])
    [-0.02286961 -0.06362756  0.57310236  0.57310236 -0.06362756 -0.02286961]

    """

    if len(freq) != len(gain):
        raise ValueError('freq and gain must be of same length.')

    if nfreqs is not None and numtaps >= nfreqs:
        raise ValueError(('ntaps must be less than nfreqs, but firwin2 was '
                          'called with ntaps=%d and nfreqs=%s') %
                         (numtaps, nfreqs))

    if freq[0] != 0 or freq[-1] != nyq:
        raise ValueError('freq must start with 0 and end with `nyq`.')
    d = np.diff(freq)
    if (d < 0).any():
        raise ValueError('The values in freq must be nondecreasing.')
    d2 = d[:-1] + d[1:]
    if (d2 == 0).any():
        raise ValueError('A value in freq must not occur more than twice.')

    if antisymmetric:
        if numtaps % 2 == 0:
            ftype = 4
        else:
            ftype = 3
    else:
        if numtaps % 2 == 0:
            ftype = 2
        else:
            ftype = 1

    if ftype == 2 and gain[-1] != 0.0:
        raise ValueError("A Type II filter must have zero gain at the "
                         "Nyquist rate.")
    elif ftype == 3 and (gain[0] != 0.0 or gain[-1] != 0.0):
        raise ValueError("A Type III filter must have zero gain at zero "
                         "and Nyquist rates.")
    elif ftype == 4 and gain[0] != 0.0:
        raise ValueError("A Type IV filter must have zero gain at zero rate.")

    if nfreqs is None:
        nfreqs = 1 + 2 ** int(ceil(log(numtaps, 2)))

    # Tweak any repeated values in freq so that interp works.
    eps = np.finfo(float).eps
    for k in range(len(freq)):
        if k < len(freq) - 1 and freq[k] == freq[k + 1]:
            freq[k] = freq[k] - eps
            freq[k + 1] = freq[k + 1] + eps

    # Linearly interpolate the desired response on a uniform mesh `x`.
    x = np.linspace(0.0, nyq, nfreqs)
    fx = np.interp(x, freq, gain)

    # Adjust the phases of the coefficients so that the first `ntaps` of the
    # inverse FFT are the desired filter coefficients.
    shift = np.exp(-(numtaps - 1) / 2. * 1.j * np.pi * x / nyq)
    if ftype > 2:
        shift *= 1j

    fx2 = fx * shift

    # Use irfft to compute the inverse FFT.
    out_full = irfft(fx2)

    if window is not None:
        # Create the window to apply to the filter coefficients.
        from .signaltools import get_window
        wind = get_window(window, numtaps, fftbins=False)
    else:
        wind = 1

    # Keep only the first `numtaps` coefficients in `out`, and multiply by
    # the window.
    out = out_full[:numtaps] * wind

    if ftype == 3:
        out[out.size // 2] = 0.0

    return out


def remez(numtaps, bands, desired, weight=None, Hz=1, type='bandpass',
          maxiter=25, grid_density=16):
    """
    Calculate the minimax optimal filter using the Remez exchange algorithm.

    Calculate the filter-coefficients for the finite impulse response
    (FIR) filter whose transfer function minimizes the maximum error
    between the desired gain and the realized gain in the specified
    frequency bands using the Remez exchange algorithm.

    Parameters
    ----------
    numtaps : int
        The desired number of taps in the filter. The number of taps is
        the number of terms in the filter, or the filter order plus one.
    bands : array_like
        A monotonic sequence containing the band edges in Hz.
        All elements must be non-negative and less than half the sampling
        frequency as given by `Hz`.
    desired : array_like
        A sequence half the size of bands containing the desired gain
        in each of the specified bands.
    weight : array_like, optional
        A relative weighting to give to each band region. The length of
        `weight` has to be half the length of `bands`.
    Hz : scalar, optional
        The sampling frequency in Hz. Default is 1.
    type : {'bandpass', 'differentiator', 'hilbert'}, optional
        The type of filter:

          'bandpass' : flat response in bands. This is the default.

          'differentiator' : frequency proportional response in bands.

          'hilbert' : filter with odd symmetry, that is, type III
                      (for even order) or type IV (for odd order)
                      linear phase filters.

    maxiter : int, optional
        Maximum number of iterations of the algorithm. Default is 25.
    grid_density : int, optional
        Grid density. The dense grid used in `remez` is of size
        ``(numtaps + 1) * grid_density``. Default is 16.

    Returns
    -------
    out : ndarray
        A rank-1 array containing the coefficients of the optimal
        (in a minimax sense) filter.

    See Also
    --------
    freqz : Compute the frequency response of a digital filter.

    References
    ----------
    .. [1] J. H. McClellan and T. W. Parks, "A unified approach to the
           design of optimum FIR linear phase digital filters",
           IEEE Trans. Circuit Theory, vol. CT-20, pp. 697-701, 1973.
    .. [2] J. H. McClellan, T. W. Parks and L. R. Rabiner, "A Computer
           Program for Designing Optimum FIR Linear Phase Digital
           Filters", IEEE Trans. Audio Electroacoust., vol. AU-21,
           pp. 506-525, 1973.

    Examples
    --------
    We want to construct a filter with a passband at 0.2-0.4 Hz, and
    stop bands at 0-0.1 Hz and 0.45-0.5 Hz. Note that this means that the
    behavior in the frequency ranges between those bands is unspecified and
    may overshoot.

    >>> from scipy import signal
    >>> bpass = signal.remez(72, [0, 0.1, 0.2, 0.4, 0.45, 0.5], [0, 1, 0])
    >>> freq, response = signal.freqz(bpass)
    >>> ampl = np.abs(response)

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax1 = fig.add_subplot(111)
    >>> ax1.semilogy(freq/(2*np.pi), ampl, 'b-')  # freq in Hz
    >>> plt.show()

    """
    # Convert type
    try:
        tnum = {'bandpass': 1, 'differentiator': 2, 'hilbert': 3}[type]
    except KeyError:
        raise ValueError("Type must be 'bandpass', 'differentiator', "
                         "or 'hilbert'")

    # Convert weight
    if weight is None:
        weight = [1] * len(desired)

    bands = np.asarray(bands).copy()
    return sigtools._remez(numtaps, bands, desired, weight, tnum, Hz,
                           maxiter, grid_density)


def oddround(x):
    """Return the nearest odd integer from x."""

    return x-np.mod(x, 2)+1


def oddceil(x):
    """Return the smallest odd integer not less than x."""

    return oddround(x+1)


def remlplen_herrmann(fp, fs, dp, ds):
    """
    Determine the length of the low pass filter with passband frequency
    fp, stopband frequency fs, passband ripple dp, and stopband ripple ds.
    fp and fs must be normalized with respect to the sampling frequency.
    Note that the filter order is one less than the filter length.

    Uses approximation algorithm described by Herrmann et al.:

    O. Herrmann, L.R. Raviner, and D.S.K. Chan, Practical Design Rules for
    Optimum Finite Impulse Response Low-Pass Digital Filters, Bell Syst. Tech.
    Jour., 52(6):769-799, Jul./Aug. 1973.
    """

    dF = fs-fp
    a = [5.309e-3, 7.114e-2, -4.761e-1, -2.66e-3, -5.941e-1, -4.278e-1]
    b = [11.01217, 0.51244]
    Dinf = np.log10(ds)*(a[0]*np.log10(dp)**2+a[1]*np.log10(dp)+a[2]) + \
        a[3]*np.log10(dp)**2+a[4]*np.log10(dp)+a[5]
    f = b[0]+b[1]*(np.log10(dp)-np.log10(ds))
    N1 = Dinf/dF-f*dF+1

    return int(oddround(N1))


def remlplen_kaiser(fp, fs, dp, ds):
    """
    Determine the length of the low pass filter with passband frequency
    fp, stopband frequency fs, passband ripple dp, and stopband ripple ds.
    fp and fs must be normalized with respect to the sampling frequency.
    Note that the filter order is one less than the filter length.

    Uses approximation algorithm described by Kaiser:

    J.F. Kaiser, Nonrecursive Digital Filter Design Using I_0-sinh Window
    function, Proc. IEEE Int. Symp. Circuits and Systems, 20-23, April 1974.
    """

    dF = fs-fp
    N2 = (-20*np.log10(np.sqrt(dp*ds))-13.0)/(14.6*dF)+1.0

    return int(oddceil(N2))


def remlplen_ichige(fp, fs, dp, ds):
    """
    Determine the length of the low pass filter with passband frequency
    fp, stopband frequency fs, passband ripple dp, and stopband ripple ds.
    fp and fs must be normalized with respect to the sampling frequency.
    Note that the filter order is one less than the filter length.

    Uses approximation algorithm described by Ichige et al.:

    K. Ichige, M. Iwaki, and R. Ishii, Accurate Estimation of Minimum
    Filter Length for Optimum FIR Digital Filters, IEEE Transactions on
    Circuits and Systems, 47(10):1008-1017, October 2000.
    """

    dF = fs-fp
    v = lambda dF, dp: 2.325*((-np.log10(dp))**-0.445)*dF**(-1.39)
    g = lambda fp, dF, d: (2.0/np.pi)*np.arctan(v(dF, dp)*(1.0/fp-1.0/(0.5-dF)))
    h = lambda fp, dF, c: (2.0/np.pi)*np.arctan((c/dF)*(1.0/fp-1.0/(0.5-dF)))
    Nc = ceil(1.0+(1.101/dF)*(-np.log10(2.0*dp))**1.1)
    Nm = (0.52/dF)*np.log10(dp/ds)*(-np.log10(dp))**0.17
    N3 = ceil(Nc*(g(fp, dF, dp)+g(0.5-dF-fp, dF, dp)+1.0)/3.0)
    DN = ceil(Nm*(h(fp, dF, 1.1)-(h(0.5-dF-fp, dF, 0.29)-1.0)/2.0))
    N4 = N3+DN

    return int(N4)


def remezord(freqs, amps, rips, Hz=1, alg='ichige'):
    """Filter parameter selection for the Remez exchange algorithm.

    Calculate the parameters required by the Remez exchange algorithm to
    construct a finite impulse response (FIR) filter that approximately
    meets the specified design.

    Paramters:
    ----------

    freqs : ndarray
        A monotonic sequence of band edges specified in Hertz. All
        elements must be non-negative and less than 1/2 the
        sampling frequency as given by the Hz parameter.
    amps : ndarray
        A sequence containing the amplitudes of the signal to be
        filtered over the various bands.
    rips : ndarray
        A sequence specifying the maximum ripples of each band.
    alg : string
        Filter length approximation algorithm. May be
        'herrmann', 'kaiser', or 'ichige'.

    Returns
    -------

    numtaps, bands, desired, weight : mixed
        See help for the remez function.

    See
    ---
    remez : Calculate the minimax optimal filter using the
        Remez exchange algorithm.

    Notes
    -----

    .. versionadded:: 0.15.0

    """

    # Make sure the parameters are floating point numpy arrays:
    freqs = np.asarray(freqs, 'd')
    amps = np.asarray(amps, 'd')
    rips = np.asarray(rips, 'd')

    # Scale ripples with respect to band amplitudes:
    rips /= (amps + (amps == 0.0))

    # Normalize input frequencies with respect to sampling frequency:
    freqs /= Hz

    # Select filter length approximation algorithm:
    if alg == 'herrmann':
        remlplen = remlplen_herrmann
    elif alg == 'kaiser':
        remlplen = remlplen_kaiser
    elif alg == 'ichige':
        remlplen = remlplen_ichige
    else:
        raise ValueError('Unknown filter length approximation algorithm.')

    # Validate inputs:
    if np.any(freqs > 0.5):
        raise ValueError('Frequency band edges must not exceed the Nyquist frequency.')
    if np.any(freqs < 0.0):
        raise ValueError('Frequency band edges must be nonnegative.')
    if np.any(rips < 0.0):
        raise ValueError('Ripples must be nonnegative.')
    if len(amps) != len(rips):
        raise ValueError('Number of amplitudes must equal number of ripples.')
    if len(freqs) != 2*(len(amps)-1):
        raise ValueError('Number of band edges must equal 2*((number of amplitudes)-1)')

    # Find the longest filter length needed to implement any of the
    # low-pass or high-pass filters with the specified edges:
    f1 = freqs[0:-1:2]
    f2 = freqs[1::2]
    L = 0
    for i in range(len(amps)-1):
        L = max((L,
                 remlplen(f1[i], f2[i], rips[i], rips[i+1]),
                 remlplen(0.5-f2[i], 0.5-f1[i], rips[i+1], rips[i])))

    # Cap the sequence of band edges with the limits of the digital frequency
    # range:
    bands = np.hstack((0.0, freqs, 0.5))

    # The filter design weights correspond to the ratios between the maximum
    # ripple and all of the other ripples:
    weight = max(rips)/rips

    return [L, bands, amps, weight]
