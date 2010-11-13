"""Functions for FIR filter design."""

import numpy
from numpy import pi, ceil, cos, asarray
from scipy.special import sinc
import sigtools

# Some notes on function parameters:
#
# `cutoff` and `width` are given as a numbers between 0 and 1.  These
# are relative frequencies, expressed as a fraction of the Nyquist rate.
# For example, if the Nyquist rate is 2KHz, then width=0.15 is a width
# of 300 Hz.
#
# The `order` of a FIR filter is one less than the number of taps.
# This is a potential source of confusion, so in the following code,
# we will always use the number of taps as the parameterization of
# the 'size' of the filter. The "number of taps" means the number
# of coefficients, which is the same as the length of the impulse
# response of the filter.


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
    if  a > 50:
        beta = 0.1102 * (a - 8.7)
    elif a > 21:
        beta = 0.5842 * (a - 21)**0.4 + 0.07886 * (a - 21)
    else:
        beta = 0.0
    return beta


def kaiser_atten(N, width):
    """Compute the attenuation of a Kaiser FIR filter.
    
    Given the number of taps `N` and the transition width `width`, compute the
    attenuation `a` in dB, given by Kaiser's formula:

        a = 2.285 * (N - 1) * pi * width + 7.95

    Parameters
    ----------
    N : int
        The number of taps in the FIR filter.
    width : float
        The desired width of the transition region between passband and stopband
        (or, in general, at any discontinuity) for the filter.

    Returns
    -------
    a : float
        The attenuation of the ripple, in dB.
        
    See Also
    --------
    kaiserord, kaiser_beta
    """
    a = 2.285 * (N - 1) * pi * width + 7.95
    return a


def kaiserord(ripple, width):
    """Design a Kaiser window to limit ripple and width of transition region.

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
    N : int
        The length of the kaiser window.
    beta :
        The beta parameter for the kaiser window.

    Notes
    -----
    There are several ways to obtain the Kaiser window:

      signal.kaiser(N, beta, sym=0)
      signal.get_window(beta, N)
      signal.get_window(('kaiser', beta), N)

    The empirical equations discovered by Kaiser are used.

    See Also
    --------
    kaiser_beta, kaiser_atten

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
    N = (A - 7.95) / 2.285 / (pi * width) + 1

    return int(ceil(N)), beta


def firwin(N, cutoff, width=None, window='hamming', pass_zero=True, scale=True):
    """
    FIR filter design using the window method.
    
    This function computes the coefficients of a finite impulse response filter.
    The filter will have linear phase; it will be Type I if `N` is odd and
    Type II if `N` is even.
    
    Type II filters always have zero response at the Nyquist rate, so a
    ValueError exception is raised if firwin is called with `N` even and
    having a passband whose right end is at the Nyquist rate.

    Parameters
    ----------
    N : int
        Length of the filter (number of coefficients, i.e. the filter
        order + 1).  `N` must be even if a passband includes the Nyquist
        frequency.

    cutoff : float or 1D array_like
        Cutoff frequency of filter (normalized so that 1 corresponds to
        Nyquist or pi radians / sample) OR an array of cutoff frequencies
        (that is, band edges). In the latter case, the frequencies in 
        `cutoff` should be positive and monotonically increasing between
        0 and 1.  The values 0 and 1 must not be included in `cutoff`.

    width : float or None
        If `width` is not None, then assume it is the approximate width of
        the transition region (normalized so that 1 corresponds to pi)
        for use in Kaiser FIR filter design.  In this case, the `window`
        argument is ignored.

    window : string or tuple of string and parameter values
        Desired window to use. See `scipy.signal.get_window` for a list of
        windows and required parameters.

    pass_zero : bool
        If True, the gain at the frequency 0 (i.e. the "DC gain") is 1.
        Otherwise the DC gain is 0.

    scale : bool
        Set to True to scale the coefficients so that the frequency
        response is exactly unity at a certain frequency.
        That frequency is either:
            0 (DC) if the first passband starts at 0 (i.e. pass_zero
                is True);
            1 (Nyquist) if the first passband ends at 1 (i.e the
                filter is a single band highpass filter);
            center of first passband otherwise.

    Returns
    -------
    h : 1D ndarray
        Coefficients of length N FIR filter.

    Raises
    ------
    ValueError
        If any value in cutoff is less than or equal to 0 or greater
        than or equal to 1, if the values in cutoff are not strictly
        monotonically increasing, or if `N` is even but a passband
        includes the Nyquist frequency.

    Examples
    --------
    
    Low-pass from 0 to f::
    
    >>> firwin(N, f)
    
    Use a specific window function::
    
    >>> firwin(N, f, window='nuttall')
    
    High-pass ('stop' from 0 to f)::
     
    >>> firwin(N, f, pass_zero=False)

    Band-pass::
    
    >>> firwin(N, [f1, f2], pass_zero=False)
    
    Band-stop::
    
    >>> firwin(N, [f1, f2]) 

    Multi-band (passbands are [0, f1], [f2, f3] and [f4, 1])::

    >>>firwin(N, [f1, f2, f3, f4])
    
    Multi-band (passbands are [f1, f2] and [f3,f4])::
    
    >>> firwin(N, [f1, f2, f3, f4], pass_zero=False)

    """

    # The major enhancements to this function added in November 2010 were
    # developed by Tom Krauss (see ticket #902).

    cutoff = numpy.atleast_1d(cutoff)
    
    # Check for invalid input.
    if cutoff.ndim > 1:
        raise ValueError("The cutoff argument must be at most one-dimensional.")
    if cutoff.size == 0:
        raise ValueError("At least one cutoff frequency must be given.")
    if cutoff.min() <= 0 or cutoff.max() >= 1:
        raise ValueError("Invalid cutoff frequency: frequencies must be greater than 0 and less than 1.")
    if numpy.any(numpy.diff(cutoff) <= 0):
        raise ValueError("Invalid cutoff frequencies: the frequencies must be strictly increasing.")

    if width is not None:
        # A width was given.  Verify that it is a float, find the beta parameter
        # of the Kaiser window, and set `window`.  This overrides the value of
        # `window` passed in.
        if isinstance(width, float):
            atten = kaiser_atten(N, width)
            beta = kaiser_beta(atten)
            window = ('kaiser', beta)
        else:
            raise ValueError("Invalid value for width: %s", width)

    pass_nyquist = bool(cutoff.size & 1) ^ pass_zero
    if pass_nyquist and N % 2 == 0:
        raise ValueError("A filter with an even number of coefficients must "
                            "have zero response at the Nyquist rate.")

    # Insert 0 and/or 1 at the ends of cutoff so that the length of cutoff is even,
    # and each pair in cutoff corresponds to passband.
    cutoff = numpy.hstack(([0.0]*pass_zero, cutoff, [1.0]*pass_nyquist))

    # `bands` is a 2D array; each row gives the left and right edges of a passband.
    bands = cutoff.reshape(-1,2)

    # Build up the coefficients.
    alpha = 0.5 * (N-1)
    m = numpy.arange(0, N) - alpha
    h = 0
    for left, right in bands:
        h += right * sinc(right * m)
        h -= left * sinc(left * m)

    # Get and apply the window function.
    from signaltools import get_window
    win = get_window(window, N, fftbins=False)
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
        c = cos(pi * m * scale_frequency)
        s = numpy.sum(h * c)
        h /= s
 
    return h


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

    >>> bpass = sp.signal.remez(72, [0, 0.1, 0.2, 0.4, 0.45, 0.5], [0, 1, 0])
    >>> freq, response = sp.signal.freqz(bpass)
    >>> ampl = np.abs(response)

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax1 = fig.add_subplot(111)
    >>> ax1.semilogy(freq/(2*np.pi), ampl, 'b-') # freq in Hz
    [<matplotlib.lines.Line2D object at 0xf486790>]
    >>> plt.show()

    """
    # Convert type
    try:
        tnum = {'bandpass':1, 'differentiator':2, 'hilbert':3}[type]
    except KeyError:
        raise ValueError("Type must be 'bandpass', 'differentiator', or 'hilbert'")

    # Convert weight
    if weight is None:
        weight = [1] * len(desired)

    bands = asarray(bands).copy()
    return sigtools._remez(numtaps, bands, desired, weight, tnum, Hz,
                           maxiter, grid_density)
