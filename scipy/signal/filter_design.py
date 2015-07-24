"""Filter design.
"""
from __future__ import division, print_function, absolute_import

import warnings

import numpy
from numpy import (atleast_1d, poly, polyval, roots, real, asarray, allclose,
                   resize, pi, absolute, logspace, r_, sqrt, tan, log10,
                   arctan, arcsinh, sin, exp, cosh, arccosh, ceil, conjugate,
                   zeros, sinh, append, concatenate, prod, ones, array)
from numpy import mintypecode
import numpy as np
from scipy import special, optimize
from scipy.special import comb

__all__ = ['findfreqs', 'freqs', 'freqz', 'tf2zpk', 'zpk2tf', 'normalize',
           'lp2lp', 'lp2hp', 'lp2bp', 'lp2bs', 'bilinear', 'iirdesign',
           'iirfilter', 'butter', 'cheby1', 'cheby2', 'ellip', 'bessel',
           'band_stop_obj', 'buttord', 'cheb1ord', 'cheb2ord', 'ellipord',
           'buttap', 'cheb1ap', 'cheb2ap', 'ellipap', 'besselap',
           'filter_dict', 'band_dict', 'BadCoefficients',
           'tf2sos', 'sos2tf', 'zpk2sos', 'sos2zpk', 'group_delay']


class BadCoefficients(UserWarning):
    pass

abs = absolute


def findfreqs(num, den, N):
    """
    Find an array of frequencies for computing the response of a filter.

    Parameters
    ----------
    num, den : array_like, 1-D
        The polynomial coefficients of the numerator and denominator of the
        transfer function of the filter or LTI system.  The coefficients are
        ordered from highest to lowest degree.
    N : int
        The length of the array to be computed.

    Returns
    -------
    w : (N,) ndarray
        A 1-D array of frequencies, logarithmically spaced.

    Examples
    --------
    Find a set of nine frequencies that span the "interesting part" of the
    frequency response for the filter with the transfer function

        H(s) = s / (s^2 + 8s + 25)

    >>> findfreqs([1, 0], [1, 8, 25], N=9)
    array([  1.00000000e-02,   3.16227766e-02,   1.00000000e-01,
             3.16227766e-01,   1.00000000e+00,   3.16227766e+00,
             1.00000000e+01,   3.16227766e+01,   1.00000000e+02])
    """
    ep = atleast_1d(roots(den)) + 0j
    tz = atleast_1d(roots(num)) + 0j

    if len(ep) == 0:
        ep = atleast_1d(-1000) + 0j

    ez = r_['-1',
            numpy.compress(ep.imag >= 0, ep, axis=-1),
            numpy.compress((abs(tz) < 1e5) & (tz.imag >= 0), tz, axis=-1)]

    integ = abs(ez) < 1e-10
    hfreq = numpy.around(numpy.log10(numpy.max(3 * abs(ez.real + integ) +
                                               1.5 * ez.imag)) + 0.5)
    lfreq = numpy.around(numpy.log10(0.1 * numpy.min(abs(real(ez + integ)) +
                                                     2 * ez.imag)) - 0.5)

    w = logspace(lfreq, hfreq, N)
    return w


def freqs(b, a, worN=None, plot=None):
    """
    Compute frequency response of analog filter.

    Given the numerator `b` and denominator `a` of a filter, compute its
    frequency response::

             b[0]*(jw)**(nb-1) + b[1]*(jw)**(nb-2) + ... + b[nb-1]
     H(w) = -------------------------------------------------------
             a[0]*(jw)**(na-1) + a[1]*(jw)**(na-2) + ... + a[na-1]

    Parameters
    ----------
    b : ndarray
        Numerator of a linear filter.
    a : ndarray
        Denominator of a linear filter.
    worN : {None, int}, optional
        If None, then compute at 200 frequencies around the interesting parts
        of the response curve (determined by pole-zero locations).  If a single
        integer, then compute at that many frequencies.  Otherwise, compute the
        response at the angular frequencies (e.g. rad/s) given in `worN`.
    plot : callable, optional
        A callable that takes two arguments. If given, the return parameters
        `w` and `h` are passed to plot. Useful for plotting the frequency
        response inside `freqs`.

    Returns
    -------
    w : ndarray
        The angular frequencies at which h was computed.
    h : ndarray
        The frequency response.

    See Also
    --------
    freqz : Compute the frequency response of a digital filter.

    Notes
    -----
    Using Matplotlib's "plot" function as the callable for `plot` produces
    unexpected results,  this plots the real part of the complex transfer
    function, not the magnitude.  Try ``lambda w, h: plot(w, abs(h))``.

    Examples
    --------
    >>> from scipy.signal import freqs, iirfilter

    >>> b, a = iirfilter(4, [1, 10], 1, 60, analog=True, ftype='cheby1')

    >>> w, h = freqs(b, a, worN=np.logspace(-1, 2, 1000))

    >>> import matplotlib.pyplot as plt
    >>> plt.semilogx(w, 20 * np.log10(abs(h)))
    >>> plt.xlabel('Frequency')
    >>> plt.ylabel('Amplitude response [dB]')
    >>> plt.grid()
    >>> plt.show()

    """
    if worN is None:
        w = findfreqs(b, a, 200)
    elif isinstance(worN, int):
        N = worN
        w = findfreqs(b, a, N)
    else:
        w = worN
    w = atleast_1d(w)
    s = 1j * w
    h = polyval(b, s) / polyval(a, s)
    if plot is not None:
        plot(w, h)

    return w, h


def freqz(b, a=1, worN=None, whole=0, plot=None):
    """
    Compute the frequency response of a digital filter.

    Given the numerator `b` and denominator `a` of a digital filter,
    compute its frequency response::

               jw               -jw            -jmw
        jw  B(e)    b[0] + b[1]e + .... + b[m]e
     H(e) = ---- = ------------------------------------
               jw               -jw            -jnw
            A(e)    a[0] + a[1]e + .... + a[n]e

    Parameters
    ----------
    b : ndarray
        numerator of a linear filter
    a : ndarray
        denominator of a linear filter
    worN : {None, int, array_like}, optional
        If None (default), then compute at 512 frequencies equally spaced
        around the unit circle.
        If a single integer, then compute at that many frequencies.
        If an array_like, compute the response at the frequencies given (in
        radians/sample).
    whole : bool, optional
        Normally, frequencies are computed from 0 to the Nyquist frequency,
        pi radians/sample (upper-half of unit-circle).  If `whole` is True,
        compute frequencies from 0 to 2*pi radians/sample.
    plot : callable
        A callable that takes two arguments. If given, the return parameters
        `w` and `h` are passed to plot. Useful for plotting the frequency
        response inside `freqz`.

    Returns
    -------
    w : ndarray
        The normalized frequencies at which h was computed, in radians/sample.
    h : ndarray
        The frequency response.

    Notes
    -----
    Using Matplotlib's "plot" function as the callable for `plot` produces
    unexpected results,  this plots the real part of the complex transfer
    function, not the magnitude.  Try ``lambda w, h: plot(w, abs(h))``.

    Examples
    --------
    >>> from scipy import signal
    >>> b = signal.firwin(80, 0.5, window=('kaiser', 8))
    >>> w, h = signal.freqz(b)

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> plt.title('Digital filter frequency response')
    >>> ax1 = fig.add_subplot(111)

    >>> plt.plot(w, 20 * np.log10(abs(h)), 'b')
    >>> plt.ylabel('Amplitude [dB]', color='b')
    >>> plt.xlabel('Frequency [rad/sample]')

    >>> ax2 = ax1.twinx()
    >>> angles = np.unwrap(np.angle(h))
    >>> plt.plot(w, angles, 'g')
    >>> plt.ylabel('Angle (radians)', color='g')
    >>> plt.grid()
    >>> plt.axis('tight')
    >>> plt.show()

    """
    b, a = map(atleast_1d, (b, a))
    if whole:
        lastpoint = 2 * pi
    else:
        lastpoint = pi
    if worN is None:
        N = 512
        w = numpy.linspace(0, lastpoint, N, endpoint=False)
    elif isinstance(worN, int):
        N = worN
        w = numpy.linspace(0, lastpoint, N, endpoint=False)
    else:
        w = worN
    w = atleast_1d(w)
    zm1 = exp(-1j * w)
    h = polyval(b[::-1], zm1) / polyval(a[::-1], zm1)
    if plot is not None:
        plot(w, h)

    return w, h


def group_delay(system, w=None, whole=False):
    r"""Compute the group delay of a digital filter.

    The group delay measures by how many samples amplitude envelopes of
    various spectral components of a signal are delayed by a filter.
    It is formally defined as the derivative of continuous (unwrapped) phase::

               d        jw
     D(w) = - -- arg H(e)
              dw

    Parameters
    ----------
    system : tuple of array_like (b, a)
        Numerator and denominator coefficients of a filter transfer function.
    w : {None, int, array-like}, optional
        If None (default), then compute at 512 frequencies equally spaced
        around the unit circle.
        If a single integer, then compute at that many frequencies.
        If array, compute the delay at the frequencies given
        (in radians/sample).
    whole : bool, optional
        Normally, frequencies are computed from 0 to the Nyquist frequency,
        pi radians/sample (upper-half of unit-circle).  If `whole` is True,
        compute frequencies from 0 to ``2*pi`` radians/sample.

    Returns
    -------
    w : ndarray
        The normalized frequencies at which the group delay was computed,
        in radians/sample.
    gd : ndarray
        The group delay.

    Notes
    -----
    The similar function in MATLAB is called `grpdelay`.

    If the transfer function :math:`H(z)` has zeros or poles on the unit
    circle, the group delay at corresponding frequencies is undefined.
    When such a case arises the warning is raised and the group delay
    is set to 0 at those frequencies.

    For the details of numerical computation of the group delay refer to [1]_.

    .. versionadded: 0.16.0

    See Also
    --------
    freqz : Frequency response of a digital filter

    References
    ----------
    .. [1] Richard G. Lyons, "Understanding Digital Signal Processing,
           3rd edition", p. 830.

    Examples
    --------
    >>> from scipy import signal
    >>> b, a = signal.iirdesign(0.1, 0.3, 5, 50, ftype='cheby1')
    >>> w, gd = signal.group_delay((b, a))

    >>> import matplotlib.pyplot as plt
    >>> plt.title('Digital filter group delay')
    >>> plt.plot(w, gd)
    >>> plt.ylabel('Group delay [samples]')
    >>> plt.xlabel('Frequency [rad/sample]')
    >>> plt.show()

    """
    if w is None:
        w = 512

    if isinstance(w, int):
        if whole:
            w = np.linspace(0, 2 * pi, w, endpoint=False)
        else:
            w = np.linspace(0, pi, w, endpoint=False)

    w = np.atleast_1d(w)
    b, a = map(np.atleast_1d, system)
    c = np.convolve(b, a[::-1])
    cr = c * np.arange(c.size)
    z = np.exp(-1j * w)
    num = np.polyval(cr[::-1], z)
    den = np.polyval(c[::-1], z)
    singular = np.absolute(den) < 10 * EPSILON
    if np.any(singular):
        warnings.warn(
            "The group delay is singular at frequencies [{0}], setting to 0".
            format(", ".join("{0:.3f}".format(ws) for ws in w[singular]))
        )

    gd = np.zeros_like(w)
    gd[~singular] = np.real(num[~singular] / den[~singular]) - a.size + 1
    return w, gd


def _cplxreal(z, tol=None):
    """
    Split into complex and real parts, combining conjugate pairs.

    The 1D input vector `z` is split up into its complex (`zc`) and real (`zr`)
    elements.  Every complex element must be part of a complex-conjugate pair,
    which are combined into a single number (with positive imaginary part) in
    the output.  Two complex numbers are considered a conjugate pair if their
    real and imaginary parts differ in magnitude by less than ``tol * abs(z)``.

    Parameters
    ----------
    z : array_like
        Vector of complex numbers to be sorted and split
    tol : float, optional
        Relative tolerance for testing realness and conjugate equality.
        Default is ``100 * spacing(1)`` of `z`'s data type (i.e. 2e-14 for
        float64)

    Returns
    -------
    zc : ndarray
        Complex elements of `z`, with each pair represented by a single value
        having positive imaginary part, sorted first by real part, and then
        by magnitude of imaginary part.  The pairs are averaged when combined
        to reduce error.
    zr : ndarray
        Real elements of `z` (those having imaginary part less than
        `tol` times their magnitude), sorted by value.

    Raises
    ------
    ValueError
        If there are any complex numbers in `z` for which a conjugate
        cannot be found.

    See Also
    --------
    _cplxpair

    Examples
    --------
    >>> a = [4, 3, 1, 2-2j, 2+2j, 2-1j, 2+1j, 2-1j, 2+1j, 1+1j, 1-1j]
    >>> zc, zr = _cplxreal(a)
    >>> print zc
    [ 1.+1.j  2.+1.j  2.+1.j  2.+2.j]
    >>> print zr
    [ 1.  3.  4.]
    """

    z = atleast_1d(z)
    if z.size == 0:
        return z, z
    elif z.ndim != 1:
        raise ValueError('_cplxreal only accepts 1D input')

    if tol is None:
        # Get tolerance from dtype of input
        tol = 100 * np.finfo((1.0 * z).dtype).eps

    # Sort by real part, magnitude of imaginary part (speed up further sorting)
    z = z[np.lexsort((abs(z.imag), z.real))]

    # Split reals from conjugate pairs
    real_indices = abs(z.imag) <= tol * abs(z)
    zr = z[real_indices].real

    if len(zr) == len(z):
        # Input is entirely real
        return array([]), zr

    # Split positive and negative halves of conjugates
    z = z[~real_indices]
    zp = z[z.imag > 0]
    zn = z[z.imag < 0]

    if len(zp) != len(zn):
        raise ValueError('Array contains complex value with no matching '
                         'conjugate.')

    # Find runs of (approximately) the same real part
    same_real = np.diff(zp.real) <= tol * abs(zp[:-1])
    diffs = numpy.diff(concatenate(([0], same_real, [0])))
    run_starts = numpy.where(diffs > 0)[0]
    run_stops = numpy.where(diffs < 0)[0]

    # Sort each run by their imaginary parts
    for i in range(len(run_starts)):
        start = run_starts[i]
        stop = run_stops[i] + 1
        for chunk in (zp[start:stop], zn[start:stop]):
            chunk[...] = chunk[np.lexsort([abs(chunk.imag)])]

    # Check that negatives match positives
    if any(abs(zp - zn.conj()) > tol * abs(zn)):
        raise ValueError('Array contains complex value with no matching '
                         'conjugate.')

    # Average out numerical inaccuracy in real vs imag parts of pairs
    zc = (zp + zn.conj()) / 2

    return zc, zr


def _cplxpair(z, tol=None):
    """
    Sort into pairs of complex conjugates.

    Complex conjugates in `z` are sorted by increasing real part.  In each
    pair, the number with negative imaginary part appears first.

    If pairs have identical real parts, they are sorted by increasing
    imaginary magnitude.

    Two complex numbers are considered a conjugate pair if their real and
    imaginary parts differ in magnitude by less than ``tol * abs(z)``.  The
    pairs are forced to be exact complex conjugates by averaging the positive
    and negative values.

    Purely real numbers are also sorted, but placed after the complex
    conjugate pairs.  A number is considered real if its imaginary part is
    smaller than `tol` times the magnitude of the number.

    Parameters
    ----------
    z : array_like
        1-dimensional input array to be sorted.
    tol : float, optional
        Relative tolerance for testing realness and conjugate equality.
        Default is ``100 * spacing(1)`` of `z`'s data type (i.e. 2e-14 for
        float64)

    Returns
    -------
    y : ndarray
        Complex conjugate pairs followed by real numbers.

    Raises
    ------
    ValueError
        If there are any complex numbers in `z` for which a conjugate
        cannot be found.

    See Also
    --------
    _cplxreal

    Examples
    --------
    >>> a = [4, 3, 1, 2-2j, 2+2j, 2-1j, 2+1j, 2-1j, 2+1j, 1+1j, 1-1j]
    >>> z = _cplxpair(a)
    >>> print(z)
    [ 1.-1.j  1.+1.j  2.-1.j  2.+1.j  2.-1.j  2.+1.j  2.-2.j  2.+2.j  1.+0.j
      3.+0.j  4.+0.j]
    """

    z = atleast_1d(z)
    if z.size == 0 or np.isrealobj(z):
        return np.sort(z)

    if z.ndim != 1:
        raise ValueError('z must be 1-dimensional')

    zc, zr = _cplxreal(z, tol)

    # Interleave complex values and their conjugates, with negative imaginary
    # parts first in each pair
    zc = np.dstack((zc.conj(), zc)).flatten()
    z = np.append(zc, zr)
    return z


def tf2zpk(b, a):
    r"""Return zero, pole, gain (z, p, k) representation from a numerator,
    denominator representation of a linear filter.

    Parameters
    ----------
    b : array_like
        Numerator polynomial coefficients.
    a : array_like
        Denominator polynomial coefficients.

    Returns
    -------
    z : ndarray
        Zeros of the transfer function.
    p : ndarray
        Poles of the transfer function.
    k : float
        System gain.

    Notes
    -----
    If some values of `b` are too close to 0, they are removed. In that case,
    a BadCoefficients warning is emitted.

    The `b` and `a` arrays are interpreted as coefficients for positive,
    descending powers of the transfer function variable.  So the inputs
    :math:`b = [b_0, b_1, ..., b_M]` and :math:`a =[a_0, a_1, ..., a_N]`
    can represent an analog filter of the form:

    .. math::

        H(s) = \frac
        {b_0 s^M + b_1 s^{(M-1)} + \cdots + b_M}
        {a_0 s^N + a_1 s^{(N-1)} + \cdots + a_N}

    or a discrete-time filter of the form:

    .. math::

        H(z) = \frac
        {b_0 z^M + b_1 z^{(M-1)} + \cdots + b_M}
        {a_0 z^N + a_1 z^{(N-1)} + \cdots + a_N}

    This "positive powers" form is found more commonly in controls
    engineering.  If `M` and `N` are equal (which is true for all filters
    generated by the bilinear transform), then this happens to be equivalent
    to the "negative powers" discrete-time form preferred in DSP:

    .. math::

        H(z) = \frac
        {b_0 + b_1 z^{-1} + \cdots + b_M z^{-M}}
        {a_0 + a_1 z^{-1} + \cdots + a_N z^{-N}}

    Although this is true for common filters, remember that this is not true
    in the general case.  If `M` and `N` are not equal, the discrete-time
    transfer function coefficients must first be converted to the "positive
    powers" form before finding the poles and zeros.

    """
    b, a = normalize(b, a)
    b = (b + 0.0) / a[0]
    a = (a + 0.0) / a[0]
    k = b[0]
    b /= b[0]
    z = roots(b)
    p = roots(a)
    return z, p, k


def zpk2tf(z, p, k):
    """
    Return polynomial transfer function representation from zeros and poles

    Parameters
    ----------
    z : array_like
        Zeros of the transfer function.
    p : array_like
        Poles of the transfer function.
    k : float
        System gain.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    """
    z = atleast_1d(z)
    k = atleast_1d(k)
    if len(z.shape) > 1:
        temp = poly(z[0])
        b = zeros((z.shape[0], z.shape[1] + 1), temp.dtype.char)
        if len(k) == 1:
            k = [k[0]] * z.shape[0]
        for i in range(z.shape[0]):
            b[i] = k[i] * poly(z[i])
    else:
        b = k * poly(z)
    a = atleast_1d(poly(p))

    # Use real output if possible.  Copied from numpy.poly, since
    # we can't depend on a specific version of numpy.
    if issubclass(b.dtype.type, numpy.complexfloating):
        # if complex roots are all complex conjugates, the roots are real.
        roots = numpy.asarray(z, complex)
        pos_roots = numpy.compress(roots.imag > 0, roots)
        neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
        if len(pos_roots) == len(neg_roots):
            if numpy.all(numpy.sort_complex(neg_roots) ==
                         numpy.sort_complex(pos_roots)):
                b = b.real.copy()

    if issubclass(a.dtype.type, numpy.complexfloating):
        # if complex roots are all complex conjugates, the roots are real.
        roots = numpy.asarray(p, complex)
        pos_roots = numpy.compress(roots.imag > 0, roots)
        neg_roots = numpy.conjugate(numpy.compress(roots.imag < 0, roots))
        if len(pos_roots) == len(neg_roots):
            if numpy.all(numpy.sort_complex(neg_roots) ==
                         numpy.sort_complex(pos_roots)):
                a = a.real.copy()

    return b, a


def tf2sos(b, a, pairing='nearest'):
    """
    Return second-order sections from transfer function representation

    Parameters
    ----------
    b : array_like
        Numerator polynomial coefficients.
    a : array_like
        Denominator polynomial coefficients.
    pairing : {'nearest', 'keep_odd'}, optional
        The method to use to combine pairs of poles and zeros into sections.
        See `zpk2sos`.

    Returns
    -------
    sos : ndarray
        Array of second-order filter coefficients, with shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    See Also
    --------
    zpk2sos, sosfilt

    Notes
    -----
    It is generally discouraged to convert from TF to SOS format, since doing
    so usually will not improve numerical precision errors. Instead, consider
    designing filters in ZPK format and converting directly to SOS. TF is
    converted to SOS by first converting to ZPK format, then converting
    ZPK to SOS.

    .. versionadded:: 0.16.0
    """
    return zpk2sos(*tf2zpk(b, a), pairing=pairing)


def sos2tf(sos):
    """
    Return a single transfer function from a series of second-order sections

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    Returns
    -------
    b : ndarray
        Numerator polynomial coefficients.
    a : ndarray
        Denominator polynomial coefficients.

    Notes
    -----
    .. versionadded:: 0.16.0
    """
    sos = np.asarray(sos)
    b = [1.]
    a = [1.]
    n_sections = sos.shape[0]
    for section in range(n_sections):
        b = np.polymul(b, sos[section, :3])
        a = np.polymul(a, sos[section, 3:])
    return b, a


def sos2zpk(sos):
    """
    Return zeros, poles, and gain of a series of second-order sections

    Parameters
    ----------
    sos : array_like
        Array of second-order filter coefficients, must have shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    Returns
    -------
    z : ndarray
        Zeros of the transfer function.
    p : ndarray
        Poles of the transfer function.
    k : float
        System gain.

    Notes
    -----
    .. versionadded:: 0.16.0
    """
    sos = np.asarray(sos)
    n_sections = sos.shape[0]
    z = np.empty(n_sections*2, np.complex128)
    p = np.empty(n_sections*2, np.complex128)
    k = 1.
    for section in range(n_sections):
        zpk = tf2zpk(sos[section, :3], sos[section, 3:])
        z[2*section:2*(section+1)] = zpk[0]
        p[2*section:2*(section+1)] = zpk[1]
        k *= zpk[2]
    return z, p, k


def _nearest_real_complex_idx(fro, to, which):
    """Get the next closest real or complex element based on distance"""
    assert which in ('real', 'complex')
    order = np.argsort(np.abs(fro - to))
    mask = np.isreal(fro[order])
    if which == 'complex':
        mask = ~mask
    return order[np.where(mask)[0][0]]


def zpk2sos(z, p, k, pairing='nearest'):
    """
    Return second-order sections from zeros, poles, and gain of a system

    Parameters
    ----------
    z : array_like
        Zeros of the transfer function.
    p : array_like
        Poles of the transfer function.
    k : float
        System gain.
    pairing : {'nearest', 'keep_odd'}, optional
        The method to use to combine pairs of poles and zeros into sections.
        See Notes below.

    Returns
    -------
    sos : ndarray
        Array of second-order filter coefficients, with shape
        ``(n_sections, 6)``. See `sosfilt` for the SOS filter format
        specification.

    See Also
    --------
    sosfilt

    Notes
    -----
    The algorithm used to convert ZPK to SOS format is designed to
    minimize errors due to numerical precision issues. The pairing
    algorithm attempts to minimize the peak gain of each biquadratic
    section. This is done by pairing poles with the nearest zeros, starting
    with the poles closest to the unit circle.

    *Algorithms*

    The current algorithms are designed specifically for use with digital
    filters. Although they can operate on analog filters, the results may
    be sub-optimal.

    The steps in the ``pairing='nearest'`` and ``pairing='keep_odd'``
    algorithms are mostly shared. The ``nearest`` algorithm attempts to
    minimize the peak gain, while ``'keep_odd'`` minimizes peak gain under
    the constraint that odd-order systems should retain one section
    as first order. The algorithm steps and are as follows:

    As a pre-processing step, add poles or zeros to the origin as
    necessary to obtain the same number of poles and zeros for pairing.
    If ``pairing == 'nearest'`` and there are an odd number of poles,
    add an additional pole and a zero at the origin.

    The following steps are then iterated over until no more poles or
    zeros remain:

    1. Take the (next remaining) pole (complex or real) closest to the
       unit circle to begin a new filter section.

    2. If the pole is real and there are no other remaining real poles [#]_,
       add the closest real zero to the section and leave it as a first
       order section. Note that after this step we are guaranteed to be
       left with an even number of real poles, complex poles, real zeros,
       and complex zeros for subsequent pairing iterations.

    3. Else:

        1. If the pole is complex and the zero is the only remaining real
           zero*, then pair the pole with the *next* closest zero
           (guaranteed to be complex). This is necessary to ensure that
           there will be a real zero remaining to eventually create a
           first-order section (thus keeping the odd order).

        2. Else pair the pole with the closest remaining zero (complex or
           real).

        3. Proceed to complete the second-order section by adding another
           pole and zero to the current pole and zero in the section:

            1. If the current pole and zero are both complex, add their
               conjugates.

            2. Else if the pole is complex and the zero is real, add the
               conjugate pole and the next closest real zero.

            3. Else if the pole is real and the zero is complex, add the
               conjugate zero and the real pole closest to those zeros.

            4. Else (we must have a real pole and real zero) add the next
               real pole closest to the unit circle, and then add the real
               zero closest to that pole.

    .. [#] This conditional can only be met for specific odd-order inputs
           with the ``pairing == 'keep_odd'`` method.

    .. versionadded:: 0.16.0

    Examples
    --------

    Design a 6th order low-pass elliptic digital filter for a system with a
    sampling rate of 8000 Hz that has a pass-band corner frequency of
    1000 Hz.  The ripple in the pass-band should not exceed 0.087 dB, and
    the attenuation in the stop-band should be at least 90 dB.

    In the following call to `signal.ellip`, we could use ``output='sos'``,
    but for this example, we'll use ``output='zpk'``, and then convert to SOS
    format with `zpk2sos`:

    >>> from scipy import signal
    >>> z, p, k = signal.ellip(6, 0.087, 90, 1000/(0.5*8000), output='zpk')

    Now convert to SOS format.

    >>> sos = signal.zpk2sos(z, p, k)

    The coefficents of the numerators of the sections:

    >>> sos[:, :3]
    array([[ 0.0014154 ,  0.00248707,  0.0014154 ],
           [ 1.        ,  0.72965193,  1.        ],
           [ 1.        ,  0.17594966,  1.        ]])

    The symmetry in the coefficients occurs because all the zeros are on the
    unit circle.

    The coefficients of the denominators of the sections:

    >>> sos[:, 3:]
    array([[ 1.        , -1.32543251,  0.46989499],
           [ 1.        , -1.26117915,  0.6262586 ],
           [ 1.        , -1.25707217,  0.86199667]])

    The next example shows the effect of the `pairing` option.  We have a
    system with three poles and three zeros, so the SOS array will have
    shape (2, 6).  The means there is, in effect, an extra pole and an extra
    zero at the origin in the SOS representation.

    >>> z1 = np.array([-1, -0.5-0.5j, -0.5+0.5j])
    >>> p1 = np.array([0.75, 0.8+0.1j, 0.8-0.1j])

    With ``pairing='nearest'`` (the default), we obtain

    >>> signal.zpk2sos(z1, p1, 1)
    array([[ 1.  ,  1.  ,  0.5 ,  1.  , -0.75,  0.  ],
           [ 1.  ,  1.  ,  0.  ,  1.  , -1.6 ,  0.65]])

    The first section has the zeros {-0.5-0.05j, -0.5+0.5j} and the poles
    {0, 0.75}, and the second section has the zeros {-1, 0} and poles
    {0.8+0.1j, 0.8-0.1j}.  Note that the extra pole and zero at the origin
    have been assigned to different sections.

    With ``pairing='keep_odd'``, we obtain:

    >>> signal.zpk2sos(z1, p1, 1, pairing='keep_odd')
    array([[ 1.  ,  1.  ,  0.  ,  1.  , -0.75,  0.  ],
           [ 1.  ,  1.  ,  0.5 ,  1.  , -1.6 ,  0.65]])

    The extra pole and zero at the origin are in the same section.
    The first section is, in effect, a first-order section.

    """
    # TODO in the near future:
    # 1. Add SOS capability to `filtfilt`, `freqz`, etc. somehow (#3259).
    # 2. Make `decimate` use `sosfilt` instead of `lfilter`.
    # 3. Make sosfilt automatically simplify sections to first order
    #    when possible. Note this might make `sosfiltfilt` a bit harder (ICs).
    # 4. Further optimizations of the section ordering / pole-zero pairing.
    # See the wiki for other potential issues.

    valid_pairings = ['nearest', 'keep_odd']
    if pairing not in valid_pairings:
        raise ValueError('pairing must be one of %s, not %s'
                         % (valid_pairings, pairing))
    if len(z) == len(p) == 0:
        return array([[k, 0., 0., 1., 0., 0.]])

    # ensure we have the same number of poles and zeros, and make copies
    p = np.concatenate((p, np.zeros(max(len(z) - len(p), 0))))
    z = np.concatenate((z, np.zeros(max(len(p) - len(z), 0))))
    n_sections = (max(len(p), len(z)) + 1) // 2
    sos = zeros((n_sections, 6))

    if len(p) % 2 == 1 and pairing == 'nearest':
        p = np.concatenate((p, [0.]))
        z = np.concatenate((z, [0.]))
    assert len(p) == len(z)

    # Ensure we have complex conjugate pairs
    # (note that _cplxreal only gives us one element of each complex pair):
    z = np.concatenate(_cplxreal(z))
    p = np.concatenate(_cplxreal(p))

    p_sos = np.zeros((n_sections, 2), np.complex128)
    z_sos = np.zeros_like(p_sos)
    for si in range(n_sections):
        # Select the next "worst" pole
        p1_idx = np.argmin(np.abs(1 - np.abs(p)))
        p1 = p[p1_idx]
        p = np.delete(p, p1_idx)

        # Pair that pole with a zero

        if np.isreal(p1) and np.isreal(p).sum() == 0:
            # Special case to set a first-order section
            z1_idx = _nearest_real_complex_idx(z, p1, 'real')
            z1 = z[z1_idx]
            z = np.delete(z, z1_idx)
            p2 = z2 = 0
        else:
            if not np.isreal(p1) and np.isreal(z).sum() == 1:
                # Special case to ensure we choose a complex zero to pair
                # with so later (setting up a first-order section)
                z1_idx = _nearest_real_complex_idx(z, p1, 'complex')
                assert not np.isreal(z[z1_idx])
            else:
                # Pair the pole with the closest zero (real or complex)
                z1_idx = np.argmin(np.abs(p1 - z))
            z1 = z[z1_idx]
            z = np.delete(z, z1_idx)

            # Now that we have p1 and z1, figure out what p2 and z2 need to be
            if not np.isreal(p1):
                if not np.isreal(z1):  # complex pole, complex zero
                    p2 = p1.conj()
                    z2 = z1.conj()
                else:  # complex pole, real zero
                    p2 = p1.conj()
                    z2_idx = _nearest_real_complex_idx(z, p1, 'real')
                    z2 = z[z2_idx]
                    assert np.isreal(z2)
                    z = np.delete(z, z2_idx)
            else:
                if not np.isreal(z1):  # real pole, complex zero
                    z2 = z1.conj()
                    p2_idx = _nearest_real_complex_idx(p, z1, 'real')
                    p2 = p[p2_idx]
                    assert np.isreal(p2)
                else:  # real pole, real zero
                    # pick the next "worst" pole to use
                    idx = np.where(np.isreal(p))[0]
                    assert len(idx) > 0
                    p2_idx = idx[np.argmin(np.abs(np.abs(p[idx]) - 1))]
                    p2 = p[p2_idx]
                    # find a real zero to match the added pole
                    assert np.isreal(p2)
                    z2_idx = _nearest_real_complex_idx(z, p2, 'real')
                    z2 = z[z2_idx]
                    assert np.isreal(z2)
                    z = np.delete(z, z2_idx)
                p = np.delete(p, p2_idx)
        p_sos[si] = [p1, p2]
        z_sos[si] = [z1, z2]
    assert len(p) == len(z) == 0  # we've consumed all poles and zeros
    del p, z

    # Construct the system, reversing order so the "worst" are last
    p_sos = np.reshape(p_sos[::-1], (n_sections, 2))
    z_sos = np.reshape(z_sos[::-1], (n_sections, 2))
    gains = np.ones(n_sections)
    gains[0] = k
    for si in range(n_sections):
        x = zpk2tf(z_sos[si], p_sos[si], gains[si])
        sos[si] = np.concatenate(x)
    return sos


def normalize(b, a):
    """Normalize polynomial representation of a transfer function.

    If values of `b` are too close to 0, they are removed. In that case, a
    BadCoefficients warning is emitted.

    """
    b, a = map(atleast_1d, (b, a))
    if len(a.shape) != 1:
        raise ValueError("Denominator polynomial must be rank-1 array.")
    if len(b.shape) > 2:
        raise ValueError("Numerator polynomial must be rank-1 or"
                         " rank-2 array.")
    if len(b.shape) == 1:
        b = asarray([b], b.dtype.char)
    while a[0] == 0.0 and len(a) > 1:
        a = a[1:]
    outb = b * (1.0) / a[0]
    outa = a * (1.0) / a[0]
    if allclose(0, outb[:, 0], atol=1e-14):
        warnings.warn("Badly conditioned filter coefficients (numerator): the "
                      "results may be meaningless", BadCoefficients)
        while allclose(0, outb[:, 0], atol=1e-14) and (outb.shape[-1] > 1):
            outb = outb[:, 1:]
    if outb.shape[0] == 1:
        outb = outb[0]
    return outb, outa


def lp2lp(b, a, wo=1.0):
    """
    Transform a lowpass filter prototype to a different frequency.

    Return an analog low-pass filter with cutoff frequency `wo`
    from an analog low-pass filter prototype with unity cutoff frequency, in
    transfer function ('ba') representation.

    """
    a, b = map(atleast_1d, (a, b))
    try:
        wo = float(wo)
    except TypeError:
        wo = float(wo[0])
    d = len(a)
    n = len(b)
    M = max((d, n))
    pwo = pow(wo, numpy.arange(M - 1, -1, -1))
    start1 = max((n - d, 0))
    start2 = max((d - n, 0))
    b = b * pwo[start1] / pwo[start2:]
    a = a * pwo[start1] / pwo[start1:]
    return normalize(b, a)


def lp2hp(b, a, wo=1.0):
    """
    Transform a lowpass filter prototype to a highpass filter.

    Return an analog high-pass filter with cutoff frequency `wo`
    from an analog low-pass filter prototype with unity cutoff frequency, in
    transfer function ('ba') representation.

    """
    a, b = map(atleast_1d, (a, b))
    try:
        wo = float(wo)
    except TypeError:
        wo = float(wo[0])
    d = len(a)
    n = len(b)
    if wo != 1:
        pwo = pow(wo, numpy.arange(max((d, n))))
    else:
        pwo = numpy.ones(max((d, n)), b.dtype.char)
    if d >= n:
        outa = a[::-1] * pwo
        outb = resize(b, (d,))
        outb[n:] = 0.0
        outb[:n] = b[::-1] * pwo[:n]
    else:
        outb = b[::-1] * pwo
        outa = resize(a, (n,))
        outa[d:] = 0.0
        outa[:d] = a[::-1] * pwo[:d]

    return normalize(outb, outa)


def lp2bp(b, a, wo=1.0, bw=1.0):
    """
    Transform a lowpass filter prototype to a bandpass filter.

    Return an analog band-pass filter with center frequency `wo` and
    bandwidth `bw` from an analog low-pass filter prototype with unity
    cutoff frequency, in transfer function ('ba') representation.

    """
    a, b = map(atleast_1d, (a, b))
    D = len(a) - 1
    N = len(b) - 1
    artype = mintypecode((a, b))
    ma = max([N, D])
    Np = N + ma
    Dp = D + ma
    bprime = numpy.zeros(Np + 1, artype)
    aprime = numpy.zeros(Dp + 1, artype)
    wosq = wo * wo
    for j in range(Np + 1):
        val = 0.0
        for i in range(0, N + 1):
            for k in range(0, i + 1):
                if ma - i + 2 * k == j:
                    val += comb(i, k) * b[N - i] * (wosq) ** (i - k) / bw ** i
        bprime[Np - j] = val
    for j in range(Dp + 1):
        val = 0.0
        for i in range(0, D + 1):
            for k in range(0, i + 1):
                if ma - i + 2 * k == j:
                    val += comb(i, k) * a[D - i] * (wosq) ** (i - k) / bw ** i
        aprime[Dp - j] = val

    return normalize(bprime, aprime)


def lp2bs(b, a, wo=1.0, bw=1.0):
    """
    Transform a lowpass filter prototype to a bandstop filter.

    Return an analog band-stop filter with center frequency `wo` and
    bandwidth `bw` from an analog low-pass filter prototype with unity
    cutoff frequency, in transfer function ('ba') representation.

    """
    a, b = map(atleast_1d, (a, b))
    D = len(a) - 1
    N = len(b) - 1
    artype = mintypecode((a, b))
    M = max([N, D])
    Np = M + M
    Dp = M + M
    bprime = numpy.zeros(Np + 1, artype)
    aprime = numpy.zeros(Dp + 1, artype)
    wosq = wo * wo
    for j in range(Np + 1):
        val = 0.0
        for i in range(0, N + 1):
            for k in range(0, M - i + 1):
                if i + 2 * k == j:
                    val += (comb(M - i, k) * b[N - i] *
                            (wosq) ** (M - i - k) * bw ** i)
        bprime[Np - j] = val
    for j in range(Dp + 1):
        val = 0.0
        for i in range(0, D + 1):
            for k in range(0, M - i + 1):
                if i + 2 * k == j:
                    val += (comb(M - i, k) * a[D - i] *
                            (wosq) ** (M - i - k) * bw ** i)
        aprime[Dp - j] = val

    return normalize(bprime, aprime)


def bilinear(b, a, fs=1.0):
    """Return a digital filter from an analog one using a bilinear transform.

    The bilinear transform substitutes ``(z-1) / (z+1)`` for ``s``.
    """
    fs = float(fs)
    a, b = map(atleast_1d, (a, b))
    D = len(a) - 1
    N = len(b) - 1
    artype = float
    M = max([N, D])
    Np = M
    Dp = M
    bprime = numpy.zeros(Np + 1, artype)
    aprime = numpy.zeros(Dp + 1, artype)
    for j in range(Np + 1):
        val = 0.0
        for i in range(N + 1):
            for k in range(i + 1):
                for l in range(M - i + 1):
                    if k + l == j:
                        val += (comb(i, k) * comb(M - i, l) * b[N - i] *
                                pow(2 * fs, i) * (-1) ** k)
        bprime[j] = real(val)
    for j in range(Dp + 1):
        val = 0.0
        for i in range(D + 1):
            for k in range(i + 1):
                for l in range(M - i + 1):
                    if k + l == j:
                        val += (comb(i, k) * comb(M - i, l) * a[D - i] *
                                pow(2 * fs, i) * (-1) ** k)
        aprime[j] = real(val)

    return normalize(bprime, aprime)


def iirdesign(wp, ws, gpass, gstop, analog=False, ftype='ellip', output='ba'):
    """Complete IIR digital and analog filter design.

    Given passband and stopband frequencies and gains, construct an analog or
    digital IIR filter of minimum order for a given basic type.  Return the
    output in numerator, denominator ('ba'), pole-zero ('zpk') or second order
    sections ('sos') form.

    Parameters
    ----------
    wp, ws : float
        Passband and stopband edge frequencies.
        For digital filters, these are normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`wp` and `ws` are thus in
        half-cycles / sample.)  For example:

            - Lowpass:   wp = 0.2,          ws = 0.3
            - Highpass:  wp = 0.3,          ws = 0.2
            - Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
            - Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]

        For analog filters, `wp` and `ws` are angular frequencies (e.g. rad/s).

    gpass : float
        The maximum loss in the passband (dB).
    gstop : float
        The minimum attenuation in the stopband (dB).
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    ftype : str, optional
        The type of IIR filter to design:

            - Butterworth   : 'butter'
            - Chebyshev I   : 'cheby1'
            - Chebyshev II  : 'cheby2'
            - Cauer/elliptic: 'ellip'
            - Bessel/Thomson: 'bessel'

    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba'.

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    See Also
    --------
    butter : Filter design using order and critical points
    cheby1, cheby2, ellip, bessel
    buttord : Find order and critical points from passband and stopband spec
    cheb1ord, cheb2ord, ellipord
    iirfilter : General filter design using order and critical frequencies

    Notes
    -----
    The ``'sos'`` output parameter was added in 0.16.0.
    """
    try:
        ordfunc = filter_dict[ftype][1]
    except KeyError:
        raise ValueError("Invalid IIR filter type: %s" % ftype)
    except IndexError:
        raise ValueError(("%s does not have order selection. Use "
                          "iirfilter function.") % ftype)

    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    band_type = 2 * (len(wp) - 1)
    band_type += 1
    if wp[0] >= ws[0]:
        band_type += 1

    btype = {1: 'lowpass', 2: 'highpass',
             3: 'bandstop', 4: 'bandpass'}[band_type]

    N, Wn = ordfunc(wp, ws, gpass, gstop, analog=analog)
    return iirfilter(N, Wn, rp=gpass, rs=gstop, analog=analog, btype=btype,
                     ftype=ftype, output=output)


def iirfilter(N, Wn, rp=None, rs=None, btype='band', analog=False,
              ftype='butter', output='ba'):
    """
    IIR digital and analog filter design given order and critical points.

    Design an Nth order digital or analog filter and return the filter
    coefficients.

    Parameters
    ----------
    N : int
        The order of the filter.
    Wn : array_like
        A scalar or length-2 sequence giving the critical frequencies.
        For digital filters, `Wn` is normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`Wn` is thus in
        half-cycles / sample.)
        For analog filters, `Wn` is an angular frequency (e.g. rad/s).
    rp : float, optional
        For Chebyshev and elliptic filters, provides the maximum ripple
        in the passband. (dB)
    rs : float, optional
        For Chebyshev and elliptic filters, provides the minimum attenuation
        in the stop band. (dB)
    btype : {'bandpass', 'lowpass', 'highpass', 'bandstop'}, optional
        The type of filter.  Default is 'bandpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    ftype : str, optional
        The type of IIR filter to design:

            - Butterworth   : 'butter'
            - Chebyshev I   : 'cheby1'
            - Chebyshev II  : 'cheby2'
            - Cauer/elliptic: 'ellip'
            - Bessel/Thomson: 'bessel'

    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba'.

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    See Also
    --------
    butter : Filter design using order and critical points
    cheby1, cheby2, ellip, bessel
    buttord : Find order and critical points from passband and stopband spec
    cheb1ord, cheb2ord, ellipord
    iirdesign : General filter design using passband and stopband spec

    Notes
    -----
    The ``'sos'`` output parameter was added in 0.16.0.

    Examples
    --------
    Generate a 17th-order Chebyshev II bandpass filter and plot the frequency
    response:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> b, a = signal.iirfilter(17, [50, 200], rs=60, btype='band',
    ...                         analog=True, ftype='cheby2')
    >>> w, h = signal.freqs(b, a, 1000)
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.semilogx(w, 20 * np.log10(abs(h)))
    >>> ax.set_title('Chebyshev Type II bandpass frequency response')
    >>> ax.set_xlabel('Frequency [radians / second]')
    >>> ax.set_ylabel('Amplitude [dB]')
    >>> ax.axis((10, 1000, -100, 10))
    >>> ax.grid(which='both', axis='both')
    >>> plt.show()

    """
    ftype, btype, output = [x.lower() for x in (ftype, btype, output)]
    Wn = asarray(Wn)
    try:
        btype = band_dict[btype]
    except KeyError:
        raise ValueError("'%s' is an invalid bandtype for filter." % btype)

    try:
        typefunc = filter_dict[ftype][0]
    except KeyError:
        raise ValueError("'%s' is not a valid basic IIR filter." % ftype)

    if output not in ['ba', 'zpk', 'sos']:
        raise ValueError("'%s' is not a valid output form." % output)

    if rp is not None and rp < 0:
        raise ValueError("passband ripple (rp) must be positive")

    if rs is not None and rs < 0:
        raise ValueError("stopband attenuation (rs) must be positive")

    # Get analog lowpass prototype
    if typefunc in [buttap, besselap]:
        z, p, k = typefunc(N)
    elif typefunc == cheb1ap:
        if rp is None:
            raise ValueError("passband ripple (rp) must be provided to "
                             "design a Chebyshev I filter.")
        z, p, k = typefunc(N, rp)
    elif typefunc == cheb2ap:
        if rs is None:
            raise ValueError("stopband attenuation (rs) must be provided to "
                             "design an Chebyshev II filter.")
        z, p, k = typefunc(N, rs)
    elif typefunc == ellipap:
        if rs is None or rp is None:
            raise ValueError("Both rp and rs must be provided to design an "
                             "elliptic filter.")
        z, p, k = typefunc(N, rp, rs)
    else:
        raise NotImplementedError("'%s' not implemented in iirfilter." % ftype)

    # Pre-warp frequencies for digital filter design
    if not analog:
        if numpy.any(Wn < 0) or numpy.any(Wn > 1):
            raise ValueError("Digital filter critical frequencies "
                             "must be 0 <= Wn <= 1")
        fs = 2.0
        warped = 2 * fs * tan(pi * Wn / fs)
    else:
        warped = Wn

    # transform to lowpass, bandpass, highpass, or bandstop
    if btype in ('lowpass', 'highpass'):
        if numpy.size(Wn) != 1:
            raise ValueError('Must specify a single critical frequency Wn')

        if btype == 'lowpass':
            z, p, k = _zpklp2lp(z, p, k, wo=warped)
        elif btype == 'highpass':
            z, p, k = _zpklp2hp(z, p, k, wo=warped)
    elif btype in ('bandpass', 'bandstop'):
        try:
            bw = warped[1] - warped[0]
            wo = sqrt(warped[0] * warped[1])
        except IndexError:
            raise ValueError('Wn must specify start and stop frequencies')

        if btype == 'bandpass':
            z, p, k = _zpklp2bp(z, p, k, wo=wo, bw=bw)
        elif btype == 'bandstop':
            z, p, k = _zpklp2bs(z, p, k, wo=wo, bw=bw)
    else:
        raise NotImplementedError("'%s' not implemented in iirfilter." % btype)

    # Find discrete equivalent if necessary
    if not analog:
        z, p, k = _zpkbilinear(z, p, k, fs=fs)

    # Transform to proper out type (pole-zero, state-space, numer-denom)
    if output == 'zpk':
        return z, p, k
    elif output == 'ba':
        return zpk2tf(z, p, k)
    elif output == 'sos':
        return zpk2sos(z, p, k)


def _relative_degree(z, p):
    """
    Return relative degree of transfer function from zeros and poles
    """
    degree = len(p) - len(z)
    if degree < 0:
        raise ValueError("Improper transfer function. "
                         "Must have at least as many poles as zeros.")
    else:
        return degree


# TODO: merge these into existing functions or make public versions

def _zpkbilinear(z, p, k, fs):
    """
    Return a digital filter from an analog one using a bilinear transform.

    Transform a set of poles and zeros from the analog s-plane to the digital
    z-plane using Tustin's method, which substitutes ``(z-1) / (z+1)`` for
    ``s``, maintaining the shape of the frequency response.

    Parameters
    ----------
    z : ndarray
        Zeros of the analog IIR filter transfer function.
    p : ndarray
        Poles of the analog IIR filter transfer function.
    k : float
        System gain of the analog IIR filter transfer function.
    fs : float
        Sample rate, as ordinary frequency (e.g. hertz). No prewarping is
        done in this function.

    Returns
    -------
    z : ndarray
        Zeros of the transformed digital filter transfer function.
    p : ndarray
        Poles of the transformed digital filter transfer function.
    k : float
        System gain of the transformed digital filter.

    """
    z = atleast_1d(z)
    p = atleast_1d(p)

    degree = _relative_degree(z, p)

    fs2 = 2*fs

    # Bilinear transform the poles and zeros
    z_z = (fs2 + z) / (fs2 - z)
    p_z = (fs2 + p) / (fs2 - p)

    # Any zeros that were at infinity get moved to the Nyquist frequency
    z_z = append(z_z, -ones(degree))

    # Compensate for gain change
    k_z = k * real(prod(fs2 - z) / prod(fs2 - p))

    return z_z, p_z, k_z


def _zpklp2lp(z, p, k, wo=1.0):
    """
    Transform a lowpass filter prototype to a different frequency.

    Return an analog low-pass filter with cutoff frequency `wo`
    from an analog low-pass filter prototype with unity cutoff frequency,
    using zeros, poles, and gain ('zpk') representation.

    Parameters
    ----------
    z : ndarray
        Zeros of the analog IIR filter transfer function.
    p : ndarray
        Poles of the analog IIR filter transfer function.
    k : float
        System gain of the analog IIR filter transfer function.
    wo : float
        Desired cutoff, as angular frequency (e.g. rad/s).
        Defaults to no change.

    Returns
    -------
    z : ndarray
        Zeros of the transformed low-pass filter transfer function.
    p : ndarray
        Poles of the transformed low-pass filter transfer function.
    k : float
        System gain of the transformed low-pass filter.

    Notes
    -----
    This is derived from the s-plane substitution

    .. math:: s \rightarrow \frac{s}{\omega_0}

    """
    z = atleast_1d(z)
    p = atleast_1d(p)
    wo = float(wo)  # Avoid int wraparound

    degree = _relative_degree(z, p)

    # Scale all points radially from origin to shift cutoff frequency
    z_lp = wo * z
    p_lp = wo * p

    # Each shifted pole decreases gain by wo, each shifted zero increases it.
    # Cancel out the net change to keep overall gain the same
    k_lp = k * wo**degree

    return z_lp, p_lp, k_lp


def _zpklp2hp(z, p, k, wo=1.0):
    """
    Transform a lowpass filter prototype to a highpass filter.

    Return an analog high-pass filter with cutoff frequency `wo`
    from an analog low-pass filter prototype with unity cutoff frequency,
    using zeros, poles, and gain ('zpk') representation.

    Parameters
    ----------
    z : ndarray
        Zeros of the analog IIR filter transfer function.
    p : ndarray
        Poles of the analog IIR filter transfer function.
    k : float
        System gain of the analog IIR filter transfer function.
    wo : float
        Desired cutoff, as angular frequency (e.g. rad/s).
        Defaults to no change.

    Returns
    -------
    z : ndarray
        Zeros of the transformed high-pass filter transfer function.
    p : ndarray
        Poles of the transformed high-pass filter transfer function.
    k : float
        System gain of the transformed high-pass filter.

    Notes
    -----
    This is derived from the s-plane substitution

    .. math:: s \rightarrow \frac{\omega_0}{s}

    This maintains symmetry of the lowpass and highpass responses on a
    logarithmic scale.

    """
    z = atleast_1d(z)
    p = atleast_1d(p)
    wo = float(wo)

    degree = _relative_degree(z, p)

    # Invert positions radially about unit circle to convert LPF to HPF
    # Scale all points radially from origin to shift cutoff frequency
    z_hp = wo / z
    p_hp = wo / p

    # If lowpass had zeros at infinity, inverting moves them to origin.
    z_hp = append(z_hp, zeros(degree))

    # Cancel out gain change caused by inversion
    k_hp = k * real(prod(-z) / prod(-p))

    return z_hp, p_hp, k_hp


def _zpklp2bp(z, p, k, wo=1.0, bw=1.0):
    """
    Transform a lowpass filter prototype to a bandpass filter.

    Return an analog band-pass filter with center frequency `wo` and
    bandwidth `bw` from an analog low-pass filter prototype with unity
    cutoff frequency, using zeros, poles, and gain ('zpk') representation.

    Parameters
    ----------
    z : ndarray
        Zeros of the analog IIR filter transfer function.
    p : ndarray
        Poles of the analog IIR filter transfer function.
    k : float
        System gain of the analog IIR filter transfer function.
    wo : float
        Desired passband center, as angular frequency (e.g. rad/s).
        Defaults to no change.
    bw : float
        Desired passband width, as angular frequency (e.g. rad/s).
        Defaults to 1.

    Returns
    -------
    z : ndarray
        Zeros of the transformed band-pass filter transfer function.
    p : ndarray
        Poles of the transformed band-pass filter transfer function.
    k : float
        System gain of the transformed band-pass filter.

    Notes
    -----
    This is derived from the s-plane substitution

    .. math:: s \rightarrow \frac{s^2 + {\omega_0}^2}{s \cdot \mathrm{BW}}

    This is the "wideband" transformation, producing a passband with
    geometric (log frequency) symmetry about `wo`.

    """
    z = atleast_1d(z)
    p = atleast_1d(p)
    wo = float(wo)
    bw = float(bw)

    degree = _relative_degree(z, p)

    # Scale poles and zeros to desired bandwidth
    z_lp = z * bw/2
    p_lp = p * bw/2

    # Square root needs to produce complex result, not NaN
    z_lp = z_lp.astype(complex)
    p_lp = p_lp.astype(complex)

    # Duplicate poles and zeros and shift from baseband to +wo and -wo
    z_bp = concatenate((z_lp + sqrt(z_lp**2 - wo**2),
                        z_lp - sqrt(z_lp**2 - wo**2)))
    p_bp = concatenate((p_lp + sqrt(p_lp**2 - wo**2),
                        p_lp - sqrt(p_lp**2 - wo**2)))

    # Move degree zeros to origin, leaving degree zeros at infinity for BPF
    z_bp = append(z_bp, zeros(degree))

    # Cancel out gain change from frequency scaling
    k_bp = k * bw**degree

    return z_bp, p_bp, k_bp


def _zpklp2bs(z, p, k, wo=1.0, bw=1.0):
    """
    Transform a lowpass filter prototype to a bandstop filter.

    Return an analog band-stop filter with center frequency `wo` and
    stopband width `bw` from an analog low-pass filter prototype with unity
    cutoff frequency, using zeros, poles, and gain ('zpk') representation.

    Parameters
    ----------
    z : ndarray
        Zeros of the analog IIR filter transfer function.
    p : ndarray
        Poles of the analog IIR filter transfer function.
    k : float
        System gain of the analog IIR filter transfer function.
    wo : float
        Desired stopband center, as angular frequency (e.g. rad/s).
        Defaults to no change.
    bw : float
        Desired stopband width, as angular frequency (e.g. rad/s).
        Defaults to 1.

    Returns
    -------
    z : ndarray
        Zeros of the transformed band-stop filter transfer function.
    p : ndarray
        Poles of the transformed band-stop filter transfer function.
    k : float
        System gain of the transformed band-stop filter.

    Notes
    -----
    This is derived from the s-plane substitution

    .. math:: s \rightarrow \frac{s \cdot \mathrm{BW}}{s^2 + {\omega_0}^2}

    This is the "wideband" transformation, producing a stopband with
    geometric (log frequency) symmetry about `wo`.

    """
    z = atleast_1d(z)
    p = atleast_1d(p)
    wo = float(wo)
    bw = float(bw)

    degree = _relative_degree(z, p)

    # Invert to a highpass filter with desired bandwidth
    z_hp = (bw/2) / z
    p_hp = (bw/2) / p

    # Square root needs to produce complex result, not NaN
    z_hp = z_hp.astype(complex)
    p_hp = p_hp.astype(complex)

    # Duplicate poles and zeros and shift from baseband to +wo and -wo
    z_bs = concatenate((z_hp + sqrt(z_hp**2 - wo**2),
                        z_hp - sqrt(z_hp**2 - wo**2)))
    p_bs = concatenate((p_hp + sqrt(p_hp**2 - wo**2),
                        p_hp - sqrt(p_hp**2 - wo**2)))

    # Move any zeros that were at infinity to the center of the stopband
    z_bs = append(z_bs, +1j*wo * ones(degree))
    z_bs = append(z_bs, -1j*wo * ones(degree))

    # Cancel out gain change caused by inversion
    k_bs = k * real(prod(-z) / prod(-p))

    return z_bs, p_bs, k_bs


def butter(N, Wn, btype='low', analog=False, output='ba'):
    """
    Butterworth digital and analog filter design.

    Design an Nth order digital or analog Butterworth filter and return
    the filter coefficients.

    Parameters
    ----------
    N : int
        The order of the filter.
    Wn : array_like
        A scalar or length-2 sequence giving the critical frequencies.
        For a Butterworth filter, this is the point at which the gain
        drops to 1/sqrt(2) that of the passband (the "-3 dB point").
        For digital filters, `Wn` is normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`Wn` is thus in
        half-cycles / sample.)
        For analog filters, `Wn` is an angular frequency (e.g. rad/s).
    btype : {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
        The type of filter.  Default is 'lowpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba'.

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    See Also
    --------
    buttord

    Notes
    -----
    The Butterworth filter has maximally flat frequency response in the
    passband.

    The ``'sos'`` output parameter was added in 0.16.0.

    Examples
    --------
    Plot the filter's frequency response, showing the critical points:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> b, a = signal.butter(4, 100, 'low', analog=True)
    >>> w, h = signal.freqs(b, a)
    >>> plt.semilogx(w, 20 * np.log10(abs(h)))
    >>> plt.title('Butterworth filter frequency response')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.margins(0, 0.1)
    >>> plt.grid(which='both', axis='both')
    >>> plt.axvline(100, color='green') # cutoff frequency
    >>> plt.show()

    """
    return iirfilter(N, Wn, btype=btype, analog=analog,
                     output=output, ftype='butter')


def cheby1(N, rp, Wn, btype='low', analog=False, output='ba'):
    """
    Chebyshev type I digital and analog filter design.

    Design an Nth order digital or analog Chebyshev type I filter and
    return the filter coefficients.

    Parameters
    ----------
    N : int
        The order of the filter.
    rp : float
        The maximum ripple allowed below unity gain in the passband.
        Specified in decibels, as a positive number.
    Wn : array_like
        A scalar or length-2 sequence giving the critical frequencies.
        For Type I filters, this is the point in the transition band at which
        the gain first drops below -`rp`.
        For digital filters, `Wn` is normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`Wn` is thus in
        half-cycles / sample.)
        For analog filters, `Wn` is an angular frequency (e.g. rad/s).
    btype : {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
        The type of filter.  Default is 'lowpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba'.

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    See Also
    --------
    cheb1ord

    Notes
    -----
    The Chebyshev type I filter maximizes the rate of cutoff between the
    frequency response's passband and stopband, at the expense of ripple in
    the passband and increased ringing in the step response.

    Type I filters roll off faster than Type II (`cheby2`), but Type II
    filters do not have any ripple in the passband.

    The equiripple passband has N maxima or minima (for example, a
    5th-order filter has 3 maxima and 2 minima).  Consequently, the DC gain is
    unity for odd-order filters, or -rp dB for even-order filters.

    The ``'sos'`` output parameter was added in 0.16.0.

    Examples
    --------
    Plot the filter's frequency response, showing the critical points:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> b, a = signal.cheby1(4, 5, 100, 'low', analog=True)
    >>> w, h = signal.freqs(b, a)
    >>> plt.semilogx(w, 20 * np.log10(abs(h)))
    >>> plt.title('Chebyshev Type I frequency response (rp=5)')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.margins(0, 0.1)
    >>> plt.grid(which='both', axis='both')
    >>> plt.axvline(100, color='green') # cutoff frequency
    >>> plt.axhline(-5, color='green') # rp
    >>> plt.show()

    """
    return iirfilter(N, Wn, rp=rp, btype=btype, analog=analog,
                     output=output, ftype='cheby1')


def cheby2(N, rs, Wn, btype='low', analog=False, output='ba'):
    """
    Chebyshev type II digital and analog filter design.

    Design an Nth order digital or analog Chebyshev type II filter and
    return the filter coefficients.

    Parameters
    ----------
    N : int
        The order of the filter.
    rs : float
        The minimum attenuation required in the stop band.
        Specified in decibels, as a positive number.
    Wn : array_like
        A scalar or length-2 sequence giving the critical frequencies.
        For Type II filters, this is the point in the transition band at which
        the gain first reaches -`rs`.
        For digital filters, `Wn` is normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`Wn` is thus in
        half-cycles / sample.)
        For analog filters, `Wn` is an angular frequency (e.g. rad/s).
    btype : {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
        The type of filter.  Default is 'lowpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba'.

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    See Also
    --------
    cheb2ord

    Notes
    -----
    The Chebyshev type II filter maximizes the rate of cutoff between the
    frequency response's passband and stopband, at the expense of ripple in
    the stopband and increased ringing in the step response.

    Type II filters do not roll off as fast as Type I (`cheby1`).

    The ``'sos'`` output parameter was added in 0.16.0.

    Examples
    --------
    Plot the filter's frequency response, showing the critical points:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> b, a = signal.cheby2(4, 40, 100, 'low', analog=True)
    >>> w, h = signal.freqs(b, a)
    >>> plt.semilogx(w, 20 * np.log10(abs(h)))
    >>> plt.title('Chebyshev Type II frequency response (rs=40)')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.margins(0, 0.1)
    >>> plt.grid(which='both', axis='both')
    >>> plt.axvline(100, color='green') # cutoff frequency
    >>> plt.axhline(-40, color='green') # rs
    >>> plt.show()

    """
    return iirfilter(N, Wn, rs=rs, btype=btype, analog=analog,
                     output=output, ftype='cheby2')


def ellip(N, rp, rs, Wn, btype='low', analog=False, output='ba'):
    """
    Elliptic (Cauer) digital and analog filter design.

    Design an Nth order digital or analog elliptic filter and return
    the filter coefficients.

    Parameters
    ----------
    N : int
        The order of the filter.
    rp : float
        The maximum ripple allowed below unity gain in the passband.
        Specified in decibels, as a positive number.
    rs : float
        The minimum attenuation required in the stop band.
        Specified in decibels, as a positive number.
    Wn : array_like
        A scalar or length-2 sequence giving the critical frequencies.
        For elliptic filters, this is the point in the transition band at
        which the gain first drops below -`rp`.
        For digital filters, `Wn` is normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`Wn` is thus in
        half-cycles / sample.)
        For analog filters, `Wn` is an angular frequency (e.g. rad/s).
    btype : {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
        The type of filter.  Default is 'lowpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba'.

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    See Also
    --------
    ellipord

    Notes
    -----
    Also known as Cauer or Zolotarev filters, the elliptical filter maximizes
    the rate of transition between the frequency response's passband and
    stopband, at the expense of ripple in both, and increased ringing in the
    step response.

    As `rp` approaches 0, the elliptical filter becomes a Chebyshev
    type II filter (`cheby2`).  As `rs` approaches 0, it becomes a Chebyshev
    type I filter (`cheby1`).  As both approach 0, it becomes a Butterworth
    filter (`butter`).

    The equiripple passband has N maxima or minima (for example, a
    5th-order filter has 3 maxima and 2 minima).  Consequently, the DC gain is
    unity for odd-order filters, or -rp dB for even-order filters.

    The ``'sos'`` output parameter was added in 0.16.0.

    Examples
    --------
    Plot the filter's frequency response, showing the critical points:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> b, a = signal.ellip(4, 5, 40, 100, 'low', analog=True)
    >>> w, h = signal.freqs(b, a)
    >>> plt.semilogx(w, 20 * np.log10(abs(h)))
    >>> plt.title('Elliptic filter frequency response (rp=5, rs=40)')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.margins(0, 0.1)
    >>> plt.grid(which='both', axis='both')
    >>> plt.axvline(100, color='green') # cutoff frequency
    >>> plt.axhline(-40, color='green') # rs
    >>> plt.axhline(-5, color='green') # rp
    >>> plt.show()

    """
    return iirfilter(N, Wn, rs=rs, rp=rp, btype=btype, analog=analog,
                     output=output, ftype='elliptic')


def bessel(N, Wn, btype='low', analog=False, output='ba'):
    """Bessel/Thomson digital and analog filter design.

    Design an Nth order digital or analog Bessel filter and return the
    filter coefficients.

    Parameters
    ----------
    N : int
        The order of the filter.
    Wn : array_like
        A scalar or length-2 sequence giving the critical frequencies.
        For a Bessel filter, this is defined as the point at which the
        asymptotes of the response are the same as a Butterworth filter of
        the same order.
        For digital filters, `Wn` is normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`Wn` is thus in
        half-cycles / sample.)
        For analog filters, `Wn` is an angular frequency (e.g. rad/s).
    btype : {'lowpass', 'highpass', 'bandpass', 'bandstop'}, optional
        The type of filter.  Default is 'lowpass'.
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.
    output : {'ba', 'zpk', 'sos'}, optional
        Type of output:  numerator/denominator ('ba'), pole-zero ('zpk'), or
        second-order sections ('sos'). Default is 'ba'.

    Returns
    -------
    b, a : ndarray, ndarray
        Numerator (`b`) and denominator (`a`) polynomials of the IIR filter.
        Only returned if ``output='ba'``.
    z, p, k : ndarray, ndarray, float
        Zeros, poles, and system gain of the IIR filter transfer
        function.  Only returned if ``output='zpk'``.
    sos : ndarray
        Second-order sections representation of the IIR filter.
        Only returned if ``output=='sos'``.

    Notes
    -----
    Also known as a Thomson filter, the analog Bessel filter has maximally
    flat group delay and maximally linear phase response, with very little
    ringing in the step response.

    As order increases, the Bessel filter approaches a Gaussian filter.

    The digital Bessel filter is generated using the bilinear
    transform, which does not preserve the phase response of the analog
    filter. As such, it is only approximately correct at frequencies
    below about fs/4.  To get maximally flat group delay at higher
    frequencies, the analog Bessel filter must be transformed using
    phase-preserving techniques.

    For a given `Wn`, the lowpass and highpass filter have the same phase vs
    frequency curves; they are "phase-matched".

    The ``'sos'`` output parameter was added in 0.16.0.

    Examples
    --------
    Plot the filter's frequency response, showing the flat group delay and
    the relationship to the Butterworth's cutoff frequency:

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> b, a = signal.butter(4, 100, 'low', analog=True)
    >>> w, h = signal.freqs(b, a)
    >>> plt.plot(w, 20 * np.log10(np.abs(h)), color='silver', ls='dashed')
    >>> b, a = signal.bessel(4, 100, 'low', analog=True)
    >>> w, h = signal.freqs(b, a)
    >>> plt.semilogx(w, 20 * np.log10(np.abs(h)))
    >>> plt.title('Bessel filter frequency response (with Butterworth)')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.margins(0, 0.1)
    >>> plt.grid(which='both', axis='both')
    >>> plt.axvline(100, color='green') # cutoff frequency
    >>> plt.show()

    >>> plt.figure()
    >>> plt.semilogx(w[1:], -np.diff(np.unwrap(np.angle(h)))/np.diff(w))
    >>> plt.title('Bessel filter group delay')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Group delay [seconds]')
    >>> plt.margins(0, 0.1)
    >>> plt.grid(which='both', axis='both')
    >>> plt.show()

    """
    return iirfilter(N, Wn, btype=btype, analog=analog,
                     output=output, ftype='bessel')


def maxflat():
    pass


def yulewalk():
    pass


def band_stop_obj(wp, ind, passb, stopb, gpass, gstop, type):
    """
    Band Stop Objective Function for order minimization.

    Returns the non-integer order for an analog band stop filter.

    Parameters
    ----------
    wp : scalar
        Edge of passband `passb`.
    ind : int, {0, 1}
        Index specifying which `passb` edge to vary (0 or 1).
    passb : ndarray
        Two element sequence of fixed passband edges.
    stopb : ndarray
        Two element sequence of fixed stopband edges.
    gstop : float
        Amount of attenuation in stopband in dB.
    gpass : float
        Amount of ripple in the passband in dB.
    type : {'butter', 'cheby', 'ellip'}
        Type of filter.

    Returns
    -------
    n : scalar
        Filter order (possibly non-integer).

    """
    passbC = passb.copy()
    passbC[ind] = wp
    nat = (stopb * (passbC[0] - passbC[1]) /
           (stopb ** 2 - passbC[0] * passbC[1]))
    nat = min(abs(nat))

    if type == 'butter':
        GSTOP = 10 ** (0.1 * abs(gstop))
        GPASS = 10 ** (0.1 * abs(gpass))
        n = (log10((GSTOP - 1.0) / (GPASS - 1.0)) / (2 * log10(nat)))
    elif type == 'cheby':
        GSTOP = 10 ** (0.1 * abs(gstop))
        GPASS = 10 ** (0.1 * abs(gpass))
        n = arccosh(sqrt((GSTOP - 1.0) / (GPASS - 1.0))) / arccosh(nat)
    elif type == 'ellip':
        GSTOP = 10 ** (0.1 * gstop)
        GPASS = 10 ** (0.1 * gpass)
        arg1 = sqrt((GPASS - 1.0) / (GSTOP - 1.0))
        arg0 = 1.0 / nat
        d0 = special.ellipk([arg0 ** 2, 1 - arg0 ** 2])
        d1 = special.ellipk([arg1 ** 2, 1 - arg1 ** 2])
        n = (d0[0] * d1[1] / (d0[1] * d1[0]))
    else:
        raise ValueError("Incorrect type: %s" % type)
    return n


def buttord(wp, ws, gpass, gstop, analog=False):
    """Butterworth filter order selection.

    Return the order of the lowest order digital or analog Butterworth filter
    that loses no more than `gpass` dB in the passband and has at least
    `gstop` dB attenuation in the stopband.

    Parameters
    ----------
    wp, ws : float
        Passband and stopband edge frequencies.
        For digital filters, these are normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`wp` and `ws` are thus in
        half-cycles / sample.)  For example:

            - Lowpass:   wp = 0.2,          ws = 0.3
            - Highpass:  wp = 0.3,          ws = 0.2
            - Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
            - Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]

        For analog filters, `wp` and `ws` are angular frequencies (e.g. rad/s).

    gpass : float
        The maximum loss in the passband (dB).
    gstop : float
        The minimum attenuation in the stopband (dB).
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.

    Returns
    -------
    ord : int
        The lowest order for a Butterworth filter which meets specs.
    wn : ndarray or float
        The Butterworth natural frequency (i.e. the "3dB frequency").  Should
        be used with `butter` to give filter results.

    See Also
    --------
    butter : Filter design using order and critical points
    cheb1ord : Find order and critical points from passband and stopband spec
    cheb2ord, ellipord
    iirfilter : General filter design using order and critical frequencies
    iirdesign : General filter design using passband and stopband spec

    Examples
    --------
    Design an analog bandpass filter with passband within 3 dB from 20 to
    50 rad/s, while rejecting at least -40 dB below 14 and above 60 rad/s.
    Plot its frequency response, showing the passband and stopband
    constraints in gray.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> N, Wn = signal.buttord([20, 50], [14, 60], 3, 40, True)
    >>> b, a = signal.butter(N, Wn, 'band', True)
    >>> w, h = signal.freqs(b, a, np.logspace(1, 2, 500))
    >>> plt.semilogx(w, 20 * np.log10(abs(h)))
    >>> plt.title('Butterworth bandpass filter fit to constraints')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.grid(which='both', axis='both')
    >>> plt.fill([1,  14,  14,   1], [-40, -40, 99, 99], '0.9', lw=0) # stop
    >>> plt.fill([20, 20,  50,  50], [-99, -3, -3, -99], '0.9', lw=0) # pass
    >>> plt.fill([60, 60, 1e9, 1e9], [99, -40, -40, 99], '0.9', lw=0) # stop
    >>> plt.axis([10, 100, -60, 3])
    >>> plt.show()

    """
    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2 * (len(wp) - 1)
    filter_type += 1
    if wp[0] >= ws[0]:
        filter_type += 1

    # Pre-warp frequencies for digital filter design
    if not analog:
        passb = tan(pi * wp / 2.0)
        stopb = tan(pi * ws / 2.0)
    else:
        passb = wp * 1.0
        stopb = ws * 1.0

    if filter_type == 1:            # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:          # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0] - 1e-12,
                                 args=(0, passb, stopb, gpass, gstop,
                                       'butter'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1] + 1e-12, passb[1],
                                 args=(1, passb, stopb, gpass, gstop,
                                       'butter'),
                                 disp=0)
        passb[1] = wp1
        nat = ((stopb * (passb[0] - passb[1])) /
               (stopb ** 2 - passb[0] * passb[1]))
    elif filter_type == 4:          # pass
        nat = ((stopb ** 2 - passb[0] * passb[1]) /
               (stopb * (passb[0] - passb[1])))

    nat = min(abs(nat))

    GSTOP = 10 ** (0.1 * abs(gstop))
    GPASS = 10 ** (0.1 * abs(gpass))
    ord = int(ceil(log10((GSTOP - 1.0) / (GPASS - 1.0)) / (2 * log10(nat))))

    # Find the Butterworth natural frequency WN (or the "3dB" frequency")
    # to give exactly gpass at passb.
    try:
        W0 = (GPASS - 1.0) ** (-1.0 / (2.0 * ord))
    except ZeroDivisionError:
        W0 = 1.0
        print("Warning, order is zero...check input parameters.")

    # now convert this frequency back from lowpass prototype
    # to the original analog filter

    if filter_type == 1:  # low
        WN = W0 * passb
    elif filter_type == 2:  # high
        WN = passb / W0
    elif filter_type == 3:  # stop
        WN = numpy.zeros(2, float)
        discr = sqrt((passb[1] - passb[0]) ** 2 +
                     4 * W0 ** 2 * passb[0] * passb[1])
        WN[0] = ((passb[1] - passb[0]) + discr) / (2 * W0)
        WN[1] = ((passb[1] - passb[0]) - discr) / (2 * W0)
        WN = numpy.sort(abs(WN))
    elif filter_type == 4:  # pass
        W0 = numpy.array([-W0, W0], float)
        WN = (-W0 * (passb[1] - passb[0]) / 2.0 +
              sqrt(W0 ** 2 / 4.0 * (passb[1] - passb[0]) ** 2 +
                   passb[0] * passb[1]))
        WN = numpy.sort(abs(WN))
    else:
        raise ValueError("Bad type: %s" % filter_type)

    if not analog:
        wn = (2.0 / pi) * arctan(WN)
    else:
        wn = WN

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn


def cheb1ord(wp, ws, gpass, gstop, analog=False):
    """Chebyshev type I filter order selection.

    Return the order of the lowest order digital or analog Chebyshev Type I
    filter that loses no more than `gpass` dB in the passband and has at
    least `gstop` dB attenuation in the stopband.

    Parameters
    ----------
    wp, ws : float
        Passband and stopband edge frequencies.
        For digital filters, these are normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`wp` and `ws` are thus in
        half-cycles / sample.)  For example:

            - Lowpass:   wp = 0.2,          ws = 0.3
            - Highpass:  wp = 0.3,          ws = 0.2
            - Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
            - Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]

        For analog filters, `wp` and `ws` are angular frequencies (e.g. rad/s).

    gpass : float
        The maximum loss in the passband (dB).
    gstop : float
        The minimum attenuation in the stopband (dB).
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.

    Returns
    -------
    ord : int
        The lowest order for a Chebyshev type I filter that meets specs.
    wn : ndarray or float
        The Chebyshev natural frequency (the "3dB frequency") for use with
        `cheby1` to give filter results.

    See Also
    --------
    cheby1 : Filter design using order and critical points
    buttord : Find order and critical points from passband and stopband spec
    cheb2ord, ellipord
    iirfilter : General filter design using order and critical frequencies
    iirdesign : General filter design using passband and stopband spec

    Examples
    --------
    Design a digital lowpass filter such that the passband is within 3 dB up
    to 0.2*(fs/2), while rejecting at least -40 dB above 0.3*(fs/2).  Plot its
    frequency response, showing the passband and stopband constraints in gray.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> N, Wn = signal.cheb1ord(0.2, 0.3, 3, 40)
    >>> b, a = signal.cheby1(N, 3, Wn, 'low')
    >>> w, h = signal.freqz(b, a)
    >>> plt.semilogx(w / np.pi, 20 * np.log10(abs(h)))
    >>> plt.title('Chebyshev I lowpass filter fit to constraints')
    >>> plt.xlabel('Normalized frequency')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.grid(which='both', axis='both')
    >>> plt.fill([.01, 0.2, 0.2, .01], [-3, -3, -99, -99], '0.9', lw=0) # stop
    >>> plt.fill([0.3, 0.3,   2,   2], [ 9, -40, -40,  9], '0.9', lw=0) # pass
    >>> plt.axis([0.08, 1, -60, 3])
    >>> plt.show()

    """
    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2 * (len(wp) - 1)
    if wp[0] < ws[0]:
        filter_type += 1
    else:
        filter_type += 2

    # Pre-warp frequencies for digital filter design
    if not analog:
        passb = tan(pi * wp / 2.0)
        stopb = tan(pi * ws / 2.0)
    else:
        passb = wp * 1.0
        stopb = ws * 1.0

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0] - 1e-12,
                                 args=(0, passb, stopb, gpass, gstop, 'cheby'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1] + 1e-12, passb[1],
                                 args=(1, passb, stopb, gpass, gstop, 'cheby'),
                                 disp=0)
        passb[1] = wp1
        nat = ((stopb * (passb[0] - passb[1])) /
               (stopb ** 2 - passb[0] * passb[1]))
    elif filter_type == 4:  # pass
        nat = ((stopb ** 2 - passb[0] * passb[1]) /
               (stopb * (passb[0] - passb[1])))

    nat = min(abs(nat))

    GSTOP = 10 ** (0.1 * abs(gstop))
    GPASS = 10 ** (0.1 * abs(gpass))
    ord = int(ceil(arccosh(sqrt((GSTOP - 1.0) / (GPASS - 1.0))) /
                   arccosh(nat)))

    # Natural frequencies are just the passband edges
    if not analog:
        wn = (2.0 / pi) * arctan(passb)
    else:
        wn = passb

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn


def cheb2ord(wp, ws, gpass, gstop, analog=False):
    """Chebyshev type II filter order selection.

    Return the order of the lowest order digital or analog Chebyshev Type II
    filter that loses no more than `gpass` dB in the passband and has at least
    `gstop` dB attenuation in the stopband.

    Parameters
    ----------
    wp, ws : float
        Passband and stopband edge frequencies.
        For digital filters, these are normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`wp` and `ws` are thus in
        half-cycles / sample.)  For example:

            - Lowpass:   wp = 0.2,          ws = 0.3
            - Highpass:  wp = 0.3,          ws = 0.2
            - Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
            - Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]

        For analog filters, `wp` and `ws` are angular frequencies (e.g. rad/s).

    gpass : float
        The maximum loss in the passband (dB).
    gstop : float
        The minimum attenuation in the stopband (dB).
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.

    Returns
    -------
    ord : int
        The lowest order for a Chebyshev type II filter that meets specs.
    wn : ndarray or float
        The Chebyshev natural frequency (the "3dB frequency") for use with
        `cheby2` to give filter results.

    See Also
    --------
    cheby2 : Filter design using order and critical points
    buttord : Find order and critical points from passband and stopband spec
    cheb1ord, ellipord
    iirfilter : General filter design using order and critical frequencies
    iirdesign : General filter design using passband and stopband spec

    Examples
    --------
    Design a digital bandstop filter which rejects -60 dB from 0.2*(fs/2) to
    0.5*(fs/2), while staying within 3 dB below 0.1*(fs/2) or above
    0.6*(fs/2).  Plot its frequency response, showing the passband and
    stopband constraints in gray.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> N, Wn = signal.cheb2ord([0.1, 0.6], [0.2, 0.5], 3, 60)
    >>> b, a = signal.cheby2(N, 60, Wn, 'stop')
    >>> w, h = signal.freqz(b, a)
    >>> plt.semilogx(w / np.pi, 20 * np.log10(abs(h)))
    >>> plt.title('Chebyshev II bandstop filter fit to constraints')
    >>> plt.xlabel('Normalized frequency')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.grid(which='both', axis='both')
    >>> plt.fill([.01, .1, .1, .01], [-3,  -3, -99, -99], '0.9', lw=0) # stop
    >>> plt.fill([.2,  .2, .5,  .5], [ 9, -60, -60,   9], '0.9', lw=0) # pass
    >>> plt.fill([.6,  .6,  2,   2], [-99, -3,  -3, -99], '0.9', lw=0) # stop
    >>> plt.axis([0.06, 1, -80, 3])
    >>> plt.show()

    """
    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2 * (len(wp) - 1)
    if wp[0] < ws[0]:
        filter_type += 1
    else:
        filter_type += 2

    # Pre-warp frequencies for digital filter design
    if not analog:
        passb = tan(pi * wp / 2.0)
        stopb = tan(pi * ws / 2.0)
    else:
        passb = wp * 1.0
        stopb = ws * 1.0

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0] - 1e-12,
                                 args=(0, passb, stopb, gpass, gstop, 'cheby'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1] + 1e-12, passb[1],
                                 args=(1, passb, stopb, gpass, gstop, 'cheby'),
                                 disp=0)
        passb[1] = wp1
        nat = ((stopb * (passb[0] - passb[1])) /
               (stopb ** 2 - passb[0] * passb[1]))
    elif filter_type == 4:  # pass
        nat = ((stopb ** 2 - passb[0] * passb[1]) /
               (stopb * (passb[0] - passb[1])))

    nat = min(abs(nat))

    GSTOP = 10 ** (0.1 * abs(gstop))
    GPASS = 10 ** (0.1 * abs(gpass))
    ord = int(ceil(arccosh(sqrt((GSTOP - 1.0) / (GPASS - 1.0))) /
                   arccosh(nat)))

    # Find frequency where analog response is -gpass dB.
    # Then convert back from low-pass prototype to the original filter.

    new_freq = cosh(1.0 / ord * arccosh(sqrt((GSTOP - 1.0) / (GPASS - 1.0))))
    new_freq = 1.0 / new_freq

    if filter_type == 1:
        nat = passb / new_freq
    elif filter_type == 2:
        nat = passb * new_freq
    elif filter_type == 3:
        nat = numpy.zeros(2, float)
        nat[0] = (new_freq / 2.0 * (passb[0] - passb[1]) +
                  sqrt(new_freq ** 2 * (passb[1] - passb[0]) ** 2 / 4.0 +
                       passb[1] * passb[0]))
        nat[1] = passb[1] * passb[0] / nat[0]
    elif filter_type == 4:
        nat = numpy.zeros(2, float)
        nat[0] = (1.0 / (2.0 * new_freq) * (passb[0] - passb[1]) +
                  sqrt((passb[1] - passb[0]) ** 2 / (4.0 * new_freq ** 2) +
                       passb[1] * passb[0]))
        nat[1] = passb[0] * passb[1] / nat[0]

    if not analog:
        wn = (2.0 / pi) * arctan(nat)
    else:
        wn = nat

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn


def ellipord(wp, ws, gpass, gstop, analog=False):
    """Elliptic (Cauer) filter order selection.

    Return the order of the lowest order digital or analog elliptic filter
    that loses no more than `gpass` dB in the passband and has at least
    `gstop` dB attenuation in the stopband.

    Parameters
    ----------
    wp, ws : float
        Passband and stopband edge frequencies.
        For digital filters, these are normalized from 0 to 1, where 1 is the
        Nyquist frequency, pi radians/sample.  (`wp` and `ws` are thus in
        half-cycles / sample.)  For example:

            - Lowpass:   wp = 0.2,          ws = 0.3
            - Highpass:  wp = 0.3,          ws = 0.2
            - Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
            - Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]

        For analog filters, `wp` and `ws` are angular frequencies (e.g. rad/s).

    gpass : float
        The maximum loss in the passband (dB).
    gstop : float
        The minimum attenuation in the stopband (dB).
    analog : bool, optional
        When True, return an analog filter, otherwise a digital filter is
        returned.

    Returns
    -------
    ord : int
        The lowest order for an Elliptic (Cauer) filter that meets specs.
    wn : ndarray or float
        The Chebyshev natural frequency (the "3dB frequency") for use with
        `ellip` to give filter results.

    See Also
    --------
    ellip : Filter design using order and critical points
    buttord : Find order and critical points from passband and stopband spec
    cheb1ord, cheb2ord
    iirfilter : General filter design using order and critical frequencies
    iirdesign : General filter design using passband and stopband spec

    Examples
    --------
    Design an analog highpass filter such that the passband is within 3 dB
    above 30 rad/s, while rejecting -60 dB at 10 rad/s.  Plot its
    frequency response, showing the passband and stopband constraints in gray.

    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> N, Wn = signal.ellipord(30, 10, 3, 60, True)
    >>> b, a = signal.ellip(N, 3, 60, Wn, 'high', True)
    >>> w, h = signal.freqs(b, a, np.logspace(0, 3, 500))
    >>> plt.semilogx(w, 20 * np.log10(abs(h)))
    >>> plt.title('Elliptical highpass filter fit to constraints')
    >>> plt.xlabel('Frequency [radians / second]')
    >>> plt.ylabel('Amplitude [dB]')
    >>> plt.grid(which='both', axis='both')
    >>> plt.fill([.1, 10,  10,  .1], [1e4, 1e4, -60, -60], '0.9', lw=0) # stop
    >>> plt.fill([30, 30, 1e9, 1e9], [-99,  -3,  -3, -99], '0.9', lw=0) # pass
    >>> plt.axis([1, 300, -80, 3])
    >>> plt.show()

    """
    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2 * (len(wp) - 1)
    filter_type += 1
    if wp[0] >= ws[0]:
        filter_type += 1

    # Pre-warp frequencies for digital filter design
    if not analog:
        passb = tan(pi * wp / 2.0)
        stopb = tan(pi * ws / 2.0)
    else:
        passb = wp * 1.0
        stopb = ws * 1.0

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0] - 1e-12,
                                 args=(0, passb, stopb, gpass, gstop, 'ellip'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1] + 1e-12, passb[1],
                                 args=(1, passb, stopb, gpass, gstop, 'ellip'),
                                 disp=0)
        passb[1] = wp1
        nat = ((stopb * (passb[0] - passb[1])) /
               (stopb ** 2 - passb[0] * passb[1]))
    elif filter_type == 4:  # pass
        nat = ((stopb ** 2 - passb[0] * passb[1]) /
               (stopb * (passb[0] - passb[1])))

    nat = min(abs(nat))

    GSTOP = 10 ** (0.1 * gstop)
    GPASS = 10 ** (0.1 * gpass)
    arg1 = sqrt((GPASS - 1.0) / (GSTOP - 1.0))
    arg0 = 1.0 / nat
    d0 = special.ellipk([arg0 ** 2, 1 - arg0 ** 2])
    d1 = special.ellipk([arg1 ** 2, 1 - arg1 ** 2])
    ord = int(ceil(d0[0] * d1[1] / (d0[1] * d1[0])))

    if not analog:
        wn = arctan(passb) * 2.0 / pi
    else:
        wn = passb

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn


def buttap(N):
    """Return (z,p,k) for analog prototype of Nth order Butterworth filter.

    The filter will have an angular (e.g. rad/s) cutoff frequency of 1.

    """
    if abs(int(N)) != N:
        raise ValueError("Filter order must be a nonnegative integer")
    z = numpy.array([])
    m = numpy.arange(-N+1, N, 2)
    # Middle value is 0 to ensure an exactly real pole
    p = -numpy.exp(1j * pi * m / (2 * N))
    k = 1
    return z, p, k


def cheb1ap(N, rp):
    """
    Return (z,p,k) for Nth order Chebyshev type I analog lowpass filter.

    The returned filter prototype has `rp` decibels of ripple in the passband.

    The filter's angular (e.g. rad/s) cutoff frequency is normalized to 1,
    defined as the point at which the gain first drops below ``-rp``.

    """
    if abs(int(N)) != N:
        raise ValueError("Filter order must be a nonnegative integer")
    elif N == 0:
        # Avoid divide-by-zero error
        # Even order filters have DC gain of -rp dB
        return numpy.array([]), numpy.array([]), 10**(-rp/20)
    z = numpy.array([])

    # Ripple factor (epsilon)
    eps = numpy.sqrt(10 ** (0.1 * rp) - 1.0)
    mu = 1.0 / N * arcsinh(1 / eps)

    # Arrange poles in an ellipse on the left half of the S-plane
    m = numpy.arange(-N+1, N, 2)
    theta = pi * m / (2*N)
    p = -sinh(mu + 1j*theta)

    k = numpy.prod(-p, axis=0).real
    if N % 2 == 0:
        k = k / sqrt((1 + eps * eps))

    return z, p, k


def cheb2ap(N, rs):
    """
    Return (z,p,k) for Nth order Chebyshev type I analog lowpass filter.

    The returned filter prototype has `rs` decibels of ripple in the stopband.

    The filter's angular (e.g. rad/s) cutoff frequency is normalized to 1,
    defined as the point at which the gain first reaches ``-rs``.

    """
    if abs(int(N)) != N:
        raise ValueError("Filter order must be a nonnegative integer")
    elif N == 0:
        # Avoid divide-by-zero warning
        return numpy.array([]), numpy.array([]), 1

    # Ripple factor (epsilon)
    de = 1.0 / sqrt(10 ** (0.1 * rs) - 1)
    mu = arcsinh(1.0 / de) / N

    if N % 2:
        m = numpy.concatenate((numpy.arange(-N+1, 0, 2),
                               numpy.arange(2, N, 2)))
    else:
        m = numpy.arange(-N+1, N, 2)

    z = -conjugate(1j / sin(m * pi / (2.0 * N)))

    # Poles around the unit circle like Butterworth
    p = -exp(1j * pi * numpy.arange(-N+1, N, 2) / (2 * N))
    # Warp into Chebyshev II
    p = sinh(mu) * p.real + 1j * cosh(mu) * p.imag
    p = 1.0 / p

    k = (numpy.prod(-p, axis=0) / numpy.prod(-z, axis=0)).real
    return z, p, k


EPSILON = 2e-16


def _vratio(u, ineps, mp):
    [s, c, d, phi] = special.ellipj(u, mp)
    ret = abs(ineps - s / c)
    return ret


def _kratio(m, k_ratio):
    m = float(m)
    if m < 0:
        m = 0.0
    if m > 1:
        m = 1.0
    if abs(m) > EPSILON and (abs(m) + EPSILON) < 1:
        k = special.ellipk([m, 1 - m])
        r = k[0] / k[1] - k_ratio
    elif abs(m) > EPSILON:
        r = -k_ratio
    else:
        r = 1e20
    return abs(r)


def ellipap(N, rp, rs):
    """Return (z,p,k) of Nth order elliptic analog lowpass filter.

    The filter is a normalized prototype that has `rp` decibels of ripple
    in the passband and a stopband `rs` decibels down.

    The filter's angular (e.g. rad/s) cutoff frequency is normalized to 1,
    defined as the point at which the gain first drops below ``-rp``.

    References
    ----------
    Lutova, Tosic, and Evans, "Filter Design for Signal Processing", Chapters 5
    and 12.

    """
    if abs(int(N)) != N:
        raise ValueError("Filter order must be a nonnegative integer")
    elif N == 0:
        # Avoid divide-by-zero warning
        # Even order filters have DC gain of -rp dB
        return numpy.array([]), numpy.array([]), 10**(-rp/20)
    elif N == 1:
        p = -sqrt(1.0 / (10 ** (0.1 * rp) - 1.0))
        k = -p
        z = []
        return asarray(z), asarray(p), k

    eps = numpy.sqrt(10 ** (0.1 * rp) - 1)
    ck1 = eps / numpy.sqrt(10 ** (0.1 * rs) - 1)
    ck1p = numpy.sqrt(1 - ck1 * ck1)
    if ck1p == 1:
        raise ValueError("Cannot design a filter with given rp and rs"
                         " specifications.")

    val = special.ellipk([ck1 * ck1, ck1p * ck1p])
    if abs(1 - ck1p * ck1p) < EPSILON:
        krat = 0
    else:
        krat = N * val[0] / val[1]

    m = optimize.fmin(_kratio, [0.5], args=(krat,), maxfun=250, maxiter=250,
                      disp=0)
    if m < 0 or m > 1:
        m = optimize.fminbound(_kratio, 0, 1, args=(krat,), maxfun=250,
                               maxiter=250, disp=0)

    capk = special.ellipk(m)

    j = numpy.arange(1 - N % 2, N, 2)
    jj = len(j)

    [s, c, d, phi] = special.ellipj(j * capk / N, m * numpy.ones(jj))
    snew = numpy.compress(abs(s) > EPSILON, s, axis=-1)
    z = 1.0 / (sqrt(m) * snew)
    z = 1j * z
    z = numpy.concatenate((z, conjugate(z)))

    r = optimize.fmin(_vratio, special.ellipk(m), args=(1. / eps, ck1p * ck1p),
                      maxfun=250, maxiter=250, disp=0)
    v0 = capk * r / (N * val[0])

    [sv, cv, dv, phi] = special.ellipj(v0, 1 - m)
    p = -(c * d * sv * cv + 1j * s * dv) / (1 - (d * sv) ** 2.0)

    if N % 2:
        newp = numpy.compress(abs(p.imag) > EPSILON *
                              numpy.sqrt(numpy.sum(p * numpy.conjugate(p),
                                                   axis=0).real),
                              p, axis=-1)
        p = numpy.concatenate((p, conjugate(newp)))
    else:
        p = numpy.concatenate((p, conjugate(p)))

    k = (numpy.prod(-p, axis=0) / numpy.prod(-z, axis=0)).real
    if N % 2 == 0:
        k = k / numpy.sqrt((1 + eps * eps))

    return z, p, k


def besselap(N):
    """Return (z,p,k) for analog prototype of an Nth order Bessel filter.

    The filter is normalized such that the filter asymptotes are the same as
    a Butterworth filter of the same order with an angular (e.g. rad/s)
    cutoff frequency of 1.

    Parameters
    ----------
    N : int
        The order of the Bessel filter to return zeros, poles and gain for.
        Values in the range 0-25 are supported.

    Returns
    -------
    z : ndarray
        Zeros. Is always an empty array.
    p : ndarray
        Poles.
    k : scalar
        Gain. Always 1.

    """
    z = []
    k = 1
    if N == 0:
        p = []
    elif N == 1:
        p = [-1]
    elif N == 2:
        p = [-.8660254037844386467637229 + .4999999999999999999999996j,
             -.8660254037844386467637229 - .4999999999999999999999996j]
    elif N == 3:
        p = [-.9416000265332067855971980,
             -.7456403858480766441810907 - .7113666249728352680992154j,
             -.7456403858480766441810907 + .7113666249728352680992154j]
    elif N == 4:
        p = [-.6572111716718829545787781 - .8301614350048733772399715j,
             -.6572111716718829545787788 + .8301614350048733772399715j,
             -.9047587967882449459642637 - .2709187330038746636700923j,
             -.9047587967882449459642624 + .2709187330038746636700926j]
    elif N == 5:
        p = [-.9264420773877602247196260,
             -.8515536193688395541722677 - .4427174639443327209850002j,
             -.8515536193688395541722677 + .4427174639443327209850002j,
             -.5905759446119191779319432 - .9072067564574549539291747j,
             -.5905759446119191779319432 + .9072067564574549539291747j]
    elif N == 6:
        p = [-.9093906830472271808050953 - .1856964396793046769246397j,
             -.9093906830472271808050953 + .1856964396793046769246397j,
             -.7996541858328288520243325 - .5621717346937317988594118j,
             -.7996541858328288520243325 + .5621717346937317988594118j,
             -.5385526816693109683073792 - .9616876881954277199245657j,
             -.5385526816693109683073792 + .9616876881954277199245657j]
    elif N == 7:
        p = [-.9194871556490290014311619,
             -.8800029341523374639772340 - .3216652762307739398381830j,
             -.8800029341523374639772340 + .3216652762307739398381830j,
             -.7527355434093214462291616 - .6504696305522550699212995j,
             -.7527355434093214462291616 + .6504696305522550699212995j,
             -.4966917256672316755024763 - 1.002508508454420401230220j,
             -.4966917256672316755024763 + 1.002508508454420401230220j]
    elif N == 8:
        p = [-.9096831546652910216327629 - .1412437976671422927888150j,
             -.9096831546652910216327629 + .1412437976671422927888150j,
             -.8473250802359334320103023 - .4259017538272934994996429j,
             -.8473250802359334320103023 + .4259017538272934994996429j,
             -.7111381808485399250796172 - .7186517314108401705762571j,
             -.7111381808485399250796172 + .7186517314108401705762571j,
             -.4621740412532122027072175 - 1.034388681126901058116589j,
             -.4621740412532122027072175 + 1.034388681126901058116589j]
    elif N == 9:
        p = [-.9154957797499037686769223,
             -.8911217017079759323183848 - .2526580934582164192308115j,
             -.8911217017079759323183848 + .2526580934582164192308115j,
             -.8148021112269012975514135 - .5085815689631499483745341j,
             -.8148021112269012975514135 + .5085815689631499483745341j,
             -.6743622686854761980403401 - .7730546212691183706919682j,
             -.6743622686854761980403401 + .7730546212691183706919682j,
             -.4331415561553618854685942 - 1.060073670135929666774323j,
             -.4331415561553618854685942 + 1.060073670135929666774323j]
    elif N == 10:
        p = [-.9091347320900502436826431 - .1139583137335511169927714j,
             -.9091347320900502436826431 + .1139583137335511169927714j,
             -.8688459641284764527921864 - .3430008233766309973110589j,
             -.8688459641284764527921864 + .3430008233766309973110589j,
             -.7837694413101441082655890 - .5759147538499947070009852j,
             -.7837694413101441082655890 + .5759147538499947070009852j,
             -.6417513866988316136190854 - .8175836167191017226233947j,
             -.6417513866988316136190854 + .8175836167191017226233947j,
             -.4083220732868861566219785 - 1.081274842819124562037210j,
             -.4083220732868861566219785 + 1.081274842819124562037210j]
    elif N == 11:
        p = [-.9129067244518981934637318,
             -.8963656705721166099815744 - .2080480375071031919692341j,
             -.8963656705721166099815744 + .2080480375071031919692341j,
             -.8453044014712962954184557 - .4178696917801248292797448j,
             -.8453044014712962954184557 + .4178696917801248292797448j,
             -.7546938934722303128102142 - .6319150050721846494520941j,
             -.7546938934722303128102142 + .6319150050721846494520941j,
             -.6126871554915194054182909 - .8547813893314764631518509j,
             -.6126871554915194054182909 + .8547813893314764631518509j,
             -.3868149510055090879155425 - 1.099117466763120928733632j,
             -.3868149510055090879155425 + 1.099117466763120928733632j]
    elif N == 12:
        p = [-.9084478234140682638817772 - 95506365213450398415258360.0e-27j,
             -.9084478234140682638817772 + 95506365213450398415258360.0e-27j,
             -.8802534342016826507901575 - .2871779503524226723615457j,
             -.8802534342016826507901575 + .2871779503524226723615457j,
             -.8217296939939077285792834 - .4810212115100676440620548j,
             -.8217296939939077285792834 + .4810212115100676440620548j,
             -.7276681615395159454547013 - .6792961178764694160048987j,
             -.7276681615395159454547013 + .6792961178764694160048987j,
             -.5866369321861477207528215 - .8863772751320727026622149j,
             -.5866369321861477207528215 + .8863772751320727026622149j,
             -.3679640085526312839425808 - 1.114373575641546257595657j,
             -.3679640085526312839425808 + 1.114373575641546257595657j]
    elif N == 13:
        p = [-.9110914665984182781070663,
             -.8991314665475196220910718 - .1768342956161043620980863j,
             -.8991314665475196220910718 + .1768342956161043620980863j,
             -.8625094198260548711573628 - .3547413731172988997754038j,
             -.8625094198260548711573628 + .3547413731172988997754038j,
             -.7987460692470972510394686 - .5350752120696801938272504j,
             -.7987460692470972510394686 + .5350752120696801938272504j,
             -.7026234675721275653944062 - .7199611890171304131266374j,
             -.7026234675721275653944062 + .7199611890171304131266374j,
             -.5631559842430199266325818 - .9135900338325109684927731j,
             -.5631559842430199266325818 + .9135900338325109684927731j,
             -.3512792323389821669401925 - 1.127591548317705678613239j,
             -.3512792323389821669401925 + 1.127591548317705678613239j]
    elif N == 14:
        p = [-.9077932138396487614720659 - 82196399419401501888968130.0e-27j,
             -.9077932138396487614720659 + 82196399419401501888968130.0e-27j,
             -.8869506674916445312089167 - .2470079178765333183201435j,
             -.8869506674916445312089167 + .2470079178765333183201435j,
             -.8441199160909851197897667 - .4131653825102692595237260j,
             -.8441199160909851197897667 + .4131653825102692595237260j,
             -.7766591387063623897344648 - .5819170677377608590492434j,
             -.7766591387063623897344648 + .5819170677377608590492434j,
             -.6794256425119233117869491 - .7552857305042033418417492j,
             -.6794256425119233117869491 + .7552857305042033418417492j,
             -.5418766775112297376541293 - .9373043683516919569183099j,
             -.5418766775112297376541293 + .9373043683516919569183099j,
             -.3363868224902037330610040 - 1.139172297839859991370924j,
             -.3363868224902037330610040 + 1.139172297839859991370924j]
    elif N == 15:
        p = [-.9097482363849064167228581,
             -.9006981694176978324932918 - .1537681197278439351298882j,
             -.9006981694176978324932918 + .1537681197278439351298882j,
             -.8731264620834984978337843 - .3082352470564267657715883j,
             -.8731264620834984978337843 + .3082352470564267657715883j,
             -.8256631452587146506294553 - .4642348752734325631275134j,
             -.8256631452587146506294553 + .4642348752734325631275134j,
             -.7556027168970728127850416 - .6229396358758267198938604j,
             -.7556027168970728127850416 + .6229396358758267198938604j,
             -.6579196593110998676999362 - .7862895503722515897065645j,
             -.6579196593110998676999362 + .7862895503722515897065645j,
             -.5224954069658330616875186 - .9581787261092526478889345j,
             -.5224954069658330616875186 + .9581787261092526478889345j,
             -.3229963059766444287113517 - 1.149416154583629539665297j,
             -.3229963059766444287113517 + 1.149416154583629539665297j]
    elif N == 16:
        p = [-.9072099595087001356491337 - 72142113041117326028823950.0e-27j,
             -.9072099595087001356491337 + 72142113041117326028823950.0e-27j,
             -.8911723070323647674780132 - .2167089659900576449410059j,
             -.8911723070323647674780132 + .2167089659900576449410059j,
             -.8584264231521330481755780 - .3621697271802065647661080j,
             -.8584264231521330481755780 + .3621697271802065647661080j,
             -.8074790293236003885306146 - .5092933751171800179676218j,
             -.8074790293236003885306146 + .5092933751171800179676218j,
             -.7356166304713115980927279 - .6591950877860393745845254j,
             -.7356166304713115980927279 + .6591950877860393745845254j,
             -.6379502514039066715773828 - .8137453537108761895522580j,
             -.6379502514039066715773828 + .8137453537108761895522580j,
             -.5047606444424766743309967 - .9767137477799090692947061j,
             -.5047606444424766743309967 + .9767137477799090692947061j,
             -.3108782755645387813283867 - 1.158552841199330479412225j,
             -.3108782755645387813283867 + 1.158552841199330479412225j]
    elif N == 17:
        p = [-.9087141161336397432860029,
             -.9016273850787285964692844 - .1360267995173024591237303j,
             -.9016273850787285964692844 + .1360267995173024591237303j,
             -.8801100704438627158492165 - .2725347156478803885651973j,
             -.8801100704438627158492165 + .2725347156478803885651973j,
             -.8433414495836129204455491 - .4100759282910021624185986j,
             -.8433414495836129204455491 + .4100759282910021624185986j,
             -.7897644147799708220288138 - .5493724405281088674296232j,
             -.7897644147799708220288138 + .5493724405281088674296232j,
             -.7166893842372349049842743 - .6914936286393609433305754j,
             -.7166893842372349049842743 + .6914936286393609433305754j,
             -.6193710717342144521602448 - .8382497252826992979368621j,
             -.6193710717342144521602448 + .8382497252826992979368621j,
             -.4884629337672704194973683 - .9932971956316781632345466j,
             -.4884629337672704194973683 + .9932971956316781632345466j,
             -.2998489459990082015466971 - 1.166761272925668786676672j,
             -.2998489459990082015466971 + 1.166761272925668786676672j]
    elif N == 18:
        p = [-.9067004324162775554189031 - 64279241063930693839360680.0e-27j,
             -.9067004324162775554189031 + 64279241063930693839360680.0e-27j,
             -.8939764278132455733032155 - .1930374640894758606940586j,
             -.8939764278132455733032155 + .1930374640894758606940586j,
             -.8681095503628830078317207 - .3224204925163257604931634j,
             -.8681095503628830078317207 + .3224204925163257604931634j,
             -.8281885016242836608829018 - .4529385697815916950149364j,
             -.8281885016242836608829018 + .4529385697815916950149364j,
             -.7726285030739558780127746 - .5852778162086640620016316j,
             -.7726285030739558780127746 + .5852778162086640620016316j,
             -.6987821445005273020051878 - .7204696509726630531663123j,
             -.6987821445005273020051878 + .7204696509726630531663123j,
             -.6020482668090644386627299 - .8602708961893664447167418j,
             -.6020482668090644386627299 + .8602708961893664447167418j,
             -.4734268069916151511140032 - 1.008234300314801077034158j,
             -.4734268069916151511140032 + 1.008234300314801077034158j,
             -.2897592029880489845789953 - 1.174183010600059128532230j,
             -.2897592029880489845789953 + 1.174183010600059128532230j]
    elif N == 19:
        p = [-.9078934217899404528985092,
             -.9021937639390660668922536 - .1219568381872026517578164j,
             -.9021937639390660668922536 + .1219568381872026517578164j,
             -.8849290585034385274001112 - .2442590757549818229026280j,
             -.8849290585034385274001112 + .2442590757549818229026280j,
             -.8555768765618421591093993 - .3672925896399872304734923j,
             -.8555768765618421591093993 + .3672925896399872304734923j,
             -.8131725551578197705476160 - .4915365035562459055630005j,
             -.8131725551578197705476160 + .4915365035562459055630005j,
             -.7561260971541629355231897 - .6176483917970178919174173j,
             -.7561260971541629355231897 + .6176483917970178919174173j,
             -.6818424412912442033411634 - .7466272357947761283262338j,
             -.6818424412912442033411634 + .7466272357947761283262338j,
             -.5858613321217832644813602 - .8801817131014566284786759j,
             -.5858613321217832644813602 + .8801817131014566284786759j,
             -.4595043449730988600785456 - 1.021768776912671221830298j,
             -.4595043449730988600785456 + 1.021768776912671221830298j,
             -.2804866851439370027628724 - 1.180931628453291873626003j,
             -.2804866851439370027628724 + 1.180931628453291873626003j]
    elif N == 20:
        p = [-.9062570115576771146523497 - 57961780277849516990208850.0e-27j,
             -.9062570115576771146523497 + 57961780277849516990208850.0e-27j,
             -.8959150941925768608568248 - .1740317175918705058595844j,
             -.8959150941925768608568248 + .1740317175918705058595844j,
             -.8749560316673332850673214 - .2905559296567908031706902j,
             -.8749560316673332850673214 + .2905559296567908031706902j,
             -.8427907479956670633544106 - .4078917326291934082132821j,
             -.8427907479956670633544106 + .4078917326291934082132821j,
             -.7984251191290606875799876 - .5264942388817132427317659j,
             -.7984251191290606875799876 + .5264942388817132427317659j,
             -.7402780309646768991232610 - .6469975237605228320268752j,
             -.7402780309646768991232610 + .6469975237605228320268752j,
             -.6658120544829934193890626 - .7703721701100763015154510j,
             -.6658120544829934193890626 + .7703721701100763015154510j,
             -.5707026806915714094398061 - .8982829066468255593407161j,
             -.5707026806915714094398061 + .8982829066468255593407161j,
             -.4465700698205149555701841 - 1.034097702560842962315411j,
             -.4465700698205149555701841 + 1.034097702560842962315411j,
             -.2719299580251652601727704 - 1.187099379810885886139638j,
             -.2719299580251652601727704 + 1.187099379810885886139638j]
    elif N == 21:
        p = [-.9072262653142957028884077,
             -.9025428073192696303995083 - .1105252572789856480992275j,
             -.9025428073192696303995083 + .1105252572789856480992275j,
             -.8883808106664449854431605 - .2213069215084350419975358j,
             -.8883808106664449854431605 + .2213069215084350419975358j,
             -.8643915813643204553970169 - .3326258512522187083009453j,
             -.8643915813643204553970169 + .3326258512522187083009453j,
             -.8299435470674444100273463 - .4448177739407956609694059j,
             -.8299435470674444100273463 + .4448177739407956609694059j,
             -.7840287980408341576100581 - .5583186348022854707564856j,
             -.7840287980408341576100581 + .5583186348022854707564856j,
             -.7250839687106612822281339 - .6737426063024382240549898j,
             -.7250839687106612822281339 + .6737426063024382240549898j,
             -.6506315378609463397807996 - .7920349342629491368548074j,
             -.6506315378609463397807996 + .7920349342629491368548074j,
             -.5564766488918562465935297 - .9148198405846724121600860j,
             -.5564766488918562465935297 + .9148198405846724121600860j,
             -.4345168906815271799687308 - 1.045382255856986531461592j,
             -.4345168906815271799687308 + 1.045382255856986531461592j,
             -.2640041595834031147954813 - 1.192762031948052470183960j,
             -.2640041595834031147954813 + 1.192762031948052470183960j]
    elif N == 22:
        p = [-.9058702269930872551848625 - 52774908289999045189007100.0e-27j,
             -.9058702269930872551848625 + 52774908289999045189007100.0e-27j,
             -.8972983138153530955952835 - .1584351912289865608659759j,
             -.8972983138153530955952835 + .1584351912289865608659759j,
             -.8799661455640176154025352 - .2644363039201535049656450j,
             -.8799661455640176154025352 + .2644363039201535049656450j,
             -.8534754036851687233084587 - .3710389319482319823405321j,
             -.8534754036851687233084587 + .3710389319482319823405321j,
             -.8171682088462720394344996 - .4785619492202780899653575j,
             -.8171682088462720394344996 + .4785619492202780899653575j,
             -.7700332930556816872932937 - .5874255426351153211965601j,
             -.7700332930556816872932937 + .5874255426351153211965601j,
             -.7105305456418785989070935 - .6982266265924524000098548j,
             -.7105305456418785989070935 + .6982266265924524000098548j,
             -.6362427683267827226840153 - .8118875040246347267248508j,
             -.6362427683267827226840153 + .8118875040246347267248508j,
             -.5430983056306302779658129 - .9299947824439872998916657j,
             -.5430983056306302779658129 + .9299947824439872998916657j,
             -.4232528745642628461715044 - 1.055755605227545931204656j,
             -.4232528745642628461715044 + 1.055755605227545931204656j,
             -.2566376987939318038016012 - 1.197982433555213008346532j,
             -.2566376987939318038016012 + 1.197982433555213008346532j]
    elif N == 23:
        p = [-.9066732476324988168207439,
             -.9027564979912504609412993 - .1010534335314045013252480j,
             -.9027564979912504609412993 + .1010534335314045013252480j,
             -.8909283242471251458653994 - .2023024699381223418195228j,
             -.8909283242471251458653994 + .2023024699381223418195228j,
             -.8709469395587416239596874 - .3039581993950041588888925j,
             -.8709469395587416239596874 + .3039581993950041588888925j,
             -.8423805948021127057054288 - .4062657948237602726779246j,
             -.8423805948021127057054288 + .4062657948237602726779246j,
             -.8045561642053176205623187 - .5095305912227258268309528j,
             -.8045561642053176205623187 + .5095305912227258268309528j,
             -.7564660146829880581478138 - .6141594859476032127216463j,
             -.7564660146829880581478138 + .6141594859476032127216463j,
             -.6965966033912705387505040 - .7207341374753046970247055j,
             -.6965966033912705387505040 + .7207341374753046970247055j,
             -.6225903228771341778273152 - .8301558302812980678845563j,
             -.6225903228771341778273152 + .8301558302812980678845563j,
             -.5304922463810191698502226 - .9439760364018300083750242j,
             -.5304922463810191698502226 + .9439760364018300083750242j,
             -.4126986617510148836149955 - 1.065328794475513585531053j,
             -.4126986617510148836149955 + 1.065328794475513585531053j,
             -.2497697202208956030229911 - 1.202813187870697831365338j,
             -.2497697202208956030229911 + 1.202813187870697831365338j]
    elif N == 24:
        p = [-.9055312363372773709269407 - 48440066540478700874836350.0e-27j,
             -.9055312363372773709269407 + 48440066540478700874836350.0e-27j,
             -.8983105104397872954053307 - .1454056133873610120105857j,
             -.8983105104397872954053307 + .1454056133873610120105857j,
             -.8837358034555706623131950 - .2426335234401383076544239j,
             -.8837358034555706623131950 + .2426335234401383076544239j,
             -.8615278304016353651120610 - .3403202112618624773397257j,
             -.8615278304016353651120610 + .3403202112618624773397257j,
             -.8312326466813240652679563 - .4386985933597305434577492j,
             -.8312326466813240652679563 + .4386985933597305434577492j,
             -.7921695462343492518845446 - .5380628490968016700338001j,
             -.7921695462343492518845446 + .5380628490968016700338001j,
             -.7433392285088529449175873 - .6388084216222567930378296j,
             -.7433392285088529449175873 + .6388084216222567930378296j,
             -.6832565803536521302816011 - .7415032695091650806797753j,
             -.6832565803536521302816011 + .7415032695091650806797753j,
             -.6096221567378335562589532 - .8470292433077202380020454j,
             -.6096221567378335562589532 + .8470292433077202380020454j,
             -.5185914574820317343536707 - .9569048385259054576937721j,
             -.5185914574820317343536707 + .9569048385259054576937721j,
             -.4027853855197518014786978 - 1.074195196518674765143729j,
             -.4027853855197518014786978 + 1.074195196518674765143729j,
             -.2433481337524869675825448 - 1.207298683731972524975429j,
             -.2433481337524869675825448 + 1.207298683731972524975429j]
    elif N == 25:
        p = [-.9062073871811708652496104,
             -.9028833390228020537142561 - 93077131185102967450643820.0e-27j,
             -.9028833390228020537142561 + 93077131185102967450643820.0e-27j,
             -.8928551459883548836774529 - .1863068969804300712287138j,
             -.8928551459883548836774529 + .1863068969804300712287138j,
             -.8759497989677857803656239 - .2798521321771408719327250j,
             -.8759497989677857803656239 + .2798521321771408719327250j,
             -.8518616886554019782346493 - .3738977875907595009446142j,
             -.8518616886554019782346493 + .3738977875907595009446142j,
             -.8201226043936880253962552 - .4686668574656966589020580j,
             -.8201226043936880253962552 + .4686668574656966589020580j,
             -.7800496278186497225905443 - .5644441210349710332887354j,
             -.7800496278186497225905443 + .5644441210349710332887354j,
             -.7306549271849967721596735 - .6616149647357748681460822j,
             -.7306549271849967721596735 + .6616149647357748681460822j,
             -.6704827128029559528610523 - .7607348858167839877987008j,
             -.6704827128029559528610523 + .7607348858167839877987008j,
             -.5972898661335557242320528 - .8626676330388028512598538j,
             -.5972898661335557242320528 + .8626676330388028512598538j,
             -.5073362861078468845461362 - .9689006305344868494672405j,
             -.5073362861078468845461362 + .9689006305344868494672405j,
             -.3934529878191079606023847 - 1.082433927173831581956863j,
             -.3934529878191079606023847 + 1.082433927173831581956863j,
             -.2373280669322028974199184 - 1.211476658382565356579418j,
             -.2373280669322028974199184 + 1.211476658382565356579418j]
    else:
        raise ValueError("Bessel Filter not supported for order %s" % N)

    return asarray(z), asarray(p), k

filter_dict = {'butter': [buttap, buttord],
               'butterworth': [buttap, buttord],

               'cauer': [ellipap, ellipord],
               'elliptic': [ellipap, ellipord],
               'ellip': [ellipap, ellipord],

               'bessel': [besselap],

               'cheby1': [cheb1ap, cheb1ord],
               'chebyshev1': [cheb1ap, cheb1ord],
               'chebyshevi': [cheb1ap, cheb1ord],

               'cheby2': [cheb2ap, cheb2ord],
               'chebyshev2': [cheb2ap, cheb2ord],
               'chebyshevii': [cheb2ap, cheb2ord],
               }

band_dict = {'band': 'bandpass',
             'bandpass': 'bandpass',
             'pass': 'bandpass',
             'bp': 'bandpass',

             'bs': 'bandstop',
             'bandstop': 'bandstop',
             'bands': 'bandstop',
             'stop': 'bandstop',

             'l': 'lowpass',
             'low': 'lowpass',
             'lowpass': 'lowpass',
             'lp': 'lowpass',

             'high': 'highpass',
             'highpass': 'highpass',
             'h': 'highpass',
             'hp': 'highpass',
             }
