# This program is public domain
# Authors: Paul Kienzle, Nadav Horesh
"""
Chirp z-transform.

We provide two interfaces to the chirp z transform: an object interface
which precalculates part of the transform and can be applied efficiently
to many different data sets, and a functional interface which is applied
only to the given data set.

Transforms
----------

CZT : callable (x, axis=-1) -> array
   define a chirp z-transform that can be applied to different signals
ZoomFFT : callable (x, axis=-1) -> array
   define a Fourier transform on a range of frequencies

Functions
---------

czt : array
   compute the chirp z-transform for a signal
zoomfft : array
   compute the Fourier transform on a range of frequencies
"""

from __future__ import division, absolute_import, print_function

import cmath

import numpy as np
from numpy import pi, arange
from scipy.fftpack import fft, ifft
from scipy.fftpack.helper import _next_opt_len

__all__ = ['czt', 'zoomfft', 'CZT', 'ZoomFFT', 'czt_points']


def _validate_sizes(n, m):
    if n < 1:
        raise ValueError("Invalid number of CZT data "
                         "points (%d) specified." % n)

    if int(n) != n:
        raise ValueError('n must be a positive integer')

    if m is None:
        m = n
    elif m < 1:
        raise ValueError("Invalid number of CZT output "
                         "points (%d) specified." % n)
    elif int(m) != m:
        raise ValueError('m must be a positive integer')

    return m


def czt_points(m, w=None, a=1+0j):
    """
    The points at which the Z-transform is computed when doing a `CZT`
    with the same arguments.

    Parameters
    ----------
    m : int, optional
        The number of points desired.
    w : complex, optional
        The ratio between points in each step.
    a : complex, optional
        The starting point in the complex plane.  The default is 1.

    Returns
    -------
    out : ndarray
        The points in the Z plane at which the CZT samples the Z-transform,
        as complex numbers.

    Examples
    --------
    Plot the points of a 16-point FFT:

    >>> import matplotlib.pyplot as plt
    >>> points = czt_points(16)
    >>> plt.plot(np.real(points), np.imag(points), 'o')
    >>> plt.margins(0.1, 0.1); plt.axis('equal')
    >>> plt.show()

    and a 91-point logarithmic spiral:

    >>> m, w, a = 91, 0.995*np.exp(-1j*np.pi*.05), 0.8*np.exp(1j*np.pi/6)
    >>> points = czt_points(m, w, a)
    >>> plt.plot(np.real(points), np.imag(points), 'o')
    >>> plt.margins(0.1, 0.1); plt.axis('equal')
    >>> plt.show()

    """
    m = _validate_sizes(1, m)

    k = arange(m)

    a = 1.0 * a  # at least float

    if w is None:
        # Nothing specified, default to FFT
        return a * np.exp(2j * pi * k / m)
    else:
        # w specified
        w = 1.0 * w  # at least float
        return a * w**-k


class CZT(object):
    """
    Create a chirp Z-transform function.

    Transform to compute the frequency response around a spiral.
    Objects of this class are callables which can compute the
    chirp Z-transform on their inputs.  This object precalculates the constant
    chirps used in the given transform.

    Parameters
    ----------
    n : int
        The size of the signal.
    m : int, optional
        The number of output points desired.  The default is `n`.
    w : complex, optional
        The ratio between points in each step.  This must be precise or the
        accumulated error will degrade the tail of the output sequence.
    a : complex, optional
        The starting point in the complex plane.  The default is 1+0j.

    Returns
    -------
    f : CZT
        callable object ``f(x, axis=-1)`` for computing the chirp z-transform
        on `x`

    See Also
    --------
    ZoomFFT : for a friendly interface to partial FFT calculations

    Notes
    -----
    If `w` does not lie on the unit circle, then the transform will be
    around a spiral with exponentially-increasing radius.  Regardless,
    angle will increase linearly.

    For transforms that do lie on the unit circle, accuracy is better when
    using `ZoomFFT`, since any numerical error in `w` is
    accumulated for long data lengths, drifting away from the unit circle.

    The chirp z-transform can be faster than an equivalent FFT with
    zero padding.  Try it with your own array sizes to see.

    However, the chirp z-transform is considerably less precise than the
    equivalent zero-padded FFT, with differences on the order of 1e-7 from the
    direct transform rather than 1e-15 as seen with zero-padding.

    Its most common application is computing large prime-length Fourier
    transforms in O(N log N) time, rather than the O(N**2) time required by
    the direct DFT calculation used in `fft`.

    References
    ----------
    .. [1] Leo I. Bluestein, "A linear filtering approach to the computation
           of the discrete Fourier transform," Northeast Electronics Research
           and Engineering Meeting Record 10, 218-219 (1968).
    .. [2] Rabiner, Schafer, and Rader, "The chirp z-transform algorithm and
           its application," Bell Syst. Tech. J. 48, 1249-1292 (1969).

    Examples
    --------
    Compute multiple prime-length FFTs:

    >>> a = np.random.rand(7)
    >>> b = np.random.rand(7)
    >>> c = np.random.rand(7)
    >>> czt_7 = CZT(n=7)
    >>> A = czt_7(a)
    >>> B = czt_7(b)
    >>> C = czt_7(c)

    Display the points at which the FFT is calculated:

    >>> czt_7.points()
    array([ 1.00000000+0.j        ,  0.62348980+0.78183148j,
           -0.22252093+0.97492791j, -0.90096887+0.43388374j,
           -0.90096887-0.43388374j, -0.22252093-0.97492791j,
            0.62348980-0.78183148j])
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(np.real(czt_7.points()), np.imag(czt_7.points()), 'o')
    >>> plt.show()

    """
    def __init__(self, n, m=None, w=None, a=1+0j):
        m = _validate_sizes(n, m)

        k = arange(max(m, n), dtype=np.min_scalar_type(-max(m, n)**2))

        if w is None:
            # Nothing specified, default to FFT-like
            w = cmath.exp(-2j*pi/m)
            wk2 = np.exp(-(1j * pi * k**2) / m)
        else:
            # w specified
            wk2 = w**(k**2/2.)

        a = 1.0 * a  # at least float

        self.w, self.a = w, a
        self.m, self.n = m, n

        nfft = _next_opt_len(n + m - 1)
        self._Awk2 = (a**-k * wk2)[:n]
        self._nfft = nfft
        self._Fwk2 = fft(1/np.hstack((wk2[n-1:0:-1], wk2[:m])), nfft)
        self._wk2 = wk2[:m]
        self._yidx = slice(n-1, n+m-1)

    def __call__(self, x, axis=-1):
        """
        Parameters
        ----------
        x : array
            The signal to transform.
        axis : int, optional
            Array dimension to operate over.  The default is the final
            dimension.

        Returns
        -------
        out : ndarray
            An array of the same dimensions as `x`, but with the length of the
            transformed axis set to `m`.  Note that this is a view on a much
            larger array.
        """
        x = np.asarray(x)
        if x.shape[axis] != self.n:
            raise ValueError("CZT defined for length %d, not %d" %
                             (self.n, x.shape[axis]))
        # Calculate transpose coordinates, to allow operation on any given axis
        trnsp = np.arange(x.ndim)
        trnsp[[axis, -1]] = [-1, axis]
        x = x.transpose(*trnsp)
        y = ifft(self._Fwk2 * fft(x*self._Awk2, self._nfft))
        y = y[..., self._yidx] * self._wk2
        return y.transpose(*trnsp)

    def points(self):
        """
        The points at which the Z-transform is computed when calling this
        `CZT`.
        """
        return czt_points(self.m, self.w, self.a)


class ZoomFFT(CZT):
    """
    Create a zoom FFT transform function.

    This is a specialization of the chirp Z-transform (`CZT`) for a set of
    equally-spaced frequencies around the unit circle.

    Parameters
    ----------
    n : int
        The size of the signal.
    fn : array_like, optional
        A length-2 sequence [f1, f2] giving the frequency range, or a scalar,
        for which the range [0, fn] is assumed.
    m : int, optional
        The number of output points desired.  The default is `n`.
    Fs : float, optional
        The sampling frequency. (default=2)

    Returns
    -------
    f : ZoomFFT
        A callable object ``f(x, axis=-1)`` for computing the zoom FFT on `x`.

    Notes
    -----
    Sampling frequency is 1/dt, the time step between samples in the
    signal `x`.  The unit circle corresponds to frequencies from 0 up
    to the sampling frequency.  The default sampling frequency of 2
    means that f1, f2 values up to the Nyquist frequency are in the
    range [0, 1). For f1, f2 values expressed in radians, a sampling
    frequency of 1/pi should be used.

    Remember that a zoom FFT can only interpolate the points of the existing
    FFT.  It cannot help to resolve two separate nearby frequencies.
    Frequency resolution can only be increased by increasing acquisition
    time.

    Examples
    --------
    To plot the transform results use something like the following:

    >>> t = np.linspace(0, 1, 1021)
    >>> x = np.cos(2*np.pi*15*t) + np.sin(2*np.pi*17*t)
    >>> f1, f2 = 5, 27
    >>> transform = ZoomFFT(len(x), [f1, f2], len(x), 1021)
    >>> X = transform(x)
    >>> f = np.linspace(f1, f2, len(x))
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(f, 20*np.log10(np.abs(X)))
    >>> plt.show()

    """
    def __init__(self, n, fn, m=None, Fs=2):
        m = _validate_sizes(n, m)

        k = arange(max(m, n), dtype=np.min_scalar_type(-max(m, n)**2))

        if np.size(fn) == 2:
            f1, f2 = fn
        elif np.size(fn) == 1:
            f1, f2 = 0.0, fn
        else:
            raise ValueError('fn must be a scalar or 2-length sequence')

        self.f1, self.f2, self.Fs = f1, f2, Fs

        scale = ((f2 - f1) * m) / (Fs * (m - 1))
        a = cmath.exp(2j * pi * f1/Fs)
        w = cmath.exp(-2j*pi/m * scale)
        wk2 = np.exp(-(1j * pi * scale * k**2) / m)

        self.w, self.a = w, a
        self.m, self.n = m, n

        nfft = _next_opt_len(n + m - 1)
        self._Awk2 = (a**-k * wk2)[:n]
        self._nfft = nfft
        self._Fwk2 = fft(1/np.hstack((wk2[n-1:0:-1], wk2[:m])), nfft)
        self._wk2 = wk2[:m]
        self._yidx = slice(n-1, n+m-1)


def czt(x, m=None, w=None, a=1+0j, axis=-1):
    """
    Compute the frequency response around a spiral in the Z plane.

    Parameters
    ----------
    x : array
        The set of data to transform.
    m : int, optional
        The number of points desired. Default is the length of the input data.
    w : complex, optional
        The ratio between points in each step.
    a : complex, optional
        The starting point in the complex plane.  Default is 1+0j.
    axis : int, optional
        Array dimension to operate over.  Default is the final dimension.

    Returns
    -------
    out : ndarray
        An array of the same dimensions as `x`, but with the length of the
        transformed axis set to `m`.  Note that this is a view on a much
        larger array.  To save space, you may want to call it as
        ``y = ascontiguousarray(czt(x))``

    See Also
    --------
    zoomfft : for a friendlier interface to partial FFT calculations

    Notes
    -----
    If the transform needs to be repeated, use `CZT` to construct a
    specialized transform function which can be reused without
    recomputing constants.

    Examples
    --------
    `czt` can quickly compute a prime-length FFT:

    >>> n = 30011  # prime
    >>> a = np.ones(n)  # DC (spectrum is [n, 0, 0, 0, ...])
    >>> f = fft(a)  # takes 1.4 seconds
    >>> c = czt(a)  # takes 21 milliseconds
    >>> np.allclose(f, c, atol=1e-6)
    True

    However, the CZT has more error:

    >>> print(abs(np.amax(f[1:])))
    1.59615876072e-10
    >>> print(abs(np.amax(c[1:])))
    2.25172229068e-09
    >>> print(np.linalg.norm(f[1:]))
    9.82334904574e-09
    >>> print(np.linalg.norm(c[1:]))
    1.08222546282e-07

    With `scale` parameter, CZT can be used to calculate even-length rfft:

    >>> from numpy.fft import rfft
    >>> a = np.random.rand(1024)
    >>> r = rfft(a)
    >>> c = czt(a, m=513, scale=0.5*513/512)
    >>> np.allclose(c, r)
    True

    """
    x = np.asarray(x)
    transform = CZT(x.shape[axis], m=m, w=w, a=a)
    return transform(x, axis=axis)


def zoomfft(x, fn, m=None, Fs=2, axis=-1):
    """
    Compute the DFT of `x` only for frequencies in range `fn`.

    Parameters
    ----------
    x : array
        The input signal.
    fn : array_like, optional
        A length-2 sequence [f1, f2] giving the frequency range, or a scalar,
        for which the range [0, fn] is assumed.
    m : int, optional
        The number of points to evaluate.  The default is the length of `x`.
    Fs : float, optional
        The sampling frequency.  With a sampling frequency of
        10kHz for example, the range f1 and f2 can be expressed in kHz.
        The default sampling frequency is 2, so f1 and f2 should be
        in the range 0, 1 to keep the transform below the Nyquist
        frequency.
    axis : int, optional
        The array dimension the transform operates over.  The default is the
        final dimension.

    Returns:
    -------
    out : ndarray
        The transformed signal.  The Fourier transform will be calculated
        at the points f1, f1+df, f1+2df, ..., f2, where df=(f2-f1)/m.

    Notes
    -----
    ``zoomfft(x, [0, 2-2./len(x)])`` is equivalent to ``fft(x)``.

    To graph the magnitude of the resulting transform, use::

        plot(linspace(f1, f2, m), abs(zoomfft(x, [f1, f2], m)))

    If the transform needs to be repeated, use `ZoomFFT` to construct
    a specialized transform function which can be reused without
    recomputing constants.
    """
    x = np.asarray(x)
    transform = ZoomFFT(x.shape[axis], fn, m=m, Fs=Fs)
    return transform(x, axis=axis)
