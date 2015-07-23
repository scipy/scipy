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
ScaledFFT : callable (x, axis=-1) -> array
   define a limited frequency FFT

Functions
---------

czt : array
   compute the chirp z-transform for a signal
zoomfft : array
   compute the Fourier transform on a range of frequencies
scaledfft : array
   compute a limited frequency FFT for a signal
"""
from __future__ import division, absolute_import, print_function

import cmath

import numpy as np
from numpy import pi, arange
from scipy.fftpack import fft, ifft, fftshift
from scipy.fftpack.helper import _next_regular

__all__ = ['czt', 'zoomfft', 'scaledfft', 'CZT', 'ZoomFFT', 'ScaledFFT',
           'czt_points']


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


def czt_points(m, w=None, a=1, scale=None):
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
    scale : float, optional
        Frequency scaling factor. This is a simpler and more precise way
        to specify FFT-like transforms along a circle.  For
        instance, when assigning scale=0.5, the resulting FT will span
        half of the frequency range that an FFT would produce, at half of
        the frequency step size.
        This is an alternative to `w`, and both cannot be specified at the
        same time.

    Returns
    -------
    out : ndarray
        The points in the Z plane at which the CZT samples the Z-transform,
        as complex numbers.

    Examples
    --------
    Plot the points of a 16-point FFT:

    >>> import matplotlib.pyplot as plt
    >>> points = czt_points(16, a=1+0j, scale=1)
    >>> plt.plot(real(points), imag(points), 'o')
    >>> plt.margins(0.1, 0.1); plt.axis('equal')
    >>> plt.show()

    and a 91-point logarithmic spiral:

    >>> points = czt_points(91, 0.995*exp(-1j*pi*.05), 0.8*exp(1j*pi/6))
    >>> plt.plot(real(points), imag(points), 'o')
    >>> plt.margins(0.1, 0.1); plt.axis('equal')
    >>> plt.show()

    """
    m = _validate_sizes(1, m)

    k = arange(0, m)

    a = 1.0 * a  # at least float

    if w is not None and scale is not None:
        raise ValueError('Only w or scale can be specified; not both.')
    elif w is None:
        if scale is None:
            # Nothing specified, default to FFT-like
            scale = 1

        return a * np.exp(2j * pi * scale * arange(m) / m)
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
    scale : float, optional
        Frequency scaling factor. This is a simpler and more accurate way
        to specify FFT-like transforms along a circle.  For
        instance, when assigning scale=0.5, the resulting FT will span
        half of the frequency range that an FFT would produce, at half of
        the frequency step size.  This is an alternative to `w`, and both
        cannot be specified at the same time.

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
    using `scale` rather than `w`, as any numerical error in `w` is
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
    array([ 1.00000000+0.j,  0.62348980+0.78183148j,
      -0.22252093+0.97492791j, -0.90096887+0.43388374j,
      -0.90096887-0.43388374j, -0.22252093-0.97492791j,
       0.62348980-0.78183148j])
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(real(czt_7.points()), imag(czt_7.points()), 'o')
    >>> plt.show()

    """
    def __init__(self, n, m=None, w=None, a=1, scale=None):
        m = _validate_sizes(n, m)

        k = arange(max(m, n), dtype=np.min_scalar_type(-max(m, n)**2))

        if w is not None and scale is not None:
            raise ValueError('Only w or scale can be specified; not both.')
        elif w is None:
            if scale is None:
                # Nothing specified, default to FFT-like
                scale = 1

            w = cmath.exp(-2j*pi/m * scale)
            wk2 = np.exp(-(1j * pi * scale * k**2) / m)
        else:
            # w specified
            wk2 = w**(k**2/2.)

        a = 1.0 * a  # at least float

        self.w, self.a = w, a
        self.m, self.n = m, n

        nfft = _next_regular(n+m-1)
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

    >>> t = linspace(0, 1, 1021)
    >>> x = cos(2*pi*15*t) + sin(2*pi*17*t)
    >>> f1, f2 = 5, 27
    >>> transform = ZoomFFT(len(x), [f1, f2], len(x), 1021)
    >>> X = transform(x)
    >>> f = linspace(f1, f2, len(x))
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(f, 20*log10(abs(X)))
    >>> plt.show()

    """
    def __init__(self, n, fn, m=None, Fs=2):
        m = _validate_sizes(n, m)

        if np.size(fn) == 2:
            f1, f2 = fn
        elif np.size(fn) == 1:
            f1, f2 = 0.0, fn
        else:
            raise ValueError('fn must be a scalar or 2-length sequence')

        scale = ((f2 - f1) * m) / (Fs * (m - 1))

        a = cmath.exp(2j * pi * f1/Fs)
        CZT.__init__(self, n, m=m, a=a, scale=scale)
        self.f1, self.f2, self.Fs = f1, f2, Fs


class ScaledFFT(CZT):
    """
    Create a scaled FFT transform function.

    Similar to FFT, where the frequency range is scaled and divided
    into m-1 equal steps.  Like the FFT, frequencies are arranged from
    0 to scale*Fs/2-delta followed by -scale*Fs/2 to -delta, where delta
    is the step size scale*Fs/m for sampling frequency Fs. The intended
    use is in a convolution of two signals, each has its own sampling step.

    This is equivalent to:

        fftshift(zoomfft(x, [-scale, scale*(m-2.)/m], m=m))

    Parameters
    ----------
    n : int
        Size of the signal
    m : int, optional
        The size of the output. Default: m=n
    scale : float, optional
        Frequency scaling factor. Default: scale=1.0

    Returns
    -------
    f : ScaledFFT
        callable function ``f(x, axis=-1)`` for computing the scaled FFT on `x`.

    Examples
    --------
    Show the relationship between scaled FFT and zoom FFT:

    >>> scale, m, n = 0.25, 200, 100
    >>> v = np.zeros(n)
    >>> v[[1, 5, 21]] = 1
    >>> shift = scipy.fftpack.fftshift
    >>> x = np.linspace(-0.5, 0.5 - 1./n, n)
    >>> xz = np.linspace(-scale * 0.5, scale * 0.5 * (m - 2.) / m, m)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, shift(abs(fft(v))), label='fft')
    >>> plt.plot(x, shift(abs(scaledfft(v))), 'g.', label='x1 scaled fft')
    >>> plt.plot(xz, abs(zoomfft(v, [-scale, scale * (m - 2.) / m], m=m)),
    ...          'b+', label='zoomfft')
    >>> plt.plot(xz, shift(abs(scaledfft(v, m=m, scale=scale))),
    ...          'rx', label='x' + str(scale) + ' scaled fft')
    >>> plt.yscale('log')
    >>> plt.legend(loc='best')
    >>> plt.show()

    """
    def __init__(self, n, m=None, scale=1.0):
        m = _validate_sizes(n, m)

        w = np.exp(-2j * pi / m * scale)
        a = w**((m+1)//2)

        CZT.__init__(self, n=n, m=m, a=a, scale=scale)
        self.scale = scale

    def __call__(self, x, axis=-1):
        return fftshift(CZT.__call__(self, x, axis), axes=(axis,))
    __call__.__doc__ = CZT.__call__.__doc__


def scaledfft(x, m=None, scale=1.0, axis=-1):
    """
    Limited frequency FFT.

    See `ScaledFFT` doc for details

    Parameters
    ----------
    x : array
        input array
    m : int, optional
        The length of the output signal
    scale : float, optional
        A frequency scaling factor
    axis : int, optional
        The array dimension to operate over.  The default is the
        final dimension.

    Returns
    -------
    out : ndarray
        An array of the same rank of `x`, but with the size of
        the `axis` dimension set to `m`
    """
    x = np.asarray(x)
    transform = ScaledFFT(x.shape[axis], m, scale)
    return transform(x, axis)


def czt(x, m=None, w=None, a=1, scale=None, axis=-1):
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
    scale : float, optional
        The frequency step scale (relative to the normal DFT frequency step).
        Cannot be specified at the same time as `w`.
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
    CZT can be used to calculate even-length rfft:

    >>> from numpy.fft import rfft
    >>> a = np.random.rand(1024)
    >>> r = rfft(a)
    >>> c = czt(a, m=513, scale=0.5*513/512)
    >>> allclose(c, r)
    True

    """
    x = np.asarray(x)
    transform = CZT(x.shape[axis], m=m, w=w, a=a, scale=scale)
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


def _test1(x, show=False, plots=[1, 2, 3, 4]):
    # Normal fft and zero-padded fft equivalent to 10x oversampling
    over = 10
    w = np.linspace(0, 2-2./len(x), len(x))
    y = fft(x)
    wover = np.linspace(0, 2-2./(over*len(x)), over*len(x))
    yover = fft(x, over*len(x))

    # Check that zoomfft is the equivalent of fft
    y1 = zoomfft(x, [0, 2-2./len(y)])

    # Check that zoomfft with oversampling is equivalent to zero padding
    y2 = zoomfft(x, [0, 2-2./len(yover)], m=len(yover))

    # Check that zoomfft works on a subrange
    f1, f2 = w[3], w[6]
    y3 = zoomfft(x, [f1, f2], m=3*over + 1)
    w3 = np.linspace(f1, f2, len(y3))
    idx3 = slice(3*over, 6*over + 1)

    if not show:
        plots = []
    if plots != []:
        import pylab
    if 0 in plots:
        pylab.figure(0)
        pylab.plot(x)
        pylab.ylabel('Intensity')
    if 1 in plots:
        pylab.figure(1)
        pylab.subplot(311)
        pylab.plot(w, abs(y), 'o', w, abs(y1))
        pylab.legend(['fft', 'zoom'])
        pylab.ylabel('Magnitude')
        pylab.title('FFT equivalent')
        pylab.subplot(312)
        pylab.plot(w, np.angle(y), 'o', w, np.angle(y1))
        pylab.legend(['fft', 'zoom'])
        pylab.ylabel('Phase (radians)')
        pylab.subplot(313)
        pylab.plot(w, abs(y)-abs(y1))  # ,w, np.angle(y)-np.angle(y1))
        # pylab.legend(['magnitude', 'phase'])
        pylab.ylabel('Residuals')
    if 2 in plots:
        pylab.figure(2)
        pylab.subplot(211)
        pylab.plot(w, abs(y), 'o', wover, abs(y2), wover, abs(yover))
        pylab.ylabel('Magnitude')
        pylab.title('Oversampled FFT')
        pylab.legend(['fft', 'zoom', 'pad'])
        pylab.subplot(212)
        pylab.plot(wover, abs(yover)-abs(y2),
                   w, abs(y)-abs(y2[0::over]), 'o',
                   w, abs(y)-abs(yover[0::over]), 'x')
        pylab.legend(['pad-zoom', 'fft-zoom', 'fft-pad'])
        pylab.ylabel('Residuals')
    if 3 in plots:
        pylab.figure(3)
        ax1 = pylab.subplot(211)
        pylab.plot(w, abs(y), 'o', w3, abs(y3), wover, abs(yover),
                   w[3:7], abs(y3[::over]), 'x')
        pylab.title('Zoomed FFT')
        pylab.ylabel('Magnitude')
        pylab.legend(['fft', 'zoom', 'pad'])
        pylab.plot(w3, abs(y3), 'x')
        ax1.set_xlim(f1, f2)
        ax2 = pylab.subplot(212)
        pylab.plot(wover[idx3], abs(yover[idx3])-abs(y3),
                   w[3:7], abs(y[3:7])-abs(y3[::over]), 'o',
                   w[3:7], abs(y[3:7])-abs(yover[3*over:6*over+1:over]), 'x')
        pylab.legend(['pad-zoom', 'fft-zoom', 'fft-pad'])
        ax2.set_xlim(f1, f2)
        pylab.ylabel('Residuals')
    if plots != []:
        pylab.show()


def test(demo=None, plots=[1, 2, 3]):
    # 0: Gauss
    t = np.linspace(-2, 2, 128)
    x = np.exp(-t**2/0.01)
    _test1(x, show=(demo == 0), plots=plots)

    # 1: Linear
    x = [1, 2, 3, 4, 5, 6, 7]
    _test1(x, show=(demo == 1), plots=plots)

    # Check near powers of two
    _test1(range(126-31), show=False)
    _test1(range(127-31), show=False)
    _test1(range(128-31), show=False)
    _test1(range(129-31), show=False)
    _test1(range(130-31), show=False)

    # Check transform on n-D array input
    x = np.reshape(np.arange(3*2*28), (3, 2, 28))

    # 2: Random (not a test condition)
    if demo == 2:
        x = np.random.rand(101)
        _test1(x, show=True, plots=plots)

    # 3: Spikes
    t = np.linspace(0, 1, 128)
    x = np.sin(2*pi*t*5)+np.sin(2*pi*t*13)
    _test1(x, show=(demo == 3), plots=plots)

    # 4: Sines
    x = np.zeros(100, dtype=complex)
    x[[1, 5, 21]] = 1
    _test1(x, show=(demo == 4), plots=plots)

    # 5: Sines plus complex component
    x += 1j*np.linspace(0, 0.5, x.shape[0])
    _test1(x, show=(demo == 5), plots=plots)

    # 6: Scaled FFT on complex sines
    x += 1j*np.linspace(0, 0.5, x.shape[0])
    if demo == 6:
        demo_scaledfft(x, 0.25, 200)
    _testscaled(x)


def demo_scaledfft(v, scale, m):
    import pylab
    shift = pylab.fftshift
    n = len(v)
    x = pylab.linspace(-0.5, 0.5 - 1./n, n)
    xz = pylab.linspace(-scale*0.5, scale*0.5*(m-2.)/m, m)
    pylab.figure()
    pylab.plot(x, shift(abs(fft(v))), label='fft')
    pylab.plot(x, shift(abs(scaledfft(v))), 'ro', label='x1 scaled fft')
    pylab.plot(xz, abs(zoomfft(v, [-scale, scale*(m-2.)/m], m=m)),
               'bo', label='zoomfft')
    pylab.plot(xz, shift(abs(scaledfft(v, m=m, scale=scale))),
               'gx', label='x'+str(scale)+' scaled fft')
    pylab.gca().set_yscale('log')
    pylab.legend()
    pylab.show()

if __name__ == "__main__":
    # Choose demo in [0,4] to show plot, or None for testing only
    test(demo=None)
