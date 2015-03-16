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
   define a chirp-z transform that can be applied to different signals
ZoomFFT : callable (x, axis=-1) -> array
   define a Fourier transform on a range of frequencies
ScaledFFT : callable (x, axis=-1) -> array
   define a limited frequency FFT

Functions
---------

czt : array
   compute the chirp-z transform for a signal
zoomfft : array
   compute the Fourier transform on a range of frequencies
scaledfft : array
   compute a limited frequency FFT for a signal
"""
from __future__ import division, absolute_import, print_function

__all__ = ['czt', 'zoomfft', 'scaledfft', 'CZT', 'ZoomFFT', 'ScaledFFT']

import cmath

import numpy as np
from numpy import pi, arange
from scipy.fftpack import fft, ifft, fftshift
from scipy.fftpack.helper import _next_regular


class CZT:
    """
    Chirp-Z Transform.

    Transform to compute the frequency response around a spiral.
    Objects of this class are callables which can compute the
    chirp-z transform on their inputs.  This object precalculates
    constants for the given transform.

    If w does not lie on the unit circle, then the transform will be
    around a spiral with exponentially increasing radius.  Regardless,
    angle will increase linearly.

    The chirp-z transform can be faster than an equivalent FFT with
    zero padding.  Try it with your own array sizes to see.  It is
    theoretically faster for large prime Fourier transforms, but not
    in practice.

    The chirp-z transform is considerably less precise than the
    equivalent zero-padded FFT, with differences on the order of 1e-7
    from the direct transform rather than the on the order of 1e-15 as
    seen with zero-padding.

    See zoomfft for a friendlier interface to partial FFT calculations.
    """
    def __init__(self, n, m=None, w=None, a=1, factor=None):
        """
        Chirp-Z transform definition.

        Parameters:
        ----------
        n : int
          The size of the signal
        m : int, optional
          The number of points desired.  The default is the length of
          the input data.
        w : complex, optional
          The ratio between points in each step.
        a : complex, optional
          The starting point in the complex plane.  The default is 1.
        factor : float, optional
          Frequency scaling factor. For instance, when assigning factor=0.5,
          the resulting FT will span half of the frequency range that an FFT
          would produce, at half of the frequency step size.  This is an
          alternative to `w`, and both cannot be specified at the same time.

        Returns:
        --------
        CZT
          callable object ``f(x, axis=-1)`` for computing the chirp-z transform
          on `x`
        """
        if n < 1:
            raise ValueError("Invalid number of CZT data "
                             "points (%d) specified." % n)

        if m is None:
            m = n

        if w is not None and factor is not None:
            raise ValueError('Only w or factor can be specified; not both.')
        elif w is None and factor is None:
            # Default to FFT
            w = cmath.exp(-2j*pi/m)
        elif w is None:
            w = cmath.exp(-2j*pi/m * factor)

        self.w, self.a = w, a
        self.m, self.n = m, n

        k = arange(max(m, n), dtype=np.min_scalar_type(-max(m, n)**2))
        wk2 = w**(k**2/2.)
        nfft = _next_regular(n+m-1)
        self._Awk2 = (a**-k * wk2)[:n]
        self._nfft = nfft
        self._Fwk2 = fft(1/np.hstack((wk2[n-1:0:-1], wk2[:m])), nfft)
        self._wk2 = wk2[:m]
        self._yidx = slice(n-1, n+m-1)

    def __call__(self, x, axis=-1):
        """
        Parameters:
        ----------
        x : array
          The signal to transform.
        axis : int, optional
          Array dimension to operate over.  The default is the final
          dimension.

        Returns:
        -------
          An array of the same dimensions as x, but with the length of the
          transformed axis set to m.  Note that this is a view on a much
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


class ZoomFFT(CZT):
    """
    Zoom FFT transform.

    This is a specialization of the chirp Z transform, CZT for a set of
    equally spaced frequencies.
    """
    def __init__(self, n, fn, m=None, Fs=2):
        """
        Zoom FFT transform.

        Defines a Fourier transform for a set of equally spaced frequencies.

        Parameters:
        ----------
        n : int
          size of the signal
        fn : array_like, optional
          A scalar or length-2 sequence giving the start and end frequencies.
          If a scalar, use 0 to f1.
        m : int, optional
          size of the output
        Fs : float, optional
          sampling frequency (default=2)

        Returns:
        -------
        A ZoomFFT instance
          A callable object f(x, axis=-1) for computing the zoom FFT on x.

        Sampling frequency is 1/dt, the time step between samples in the
        signal x.  The unit circle corresponds to frequencies from 0 up
        to the sampling frequency.  The default sampling frequency of 2
        means that f1, f2 values up to the Nyquist frequency are in the
        range [0, 1). For f1, f2 values expressed in radians, a sampling
        frequency of 1/pi should be used.

        To plot the transform results use something like the following:

            t = transform(len(x), f1, f2, m)
            f = linspace(f1, f2, m)
            y = t(x)
            plot(f, y)

        """
        if n < 1:
            raise ValueError("Invalid number of CZT data "
                             "points (%d) specified." % n)

        if m is None:
            m = n
        if np.size(fn) == 2:
            f1, f2 = fn
        elif np.size(fn) == 1:
            f1, f2 = 0.0, fn
        else:
            raise ValueError('fn must be a scalar or 2-length sequence')

        w = cmath.exp(-2j * pi * (f2-f1) / ((m-1)*Fs))
        a = cmath.exp(2j * pi * f1/Fs)
        CZT.__init__(self, n, m=m, w=w, a=a)
        self.f1, self.f2, self.Fs = f1, f2, Fs


class ScaledFFT(CZT):

    def __init__(self, n, m=None, scale=1.0):
        """
        Scaled FFT transform.

        Similar to FFT, where the frequency range is scaled and divided
        into m-1 equal steps.  Like the FFT, frequencies are arranged from
        0 to scale*Fs/2-delta followed by -scale*Fs/2 to -delta, where delta
        is the step size scale*Fs/m for sampling frequency Fs. The intended
        use is in a convolution of two signals, each has its own sampling step.

        This is equivalent to:

            fftshift(zoomfft(x, [-scale, scale*(m-2.)/m], m=m))

        For example:

            m, n = 10, len(x)
            sf = ScaledFFT(n, m=m, scale=0.25)
            X = fftshift(fft(x))
            W = linspace(-8, 8*(n-2.)/n, n)
            SX = fftshift(sf(x))
            SW = linspace(-2, 2*(m-2.)/m, m)
            plot(X, W, SX, SW)

        Parameters:
        ----------
        n : int
          Size of the signal
        m : int, optional
          The size of the output.
          Default: m=n
        scale : float, optional
          Frequency scaling factor.
          Default: scale=1.0

        Returns:
        -------
        callable f(x, axis=-1)
          function for computing the scaled FFT on x.
        """
        if n < 1:
            raise ValueError("Invalid number of CZT data "
                             "points (%d) specified." % n)

        if m is None:
            m = n
        w = np.exp(-2j * pi / m * scale)
        a = w**((m+1)//2)
        CZT.__init__(self, n=n, m=m, a=a, w=w)
        self.scale = scale

    def __call__(self, x, axis=-1):
        return fftshift(CZT.__call__(self, x, axis), axes=(axis,))
    __call__.__doc__ = CZT.__call__.__doc__


def scaledfft(x, m=None, scale=1.0, axis=-1):
    """
    Limited frequency FFT.

    See ScaledFFT doc for details

    Parameters:
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

    Returns:
    -------
        An array of the same rank of 'x', but with the size if
        the 'axis' dimension set to 'm'
    """
    x = np.asarray(x)
    transform = ScaledFFT(x.shape[axis], m, scale)
    return transform(x, axis)


def czt(x, m=None, w=None, a=1, factor=None, axis=-1):
    """
    Compute the frequency response around a spiral.

    Parameters:
    ----------
    x : array
        The set of data to transform.
    m : int, optional
        The number of points desired. Default is the length of the input data.
    w : complex, optional
        The ratio between points in each step.
    a : complex, optional
        The starting point in the complex plane.  Default is 1.
    factor : float, optional
        The frequency step scale (relative to the normal DFT frequency step).
        Cannot be specified at the same time as `w`.
    axis : int, optional
        Array dimension to operate over.  Default is the final dimension.

    Returns:
    -------
    An array of the same dimensions as x, but with the length of the
    transformed axis set to m.  Note that this is a view on a much
    larger array.  To save space, you may want to call it as
    y = ascontiguousarray(czt(x))

    See zoomfft for a friendlier interface to partial FFT calculations.

    If the transform needs to be repeated, use CZT to construct a
    specialized transform function which can be reused without
    recomputing constants.
    """
    x = np.asarray(x)
    transform = CZT(x.shape[axis], m=m, w=w, a=a, factor=factor)
    return transform(x, axis=axis)


def zoomfft(x, fn, m=None, Fs=2, axis=-1):
    """
    Compute the Fourier transform of x for frequencies in [f1, f2].

    Parameters:
    ----------
    x : array
        The input signal.
    fn : array_like, optional
        A scalar or length-2 sequence giving the frequency range. If a scalar,
        the range 0-fn is assumed.
    m : int, optional
        The number of points to evaluate.  The default is the length of x.
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
    array
        The transformed signal.  The Fourier transform will be calculate
        at the points f1, f1+df, f1+2df, ..., f2, where df=(f2-f1)/m.

    Notes
    -----
    ``zoomfft(x, [0, 2-2./len(x)])`` is equivalent to ``fft(x)``.

    To graph the magnitude of the resulting transform, use::

        plot(linspace(f1, f2, m), abs(zoomfft(x, [f1, f2], m)))

    If the transform needs to be repeated, use ZoomFFT to construct
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
