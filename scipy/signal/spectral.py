"""Tools for spectral analysis.
"""

from __future__ import division, print_function, absolute_import

import numpy as np
from scipy import fftpack
from . import signaltools
from .windows import get_window
from ._spectral import lombscargle
import warnings

from scipy.lib.six import string_types

__all__ = ['periodogram', 'welch', 'lombscargle']


def periodogram(x, fs=1.0, window=None, nfft=None, detrend='constant',
                return_onesided=True, scaling='density', axis=-1):
    """
    Estimate power spectral density using a periodogram.

    Parameters
    ----------
    x : array_like
        Time series of measurement values
    fs : float, optional
        Sampling frequency of the `x` time series in units of Hz. Defaults
        to 1.0.
    window : str or tuple or array_like, optional
        Desired window to use. See `get_window` for a list of windows and
        required parameters. If `window` is an array it will be used
        directly as the window. Defaults to None; equivalent to 'boxcar'.
    nfft : int, optional
        Length of the FFT used. If None the length of `x` will be used.
    detrend : str or function, optional
        Specifies how to detrend `x` prior to computing the spectrum. If
        `detrend` is a string, it is passed as the ``type`` argument to
        `detrend`. If it is a function, it should return a detrended array.
        Defaults to 'constant'.
    return_onesided : bool, optional
        If True, return a one-sided spectrum for real data. If False return
        a two-sided spectrum. Note that for complex data, a two-sided
        spectrum is always returned.
    scaling : { 'density', 'spectrum' }, optional
        Selects between computing the power spectral density ('density')
        where `Pxx` has units of V**2/Hz if `x` is measured in V and computing
        the power spectrum ('spectrum') where `Pxx` has units of V**2 if `x` is
        measured in V. Defaults to 'density'
    axis : int, optional
        Axis along which the periodogram is computed; the default is over
        the last axis (i.e. ``axis=-1``).

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxx : ndarray
        Power spectral density or power spectrum of `x`.

    Notes
    -----
    .. versionadded:: 0.12.0

    See Also
    --------
    welch: Estimate power spectral density using Welch's method
    lombscargle: Lomb-Scargle periodogram for unevenly sampled data

    Examples
    --------
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    Generate a test signal, a 2 Vrms sine wave at 1234 Hz, corrupted by
    0.001 V**2/Hz of white noise sampled at 10 kHz.

    >>> fs = 10e3
    >>> N = 1e5
    >>> amp = 2*np.sqrt(2)
    >>> freq = 1234.0
    >>> noise_power = 0.001 * fs / 2
    >>> time = np.arange(N) / fs
    >>> x = amp*np.sin(2*np.pi*freq*time)
    >>> x += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)

    Compute and plot the power spectral density.

    >>> f, Pxx_den = signal.periodogram(x, fs)
    >>> plt.semilogy(f, Pxx_den)
    >>> plt.xlabel('frequency [Hz]')
    >>> plt.ylabel('PSD [V**2/Hz]')
    >>> plt.show()

    If we average the last half of the spectral density, to exclude the
    peak, we can recover the noise power on the signal.

    >>> np.mean(Pxx_den[256:])
    0.0009924865443739191

    Now compute and plot the power spectrum.

    >>> f, Pxx_spec = signal.periodogram(x, fs, 'flattop', scaling='spectrum')
    >>> plt.figure()
    >>> plt.semilogy(f, np.sqrt(Pxx_spec))
    >>> plt.xlabel('frequency [Hz]')
    >>> plt.ylabel('Linear spectrum [V RMS]')
    >>> plt.show()

    The peak height in the power spectrum is an estimate of the RMS amplitude.

    >>> np.sqrt(Pxx_spec.max())
    2.0077340678640727

    """
    x = np.asarray(x)

    if x.size == 0:
        return np.empty(x.shape), np.empty(x.shape)

    if window is None:
        window = 'boxcar'

    if nfft is None:
        nperseg = x.shape[axis]
    elif nfft == x.shape[axis]:
        nperseg = nfft
    elif nfft > x.shape[axis]:
        nperseg = x.shape[axis]
    elif nfft < x.shape[axis]:
        s = [np.s_[:]]*len(x.shape)
        s[axis] = np.s_[:nfft]
        x = x[s]
        nperseg = nfft
        nfft = None

    return welch(x, fs, window, nperseg, 0, nfft, detrend, return_onesided,
                 scaling, axis)


def welch(x, fs=1.0, window='hanning', nperseg=256, noverlap=None, nfft=None,
          detrend='constant', return_onesided=True, scaling='density', axis=-1):
    """
    Estimate power spectral density using Welch's method.

    Welch's method [1]_ computes an estimate of the power spectral density
    by dividing the data into overlapping segments, computing a modified
    periodogram for each segment and averaging the periodograms.

    Parameters
    ----------
    x : array_like
        Time series of measurement values
    fs : float, optional
        Sampling frequency of the `x` time series in units of Hz. Defaults
        to 1.0.
    window : str or tuple or array_like, optional
        Desired window to use. See `get_window` for a list of windows and
        required parameters. If `window` is array_like it will be used
        directly as the window and its length will be used for nperseg.
        Defaults to 'hanning'.
    nperseg : int, optional
        Length of each segment.  Defaults to 256.
    noverlap: int, optional
        Number of points to overlap between segments. If None,
        ``noverlap = nperseg / 2``.  Defaults to None.
    nfft : int, optional
        Length of the FFT used, if a zero padded FFT is desired.  If None,
        the FFT length is `nperseg`. Defaults to None.
    detrend : str or function, optional
        Specifies how to detrend each segment. If `detrend` is a string,
        it is passed as the ``type`` argument to `detrend`. If it is a
        function, it takes a segment and returns a detrended segment.
        Defaults to 'constant'.
    return_onesided : bool, optional
        If True, return a one-sided spectrum for real data. If False return
        a two-sided spectrum. Note that for complex data, a two-sided
        spectrum is always returned.
    scaling : { 'density', 'spectrum' }, optional
        Selects between computing the power spectral density ('density')
        where Pxx has units of V**2/Hz if x is measured in V and computing
        the power spectrum ('spectrum') where Pxx has units of V**2 if x is
        measured in V. Defaults to 'density'.
    axis : int, optional
        Axis along which the periodogram is computed; the default is over
        the last axis (i.e. ``axis=-1``).

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxx : ndarray
        Power spectral density or power spectrum of x.

    See Also
    --------
    periodogram: Simple, optionally modified periodogram
    lombscargle: Lomb-Scargle periodogram for unevenly sampled data

    Notes
    -----
    An appropriate amount of overlap will depend on the choice of window
    and on your requirements.  For the default 'hanning' window an
    overlap of 50% is a reasonable trade off between accurately estimating
    the signal power, while not over counting any of the data.  Narrower
    windows may require a larger overlap.

    If `noverlap` is 0, this method is equivalent to Bartlett's method [2]_.

    .. versionadded:: 0.12.0

    References
    ----------
    .. [1] P. Welch, "The use of the fast Fourier transform for the
           estimation of power spectra: A method based on time averaging
           over short, modified periodograms", IEEE Trans. Audio
           Electroacoust. vol. 15, pp. 70-73, 1967.
    .. [2] M.S. Bartlett, "Periodogram Analysis and Continuous Spectra",
           Biometrika, vol. 37, pp. 1-16, 1950.

    Examples
    --------
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    Generate a test signal, a 2 Vrms sine wave at 1234 Hz, corrupted by
    0.001 V**2/Hz of white noise sampled at 10 kHz.

    >>> fs = 10e3
    >>> N = 1e5
    >>> amp = 2*np.sqrt(2)
    >>> freq = 1234.0
    >>> noise_power = 0.001 * fs / 2
    >>> time = np.arange(N) / fs
    >>> x = amp*np.sin(2*np.pi*freq*time)
    >>> x += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)

    Compute and plot the power spectral density.

    >>> f, Pxx_den = signal.welch(x, fs, 1024)
    >>> plt.semilogy(f, Pxx_den)
    >>> plt.xlabel('frequency [Hz]')
    >>> plt.ylabel('PSD [V**2/Hz]')
    >>> plt.show()

    If we average the last half of the spectral density, to exclude the
    peak, we can recover the noise power on the signal.

    >>> np.mean(Pxx_den[256:])
    0.0009924865443739191

    Now compute and plot the power spectrum.

    >>> f, Pxx_spec = signal.welch(x, fs, 'flattop', 1024, scaling='spectrum')
    >>> plt.figure()
    >>> plt.semilogy(f, np.sqrt(Pxx_spec))
    >>> plt.xlabel('frequency [Hz]')
    >>> plt.ylabel('Linear spectrum [V RMS]')
    >>> plt.show()

    The peak height in the power spectrum is an estimate of the RMS amplitude.

    >>> np.sqrt(Pxx_spec.max())
    2.0077340678640727

    """
    x = np.asarray(x)

    if x.size == 0:
        return np.empty(x.shape), np.empty(x.shape)

    if axis != -1:
        x = np.rollaxis(x, axis, len(x.shape))

    if x.shape[-1] < nperseg:
        warnings.warn('nperseg = %d, is greater than x.shape[%d] = %d, using '
                      'nperseg = x.shape[%d]'
                      % (nperseg, axis, x.shape[axis], axis))
        nperseg = x.shape[-1]

    if isinstance(window, string_types) or type(window) is tuple:
        win = get_window(window, nperseg)
    else:
        win = np.asarray(window)
        if len(win.shape) != 1:
            raise ValueError('window must be 1-D')
        if  win.shape[0] > x.shape[-1]:
            raise ValueError('window is longer than x.')
        nperseg = win.shape[0]

    if scaling == 'density':
        scale = 1.0 / (fs * (win*win).sum())
    elif scaling == 'spectrum':
        scale = 1.0 / win.sum()**2
    else:
        raise ValueError('Unknown scaling: %r' % scaling)

    if noverlap is None:
        noverlap = nperseg // 2
    elif noverlap >= nperseg:
        raise ValueError('noverlap must be less than nperseg.')

    if nfft is None:
        nfft = nperseg
    elif nfft < nperseg:
        raise ValueError('nfft must be greater than or equal to nperseg.')

    if not hasattr(detrend, '__call__'):
        detrend_func = lambda seg: signaltools.detrend(seg, type=detrend)
    elif axis != -1:
        # Wrap this function so that it receives a shape that it could
        # reasonably expect to receive.
        def detrend_func(seg):
            seg = np.rollaxis(seg, -1, axis)
            seg = detrend(seg)
            return np.rollaxis(seg, axis, len(seg.shape))
    else:
        detrend_func = detrend

    step = nperseg - noverlap
    indices = np.arange(0, x.shape[-1]-nperseg+1, step)

    if np.isrealobj(x) and return_onesided:
        outshape = list(x.shape)
        if nfft % 2 == 0: # even
            outshape[-1] = nfft // 2 + 1
            Pxx = np.empty(outshape, x.dtype)
            for k, ind in enumerate(indices):
                x_dt = detrend_func(x[..., ind:ind+nperseg])
                xft = fftpack.rfft(x_dt*win, nfft)
                # fftpack.rfft returns the positive frequency part of the fft
                # as real values, packed r r i r i r i ...
                # this indexing is to extract the matching real and imaginary
                # parts, while also handling the pure real zero and nyquist
                # frequencies.
                if k == 0:
                    Pxx[..., (0,-1)] = xft[..., (0,-1)]**2
                    Pxx[..., 1:-1] = xft[..., 1:-1:2]**2 + xft[..., 2::2]**2
                else:
                    Pxx *= k/(k+1.0)
                    Pxx[..., (0,-1)] += xft[..., (0,-1)]**2 / (k+1.0)
                    Pxx[..., 1:-1] += (xft[..., 1:-1:2]**2 + xft[..., 2::2]**2) \
                                    / (k+1.0)
        else: # odd
            outshape[-1] = (nfft+1) // 2
            Pxx = np.empty(outshape, x.dtype)
            for k, ind in enumerate(indices):
                x_dt = detrend_func(x[..., ind:ind+nperseg])
                xft = fftpack.rfft(x_dt*win, nfft)
                if k == 0:
                    Pxx[..., 0] = xft[..., 0]**2
                    Pxx[..., 1:] = xft[..., 1::2]**2 + xft[..., 2::2]**2
                else:
                    Pxx *= k/(k+1.0)
                    Pxx[..., 0] += xft[..., 0]**2 / (k+1)
                    Pxx[..., 1:] += (xft[..., 1::2]**2 + xft[..., 2::2]**2) \
                                  / (k+1.0)

        Pxx[..., 1:-1] *= 2*scale
        Pxx[..., (0,-1)] *= scale
        f = np.arange(Pxx.shape[-1]) * (fs/nfft)
    else:
        for k, ind in enumerate(indices):
            x_dt = detrend_func(x[..., ind:ind+nperseg])
            xft = fftpack.fft(x_dt*win, nfft)
            if k == 0:
                Pxx = (xft * xft.conj()).real
            else:
                Pxx *= k/(k+1.0)
                Pxx += (xft * xft.conj()).real / (k+1.0)
        Pxx *= scale
        f = fftpack.fftfreq(nfft, 1.0/fs)

    if axis != -1:
        Pxx = np.rollaxis(Pxx, -1, axis)

    return f, Pxx

