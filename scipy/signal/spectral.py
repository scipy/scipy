"""Tools for spectral analysis.
"""

import numpy as np
from scipy import fftpack
from windows import get_window
from _spectral import *

__all__ = ['periodogram', 'welch', 'lombscargle']


def periodogram(x, fs=1.0, window=None, nfft=None, sides='default', scaling='density', axis=-1):
    """
    Estimate power spectal density using a periodogram.  

    Parameters
    ----------
    x : array_like
        Time series of measurement values
        
    fs : float, optional 
        Sampling frequency of the `x` time series in units of Hz. Defaults
        to 1.0.

   window : string or tuple of string and parameter values or array_like
        Desired window to use. See scipy.signal.get_window for a list of
        windows and required parameters. If `window` is an array it will
        be used directly as the window. Defaults to 'None'; equivalent to 
        'boxcar'.

    nfft : int or None, optional 
        Length of the FFT used. If 'None' the length of `x` will be used.

    sides : { 'default', 'onesided', 'twosided' }, optional
        Selects which sides of the periodogram to return. 'Default'
        computes a one-sided spectrum for real data and a two-sided
        spectrum for complex data. 'onesided' or 'twosided' will force
        a one- or two-sided spectrum.  'onesided' cannot be specified with
        complex data.

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

    See Also
    --------
    welch: Estimate power spectral density using Welch's method 
    lombscargle: Lomb-Scargle periodogram for unevenly sampled data

    Examples
    --------
    >>> from scipy import signal

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
    >>> plt.ylabel('PSD [V/sqrt(Hz)]')

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

    The peak height in the power spectrum is an estimate of the RMS amplitude.

    >>> np.sqrt(Pxx_spec.max())
    2.0077340678640727 
    """
    x = np.asarray(x)

    if axis != -1:
        x = np.rollaxis(x, axis, len(x.shape))

    if nfft is None:
        nfft = x.shape[-1]

    if window is not None:
        if type(window) is str or type(window) is tuple:
            win = get_window(window, x.shape[-1])
        else:
            win = np.asarray(window)
            if len(win.shape) != 1:
                raise ValueError('window must be 1-D')
            elif win.shape[0] != x.shape[-1]:
                raise ValueError('window length must match data length.')
        x = win*x
        s1 = win.sum()**2
        s2 = (win*win).sum()
    else:
        s1 = nfft**2
        s2 = nfft

    if scaling == 'density':
        scale = 1.0 / (fs * s2)
    elif scaling == 'spectrum':
        scale = 1.0 / s1
    else:
        raise ValueError('Unknown scaling "{2}".'.format(scaling))

    if np.isrealobj(x) and (sides == 'default' or sides == 'onesided'):
        x = fftpack.rfft(x, nfft)
        outshape = list(x.shape)
        if nfft % 2 == 0: # even
            outshape[-1] = nfft/2+1
            Pxx = np.empty(outshape, x.dtype)
            Pxx[..., (0,-1)] = x[..., (0,-1)]**2
            Pxx[..., 1:-1] = x[..., 1:-1:2]**2 + x[..., 2::2]**2
        else: # odd
            outshape[-1] = (nfft+1)/2
            Pxx = np.empty(outshape, x.dtype)
            Pxx[..., 0] = x[..., 0]**2
            Pxx[..., 1:] = x[..., 1::2]**2 + x[..., 2::2]**2
        Pxx[..., 1:-1] *= 2*scale
        Pxx[..., (0,-1)] *= scale
        f = np.arange(Pxx.shape[-1])*(fs/nfft)
    elif (np.iscomplexobj(x) and sides == 'default') or sides == 'twosided':
        x = fftpack.fft(x, nfft)
        Pxx = (x * x.conj()).real
        f = fftpack.fftfreq(nfft, 1.0/fs)
        Pxx *= scale
    elif np.iscomplexobj(x) and sides == 'onesided':
        raise ValueError('Cannot compute onsided spectrum of complex data')
    else: 
        raise ValueError('Unknown sides: "{0}".'.format(sides))

    if axis != -1:
        Pxx = np.rollaxis(Pxx, len(Pxx.shape)-1, axis)
    return f, Pxx


def welch(x, fs=1.0, window='hanning', nfft=256, noverlap=None, sides='default', scaling='density', axis=-1):
    """
    Estimate power spectral density using Welch's method.  

    Welch's method [1]_ computes an estimate of the power spectral desnsity 
    by dividing the data into overlaping segments, computing a modified
    periodogram for each segment and averaging the periodograms.

    Parameters
    ----------
    x : array_like
        Time series of measurement values
        
    fs : float, optional 
        Sampling frequency of the `x` time series in units of Hz. Defaults
        to 1.0.

   window : string or tuple of string and parameter values or array_like
        Desired window to use. See scipy.signal.get_window for a list of
        windows and required parameters. If 'window' is array_like it will
        be used directly as the window and its length will be used for
        nfft. Defaults to 'hanning'.

    nfft : int,  optional 
        Length of the FFT used. Defaults to 256.

    noverlap: int or None, optional 
        Number of points to overlap between segments. If None, `noverlap`
        = nfft / 2.  Defaults to None.

    sides : { 'default', 'onesided', 'twosided' }, optional
        Selects which sides of the periodogram to return. 'default'
        computes a one-sided spectrum for real data and a two-sided
        spectrum for complex data. 'onesided' or 'twosided' will force
        a one- or two-sided spectrum.  'onesided' cannot be specified with
        complex data.

    scaling : { 'density', 'spectrum' }, optional 
        Selects between computing the power spectral density ('density')
        where Pxx has units of V**2/Hz if x is measured in V and computing
        the power spectrum ('spectrum') where Pxx has units of V**2 if x is
        measured in V. Defaults to 'density'

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
    the signal power, while not overcounting any of the data.  Narrower 
    windows may require a larger overlap.

    If `noverlap` is 0, this method is equivalent to Bartlett's method [2]_.

    References
    ----------
    .. [1] Welch, P. "The use of the fast Fourier transform for the
           estimation of power spectra: A method based on time averaging
           over short, modified periodograms", IEEE Trans. Audio 
           Electroacoust. vol. 15, pp. 70-73, 1967.
    .. [2] Bartlett, M. S. "Periodogram Analysis and Continuous Spectra",
           Biometrika, vol. 37, pp. 1-16, 1950.

    Examples
    --------
    >>> from scipy import signal

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

    >>> f, Pxx_den = signal.welch(x, fs, nfft=1024)
    >>> plt.semilogy(f, Pxx_den)
    >>> plt.xlabel('frequency [Hz]')
    >>> plt.ylabel('PSD [V/sqrt(Hz)]')

    If we average the last half of the spectral density, to exclude the
    peak, we can recover the noise power on the signal. 

    >>> np.mean(Pxx_den[256:])
    0.0009924865443739191

    Now compute and plot the power spectrum.

    >>> f, Pxx_spec = signal.welch(x, fs, 'flattop', nfft=1024, scaling='spectrum')
    >>> plt.figure()
    >>> plt.semilogy(f, np.sqrt(Pxx_spec))
    >>> plt.xlabel('frequency [Hz]')
    >>> plt.ylabel('Linear spectrum [V RMS]')

    The peak height in the power spectrum is an estimate of the RMS amplitude.

    >>> np.sqrt(Pxx_spec.max())
    2.0077340678640727 
    """
    x = np.asarray(x)

    if axis != -1:
        x = np.rollaxis(x, axis, len(x.shape))

    if type(window) is str or type(window) is tuple:
        win = get_window(window, nfft)
    else:
        win = np.asarray(window)
        if len(win.shape) != 1:
            raise ValueError('window must be 1-D')
        nfft = win.shape[0]

    if scaling == 'density':
        scale = 1.0 / (fs * (win*win).sum())
    elif scaling == 'spectrum':
        scale = 1.0 / win.sum()**2
    else:
        raise  ValueError('Unknown scaling "{0}".'.format(scaling))

    if noverlap is None:
        noverlap = nfft // 2

    step = nfft - noverlap
    indices = np.arange(0, x.shape[-1]-nfft+1, step)

    if np.isrealobj(x) and (sides == 'default' or sides == 'onesided'):
        outshape = list(x.shape)
        if nfft % 2 == 0: # even
            outshape[-1] = nfft/2+1
            Pxx = np.empty(outshape, x.dtype)
            for k, ind in enumerate(indices):
                xft = fftpack.rfft(x[..., ind:ind+nfft]*win, nfft)
                if k == 0:
                    Pxx[..., (0,-1)] = xft[..., (0,-1)]**2
                    Pxx[..., 1:-1] = xft[..., 1:-1:2]**2 + xft[..., 2::2]**2
                else:
                    Pxx *= k/(k+1.0)
                    Pxx[..., (0,-1)] += xft[..., (0,-1)]**2 / (k+1.0)
                    Pxx[..., 1:-1] += (xft[..., 1:-1:2]**2 + xft[..., 2::2]**2) / (k+1.0)
        else: # odd
            outshape[-1] = (nfft+1)/2
            Pxx = np.empty(outshape, x.dtype)
            for k, ind in enumerate(indices):
                xft = fftpack.rfft(x[..., ind:ind+nfft]*win, nfft)
                if k == 0:
                    Pxx[..., 0] = xft[..., 0]**2
                    Pxx[..., 1:] = xft[..., 1::2]**2 + xft[..., 2::2]**2
                else:
                    Pxx *= k/(k+1.0)
                    Pxx[..., 0] += xft[..., 0]**2 / (k+1)
                    Pxx[..., 1:] += (xft[..., 1::2]**2 + xft[..., 2::2]**2) / (k+1.0)
        Pxx[..., 1:-1] *= 2*scale
        Pxx[..., (0,-1)] *= scale
        f = np.arange(Pxx.shape[-1])*(fs/nfft)
    elif (np.iscomplexobj(x) and sides == 'default') or sides == 'twosided':
        for k, ind in enumerate(indices):
            xft = fftpack.fft(x[..., ind:ind+nfft]*win, nfft)
            if k == 0:
                Pxx = (xft * xft.conj()).real
            else:
                Pxx *= k/(k+1.0)
                Pxx += (xft * xft.conj()).real / (k+1.0)
        Pxx *= scale
        f = fftpack.fftfreq(nfft, 1.0/fs)
    elif np.iscomplexobj(x) and sides == 'onesided':
        raise ValueError('Cannot compute onsided spectrum of complex data')
    else: 
        raise ValueError('Unknown sides: "{0}".'.format(sides))

    if axis != -1:
        Pxx = np.rollaxis(Pxx, len(Pxx.shape)-1, axis)
    return f, Pxx

