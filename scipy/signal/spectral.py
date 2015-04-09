"""Tools for spectral analysis.
"""

from __future__ import division, print_function, absolute_import

import numpy as np
from scipy import fftpack
from . import signaltools
from .windows import get_window
from ._spectral import lombscargle
import warnings

from scipy._lib.six import string_types

__all__ = ['periodogram', 'welch', 'lombscargle', 'csd', 'coherence']


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
    detrend : str or function or False, optional
        Specifies how to detrend `x` prior to computing the spectrum. If
        `detrend` is a string, it is passed as the ``type`` argument to
        `detrend`.  If it is a function, it should return a detrended array.
        If `detrend` is False, no detrending is done.  Defaults to 'constant'.
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
    >>> plt.ylim([1e-7, 1e2])
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
    >>> plt.ylim([1e-4, 1e1])
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
    noverlap : int, optional
        Number of points to overlap between segments. If None,
        ``noverlap = nperseg // 2``.  Defaults to None.
    nfft : int, optional
        Length of the FFT used, if a zero padded FFT is desired.  If None,
        the FFT length is `nperseg`. Defaults to None.
    detrend : str or function or False, optional
        Specifies how to detrend each segment. If `detrend` is a string,
        it is passed as the ``type`` argument to `detrend`.  If it is a
        function, it takes a segment and returns a detrended segment.
        If `detrend` is False, no detrending is done.  Defaults to 'constant'.
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

    >>> f, Pxx_den = signal.welch(x, fs, nperseg=1024)
    >>> plt.semilogy(f, Pxx_den)
    >>> plt.ylim([0.5e-3, 1])
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

    freqs, Pxx = csd(x, x, fs, window, nperseg, noverlap, nfft, detrend,
                     return_onesided, scaling, axis)

    return freqs, Pxx.real


def csd(x, y, fs=1.0, window='hanning', nperseg=256, noverlap=None, nfft=None,
        detrend='constant', return_onesided=True, scaling='density', axis=-1):
    """
    Estimate the cross power spectral density, Pxy, using Welch's method.

    Parameters
    ----------
    x : array_like
        Time series of measurement values
    y : array_like
        Time series of measurement values
    fs : float, optional
        Sampling frequency of the `x` and `y` time series in units of Hz.
        Defaults to 1.0.
    window : str or tuple or array_like, optional
        Desired window to use. See `get_window` for a list of windows and
        required parameters. If `window` is array_like it will be used
        directly as the window and its length will be used for nperseg.
        Defaults to 'hanning'.
    nperseg : int, optional
        Length of each segment.  Defaults to 256.
    noverlap: int, optional
        Number of points to overlap between segments. If None,
        ``noverlap = nperseg // 2``.  Defaults to None.
    nfft : int, optional
        Length of the FFT used, if a zero padded FFT is desired.  If None,
        the FFT length is `nperseg`. Defaults to None.
    detrend : str or function or False, optional
        Specifies how to detrend each segment. If `detrend` is a string,
        it is passed as the ``type`` argument to `detrend`.  If it is a
        function, it takes a segment and returns a detrended segment.
        If `detrend` is False, no detrending is done.  Defaults to 'constant'.
    return_onesided : bool, optional
        If True, return a one-sided spectrum for real data. If False return
        a two-sided spectrum. Note that for complex data, a two-sided
        spectrum is always returned.
    scaling : { 'density', 'spectrum' }, optional
        Selects between computing the power spectral density ('density')
        where Pxy has units of V**2/Hz if x and y are measured in V and
        computing the power spectrum ('spectrum') where Pxy has units of V**2
        if x and y are measured in V. Defaults to 'density'.
    axis : int, optional
        Axis along which the CSD is computed for both inputs; the default is
        over the last axis (i.e. ``axis=-1``).

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Pxy : ndarray
        Cross spectral density or cross power spectrum of x,y.

    See Also
    --------
    periodogram: Simple, optionally modified periodogram
    lombscargle: Lomb-Scargle periodogram for unevenly sampled data
    welch: Power spectral density by Welch's method. [Equivalent to csd(x,x)]
    coherence: Magnitude squared coherence by Welch's method.

    Notes
    --------
    By convention, Pxy is computed with the conjugate FFT of X multiplied by
    the FFT of Y.

    If the input series differ in length, the shorter series will be
    zero-padded to match.

    An appropriate amount of overlap will depend on the choice of window
    and on your requirements.  For the default 'hanning' window an
    overlap of 50\% is a reasonable trade off between accurately estimating
    the signal power, while not over counting any of the data.  Narrower
    windows may require a larger overlap.

    If `noverlap` is 0, this method is equivalent to Bartlett's method [2]_.

    References
    ----------
    .. [1] P. Welch, "The use of the fast Fourier transform for the
           estimation of power spectra: A method based on time averaging
           over short, modified periodograms", IEEE Trans. Audio
           Electroacoust. vol. 15, pp. 70-73, 1967.
    .. [2] M.S. Bartlett, "Periodogram Analysis and Continuous Spectra",
           Biometrika, vol. 37, pp. 1-16, 1950.
    """

    freqs, Pxy, _ = _spectral_helper(x, y, fs, window, nperseg, noverlap, nfft,
                                     detrend, return_onesided, scaling, axis,
                                     mode='psd')

    # Average over windows.
    if len(Pxy.shape) >= 2 and Pxy.size > 0:
        if Pxy.shape[-1] > 1:
            Pxy = Pxy.mean(axis=-1)
        else:
            Pxy = np.squeeze(Pxy, axis=-1)

    return freqs, Pxy


def coherence(x, y, fs=1.0, window='hanning', nperseg=256, noverlap=None,
              nfft=None, detrend='constant', axis=-1):
    """
    Estimate the magnitude squared coherence estimate, Cxy, of discrete-time
    signals X and Y using Welch's method.

    Cxy = abs(Pxy)**2/(Pxx*Pyy), where Pxx and Pyy are power spectral density
    estimates of X and Y, and Pxy is the cross spectral density estimate of X
    and Y.

    Parameters
    ----------
    x : array_like
        Time series of measurement values
    y : array_like
        Time series of measurement values
    fs : float, optional
        Sampling frequency of the `x` and `y` time series in units of Hz.
        Defaults to 1.0.
    window : str or tuple or array_like, optional
        Desired window to use. See `get_window` for a list of windows and
        required parameters. If `window` is array_like it will be used
        directly as the window and its length will be used for nperseg.
        Defaults to 'hanning'.
    nperseg : int, optional
        Length of each segment.  Defaults to 256.
    noverlap: int, optional
        Number of points to overlap between segments. If None,
        ``noverlap = nperseg // 2``.  Defaults to None.
    nfft : int, optional
        Length of the FFT used, if a zero padded FFT is desired.  If None,
        the FFT length is `nperseg`. Defaults to None.
    detrend : str or function or False, optional
        Specifies how to detrend each segment. If `detrend` is a string,
        it is passed as the ``type`` argument to `detrend`.  If it is a
        function, it takes a segment and returns a detrended segment.
        If `detrend` is False, no detrending is done.  Defaults to 'constant'.
    axis : int, optional
        Axis along which the CSD is computed for both inputs; the default is
        over the last axis (i.e. ``axis=-1``).

    Returns
    -------
    f : ndarray
        Array of sample frequencies.
    Cxy : ndarray
        Magnitude squared coherence of x and y.

    See Also
    --------
    periodogram: Simple, optionally modified periodogram
    lombscargle: Lomb-Scargle periodogram for unevenly sampled data
    welch: Power spectral density by Welch's method.
    csd: Cross spectral density by Welch's method.

    Notes
    --------
    An appropriate amount of overlap will depend on the choice of window
    and on your requirements.  For the default 'hanning' window an
    overlap of 50\% is a reasonable trade off between accurately estimating
    the signal power, while not over counting any of the data.  Narrower
    windows may require a larger overlap.
    """
    freqs, Pxx = welch(x, fs, window, nperseg, noverlap, nfft, detrend,
                       axis=axis)
    _, Pyy = welch(y, fs, window, nperseg, noverlap, nfft, detrend, axis=axis)
    _, Pxy = csd(x, y, fs, window, nperseg, noverlap, nfft, detrend, axis=axis)

    Cxy = np.abs(Pxy)**2 / Pxx / Pyy

    return freqs, Cxy


def _spectral_helper(x, y, fs=1.0, window='hanning', nperseg=256,
                    noverlap=None, nfft=None, detrend='constant',
                    return_onesided=True, scaling='spectrum', axis=-1,
                    mode='psd'):
    '''
    Calculate various forms of windowed FFTs for PSD, CSD, etc.

    This is a helper function that implements the commonality between the
    psd, csd, and spectrogram functions. It is not designed to be called
    externally. The windows are not averaged over; the result from each window
    is returned.

    Parameters
    ---------
    x : array_like
        Array or sequence containing the data to be analyzed.
    y : array_like
        Array or sequence containing the data to be analyzed. If this is
        the same object in memoery as x (i.e. _spectral_helper(x, x, ...)),
        the extra computations are spared.
    fs : float, optional
        Sampling frequency of the time series in units of Hz. Defaults
        to 1.0.
    window : str or tuple or array_like, optional
        Desired window to use. See `get_window` for a list of windows and
        required parameters. If `window` is array_like it will be used
        directly as the window and its length will be used for nperseg.
        Defaults to 'hanning'.
    nperseg : int, optional
        Length of each segment.  Defaults to 256.
    noverlap : int, optional
        Number of points to overlap between segments. If None,
        ``noverlap = nperseg // 2``.  Defaults to None.
    nfft : int, optional
        Length of the FFT used, if a zero padded FFT is desired.  If None,
        the FFT length is `nperseg`. Defaults to None.
    detrend : str or function or False, optional
        Specifies how to detrend each segment. If `detrend` is a string,
        it is passed as the ``type`` argument to `detrend`.  If it is a
        function, it takes a segment and returns a detrended segment.
        If `detrend` is False, no detrending is done.  Defaults to 'constant'.
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
    mode : str, optional
        Defines what kind of return values are expected. Options are ['psd',
        'complex', 'magnitude', 'angle', 'phase'].

    Returns
    -------
    freqs : ndarray
        Array of sample frequencies.
    result : ndarray
        Array of output data, contents dependant on *mode* kwarg.
    t : ndarray
        Array of times corresponding to each data segment

    References
    ----------
    stackoverflow: Rolling window for 1D arrays in Numpy?
    <http://stackoverflow.com/a/6811241>
    stackoverflow: Using strides for an efficient moving average filter
    <http://stackoverflow.com/a/4947453>

    Notes
    -----
    Adapted from matplotlib.mlab
    '''
    if mode not in ['psd', 'complex', 'magnitude', 'angle', 'phase']:
        raise ValueError("Unknown value for mode %s, must be one of: "
                         "'default', 'psd', 'complex', "
                         "'magnitude', 'angle', 'phase'" % mode)

    # If x and y are the same object we can save ourselves some computation.
    same_data = y is x

    if not same_data and mode != 'psd':
        raise ValueError("x and y must be equal if mode is not 'psd'")

    axis = int(axis)

    # Ensure we have np.arrays, get outdtype
    x = np.asarray(x)
    if not same_data:
        y = np.asarray(y)
        outdtype = np.result_type(x,y,np.complex64)
    else:
        outdtype = np.result_type(x,np.complex64)

    if not same_data:
        # Check if we can broadcast the outer axes together
        xouter = list(x.shape)
        youter = list(y.shape)
        xouter.pop(axis)
        youter.pop(axis)
        try:
            outershape = np.broadcast(np.empty(xouter), np.empty(youter)).shape
        except ValueError:
            raise ValueError('x and y cannot be broadcast together.')

    if same_data:
        if x.size == 0:
            return np.empty(x.shape), np.empty(x.shape), np.empty(x.shape)
    else:
        if x.size == 0 or y.size == 0:
            outshape = outershape + (min([x.shape[axis], y.shape[axis]]),)
            emptyout = np.rollaxis(np.empty(outshape), -1, axis)
            return emptyout, emptyout, emptyout

    if x.ndim > 1:
        if axis != -1:
            x = np.rollaxis(x, axis, len(x.shape))
            if not same_data and y.ndim > 1:
                y = np.rollaxis(y, axis, len(y.shape))

    # Check if x and y are the same length, zero-pad if neccesary
    if not same_data:
        if x.shape[-1] != y.shape[-1]:
            if x.shape[-1] < y.shape[-1]:
                pad_shape = list(x.shape)
                pad_shape[-1] = y.shape[-1] - x.shape[-1]
                x = np.concatenate((x, np.zeros(pad_shape)), -1)
            else:
                pad_shape = list(y.shape)
                pad_shape[-1] = x.shape[-1] - y.shape[-1]
                y = np.concatenate((y, np.zeros(pad_shape)), -1)

    # X and Y are same length now, can test nperseg with either
    if x.shape[-1] < nperseg:
        warnings.warn('nperseg = {0:d}, is greater than input length = {1:d}, '
                      'using nperseg = {1:d}'.format(nperseg, x.shape[-1]))
        nperseg = x.shape[-1]

    nperseg = int(nperseg)
    if nperseg < 1:
        raise ValueError('nperseg must be a positive integer')

    if nfft is None:
        nfft = nperseg
    elif nfft < nperseg:
        raise ValueError('nfft must be greater than or equal to nperseg.')
    else:
        nfft = int(nfft)

    if noverlap is None:
        noverlap = nperseg//2
    elif noverlap >= nperseg:
        raise ValueError('noverlap must be less than nperseg.')
    else:
        noverlap = int(noverlap)

    # Handle detrending and window functions
    if not detrend:
        detrend_func = lambda d: d
    elif not hasattr(detrend, '__call__'):
        detrend_func = lambda d: signaltools.detrend(d, type=detrend, axis=-1)
    elif axis != -1:
        # Wrap this function so that it receives a shape that it could
        # reasonably expect to receive.
        def detrend_func(d):
            d = np.rollaxis(d, -1, axis)
            d = detrend(d)
            return np.rollaxis(d, axis, len(d.shape))
    else:
        detrend_func = detrend

    if isinstance(window, string_types) or type(window) is tuple:
        win = get_window(window, nperseg)
    else:
        win = np.asarray(window)
        if len(win.shape) != 1:
            raise ValueError('window must be 1-D')
        if win.shape[0] != nperseg:
            raise ValueError('window must have length of nperseg')

    if np.result_type(win,np.complex64) != outdtype:
        win = win.astype(outdtype)

    if mode == 'psd':
        if scaling == 'density':
            scale = 1.0 / (fs * (win*win).sum())
        elif scaling == 'spectrum':
            scale = 1.0 / win.sum()**2
        else:
            raise ValueError('Unknown scaling: %r' % scaling)
    else:
        scale = 1

    if return_onesided is True:
        if np.iscomplexobj(x):
            sides = 'twosided'
        else:
            sides = 'onesided'
            if not same_data:
                if np.iscomplexobj(y):
                    sides = 'twosided'
    else:
        sides = 'twosided'

    if sides == 'twosided':
        num_freqs = nfft
    elif sides == 'onesided':
        if nperseg % 2:
            num_freqs = (nfft + 1)//2
        else:
            num_freqs = nfft//2 + 1

    # Perform the windowed FFTs
    result = _fft_helper(x, win, detrend_func, nperseg, noverlap, nfft)
    result = result[..., :num_freqs]
    freqs = fftpack.fftfreq(nfft, 1/fs)[:num_freqs]

    if not same_data:
        # All the same operations on the y data
        result_y = _fft_helper(y, win, detrend_func, nperseg, noverlap, nfft)
        result_y = result_y[..., :num_freqs]
        result = np.conjugate(result) * result_y
    elif mode == 'psd':
        result = np.conjugate(result) * result
    elif mode == 'magnitude':
        result = np.absolute(result)
    elif mode == 'angle' or mode == 'phase':
        result = np.angle(result)
    elif mode == 'complex':
        pass

    result *= scale
    if mode == 'psd' and sides == 'onesided':
        result[...,1:-1] *= 2

    t = np.arange(nfft/2, x.shape[-1] - nfft/2 + 1, nfft - noverlap)/float(fs)

    if sides != 'twosided' and not nperseg % 2:
        # get the last value correctly, it is negative otherwise
        freqs[-1] *= -1

    # we unwrap the phase here to handle the onesided vs. twosided case
    if mode == 'phase':
        result = np.unwrap(result, axis=-1)

    result = result.astype(outdtype)

    # All imaginary parts are zero anyways
    if same_data:
        result = result.real

    # Output is going to have new last axis for window index
    if axis != -1:
        # Specify as positive axis index
        if axis < 0:
            axis = len(result.shape)-1-axis

        # Roll frequency axis back to axis where the data came from
        result = np.rollaxis(result, -1, axis)
    else:
        # Make sure window/time index is last axis
        result = np.rollaxis(result, -1, -2)

    return freqs, result, t


def _fft_helper(x, win, detrend_func, nperseg, noverlap, nfft):
    '''
    Calculate windowed FFT, for internal use by scipy.signal._spectral_helper

    This is a helper function that does the main FFT calculation for
    _spectral helper. All input valdiation is performed there, and the data
    axis is assumed to be the last axis of x. It is not designed to be called
    externally. The windows are not averaged over; the result from each window
    is returned.

    Returns
    -------
    result : ndarray
        Array of FFT data

    References
    ----------
    stackoverflow: Repeat NumPy array without replicating data?
    <http://stackoverflow.com/a/5568169>

    Notes
    -----
    Adapted from matplotlib.mlab
    '''
    # Created strided array of data segments
    if nperseg == 1 and noverlap == 0:
        result = x[..., np.newaxis]
    else:
        step = nperseg - noverlap
        shape = x.shape[:-1]+((x.shape[-1]-noverlap)//step, nperseg)
        strides = x.strides[:-1]+(step*x.strides[-1], x.strides[-1])
        result = np.lib.stride_tricks.as_strided(x, shape=shape,
                                                 strides=strides)

    # Detrend each data segment individually
    result = detrend_func(result)

    # Apply window by multiplication
    result = win * result

    # Perform the fft. Acts on last axis by default. Zero-pads automatically
    result = fftpack.fft(result, n=nfft)

    return result
