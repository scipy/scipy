"""Tools based on the Wigner quasiprobability distributions.
"""

from __future__ import annotations
from typing import Tuple, TYPE_CHECKING

import numpy as np
from scipy.fft import fft

if TYPE_CHECKING:
    import numpy.typing as npt


__all__ = ['wigner_ville']


def wigner_ville(
        x: npt.ArrayLike,
        fs: float = 1.0,
        *,
        resolution: int = 1,
        window_size: None | int = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Wigner-Ville distribution.

    The Wigner–Ville distribution is a method of signal analysis that yields
    (under ideal conditions) the signal energy distributed over time and
    frequency. It provides the highest temporal vs frequency resolution
    possible within the limitations of the uncertainty principle.

    Parameters
    ----------
    x : array_like
        Time series of measurement values.
    fs : float, optional
        Sampling frequency of the `x` time series. Defaults to 1.0.
    resolution : int, optional
        Resolution, number of samples between windows. Defaults to 1.
    window_size : int, optional
        Window size, length of frequency axis. Defaults to ``len(x) // 2``.

    Returns
    -------
    t : ndarray
        Array of time samples.
    f : ndarray
        Array of frequency samples.
    wv : ndarray
        Wigner-Ville distribution of `x`.

    References
    ----------
    .. [1] S. Mallat, "A Wavelet Tour of Signal Processing (3rd Edition)",
        Academic Press, 2009.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from scipy import signal

    Generate a test signal with two superimposed sine wave at 1 Hz and 5 Hz.

    >>> f_s = 100
    >>> t = np.arange(0, 1024) / 100
    >>> x = np.sin(2*np.pi*1*t) + np.sin(2*np.pi*5*t)

    Compute the Wigner-Ville distribution. Note that the analystic signal is
    supplied to the function.

    >>> t, f, wv = signal.wigner_ville(signal.hilbert(x), f_s)

    Plot the Wigner-Ville distribution with expected cross-terms clearly
    visible at 3 Hz.

    >>> plt.figure()
    >>> plt.pcolormesh(t, f, wv, shading='nearest')
    >>> plt.xlabel('Time $t$ / s')
    >>> plt.ylabel('Frequency $f$ / Hz')
    >>> plt.ylim([0, 6])
    >>> plt.show()
    """
    x = np.asarray(x)
    if window_size is None:
        window_size = len(x) // 2

    n_pts = int(np.floor(np.floor(len(x) / resolution) / 2) * 2)
    odd_win_size = int((np.floor((window_size - 1) / 2) * 2)) + 1
    half_win_size = (odd_win_size + 1) // 2 - 1
    x_padded = np.concatenate((np.zeros(odd_win_size - 1), x, np.zeros(
        odd_win_size - 1)))

    wv = np.zeros((window_size, n_pts), dtype=float)
    r = np.zeros(window_size, dtype=complex)
    idx = np.arange(1, half_win_size + 1, dtype=int)

    for n in range(0, n_pts // 2):
        idy = 2 * n * resolution + odd_win_size - 1
        r[0] = x_padded[idy] * np.conj(x_padded[idy]) + 1j * x_padded[
            idy + resolution] * np.conj(x_padded[idy + resolution])
        v1 = x_padded[idy + idx] * np.conj(x_padded[idy - idx])
        v2 = x_padded[idy + resolution + idx] * np.conj(
            x_padded[idy + resolution - idx])
        r[idx] = v1 + 1j * v2
        r[window_size - idx] = np.conj(v1) + 1j * np.conj(v2)
        r_fft = fft(r, window_size)
        wv[:, 2 * n] = np.real(r_fft)
        wv[:, 2 * n + 1] = np.imag(r_fft)

    t = np.arange(0, n_pts) * resolution / fs
    f = np.arange(0, window_size) * (fs / 2) / (window_size - 1)
    return t, f, wv
