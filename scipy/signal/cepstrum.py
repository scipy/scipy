"""Cepstrum.
"""

import numpy as np

__all__ = ['complex_cepstrum', 'real_cepstrum', 'inverse_complex_cepstrum',
           'minimum_phase']
 
def complex_cepstrum(x, N=None):
    """Compute the complex cepstrum of a real sequence.
        
    Parameters
    ----------
    x : ndarray
        Real sequence to compute complex cepstrum of.
    N : {None, int}, optional
        Length of the Fourier transform.
    
    Returns
    -------
    ceps : ndarray
        The complex cepstrum of the real data sequence `x` computed using the 
        Fourier transform.
    ndelay : int
        The amount of samples of circular delay added to `x`.
    
    The complex cepstrum is given by
    
    .. math:: c[n] = F^{-1}\\left{\\log_{10}{\\left(F{x[n}\\right)}\\right}
    
    where :math:`x_[n]` is the input signal and :math:`F` and :math:`F_{-1} 
    are respectively the forward and backward Fourier transform.

    See Also
    --------
    real_cepstrum: Compute the real cepstrum.
    inverse_complex_cepstrum: Compute the inverse complex cepstrum of a real sequence.

    
    Examples
    --------
    In the following example we use the cepstrum to determine the fundamental 
    frequency of a set of harmonics. There is a distinct peak at the quefrency 
    corresponding to the fundamental frequency. To be more precise, the peak 
    corresponds to the spacing between the harmonics.
    
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import complex_cepstrum
    
    >>> duration = 5.0
    >>> fs = 8000.0
    >>> samples = int(fs*duration)
    >>> t = np.arange(samples) / fs
    
    >>> fundamental = 100.0
    >>> harmonics = np.arange(1, 30) * fundamental
    >>> signal = np.sin(2.0*np.pi*harmonics[:,None]*t).sum(axis=0)
    >>> ceps, _ = complex_cepstrum(signal)
    
    >>> fig = plt.figure()
    >>> ax0 = fig.add_subplot(211)
    >>> ax0.plot(t, signal)
    >>> ax0.set_xlabel('time in seconds')
    >>> ax0.set_xlim(0.0, 0.05)
    >>> ax1 = fig.add_subplot(212)
    >>> ax1.plot(t, ceps)
    >>> ax1.set_xlabel('quefrency in seconds')
    >>> ax1.set_xlim(0.005, 0.015)
    >>> ax1.set_ylim(-5., +10.)
    
    References
    ----------
    .. [1] Wikipedia, "Cepstrum".
           http://en.wikipedia.org/wiki/Cepstrum
    .. [2] M.P. Norton and D.G. Karczub, D.G., 
           "Fundamentals of Noise and Vibration Analysis for Engineers", 2003.
    .. [3] B. P. Bogert, M. J. R. Healy, and J. W. Tukey: 
           "The Quefrency Alanysis of Time Series for Echoes: Cepstrum, Pseudo 
           Autocovariance, Cross-Cepstrum and Saphe Cracking". 
           Proceedings of the Symposium on Time Series Analysis
           Chapter 15, 209-243. New York: Wiley, 1963.
    
    """
    def _unwrap(phase):
        samples = len(phase)
        unwrapped = np.unwrap(phase)
        center = (samples+1)//2
        if samples == 1: 
            center = 0  
        ndelay = np.round(unwrapped[center]/np.pi)
        unwrapped -= np.pi * ndelay * np.arange(samples) / center
        return unwrapped, ndelay
        
    N = N if N is not None else x.shape[-1]
    spectrum = np.fft.fft(x, n=N)
    unwrapped_phase, ndelay = _unwrap(np.angle(spectrum))
    log_spectrum = np.log(np.abs(spectrum)) + 1j*unwrapped_phase
    ceps = np.fft.ifft(log_spectrum).real
    
    return ceps, ndelay


def real_cepstrum(x, N=None):
    """Compute the real cepstrum of a real sequence. 
    
    x : ndarray
        Real sequence to compute real cepstrum of.
    N : {None, int}, optional
        Length of the Fourier transform.
    
    Returns
    -------
    ceps: ndarray
        The real cepstrum.
    
    The real cepstrum is given by
    
    .. math:: c[n] = F^{-1}\\left{\\log_{10}{\\left|F{x[n}\\right|}\\right}
    
    where :math:`x_[n]` is the input signal and :math:`F` and :math:`F_{-1} 
    are respectively the forward and backward Fourier transform. Note that 
    contrary to the complex cepstrum the magnitude is taken of the spectrum.
    
    
    See Also
    --------
    complex_cepstrum: Compute the complex cepstrum of a real sequence.
    inverse_complex_cepstrum: Compute the inverse complex cepstrum of a real sequence.
    
    Examples
    --------
    >>> from scipy.signal import real_cepstrum
    
    
    References
    ----------
    .. [1] Wikipedia, "Cepstrum".
           http://en.wikipedia.org/wiki/Cepstrum
    
    """
    N = N if N is not None else x.shape[-1]
    spectrum = np.fft.fft(x, n=N)
    ceps = np.fft.ifft(np.log(np.abs(spectrum))).real
    
    return ceps


def inverse_complex_cepstrum(ceps, ndelay):
    """Compute the inverse complex cepstrum of a real sequence.
    
    ceps : ndarray
        Real sequence to compute inverse complex cepstrum of.
    ndelay: int
        The amount of samples of circular delay added to `x`.
    
    Returns
    -------
    x : ndarray
        The inverse complex cepstrum of the real sequence `ceps`.
    
    The inverse complex cepstrum is given by
    
    .. math:: x[n] = F^{-1}\\left{\\exp(F(c[n]))\\right}
    
    where :math:`c_[n]` is the input signal and :math:`F` and :math:`F_{-1} 
    are respectively the forward and backward Fourier transform.
    
    See Also
    --------
    complex_cepstrum: Compute the complex cepstrum of a real sequence.
    real_cepstrum: Compute the real cepstrum of a real sequence.
      
    Examples
    --------
    Taking the complex cepstrum and then the inverse complex cepstrum results 
    in the original sequence.
    
    >>> import numpy as np
    >>> from scipy.signal import inverse_complex_cepstrum
    >>> x = np.arange(10)
    >>> ceps, ndelay = complex_cepstrum(x)
    >>> y = inverse_complex_cepstrum(ceps, ndelay)
    >>> print(x)
    >>> print(y)

    References
    ----------
    .. [1] Wikipedia, "Cepstrum".
           http://en.wikipedia.org/wiki/Cepstrum
      
    """
    def _wrap(phase, ndelay):
        samples = phase.shape[-1]
        center = (samples+1)//2
        wrapped = phase + np.pi * ndelay * np.arange(samples) / center
        return wrapped
    
    log_spectrum = np.fft.fft(ceps)
    spectrum = np.exp(log_spectrum.real + 1j * _wrap(log_spectrum.imag, ndelay))
    x = np.fft.ifft(spectrum).real
    return x


def minimum_phase(x, N=None):
    """Compute the minimum phase reconstruction of a real sequence.
    
    x : ndarray
        Real sequence to compute the minimum phase reconstruction of.    
    N : {None, int}, optional
        Length of the Fourier transform.
    
    Compute the minimum phase reconstruction of a real sequence using the 
    real cepstrum.
    
    Returns
    -------
    m : ndarray
        The minimum phase reconstruction of the real sequence `x`.
    
    See Also
    --------
    real_cepstrum: Compute the real cepstrum.
    
    Examples
    --------
    >>> from scipy.signal import minimum_phase
    
    
    References
    ----------
    .. [1] Soo-Chang Pei, Huei-Shan Lin. Minimum-Phase FIR Filter Design Using 
           Real Cepstrum. IEEE TRANSACTIONS ON CIRCUITS AND SYSTEMS-II: 
           EXPRESS BRIEFS, VOL. 53, NO. 10, OCTOBER 2006
    
    """
    N = N if N is not None else x.shape[-1]
    ceps = real_cepstrum(x, N=N)
    odd = N % 2 
    window = np.concatenate(([1.0], 2.0*np.ones((N+odd)/2-1), 
                             np.ones(1-odd), np.zeros((N+odd)/2-1)))
    
    m = np.fft.ifft(np.exp(np.fft.fft(window*ceps))).real
    
    return m
