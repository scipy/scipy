"""The suite of window functions."""

import numpy as np
from scipy import special, linalg
from scipy.fftpack import fft

__all__ = ['boxcar', 'triang', 'parzen', 'bohman', 'blackman', 'nuttall',
           'blackmanharris', 'flattop', 'bartlett', 'hanning', 'barthann',
           'hamming', 'kaiser', 'gaussian', 'general_gaussian', 'chebwin',
           'slepian', 'hann', 'get_window']


def boxcar(M, sym=True):
    """The M-point boxcar window.

    """
    return np.ones(M, float)


def triang(M, sym=True):
    """The M-point triangular window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(1, int((M + 1) / 2) + 1)
    if M % 2 == 0:
        w = (2 * n - 1.0) / M
        w = np.r_[w, w[::-1]]
    else:
        w = 2 * n / (M + 1.0)
        w = np.r_[w, w[-2::-1]]

    if not sym and not odd:
        w = w[:-1]
    return w


def parzen(M, sym=True):
    """The M-point Parzen window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(-(M - 1) / 2.0, (M - 1) / 2.0 + 0.5, 1.0)
    na = np.extract(n < -(M - 1) / 4.0, n)
    nb = np.extract(abs(n) <= (M - 1) / 4.0, n)
    wa = 2 * (1 - np.abs(na) / (M / 2.0)) ** 3.0
    wb = (1 - 6 * (np.abs(nb) / (M / 2.0)) ** 2.0 +
          6 * (np.abs(nb) / (M / 2.0)) ** 3.0)
    w = np.r_[wa, wb, wa[::-1]]
    if not sym and not odd:
        w = w[:-1]
    return w


def bohman(M, sym=True):
    """The M-point Bohman window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    fac = np.abs(np.linspace(-1, 1, M)[1:-1])
    w = (1 - fac) * np.cos(np.pi * fac) + 1.0 / np.pi * np.sin(np.pi * fac)
    w = np.r_[0, w, 0]
    if not sym and not odd:
        w = w[:-1]
    return w


def blackman(M, sym=True):
    """The M-point Blackman window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M)
    w = (0.42 - 0.5 * np.cos(2.0 * np.pi * n / (M - 1)) +
         0.08 * np.cos(4.0 * np.pi * n / (M - 1)))
    if not sym and not odd:
        w = w[:-1]
    return w


def nuttall(M, sym=True):
    """A minimum 4-term Blackman-Harris window according to Nuttall.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    a = [0.3635819, 0.4891775, 0.1365995, 0.0106411]
    n = np.arange(0, M)
    fac = n * 2 * np.pi / (M - 1.0)
    w = (a[0] - a[1] * np.cos(fac) +
         a[2] * np.cos(2 * fac) - a[3] * np.cos(3 * fac))
    if not sym and not odd:
        w = w[:-1]
    return w


def blackmanharris(M, sym=True):
    """The M-point minimum 4-term Blackman-Harris window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    a = [0.35875, 0.48829, 0.14128, 0.01168]
    n = np.arange(0, M)
    fac = n * 2 * np.pi / (M - 1.0)
    w = (a[0] - a[1] * np.cos(fac) +
         a[2] * np.cos(2 * fac) - a[3] * np.cos(3 * fac))
    if not sym and not odd:
        w = w[:-1]
    return w


def flattop(M, sym=True):
    """The M-point Flat top window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    a = [0.2156, 0.4160, 0.2781, 0.0836, 0.0069]
    n = np.arange(0, M)
    fac = n * 2 * np.pi / (M - 1.0)
    w = (a[0] - a[1] * np.cos(fac) +
         a[2] * np.cos(2 * fac) - a[3] * np.cos(3 * fac) +
         a[4] * np.cos(4 * fac))
    if not sym and not odd:
        w = w[:-1]
    return w


def bartlett(M, sym=True):
    """The M-point Bartlett window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M)
    w = np.where(np.less_equal(n, (M - 1) / 2.0),
                 2.0 * n / (M - 1), 2.0 - 2.0 * n / (M - 1))
    if not sym and not odd:
        w = w[:-1]
    return w


def hanning(M, sym=True):
    """The M-point Hanning window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M)
    w = 0.5 - 0.5 * np.cos(2.0 * np.pi * n / (M - 1))
    if not sym and not odd:
        w = w[:-1]
    return w

hann = hanning


def barthann(M, sym=True):
    """Return the M-point modified Bartlett-Hann window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M)
    fac = np.abs(n / (M - 1.0) - 0.5)
    w = 0.62 - 0.48 * fac + 0.38 * np.cos(2 * np.pi * fac)
    if not sym and not odd:
        w = w[:-1]
    return w


def hamming(M, sym=True):
    """The M-point Hamming window.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M)
    w = 0.54 - 0.46 * np.cos(2.0 * np.pi * n / (M - 1))
    if not sym and not odd:
        w = w[:-1]
    return w


def kaiser(M, beta, sym=True):
    """Return a Kaiser window of length M with shape parameter beta.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M)
    alpha = (M - 1) / 2.0
    w = (special.i0(beta * np.sqrt(1 - ((n - alpha) / alpha) ** 2.0)) /
         special.i0(beta))
    if not sym and not odd:
        w = w[:-1]
    return w


def gaussian(M, std, sym=True):
    """Return a Gaussian window of length M with standard-deviation std.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M) - (M - 1.0) / 2.0
    sig2 = 2 * std * std
    w = np.exp(-n ** 2 / sig2)
    if not sym and not odd:
        w = w[:-1]
    return w


def general_gaussian(M, p, sig, sym=True):
    """Return a window with a generalized Gaussian shape.

    The Gaussian shape is defined as ``exp(-0.5*(x/sig)**(2*p))``, the
    half-power point is at ``(2*log(2)))**(1/(2*p)) * sig``.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1
    n = np.arange(0, M) - (M - 1.0) / 2.0
    w = np.exp(-0.5 * (n / sig) ** (2 * p))
    if not sym and not odd:
        w = w[:-1]
    return w


# `chebwin` contributed by Kumar Appaiah.

def chebwin(M, at, sym=True):
    """Dolph-Chebyshev window.

    Parameters
    ----------
    M : int
        Window size.
    at : float
        Attenuation (in dB).
    sym : bool
        Generates symmetric window if True.

    """
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')

    odd = M % 2
    if not sym and not odd:
        M = M + 1

    # compute the parameter beta
    order = M - 1.0
    beta = np.cosh(1.0 / order * np.arccosh(10 ** (np.abs(at) / 20.)))
    k = np.r_[0:M] * 1.0
    x = beta * np.cos(np.pi * k / M)
    # Find the window's DFT coefficients
    # Use analytic definition of Chebyshev polynomial instead of expansion
    # from scipy.special. Using the expansion in scipy.special leads to errors.
    p = np.zeros(x.shape)
    p[x > 1] = np.cosh(order * np.arccosh(x[x > 1]))
    p[x < -1] = (1 - 2 * (order % 2)) * np.cosh(order * np.arccosh(-x[x < -1]))
    p[np.abs(x) <= 1] = np.cos(order * np.arccos(x[np.abs(x) <= 1]))

    # Appropriate IDFT and filling up
    # depending on even/odd M
    if M % 2:
        w = np.real(fft(p))
        n = (M + 1) / 2
        w = w[:n] / w[0]
        w = np.concatenate((w[n - 1:0:-1], w))
    else:
        p = p * np.exp(1.j * np.pi / M * np.r_[0:M])
        w = np.real(fft(p))
        n = M / 2 + 1
        w = w / w[1]
        w = np.concatenate((w[n - 1:0:-1], w[1:n]))
    if not sym and not odd:
        w = w[:-1]
    return w


def slepian(M, width, sym=True):
    """Return the M-point slepian window.

    """
    if (M * width > 27.38):
        raise ValueError("Cannot reliably obtain slepian sequences for"
              " M*width > 27.38.")
    if M < 1:
        return np.array([])
    if M == 1:
        return np.ones(1, 'd')
    odd = M % 2
    if not sym and not odd:
        M = M + 1

    twoF = width / 2.0
    alpha = (M - 1) / 2.0
    m = np.arange(0, M) - alpha
    n = m[:, np.newaxis]
    k = m[np.newaxis, :]
    AF = twoF * special.sinc(twoF * (n - k))
    [lam, vec] = linalg.eig(AF)
    ind = np.argmax(abs(lam), axis=-1)
    w = np.abs(vec[:, ind])
    w = w / max(w)

    if not sym and not odd:
        w = w[:-1]
    return w


def get_window(window, Nx, fftbins=True):
    """
    Return a window of length `Nx` and type `window`.

    Parameters
    ----------
    window : string, float, or tuple
        The type of window to create. See below for more details.
    Nx : int
        The number of samples in the window.
    fftbins : bool, optional
        If True, create a "periodic" window ready to use with ifftshift
        and be multiplied by the result of an fft (SEE ALSO fftfreq).

    Notes
    -----
    Window types:

        boxcar, triang, blackman, hamming, hanning, bartlett,
        parzen, bohman, blackmanharris, nuttall, barthann,
        kaiser (needs beta), gaussian (needs std),
        general_gaussian (needs power, width),
        slepian (needs width), chebwin (needs attenuation)


    If the window requires no parameters, then `window` can be a string.

    If the window requires parameters, then `window` must be a tuple
    with the first argument the string name of the window, and the next
    arguments the needed parameters.

    If `window` is a floating point number, it is interpreted as the beta
    parameter of the kaiser window.

    Each of the window types listed above is also the name of
    a function that can be called directly to create a window of
    that type.

    Examples
    --------
    >>> get_window('triang', 7)
    array([ 0.25,  0.5 ,  0.75,  1.  ,  0.75,  0.5 ,  0.25])
    >>> get_window(('kaiser', 4.0), 9)
    array([ 0.08848053,  0.32578323,  0.63343178,  0.89640418,  1.        ,
            0.89640418,  0.63343178,  0.32578323,  0.08848053])
    >>> get_window(4.0, 9)
    array([ 0.08848053,  0.32578323,  0.63343178,  0.89640418,  1.        ,
            0.89640418,  0.63343178,  0.32578323,  0.08848053])

    """

    sym = not fftbins
    try:
        beta = float(window)
    except (TypeError, ValueError):
        args = ()
        if isinstance(window, tuple):
            winstr = window[0]
            if len(window) > 1:
                args = window[1:]
        elif isinstance(window, str):
            if window in ['kaiser', 'ksr', 'gaussian', 'gauss', 'gss',
                        'general gaussian', 'general_gaussian',
                        'general gauss', 'general_gauss', 'ggs',
                        'slepian', 'optimal', 'slep', 'dss',
                        'chebwin', 'cheb']:
                raise ValueError("The '" + window + "' window needs one or "
                                    "more parameters  -- pass a tuple.")
            else:
                winstr = window

        if winstr in ['blackman', 'black', 'blk']:
            winfunc = blackman
        elif winstr in ['triangle', 'triang', 'tri']:
            winfunc = triang
        elif winstr in ['hamming', 'hamm', 'ham']:
            winfunc = hamming
        elif winstr in ['bartlett', 'bart', 'brt']:
            winfunc = bartlett
        elif winstr in ['hanning', 'hann', 'han']:
            winfunc = hanning
        elif winstr in ['blackmanharris', 'blackharr', 'bkh']:
            winfunc = blackmanharris
        elif winstr in ['parzen', 'parz', 'par']:
            winfunc = parzen
        elif winstr in ['bohman', 'bman', 'bmn']:
            winfunc = bohman
        elif winstr in ['nuttall', 'nutl', 'nut']:
            winfunc = nuttall
        elif winstr in ['barthann', 'brthan', 'bth']:
            winfunc = barthann
        elif winstr in ['flattop', 'flat', 'flt']:
            winfunc = flattop
        elif winstr in ['kaiser', 'ksr']:
            winfunc = kaiser
        elif winstr in ['gaussian', 'gauss', 'gss']:
            winfunc = gaussian
        elif winstr in ['general gaussian', 'general_gaussian',
                        'general gauss', 'general_gauss', 'ggs']:
            winfunc = general_gaussian
        elif winstr in ['boxcar', 'box', 'ones']:
            winfunc = boxcar
        elif winstr in ['slepian', 'slep', 'optimal', 'dss']:
            winfunc = slepian
        elif winstr in ['chebwin', 'cheb']:
            winfunc = chebwin
        else:
            raise ValueError("Unknown window type.")

        params = (Nx,) + args + (sym,)
    else:
        winfunc = kaiser
        params = (Nx, beta, sym)

    return winfunc(*params)
