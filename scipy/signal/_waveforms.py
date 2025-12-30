# Author: Travis Oliphant
# 2003
#
# Feb. 2010: Updated by Warren Weckesser:
#   Rewrote much of chirp()
#   Added sweep_poly()
from typing import Literal

import numpy as np
from numpy import asarray, zeros, place, nan, mod, pi, extract, log, sqrt, \
    exp, cos, sin, polyval, polyint

from scipy._lib._array_api import array_namespace
from scipy.fft import rfft

__all__ = ['sawtooth', 'sawtooth_rfft', 'square', 'square_rfft', 'gausspulse', 'chirp',
           'sweep_poly', 'unit_impulse']


def sawtooth(t, width=1):
    """
    Return a periodic sawtooth or triangle waveform.

    The sawtooth waveform has a period ``2*pi``, rises from -1 to 1 on the
    interval 0 to ``width*2*pi``, then drops from 1 to -1 on the interval
    ``width*2*pi`` to ``2*pi``. `width` must be in the interval [0, 1].

    Note that the created signal is not band-limited. Consult the example of
    `sawtooth_rfft` on how to create a waveform bandlimited by the Nyquist frequency.

    Parameters
    ----------
    t : array_like
        Time.
    width : array_like, optional
        Width of the rising ramp as a proportion of the total cycle.
        Default is 1, producing a rising ramp, while 0 produces a falling
        ramp.  `width` = 0.5 produces a triangle wave.
        If an array, causes wave shape to change over time, and must be the
        same length as t.

    Returns
    -------
    y : ndarray
        Output array containing the sawtooth waveform.

    See Also
    --------
    sawtooth_rfft: Onesided FFT of a sawtooth or triangle waveform.
    square: Create periodic square-wave waveform.
    square_rfft: Onesided FFT of a square-wave waveform.

    Examples
    --------
    The following example depicts 5 Hz waveforms with various `width` parameters:

    >>> import numpy as np
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    ...
    >>> t = np.linspace(0, 1, 500)
    >>> fg0, axx = plt.subplots(5, 1, sharex='all', figsize=(5, 5), tight_layout=True)
    >>> axx[0].set_title("Triangle signals with various width parameters")
    >>> for i_, ax_ in enumerate(axx):
    ...     w_ = i_ / (len(axx) - 1)
    ...     y_ = signal.sawtooth(2 * np.pi * 5 * t, width=w_)
    ...     ax_.plot(t, y_, f'C{i_}', label=f"width={w_}")
    ...     ax_.legend(loc='upper right', framealpha=1)
    ...     ax_.grid(True)
    ...     ax_.set_ylabel("Ampl.")
    >>> axx[-1].set(xlabel="Time in seconds")
    >>> plt.show()
    """
    t, w = asarray(t), asarray(width)
    w = asarray(w + (t - t))
    t = asarray(t + (w - w))
    y = zeros(t.shape, dtype="d")

    # width must be between 0 and 1 inclusive
    mask1 = (w > 1) | (w < 0)
    place(y, mask1, nan)

    # take t modulo 2*pi
    tmod = mod(t, 2 * pi)

    # on the interval 0 to width*2*pi function is
    #  tmod / (pi*w) - 1
    mask2 = (1 - mask1) & (tmod < w * 2 * pi)
    tsub = extract(mask2, tmod)
    wsub = extract(mask2, w)
    place(y, mask2, tsub / (pi * wsub) - 1)

    # on the interval width*2*pi to 2*pi function is
    #  (pi*(w+1)-tmod) / (pi*(1-w))

    mask3 = (1 - mask1) & (1 - mask2)
    tsub = extract(mask3, tmod)
    wsub = extract(mask3, w)
    place(y, mask3, (pi * (wsub + 1) - tsub) / (pi * (1 - wsub)))
    return y


def sawtooth_rfft(n: int, m_cyc: int, duty: float = 0.5, *,
                  norm: Literal['backward', 'ortho', 'forward'] = 'backward'):
    r"""Onesided FFT of a sawtooth or triangle wave which is band-limited by the
    Nyquist frequency.

    This function returns the result of applying the :func:`~scipy.fft.rfft` function
    to a virtual sawtooth wave signal of `n` samples which had an ideal low-pass filter
    applied to remove spectral components above the Nyquist frequency ``0.5/T``. Here,
    ``T`` is assumed to be the sampling interval.

    Parameters
    ----------
    n : int
        Number of samples of the virtual input signal.
    m_cyc : int
        Number of full cycles for `n` samples (``0 < m_cyc < n // 2`` must hold).
        The frequency of the sawtooth wave is ``m_cyc / (n*T)``.
    duty : float
        Duty cycle, i.e., the ratio of the duration of the first part with positive
        slope divided by the duration of a complete cycle. Note that ``0 <= duty <= 1``
        must hold. Default: ``1``
    norm :  Literal['backward', 'ortho', 'forward']
        Scaling factor. Same as in :func:`~scipy.fft.rfft`. Default: ``'backward'``

    Returns
    -------
    X: np.ndarray
        A complex-valued array of length ``(n+1)//2`` containing the spectral values.
        ``scipy.fft.rfftfreq(n, T)`` determines the corresponding frequency values of
        `X`. The signal is calculated with ``scipy.fft.irfft(X, n, norm=norm)``.

    See Also
    --------
    sawtooth: Create a periodic sawtooth or triangle waveform.
    square: Create periodic square-wave waveform.
    square_rfft: Onesided FFT of a square-wave waveform.

    Notes
    -----
    A sawtooth wave :math:`x(t)` of frequency :math:`1/\tau` and duty cycle `d` can be
    expressed as the Fourier series

    .. math::

        x(t) = \sum_{l=-\infty}^\infty c[l] \exp\{\mathbb{j}2\pi l t /\tau\}\ ,\quad
        c[l] = \begin{cases}
               0 &\text{ for } l = 0\ ,\\
               -\frac{\mathbb{j}\tau}{\pi l} &\text{ for } l \neq 0 \wedge d = 0 \ ,\\
               +\frac{\mathbb{j}\tau}{\pi l} &\text{ for } l \neq 0 \wedge d = 1 \ ,\\
               \frac{\tau \left(1 - \operatorname{e}^{\mathbb{j}2\pi d l}\right)}{
                     2\pi^2 l^2 d (d-1)}
               &\text{ for } l \neq 0 \wedge 0 < d < 1 \ .
            \end{cases}

    The close relationship between Fourier series and the discrete Fourier transform
    allows to utilize :math:`c[l]` to obtain a closed-form expression for `X`. I.e.,
    ``X[m_cyc*l] = s * c[l] / tau``, with ``tau = n * T / m_cyc`` being the duration of
    one cycle. The scaling factor ``s`` has a value of ``n`` for parameter
    `norm` = ``'backward'``, ``sqrt(n)`` for `norm` = ``'ortho'``, or ``1`` for
    `norm` = ``'forward'``. All other values of `X` are zero.

    Note that the magnitude :math:`|c[l]|` decreases proportionally to :math:`1/l^2`
    if :math:`0 < d < 1` holds, else proportionally to :math:`1/|l|`. Hence, an ideal
    sawtooth wave is not band-limited. The function :func:`~scipy.signal.sawtooth` can
    be used to sample non-band-limited sawtooth waves.

    Examples
    --------
    The following code generates a plot of 1 Hz square wave and its Magnitude spectrum
    (Consult the `square_rfft` examples on creating a signal with varying frequency):

    >>> import numpy as np
    >>> from matplotlib import pyplot as plt
    >>> from scipy.fft import irfft, rfftfreq
    >>> from scipy.signal import sawtooth_rfft
    ...
    >>> n, T = 500, 1e-2  # number of samples and sampling interval
    >>> m_cyc, d = 5, 1  # number of full cycles and duty cycle
    ...
    >>> X = sawtooth_rfft(n, m_cyc, duty=d)
    >>> x, t = irfft(X, n), np.arange(n) * T  # signal and time stamps
    ...
    >>> f = rfftfreq(n, T)  # frequencies of X
    >>> df = f[1] - f[0]  # frequency resolution
    ...
    >>> fg, axx = plt.subplots(2, 1, tight_layout=True)
    >>> axx[0].set_title(f"{m_cyc*df:g} Hz Sawtooth Wave ({100*d:g}% duty cycle)")
    >>> axx[0].set(xlabel=rf"Time $t$ in seconds ($T={1e3*T:g}\,$ms, ${n=}$ samples)",
    ...            ylabel="Amplitude $x(t)$", xlim=t[[0, -1]])
    >>> axx[1].set_title(f"Magnitude Spectrum of $x(t)$")
    >>> axx[1].set(ylabel="Magnitude $|X(f)| / (nT)$", xlim=f[[0, -1]],
    ...            xlabel=rf"Frequency $f$ in hertz ($\Delta f={df:g}\,$Hz)")
    >>> axx[0].plot(t, x, 'C0-')
    >>> axx[1].plot(f, abs(X)/n, 'C1-')
    ...
    >>> for ax_ in axx:
    ...    ax_.grid()
    >>> plt.show()
    """
    xp = array_namespace()
    if not (xp.result_type(n, xp.int64) == xp.int64 and n > 0):
        raise ValueError(f"Parameter {n=} is not a positive integer!")
    if not (xp.result_type(m_cyc, xp.int64) == xp.int64 and 0 < m_cyc < n // 2):
        raise ValueError(f"Parameter {m_cyc=} is not a positive integer < {n//2=}!")
    if not (0 <= duty <= 1):
        raise ValueError(f"0 <= duty <= 1 does not hold for parameter {duty=}!")

    scale_factors = {'backward': n, 'ortho': xp.sqrt(n), 'forward': 1}
    if norm not in scale_factors:
        raise ValueError(f"Parameter {norm=} not in {list(scale_factors)}!")

    n_X = n // 2 + 1  # number of bins produced by rfft for signal of n samples
    X_dtype = rfft([0., ], norm=norm).dtype  # Use same dtype as rfft produces
    X = xp.zeros((n_X,), dtype=X_dtype)

    ll = xp.arange(1, (n_X + m_cyc - 1) // m_cyc) * xp.pi
    if duty == 0:
        X[m_cyc::m_cyc] = -scale_factors[norm] * 1j / ll
    elif duty == 1:
        X[m_cyc::m_cyc] = scale_factors[norm] * 1j / ll
    else:
        exponent, denominator = xp.exp(-2j * duty * ll), 2 * ll**2 * duty * (duty - 1)
        X[m_cyc::m_cyc] = scale_factors[norm] * (1 - exponent) / denominator
    return X

def square(t, duty=0.5):
    """
    Return a periodic square-wave waveform.

    The square wave has a period ``2*pi``, has value +1 from 0 to
    ``2*pi*duty`` and -1 from ``2*pi*duty`` to ``2*pi``. `duty` must be in
    the interval [0,1].

    Note that the created signal is not band-limited. Consult the example of
    `square_rfft` on how to create a waveform bandlimited by the Nyquist frequency.

    Parameters
    ----------
    t : array_like
        The input time array.
    duty : array_like, optional
        Duty cycle.  Default is 0.5 (50% duty cycle).
        If an array, causes wave shape to change over time, and must be the
        same length as t.

    Returns
    -------
    y : ndarray
        Output array containing the square waveform.

    See Also
    --------
    square_rfft: Onesided FFT of a square-wave waveform.
    sawtooth: Create a periodic sawtooth or triangle waveform.
    sawtooth_rfft: Onesided FFT of a sawtooth or triangle waveform.

    Examples
    --------
    The following example shows a 5 Hz waveform sampled at 500 Hz for 1 second:

    >>> import numpy as np
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    ...
    >>> t = np.linspace(0, 1, 500, endpoint=False)
    >>> x = signal.square(2 * np.pi * 5 * t)
    >>>
    ...
    >>> fg0, ax0 = plt.subplots(tight_layout=True)
    >>> ax0.set_title("5 Hz Square Wave signal")
    >>> ax0.set(xlabel="Time in seconds", ylabel="Amplitude")
    >>> ax0.plot(t, x)
    >>> plt.show()


    The second example shows how to generate a pulse-width modulated square wave
    signal:

    >>> import numpy as np
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    ...
    >>> t = np.linspace(0, 1, 500, endpoint=False)
    >>> sig = np.sin(2 * np.pi * t)
    >>> pwm = signal.square(2 * np.pi * 30 * t, duty=(sig + 1)/2)
    ...
    >>> fg1, (ax10, ax11) = plt.subplots(2, 1, sharex='all', tight_layout=True)
    >>> ax10.set_title("1 Hz Sine wave")
    >>> ax10.set_ylabel("Amplitude")
    >>> ax10.plot(t, sig, 'C1-')
    >>> ax11.set_title("Pulse-width modulated Square Wave")
    >>> ax11.set(xlabel="Time in seconds", ylabel="Amplitude")
    >>> ax11.plot(t, pwm, 'C0-')
    >>> plt.show()
    """
    t, w = asarray(t), asarray(duty)
    w = asarray(w + (t - t))
    t = asarray(t + (w - w))
    y = zeros(t.shape, dtype="d")

    # width must be between 0 and 1 inclusive
    mask1 = (w > 1) | (w < 0)
    place(y, mask1, nan)

    # on the interval 0 to duty*2*pi function is 1
    tmod = mod(t, 2 * pi)
    mask2 = (1 - mask1) & (tmod < w * 2 * pi)
    place(y, mask2, 1)

    # on the interval duty*2*pi to 2*pi function is
    #  (pi*(w+1)-tmod) / (pi*(1-w))
    mask3 = (1 - mask1) & (1 - mask2)
    place(y, mask3, -1)
    return y


def square_rfft(n: int, m_cyc: int, duty: float = 0.5,  *,
                norm: Literal['backward', 'ortho', 'forward'] = 'backward'):
    r"""Onesided FFT of a square wave which is band-limited by the Nyquist frequency.

    This function returns the result of applying the :func:`~scipy.fft.rfft` function
    to a virtual square wave signal of `n` samples which had an ideal low-pass filter
    applied to remove spectral components above the Nyquist frequency ``0.5/T``. Here,
    ``T`` is assumed to be the sampling interval.

    Parameters
    ----------
    n : int
        Number of samples of the virtual input signal.
    m_cyc : int
        Number of full cycles for `n` samples (``0 < m_cyc < n // 2`` must hold).
        The frequency of the square wave is ``m_cyc / (n*T)``.
    duty : float
        Duty cycle, i.e., the ratio of the duration of the positive part divided by the
        duration of a complete cycle. Note that ``0 < duty < 1`` must hold.
        Default: ``0.5``
    norm :  Literal['backward', 'ortho', 'forward']
        Scaling factor. Same as in :func:`~scipy.fft.rfft`. Default: ``'backward'``

    Returns
    -------
    X: np.ndarray
        A complex-valued array of length ``(n+1)//2`` containing the spectral values.
        ``scipy.fft.rfftfreq(n, T)`` determines the corresponding frequency values of
        `X`. The signal is calculated with ``scipy.fft.irfft(X, n, norm=norm)``.

    .. versionadded:: 1.16.0

    See Also
    --------
    square: Create periodic square-wave waveform.
    sawtooth: Create a periodic sawtooth or triangle waveform.
    sawtooth_rfft: Onesided FFT of a sawtooth or triangle waveform.


    Notes
    -----
    A square wave :math:`x(t)` of frequency :math:`1/\tau` and duty cycle `d` can be
    expressed as the Fourier series

    .. math::

        x(t) = \sum_{l=-\infty}^\infty c[l] \exp\{\mathbb{j}2\pi l t /\tau\}\ ,\quad
        c[l] = \begin{cases}
               \tau (2 d - 1) &\text{ for } l = 0\ ,\\
               \frac{\mathbb{j}\tau}{\pi l}
               \left(\operatorname{e}^{\mathbb{j}2\pi d l} - 1 \right)
               &\text{ for } l \neq 0\ .
            \end{cases}

    The close relationship between Fourier series and the discrete Fourier transform
    allows to utilize :math:`c[l]` to obtain a closed-form expression for `X`. I.e.,
    ``X[m_cyc*l] = s * c[l] / tau``, with ``tau = n * T / m_cyc`` being the duration of
    one cycle. The scaling factor ``s`` has a value of ``n`` for parameter
    `norm` = ``'backward'``, ``sqrt(n)`` for `norm` = ``'ortho'``, or ``1`` for
    `norm` = ``'forward'``. All other values of `X` are zero.

    Note that the magnitude :math:`|c[l]|` decreases proportionally to :math:`1/|l|`.
    Hence, an ideal square wave is not band-limited. The function
    :func:`~scipy.signal.square` can be used to sample non-band-limited square waves.

    Examples
    --------
    The following code generates a plot of 1 Hz square wave and its Magnitude
    spectrum:

    >>> import numpy as np
    >>> from matplotlib import pyplot as plt
    >>> from scipy.fft import irfft, rfftfreq
    >>> from scipy.signal import square_rfft
    ...
    >>> n, T = 500, 1e-2  # number of samples and sampling interval
    >>> m_cyc, d = 5, 0.5  # number of full cycles and duty cycle
    ...
    >>> X = square_rfft(n, m_cyc, duty=d)
    >>> x, t = irfft(X, n), np.arange(n) * T  # signal and time stamps
    ...
    >>> f = rfftfreq(n, T)  # frequencies of X
    >>> df = f[1] - f[0]  # frequency resolution
    ...
    >>> fg, axx = plt.subplots(2, 1, tight_layout=True)
    >>> axx[0].set_title(f"{m_cyc*df:g} Hz Square Wave ({100*d:g}% duty cycle)")
    >>> axx[0].set(xlabel=rf"Time $t$ in seconds ($T={1e3*T:g}\,$ms, ${n=}$ samples)",
    ...            ylabel="Amplitude", xlim=t[[0, -1]])
    >>> axx[1].set_title(f"Magnitude Spectrum")
    >>> axx[1].set(ylabel="Magnitude", xlim=f[[0, -1]],
    ...            xlabel=rf"Frequency $f$ in hertz ($\Delta f={df:g}\,$Hz)")
    >>> axx[0].plot(t, x, 'C0-')
    >>> axx[1].plot(f, abs(X)/n, 'C1-')
    ...
    >>> for ax_ in axx:
    ...     ax_.grid()
    >>> plt.show()

    The second example illustrates how to create a square wave with varying frequency.
    This is achieved by creating the signal directly from the Fourier series with the
    caveat that the complex exponential :math:`\exp\{\mathbb{j}2\pi l t /\tau\}` is
    replaced with a sinuoisod function with varying in frequency.

    >>> from matplotlib import pyplot as plt
    >>> from scipy.signal import chirp
    >>> from scipy.signal import square_rfft
    ...
    >>> n, T = 500, 1/500  # number of samples and sampling interval
    >>> m = 20  # number of Fourier coefficients
    >>> f0, f1 = 3, 10  # start and stop frequency of sweep
    ...
    >>> X = square_rfft(m, 1, duty=0.5, norm='forward') * n*T  # Fourier coefficients
    >>> t = np.arange(n) * T  # timestamps
    >>> x = np.ones(n) * X[0].real  # signal of zeroth Fourier coefficient
    >>> for l_, X_ in enumerate(X[1:], start=1):
    ...     x_ =  X_* chirp(t, f0*l_, n*T, f1*l_, method='linear', complex=True)
    ...     x += (x_ + np.conj(x_)).real  # since x is real, X[-l] == conj(X[l]) holds
    ...
    >>> fg0, ax0 = plt.subplots()
    >>> ax0.set_title(f"{m} Coefficient Fourier Series sweeping from " +
    ...               rf"${f0}\,$Hz to ${f1}\,$Hz ")
    >>> ax0.set(xlabel=rf"Time $t$ in seconds ({n} samples, $T={T*1e3:g}\,$ms)",
    ...         ylabel="Amplitude", xlim=(0, t[-1]))
    >>> ax0.plot(np.arange(n) * T, x)
    >>> ax0.grid()
    >>> plt.show()
    """
    xp = array_namespace()
    if not (xp.result_type(n, xp.int64) == xp.int64 and n > 0):
        raise ValueError(f"Parameter {n=} is not a positive integer!")
    if not (xp.result_type(m_cyc, xp.int64) == xp.int64 and 0 < m_cyc < n // 2):
        raise ValueError(f"Parameter {m_cyc=} is not a positive integer < {n//2=}!")
    if not (0 < duty < 1):
        raise ValueError(f"0 < duty < 1 does not hold for parameter {duty=}!")

    scale_factors = {'backward': n, 'ortho': xp.sqrt(n), 'forward': 1}
    if norm not in scale_factors:
        raise ValueError(f"Parameter {norm=} not in {list(scale_factors)}!")

    n_X = n // 2 + 1  # number of bins produced by rfft for signal of n samples
    X_dtype = rfft([0., ], norm=norm).dtype  # Use same dtype as rfft produces
    X = xp.zeros((n_X,), dtype=X_dtype)

    fac0, fac1 = scale_factors[norm], scale_factors[norm] * 1j / xp.pi
    X[0] = fac0 * (2 * duty - 1)  # zeroth Fourier coefficient
    ll = xp.arange(1, (n_X + m_cyc - 1) // m_cyc)
    X[m_cyc::m_cyc] = fac1 * (xp.exp(-2j * xp.pi * duty * ll) - 1) / ll
    return X



def gausspulse(t, fc=1000, bw=0.5, bwr=-6, tpr=-60, retquad=False,
               retenv=False):
    """
    Return a Gaussian modulated sinusoid:

        ``exp(-a t^2) exp(1j*2*pi*fc*t).``

    If `retquad` is True, then return the real and imaginary parts
    (in-phase and quadrature).
    If `retenv` is True, then return the envelope (unmodulated signal).
    Otherwise, return the real part of the modulated sinusoid.

    Parameters
    ----------
    t : ndarray or the string 'cutoff'
        Input array.
    fc : float, optional
        Center frequency (e.g. Hz).  Default is 1000.
    bw : float, optional
        Fractional bandwidth in frequency domain of pulse (e.g. Hz).
        Default is 0.5.
    bwr : float, optional
        Reference level at which fractional bandwidth is calculated (dB).
        Default is -6.
    tpr : float, optional
        If `t` is 'cutoff', then the function returns the cutoff
        time for when the pulse amplitude falls below `tpr` (in dB).
        Default is -60.
    retquad : bool, optional
        If True, return the quadrature (imaginary) as well as the real part
        of the signal.  Default is False.
    retenv : bool, optional
        If True, return the envelope of the signal.  Default is False.

    Returns
    -------
    yI : ndarray
        Real part of signal.  Always returned.
    yQ : ndarray
        Imaginary part of signal.  Only returned if `retquad` is True.
    yenv : ndarray
        Envelope of signal.  Only returned if `retenv` is True.

    Examples
    --------
    Plot real component, imaginary component, and envelope for a 5 Hz pulse,
    sampled at 100 Hz for 2 seconds:

    >>> import numpy as np
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    >>> t = np.linspace(-1, 1, 2 * 100, endpoint=False)
    >>> i, q, e = signal.gausspulse(t, fc=5, retquad=True, retenv=True)
    >>> plt.plot(t, i, t, q, t, e, '--')

    """
    if fc < 0:
        raise ValueError(f"Center frequency (fc={fc:.2f}) must be >=0.")
    if bw <= 0:
        raise ValueError(f"Fractional bandwidth (bw={bw:.2f}) must be > 0.")
    if bwr >= 0:
        raise ValueError(f"Reference level for bandwidth (bwr={bwr:.2f}) "
                         "must be < 0 dB")

    # exp(-a t^2) <->  sqrt(pi/a) exp(-pi^2/a * f^2)  = g(f)

    ref = pow(10.0, bwr / 20.0)
    # fdel = fc*bw/2:  g(fdel) = ref --- solve this for a
    #
    # pi^2/a * fc^2 * bw^2 /4=-log(ref)
    a = -(pi * fc * bw) ** 2 / (4.0 * log(ref))

    if isinstance(t, str):
        if t == 'cutoff':  # compute cut_off point
            #  Solve exp(-a tc**2) = tref  for tc
            #   tc = sqrt(-log(tref) / a) where tref = 10^(tpr/20)
            if tpr >= 0:
                raise ValueError("Reference level for time cutoff must "
                                 "be < 0 dB")
            tref = pow(10.0, tpr / 20.0)
            return sqrt(-log(tref) / a)
        else:
            raise ValueError("If `t` is a string, it must be 'cutoff'")

    yenv = exp(-a * t * t)
    yI = yenv * cos(2 * pi * fc * t)
    yQ = yenv * sin(2 * pi * fc * t)
    if not retquad and not retenv:
        return yI
    if not retquad and retenv:
        return yI, yenv
    if retquad and not retenv:
        return yI, yQ
    return yI, yQ, yenv  # if retquad and retenv


def chirp(t, f0, t1, f1, method='linear', phi=0, vertex_zero=True, *,
          complex=False):
    r"""Frequency-swept cosine generator.

    In the following, 'Hz' should be interpreted as 'cycles per unit';
    there is no requirement here that the unit is one second.  The
    important distinction is that the units of rotation are cycles, not
    radians. Likewise, `t` could be a measurement of space instead of time.

    Parameters
    ----------
    t : array_like
        Times at which to evaluate the waveform.
    f0 : float
        Frequency (e.g. Hz) at time t=0.
    t1 : float
        Time at which `f1` is specified.
    f1 : float
        Frequency (e.g. Hz) of the waveform at time `t1`.
    method : {'linear', 'quadratic', 'logarithmic', 'hyperbolic'}, optional
        Kind of frequency sweep.  If not given, `linear` is assumed.  See
        Notes below for more details.
    phi : float, optional
        Phase offset, in degrees. Default is 0.
    vertex_zero : bool, optional
        This parameter is only used when `method` is 'quadratic'.
        It determines whether the vertex of the parabola that is the graph
        of the frequency is at t=0 or t=t1.
    complex : bool, optional
        This parameter creates a complex-valued analytic signal instead of a
        real-valued signal. It allows the use of complex baseband (in communications
        domain). Default is False.

        .. versionadded:: 1.15.0

    Returns
    -------
    y : ndarray
        A numpy array containing the signal evaluated at `t` with the requested
        time-varying frequency.  More precisely, the function returns
        ``exp(1j*phase + 1j*(pi/180)*phi) if complex else cos(phase + (pi/180)*phi)``
        where `phase` is the integral (from 0 to `t`) of ``2*pi*f(t)``.
        The instantaneous frequency ``f(t)`` is defined below.

    See Also
    --------
    sweep_poly

    Notes
    -----
    There are four possible options for the parameter `method`, which have a (long)
    standard form and some allowed abbreviations. The formulas for the instantaneous
    frequency :math:`f(t)` of the generated signal are as follows:

    1. Parameter `method` in ``('linear', 'lin', 'li')``:

       .. math::
           f(t) = f_0 + \beta\, t           \quad\text{with}\quad
           \beta = \frac{f_1 - f_0}{t_1}

       Frequency :math:`f(t)` varies linearly over time with a constant rate
       :math:`\beta`.

    2. Parameter `method` in ``('quadratic', 'quad', 'q')``:

       .. math::
            f(t) =
            \begin{cases}
              f_0 + \beta\, t^2          & \text{if vertex_zero is True,}\\
              f_1 + \beta\, (t_1 - t)^2  & \text{otherwise,}
            \end{cases}
            \quad\text{with}\quad
            \beta = \frac{f_1 - f_0}{t_1^2}

       The graph of the frequency f(t) is a parabola through :math:`(0, f_0)` and
       :math:`(t_1, f_1)`.  By default, the vertex of the parabola is at
       :math:`(0, f_0)`. If `vertex_zero` is ``False``, then the vertex is at
       :math:`(t_1, f_1)`.
       To use a more general quadratic function, or an arbitrary
       polynomial, use the function `scipy.signal.sweep_poly`.

    3. Parameter `method` in ``('logarithmic', 'log', 'lo')``:

       .. math::
            f(t) = f_0  \left(\frac{f_1}{f_0}\right)^{t/t_1}

       :math:`f_0` and :math:`f_1` must be nonzero and have the same sign.
       This signal is also known as a geometric or exponential chirp.

    4. Parameter `method` in ``('hyperbolic', 'hyp')``:

       .. math::
              f(t) = \frac{\alpha}{\beta\, t + \gamma} \quad\text{with}\quad
              \alpha = f_0 f_1 t_1, \ \beta = f_0 - f_1, \ \gamma = f_1 t_1

       :math:`f_0` and :math:`f_1` must be nonzero.


    Examples
    --------
    For the first example, a linear chirp ranging from 6 Hz to 1 Hz over 10 seconds is
    plotted:

    >>> import numpy as np
    >>> from matplotlib.pyplot import tight_layout
    >>> from scipy.signal import chirp, square, ShortTimeFFT
    >>> from scipy.signal.windows import gaussian
    >>> import matplotlib.pyplot as plt
    ...
    >>> N, T = 1000, 0.01  # number of samples and sampling interval for 10 s signal
    >>> t = np.arange(N) * T  # timestamps
    ...
    >>> x_lin = chirp(t, f0=6, f1=1, t1=10, method='linear')
    ...
    >>> fg0, ax0 = plt.subplots()
    >>> ax0.set_title(r"Linear Chirp from $f(0)=6\,$Hz to $f(10)=1\,$Hz")
    >>> ax0.set(xlabel="Time $t$ in Seconds", ylabel=r"Amplitude $x_\text{lin}(t)$")
    >>> ax0.plot(t, x_lin)
    >>> plt.show()

    The following four plots each show the short-time Fourier transform of a chirp
    ranging from 45 Hz to 5 Hz with different values for the parameter `method`
    (and `vertex_zero`):

    >>> x_qu0 = chirp(t, f0=45, f1=5, t1=N*T, method='quadratic', vertex_zero=True)
    >>> x_qu1 = chirp(t, f0=45, f1=5, t1=N*T, method='quadratic', vertex_zero=False)
    >>> x_log = chirp(t, f0=45, f1=5, t1=N*T, method='logarithmic')
    >>> x_hyp = chirp(t, f0=45, f1=5, t1=N*T, method='hyperbolic')
    ...
    >>> win = gaussian(50, std=12, sym=True)
    >>> SFT = ShortTimeFFT(win, hop=2, fs=1/T, mfft=800, scale_to='magnitude')
    >>> ts = ("'quadratic', vertex_zero=True", "'quadratic', vertex_zero=False",
    ...       "'logarithmic'", "'hyperbolic'")
    >>> fg1, ax1s = plt.subplots(2, 2, sharex='all', sharey='all',
    ...                          figsize=(6, 5),  layout="constrained")
    >>> for x_, ax_, t_ in zip([x_qu0, x_qu1, x_log, x_hyp], ax1s.ravel(), ts):
    ...     aSx = abs(SFT.stft(x_))
    ...     im_ = ax_.imshow(aSx, origin='lower', aspect='auto', extent=SFT.extent(N),
    ...                      cmap='plasma')
    ...     ax_.set_title(t_)
    ...     if t_ == "'hyperbolic'":
    ...         fg1.colorbar(im_, ax=ax1s, label='Magnitude $|S_z(t,f)|$')
    >>> _ = fg1.supxlabel("Time $t$ in Seconds")  # `_ =` is needed to pass doctests
    >>> _ = fg1.supylabel("Frequency $f$ in Hertz")
    >>> plt.show()

    Finally, the short-time Fourier transform of a complex-valued linear chirp
    ranging from -30 Hz to 30 Hz is depicted:

    >>> z_lin = chirp(t, f0=-30, f1=30, t1=N*T, method="linear", complex=True)
    >>> SFT.fft_mode = 'centered'  # needed to work with complex signals
    >>> aSz = abs(SFT.stft(z_lin))
    ...
    >>> fg2, ax2 = plt.subplots()
    >>> ax2.set_title(r"Linear Chirp from $-30\,$Hz to $30\,$Hz")
    >>> ax2.set(xlabel="Time $t$ in Seconds", ylabel="Frequency $f$ in Hertz")
    >>> im2 = ax2.imshow(aSz, origin='lower', aspect='auto',
    ...                  extent=SFT.extent(N), cmap='viridis')
    >>> fg2.colorbar(im2, label='Magnitude $|S_z(t,f)|$')
    >>> plt.show()

    Note that using negative frequencies makes only sense with complex-valued signals.
    Furthermore, the magnitude of the complex exponential function is one whereas the
    magnitude of the real-valued cosine function is only 1/2.
    """
    # 'phase' is computed in _chirp_phase, to make testing easier.
    phase = _chirp_phase(t, f0, t1, f1, method, vertex_zero) + np.deg2rad(phi)
    return np.exp(1j*phase) if complex else np.cos(phase)


def _chirp_phase(t, f0, t1, f1, method='linear', vertex_zero=True):
    """
    Calculate the phase used by `chirp` to generate its output.

    See `chirp` for a description of the arguments.

    """
    t = asarray(t)
    f0 = float(f0)
    t1 = float(t1)
    f1 = float(f1)
    if method in ['linear', 'lin', 'li']:
        beta = (f1 - f0) / t1
        phase = 2 * pi * (f0 * t + 0.5 * beta * t * t)

    elif method in ['quadratic', 'quad', 'q']:
        beta = (f1 - f0) / (t1 ** 2)
        if vertex_zero:
            phase = 2 * pi * (f0 * t + beta * t ** 3 / 3)
        else:
            phase = 2 * pi * (f1 * t + beta * ((t1 - t) ** 3 - t1 ** 3) / 3)

    elif method in ['logarithmic', 'log', 'lo']:
        if f0 * f1 <= 0.0:
            raise ValueError("For a logarithmic chirp, f0 and f1 must be "
                             "nonzero and have the same sign.")
        if f0 == f1:
            phase = 2 * pi * f0 * t
        else:
            beta = t1 / log(f1 / f0)
            phase = 2 * pi * beta * f0 * (pow(f1 / f0, t / t1) - 1.0)

    elif method in ['hyperbolic', 'hyp']:
        if f0 == 0 or f1 == 0:
            raise ValueError("For a hyperbolic chirp, f0 and f1 must be "
                             "nonzero.")
        if f0 == f1:
            # Degenerate case: constant frequency.
            phase = 2 * pi * f0 * t
        else:
            # Singular point: the instantaneous frequency blows up
            # when t == sing.
            sing = -f1 * t1 / (f0 - f1)
            phase = 2 * pi * (-sing * f0) * log(np.abs(1 - t/sing))

    else:
        raise ValueError("method must be 'linear', 'quadratic', 'logarithmic', "
                         f"or 'hyperbolic', but a value of {method!r} was given.")

    return phase


def sweep_poly(t, poly, phi=0):
    """
    Frequency-swept cosine generator, with a time-dependent frequency.

    This function generates a sinusoidal function whose instantaneous
    frequency varies with time.  The frequency at time `t` is given by
    the polynomial `poly`.

    Parameters
    ----------
    t : ndarray
        Times at which to evaluate the waveform.
    poly : 1-D array_like or instance of numpy.poly1d
        The desired frequency expressed as a polynomial.  If `poly` is
        a list or ndarray of length n, then the elements of `poly` are
        the coefficients of the polynomial, and the instantaneous
        frequency is

          ``f(t) = poly[0]*t**(n-1) + poly[1]*t**(n-2) + ... + poly[n-1]``

        If `poly` is an instance of numpy.poly1d, then the
        instantaneous frequency is

          ``f(t) = poly(t)``

    phi : float, optional
        Phase offset, in degrees, Default: 0.

    Returns
    -------
    sweep_poly : ndarray
        A numpy array containing the signal evaluated at `t` with the
        requested time-varying frequency.  More precisely, the function
        returns ``cos(phase + (pi/180)*phi)``, where `phase` is the integral
        (from 0 to t) of ``2 * pi * f(t)``; ``f(t)`` is defined above.

    See Also
    --------
    chirp

    Notes
    -----
    .. versionadded:: 0.8.0

    If `poly` is a list or ndarray of length `n`, then the elements of
    `poly` are the coefficients of the polynomial, and the instantaneous
    frequency is:

        ``f(t) = poly[0]*t**(n-1) + poly[1]*t**(n-2) + ... + poly[n-1]``

    If `poly` is an instance of `numpy.poly1d`, then the instantaneous
    frequency is:

          ``f(t) = poly(t)``

    Finally, the output `s` is:

        ``cos(phase + (pi/180)*phi)``

    where `phase` is the integral from 0 to `t` of ``2 * pi * f(t)``,
    ``f(t)`` as defined above.

    Examples
    --------
    Compute the waveform with instantaneous frequency::

        f(t) = 0.025*t**3 - 0.36*t**2 + 1.25*t + 2

    over the interval 0 <= t <= 10.

    >>> import numpy as np
    >>> from scipy.signal import sweep_poly
    >>> p = np.poly1d([0.025, -0.36, 1.25, 2.0])
    >>> t = np.linspace(0, 10, 5001)
    >>> w = sweep_poly(t, p)

    Plot it:

    >>> import matplotlib.pyplot as plt
    >>> plt.subplot(2, 1, 1)
    >>> plt.plot(t, w)
    >>> plt.title("Sweep Poly\\nwith frequency " +
    ...           "$f(t) = 0.025t^3 - 0.36t^2 + 1.25t + 2$")
    >>> plt.subplot(2, 1, 2)
    >>> plt.plot(t, p(t), 'r', label='f(t)')
    >>> plt.legend()
    >>> plt.xlabel('t')
    >>> plt.tight_layout()
    >>> plt.show()

    """
    # 'phase' is computed in _sweep_poly_phase, to make testing easier.
    phase = _sweep_poly_phase(t, poly)
    # Convert to radians.
    phi *= pi / 180
    return cos(phase + phi)


def _sweep_poly_phase(t, poly):
    """
    Calculate the phase used by sweep_poly to generate its output.

    See `sweep_poly` for a description of the arguments.

    """
    # polyint handles lists, ndarrays and instances of poly1d automatically.
    intpoly = polyint(poly)
    phase = 2 * pi * polyval(intpoly, t)
    return phase


def unit_impulse(shape, idx=None, dtype=float):
    r"""
    Unit impulse signal (discrete delta function) or unit basis vector.

    Parameters
    ----------
    shape : int or tuple of int
        Number of samples in the output (1-D), or a tuple that represents the
        shape of the output (N-D).
    idx : None or int or tuple of int or 'mid', optional
        Index at which the value is 1.  If None, defaults to the 0th element.
        If ``idx='mid'``, the impulse will be centered at ``shape // 2`` in
        all dimensions.  If an int, the impulse will be at `idx` in all
        dimensions.
    dtype : data-type, optional
        The desired data-type for the array, e.g., ``numpy.int8``.  Default is
        ``numpy.float64``.

    Returns
    -------
    y : ndarray
        Output array containing an impulse signal.

    Notes
    -----
    In digital signal processing literature the unit impulse signal is often
    represented by the Kronecker delta. [1]_ I.e., a signal :math:`u_k[n]`,
    which is zero everywhere except being one at the :math:`k`-th sample,
    can be expressed as

    .. math::

        u_k[n] = \delta[n-k] \equiv \delta_{n,k}\ .

    Furthermore, the unit impulse is frequently interpreted as the discrete-time
    version of the continuous-time Dirac distribution. [2]_

    References
    ----------
    .. [1] "Kronecker delta", *Wikipedia*,
           https://en.wikipedia.org/wiki/Kronecker_delta#Digital_signal_processing
    .. [2] "Dirac delta function" *Wikipedia*,
           https://en.wikipedia.org/wiki/Dirac_delta_function#Relationship_to_the_Kronecker_delta

    .. versionadded:: 0.19.0

    Examples
    --------
    An impulse at the 0th element (:math:`\\delta[n]`):

    >>> from scipy import signal
    >>> signal.unit_impulse(8)
    array([ 1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

    Impulse offset by 2 samples (:math:`\\delta[n-2]`):

    >>> signal.unit_impulse(7, 2)
    array([ 0.,  0.,  1.,  0.,  0.,  0.,  0.])

    2-dimensional impulse, centered:

    >>> signal.unit_impulse((3, 3), 'mid')
    array([[ 0.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  0.]])

    Impulse at (2, 2), using broadcasting:

    >>> signal.unit_impulse((4, 4), 2)
    array([[ 0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.],
           [ 0.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  0.]])

    Plot the impulse response of a 4th-order Butterworth lowpass filter:

    >>> imp = signal.unit_impulse(100, 'mid')
    >>> b, a = signal.butter(4, 0.2)
    >>> response = signal.lfilter(b, a, imp)

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(np.arange(-50, 50), imp)
    >>> plt.plot(np.arange(-50, 50), response)
    >>> plt.margins(0.1, 0.1)
    >>> plt.xlabel('Time [samples]')
    >>> plt.ylabel('Amplitude')
    >>> plt.grid(True)
    >>> plt.show()

    """
    out = zeros(shape, dtype)

    shape = np.atleast_1d(shape)

    if idx is None:
        idx = (0,) * len(shape)
    elif idx == 'mid':
        idx = tuple(shape // 2)
    elif not hasattr(idx, "__iter__"):
        idx = (idx,) * len(shape)

    out[idx] = 1
    return out
