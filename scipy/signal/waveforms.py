# Author: Travis Oliphant
# 2003
#
# Feb. 2010: Updated by Warren Weckesser:
#   Rewrote much of chirp()
#   Added sweep_poly()

from numpy import asarray, zeros, place, nan, mod, pi, extract, log, sqrt, \
     exp, cos, sin, polyval, polyint

def sawtooth(t,width=1):
    """Returns a periodic sawtooth waveform with period 2*pi
    which rises from -1 to 1 on the interval 0 to width*2*pi
    and drops from 1 to -1 on the interval width*2*pi to 2*pi
    width must be in the interval [0,1]

    """
    t,w = asarray(t), asarray(width)
    w = asarray(w + (t-t))
    t = asarray(t + (w-w))
    if t.dtype.char in ['fFdD']:
        ytype = t.dtype.char
    else:
        ytype = 'd'
    y = zeros(t.shape,ytype)

    # width must be between 0 and 1 inclusive
    mask1 = (w > 1) | (w < 0)
    place(y,mask1,nan)

    # take t modulo 2*pi
    tmod = mod(t,2*pi)

    # on the interval 0 to width*2*pi function is
    #  tmod / (pi*w) - 1
    mask2 = (1-mask1) & (tmod < w*2*pi)
    tsub = extract(mask2,tmod)
    wsub = extract(mask2,w)
    place(y,mask2,tsub / (pi*wsub) - 1)

    # on the interval width*2*pi to 2*pi function is
    #  (pi*(w+1)-tmod) / (pi*(1-w))

    mask3 = (1-mask1) & (1-mask2)
    tsub = extract(mask3,tmod)
    wsub = extract(mask3,w)
    place(y,mask3, (pi*(wsub+1)-tsub)/(pi*(1-wsub)))
    return y


def square(t,duty=0.5):
    """Returns a periodic square-wave waveform with period 2*pi
    which is +1 from 0 to 2*pi*duty and -1 from 2*pi*duty to 2*pi
    duty must be in the interval [0,1]

    """
    t,w = asarray(t), asarray(duty)
    w = asarray(w + (t-t))
    t = asarray(t + (w-w))
    if t.dtype.char in ['fFdD']:
        ytype = t.dtype.char
    else:
        ytype = 'd'
    y = zeros(t.shape,ytype)

    # width must be between 0 and 1 inclusive
    mask1 = (w > 1) | (w < 0)
    place(y,mask1,nan)

    # take t modulo 2*pi
    tmod = mod(t,2*pi)

    # on the interval 0 to duty*2*pi function is
    #  1
    mask2 = (1-mask1) & (tmod < w*2*pi)
    tsub = extract(mask2,tmod)
    wsub = extract(mask2,w)
    place(y,mask2,1)

    # on the interval duty*2*pi to 2*pi function is
    #  (pi*(w+1)-tmod) / (pi*(1-w))

    mask3 = (1-mask1) & (1-mask2)
    tsub = extract(mask3,tmod)
    wsub = extract(mask3,w)
    place(y,mask3,-1)
    return y

def gausspulse(t,fc=1000,bw=0.5,bwr=-6,tpr=-60,retquad=0,retenv=0):
    """Return a gaussian modulated sinusoid:  exp(-a t^2) exp(1j*2*pi*fc)

    If retquad is non-zero, then return the real and imaginary parts
       (inphase and quadrature)
    If retenv is non-zero, then return the envelope (unmodulated signal).
    Otherwise, return the real part of the modulated sinusoid.

    Inputs:

       t   --  Input array.
       fc  --  Center frequency (Hz).
       bw  --  Fractional bandwidth in frequency domain of pulse (Hz).
       bwr --  Reference level at which fractional bandwidth is calculated (dB).
       tpr --  If t is 'cutoff', then the function returns the cutoff time for when the
                  pulse amplitude falls below tpr (in dB).
       retquad -- Return the quadrature (imaginary) as well as the real part of the signal
       retenv  -- Return the envelope of th signal.

    """
    if fc < 0:
        raise ValueError, "Center frequency (fc=%.2f) must be >=0." % fc
    if bw <= 0:
        raise ValueError, "Fractional bandwidth (bw=%.2f) must be > 0." % bw
    if bwr >= 0:
        raise ValueError, "Reference level for bandwidth (bwr=%.2f) must " \
              "be < 0 dB" % bwr

    # exp(-a t^2) <->  sqrt(pi/a) exp(-pi^2/a * f^2)  = g(f)

    ref = pow(10, bwr/ 20)
    # fdel = fc*bw/2:  g(fdel) = ref --- solve this for a
    #
    # pi^2/a * fc^2 * bw^2 /4=-log(ref)
    a = -(pi*fc*bw)**2 / (4*log(ref))

    if t == 'cutoff': # compute cut_off point
        #  Solve exp(-a tc**2) = tref  for tc
        #   tc = sqrt(-log(tref) / a) where tref = 10^(tpr/20)
        if tpr >= 0:
            raise ValueError, "Reference level for time cutoff must be < 0 dB"
        tref = pow(10, tpr / 20)
        return sqrt(-log(tref)/a)

    yenv = exp(-a*t*t)
    yI = yenv * cos(2*pi*fc*t)
    yQ = yenv * sin(2*pi*fc*t)
    if not retquad and not retenv:
        return yI
    if not retquad and retenv:
        return yI, yenv
    if retquad and not retenv:
        return yI, yQ
    if retquad and retenv:
        return yI, yQ, yenv


def chirp(t, f0, t1, f1, method='linear', phi=0, vertex_zero=True):
    """Frequency-swept cosine generator.

    In the following, 'Hz' should be interpreted as 'cycles per time unit'; there is
    no assumption here that the time unit is one second.  The important distinction
    is that the units of rotation are cycles, not radians.

    Parameters
    ----------
    t : ndarray
        Times at which to evaluate the waveform.
    f0 : float
        Frequency (in Hz) at time t=0.
    t1 : float
        Time at which `f1` is specified.
    f1 : float
        Frequency (in Hz) of the waveform at time `t1`.
    method : {'linear', 'quadratic', 'logarithmic', 'hyperbolic'}, optional
        Kind of frequency sweep.  If not given, `linear` is assumed.  See Notes below
        for more details.
    phi : float, optional
        Phase offset, in degrees. Default is 0.
    vertex_zero : bool, optional
        This parameter is only used when `method` is 'quadratic'.
        It determines whether the vertex of the parabola that is the graph of the
        frequency is at t=0 or t=t1.

    Returns
    -------
    A numpy array containing the signal evaluated at 't' with the requested
    time-varying frequency.  More precisely, the function returns
        ``cos(phase + (pi/180)*phi)``
    where `phase` is the integral (from 0 to t) of ``2*pi*f(t)``.
    ``f(t)`` is defined below.

    Notes
    -----
    There are four options for the `method`.  The following formulas give the
    instantaneous frequency (in Hz) of the signal generated by `chirp()`.
    For convenience, the shorter names shown below may also be used.

    linear, lin, li:

        ``f(t) = f0 + (f1 - f0) * t / t1``

    quadratic, quad, q:

        The graph of the frequency f(t) is a parabola through (0, f0) and (t1, f1).
        By default, the vertex of the parabola is at (0, f0).  If `vertex_zero`
        is False, then the vertex is at (t1, f1).  The formula is:

        if vertex_zero is True:

            ``f(t) = f0 + (f1 - f0) * t**2 / t1**2``

        else

            ``f(t) = f1 - (f1 - f0) * (t1 - t)**2 / t1**2``

        To use a more general quadratic function, or an arbitrary polynomial,
        use the function `scipy.signal.waveforms.sweep_poly`.

    logarithmic, log, lo:

        ``f(t) = f0 * (f1/f0)**(t/t1)``

        f0 and f1 must be nonzero and have the same sign.

        This signal is also known as a geometric or exponential chirp.

    hyperbolic, hyp:

        ``f(t) = f0*f1*t1 / ((f0 - f1)*t + f1*t1)``

        f1 must be positive, and f0 must be greater than f1.

    See Also
    --------
    scipy.signal.waveforms.sweep_poly

    """
    # 'phase' is computed in _chirp_phase, to make testing easier.
    phase = _chirp_phase(t, f0, t1, f1, method, vertex_zero)
    # Convert  phi to radians.
    phi *= pi / 180
    return cos(phase + phi)


def _chirp_phase(t, f0, t1, f1, method='linear', vertex_zero=True):
    """
    Calculate the phase used by chirp_phase to generate its output.  See
    chirp_phase for a description of the arguments.

    """
    if method in ['linear', 'lin', 'li']:
        beta = (f1 - f0) / t1
        phase = 2*pi * (f0*t + 0.5*beta*t*t)

    elif method in ['quadratic','quad','q']:
        beta = (f1 - f0)/(t1**2)
        if vertex_zero:
            phase = 2*pi * (f0*t + beta * t**3/3)
        else:
            phase = 2*pi * (f1*t + beta * ((t1 - t)**3 - t1**3)/3)

    elif method in ['logarithmic', 'log', 'lo']:
        if f0*f1 <= 0.0:
            raise ValueError("For a geometric chirp, f0 and f1 must be nonzero " \
                                "and have the same sign.")
        if f0 == f1:
            phase = 2*pi * f0 * t
        else:
            beta = t1 / log(f1/f0)
            phase = 2*pi * beta * f0 * (pow(f1/f0, t/t1) - 1.0)

    elif method in ['hyperbolic', 'hyp']:
        if f1 <= 0.0 or f0 <= f1:
            raise ValueError("hyperbolic chirp requires f0 > f1 > 0.0.")
        c = f1*t1
        df = f0 - f1
        phase = 2*pi * (f0 * c / df) * log((df*t + c)/c)

    else:
        raise ValueError("method must be 'linear', 'quadratic', 'logarithmic', "
                "or 'hyperbolic', but a value of %r was given." % method)

    return phase


def sweep_poly(t, poly, phi=0):
    """Frequency-swept cosine generator, with a time-dependent frequency specified
    as a polynomial.

    This function generates a sinusoidal function whose instantaneous frequency
    varies with time.  The frequency at time `t` is given by the polynomial `poly`.

    Parameters
    ----------
    t : ndarray
        Times at which to evaluate the waveform.
    poly : 1D ndarray (or array-like), or instance of numpy.poly1d
        The desired frequency expressed as a polynomial.  If `poly` is
        a list or ndarray of length n, then the elements of `poly` are
        the coefficients of the polynomial, and the instantaneous
        frequency is
            ``f(t) = poly[0]*t**(n-1) + poly[1]*t**(n-2) + ... + poly[n-1]``
        If `poly` is an instance of numpy.poly1d, then the
        instantaneous frequency is
            ``f(t) = poly(t)``
    phi : float, optional
        Phase offset, in degrees. Default is 0.

    Returns
    -------
    A numpy array containing the signal evaluated at 't' with the requested
    time-varying frequency.  More precisely, the function returns
        ``cos(phase + (pi/180)*phi)``
    where `phase` is the integral (from 0 to t) of ``2 * pi * f(t)``;
    ``f(t)`` is defined above.

    See Also
    --------
    scipy.signal.waveforms.chirp

    """
    # 'phase' is computed in _sweep_poly_phase, to make testing easier.
    phase = _sweep_poly_phase(t, poly)
    # Convert to radians.
    phi *= pi / 180
    return cos(phase + phi)

def _sweep_poly_phase(t, poly):
    """
    Calculate the phase used by sweep_poly to generate its output.  See
    sweep_poly for a description of the arguments.

    """
    # polyint handles lists, ndarrays and instances of poly1d automatically.
    intpoly = polyint(poly)
    phase = 2*pi * polyval(intpoly, t)
    return phase
