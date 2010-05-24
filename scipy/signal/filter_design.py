"""Filter design.
"""

import types
import warnings

import numpy
from numpy import atleast_1d, poly, polyval, roots, real, asarray, allclose, \
    resize, pi, absolute, logspace, r_, sqrt, tan, log10, arctan, arcsinh, \
    cos, exp, cosh, arccosh, ceil, conjugate, zeros, sinh
from numpy import mintypecode
from scipy import special, optimize
from scipy.misc import comb

class BadCoefficients(UserWarning):
    pass

abs = absolute

def findfreqs(num, den, N):
    ep = atleast_1d(roots(den))+0j
    tz = atleast_1d(roots(num))+0j

    if len(ep) == 0:
        ep = atleast_1d(-1000)+0j

    ez = r_['-1',numpy.compress(ep.imag >=0, ep,axis=-1), numpy.compress((abs(tz) < 1e5) & (tz.imag >=0),tz,axis=-1)]

    integ = abs(ez) < 1e-10
    hfreq = numpy.around(numpy.log10(numpy.max(3*abs(ez.real + integ)+1.5*ez.imag))+0.5)
    lfreq = numpy.around(numpy.log10(0.1*numpy.min(abs(real(ez+integ))+2*ez.imag))-0.5)

    w = logspace(lfreq, hfreq, N)
    return w

def freqs(b,a,worN=None,plot=None):
    """Compute frequency response of analog filter.

    Given the numerator (b) and denominator (a) of a filter compute its
    frequency response.

            b[0]*(jw)**(nb-1) + b[1]*(jw)**(nb-2) + ... + b[nb-1]
    H(w) = --------------------------------------------------------
            a[0]*(jw)**(na-1) + a[1]*(jw)**(na-2) + ... + a[na-1]

    Parameters
    ----------
    b : ndarray
        numerator of a linear filter
    a : ndarray
        numerator of a linear filter
    worN : {None, int}, optional
        If None, then compute at 200 frequencies around the interesting parts
        of the response curve (determined by pole-zero locations).  If a single
        integer, the compute at that many frequencies.  Otherwise, compute the
        response at frequencies given in worN.

    Returns
    -------
    w : ndarray
        The frequencies at which h was computed.
    h : ndarray
        The frequency response.
    """
    if worN is None:
        w = findfreqs(b,a,200)
    elif isinstance(worN, types.IntType):
        N = worN
        w = findfreqs(b,a,N)
    else:
        w = worN
    w = atleast_1d(w)
    s = 1j*w
    h = polyval(b, s) / polyval(a, s)
    if not plot is None:
        plot(w, h)
    return w, h

def freqz(b, a=1, worN=None, whole=0, plot=None):
    """
    Compute the frequency response of a digital filter.

    Given the numerator ``b`` and denominator ``a`` of a digital filter compute
    its frequency response::

               jw               -jw            -jmw
        jw  B(e)    b[0] + b[1]e + .... + b[m]e
     H(e) = ---- = ------------------------------------
               jw               -jw            -jnw
            A(e)    a[0] + a[1]e + .... + a[n]e

    Parameters
    ----------
    b : ndarray
        numerator of a linear filter
    a : ndarray
        numerator of a linear filter
    worN : {None, int}, optional
        If None, then compute at 512 frequencies around the unit circle.
        If a single integer, the compute at that many frequencies.
        Otherwise, compute the response at frequencies given in worN
    whole : {0,1}, optional
        Normally, frequencies are computed from 0 to pi (upper-half of
        unit-circle.  If whole is non-zero compute frequencies from 0 to 2*pi.

    Returns
    -------
    w : ndarray
        The frequencies at which h was computed.
    h : ndarray
        The frequency response.

    Examples
    --------

    >>> b = firwin(80, 0.5, window=('kaiser', 8))
    >>> h, w = freqz(b)

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> plt.title('Digital filter frequency response')
    >>> ax1 = fig.add_subplot(111)

    >>> plt.semilogy(h, np.abs(w), 'b')
    >>> plt.ylabel('Amplitude (dB)', color='b')
    >>> plt.xlabel('Frequency (rad/sample)')
    >>> plt.grid()
    >>> plt.legend()

    >>> ax2 = ax1.twinx()
    >>> angles = np.unwrap(np.angle(w))
    >>> plt.plot(h, angles, 'g')
    >>> plt.ylabel('Angle (radians)', color='g')
    >>> plt.show()

    """
    b, a = map(atleast_1d, (b,a))
    if whole:
        lastpoint = 2*pi
    else:
        lastpoint = pi
    if worN is None:
        N = 512
        w = numpy.arange(0,lastpoint,lastpoint/N)
    elif isinstance(worN, types.IntType):
        N = worN
        w = numpy.arange(0,lastpoint,lastpoint/N)
    else:
        w = worN
    w = atleast_1d(w)
    zm1 = exp(-1j*w)
    h = polyval(b[::-1], zm1) / polyval(a[::-1], zm1)
    if not plot is None:
        plot(w, h)
    return w, h

def tf2zpk(b,a):
    """Return zero, pole, gain (z,p,k) representation from a numerator,
    denominator representation of a linear filter.

    Parameters
    ----------
    b : ndarray
        numerator polynomial.
    a : ndarray
        numerator and denominator polynomials.

    Returns
    -------
    z : ndarray
        zeros of the transfer function.
    p : ndarray
        poles of the transfer function.
    k : float
        system gain.

    If some values of b are too close to 0, they are removed. In that case, a
    BadCoefficients warning is emitted.
    """
    b,a = normalize(b,a)
    b = (b+0.0) / a[0]
    a = (a+0.0) / a[0]
    k = b[0]
    b /= b[0]
    z = roots(b)
    p = roots(a)
    return z, p, k

def zpk2tf(z,p,k):
    """Return polynomial transfer function representation from zeros
    and poles

    Parameters
    ----------
    z : ndarray
        zeros of the transfer function.
    p : ndarray
        poles of the transfer function.
    k : float
        system gain.

    Returns
    -------
    b : ndarray
        numerator polynomial.
    a : ndarray
        numerator and denominator polynomials.

    Note
    ----
    If some values of b are too close to 0, they are removed. In that case, a
    BadCoefficients warning is emitted.
    """
    z = atleast_1d(z)
    k = atleast_1d(k)
    if len(z.shape) > 1:
        temp = poly(z[0])
        b = zeros((z.shape[0], z.shape[1]+1), temp.dtype.char)
        if len(k) == 1:
            k = [k[0]]*z.shape[0]
        for i in range(z.shape[0]):
            b[i] = k[i] * poly(z[i])
    else:
        b = k * poly(z)
    a = poly(p)
    return b, a

def normalize(b,a):
    """Normalize polynomial representation of a transfer function.

    If values of b are too close to 0, they are removed. In that case, a
    BadCoefficients warning is emitted.
    """
    b,a = map(atleast_1d,(b,a))
    if len(a.shape) != 1:
        raise ValueError, "Denominator polynomial must be rank-1 array."
    if len(b.shape) > 2:
        raise ValueError, "Numerator polynomial must be rank-1 or rank-2 array."
    if len(b.shape) == 1:
        b = asarray([b],b.dtype.char)
    while a[0] == 0.0 and len(a) > 1:
        a = a[1:]
    outb = b * (1.0) / a[0]
    outa = a * (1.0) / a[0]
    if allclose(outb[:,0], 0, rtol=1e-14):
        warnings.warn("Badly conditioned filter coefficients (numerator): the "
                      "results may be meaningless", BadCoefficients)
        while allclose(outb[:,0], 0, rtol=1e-14) and (outb.shape[-1] > 1):
            outb = outb[:,1:]
    if outb.shape[0] == 1:
        outb = outb[0]
    return outb, outa


def lp2lp(b,a,wo=1.0):
    """Return a low-pass filter with cuttoff frequency wo
    from a low-pass filter prototype with unity cutoff frequency.
    """
    a,b = map(atleast_1d,(a,b))
    try:
        wo = float(wo)
    except TypeError:
        wo = float(wo[0])
    d = len(a)
    n = len(b)
    M = max((d,n))
    pwo = pow(wo,numpy.arange(M-1,-1,-1))
    start1 = max((n-d,0))
    start2 = max((d-n,0))
    b = b * pwo[start1]/pwo[start2:]
    a = a * pwo[start1]/pwo[start1:]
    return normalize(b, a)

def lp2hp(b,a,wo=1.0):
    """Return a high-pass filter with cuttoff frequency wo
    from a low-pass filter prototype with unity cutoff frequency.
    """
    a,b = map(atleast_1d,(a,b))
    try:
        wo = float(wo)
    except TypeError:
        wo = float(wo[0])
    d = len(a)
    n = len(b)
    if wo != 1:
        pwo = pow(wo,numpy.arange(max((d,n))))
    else:
        pwo = numpy.ones(max((d,n)),b.dtype.char)
    if d >= n:
        outa = a[::-1] * pwo
        outb = resize(b,(d,))
        outb[n:] = 0.0
        outb[:n] = b[::-1] * pwo[:n]
    else:
        outb = b[::-1] * pwo
        outa = resize(a,(n,))
        outa[d:] = 0.0
        outa[:d] = a[::-1] * pwo[:d]

    return normalize(outb, outa)

def lp2bp(b,a,wo=1.0, bw=1.0):
    """Return a band-pass filter with center frequency wo and bandwidth bw
    from a low-pass filter prototype with unity cutoff frequency.
    """
    a,b = map(atleast_1d,(a,b))
    D = len(a) - 1
    N = len(b) - 1
    artype = mintypecode((a,b))
    ma = max([N,D])
    Np = N + ma
    Dp = D + ma
    bprime = numpy.zeros(Np+1,artype)
    aprime = numpy.zeros(Dp+1,artype)
    wosq = wo*wo
    for j in range(Np+1):
        val = 0.0
        for i in range(0,N+1):
            for k in range(0,i+1):
                if ma-i+2*k == j:
                    val += comb(i,k)*b[N-i]*(wosq)**(i-k) / bw**i
        bprime[Np-j] = val
    for j in range(Dp+1):
        val = 0.0
        for i in range(0,D+1):
            for k in range(0,i+1):
                if ma-i+2*k == j:
                    val += comb(i,k)*a[D-i]*(wosq)**(i-k) / bw**i
        aprime[Dp-j] = val

    return normalize(bprime, aprime)

def lp2bs(b,a,wo=1,bw=1):
    """Return a band-stop filter with center frequency wo and bandwidth bw
    from a low-pass filter prototype with unity cutoff frequency.
    """
    a,b = map(atleast_1d,(a,b))
    D = len(a) - 1
    N = len(b) - 1
    artype = mintypecode((a,b))
    M = max([N,D])
    Np = M + M
    Dp = M + M
    bprime = numpy.zeros(Np+1,artype)
    aprime = numpy.zeros(Dp+1,artype)
    wosq = wo*wo
    for j in range(Np+1):
        val = 0.0
        for i in range(0,N+1):
            for k in range(0,M-i+1):
                if i+2*k == j:
                    val += comb(M-i,k)*b[N-i]*(wosq)**(M-i-k) * bw**i
        bprime[Np-j] = val
    for j in range(Dp+1):
        val = 0.0
        for i in range(0,D+1):
            for k in range(0,M-i+1):
                if i+2*k == j:
                    val += comb(M-i,k)*a[D-i]*(wosq)**(M-i-k) * bw**i
        aprime[Dp-j] = val

    return normalize(bprime, aprime)

def bilinear(b,a,fs=1.0):
    """Return a digital filter from an analog filter using the bilinear transform.

    The bilinear transform substitutes (z-1) / (z+1) for s
    """
    fs =float(fs)
    a,b = map(atleast_1d,(a,b))
    D = len(a) - 1
    N = len(b) - 1
    artype = float
    M = max([N,D])
    Np = M
    Dp = M
    bprime = numpy.zeros(Np+1,artype)
    aprime = numpy.zeros(Dp+1,artype)
    for j in range(Np+1):
        val = 0.0
        for i in range(N+1):
            for k in range(i+1):
                for l in range(M-i+1):
                    if k+l == j:
                        val += comb(i,k)*comb(M-i,l)*b[N-i]*pow(2*fs,i)*(-1)**k
        bprime[j] = real(val)
    for j in range(Dp+1):
        val = 0.0
        for i in range(D+1):
            for k in range(i+1):
                for l in range(M-i+1):
                    if k+l == j:
                        val += comb(i,k)*comb(M-i,l)*a[D-i]*pow(2*fs,i)*(-1)**k
        aprime[j] = real(val)

    return normalize(bprime, aprime)

def iirdesign(wp, ws, gpass, gstop, analog=0, ftype='ellip', output='ba'):
    """Complete IIR digital and analog filter design.

    Given passband and stopband frequencies and gains construct an analog or
    digital IIR filter of minimum order for a given basic type.  Return the
    output in numerator, denominator ('ba') or pole-zero ('zpk') form.

    Parameters
    ----------
    wp, ws -- Passband and stopband edge frequencies, normalized from 0
              to 1 (1 corresponds to pi radians / sample).  For example:
                 Lowpass:   wp = 0.2,          ws = 0.3
                 Highpass:  wp = 0.3,          ws = 0.2
                 Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
                 Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]
    gpass -- The maximum loss in the passband (dB).
    gstop -- The minimum attenuation in the stopband (dB).
    analog -- Non-zero to design an analog filter (in this case wp and
              ws are in radians / second).
    ftype -- The type of iir filter to design:
               elliptic    : 'ellip'
               Butterworth : 'butter',
               Chebyshev I : 'cheby1',
               Chebyshev II: 'cheby2',
               Bessel :      'bessel'
    output -- Type of output:  numerator/denominator ('ba') or pole-zero ('zpk')

    Returns
    -------
      b,a -- Numerator and denominator of the iir filter.
      z,p,k -- Zeros, poles, and gain of the iir filter.
    """

    try:
        ordfunc = filter_dict[ftype][1]
    except KeyError:
        raise ValueError, "Invalid IIR filter type."
    except IndexError:
        raise ValueError, "%s does not have order selection use iirfilter function." % ftype

    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    band_type = 2*(len(wp)-1)
    band_type +=1
    if wp[0] >= ws[0]:
        band_type += 1

    btype = {1:'lowpass', 2:'highpass', 3:'bandstop', 4:'bandpass'}[band_type]

    N, Wn = ordfunc(wp, ws, gpass, gstop, analog=analog)
    return iirfilter(N, Wn, rp=gpass, rs=gstop, analog=analog, btype=btype, ftype=ftype, output=output)


def iirfilter(N, Wn, rp=None, rs=None, btype='band', analog=0, ftype='butter', output='ba'):
    """IIR digital and analog filter design given order and critical points.

    Design an Nth order lowpass digital or analog filter and return the filter
    coefficients in (B,A) (numerator, denominator) or (Z,P,K) form.

    Parameters
    ----------
    N -- the order of the filter.
    Wn -- a scalar or length-2 sequence giving the critical frequencies.
    rp, rs -- For chebyshev and elliptic filters provides the maximum ripple
              in the passband and the minimum attenuation in the stop band.
    btype -- the type of filter (lowpass, highpass, bandpass, or bandstop).
    analog -- non-zero to return an analog filter, otherwise
              a digital filter is returned.
    ftype -- the type of IIR filter (Butterworth, Cauer (Elliptic),
             Bessel, Chebyshev1, Chebyshev2)
    output -- 'ba' for (b,a) output, 'zpk' for (z,p,k) output.

    SEE ALSO butterord, cheb1ord, cheb2ord, ellipord
    """

    ftype, btype, output = [x.lower() for x in (ftype, btype, output)]
    Wn = asarray(Wn)
    try:
        btype = band_dict[btype]
    except KeyError:
        raise ValueError, "%s is an invalid bandtype for filter." % btype

    try:
        typefunc = filter_dict[ftype][0]
    except KeyError:
        raise ValueError, "%s is not a valid basic iir filter." % ftype

    if output not in ['ba', 'zpk']:
        raise ValueError, "%s is not a valid output form." % output

    #pre-warp frequencies for digital filter design
    if not analog:
        fs = 2.0
        warped = 2*fs*tan(pi*Wn/fs)
    else:
        warped = Wn

    # convert to low-pass prototype
    if btype in ['lowpass', 'highpass']:
        wo = warped
    else:
        bw = warped[1] - warped[0]
        wo = sqrt(warped[0]*warped[1])

    # Get analog lowpass prototype
    if typefunc in [buttap, besselap]:
        z, p, k = typefunc(N)
    elif typefunc == cheb1ap:
        if rp is None:
            raise ValueError, "passband ripple (rp) must be provided to design a Chebyshev I filter."
        z, p, k = typefunc(N, rp)
    elif typefunc == cheb2ap:
        if rs is None:
            raise ValueError, "stopband atteunatuion (rs) must be provided to design an Chebyshev II filter."
        z, p, k = typefunc(N, rs)
    else:  # Elliptic filters
        if rs is None or rp is None:
            raise ValueError, "Both rp and rs must be provided to design an elliptic filter."
        z, p, k = typefunc(N, rp, rs)

    b, a = zpk2tf(z,p,k)

    # transform to lowpass, bandpass, highpass, or bandstop
    if btype == 'lowpass':
        b, a = lp2lp(b,a,wo=wo)
    elif btype == 'highpass':
        b, a = lp2hp(b,a,wo=wo)
    elif btype == 'bandpass':
        b, a = lp2bp(b,a,wo=wo,bw=bw)
    else: # 'bandstop'
        b, a = lp2bs(b,a,wo=wo,bw=bw)


    # Find discrete equivalent if necessary
    if not analog:
        b, a = bilinear(b, a, fs=fs)

    # Transform to proper out type (pole-zero, state-space, numer-denom)
    if output == 'zpk':
        return tf2zpk(b,a)
    else:
        return b,a


def butter(N, Wn, btype='low', analog=0, output='ba'):
    """Butterworth digital and analog filter design.

    Description:

      Design an Nth order lowpass digital or analog Butterworth filter
      and return the filter coefficients in (B,A) or (Z,P,K) form.

    See also buttord.
    """
    return iirfilter(N, Wn, btype=btype, analog=analog, output=output, ftype='butter')

def cheby1(N, rp, Wn, btype='low', analog=0, output='ba'):
    """Chebyshev type I digital and analog filter design.

    Description:

      Design an Nth order lowpass digital or analog Chebyshev type I filter
      and return the filter coefficients in (B,A) or (Z,P,K) form.

    See also cheb1ord.
    """
    return iirfilter(N, Wn, rp=rp, btype=btype, analog=analog, output=output, ftype='cheby1')

def cheby2(N, rs, Wn, btype='low', analog=0, output='ba'):
    """Chebyshev type I digital and analog filter design.

    Description:

      Design an Nth order lowpass digital or analog Chebyshev type I filter
      and return the filter coefficients in (B,A) or (Z,P,K) form.

    See also cheb2ord.
    """
    return iirfilter(N, Wn, rs=rs, btype=btype, analog=analog, output=output, ftype='cheby2')

def ellip(N, rp, rs, Wn, btype='low', analog=0, output='ba'):
    """Elliptic (Cauer) digital and analog filter design.

    Description:

      Design an Nth order lowpass digital or analog elliptic filter
      and return the filter coefficients in (B,A) or (Z,P,K) form.

    See also ellipord.
    """
    return iirfilter(N, Wn, rs=rs, rp=rp, btype=btype, analog=analog, output=output, ftype='elliptic')

def bessel(N, Wn, btype='low', analog=0, output='ba'):
    """Bessel digital and analog filter design.

    Description:

      Design an Nth order lowpass digital or analog Bessel filter
      and return the filter coefficients in (B,A) or (Z,P,K) form.

    """
    return iirfilter(N, Wn, btype=btype, analog=analog, output=output, ftype='bessel')


def maxflat():
    pass

def yulewalk():
    pass


def band_stop_obj(wp, ind, passb, stopb, gpass, gstop, type):
    """Band Stop Objective Function for order minimization

    Description:

      Returns the non-integer order for an analog band stop filter.

    Parameters
    ----------
    wp -- passb edge
    ind -- index specifying which passb edge to vary (0 or 1).
    passb -- two element vector of fixed passband edges.
    stopb -- two element vector of fixed stopband edges.
    gstop -- amount in dB of attenuation in stopband.
    gpass -- amount in dB of ripple in the passband.
    type -- 'butter', 'cheby', or 'ellip':

    Returns
    -------
    n -- filter order (possibly non-integer)
    """

    passbC = passb.copy()
    passbC[ind] = wp
    nat = stopb*(passbC[0]-passbC[1]) / (stopb**2 - passbC[0]*passbC[1])
    nat = min(abs(nat))

    if type == 'butter':
        GSTOP = 10**(0.1*abs(gstop))
        GPASS = 10**(0.1*abs(gpass))
        n = (log10((GSTOP-1.0)/(GPASS-1.0)) / (2*log10(nat)))
    elif type == 'cheby':
        GSTOP = 10**(0.1*abs(gstop))
        GPASS = 10**(0.1*abs(gpass))
        n = arccosh(sqrt((GSTOP-1.0)/(GPASS-1.0))) / arccosh(nat)
    elif type == 'ellip':
        GSTOP = 10**(0.1*gstop)
        GPASS = 10**(0.1*gpass)
        arg1 = sqrt( (GPASS-1.0) / (GSTOP-1.0) )
        arg0 = 1.0 / nat
        d0 = special.ellipk([arg0**2, 1-arg0**2])
        d1 = special.ellipk([arg1**2, 1-arg1**2])
        n = (d0[0]*d1[1] / (d0[1]*d1[0]))
    else:
        raise ValueError, "Incorrect type: ", type
    return n

def buttord(wp, ws, gpass, gstop, analog=0):
    """Butterworth filter order selection.

    Return the order of the lowest order digital Butterworth filter that loses
    no more than gpass dB in the passband and has at least gstop dB attenuation
    in the stopband.

    Parameters
    ----------
    wp, ws -- Passband and stopband edge frequencies, normalized from 0
              to 1 (1 corresponds to pi radians / sample).  For example:
                 Lowpass:   wp = 0.2,          ws = 0.3
                 Highpass:  wp = 0.3,          ws = 0.2
                 Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
                 Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]
    gpass -- The maximum loss in the passband (dB).
    gstop -- The minimum attenuation in the stopband (dB).
    analog -- Non-zero to design an analog filter (in this case wp and
              ws are in radians / second).

    Returns
    -------
    ord -- The lowest order for a Butterworth filter which meets specs.
    Wn -- The Butterworth natural frequency (i.e. the "3dB frequency").
          Should be used with scipy.signal.butter to give filter results.

    """

    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2*(len(wp)-1)
    filter_type +=1
    if wp[0] >= ws[0]:
        filter_type += 1

    # Pre-warp frequencies
    if not analog:
        passb = tan(wp*pi/2.0)
        stopb = tan(ws*pi/2.0)
    else:
        passb = wp*1.0
        stopb = ws*1.0

    if filter_type == 1:            # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:          # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                 args=(0,passb,stopb,gpass,gstop,'butter'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
                                 args=(1,passb,stopb,gpass,gstop,'butter'),
                                 disp=0)
        passb[1] = wp1
        nat = (stopb * (passb[0]-passb[1])) / (stopb**2 - passb[0]*passb[1])
    elif filter_type == 4:          # pass
        nat = (stopb**2 - passb[0]*passb[1]) / (stopb* (passb[0]-passb[1]))

    nat = min(abs(nat))

    GSTOP = 10**(0.1*abs(gstop))
    GPASS = 10**(0.1*abs(gpass))
    ord = int(ceil( log10((GSTOP-1.0)/(GPASS-1.0)) / (2*log10(nat))))

    # Find the butterworth natural frequency W0 (or the "3dB" frequency")
    # to give exactly gstop at nat. W0 will be between 1 and nat
    try:
        W0 = nat / ( ( 10**(0.1*abs(gstop))-1)**(1.0/(2.0*ord)))
    except ZeroDivisionError:
        W0 = nat
        print "Warning, order is zero...check input parametegstop."

    # now convert this frequency back from lowpass prototype
    # to the original analog filter

    if filter_type == 1:  # low
        WN = W0*passb
    elif filter_type == 2: # high
        WN = passb / W0
    elif filter_type == 3:  # stop
        WN = numpy.zeros(2,float)
        WN[0] = ((passb[1] - passb[0]) + sqrt((passb[1] - passb[0])**2 + \
                                        4*W0**2 * passb[0] * passb[1])) / (2*W0)
        WN[1] = ((passb[1] - passb[0]) - sqrt((passb[1] - passb[0])**2 + \
                                        4*W0**2 * passb[0] * passb[1])) / (2*W0)
        WN = numpy.sort(abs(WN))
    elif filter_type == 4: # pass
        W0 = numpy.array([-W0, W0],float)
        WN = -W0 * (passb[1]-passb[0]) / 2.0 + sqrt(W0**2 / 4.0 * \
                                              (passb[1]-passb[0])**2 + \
                                              passb[0]*passb[1])
        WN = numpy.sort(abs(WN))
    else:
        raise ValueError, "Bad type."

    if not analog:
        wn = (2.0/pi)*arctan(WN)
    else:
        wn = WN

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn

def cheb1ord(wp, ws, gpass, gstop, analog=0):
    """Chebyshev type I filter order selection.

    Return the order of the lowest order digital Chebyshev Type I filter that
    loses no more than gpass dB in the passband and has at least gstop dB
    attenuation in the stopband.

    Parameters
    ----------
    wp, ws -- Passband and stopband edge frequencies, normalized from 0
              to 1 (1 corresponds to pi radians / sample).  For example:
                 Lowpass:   wp = 0.2,          ws = 0.3
                 Highpass:  wp = 0.3,          ws = 0.2
                 Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
                 Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]
    gpass -- The maximum loss in the passband (dB).
    gstop -- The minimum attenuation in the stopband (dB).
    analog -- Non-zero to design an analog filter (in this case wp and
              ws are in radians / second).

    Returns
    -------
    ord -- The lowest order for a Chebyshev type I filter that meets specs.
    Wn -- The Chebyshev natural frequency (the "3dB frequency") for
          use with scipy.signal.cheby1 to give filter results.

    """
    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2*(len(wp)-1)
    if wp[0] < ws[0]:
        filter_type += 1
    else:
        filter_type += 2

    # Pre-wagpass frequencies
    if not analog:
        passb = tan(pi*wp/2.)
        stopb = tan(pi*ws/2.)
    else:
        passb = wp*1.0
        stopb = ws*1.0

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                 args=(0,passb,stopb,gpass,gstop,'cheby'), disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
                                 args=(1,passb,stopb,gpass,gstop,'cheby'), disp=0)
        passb[1] = wp1
        nat = (stopb * (passb[0]-passb[1])) / (stopb**2 - passb[0]*passb[1])
    elif filter_type == 4:  # pass
        nat = (stopb**2 - passb[0]*passb[1]) / (stopb* (passb[0]-passb[1]))

    nat = min(abs(nat))

    GSTOP = 10**(0.1*abs(gstop))
    GPASS = 10**(0.1*abs(gpass))
    ord = int(ceil(arccosh(sqrt((GSTOP-1.0) / (GPASS-1.0))) / arccosh(nat)))

    # Natural frequencies are just the passband edges
    if not analog:
        wn = (2.0/pi)*arctan(passb)
    else:
        wn = passb

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn


def cheb2ord(wp, ws, gpass, gstop, analog=0):
    """Chebyshev type II filter order selection.

    Description:

      Return the order of the lowest order digital Chebyshev Type II filter
      that loses no more than gpass dB in the passband and has at least gstop dB
      attenuation in the stopband.

    Parameters
    ----------
    wp, ws -- Passband and stopband edge frequencies, normalized from 0
              to 1 (1 corresponds to pi radians / sample).  For example:
                 Lowpass:   wp = 0.2,          ws = 0.3
                 Highpass:  wp = 0.3,          ws = 0.2
                 Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
                 Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]
    gpass -- The maximum loss in the passband (dB).
    gstop -- The minimum attenuation in the stopband (dB).
    analog -- Non-zero to design an analog filter (in this case wp and
              ws are in radians / second).

    Returns
    -------
    ord -- The lowest order for a Chebyshev type II filter that meets specs.
    Wn -- The Chebyshev natural frequency for
          use with scipy.signal.cheby2 to give the filter.

    """
    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2*(len(wp)-1)
    if wp[0] < ws[0]:
        filter_type += 1
    else:
        filter_type += 2

    # Pre-wagpass frequencies
    if not analog:
        passb = tan(pi*wp/2.0)
        stopb = tan(pi*ws/2.0)
    else:
        passb = wp*1.0
        stopb = ws*1.0

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                 args=(0,passb,stopb,gpass,gstop,'cheby'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
                                 args=(1,passb,stopb,gpass,gstop,'cheby'),
                                 disp=0)
        passb[1] = wp1
        nat = (stopb * (passb[0]-passb[1])) / (stopb**2 - passb[0]*passb[1])
    elif filter_type == 4:  # pass
        nat = (stopb**2 - passb[0]*passb[1]) / (stopb* (passb[0]-passb[1]))

    nat = min(abs(nat))

    GSTOP = 10**(0.1*abs(gstop))
    GPASS = 10**(0.1*abs(gpass))
    ord = int(ceil(arccosh(sqrt((GSTOP-1.0) / (GPASS-1.0))) / arccosh(nat)))

    # Find frequency where analog response is -gpass dB.
    # Then convert back from low-pass prototype to the original filter.

    new_freq = cosh(1.0/ord * arccosh(sqrt((GSTOP-1.0)/(GPASS-1.0))))
    new_freq = 1.0 / new_freq

    if filter_type == 1:
        nat = passb / new_freq
    elif filter_type == 2:
        nat = passb * new_freq
    elif filter_type == 3:
        nat = numpy.zeros(2,float)
        nat[0] = new_freq / 2.0 * (passb[0]-passb[1]) + \
                 sqrt(new_freq**2 * (passb[1]-passb[0])**2 / 4.0 + \
                      passb[1] * passb[0])
        nat[1] = passb[1] * passb[0] / nat[0]
    elif filter_type == 4:
        nat = numpy.zeros(2,float)
        nat[0] = 1.0/(2.0*new_freq) * (passb[0] - passb[1]) + \
                 sqrt((passb[1]-passb[0])**2 / (4.0*new_freq**2) + \
                      passb[1] * passb[0])
        nat[1] = passb[0] * passb[1] / nat[0]

    if not analog:
        wn = (2.0/pi)*arctan(nat)
    else:
        wn = nat

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn


def ellipord(wp, ws, gpass, gstop, analog=0):
    """Elliptic (Cauer) filter order selection.

    Return the order of the lowest order digital elliptic filter that loses no
    more than gpass dB in the passband and has at least gstop dB attenuation in
    the stopband.

    Parameters
    ----------
    wp, ws -- Passband and stopband edge frequencies, normalized from 0
              to 1 (1 corresponds to pi radians / sample).  For example:
                 Lowpass:   wp = 0.2,          ws = 0.3
                 Highpass:  wp = 0.3,          ws = 0.2
                 Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
                 Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]
    gpass -- The maximum loss in the passband (dB).
    gstop -- The minimum attenuation in the stopband (dB).
    analog -- Non-zero to design an analog filter (in this case wp and
              ws are in radians / second).

    Returns
    -------
    ord -- The lowest order for an Elliptic (Cauer) filter that meets specs.
    Wn  -- The natural frequency for use with scipy.signal.ellip
             to give the filter.

    """
    wp = atleast_1d(wp)
    ws = atleast_1d(ws)
    filter_type = 2*(len(wp)-1)
    filter_type += 1
    if wp[0] >= ws[0]:
        filter_type += 1

    # Pre-wagpass frequencies
    if analog:
        passb = wp*1.0
        stopb = ws*1.0
    else:
        passb = tan(wp*pi/2.0)
        stopb = tan(ws*pi/2.0)

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                 args=(0,passb,stopb,gpass,gstop,'ellip'),
                                 disp=0)
        passb[0] = wp0
        wp1 = optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
                                 args=(1,passb,stopb,gpass,gstop,'ellip'),
                                 disp=0)
        passb[1] = wp1
        nat = (stopb * (passb[0]-passb[1])) / (stopb**2 - passb[0]*passb[1])
    elif filter_type == 4:  # pass
        nat = (stopb**2 - passb[0]*passb[1]) / (stopb* (passb[0]-passb[1]))

    nat = min(abs(nat))

    GSTOP = 10**(0.1*gstop)
    GPASS = 10**(0.1*gpass)
    arg1 = sqrt( (GPASS-1.0) / (GSTOP-1.0) )
    arg0 = 1.0 / nat
    d0 = special.ellipk([arg0**2, 1-arg0**2])
    d1 = special.ellipk([arg1**2, 1-arg1**2])
    ord = int(ceil(d0[0]*d1[1] / (d0[1]*d1[0])))

    if not analog:
        wn = arctan(passb)*2.0/pi
    else:
        wn = passb

    if len(wn) == 1:
        wn = wn[0]
    return ord, wn

def buttap(N):
    """Return (z,p,k) zero, pole, gain for analog prototype of an Nth
    order Butterworth filter."""
    z = []
    n = numpy.arange(1,N+1)
    p = numpy.exp(1j*(2*n-1)/(2.0*N)*pi)*1j
    k = 1
    return z, p, k

def cheb1ap(N,rp):
    """Return (z,p,k) zero, pole, gain for Nth order Chebyshev type I
    lowpass analog filter prototype with rp decibels of ripple
    in the passband.
    """
    z = []
    eps = numpy.sqrt(10**(0.1*rp)-1.0)
    n = numpy.arange(1,N+1)
    mu = 1.0/N * numpy.log((1.0+numpy.sqrt(1+eps*eps)) / eps)
    theta = pi/2.0 * (2*n-1.0)/N
    p = -numpy.sinh(mu)*numpy.sin(theta) + 1j*numpy.cosh(mu)*numpy.cos(theta)
    k = numpy.prod(-p,axis=0).real
    if N % 2 == 0:
        k = k / sqrt((1+eps*eps))
    return z, p, k
    pass

def cheb2ap(N,rs):
    """Return (z,p,k) zero, pole, gain for Nth order Chebyshev type II
    lowpass analog filter prototype with rs decibels of ripple
    in the stopband.
    """
    de = 1.0/sqrt(10**(0.1*rs)-1)
    mu = arcsinh(1.0/de)/N

    if N % 2:
        m = N - 1
        n = numpy.concatenate((numpy.arange(1,N-1,2),numpy.arange(N+2,2*N,2)))
    else:
        m = N
        n = numpy.arange(1,2*N,2)

    z = conjugate(1j / cos(n*pi/(2.0*N)))
    p = exp(1j*(pi*numpy.arange(1,2*N,2)/(2.0*N) + pi/2.0))
    p = sinh(mu) * p.real + 1j*cosh(mu)*p.imag
    p = 1.0 / p
    k = (numpy.prod(-p,axis=0)/numpy.prod(-z,axis=0)).real
    return z, p, k


EPSILON = 2e-16

def vratio(u, ineps, mp):
    [s,c,d,phi] = special.ellipj(u,mp)
    ret = abs(ineps - s/c)
    return ret

def kratio(m, k_ratio):
    m = float(m)
    if m < 0:
        m = 0.0
    if m > 1:
        m = 1.0
    if abs(m) > EPSILON and (abs(m) + EPSILON) < 1:
        k = special.ellipk([m,1-m])
        r = k[0] / k[1] - k_ratio
    elif abs(m) > EPSILON:
        r = -k_ratio
    else:
        r = 1e20
    return abs(r)

def ellipap(N,rp,rs):
    """Return (z,p,k) zeros, poles, and gain of an Nth order normalized
    prototype elliptic analog lowpass filter with rp decibels of ripple
    in the passband and a stopband rs decibels down.

    See Chapter 12 and Chapter 5 of "Filter Design for Signal Processing",
    by Lutova, Tosic, and Evans.  This is
    """
    if N == 1:
        p = -sqrt(1.0/(10**(0.1*rp)-1.0))
        k = -p
        z = []
        return z, p, k

    eps = numpy.sqrt(10**(0.1*rp)-1)
    ck1 = eps / numpy.sqrt(10**(0.1*rs)-1)
    ck1p = numpy.sqrt(1-ck1*ck1)
    if ck1p == 1:
        raise ValueError, "Cannot design a filter with given rp and rs specifications."

    wp = 1
    val = special.ellipk([ck1*ck1,ck1p*ck1p])
    if abs(1-ck1p*ck1p) < EPSILON:
        krat = 0
    else:
        krat = N*val[0] / val[1]

    m = optimize.fmin(kratio, [0.5], args=(krat,), maxfun=250, maxiter=250,
                      disp=0)
    if m < 0 or m > 1:
        m = optimize.fminbound(kratio, 0, 1, args=(krat,), maxfun=250,
                               maxiter=250, disp=0)

    capk = special.ellipk(m)
    ws = wp / sqrt(m)
    m1 = 1-m

    j = numpy.arange(1-N%2,N,2)
    jj = len(j)

    [s,c,d,phi] = special.ellipj(j*capk/N,m*numpy.ones(jj))
    snew = numpy.compress(abs(s) > EPSILON, s,axis=-1)
    z = 1.0 / (sqrt(m)*snew)
    z = 1j*z
    z = numpy.concatenate((z,conjugate(z)))

    r = optimize.fmin(vratio, special.ellipk(m), args=(1./eps, ck1p*ck1p),
                      maxfun=250, maxiter=250, disp=0)
    v0 = capk * r / (N*val[0])

    [sv,cv,dv,phi] = special.ellipj(v0,1-m)
    p = -(c*d*sv*cv + 1j*s*dv) / (1-(d*sv)**2.0)

    if N % 2:
        newp = numpy.compress(abs(p.imag) > EPSILON*numpy.sqrt(numpy.sum(p*numpy.conjugate(p),axis=0).real), p,axis=-1)
        p = numpy.concatenate((p,conjugate(newp)))
    else:
        p = numpy.concatenate((p,conjugate(p)))

    k = (numpy.prod(-p,axis=0) / numpy.prod(-z,axis=0)).real
    if N % 2 == 0:
        k = k / numpy.sqrt((1+eps*eps))

    return z, p, k


def besselap(N):
    """Return (z,p,k) zero, pole, gain for analog prototype of an Nth
    order Bessel filter."""
    z = []
    k = 1
    if N == 0:
        p = [];
    elif N == 1:
        p = [-1]
    elif N == 2:
        p = [-.8660254037844386467637229+.4999999999999999999999996*1j,
             -.8660254037844386467637229-.4999999999999999999999996*1j]
    elif N == 3:
        p = [-.9416000265332067855971980,
             -.7456403858480766441810907-.7113666249728352680992154*1j,
             -.7456403858480766441810907+.7113666249728352680992154*1j]
    elif N == 4:
        p = [-.6572111716718829545787781-.8301614350048733772399715*1j,
             -.6572111716718829545787788+.8301614350048733772399715*1j,
             -.9047587967882449459642637-.2709187330038746636700923*1j,
             -.9047587967882449459642624+.2709187330038746636700926*1j]
    elif N == 5:
        p = [-.9264420773877602247196260,
             -.8515536193688395541722677-.4427174639443327209850002*1j,
             -.8515536193688395541722677+.4427174639443327209850002*1j,
             -.5905759446119191779319432-.9072067564574549539291747*1j,
             -.5905759446119191779319432+.9072067564574549539291747*1j]
    elif N == 6:
        p = [-.9093906830472271808050953-.1856964396793046769246397*1j,
             -.9093906830472271808050953+.1856964396793046769246397*1j,
             -.7996541858328288520243325-.5621717346937317988594118*1j,
             -.7996541858328288520243325+.5621717346937317988594118*1j,
             -.5385526816693109683073792-.9616876881954277199245657*1j,
             -.5385526816693109683073792+.9616876881954277199245657*1j]
    elif N == 7:
        p = [-.9194871556490290014311619,
             -.8800029341523374639772340-.3216652762307739398381830*1j,
             -.8800029341523374639772340+.3216652762307739398381830*1j,
             -.7527355434093214462291616-.6504696305522550699212995*1j,
             -.7527355434093214462291616+.6504696305522550699212995*1j,
             -.4966917256672316755024763-1.002508508454420401230220*1j,
             -.4966917256672316755024763+1.002508508454420401230220*1j]
    elif N == 8:
        p = [-.9096831546652910216327629-.1412437976671422927888150*1j,
             -.9096831546652910216327629+.1412437976671422927888150*1j,
             -.8473250802359334320103023-.4259017538272934994996429*1j,
             -.8473250802359334320103023+.4259017538272934994996429*1j,
             -.7111381808485399250796172-.7186517314108401705762571*1j,
             -.7111381808485399250796172+.7186517314108401705762571*1j,
             -.4621740412532122027072175-1.034388681126901058116589*1j,
             -.4621740412532122027072175+1.034388681126901058116589*1j]
    elif N == 9:
        p = [-.9154957797499037686769223,
             -.8911217017079759323183848-.2526580934582164192308115*1j,
             -.8911217017079759323183848+.2526580934582164192308115*1j,
             -.8148021112269012975514135-.5085815689631499483745341*1j,
             -.8148021112269012975514135+.5085815689631499483745341*1j,
             -.6743622686854761980403401-.7730546212691183706919682*1j,
             -.6743622686854761980403401+.7730546212691183706919682*1j,
             -.4331415561553618854685942-1.060073670135929666774323*1j,
             -.4331415561553618854685942+1.060073670135929666774323*1j]
    elif N == 10:
        p = [-.9091347320900502436826431-.1139583137335511169927714*1j,
             -.9091347320900502436826431+.1139583137335511169927714*1j,
             -.8688459641284764527921864-.3430008233766309973110589*1j,
             -.8688459641284764527921864+.3430008233766309973110589*1j,
             -.7837694413101441082655890-.5759147538499947070009852*1j,
             -.7837694413101441082655890+.5759147538499947070009852*1j,
             -.6417513866988316136190854-.8175836167191017226233947*1j,
             -.6417513866988316136190854+.8175836167191017226233947*1j,
             -.4083220732868861566219785-1.081274842819124562037210*1j,
             -.4083220732868861566219785+1.081274842819124562037210*1j]
    elif N == 11:
        p = [-.9129067244518981934637318,
             -.8963656705721166099815744-.2080480375071031919692341*1j
             -.8963656705721166099815744+.2080480375071031919692341*1j,
             -.8453044014712962954184557-.4178696917801248292797448*1j,
             -.8453044014712962954184557+.4178696917801248292797448*1j,
             -.7546938934722303128102142-.6319150050721846494520941*1j,
             -.7546938934722303128102142+.6319150050721846494520941*1j,
             -.6126871554915194054182909-.8547813893314764631518509*1j,
             -.6126871554915194054182909+.8547813893314764631518509*1j,
             -.3868149510055090879155425-1.099117466763120928733632*1j,
             -.3868149510055090879155425+1.099117466763120928733632*1j]
    elif N == 12:
        p = [-.9084478234140682638817772-95506365213450398415258360.0e-27*1j,
             -.9084478234140682638817772+95506365213450398415258360.0e-27*1j,
             -.8802534342016826507901575-.2871779503524226723615457*1j,
             -.8802534342016826507901575+.2871779503524226723615457*1j,
             -.8217296939939077285792834-.4810212115100676440620548*1j,
             -.8217296939939077285792834+.4810212115100676440620548*1j,
             -.7276681615395159454547013-.6792961178764694160048987*1j,
             -.7276681615395159454547013+.6792961178764694160048987*1j,
             -.5866369321861477207528215-.8863772751320727026622149*1j,
             -.5866369321861477207528215+.8863772751320727026622149*1j,
             -.3679640085526312839425808-1.114373575641546257595657*1j,
             -.3679640085526312839425808+1.114373575641546257595657*1j]
    elif N == 13:
        p = [-.9110914665984182781070663,
             -.8991314665475196220910718-.1768342956161043620980863*1j,
             -.8991314665475196220910718+.1768342956161043620980863*1j,
             -.8625094198260548711573628-.3547413731172988997754038*1j,
             -.8625094198260548711573628+.3547413731172988997754038*1j,
             -.7987460692470972510394686-.5350752120696801938272504*1j,
             -.7987460692470972510394686+.5350752120696801938272504*1j,
             -.7026234675721275653944062-.7199611890171304131266374*1j,
             -.7026234675721275653944062+.7199611890171304131266374*1j,
             -.5631559842430199266325818-.9135900338325109684927731*1j,
             -.5631559842430199266325818+.9135900338325109684927731*1j,
             -.3512792323389821669401925-1.127591548317705678613239*1j,
             -.3512792323389821669401925+1.127591548317705678613239*1j]
    elif N == 14:
        p = [-.9077932138396487614720659-82196399419401501888968130.0e-27*1j,
             -.9077932138396487614720659+82196399419401501888968130.0e-27*1j,
             -.8869506674916445312089167-.2470079178765333183201435*1j,
             -.8869506674916445312089167+.2470079178765333183201435*1j,
             -.8441199160909851197897667-.4131653825102692595237260*1j,
             -.8441199160909851197897667+.4131653825102692595237260*1j,
             -.7766591387063623897344648-.5819170677377608590492434*1j,
             -.7766591387063623897344648+.5819170677377608590492434*1j,
             -.6794256425119233117869491-.7552857305042033418417492*1j,
             -.6794256425119233117869491+.7552857305042033418417492*1j,
             -.5418766775112297376541293-.9373043683516919569183099*1j,
             -.5418766775112297376541293+.9373043683516919569183099*1j,
             -.3363868224902037330610040-1.139172297839859991370924*1j,
             -.3363868224902037330610040+1.139172297839859991370924*1j]
    elif N == 15:
        p = [-.9097482363849064167228581,
             -.9006981694176978324932918-.1537681197278439351298882*1j,
             -.9006981694176978324932918+.1537681197278439351298882*1j,
             -.8731264620834984978337843-.3082352470564267657715883*1j,
             -.8731264620834984978337843+.3082352470564267657715883*1j,
             -.8256631452587146506294553-.4642348752734325631275134*1j,
             -.8256631452587146506294553+.4642348752734325631275134*1j,
             -.7556027168970728127850416-.6229396358758267198938604*1j,
             -.7556027168970728127850416+.6229396358758267198938604*1j,
             -.6579196593110998676999362-.7862895503722515897065645*1j,
             -.6579196593110998676999362+.7862895503722515897065645*1j,
             -.5224954069658330616875186-.9581787261092526478889345*1j,
             -.5224954069658330616875186+.9581787261092526478889345*1j,
             -.3229963059766444287113517-1.149416154583629539665297*1j,
             -.3229963059766444287113517+1.149416154583629539665297*1j]
    elif N == 16:
        p = [-.9072099595087001356491337-72142113041117326028823950.0e-27*1j,
             -.9072099595087001356491337+72142113041117326028823950.0e-27*1j,
             -.8911723070323647674780132-.2167089659900576449410059*1j,
             -.8911723070323647674780132+.2167089659900576449410059*1j,
             -.8584264231521330481755780-.3621697271802065647661080*1j,
             -.8584264231521330481755780+.3621697271802065647661080*1j,
             -.8074790293236003885306146-.5092933751171800179676218*1j,
             -.8074790293236003885306146+.5092933751171800179676218*1j,
             -.7356166304713115980927279-.6591950877860393745845254*1j,
             -.7356166304713115980927279+.6591950877860393745845254*1j,
             -.6379502514039066715773828-.8137453537108761895522580*1j,
             -.6379502514039066715773828+.8137453537108761895522580*1j,
             -.5047606444424766743309967-.9767137477799090692947061*1j,
             -.5047606444424766743309967+.9767137477799090692947061*1j,
             -.3108782755645387813283867-1.158552841199330479412225*1j,
             -.3108782755645387813283867+1.158552841199330479412225*1j]
    elif N == 17:
        p = [-.9087141161336397432860029,
             -.9016273850787285964692844-.1360267995173024591237303*1j,
             -.9016273850787285964692844+.1360267995173024591237303*1j,
             -.8801100704438627158492165-.2725347156478803885651973*1j,
             -.8801100704438627158492165+.2725347156478803885651973*1j,
             -.8433414495836129204455491-.4100759282910021624185986*1j,
             -.8433414495836129204455491+.4100759282910021624185986*1j,
             -.7897644147799708220288138-.5493724405281088674296232*1j,
             -.7897644147799708220288138+.5493724405281088674296232*1j,
             -.7166893842372349049842743-.6914936286393609433305754*1j,
             -.7166893842372349049842743+.6914936286393609433305754*1j,
             -.6193710717342144521602448-.8382497252826992979368621*1j,
             -.6193710717342144521602448+.8382497252826992979368621*1j,
             -.4884629337672704194973683-.9932971956316781632345466*1j,
             -.4884629337672704194973683+.9932971956316781632345466*1j,
             -.2998489459990082015466971-1.166761272925668786676672*1j,
             -.2998489459990082015466971+1.166761272925668786676672*1j]
    elif N == 18:
        p = [-.9067004324162775554189031-64279241063930693839360680.0e-27*1j,
             -.9067004324162775554189031+64279241063930693839360680.0e-27*1j,
             -.8939764278132455733032155-.1930374640894758606940586*1j,
             -.8939764278132455733032155+.1930374640894758606940586*1j,
             -.8681095503628830078317207-.3224204925163257604931634*1j,
             -.8681095503628830078317207+.3224204925163257604931634*1j,
             -.8281885016242836608829018-.4529385697815916950149364*1j,
             -.8281885016242836608829018+.4529385697815916950149364*1j,
             -.7726285030739558780127746-.5852778162086640620016316*1j,
             -.7726285030739558780127746+.5852778162086640620016316*1j,
             -.6987821445005273020051878-.7204696509726630531663123*1j,
             -.6987821445005273020051878+.7204696509726630531663123*1j,
             -.6020482668090644386627299-.8602708961893664447167418*1j,
             -.6020482668090644386627299+.8602708961893664447167418*1j,
             -.4734268069916151511140032-1.008234300314801077034158*1j,
             -.4734268069916151511140032+1.008234300314801077034158*1j,
             -.2897592029880489845789953-1.174183010600059128532230*1j,
             -.2897592029880489845789953+1.174183010600059128532230*1j]
    elif N == 19:
        p = [-.9078934217899404528985092,
             -.9021937639390660668922536-.1219568381872026517578164*1j,
             -.9021937639390660668922536+.1219568381872026517578164*1j,
             -.8849290585034385274001112-.2442590757549818229026280*1j,
             -.8849290585034385274001112+.2442590757549818229026280*1j,
             -.8555768765618421591093993-.3672925896399872304734923*1j,
             -.8555768765618421591093993+.3672925896399872304734923*1j,
             -.8131725551578197705476160-.4915365035562459055630005*1j,
             -.8131725551578197705476160+.4915365035562459055630005*1j,
             -.7561260971541629355231897-.6176483917970178919174173*1j,
             -.7561260971541629355231897+.6176483917970178919174173*1j,
             -.6818424412912442033411634-.7466272357947761283262338*1j,
             -.6818424412912442033411634+.7466272357947761283262338*1j,
             -.5858613321217832644813602-.8801817131014566284786759*1j,
             -.5858613321217832644813602+.8801817131014566284786759*1j,
             -.4595043449730988600785456-1.021768776912671221830298*1j,
             -.4595043449730988600785456+1.021768776912671221830298*1j,
             -.2804866851439370027628724-1.180931628453291873626003*1j,
             -.2804866851439370027628724+1.180931628453291873626003*1j]
    elif N == 20:
        p = [-.9062570115576771146523497-57961780277849516990208850.0e-27*1j,
             -.9062570115576771146523497+57961780277849516990208850.0e-27*1j,
             -.8959150941925768608568248-.1740317175918705058595844*1j,
             -.8959150941925768608568248+.1740317175918705058595844*1j,
             -.8749560316673332850673214-.2905559296567908031706902*1j,
             -.8749560316673332850673214+.2905559296567908031706902*1j,
             -.8427907479956670633544106-.4078917326291934082132821*1j,
             -.8427907479956670633544106+.4078917326291934082132821*1j,
             -.7984251191290606875799876-.5264942388817132427317659*1j,
             -.7984251191290606875799876+.5264942388817132427317659*1j,
             -.7402780309646768991232610-.6469975237605228320268752*1j,
             -.7402780309646768991232610+.6469975237605228320268752*1j,
             -.6658120544829934193890626-.7703721701100763015154510*1j,
             -.6658120544829934193890626+.7703721701100763015154510*1j,
             -.5707026806915714094398061-.8982829066468255593407161*1j,
             -.5707026806915714094398061+.8982829066468255593407161*1j,
             -.4465700698205149555701841-1.034097702560842962315411*1j,
             -.4465700698205149555701841+1.034097702560842962315411*1j,
             -.2719299580251652601727704-1.187099379810885886139638*1j,
             -.2719299580251652601727704+1.187099379810885886139638*1j]
    elif N == 21:
        p = [-.9072262653142957028884077,
             -.9025428073192696303995083-.1105252572789856480992275*1j,
             -.9025428073192696303995083+.1105252572789856480992275*1j,
             -.8883808106664449854431605-.2213069215084350419975358*1j,
             -.8883808106664449854431605+.2213069215084350419975358*1j,
             -.8643915813643204553970169-.3326258512522187083009453*1j,
             -.8643915813643204553970169+.3326258512522187083009453*1j,
             -.8299435470674444100273463-.4448177739407956609694059*1j,
             -.8299435470674444100273463+.4448177739407956609694059*1j,
             -.7840287980408341576100581-.5583186348022854707564856*1j,
             -.7840287980408341576100581+.5583186348022854707564856*1j,
             -.7250839687106612822281339-.6737426063024382240549898*1j,
             -.7250839687106612822281339+.6737426063024382240549898*1j,
             -.6506315378609463397807996-.7920349342629491368548074*1j,
             -.6506315378609463397807996+.7920349342629491368548074*1j,
             -.5564766488918562465935297-.9148198405846724121600860*1j,
             -.5564766488918562465935297+.9148198405846724121600860*1j,
             -.4345168906815271799687308-1.045382255856986531461592*1j,
             -.4345168906815271799687308+1.045382255856986531461592*1j,
             -.2640041595834031147954813-1.192762031948052470183960*1j,
             -.2640041595834031147954813+1.192762031948052470183960*1j]
    elif N == 22:
        p = [-.9058702269930872551848625-52774908289999045189007100.0e-27*1j,
             -.9058702269930872551848625+52774908289999045189007100.0e-27*1j,
             -.8972983138153530955952835-.1584351912289865608659759*1j,
             -.8972983138153530955952835+.1584351912289865608659759*1j,
             -.8799661455640176154025352-.2644363039201535049656450*1j,
             -.8799661455640176154025352+.2644363039201535049656450*1j,
             -.8534754036851687233084587-.3710389319482319823405321*1j,
             -.8534754036851687233084587+.3710389319482319823405321*1j,
             -.8171682088462720394344996-.4785619492202780899653575*1j,
             -.8171682088462720394344996+.4785619492202780899653575*1j,
             -.7700332930556816872932937-.5874255426351153211965601*1j,
             -.7700332930556816872932937+.5874255426351153211965601*1j,
             -.7105305456418785989070935-.6982266265924524000098548*1j,
             -.7105305456418785989070935+.6982266265924524000098548*1j,
             -.6362427683267827226840153-.8118875040246347267248508*1j,
             -.6362427683267827226840153+.8118875040246347267248508*1j,
             -.5430983056306302779658129-.9299947824439872998916657*1j,
             -.5430983056306302779658129+.9299947824439872998916657*1j,
             -.4232528745642628461715044-1.055755605227545931204656*1j,
             -.4232528745642628461715044+1.055755605227545931204656*1j,
             -.2566376987939318038016012-1.197982433555213008346532*1j,
             -.2566376987939318038016012+1.197982433555213008346532*1j]
    elif N == 23:
        p = [-.9066732476324988168207439,
             -.9027564979912504609412993-.1010534335314045013252480*1j,
             -.9027564979912504609412993+.1010534335314045013252480*1j,
             -.8909283242471251458653994-.2023024699381223418195228*1j,
             -.8909283242471251458653994+.2023024699381223418195228*1j,
             -.8709469395587416239596874-.3039581993950041588888925*1j,
             -.8709469395587416239596874+.3039581993950041588888925*1j,
             -.8423805948021127057054288-.4062657948237602726779246*1j,
             -.8423805948021127057054288+.4062657948237602726779246*1j,
             -.8045561642053176205623187-.5095305912227258268309528*1j,
             -.8045561642053176205623187+.5095305912227258268309528*1j,
             -.7564660146829880581478138-.6141594859476032127216463*1j,
             -.7564660146829880581478138+.6141594859476032127216463*1j,
             -.6965966033912705387505040-.7207341374753046970247055*1j,
             -.6965966033912705387505040+.7207341374753046970247055*1j,
             -.6225903228771341778273152-.8301558302812980678845563*1j,
             -.6225903228771341778273152+.8301558302812980678845563*1j,
             -.5304922463810191698502226-.9439760364018300083750242*1j,
             -.5304922463810191698502226+.9439760364018300083750242*1j,
             -.4126986617510148836149955-1.065328794475513585531053*1j,
             -.4126986617510148836149955+1.065328794475513585531053*1j,
             -.2497697202208956030229911-1.202813187870697831365338*1j,
             -.2497697202208956030229911+1.202813187870697831365338*1j]
    elif N == 24:
        p = [-.9055312363372773709269407-48440066540478700874836350.0e-27*1j,
             -.9055312363372773709269407+48440066540478700874836350.0e-27*1j,
             -.8983105104397872954053307-.1454056133873610120105857*1j,
             -.8983105104397872954053307+.1454056133873610120105857*1j,
             -.8837358034555706623131950-.2426335234401383076544239*1j,
             -.8837358034555706623131950+.2426335234401383076544239*1j,
             -.8615278304016353651120610-.3403202112618624773397257*1j,
             -.8615278304016353651120610+.3403202112618624773397257*1j,
             -.8312326466813240652679563-.4386985933597305434577492*1j,
             -.8312326466813240652679563+.4386985933597305434577492*1j,
             -.7921695462343492518845446-.5380628490968016700338001*1j,
             -.7921695462343492518845446+.5380628490968016700338001*1j,
             -.7433392285088529449175873-.6388084216222567930378296*1j,
             -.7433392285088529449175873+.6388084216222567930378296*1j,
             -.6832565803536521302816011-.7415032695091650806797753*1j,
             -.6832565803536521302816011+.7415032695091650806797753*1j,
             -.6096221567378335562589532-.8470292433077202380020454*1j,
             -.6096221567378335562589532+.8470292433077202380020454*1j,
             -.5185914574820317343536707-.9569048385259054576937721*1j,
             -.5185914574820317343536707+.9569048385259054576937721*1j,
             -.4027853855197518014786978-1.074195196518674765143729*1j,
             -.4027853855197518014786978+1.074195196518674765143729*1j,
             -.2433481337524869675825448-1.207298683731972524975429*1j,
             -.2433481337524869675825448+1.207298683731972524975429*1j]
    elif N == 25:
        p = [-.9062073871811708652496104,
             -.9028833390228020537142561-93077131185102967450643820.0e-27*1j,
             -.9028833390228020537142561+93077131185102967450643820.0e-27*1j,
             -.8928551459883548836774529-.1863068969804300712287138*1j,
             -.8928551459883548836774529+.1863068969804300712287138*1j,
             -.8759497989677857803656239-.2798521321771408719327250*1j,
             -.8759497989677857803656239+.2798521321771408719327250*1j,
             -.8518616886554019782346493-.3738977875907595009446142*1j,
             -.8518616886554019782346493+.3738977875907595009446142*1j,
             -.8201226043936880253962552-.4686668574656966589020580*1j,
             -.8201226043936880253962552+.4686668574656966589020580*1j,
             -.7800496278186497225905443-.5644441210349710332887354*1j,
             -.7800496278186497225905443+.5644441210349710332887354*1j,
             -.7306549271849967721596735-.6616149647357748681460822*1j,
             -.7306549271849967721596735+.6616149647357748681460822*1j,
             -.6704827128029559528610523-.7607348858167839877987008*1j,
             -.6704827128029559528610523+.7607348858167839877987008*1j,
             -.5972898661335557242320528-.8626676330388028512598538*1j,
             -.5972898661335557242320528+.8626676330388028512598538*1j,
             -.5073362861078468845461362-.9689006305344868494672405*1j,
             -.5073362861078468845461362+.9689006305344868494672405*1j,
             -.3934529878191079606023847-1.082433927173831581956863*1j,
             -.3934529878191079606023847+1.082433927173831581956863*1j,
             -.2373280669322028974199184-1.211476658382565356579418*1j,
             -.2373280669322028974199184+1.211476658382565356579418*1j]
    else:
        raise ValueError, "Bessel Filter not supported for order %d" % N

    return z, p, k

filter_dict = {'butter': [buttap,buttord],
               'butterworth' : [buttap,buttord],
               'cauer' : [ellipap,ellipord],
               'elliptic' : [ellipap,ellipord],
               'ellip' : [ellipap,ellipord],
               'bessel' : [besselap],
               'cheby1' : [cheb1ap, cheb1ord],
               'chebyshev1' : [cheb1ap, cheb1ord],
               'chebyshevi' : [cheb1ap, cheb1ord],
               'cheby2' : [cheb2ap, cheb2ord],
               'chebyshev2' : [cheb2ap, cheb2ord],
               'chebyshevii' : [cheb2ap, cheb2ord]
               }

band_dict = {'band':'bandpass',
             'bandpass':'bandpass',
             'pass' : 'bandpass',
             'bp':'bandpass',
             'bs':'bandstop',
             'bandstop':'bandstop',
             'bands' : 'bandstop',
             'stop' : 'bandstop',
             'l' : 'lowpass',
             'low': 'lowpass',
             'lowpass' : 'lowpass',
             'high' : 'highpass',
             'highpass' : 'highpass',
             'h' : 'highpass'
             }

def kaiserord(ripple, width):
    """Design a Kaiser window to limit ripple and width of transition region.

    Parameters
    ----------
    ripple -- positive number specifying maximum ripple in passband (dB)
                and minimum ripple in stopband
    width  -- width of transition region (normalized so that 1 corresponds
                to pi radians / sample)

    Returns
    -------
    N, beta -- the order and beta parameter for the kaiser window.

               signal.kaiser(N,beta,sym=0) returns the window as does
               signal.get_window(beta,N)
               signal.get_window(('kaiser',beta),N)

    Uses the empirical equations discovered by Kaiser.

    Oppenheim, Schafer, "Discrete-Time Signal Processing,", p.475-476.
    """
    A = abs(ripple)  # in case somebody is confused as to what's meant
    if (A>50):
        beta = 0.1102*(A-8.7)
    elif (A>21):
        beta = 0.5842*(A-21)**0.4 + 0.07886*(A-21)
    else:
        beta = 0.0
    N = (A-8)/2.285/(pi*width)
    return ceil(N), beta

def firwin(N, cutoff, width=None, window='hamming'):
    """
    FIR Filter Design using windowed ideal filter method.

    Parameters
    ----------
    N      -- order of filter (number of taps)
    cutoff -- cutoff frequency of filter (normalized so that 1 corresponds to
              Nyquist or pi radians / sample)

    width  -- if width is not None, then assume it is the approximate width of
              the transition region (normalized so that 1 corresonds to pi)
              for use in kaiser FIR filter design.
    window -- desired window to use. See get_window for a list
              of windows and required parameters.

    Returns
    -------
    h      -- coefficients of length N fir filter.

    """

    from signaltools import get_window
    if isinstance(width,float):
        A = 2.285*N*width + 8
        if (A < 21): beta = 0.0
        elif (A <= 50): beta = 0.5842*(A-21)**0.4 + 0.07886*(A-21)
        else: beta = 0.1102*(A-8.7)
        window=('kaiser',beta)

    win = get_window(window,N,fftbins=1)
    alpha = N//2
    m = numpy.arange(0,N)
    h = win*special.sinc(cutoff*(m-alpha))
    return h / numpy.sum(h,axis=0)

warnings.simplefilter("always", BadCoefficients)
