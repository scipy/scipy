import MLab
import scipy
from scipy.fastumath import *
Num = MLab
abs = absolute
pi = Numeric.pi
import scipy
from scipy import r1array, poly

def factorial(n):
    return scipy.special.gamma(n+1)

def comb(N,k):
    lgam = scipy.special.lgam
    return exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))
    
def tf2zpk(b,a):
    b = (b+0.0) / a[0]
    a = (a+0.0) / a[0]
    k = b[0]
    b /= b[0]
    z = roots(b)
    p = roots(a)
    return z, p, k
    
def zpk2tf(z,p,k):
    """Return polynomial transfer function representation from zeros and poles

    Inputs:

      z, p --- sequences representing the zeros and poles.
      k --- system gain.

    Outputs: (b,a)

      b, a --- numerator and denominator polynomials.
    """
    b = k * poly(z)
    a = poly(p)
    return b, a

def normalize(b,a):
    if a[0] != 0.0:
        outb = b / a[0]
        outa = a / a[0]
    else:
        outb, outa = b, a
    return outb, outa


def lp2lp(b,a,wo=1.0):
    """Return a low-pass filter with cuttoff frequency wo
    from a low-pass filter prototype with unity cutoff frequency.
    """
    a,b = map(r1array,(a,b))
    wo = float(wo)
    d = len(a)
    n = len(b)
    M = max((d,n))
    pwo = pow(wo,arange(M-1,-1,-1))
    start1 = max((n-d,0))
    start2 = max((d-n,0))
    b = b / pwo[start2:]
    a = a / pwo[start1:]
    return normalize(b, a)

def lp2hp(b,a,wo=1):
    """Return a high-pass filter with cuttoff frequency wo
    from a low-pass filter prototype with unity cutoff frequency.
    """
    a,b = map(r1array,(a,b))
    d = len(a)
    n = len(b)
    if wo != 1:
        pwo = pow(wo,arange(max((d,n))))
    else:
        pwo = ones(max((d,n)),b.typecode())
    if d >= n:
        outa = a[::-1] * pwo
        outb = Numeric.resize(b,d)
        outb[:n] = b[::-1] * pwo[:n]
    else:
        outb = b[::-1] * pwo
        outa = Numeric.resize(a,n)
        outa[:d] = a[::-1] * pwo[:d]

    return normalize(outb, outa)

def lp2bp(b,a,wo=1.0, bw=1.0):
    """Return a band-pass filter with center frequency wo and bandwidth bw
    from a low-pass filter prototype with unity cutoff frequency.
    """
    a,b = map(r1array,(a,b))
    D = len(a) - 1
    N = len(b) - 1
    artype = b.typecode()
    if artype not in ['F','D','f','d']:
        artype = Float
    ma = max([N,D])
    Np = N + ma
    Dp = D + ma
    bprime = zeros(Np+1,artype)
    aprime = zeros(Dp+1,artype)
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
    a,b = map(r1array,(a,b))
    D = len(a) - 1
    N = len(b) - 1
    artype = b.typecode()
    if artype not in ['F','D','f','d']:
        artype = Float
    M = max([N,D])
    Np = M + M
    Dp = M + M
    bprime = zeros(Np+1,artype)
    aprime = zeros(Dp+1,artype)
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
    a,b = map(r1array,(a,b))
    D = len(a) - 1
    N = len(b) - 1
    artype = Float
    M = max([N,D])
    Np = M
    Dp = M
    bprime = zeros(Np+1,artype)
    aprime = zeros(Dp+1,artype)
    for j in range(Np+1):
        val = 0.0
        for i in range(N+1):
             for k in range(i+1):
                for l in range(M-i+1):
                    if k+l == j:
                        val += comb(i,k)*comb(M-i,l)*b[N-i]*pow(2*fs,i)*(-1)**k
        bprime[j] = val
    for j in range(Dp+1):
        val = 0.0
        for i in range(D+1):
            for k in range(i+1):
                for l in range(M-i+1):
                    if k+l == j:
                        val += comb(i,k)*comb(M-i,l)*a[D-i]*pow(2*fs,i)*(-1)**k
        aprime[j] = val
        
    return normalize(bprime, aprime)


def butter(N, Wn, bandtype='band', analog=0, output=''):
    """Butterworth digital and analog filter design.

    Description:

      Design an Nth order lowpass digital Butterworth filter and return the
      filter coefficients in (B,A) form.

    Inputs:
      
    """

    #pre-warp frequencies for digital filter design
    if not analog:
        warped = 2*fs*tan(pi*Wn/fs)
    else:
        warped = Wn

    # convert to low-pass prototype

    # Get analog lowpass prototype

    # transform to lowpass, bandpass, highpass, or bandstop

    # Find discrete equivalent if necessary
    if not analog:
        pass

    # Transform to proper out type (pole-zero, state-space, numer-denom)    
    pass

def cheby1():
    pass

def cheby2():
    pass

def ellip():
    pass

def besself():
    pass

def maxflat():
    pass

def yulewalk():
    pass


def band_stop_obj(wp, ind, passb, stopb, gpass, gstop, type):
    """Band Stop Objective Function for order minimization

    Description:

      Returns the non-integer order for an analog band stop filter.

    Inputs:

      wp -- passb edge
      ind -- index specifying which passb edge to vary (0 or 1).
      passb -- two element vector of fixed passband edges.
      stopb -- two element vector of fixed stopband edges.
      gstop -- amount in dB of attenuation in stopband.
      gpass -- amount in dB of ripple in the passband.
      type -- 'butter', 'cheby', or 'ellip':

    Outputs: (n,)

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
        d0 = scipy.special.ellpk([1-arg0**2, arg0**2])
        d1 = scipy.special.ellpk([1-arg1**2, arg1**2])
        n = (d0[0]*d1[1] / (d0[1]*d1[0]))
    else:
        raise ValueError, "Incorrect type: ", type
    return n

def buttord(wp, ws, gpass, gstop, analog=0):
    """Butterworth filter order selection.

    Description:

      Return the order of the lowest order digital Butterworth filter that
      loses no more than gpass dB in the passband and has at least gstop dB
      attenuation in the stopband.

    Inputs:

      wp, ws -- Passb and stopb edge frequencies, normalized from 0
                to 1 (1 corresponds to pi radians / sample).  For example:
                   Lowpass:   wp = 0.2,          ws = 0.3
                   Highpass:  wp = 0.3,          ws = 0.2
                   Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
                   Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]
      gpass -- The maximum loss in the passband (dB).
      gstop -- The minimum attenuation in the stopband (dB).
      analog -- Non-zero to design an analog filter (in this case wp and
                ws are in radians / second).

    Outputs: (ord, Wn)

      ord -- The lowest order for a Butterworth filter which meets specs.
      Wn -- The Butterworth natural frequency (i.e. the "3dB frequency"). 
            Should be used with scipy.signal.butter to give filter results.

    """

    wp = r1array(wp)
    ws = r1array(ws)
    filter_type = 2*(len(wp)-1)
    filter_type +=1 
    if wp[0] >= ws[0]:
        filter_type += 1

    # Pre-warp frequencies
    if not analog:
        passb = tan(wp*pi/2.0)
        stopb = tan(ws*pi/2.0)
    else:
        passb = wp
        stopb = ws

    if filter_type == 1:            # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:          # stop
        wp0 = scipy.optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                       args=(0,passb,stopb,gpass,gstop,'butter'),
                                       disp=0)
        passb[0] = wp0
        wp1 = scipy.optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
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
        WN = zeros(2,Numeric.Float)
        WN[0] = ((passb[1] - passb[0]) + sqrt((passb[1] - passb[0])**2 + \
                                        4*W0**2 * passb[0] * passb[1])) / (2*W0)
        WN[1] = ((passb[1] - passb[0]) - sqrt((passb[1] - passb[0])**2 + \
                                        4*W0**2 * passb[0] * passb[1])) / (2*W0)
        WN = Num.sort(abs(WN))
    elif filter_type == 4: # pass
        W0 = array([-W0, W0],Numeric.Float)
        WN = -W0 * (passb[1]-passb[0]) / 2.0 + sqrt(W0**2 / 4.0 * \
                                              (passb[1]-passb[0])**2 + \
                                              passb[0]*passb[1])
        WN = Num.sort(abs(WN))
    else:
        raise ValueError, "Bad type."

    if not analog:
        wn = (2.0/pi)*arctan(WN)
    else:
        wn = WN

    return ord, wn

def cheb1ord(wp, ws, gpass, gstop, analog=0):
    """Chebyshev type I filter order selection.

    Description:

      Return the order of the lowest order digital Chebyshev Type I filter
      that loses no more than gpass dB in the passband and has at least gstop dB
      attenuation in the stopband.

    Inputs:

      wp, ws -- Passb and stopb edge frequencies, normalized from 0
                to 1 (1 corresponds to pi radians / sample).  For example:
                   Lowpass:   Wp = 0.2,          Ws = 0.3
                   Highpass:  Wp = 0.3,          Ws = 0.2
                   Bandpass:  Wp = [0.2, 0.5],   Ws = [0.1, 0.6]
                   Bandstop:  Wp = [0.1, 0.6],   Ws = [0.2, 0.5]
      gpass -- The maximum loss in the passband (dB).
      gstop -- The minimum attenuation in the stopband (dB).
      analog -- Non-zero to design an analog filter (in this case wp and
                ws are in radians / second).

    Outputs: (ord, Wn)

      ord -- The lowest order for a Chebyshev type I filter that meets specs.
      Wn -- The Chebyshev natural frequency (the "3dB frequency") for
            use with scipy.signal.cheby1 to give filter results.

    """
    wp = r1array(wp)
    ws = r1array(ws)
    filter_type = 2*(len(wp)-1)
    if wp[0] < ws[0]:
        filter_type += 1
    else:
        filter_type += 2

    # Pre-wagpass frequencies
    if not analog:
        passb = tan(pi*wp/2)
        stopb = tan(pi*ws/2)
    else:
        passb = wp
        stopb = ws

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = scipy.optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                       args=(0,passb,stopb,gpass,gstop,'cheby'), disp=0)
        passb[0] = wp0
        wp1 = scipy.optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
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

    return ord, wn
    

def cheb2ord(wp, ws, gpass, gstop, analog=0):
    """Chebyshev type II filter order selection.

    Description:

      Return the order of the lowest order digital Chebyshev Type II filter
      that loses no more than gpass dB in the passband and has at least gstop dB
      attenuation in the stopband.

    Inputs:

      wp, ws -- Passb and stopb edge frequencies, normalized from 0
                to 1 (1 corresponds to pi radians / sample).  For example:
                   Lowpass:   Wp = 0.2,          Ws = 0.3
                   Highpass:  Wp = 0.3,          Ws = 0.2
                   Bandpass:  Wp = [0.2, 0.5],   Ws = [0.1, 0.6]
                   Bandstop:  Wp = [0.1, 0.6],   Ws = [0.2, 0.5]
      gpass -- The maximum loss in the passband (dB).
      gstop -- The minimum attenuation in the stopband (dB).
      analog -- Non-zero to design an analog filter (in this case wp and
                ws are in radians / second).

    Outputs: (ord, Wn)

      ord -- The lowest order for a Chebyshev type II filter that meets specs.
      Wn -- The Chebyshev natural frequency for
            use with scipy.signal.cheby2 to give the filter.

    """
    wp = r1array(wp)
    ws = r1array(ws)
    filter_type = 2*(len(wp)-1)
    if wp[0] < ws[0]:
        filter_type += 1
    else:
        filter_type += 2

    # Pre-wagpass frequencies
    if not analog:
        passb = tan(pi*wp/2)
        stopb = tan(pi*ws/2)
    else:
        passb = wp
        stopb = ws

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = scipy.optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                       args=(0,passb,stopb,gpass,gstop,'cheby'),
                                       disp=0)
        passb[0] = wp0
        wp1 = scipy.optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
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
        nat = zeros(2,Num.Float)
        nat[0] = new_freq / 2.0 * (passb[0]-passb[1]) + \
                 sqrt(new_freq**2 * (passb[1]-passb[0])**2 / 4.0 + \
                      passb[1] * passb[0])
        nat[1] = passb[1] * passb[0] / nat[0]
    elif filter_type == 4:
        nat = zeros(2,Num.Float)
        nat[0] = 1.0/(2.0*new_freq) * (passb[0] - passb[1]) + \
                 sqrt((passb[1]-passb[0])**2 / (4.0*new_freq**2) + \
                      passb[1] * passb[0])
        nat[1] = passb[0] * passb[1] / nat[0]        

    if not analog:
        wn = (2.0/pi)*arctan(nat)
    else:
        wn = nat

    return ord, wn


def ellipord(wp, ws, gpass, gstop, analog=0):
    """Elliptic (Cauer) filter order selection.

    Description:

      Return the order of the lowest order digital elliptic filter
      that loses no more than gpass dB in the passband and has at least gstop dB
      attenuation in the stopband.

    Inputs:

      wp, ws -- Passb and stopb edge frequencies, normalized from 0
                to 1 (1 corresponds to pi radians / sample).  For example:
                   Lowpass:   Wp = 0.2,          Ws = 0.3
                   Highpass:  Wp = 0.3,          Ws = 0.2
                   Bandpass:  Wp = [0.2, 0.5],   Ws = [0.1, 0.6]
                   Bandstop:  Wp = [0.1, 0.6],   Ws = [0.2, 0.5]
      gpass -- The maximum loss in the passband (dB).
      gstop -- The minimum attenuation in the stopband (dB).
      analog -- Non-zero to design an analog filter (in this case wp and
                ws are in radians / second).

    Outputs: (ord, Wn)

      ord -- The lowest order for a Chebyshev type II filter that meets specs.
      Wn -- The Chebyshev natural frequency for
            use with scipy.signal.cheby2 to give the filter.

    """
    wp = r1array(wp)
    ws = r1array(ws)
    filter_type = 2*(len(wp)-1)
    filter_type += 1
    if wp[0] >= ws[0]:
        filter_type += 1

    # Pre-wagpass frequencies
    if analog:
        passb = wp
        stopb = ws
    else:
        passb = tan(wp*pi/2.0)
        stopb = tan(ws*pi/2.0)

    if filter_type == 1:           # low
        nat = stopb / passb
    elif filter_type == 2:          # high
        nat = passb / stopb
    elif filter_type == 3:     # stop
        wp0 = scipy.optimize.fminbound(band_stop_obj, passb[0], stopb[0]-1e-12,
                                       args=(0,passb,stopb,gpass,gstop,'ellip'),
                                       disp=0)
        passb[0] = wp0
        wp1 = scipy.optimize.fminbound(band_stop_obj, stopb[1]+1e-12, passb[1],
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
    d0 = scipy.special.ellpk([1-arg0**2, arg0**2])
    d1 = scipy.special.ellpk([1-arg1**2, arg1**2])
    ord = int(ceil(d0[0]*d1[1] / (d0[1]*d1[0])))

    if not analog:
        wn = arctan(passb)*2.0/pi
    else:
        wn = passb

    return ord, wn
    
def buttap(N):
    """Return (z,p,k) zero, pole, gain for analog prototype of an Nth
    order Butterworth filter."""
    z = []
    n = Num.arange(1,N+1)
    p = Num.exp(1j*(2*n-1)/(2.0*N)*pi)*1j
    k = 1
    return z, p, k

def cheb1ap(N,rp):
    """Return (z,p,k) zero, pole, gain for Nth order Chebyshev type I
    lowpass analog filter prototype with rp decibels of ripple
    in the passband.
    """
    z = []
    eps = Num.sqrt(10**(0.1*rp)-1.0)
    n = Num.arange(1,N+1)
    mu = 1.0/N * Num.log((1.0+Num.sqrt(1+eps*eps)) / eps)
    theta = pi/2.0 * (2*n-1.0)/N
    p = -Num.sinh(mu)*Num.sin(theta) + 1j*Num.cosh(mu)*Num.cos(theta)
    k = MLab.prod(-p).real
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
        n = Num.concatenate((arange(1,N-1,2),arange(N+2,2*N,2)))
    else:
        m = N
        n = arange(1,2*N,2)
        
    z = conjugate(1j / cos(n*pi/(2.0*N)))
    p = exp(1j*(pi*arange(1,2*N,2)/(2.0*N) + pi/2.0))
    p = sinh(mu) * p.real + 1j*cosh(mu)*p.imag
    p = 1.0 / p
    k = (MLab.prod(-p)/MLab.prod(-z)).real
    return z, p, k
    

EPSILON = 2e-16

def vratio(u, ineps, mp):
    [s,c,d,phi] = special.ellpj(u,mp)
    ret = abs(ineps - s/c)
    try:
        return ret[0]
    except IndexError, TypeError:
        return ret

def kratio(m, k_ratio):
    if m < 0:
        m = 0.0
    if m > 1:
        m = 1.0
    if abs(m) > EPSILON and (abs(m) + EPSILON) < 1:
        k = special.ellpk([1-m,m])
        r = k[0] / k[1] - k_ratio
    elif abs(m) > EPSILON:
        r = -k_ratio
    else:
        r = 1e20
    try:
        return abs(r)[0]
    except TypeError:
        return abs(r)

def ellipap(N,rp,rs):
    """Return (z,p,k) zeros, poles, and gain of an Nth order normalized
    prototype elliptic analog lowpass filter with rp decibels of ripple
    in the passband and a stopband rs decibels down.  Broken...
    """
    if N == 1:
        p = -sqrt(1.0/(10**(0.1*rp)-1.0))
        k = -p
        z = []
        return z, p, k

    eps = Num.sqrt(10**(0.1*rp)-1)
    ck1 = eps / Num.sqrt(10**(0.1*rs)-1)
    ck1p = Num.sqrt(1-ck1*ck1)
    if ck1p == 1:
        raise ValueError, "Cannot design a filter with given rp and rs specifications."

    wp = 1
    val = special.ellpk([1-ck1*ck1,1-ck1p*ck1p])
    if abs(1-ck1p*ck1p) < EPSILON:
        krat = 0
    else:
        krat = N*val[0] / val[1]

    m = optimize.fmin(kratio, 0.5, args=(krat,), maxfun=250, maxiter=250,
                      disp=0)
    if m < 0 or m > 1:
        m = optimize.fminbound(kratio, 0, 1, args=(krat,), maxfun=250,
                               maxiter=250, disp=0)
    
    capk = special.ellpk(1-m)
    ws = wp / sqrt(m)
    m1 = 1-m

    j = arange(1-N%2,N,2)
    jj = len(j)

    [s,c,d,phi] = special.ellpj(j*capk/N,m*ones(jj))
    snew = Num.compress(abs(s) > EPSILON, s)
    z = 1.0 / (sqrt(m)*snew)
    z = 1j*z
    z = Num.concatenate((z,conjugate(z)))

    r = optimize.fmin(vratio, special.ellpk(1-m), args=(1/eps, ck1p*ck1p),
                      maxfun=250, maxiter=250, printmessg=0)
    v0 = capk * r / (N*capck1[0])

    [sv,cv,dv,phi] = special.ellpj(v0,1-m)
    p = -(c*d*sv*cv + 1j*s*dv) / (1-(d*sv)**2.0)

    if N % 2:
        newp = Num.compress(abs(p.imag) > EPSILON*Num.sqrt(MLab.sum(p*Num.conjugate(p)).real), p)
        p = Num.concatenate((p,conjugate(newp)))
    else:
        p = Num.concatenate((p,conjugate(p)))

    k = (MLab.prod(-p) / MLab.prod(-z)).real
    if N % 2 == 0:
        k = k / Num.sqrt((1+eps*eps))

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
    else:
        raise ValueError, "Bessel Filter not supported for given order."

    return z, p, k
