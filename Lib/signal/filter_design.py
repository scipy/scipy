import Numeric
from scipy.fastumath import *
Num = Numeric
abs = absolute
pi = Numeric.pi
import scipy
from scipy import asarray_1d


def butter(N, Wn, bandtype='band', analog=0, output='ba'):
    """Butterworth digital and analog filter design.

    Description:

      Design an Nth order lowpass digital Butterworth filter and return the
      filter coeeficients in (B,A), zero-pole, or state-space form.

    Inputs:

      
    """
    pass

def cheby1():
    pass

def cheby2():
    pass

def ellip():
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

    wp = asarray_1d(wp)
    ws = asarray_1d(ws)
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
        W0 = nat / ( ( 10**(0.1*abs(gstop))-1)**(1.0/(2.0*or)))
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
    wp = asarray_1d(wp)
    ws = asarray_1d(ws)
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
    wp = asarray_1d(wp)
    ws = asarray_1d(ws)
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
    wp = asarray_1d(wp)
    ws = asarray_1d(ws)
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
    

def besselap():
    pass

def buttap():
    pass

def cheb1ap():
    pass

def cheb2ap():
    pass

def ellipap():
    pass

def besself():
    pass



