import Numeric
Num = Numeric
pi = Numeric.pi
import scipy.special

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


def band_stop_cost(wp, ind, WP, WS, rs, rp, type):
    """Band Stop Objective Function for order minimization

    Description:

      Returns the non-integer order for an analog band stop filter.

    Inputs:

      wp -- passband edge
      ind -- index specifying which passband edge to vary (0 or 1).
      WP -- two element vector of fixed passband edges.
      WS -- two element vector of fixed stopband edges.
      rs -- amount in dB of attenuation in stopband.
      rp -- amount in dB of ripple in the passband.
      type -- 'butter', 'cheby', or 'ellip':

    Outputs: (n,)

      n -- filter order (non-integer)
    """

    WP[ind] = wp
    WSflip = array([WS(0), -WS(1)])
    WA = WS*(WP(0)-WP(1)) / (WS**2 - WP(0)*WP(1))
    WA = min(abs(WA))

    if type == 'butter':
        n = (Num.log10 ( (10**(0.1*abs(rs)) - 1.0) / \
                         (10**(0.1*abs(rp))-1)) / (2*Num.log10(WA)))
    elif type == 'cheby':
        n = Num.acosh(Num.sqrt((10**(0.1*abs(rs))-1) / \
                               (10**(0.1*abs(rp))-1))) / acosh(WA)
    elif type == 'ellip':
        epsilon = sqrt(10**(0.1*rp)-1)
        k1 = epsilon / sqrt(10**(0.1*rs)-1)
        k = 1.0 / WA
        capk = scipy.special.ellpk([k**2, 1-k**2])
        capk1 = scipy.special.ellpk([k1**2, 1-k1**2])
        n = (capk[0]*capk1[1] / (capk[1]*capk1[0]))
    else:
        raise ValueError, "Incorrect type: ", type

def buttord(wp, ws, rp, rs, analog=0):
    """Butterworth filter order selection.

    """

    filter_type = 2*(len(wp)-1)
    if wp(0) < ws(0):
        filter_type += 1
    else:
        filter_type += 2

    # Pre-warp frequencies
    if not analog:
        WP = tan(pi*wp/2)
        WS = tan(pi*ws/2)
    else:
        WP = wp
        WS = ws

    if ftype == 1:           # low
        WA = WS / WP
    elif ftype == 2:          # high
        WA = WP / WS
    elif ftype == 3:     # stop
        wp0 = scipy.optimize.fminbound(band_stop_cost, WP[0], WS[1]-1e-12,
                                       args=(0,WP,WS,rs,rp,'butter'), disp=0)
        WP[0] = wp0
        wp1 = scipy.optimize.fminbound(band_stop_cost, WS[1]+1e-12, WP[1],
                                       args=(1,WP,WS,rs,rp,'butter'), disp=0)
        WP[1] = wp1
        WA = (WS * (WP[0]-WP[1])) / (WS**2 - WP[0]*WP[1])
    elif ftype == 4:  # pass
        WA = (WS**2 - WP[0]*WP[1]) / (WS* (WP[0]-WP[1]))

    WA = min(abs(WA))

    order = Num.ceil( Num.log10( (10**(0.1*abs(rs)) - 1) / \
                                 (10**(0.1*abs(rp)) - 1) ) / \
                      (2*Num.log10(WA)))

    # Find the butterworth natural frequency Wo (or the "3dB" frequency")
    # to give exactly rs at WA. W0 will be between 1 and WA
    W0 = WA / ( ( 10**(0.1*abs(rs))-1)**(1.0/(2*order)))

    # now convert this frequency back from lowpass prototype
    # to the original analog filter

    if ftype == 1:  # low
        WN = WO*WP
    elif ftype == 2: # high
        WN = WP / W0
    elif ftype == 3:  # stop
        WN = zeros(2,Numeric.Float)
        WN[0] = ((WP[1] - WP[0]) + sqrt((WP[1] - WP[0])**2 + \
                                        4*W0**2 * WP[0] * WP[1])) / (2*W0)
        WN[1] = ((WP[1] - WP[0]) - sqrt((WP[1] - WP[0])**2 + \
                                        4*W0**2 * WP[0] * WP[1])) / (2*W0)
        WN = Num.sort(Num.abs(WN))
    elif ftype == 4: # pass
        W0 = array([-W0, W0],Numeric.Float)
        WN = -W0 * (WP[1]-WP[0]) / 2.0 + sqrt(W0**2 / 4.0 * (WP[1]-WP[0])**2 + \
                                              WP[0]*WP[1])
        WN = Num.sort(Num.abs(WN))
    else:
        raise ValueError, "Bad type."

    if not analog:
        wn = (2.0/pi)*Num.atan(WN)
    else:
        wn = WN

    return order, wn

def cheb1ord():
    pass

def cheb2ord():
    pass

def ellipord():
    pass

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



