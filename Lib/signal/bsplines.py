import scipy.special
from scipy_base import *
from scipy_base.fastumath import sqrt, exp, greater, equal, cos, add, sin
from spline import *      # C-modules

gamma = scipy.special.gamma
def factorial(n):
    return gamma(n+1)

def spline_filter(Iin, lmbda=5.0):
    """Smoothing spline (cubic) filtering of a rank-2 array.

    Filter an input data set, Iin, using a (cubic) smoothing spline of
    fall-off lmbda.
    """
    intype = Iin.typecode()
    hcol = sarray([1.0,4.0,1.0],'f')/6.0
    if intype in ['F','D']:
        Iin = Iin.astype('F')
        ckr = cspline2d(Iin.real,lmbda)
        cki = cspline2d(Iin.imag,lmbda)
        outr = sepfir2d(ckr,hcol,hcol)
        outi = sepfir2d(cki,hcol,hcol)
        out = (outr + 1j*outi).astype(intype)
    elif intype in ['f','d']:
        ckr = cspline2d(Iin,lmbda)
        out = sepfir2d(ckr, hcol, hcol)
        out = out.astype(intype)
    else:
        raise TypeError;
    return out            

def _bspline(x,n):
    """bspline(x,n) -> y:  B-spline basis function of order n.
    """
    jlist = arange(n+2)
    val = 0.0
    baseval = x + (n+1)/2.0
    for j in jlist:
        xval = baseval - j
        if xval >=0:
            if j % 2 == 0:
                fact = 1
            else:
                fact = -1                
            term = fact * (n+1) * xval**n / factorial(j) / factorial(n+1-j)
            val = val + term
    return val

bspline = vectorize(_bspline)

def gauss_spline(x,n):
    """Gaussian approximation to B-spline basis function of order n.
    """
    signsq = (n+1) / 12.0
    return 1/sqrt(2*pi*signsq) * exp(-x**2 / 2 / signsq)


def cubic(x):
    ax = abs(x)
    f1 = 2.0/3 - 1.0/2*ax**2 * (2-ax)
    f2 = 1.0/6*(2-ax)**3
    f3 = where(greater(ax,1),f2,f1)
    return where(greater(ax,2),0,f3)

def quintic(n):
    n = asarray(n)
    an = abs(n)
    f1 = where(equal(an,0),66.0/120.0,0)
    f2 = where(equal(an,1),26.0/120.0,0)
    f3 = where(equal(an,2),1.0/120.0,0)
    return f1 + f2 + f3

def quadratic(x):
    ax = abs(x)
    f1 = 0.75-ax**2
    f2 = (ax-1.5)**2 / 2.0
    f3 = where(greater(abs(x),0.5),f2,f1)
    return where(greater(abs(x),1.5),0,f3)



def c0_P(order):
    # values taken from Unser, et.al. 1993 IEEE
    if order == 0:
        c0 = 1
        P = array([1])
    elif order == 1:
        c0 = 1
        P = array([0,1])
    elif order == 2:
        c0 = 8
        P = array([1,6,1])
    elif order == 3:
        c0 = 6
        P = array([1,4,1])
    elif order == 4:
        c0 = 384
        P = array([1,76,230,76,1])
    elif order == 5:
        c0 = 120
        P = array([1,26,66,26,1])
    elif order == 6:
        c0 = 46080
        P = array([1,722,10543,23548, 10543, 722, 1])
    elif order == 7:
        c0 = 5040
        P = array([1,120,1191,2416,1191, 120, 1])
    else:
        raise ValueError, "Unknown order."

def _coeff_smooth(lam):
    xi = 1 - 96*lam + 24*lam * sqrt(3 + 144*lam)
    omeg = arctan2(sqrt(144*lam-1),sqrt(xi))
    rho = (24*lam - 1 - sqrt(xi)) / (24*lam)
    rho = rho * sqrt((48*lam + 24*lam * sqrt(3+144*lam))/xi)
    return rho,omeg


def _cubic_smooth_coeff(signal,lamb):
    rho, omega = _coeff_smooth(lamb)
    cs = 1-2*rho*cos(omega) + rho*rho
    K = len(signal)
    yp = zeros((K,),signal.typecode())
    k = arange(K)
    yp[0] = hc(0,cs,rho,omega)*signal[0] + \
            add.reduce(hc(k+1,cs,rho,omega)*signal)

    yp[1] = hc(0,cs,rho,omega)*signal[0] + \
            hc(1,cs,rho,omega)*signal[1] + \
            add.reduce(hc(k+2,cs,rho,omega)*signal)

    for n in range(2,K):
        yp[n] = cs * signal[n] + 2*rho*cos(omega)*yp[n-1] - rho*rho*yp[n-2]
        
    y = zeros((K,),signal.typecode())

    y[K-1] = add.reduce((hs(k,cs,rho,omega) + hs(k+1,cs,rho,omega))*signal[::-1])
    y[K-2] = add.reduce((hs(k-1,cs,rho,omega) + hs(k+2,cs,rho,omega))*signal[::-1])

    for n in range(K-3,-1,-1):
        y[n] = cs*yp[n] + 2*rho*cos(omega)*y[n+1] - rho*rho*y[n+2]

    return y

def _cubic_coeff(signal):
    zi = -2 + sqrt(3)
    K = len(signal)
    yplus = zeros((K,),signal.typecode())
    powers = zi**arange(K)
    yplus[0] = signal[0] + zi*add.reduce(powers*signal)
    for k in range(1,K):
        yplus[k] = signal[k] + zi*yplus[k-1]
    output = zeros((K,),signal.typecode())
    output[K-1] = zi / (zi-1)*yplus[K-1]
    for k in range(K-2,-1,-1):
        output[k] = zi*(output[k+1]-yplus[k])
    return output*6.0

def _quadratic_coeff(signal):
    zi = -3 + 2*sqrt(2.0)    
    K = len(signal)
    yplus = zeros((K,),signal.typecode())
    powers = zi**arange(K)
    yplus[0] = signal[0] + zi*add.reduce(powers*signal)    
    for k in range(1,K):
        yplus[k] = signal[k] + zi*yplus[k-1]
    output = zeros((K,),signal.typecode())
    output[K-1] = zi / (zi-1)*yplus[K-1]
    for k in range(K-2,-1,-1):
        output[k] = zi*(output[k+1]-yplus[k])
    return output*8.0

def cspline1d(signal,lamb=0.0):
    """Compute cubic spline coefficients for rank-1 array.

    Description:

      Find the cubic spline coefficients for a 1-D signal assuming
      mirror-symmetric boundary conditions.   To obtain the signal back from
      the spline representation mirror-symmetric-convolve these coefficients
      with a length 3 FIR window [1.0, 4.0, 1.0]/ 6.0 .

    Inputs:

      signal -- a rank-1 array representing samples of a signal.
      lamb -- smoothing coefficient (default = 0.0)

    Output:

      c -- cubic spline coefficients.
    """
    if lamb != 0.0:
        return _cubic_smooth_coeff(signal,lamb)
    else:
        return _cubic_coeff(signal)


def qspline1d(signal,lamb=0.0):
    """Compute quadratic spline coefficients for rank-1 array.

    Description:

      Find the quadratic spline coefficients for a 1-D signal assuming
      mirror-symmetric boundary conditions.   To obtain the signal back from
      the spline representation mirror-symmetric-convolve these coefficients
      with a length 3 FIR window [1.0, 6.0, 1.0]/ 8.0 .

    Inputs:

      signal -- a rank-1 array representing samples of a signal.
      lamb -- smoothing coefficient (must be zero for now.)

    Output:

      c -- cubic spline coefficients.
    """
    if lamb != 0.0:
        raise ValueError, "Smoothing quadratic splines not supported yet."
    else:
        return _quadratic_coeff(signal)

    
def hc(k,cs,rho,omega):
    return cs / sin(omega) * (rho**k)*sin(omega*(k+1))*(greater(k,-1))

def hs(k,cs,rho,omega):
    c0 = cs*cs * (1 + rho*rho) / (1 - rho*rho) / (1-2*rho*rho*cos(2*omega) + rho**4)
    gamma = (1-rho*rho) / (1+rho*rho) / tan(omega)
    ak = abs(k)
    return c0 * rho**ak * (cos(omega*ak) + gamma*sin(omega*ak))
    












