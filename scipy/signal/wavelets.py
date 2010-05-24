__all__ = ['daub','qmf','cascade','morlet']

import numpy as np
from numpy.dual import eig
from scipy.misc import comb
from scipy import linspace, pi, exp

def daub(p):
    """
    The coefficients for the FIR low-pass filter producing Daubechies wavelets.

    p>=1 gives the order of the zero at f=1/2.
    There are 2p filter coefficients.

    Parameters
    ----------
    p : int
        Order of the zero at f=1/2, can have values from 1 to 34.

    """
    sqrt = np.sqrt
    assert(p>=1)
    if p==1:
        c = 1/sqrt(2)
        return np.array([c,c])
    elif p==2:
        f = sqrt(2)/8
        c = sqrt(3)
        return f*np.array([1+c,3+c,3-c,1-c])
    elif p==3:
        tmp  = 12*sqrt(10)
        z1 = 1.5 + sqrt(15+tmp)/6 - 1j*(sqrt(15)+sqrt(tmp-15))/6
        z1c = np.conj(z1)
        f = sqrt(2)/8
        d0 = np.real((1-z1)*(1-z1c))
        a0 = np.real(z1*z1c)
        a1 = 2*np.real(z1)
        return f/d0*np.array([a0, 3*a0-a1, 3*a0-3*a1+1, a0-3*a1+3, 3-a1, 1])
    elif p<35:
        # construct polynomial and factor it
        if p<35:
            P = [comb(p-1+k,k,exact=1) for k in range(p)][::-1]
            yj = np.roots(P)
        else:  # try different polynomial --- needs work
            P = [comb(p-1+k,k,exact=1)/4.0**k for k in range(p)][::-1]
            yj = np.roots(P) / 4
        # for each root, compute two z roots, select the one with |z|>1
        # Build up final polynomial
        c = np.poly1d([1,1])**p
        q = np.poly1d([1])
        for k in range(p-1):
            yval = yj[k]
            part = 2*sqrt(yval*(yval-1))
            const = 1-2*yval
            z1 = const + part
            if (abs(z1)) < 1:
                z1 = const - part
            q = q * [1,-z1]

        q = c * np.real(q)
        # Normalize result
        q = q / np.sum(q) * sqrt(2)
        return q.c[::-1]
    else:
        raise ValueError, "Polynomial factorization does not work "\
              "well for p too large."

def qmf(hk):
    """Return high-pass qmf filter from low-pass
    """
    N = len(hk)-1
    asgn = [{0:1,1:-1}[k%2] for k in range(N+1)]
    return hk[::-1]*np.array(asgn)

def wavedec(amn,hk):
    gk = qmf(hk)
    return NotImplemented

def cascade(hk,J=7):
    """(x,phi,psi) at dyadic points K/2**J from filter coefficients.

    Inputs:
      hk  -- coefficients of low-pass filter
      J   -- values will be computed at grid points $K/2^J$

    Outputs:
      x   -- the dyadic points $K/2^J$ for $K=0...N*(2^J)-1$
              where len(hk)=len(gk)=N+1
      phi -- the scaling function phi(x) at x
               $\phi(x) = \sum_{k=0}^{N} h_k \phi(2x-k)$
      psi -- the wavelet function psi(x) at x
               $\psi(x) = \sum_{k=0}^N g_k \phi(2x-k)$
             Only returned if gk is not None

    Algorithm:

      Uses the vector cascade algorithm described by Strang and Nguyen in
      "Wavelets and Filter Banks"

      Builds a dictionary of values and slices for quick reuse.
      Then inserts vectors into final vector at then end

    """

    N = len(hk)-1

    if (J > 30 - np.log2(N+1)):
        raise ValueError, "Too many levels."
    if (J < 1):
        raise ValueError, "Too few levels."


    # construct matrices needed
    nn,kk = np.ogrid[:N,:N]
    s2 = np.sqrt(2)
    # append a zero so that take works
    thk = np.r_[hk,0]
    gk = qmf(hk)
    tgk = np.r_[gk,0]

    indx1 = np.clip(2*nn-kk,-1,N+1)
    indx2 = np.clip(2*nn-kk+1,-1,N+1)
    m = np.zeros((2,2,N,N),'d')
    m[0,0] = np.take(thk,indx1,0)
    m[0,1] = np.take(thk,indx2,0)
    m[1,0] = np.take(tgk,indx1,0)
    m[1,1] = np.take(tgk,indx2,0)
    m *= s2

    # construct the grid of points
    x = np.arange(0,N*(1<<J),dtype=np.float) / (1<<J)
    phi = 0*x

    psi = 0*x

    # find phi0, and phi1
    lam, v = eig(m[0,0])
    ind = np.argmin(np.absolute(lam-1))
    # a dictionary with a binary representation of the
    #   evaluation points x < 1 -- i.e. position is 0.xxxx
    v = np.real(v[:,ind])
    # need scaling function to integrate to 1 so find
    #  eigenvector normalized to sum(v,axis=0)=1
    sm = np.sum(v)
    if sm < 0:  # need scaling function to integrate to 1
        v = -v
        sm = -sm
    bitdic = {}
    bitdic['0'] = v / sm
    bitdic['1'] = np.dot(m[0,1],bitdic['0'])
    step = 1<<J
    phi[::step] = bitdic['0']
    phi[(1<<(J-1))::step] = bitdic['1']
    psi[::step] = np.dot(m[1,0],bitdic['0'])
    psi[(1<<(J-1))::step] = np.dot(m[1,1],bitdic['0'])
    # descend down the levels inserting more and more values
    #  into bitdic -- store the values in the correct location once we
    #  have computed them -- stored in the dictionary
    #  for quicker use later.
    prevkeys = ['1']
    for level in range(2,J+1):
        newkeys = ['%d%s' % (xx,yy) for xx in [0,1] for yy in prevkeys]
        fac = 1<<(J-level)
        for key in newkeys:
            # convert key to number
            num = 0
            for pos in range(level):
                if key[pos] == '1':
                    num += (1<<(level-1-pos))
            pastphi = bitdic[key[1:]]
            ii = int(key[0])
            temp = np.dot(m[0,ii],pastphi)
            bitdic[key] = temp
            phi[num*fac::step] = temp
            psi[num*fac::step] = np.dot(m[1,ii],pastphi)
        prevkeys = newkeys

    return x, phi, psi

def morlet(M, w=5.0, s=1.0, complete=True):
    """
    Complex Morlet wavelet.

    Parameters
    ----------
    M : int
        Length of the wavelet.
    w : float
        Omega0
    s : float
        Scaling factor, windowed from -s*2*pi to +s*2*pi.
    complete : bool
        Whether to use the complete or the standard version.

    Notes
    -----
    The standard version:
        pi**-0.25 * exp(1j*w*x) * exp(-0.5*(x**2))

        This commonly used wavelet is often referred to simply as the
        Morlet wavelet.  Note that, this simplified version can cause
        admissibility problems at low values of w.

    The complete version:
        pi**-0.25 * (exp(1j*w*x) - exp(-0.5*(w**2))) * exp(-0.5*(x**2))

        The complete version of the Morlet wavelet, with a correction
        term to improve admissibility. For w greater than 5, the
        correction term is negligible.

    Note that the energy of the return wavelet is not normalised
    according to s.

    The fundamental frequency of this wavelet in Hz is given
    by f = 2*s*w*r / M where r is the sampling rate.

    """
    x = linspace(-s*2*pi,s*2*pi,M)
    output = exp(1j*w*x)

    if complete:
        output -= exp(-0.5*(w**2))

    output *= exp(-0.5*(x**2)) * pi**(-0.25)

    return output
