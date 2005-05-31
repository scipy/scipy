
import scipy_base as sb
import scipy as s

def daub(p):
    """The coefficients for the FIR low-pass filter producing Daubechies wavelets.

    p>=1 gives the order of the zero at f=1/2.  There are 2p filter coefficients.
    """
    sqrt = sb.sqrt
    assert(p>=1)
    if p==1:
        c = 1/sqrt(2)
        return sb.array([c,c])
    elif p==2:
        f = sqrt(2)/8
        c = sqrt(3)
        return f*sb.array([1+c,3+c,3-c,1-c])
    elif p==3:
        tmp  = 12*sqrt(10)
        z1 = 1.5 + sqrt(15+tmp)/6 - 1j*(sqrt(15)+sqrt(tmp-15))/6
        z1c = sb.conj(z1)
        f = sqrt(2)/8        
        d0 = sb.real((1-z1)*(1-z1c))
        a0 = sb.real(z1*z1c)
        a1 = 2*sb.real(z1)
        return f/d0*sb.array([a0, 3*a0-a1, 3*a0-3*a1+1, a0-3*a1+3, 3-a1, 1])
    elif p<81:
        # construct polynomial and factor it
        if p<35:
            P = [s.comb(p-1+k,k,exact=1) for k in range(p)][::-1]
            yj = sb.roots(P)
        else:
            raise ValueError, "Cannot factor such large polynomial well."
            k = sb.r_[0:p]
            P = s.comb(p-1+k,k)/4.0**k
            yj = sb.roots(P) / 4
        # for each root, compute two z roots, select the one with |z|>1
        # Build up final polynomial
        c = sb.poly1d([1,1])**p
        q = sb.poly1d([1])
        for k in range(p-1):
            yval = yj[k]
            part = 2*sqrt(yval*(yval-1))
            const = 1-2*yval
            z1 = const + part
            if (abs(z1)) < 1:
                z1 = const - part
            q = q * [1,-z1]
        q = sb.real(q) * c
        # Normalize result
        q = q / sb.sum(q) * sqrt(2)
        return q.c[::-1]        
    else:
        raise ValueError, "Polynomial factorization does not work "\
              "well for p too large."

def qmf(hk):
    """Return high-pass qmf filter from low-pass
    """
    N = len(hk)-1
    asgn = [{0:1,1:-1}[k%2] for k in range(N+1)]
    return hk[::-1]*sb.array(asgn)

def wavedec(amn,hk):
    gk = qmf(hk)
    return NotImplemented

def cascade(hk,J=7):
    """(x,phi,psi) at dyadic points K/2**J from filter coefficients.

    Inputs:
      hk  -- coefficients of low-pass filter
      J   -- values will be computed at grid points $K/2^J$

    Outputs:
      x   -- the dyadic points $K/2^J$ for $K=0...N*2^J-1$
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

    if (J > 30 - sb.log2(N+1)):
        raise ValueError, "Too many levels."
    if (J < 1):
        raise ValueError, "Too few levels."
    

    # construct matrices needed
    nn,kk = sb.ogrid[:N,:N]
    s2 = sqrt(2)    
    # append a zero so that take works
    thk = sb.r_[hk,0]
    gk = qmf(hk)
    tgk = sb.r_[gk,0]    
    
    indx1 = clip(2*nn-kk,-1,N+1)
    indx2 = clip(2*nn-kk+1,-1,N+1)
    m = zeros((2,2,N,N),'d')
    m[0,0] = take(thk,indx1)
    m[0,1] = take(thk,indx2)
    m[1,0] = take(tgk,indx1)
    m[1,1] = take(tgk,indx2)        
    m *= s2

    # construct the grid of points
    x = arange(0,N*(1<<J),typecode=sb.Float) / (1<<J)
    phi = 0*x

    psi = 0*x

    # find phi0, and phi1
    lam, v = s.linalg.eig(m[0,0])
    ind = argmin(abs(lam-1))
    # a dictionary with a binary representation of the
    #   evaluation points x < 1 -- i.e. position is 0.xxxx
    v = v[:,ind]
    # need scaling function to integrate to 1 so find
    #  eigenvector normalized to sum(v)=1
    sm = sum(v)
    if sm < 0:  # need scaling function to integrate to 1 
        v = -v
        sm = -sm
    bitdic = {}    
    bitdic['0'] = v / sm
    bitdic['1'] = dot(m[0,1],bitdic['0'])
    step = 1<<J
    phi[::step] = bitdic['0']
    phi[(1<<(J-1))::step] = bitdic['1']
    psi[::step] = dot(m[1,0],bitdic['0'])
    psi[(1<<(J-1))::step] = dot(m[1,1],bitdic['0'])
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
            temp = dot(m[0,ii],pastphi)
            bitdic[key] = temp
            phi[num*fac::step] = temp
            psi[num*fac::step] = dot(m[1,ii],pastphi)
        prevkeys = newkeys

    return x, phi, psi

    

        
    

    
    
        
