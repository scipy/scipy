import sigtools
import scipy.special as special
from scipy import iscomplex, fft, ifft
from scipy import polyadd, polymul, polydiv, polysub, \
                  roots, poly, polyval, polyder
import scipy
import Numeric
from Numeric import array, asarray

_modedict = {'valid':0, 'same':1, 'full':2}
_boundarydict = {'fill':0, 'pad':0, 'wrap':2, 'circular':2, 'symm':1, 'symmetric':1, 'reflect':4}
                                                                            
def _valfrommode(mode):
    try:
        val = _modedict[mode]
    except KeyError:
        if mode not in [0,1,2]:
            raise ValueError, "Acceptable mode flags are 'valid' (0), 'same' (1), or 'full' (2)."
        val = mode
    return val

def _bvalfromboundary(boundary):
    try:
        val = _boundarydict[boundary] << 2
    except KeyError:
        if val not in [0,1,2] :
            raise ValueError, "Acceptable boundary flags are 'fill', 'wrap' (or 'circular'), \n  and 'symm' (or 'symmetric')."
        val = boundary << 2
    return val


def correlate(in1, in2, mode='full'):
    """Cross-correlate two N-dimensional arrays.

  Description:

     Cross-correlate in1 and in2 with the output size determined by mode.

  Inputs:

    in1 -- an N-dimensional array.
    in2 -- an array with the same number of dimensions as in1.
    mode -- a flag indicating the size of the output
            'valid'  (0): The output consists only of those elements that
                            do not rely on the zero-padding.
            'same'   (1): The output is the same size as the input centered
                            with respect to the 'full' output.
            'full'   (2): The output is the full discrete linear
                            cross-correlation of the inputs. (Default)

  Outputs:  (out,)

    out -- an N-dimensional array containing a subset of the discrete linear
           cross-correlation of in1 with in2.
 
    """
    # Code is faster if kernel is smallest array.
    volume = Numeric.asarray(in1)
    kernel = Numeric.asarray(in2)
    if (Numeric.product(kernel.shape) > Numeric.product(volume.shape)):
        temp = kernel
        kernel = volume
        volume = temp
        del temp

    val = _valfrommode(mode)

    return sigtools._correlateND(volume, kernel, val)

def convolve(in1, in2, mode='full'):
    """Convolve two N-dimensional arrays.

  Description:

     Convolve in1 and in2 with output size determined by mode.

  Inputs:

    in1 -- an N-dimensional array.
    in2 -- an array with the same number of dimensions as in1.
    mode -- a flag indicating the size of the output
            'valid'  (0): The output consists only of those elements that
                            do not rely on the zero-padding.
            'same'   (1): The output is the same size as the input centered
                            with respect to the 'full' output.
            'full'   (2): The output is the full discrete linear convolution
                            of the inputs. (Default)

  Outputs:  (out,)

    out -- an N-dimensional array containing a subset of the discrete linear
           convolution of in1 with in2.

    """
    volume = Numeric.asarray(in1)
    kernel = Numeric.asarray(in2)
    if (Numeric.product(kernel.shape) > Numeric.product(volume.shape)):
        temp = kernel
        kernel = volume
        volume = temp
        del temp

    slice_obj = [slice(None,None,-1)]*len(kernel.shape)
    val = _valfrommode(mode)
    
    return sigtools._correlateND(volume,kernel[slice_obj],val)

def order_filter(a, domain, order):
    """Perform an order filter on an N-dimensional array.
    
  Description:

    Perform an order filter on the array in.  The domain argument acts as a
    mask centered over each pixel.  The non-zero elements of domain are
    used to select elements surrounding each input pixel which are placed
    in a list.   The list is sorted, and the output for that pixel is the
    element corresponding to rank in the sorted list.
    
  Inputs:

    in -- an N-dimensional input array.
    domain -- a mask array with the same number of dimensions as in.  Each
              dimension should have an odd number of elements.
    rank -- an non-negative integer which selects the element from the
            sorted list (0 corresponds to the largest element, 1 is the
            next largest element, etc.)

  Output: (out,)

    out -- the results of the order filter in an array with the same
           shape as in.
          
    """
    domain = Numeric.asarray(domain)
    size = domain.shape
    for k in range(len(size)):
        if (size[k] % 2) != 1:
            raise ValueError, "Each dimension of domain argument should have an odd number of elements."
    return sigtools._orderfilterND(a, domain, rank)
   

def medfilt(volume,kernel_size=None):
    """Perform a median filter on an N-dimensional array.

  Description:

    Apply a median filter to the input array using a local window-size
    given by kernel_size.

  Inputs:

    in -- An N-dimensional input array.
    kernel_size -- A scalar or an N-length list giving the size of the
                   median filter window in each dimension.  Elements of
                   kernel_size should be odd.  If kernel_size is a scalar,
                   then this scalar is used as the size in each dimension.

  Outputs: (out,)

    out -- An array the same size as input containing the median filtered
           result.
  
    """
    volume = Numeric.asarray(volume)
    if kernel_size == None:
        kernel_size = [3] * len(volume.shape)
    kernel_size = Numeric.asarray(kernel_size)
    if len(kernel_size.shape) == 0:
        kernel_size = [kernel_size.toscalar()] * len(volume.shape)
    kernel_size = Numeric.asarray(kernel_size)

    for k in range(len(volume.shape)):
        if (kernel_size[k] % 2) != 1:
            raise ValueError, "Each element of kernel_size should be odd." 

    domain = Numeric.ones(kernel_size)

    numels = Numeric.product(kernel_size)
    order = numels/2
    return sigtools._order_filterND(volume,domain,order)


def wiener(im,mysize=None,noise=None):
    """Perform a wiener filter on an N-dimensional array.

  Description:

    Apply a wiener filter to the N-dimensional array in.

  Inputs:

    in -- an N-dimensional array.
    kernel_size -- A scalar or an N-length list giving the size of the
                   median filter window in each dimension.  Elements of
                   kernel_size should be odd.  If kernel_size is a scalar,
                   then this scalar is used as the size in each dimension.
    noise -- The noise-power to use.  If None, then noise is estimated as
             the average of the local variance of the input.

  Outputs: (out,)

    out -- Wiener filtered result with the same shape as in.

    """
    im = Numeric.asarray(im)
    if mysize == None:
        mysize = [3] * len(im.shape)
    mysize = Numeric.asarray(mysize);

    # Estimate the local mean
    lMean = correlate(im,Numeric.ones(mysize),1) / Numeric.product(mysize)

    # Estimate the local variance
    lVar = correlate(im**2,Numeric.ones(mysize),1) / Numeric.product(mysize) - lMean**2

    # Estimate the noise power if needed.
    if noise==None:
        noise = Numeric.mean(Numeric.ravel(lVar))

    # Compute result
    # f = lMean + (maximum(0, lVar - noise) ./
    #               maximum(lVar, noise)) * (im - lMean) 
    #
    out = im - lMean
    im = lVar - noise
    im = Numeric.maximum(im,0)
    lVar = Numeric.maximum(lVar,noise)
    out = out / lVar
    out = out * im
    out = out + lMean

    return out

def convolve2d(in1, in2, mode='full', boundary='fill', fillvalue=0):
    """Conolve two 2-dimensional arrays.

  Description:

     Convolve in1 and in2 with output size determined by mode and boundary
     conditions determined by boundary and fillvalue.

  Inputs:

    in1 -- a 2-dimensional array.
    in2 -- a 2-dimensional array.
    mode -- a flag indicating the size of the output
            'valid'  (0): The output consists only of those elements that
                            do not rely on the zero-padding.
            'same'   (1): The output is the same size as the input centered
                            with respect to the 'full' output.
            'full'   (2): The output is the full discrete linear convolution
                            of the inputs. (*Default*)
    boundary -- a flag indicating how to handle boundaries
                'fill' : pad input arrays with fillvalue. (*Default*)
                'wrap' : circular boundary conditions.
                'symm' : symmetrical boundary conditions.
    fillvalue -- value to fill pad input arrays with (*Default* = 0)

  Outputs:  (out,)

    out -- a 2-dimensional array containing a subset of the discrete linear
           convolution of in1 with in2.

    """
    val = _valfrommode(mode)
    bval = _bvalfromboundary(boundary)
        
    return sigtools._convolve2d(in1,in2,1,val,bval,fillvalue)

def correlate2d(in1, in2, mode='full', boundary='fill', fillvalue=0):
    """Cross-correlate two 2-dimensional arrays.

  Description:

     Cross correlate in1 and in2 with output size determined by mode
     and boundary conditions determined by boundary and fillvalue.

  Inputs:

    in1 -- a 2-dimensional array.
    in2 -- a 2-dimensional array.
    mode -- a flag indicating the size of the output
            'valid'  (0): The output consists only of those elements that
                            do not rely on the zero-padding.
            'same'   (1): The output is the same size as the input centered
                            with respect to the 'full' output.
            'full'   (2): The output is the full discrete linear convolution
                            of the inputs. (*Default*)
    boundary -- a flag indicating how to handle boundaries
                'fill' : pad input arrays with fillvalue. (*Default*)
                'wrap' : circular boundary conditions.
                'symm' : symmetrical boundary conditions.
    fillvalue -- value to fill pad input arrays with (*Default* = 0)

  Outputs:  (out,)

    out -- a 2-dimensional array containing a subset of the discrete linear
           cross-correlation of in1 with in2.

    """
    val = _valfrommode(mode)
    bval = _bvalfromboundary(boundary)
        
    return sigtools._convolve2d(in1, in2, 0,val,bval,fillvalue)

def medfilt2d(input, kernel_size=3):
    """Median filter two 2-dimensional arrays.

  Description:

    Apply a median filter to the input array using a local window-size
    given by kernel_size (must be odd).

  Inputs:

    in -- An 2 dimensional input array.
    kernel_size -- A scalar or an length-2 list giving the size of the
                   median filter window in each dimension.  Elements of
                   kernel_size should be odd.  If kernel_size is a scalar,
                   then this scalar is used as the size in each dimension.

  Outputs: (out,)

    out -- An array the same size as input containing the median filtered
           result.
    """
    image = Numeric.asarray(input)
    if kernel_size == None:
        kernel_size = [3] * 2
    kernel_size = Numeric.asarray(kernel_size)
    if len(kernel_size.shape) == 0:
        kernel_size = [kernel_size.toscalar()] * 2
    kernel_size = Numeric.asarray(kernel_size)

    for size in kernel_size:
        if (size % 2) != 1:
            raise ValueError, "Each element of kernel_size should be odd." 

    return sigtools._medfilt2d(image, kernel_size)

def remez(numtaps, bands, desired, weight=None, Hz=1, type='bandpass',
          maxiter=25, grid_density=16):
    """Calculate the minimax optimal filter using Remez exchange algorithm.
    
  Description:

    Calculate the filter-coefficients for the finite impulse response
    (FIR) filter whose transfer function minimizes the maximum error
    between the desired gain and the realized gain in the specified bands
    using the remez exchange algorithm.

  Inputs:

    numtaps -- The desired number of taps in the filter.
    bands -- A montonic sequence containing the band edges.  All elements
             must be non-negative and less than 1/2 the sampling frequency
             as given by Hz.
    desired -- A sequency half the size of bands containing the desired gain
               in each of the specified bands
    weight -- A relative weighting to give to each band region.
    type --- The type of filter:
             'bandpass' : flat response in bands.
             'differentiator' : frequency proportional response in bands.

  Outputs: (out,)

    out -- A rank-1 array containing the coefficients of the optimal
           (in a minimax sense) filter.
    
    """
    # Convert type
    try:
        tnum = {'bandpass':1, 'differentiator':2}[type]
    except KeyError:
        raise ValueError, "Type must be 'bandpass', or 'differentiator'"

    # Convert weight
    if weight is None:
        weight = [1] * len(desired)

    return sigtools._remez(numtaps, bands, desired, weight, tnum, Hz,
                           maxiter, grid_density)

def lfilter(b, a, x, axis=-1, zi=None):
    """Filter data along one-dimension with an IIR or FIR filter.

  Description

    Filter a data sequence, x, using a digital filter.  This works for many
    fundamental data types (including Object type).  The filter is a direct
    form II transposed implementation of the standard difference equation
    (see "Algorithm"). 

  Inputs:

    b -- The numerator coefficient vector in a 1-D sequence.
    a -- The denominator coefficient vector in a 1-D sequence.  If a[0]
         is not 1, then both a and b are normalized by a[0].
    x -- An N-dimensional input array.
    axis -- The axis of the input data array along which to apply the
            linear filter. The filter is applied to each subarray along
            this axis (*Default* = -1)
    zi -- Initial conditions for the filter delays.  It is a vector
          (or array of vectors) of length max(len(a),len(b))-1.  If
          zi=None or is not given then initial rest is assumed.

  Outputs: (y, {zf})

    y -- The output of the digital filter.
    zf -- If zi is None, this is not returned, otherwise, zf holds the
          final filter delay values.

  Algorithm:

    The filter function is implemented as a direct II transposed structure.
    This means that the filter implements

    y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[nb]*x[n-nb]
                     - a[1]*y[n-1] + ... + a[na]*y[n-na]

    using the following difference equations:

    y[m] = b[0]*x[m] + z[0,m-1]
    z[0,m] = b[1]*x[m] + z[1,m-1] - a[1]*y[m]
    ...
    z[n-3,m] = b[n-2]*x[m] + z[n-2,m-1] - a[n-2]*y[m]
    z[n-2,m] = b[n-1]*x[m] - a[n-1]*y[m]

    where m is the output sample number and n=max(len(a),len(b)) is the
    model order.

    The rational transfer function describing this filter in the
    z-transform domain is
                                -1               -nb
                    b[0] + b[1]z  + ... + b[nb] z
            Y(z) = ---------------------------------- X(z)
                                -1               -na
                    a[0] + a[1]z  + ... + a[na] z
                    
    """
    if zi is None:
        return sigtools._linear_filter(b, a, x, axis)
    else:
        return sigtools._linear_filter(b, a, x, axis, zi)

def blackman(M):
    """blackman(M) returns the M-point Blackman window.
    """
    n = arange(0,M)
    return 0.42-0.5*cos(2.0*pi*n/(M-1)) + 0.08*cos(4.0*pi*n/(M-1))

def bartlett(M):
    """bartlett(M) returns the M-point Bartlett window.
    """
    n = arange(0,M)
    return where(less_equal(n,(M-1)/2.0),2.0*n/(M-1),2.0-2.0*n/(M-1))

def hanning(M):
    """hanning(M) returns the M-point Hanning window.
    """
    n = arange(0,M)
    return 0.5-0.5*cos(2.0*pi*n/(M-1))

def hamming(M):
    """hamming(M) returns the M-point Hamming window.
    """
    n = arange(0,M)
    return 0.54-0.46*cos(2.0*pi*n/(M-1))

def kaiser(M,beta):
    """kaiser(M, beta) returns a Kaiser window of length M with shape
    parameter beta. It depends on the cephes module for the modified bessel
    function i0.
    """
    n = arange(0,M)
    alpha = (M-1)/2.0
    return special.i0(beta * sqrt(1-((n-alpha)/alpha)**2.0))/special.i0(beta)

def hilbert(x, N=None):
    """Return the hilbert transform of x of length N.
    """
    x = Numeric.asarray(x)
    if N is None:
        N = len(x)
    if N <=0:
        raise ValueError, "N must be positive."
    if iscomplex(x):
        print "Warning: imaginary part of x ignored."
        x = real(x)
    Xf = fft(x,N,axis=0)
    h = Numeric.zeros(N)
    if N % 2 == 0:
        h[0] = h[N/2] = 1
        h[1:N/2] = 2
    else:
        h[0] = 1
        h[1:(N+1)/2] = 2

    if len(x.shape) > 1:
        h = h[:,Numeric.NewAxis]
    x = ifft(Xf*h)
    return x

def cmplx_sort(p):
    "sort roots based on magnitude."
    p = asarray(p)
    if scipy.iscomplex(p):
        indx = Numeric.argsort(abs(p))
    else:
        indx = Numeric.argsort(p)
    return Numeric.take(p,indx), indx

def unique_roots(p,tol=1e-3,rtype='min'):
    """Determine the unique roots and their multiplicities in two lists

    Inputs:

      p -- The list of roots
      tol --- The tolerance for two roots to be considered equal.
      rtype --- How to determine the returned root from the close 
                  ones:  'max': pick the maximum
                         'min': pick the minimum
                         'avg': average roots
    Outputs: (pout, mult)

      pout -- The list of sorted roots
      mult -- The multiplicity of each root
    """
    if rtype in ['max','maximum']:
        comproot = scipy.max
    elif rtype in ['min','minimum']:
        comproot = scipy.min
    elif rtype in ['avg','mean']:
        comproot = scipy.mean
    p = asarray(p)*1.0
    tol = abs(tol)
    p, indx = cmplx_sort(p)
    pout = []
    mult = []
    indx = -1
    curp = p[0] + 5*tol
    sameroots = []
    for k in range(len(p)):
        tr = p[k]
        if abs(tr-curp) < tol:
            sameroots.append(tr)
            curp = comproot(sameroots)
            pout[indx] = curp
            mult[indx] += 1
        else:
            pout.append(tr)
            curp = tr
            sameroots = [tr]            
            indx += 1
            mult.append(1)
    return array(pout), array(mult)

from scipy import real_if_close
def invres(r,p,k,tol=1e-3):
    """Compute b(s) and a(s) from partial fraction expansion: r,p,k

    If M = len(b) and N = len(a)
    
            b(s)     b[0] x**(M-1) + b[1] x**(M-2) + ... + b[M-1] 
    H(s) = ------ = ----------------------------------------------
            a(s)     b[0] x**(M-1) + b[1] x**(M-2) + ... + b[M-1]

             r[0]       r[1]             r[-1]
         = -------- + -------- + ... + --------- + k(s)
           (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like

            r[i]      r[i+1]              r[i+n-1]
          -------- + ----------- + ... + ----------- 
          (s-p[i])  (s-p[i])**2          (s-p[i])**n

    See also:  residue, poly, polyval
    """
    extra = k
    p, indx = cmplx_sort(p)
    r = Numeric.take(r,indx)
    pout, mult = unique_roots(p,tol=tol,rtype='avg')
    p = []
    for k in range(len(pout)):
        p.extend([pout[k]]*mult[k])
    a = poly(p)
    if len(extra) > 0:
        b = polymul(extra,a)
    else:
        b = [0]
    indx = 0
    for k in range(len(pout)):
        temp = []
        for l in range(len(pout)):
            if l != k:
                temp.extend([pout[l]]*mult[l])
        for m in range(mult[k]):
            t2 = temp[:]
            t2.extend([pout[k]]*(mult[k]-m-1))
            b = polyadd(b,r[indx]*poly(t2))
            indx += 1
    b = real_if_close(b)
    while Numeric.allclose(b[0], 0, rtol=1e-14) and (b.shape[-1] > 1):
        b = b[1:]
    return b, a

from scipy import factorial
def residue(b,a,tol=1e-3):
    """Compute partial-fraction expansion of b(s) / a(s).

    If M = len(b) and N = len(a)
    
            b(s)     b[0] x**(M-1) + b[1] x**(M-2) + ... + b[M-1] 
    H(s) = ------ = ----------------------------------------------
            a(s)     b[0] x**(M-1) + b[1] x**(M-2) + ... + b[M-1]

             r[0]       r[1]             r[-1]
         = -------- + -------- + ... + --------- + k(s)
           (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like 

            r[i]      r[i+1]              r[i+n-1]
          -------- + ----------- + ... + ----------- 
          (s-p[i])  (s-p[i])**2          (s-p[i])**n

    See also:  invres, poly, polyval
    """

    b,a = map(asarray,(b,a))
    k,b = polydiv(b,a)
    p = roots(a)
    r = p*0.0
    pout, mult = unique_roots(p,tol=tol,rtype='avg')
    p = []
    for n in range(len(pout)):
        p.extend([pout[n]]*mult[n])
    p = asarray(p)
    # Compute the residue from the general formula
    indx = 0
    for n in range(len(pout)):
        bn = b.copy()
        pn = []
        for l in range(len(pout)):
            if l != n:
                pn.extend([pout[l]]*mult[l])
        an = poly(pn)
        # bn(s) / an(s) is (s-po[n])**Nn * b(s) / a(s) where Nn is
        # multiplicity of pole at po[n]
        sig = mult[n]
        for m in range(sig,0,-1):
            if sig > m:
                # compute next derivative of bn(s) / an(s)
                term1 = polymul(polyder(bn,1),an)
                term2 = polymul(bn,polyder(dn))
                bn = polysub(term1,term2)
                an = polymul(an,an)                
            r[indx] = polyval(bn,pout[n]) / polyval(an,pout[n]) \
                      / factorial(sig-m)
            indx += 1
    return r, p, k


def test():
    a = [3,4,5,6,5,4]
    b = [1,2,3]
    c = convolve(a,b)
    if (Numeric.product(equal(c,[3,10,22,28,32,32,23,12]))==0):
        print "Error in convolve."

    f = [[3,4,5],[2,3,4],[1,2,5]]
    d = medfilt(f)
    if (Numeric.product(ravel(equal(d,[[0,3,0],[2,3,3],[0,2,0]])))==0):
        print "Error in medfilt."

    g = Numeric.array([[5,6,4,3],[3,5,6,2],[2,3,5,6],[1,6,9,7]],'d')
    correct = Numeric.array([[2.16374269,3.2222222222, 2.8888888889, 1.6666666667],[2.666666667, 4.33333333333, 4.44444444444, 2.8888888888],[2.222222222, 4.4444444444, 5.4444444444, 4.801066874837],[1.33333333333, 3.92735042735, 6.0712560386, 5.0404040404]])
    h = wiener(g)
    if (Numeric.abs(Numeric.product(Numeric.ravel(h-correct)))> 1e-7):
        print "Error in wiener."

    return

if __name__ == "__main__":
    test()


