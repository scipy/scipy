# Author: Travis Oliphant
# 1999 -- 2002

import warnings

import sigtools
from scipy import linalg
from scipy.fftpack import fft, ifft, ifftshift, fft2, ifft2, fftn, \
        ifftn, fftfreq
from numpy import polyadd, polymul, polydiv, polysub, roots, \
        poly, polyval, polyder, cast, asarray, isscalar, atleast_1d, \
        ones, real, real_if_close, zeros, array, arange, where, rank, \
        newaxis, product, ravel, sum, r_, iscomplexobj, take, \
        argsort, allclose, expand_dims, unique, prod, sort, reshape, \
        transpose, dot, any, mean, flipud, ndarray
import numpy as np
from scipy.misc import factorial
from windows import get_window

_modedict = {'valid':0, 'same':1, 'full':2}

_boundarydict = {'fill':0, 'pad':0, 'wrap':2, 'circular':2, 'symm':1,
                 'symmetric':1, 'reflect':4}

_SWAP_INPUTS_DEPRECATION_MSG = """\
Current default behavior of convolve and correlate functions is deprecated.

Convolve and corelate currently swap their arguments if the second argument
has dimensions larger than the first one, and the mode is relative to the input
with the largest dimension. The new behavior is to never swap the inputs, which
is what most people expects, and is how correlation is usually defined.

You can control the behavior with the old_behavior flag - the flag will
disappear in scipy 0.9.0, and the functions will then implement the new
behavior only."""

def _valfrommode(mode):
    try:
        val = _modedict[mode]
    except KeyError:
        if mode not in [0,1,2]:
            raise ValueError, "Acceptable mode flags are 'valid' (0)," \
                  "'same' (1), or 'full' (2)."
        val = mode
    return val

def _bvalfromboundary(boundary):
    try:
        val = _boundarydict[boundary] << 2
    except KeyError:
        if val not in [0,1,2] :
            raise ValueError, "Acceptable boundary flags are 'fill', 'wrap'" \
                  " (or 'circular'), \n  and 'symm' (or 'symmetric')."
        val = boundary << 2
    return val


def correlate(in1, in2, mode='full', old_behavior=True):
    """
    Cross-correlate two N-dimensional arrays.

    Cross-correlate in1 and in2 with the output size determined by the mode
    argument.

    Parameters
    ----------
    in1: array
        first input.
    in2: array
        second input. Should have the same number of dimensions as in1.
    mode: str {'valid', 'same', 'full'}
        a string indicating the size of the output:
            - 'valid': the output consists only of those elements that do not
            rely on the zero-padding.
            - 'same': the output is the same size as the largest input centered
              with respect to the 'full' output.
            - 'full': the output is the full discrete linear cross-correlation
              of the inputs. (Default)
    old_behavior: bool
        If True (default), the old behavior of correlate is implemented:
            - if in1.size < in2.size, in1 and in2 are swapped (correlate(in1,
              in2) == correlate(in2, in1))
            - For complex inputs, the conjugate is not taken for in2
        If False, the new, conventional definition of correlate is implemented.

    Returns
    -------
    out: array
        an N-dimensional array containing a subset of the discrete linear
        cross-correlation of in1 with in2.

    Notes
    -----
    The correlation z of two arrays x and y of rank d is defined as

      z[...,k,...] = sum[..., i_l, ...]
            x[..., i_l,...] * conj(y[..., i_l + k,...])

    """
    val = _valfrommode(mode)

    if old_behavior:
        warnings.warn(DeprecationWarning(_SWAP_INPUTS_DEPRECATION_MSG))
        if np.iscomplexobj(in2):
            in2 = in2.conjugate()
        if in1.size < in2.size:
            swp = in2
            in2 = in1
            in1 = swp

    if mode == 'valid':
        ps = [i - j + 1 for i, j in zip(in1.shape, in2.shape)]
        out = np.empty(ps, in1.dtype)
        for i in range(len(ps)):
            if ps[i] <= 0:
                raise ValueError("Dimension of x(%d) < y(%d) " \
                                 "not compatible with valid mode" % \
                                 (in1.shape[i], in2.shape[i]))

        z = sigtools._correlateND(in1, in2, out, val)
    else:
        ps = [i + j - 1 for i, j in zip(in1.shape, in2.shape)]
        # zero pad input
        in1zpadded = np.zeros(ps, in1.dtype)
        sc = [slice(0, i) for i in in1.shape]
        in1zpadded[sc] = in1.copy()

        if mode == 'full':
            out = np.empty(ps, in1.dtype)
            z = sigtools._correlateND(in1zpadded, in2, out, val)
        elif mode == 'same':
            out = np.empty(in1.shape, in1.dtype)

            z = sigtools._correlateND(in1zpadded, in2, out, val)
        else:
            raise ValueError("Uknown mode %s" % mode)

    return z

def _centered(arr, newsize):
    # Return the center newsize portion of the array.
    newsize = asarray(newsize)
    currsize = array(arr.shape)
    startind = (currsize - newsize) / 2
    endind = startind + newsize
    myslice = [slice(startind[k], endind[k]) for k in range(len(endind))]
    return arr[tuple(myslice)]

def fftconvolve(in1, in2, mode="full"):
    """Convolve two N-dimensional arrays using FFT. See convolve.

    """
    s1 = array(in1.shape)
    s2 = array(in2.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex) or
                      np.issubdtype(in2.dtype, np.complex))
    size = s1+s2-1

    # Always use 2**n-sized FFT
    fsize = 2**np.ceil(np.log2(size))
    IN1 = fftn(in1,fsize)
    IN1 *= fftn(in2,fsize)
    fslice = tuple([slice(0, int(sz)) for sz in size])
    ret = ifftn(IN1)[fslice].copy()
    del IN1
    if not complex_result:
        ret = ret.real
    if mode == "full":
        return ret
    elif mode == "same":
        if product(s1,axis=0) > product(s2,axis=0):
            osize = s1
        else:
            osize = s2
        return _centered(ret,osize)
    elif mode == "valid":
        return _centered(ret,abs(s2-s1)+1)


def convolve(in1, in2, mode='full', old_behavior=True):
    """
    Convolve two N-dimensional arrays.

    Convolve in1 and in2 with output size determined by mode.

    Parameters
    ----------
    in1: array
        first input.
    in2: array
        second input. Should have the same number of dimensions as in1.
    mode: str {'valid', 'same', 'full'}
        a string indicating the size of the output:

        ``valid`` : the output consists only of those elements that do not
           rely on the zero-padding.

        ``same`` : the output is the same size as the largest input centered
           with respect to the 'full' output.

        ``full`` : the output is the full discrete linear cross-correlation
           of the inputs. (Default)


    Returns
    -------
    out: array
        an N-dimensional array containing a subset of the discrete linear
        cross-correlation of in1 with in2.

    """
    volume = asarray(in1)
    kernel = asarray(in2)

    if rank(volume) == rank(kernel) == 0:
        return volume*kernel
    elif not volume.ndim == kernel.ndim:
        raise ValueError("in1 and in2 should have the same rank")

    slice_obj = [slice(None,None,-1)]*len(kernel.shape)

    if old_behavior:
        warnings.warn(DeprecationWarning(_SWAP_INPUTS_DEPRECATION_MSG))
        if (product(kernel.shape,axis=0) > product(volume.shape,axis=0)):
            temp = kernel
            kernel = volume
            volume = temp
            del temp

        return correlate(volume, kernel[slice_obj], mode, old_behavior=True)
    else:
        if mode == 'valid':
            for d1, d2 in zip(volume.shape, kernel.shape):
                if not d1 >= d2:
                    raise ValueError(
                        "in1 should have at least as many items as in2 in " \
                        "every dimension for valid mode.")
        if np.iscomplexobj(kernel):
            return correlate(volume, kernel[slice_obj].conj(), mode, old_behavior=False)
        else:
            return correlate(volume, kernel[slice_obj], mode, old_behavior=False)

def order_filter(a, domain, rank):
    """
    Perform an order filter on an N-dimensional array.

    Description:

      Perform an order filter on the array in.  The domain argument acts as a
      mask centered over each pixel.  The non-zero elements of domain are
      used to select elements surrounding each input pixel which are placed
      in a list.   The list is sorted, and the output for that pixel is the
      element corresponding to rank in the sorted list.

    Parameters
    ----------
      in -- an N-dimensional input array.
      domain -- a mask array with the same number of dimensions as in.  Each
                dimension should have an odd number of elements.
      rank -- an non-negative integer which selects the element from the
              sorted list (0 corresponds to the largest element, 1 is the
              next largest element, etc.)

    Returns
    -------
      out -- the results of the order filter in an array with the same
             shape as in.

    """
    domain = asarray(domain)
    size = domain.shape
    for k in range(len(size)):
        if (size[k] % 2) != 1:
            raise ValueError, "Each dimension of domain argument " \
                  "should have an odd number of elements."
    return sigtools._order_filterND(a, domain, rank)


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
    volume = atleast_1d(volume)
    if kernel_size is None:
        kernel_size = [3] * len(volume.shape)
    kernel_size = asarray(kernel_size)
    if len(kernel_size.shape) == 0:
        kernel_size = [kernel_size.item()] * len(volume.shape)
    kernel_size = asarray(kernel_size)

    for k in range(len(volume.shape)):
        if (kernel_size[k] % 2) != 1:
            raise ValueError, "Each element of kernel_size should be odd."

    domain = ones(kernel_size)

    numels = product(kernel_size,axis=0)
    order = int(numels/2)
    return sigtools._order_filterND(volume,domain,order)


def wiener(im,mysize=None,noise=None):
    """
    Perform a Wiener filter on an N-dimensional array.

    Description:

      Apply a Wiener filter to the N-dimensional array in.

    Inputs:

      in -- an N-dimensional array.
      kernel_size -- A scalar or an N-length list giving the size of the
                     Wiener filter window in each dimension.  Elements of
                     kernel_size should be odd.  If kernel_size is a scalar,
                     then this scalar is used as the size in each dimension.
      noise -- The noise-power to use.  If None, then noise is estimated as
               the average of the local variance of the input.

    Outputs: (out,)

      out -- Wiener filtered result with the same shape as in.

    """
    im = asarray(im)
    if mysize is None:
        mysize = [3] * len(im.shape)
    mysize = asarray(mysize);

    # Estimate the local mean
    lMean = correlate(im,ones(mysize), 'same', old_behavior=False) / product(mysize,axis=0)

    # Estimate the local variance
    lVar = correlate(im**2,ones(mysize), 'same', old_behavior=False) / product(mysize,axis=0) - lMean**2

    # Estimate the noise power if needed.
    if noise==None:
        noise = mean(ravel(lVar),axis=0)

    res = (im - lMean)
    res *= (1-noise / lVar)
    res += lMean
    out = where(lVar < noise, lMean, res)

    return out


def convolve2d(in1, in2, mode='full', boundary='fill', fillvalue=0, old_behavior=True):
    """Convolve two 2-dimensional arrays.

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
    if old_behavior:
        warnings.warn(DeprecationWarning(_SWAP_INPUTS_DEPRECATION_MSG))

    if old_behavior:
        warnings.warn(DeprecationWarning(_SWAP_INPUTS_DEPRECATION_MSG))
        if (product(np.shape(in2),axis=0) > product(np.shape(in1),axis=0)):
            temp = in1
            in1 = in2
            in2 = temp
            del temp
    else:
        if mode == 'valid':
            for d1, d2 in zip(np.shape(in1), np.shape(in2)):
                if not d1 >= d2:
                    raise ValueError(
                        "in1 should have at least as many items as in2 in " \
                        "every dimension for valid mode.")

    val = _valfrommode(mode)
    bval = _bvalfromboundary(boundary)

    return sigtools._convolve2d(in1,in2,1,val,bval,fillvalue)

def correlate2d(in1, in2, mode='full', boundary='fill', fillvalue=0, old_behavior=True):
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
    if old_behavior:
        warnings.warn(DeprecationWarning(_SWAP_INPUTS_DEPRECATION_MSG))
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
    image = asarray(input)
    if kernel_size is None:
        kernel_size = [3] * 2
    kernel_size = asarray(kernel_size)
    if len(kernel_size.shape) == 0:
        kernel_size = [kernel_size.item()] * 2
    kernel_size = asarray(kernel_size)

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

    bands = asarray(bands).copy()
    return sigtools._remez(numtaps, bands, desired, weight, tnum, Hz,
                           maxiter, grid_density)

def lfilter(b, a, x, axis=-1, zi=None):
    """
    Filter data along one-dimension with an IIR or FIR filter.

    Filter a data sequence, x, using a digital filter.  This works for many
    fundamental data types (including Object type).  The filter is a direct
    form II transposed implementation of the standard difference equation
    (see Notes).

    Parameters
    ----------
    b : array_like
        The numerator coefficient vector in a 1-D sequence.
    a : array_like
        The denominator coefficient vector in a 1-D sequence.  If a[0]
        is not 1, then both a and b are normalized by a[0].
    x : array_like
        An N-dimensional input array.
    axis : int
        The axis of the input data array along which to apply the
        linear filter. The filter is applied to each subarray along
        this axis (*Default* = -1)
    zi : array_like (optional)
        Initial conditions for the filter delays.  It is a vector
        (or array of vectors for an N-dimensional input) of length
        max(len(a),len(b))-1.  If zi=None or is not given then initial
        rest is assumed.  SEE signal.lfiltic for more information.

    Returns
    -------
    y : array
        The output of the digital filter.
    zf : array (optional)
        If zi is None, this is not returned, otherwise, zf holds the
        final filter delay values.

    Notes
    -----
    The filter function is implemented as a direct II transposed structure.
    This means that the filter implements

    ::

       a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[nb]*x[n-nb]
                               - a[1]*y[n-1] - ... - a[na]*y[n-na]

    using the following difference equations::

         y[m] = b[0]*x[m] + z[0,m-1]
         z[0,m] = b[1]*x[m] + z[1,m-1] - a[1]*y[m]
         ...
         z[n-3,m] = b[n-2]*x[m] + z[n-2,m-1] - a[n-2]*y[m]
         z[n-2,m] = b[n-1]*x[m] - a[n-1]*y[m]

    where m is the output sample number and n=max(len(a),len(b)) is the
    model order.

    The rational transfer function describing this filter in the
    z-transform domain is::

                             -1               -nb
                 b[0] + b[1]z  + ... + b[nb] z
         Y(z) = ---------------------------------- X(z)
                             -1               -na
                 a[0] + a[1]z  + ... + a[na] z

    """
    if isscalar(a):
        a = [a]
    if zi is None:
        return sigtools._linear_filter(b, a, x, axis)
    else:
        return sigtools._linear_filter(b, a, x, axis, zi)

def lfiltic(b,a,y,x=None):
    """
    Construct initial conditions for lfilter

    Given a linear filter (b,a) and initial conditions on the output y
    and the input x, return the inital conditions on the state vector zi
    which is used by lfilter to generate the output given the input.

    If M=len(b)-1 and N=len(a)-1.  Then, the initial conditions are given
    in the vectors x and y as::

     x = {x[-1],x[-2],...,x[-M]}
     y = {y[-1],y[-2],...,y[-N]}

    If x is not given, its inital conditions are assumed zero.
    If either vector is too short, then zeros are added
    to achieve the proper length.

    The output vector zi contains::

     zi = {z_0[-1], z_1[-1], ..., z_K-1[-1]}  where K=max(M,N).

    """
    N = np.size(a)-1
    M = np.size(b)-1
    K = max(M,N)
    y = asarray(y)
    zi = zeros(K,y.dtype.char)
    if x is None:
        x = zeros(M,y.dtype.char)
    else:
        x = asarray(x)
        L = np.size(x)
        if L < M:
            x = r_[x,zeros(M-L)]
    L = np.size(y)
    if L < N:
        y = r_[y,zeros(N-L)]

    for m in range(M):
        zi[m] = sum(b[m+1:]*x[:M-m],axis=0)

    for m in range(N):
        zi[m] -= sum(a[m+1:]*y[:N-m],axis=0)

    return zi

def deconvolve(signal, divisor):
    """Deconvolves divisor out of signal.

    """
    num = atleast_1d(signal)
    den = atleast_1d(divisor)
    N = len(num)
    D = len(den)
    if D > N:
        quot = [];
        rem = num;
    else:
        input = ones(N-D+1, float)
        input[1:] = 0
        quot = lfilter(num, den, input)
        rem = num - convolve(den, quot, mode='full')
    return quot, rem



def hilbert(x, N=None, axis=-1):
    """Compute the analytic signal.

    The transformation is done along the last axis by default.

    Parameters
    ----------
    x : array-like
        Signal data
    N : int, optional
        Number of Fourier components. Default: ``x.shape[axis]``
    axis : int, optional
        

    Returns
    -------
    xa : ndarray
        Analytic signal of `x`, of each 1d array along axis

    Notes
    -----
    The analytic signal `x_a(t)` of `x(t)` is::

        x_a = F^{-1}(F(x) 2U) = x + i y

    where ``F`` is the Fourier transform, ``U`` the unit step function,
    and ``y`` the Hilbert transform of ``x``. [1]

    changes in scipy 0.8.0: new axis argument, new default axis=-1

    References
    ----------
    .. [1] Wikipedia, "Analytic signal".
           http://en.wikipedia.org/wiki/Analytic_signal

    """
    x = asarray(x)
    if N is None:
        N = x.shape[axis]
    if N <=0:
        raise ValueError, "N must be positive."
    if iscomplexobj(x):
        print "Warning: imaginary part of x ignored."
        x = real(x)
    Xf = fft(x, N, axis=axis)
    h = zeros(N)
    if N % 2 == 0:
        h[0] = h[N/2] = 1
        h[1:N/2] = 2
    else:
        h[0] = 1
        h[1:(N+1)/2] = 2

    if len(x.shape) > 1:
        ind = [newaxis]*x.ndim
        ind[axis] = slice(None)
        h = h[ind]
    x = ifft(Xf*h, axis=axis)
    return x

def hilbert2(x,N=None):
    """
    Compute the '2-D' analytic signal of `x`


    Parameters
    ----------
    x : array_like
        2-D signal data.
    N : int, optional
        Number of Fourier components. Default is ``x.shape``

    Returns
    -------
    xa : ndarray
        Analytic signal of `x` taken along axes (0,1).

    References
    ----------
    .. [1] Wikipedia, "Analytic signal",
        http://en.wikipedia.org/wiki/Analytic_signal

    """
    x = asarray(x)
    x = asarray(x)
    if N is None:
        N = x.shape
    if len(N) < 2:
        if N <=0:
            raise ValueError, "N must be positive."
        N = (N,N)
    if iscomplexobj(x):
        print "Warning: imaginary part of x ignored."
        x = real(x)
    Xf = fft2(x,N,axes=(0,1))
    h1 = zeros(N[0],'d')
    h2 = zeros(N[1],'d')
    for p in range(2):
        h = eval("h%d"%(p+1))
        N1 = N[p]
        if N1 % 2 == 0:
            h[0] = h[N1/2] = 1
            h[1:N1/2] = 2
        else:
            h[0] = 1
            h[1:(N1+1)/2] = 2
        exec("h%d = h" % (p+1), globals(), locals())

    h = h1[:,newaxis] * h2[newaxis,:]
    k = len(x.shape)
    while k > 2:
        h = h[:, newaxis]
        k -= 1
    x = ifft2(Xf*h,axes=(0,1))
    return x


def cmplx_sort(p):
    "sort roots based on magnitude."
    p = asarray(p)
    if iscomplexobj(p):
        indx = argsort(abs(p))
    else:
        indx = argsort(p)
    return take(p,indx,0), indx

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
        comproot = np.maximum
    elif rtype in ['min','minimum']:
        comproot = np.minimum
    elif rtype in ['avg','mean']:
        comproot = np.mean
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


def invres(r,p,k,tol=1e-3,rtype='avg'):
    """Compute b(s) and a(s) from partial fraction expansion: r,p,k

    If M = len(b) and N = len(a)

            b(s)     b[0] x**(M-1) + b[1] x**(M-2) + ... + b[M-1]
    H(s) = ------ = ----------------------------------------------
            a(s)     a[0] x**(N-1) + a[1] x**(N-2) + ... + a[N-1]

             r[0]       r[1]             r[-1]
         = -------- + -------- + ... + --------- + k(s)
           (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like

            r[i]      r[i+1]              r[i+n-1]
          -------- + ----------- + ... + -----------
          (s-p[i])  (s-p[i])**2          (s-p[i])**n

    See Also
    --------
    residue, poly, polyval, unique_roots

    """
    extra = k
    p, indx = cmplx_sort(p)
    r = take(r,indx,0)
    pout, mult = unique_roots(p,tol=tol,rtype=rtype)
    p = []
    for k in range(len(pout)):
        p.extend([pout[k]]*mult[k])
    a = atleast_1d(poly(p))
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
    while allclose(b[0], 0, rtol=1e-14) and (b.shape[-1] > 1):
        b = b[1:]
    return b, a

def residue(b,a,tol=1e-3,rtype='avg'):
    """Compute partial-fraction expansion of b(s) / a(s).

    If M = len(b) and N = len(a)

            b(s)     b[0] s**(M-1) + b[1] s**(M-2) + ... + b[M-1]
    H(s) = ------ = ----------------------------------------------
            a(s)     a[0] s**(N-1) + a[1] s**(N-2) + ... + a[N-1]

             r[0]       r[1]             r[-1]
         = -------- + -------- + ... + --------- + k(s)
           (s-p[0])   (s-p[1])         (s-p[-1])

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like

            r[i]      r[i+1]              r[i+n-1]
          -------- + ----------- + ... + -----------
          (s-p[i])  (s-p[i])**2          (s-p[i])**n

    Returns
    -------
    r : ndarray
        Residues
    p : ndarray
        Poles
    k : ndarray
        Coefficients of the direct polynomial term.

    See Also
    --------
    invres, poly, polyval, unique_roots

    """

    b,a = map(asarray,(b,a))
    rscale = a[0]
    k,b = polydiv(b,a)
    p = roots(a)
    r = p*0.0
    pout, mult = unique_roots(p,tol=tol,rtype=rtype)
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
        an = atleast_1d(poly(pn))
        # bn(s) / an(s) is (s-po[n])**Nn * b(s) / a(s) where Nn is
        # multiplicity of pole at po[n]
        sig = mult[n]
        for m in range(sig,0,-1):
            if sig > m:
                # compute next derivative of bn(s) / an(s)
                term1 = polymul(polyder(bn,1),an)
                term2 = polymul(bn,polyder(an,1))
                bn = polysub(term1,term2)
                an = polymul(an,an)
            r[indx+m-1] = polyval(bn,pout[n]) / polyval(an,pout[n]) \
                          / factorial(sig-m)
        indx += sig
    return r/rscale, p, k

def residuez(b,a,tol=1e-3,rtype='avg'):
    """Compute partial-fraction expansion of b(z) / a(z).

    If M = len(b) and N = len(a)

            b(z)     b[0] + b[1] z**(-1) + ... + b[M-1] z**(-M+1)
    H(z) = ------ = ----------------------------------------------
            a(z)     a[0] + a[1] z**(-1) + ... + a[N-1] z**(-N+1)

                 r[0]                   r[-1]
         = --------------- + ... + ---------------- + k[0] + k[1]z**(-1) ...
           (1-p[0]z**(-1))         (1-p[-1]z**(-1))

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like

               r[i]              r[i+1]                    r[i+n-1]
          -------------- + ------------------ + ... + ------------------
          (1-p[i]z**(-1))  (1-p[i]z**(-1))**2         (1-p[i]z**(-1))**n

    See also:  invresz, poly, polyval, unique_roots
    """
    b,a = map(asarray,(b,a))
    gain = a[0]
    brev, arev = b[::-1],a[::-1]
    krev,brev = polydiv(brev,arev)
    if krev == []:
        k = []
    else:
        k = krev[::-1]
    b = brev[::-1]
    p = roots(a)
    r = p*0.0
    pout, mult = unique_roots(p,tol=tol,rtype=rtype)
    p = []
    for n in range(len(pout)):
        p.extend([pout[n]]*mult[n])
    p = asarray(p)
    # Compute the residue from the general formula (for discrete-time)
    #  the polynomial is in z**(-1) and the multiplication is by terms
    #  like this (1-p[i] z**(-1))**mult[i].  After differentiation,
    #  we must divide by (-p[i])**(m-k) as well as (m-k)!
    indx = 0
    for n in range(len(pout)):
        bn = brev.copy()
        pn = []
        for l in range(len(pout)):
            if l != n:
                pn.extend([pout[l]]*mult[l])
        an = atleast_1d(poly(pn))[::-1]
        # bn(z) / an(z) is (1-po[n] z**(-1))**Nn * b(z) / a(z) where Nn is
        # multiplicity of pole at po[n] and b(z) and a(z) are polynomials.
        sig = mult[n]
        for m in range(sig,0,-1):
            if sig > m:
                # compute next derivative of bn(s) / an(s)
                term1 = polymul(polyder(bn,1),an)
                term2 = polymul(bn,polyder(an,1))
                bn = polysub(term1,term2)
                an = polymul(an,an)
            r[indx+m-1] = polyval(bn,1.0/pout[n]) / polyval(an,1.0/pout[n]) \
                          / factorial(sig-m) / (-pout[n])**(sig-m)
        indx += sig
    return r/gain, p, k

def invresz(r,p,k,tol=1e-3,rtype='avg'):
    """Compute b(z) and a(z) from partial fraction expansion: r,p,k

    If M = len(b) and N = len(a)

            b(z)     b[0] + b[1] z**(-1) + ... + b[M-1] z**(-M+1)
    H(z) = ------ = ----------------------------------------------
            a(z)     a[0] + a[1] z**(-1) + ... + a[N-1] z**(-N+1)

                 r[0]                   r[-1]
         = --------------- + ... + ---------------- + k[0] + k[1]z**(-1) ...
           (1-p[0]z**(-1))         (1-p[-1]z**(-1))

    If there are any repeated roots (closer than tol), then the partial
    fraction expansion has terms like

               r[i]              r[i+1]                    r[i+n-1]
          -------------- + ------------------ + ... + ------------------
          (1-p[i]z**(-1))  (1-p[i]z**(-1))**2         (1-p[i]z**(-1))**n

    See also:  residuez, poly, polyval, unique_roots
    """
    extra = asarray(k)
    p, indx = cmplx_sort(p)
    r = take(r,indx,0)
    pout, mult = unique_roots(p,tol=tol,rtype=rtype)
    p = []
    for k in range(len(pout)):
        p.extend([pout[k]]*mult[k])
    a = atleast_1d(poly(p))
    if len(extra) > 0:
        b = polymul(extra,a)
    else:
        b = [0]
    indx = 0
    brev = asarray(b)[::-1]
    for k in range(len(pout)):
        temp = []
        # Construct polynomial which does not include any of this root
        for l in range(len(pout)):
            if l != k:
                temp.extend([pout[l]]*mult[l])
        for m in range(mult[k]):
            t2 = temp[:]
            t2.extend([pout[k]]*(mult[k]-m-1))
            brev = polyadd(brev,(r[indx]*poly(t2))[::-1])
            indx += 1
    b = real_if_close(brev[::-1])
    return b, a


def resample(x,num,t=None,axis=0,window=None):
    """Resample to num samples using Fourier method along the given axis.

    The resampled signal starts at the same value of x but is sampled
    with a spacing of len(x) / num * (spacing of x).  Because a
    Fourier method is used, the signal is assumed periodic.

    Window controls a Fourier-domain window that tapers the Fourier
    spectrum before zero-padding to alleviate ringing in the resampled
    values for sampled signals you didn't intend to be interpreted as
    band-limited.

    If window is a function, then it is called with a vector of inputs
    indicating the frequency bins (i.e. fftfreq(x.shape[axis]) )

    If window is an array of the same length as x.shape[axis] it is
    assumed to be the window to be applied directly in the Fourier
    domain (with dc and low-frequency first).

    If window is a string then use the named window.  If window is a
    float, then it represents a value of beta for a kaiser window.  If
    window is a tuple, then the first component is a string
    representing the window, and the next arguments are parameters for
    that window.  
    
    Possible windows are:
           'flattop'        -- 'flat', 'flt'
           'boxcar'         -- 'ones', 'box'
           'triang'         -- 'traing', 'tri'
           'parzen'         -- 'parz', 'par'
           'bohman'         -- 'bman', 'bmn'
           'blackmanharris' -- 'blackharr', 'bkh'
           'nuttall',       -- 'nutl', 'nut'
           'barthann'       -- 'brthan', 'bth'
           'blackman'       -- 'black',   'blk'
           'hamming'        -- 'hamm',    'ham'
           'bartlett'       -- 'bart',    'brt'
           'hanning'        -- 'hann',    'han'
           ('kaiser', beta)                 -- 'ksr'
           ('gaussian', std)                -- 'gauss',   'gss' 
           ('general gauss', power, width)  -- 'general', 'ggs'
           ('slepian', width)               -- 'slep', 'optimal', 'dss'

    The first sample of the returned vector is the same as the first
    sample of the input vector, the spacing between samples is changed
    from dx to

        dx * len(x) / num

    If t is not None, then it represents the old sample positions, and the new
    sample positions will be returned as well as the new samples.
    """
    x = asarray(x)
    X = fft(x,axis=axis)
    Nx = x.shape[axis]
    if window is not None:
        if callable(window):
            W = window(fftfreq(Nx))
        elif isinstance(window, ndarray) and window.shape == (Nx,):
            W = window
        else:
            W = ifftshift(get_window(window,Nx))
        newshape = ones(len(x.shape))
        newshape[axis] = len(W)
        W.shape = newshape
        X = X*W
    sl = [slice(None)]*len(x.shape)
    newshape = list(x.shape)
    newshape[axis] = num
    N = int(np.minimum(num,Nx))
    Y = zeros(newshape,'D')
    sl[axis] = slice(0,(N+1)/2)
    Y[sl] = X[sl]
    sl[axis] = slice(-(N-1)/2,None)
    Y[sl] = X[sl]
    y = ifft(Y,axis=axis)*(float(num)/float(Nx))

    if x.dtype.char not in ['F','D']:
        y = y.real

    if t is None:
        return y
    else:
        new_t = arange(0,num)*(t[1]-t[0])* Nx / float(num) + t[0]
        return y, new_t

def detrend(data, axis=-1, type='linear', bp=0):
    """Remove linear trend along axis from data.

    If type is 'constant' then remove mean only.

    If bp is given, then it is a sequence of points at which to
       break a piecewise-linear fit to the data.

    """
    if type not in ['linear','l','constant','c']:
        raise ValueError, "Trend type must be linear or constant"
    data = asarray(data)
    dtype = data.dtype.char
    if dtype not in 'dfDF':
        dtype = 'd'
    if type in ['constant','c']:
        ret = data - expand_dims(mean(data,axis),axis)
        return ret
    else:
        dshape = data.shape
        N = dshape[axis]
        bp = sort(unique(r_[0,bp,N]))
        if any(bp > N):
            raise ValueError, "Breakpoints must be less than length " \
                  "of data along given axis."
        Nreg = len(bp) - 1
        # Restructure data so that axis is along first dimension and
        #  all other dimensions are collapsed into second dimension
        rnk = len(dshape)
        if axis < 0: axis = axis + rnk
        newdims = r_[axis,0:axis,axis+1:rnk]
        newdata = reshape(transpose(data, tuple(newdims)),
                          (N, prod(dshape, axis=0)/N))
        newdata = newdata.copy()  # make sure we have a copy
        if newdata.dtype.char not in 'dfDF':
            newdata = newdata.astype(dtype)
        # Find leastsq fit and remove it for each piece
        for m in range(Nreg):
            Npts = bp[m+1] - bp[m]
            A = ones((Npts,2),dtype)
            A[:,0] = cast[dtype](arange(1,Npts+1)*1.0/Npts)
            sl = slice(bp[m],bp[m+1])
            coef,resids,rank,s = linalg.lstsq(A,newdata[sl])
            newdata[sl] = newdata[sl] - dot(A,coef)
        # Put data back in original shape.
        tdshape = take(dshape,newdims,0)
        ret = reshape(newdata,tuple(tdshape))
        vals = range(1,rnk)
        olddims = vals[:axis] + [0] + vals[axis:]
        ret = transpose(ret,tuple(olddims))
        return ret

def lfilter_zi(b,a):
    #compute the zi state from the filter parameters. see [Gust96].

    #Based on:
    # [Gust96] Fredrik Gustafsson, Determining the initial states in
    #          forward-backward filtering, IEEE Transactions on
    #          Signal Processing, pp. 988--992, April 1996,
    #          Volume 44, Issue 4

    n=max(len(a),len(b))

    zin = (np.eye(n-1) - np.hstack((-a[1:n,newaxis],
                                    np.vstack((np.eye(n-2),zeros(n-2))))))

    zid=  b[1:n] - a[1:n]*b[0]

    zi_matrix=linalg.inv(zin)*(np.matrix(zid).transpose())
    zi_return=[]

    #convert the result into a regular array (not a matrix)
    for i in range(len(zi_matrix)):
        zi_return.append(float(zi_matrix[i][0]))

    return array(zi_return)

def filtfilt(b,a,x):
    b, a, x = map(asarray, [b, a, x])
    # FIXME:  For now only accepting 1d arrays
    ntaps=max(len(a),len(b))
    edge=ntaps*3

    if x.ndim != 1:
        raise ValueError, "filtfilt only accepts 1-d arrays."

    #x must be bigger than edge
    if x.size < edge:
        raise ValueError, "Input vector needs to be bigger than " \
              "3 * max(len(a),len(b)."

    if len(a) < ntaps:
        a=r_[a,zeros(len(b)-len(a))]

    if len(b) < ntaps:
        b=r_[b,zeros(len(a)-len(b))]

    zi = lfilter_zi(b,a)

    #Grow the signal to have edges for stabilizing
    #the filter with inverted replicas of the signal
    s=r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
    #in the case of one go we only need one of the extrems
    # both are needed for filtfilt

    (y,zf)=lfilter(b,a,s,-1,zi*s[0])

    (y,zf)=lfilter(b,a,flipud(y),-1,zi*y[-1])

    return flipud(y[edge-1:-edge+1])


from scipy.signal.filter_design import cheby1, firwin

def decimate(x, q, n=None, ftype='iir', axis=-1):
    """downsample the signal x by an integer factor q, using an order n filter

    By default an order 8 Chebyshev type I filter is used or a 30 point FIR
    filter with hamming window if ftype is 'fir'.

    Parameters
    ----------
    x : N-d array
      the signal to be downsampled
    q : int
      the downsampling factor
    n : int or None
      the order of the filter (1 less than the length for 'fir')
    ftype : {'iir' or 'fir'}
      the type of the lowpass filter
    axis : int
      the axis along which to decimate

    Returns
    -------
    y : N-d array
      the down-sampled signal

    See also:  resample
    """

    if not isinstance(q, int):
        raise TypeError, "q must be an integer"

    if n is None:
        if ftype == 'fir':
            n = 30
        else:
            n = 8

    if ftype == 'fir':
        b = firwin(n+1, 1./q, window='hamming')
        a = 1.
    else:
        b, a = cheby1(n, 0.05, 0.8/q)

    y = lfilter(b, a, x, axis=axis)

    sl = [None]*y.ndim
    sl[axis] = slice(None, None, q)
    return y[sl]
