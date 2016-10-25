"""
Functions which are common and require SciPy Base and Level 1 SciPy
(special, linalg)
"""

from __future__ import division, print_function, absolute_import

import warnings

import numpy
import numpy as np
from numpy import (exp, log, asarray, arange, newaxis, hstack, product, array,
                   zeros, eye, poly1d, r_, fromstring, isfinite,
                   squeeze, amax, reshape, sign, broadcast_arrays)

from scipy._lib._util import _asarray_validated
from scipy.optimize import minimize_scalar

__all__ = ['logsumexp', 'inversefunc', 'central_diff_weights', 'derivative', 'pade', 'lena',
           'ascent', 'face']


def logsumexp(a, axis=None, b=None, keepdims=False, return_sign=False):
    """Compute the log of the sum of exponentials of input elements.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : None or int or tuple of ints, optional
        Axis or axes over which the sum is taken. By default `axis` is None,
        and all elements are summed.

        .. versionadded:: 0.11.0
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left in the
        result as dimensions with size one. With this option, the result
        will broadcast correctly against the original array.

        .. versionadded:: 0.15.0
    b : array-like, optional
        Scaling factor for exp(`a`) must be of the same shape as `a` or
        broadcastable to `a`. These values may be negative in order to
        implement subtraction.

        .. versionadded:: 0.12.0
    return_sign : bool, optional
        If this is set to True, the result will be a pair containing sign
        information; if False, results that are negative will be returned
        as NaN. Default is False (no sign information).

        .. versionadded:: 0.16.0
    Returns
    -------
    res : ndarray
        The result, ``np.log(np.sum(np.exp(a)))`` calculated in a numerically
        more stable way. If `b` is given then ``np.log(np.sum(b*np.exp(a)))``
        is returned.
    sgn : ndarray
        If return_sign is True, this will be an array of floating-point
        numbers matching res and +1, 0, or -1 depending on the sign
        of the result. If False, only one result is returned.

    See Also
    --------
    numpy.logaddexp, numpy.logaddexp2

    Notes
    -----
    Numpy has a logaddexp function which is very similar to `logsumexp`, but
    only handles two arguments. `logaddexp.reduce` is similar to this
    function, but may be less stable.

    Examples
    --------
    >>> from scipy.misc import logsumexp
    >>> a = np.arange(10)
    >>> np.log(np.sum(np.exp(a)))
    9.4586297444267107
    >>> logsumexp(a)
    9.4586297444267107

    With weights

    >>> a = np.arange(10)
    >>> b = np.arange(10, 0, -1)
    >>> logsumexp(a, b=b)
    9.9170178533034665
    >>> np.log(np.sum(b*np.exp(a)))
    9.9170178533034647

    Returning a sign flag

    >>> logsumexp([1,2],b=[1,-1],return_sign=True)
    (1.5413248546129181, -1.0)

    Notice that `logsumexp` does not directly support masked arrays. To use it
    on a masked array, convert the mask into zero weights:

    >>> a = np.ma.array([np.log(2), 2, np.log(3)],
    ...                  mask=[False, True, False])
    >>> b = (~a.mask).astype(int)
    >>> logsumexp(a.data, b=b), np.log(5)
    1.6094379124341005, 1.6094379124341005

    """
    a = _asarray_validated(a, check_finite=False)
    if b is not None:
        a, b = broadcast_arrays(a,b)
        if np.any(b == 0):
            a = a + 0.  # promote to at least float
            a[b == 0] = -np.inf

    a_max = amax(a, axis=axis, keepdims=True)

    if a_max.ndim > 0:
        a_max[~isfinite(a_max)] = 0
    elif not isfinite(a_max):
        a_max = 0

    if b is not None:
        b = asarray(b)
        tmp = b * exp(a - a_max)
    else:
        tmp = exp(a - a_max)

    # suppress warnings about log of zero
    with np.errstate(divide='ignore'):
        s = np.sum(tmp, axis=axis, keepdims=keepdims)
        if return_sign:
            sgn = sign(s)
            s *= sgn  # /= makes more sense but we need zero -> zero
        out = log(s)

    if not keepdims:
        a_max = squeeze(a_max, axis=axis)
    out += a_max

    if return_sign:
        return out, sgn
    else:
        return out


def central_diff_weights(Np, ndiv=1):
    """
    Return weights for an Np-point central derivative.

    Assumes equally-spaced function points.

    If weights are in the vector w, then
    derivative is w[0] * f(x-ho*dx) + ... + w[-1] * f(x+h0*dx)

    Parameters
    ----------
    Np : int
        Number of points for the central derivative.
    ndiv : int, optional
        Number of divisions.  Default is 1.

    Notes
    -----
    Can be inaccurate for large number of points.

    """
    if Np < ndiv + 1:
        raise ValueError("Number of points must be at least the derivative order + 1.")
    if Np % 2 == 0:
        raise ValueError("The number of points must be odd.")
    from scipy import linalg
    ho = Np >> 1
    x = arange(-ho,ho+1.0)
    x = x[:,newaxis]
    X = x**0.0
    for k in range(1,Np):
        X = hstack([X,x**k])
    w = product(arange(1,ndiv+1),axis=0)*linalg.inv(X)[ndiv]
    return w


def derivative(func, x0, dx=1.0, n=1, args=(), order=3):
    """
    Find the n-th derivative of a function at a point.

    Given a function, use a central difference formula with spacing `dx` to
    compute the `n`-th derivative at `x0`.

    Parameters
    ----------
    func : function
        Input function.
    x0 : float
        The point at which `n`-th derivative is found.
    dx : float, optional
        Spacing.
    n : int, optional
        Order of the derivative. Default is 1.
    args : tuple, optional
        Arguments
    order : int, optional
        Number of points to use, must be odd.

    Notes
    -----
    Decreasing the step size too small can result in round-off error.

    Examples
    --------
    >>> from scipy.misc import derivative
    >>> def f(x):
    ...     return x**3 + x**2
    >>> derivative(f, 1.0, dx=1e-6)
    4.9999999999217337

    """
    if order < n + 1:
        raise ValueError("'order' (the number of points used to compute the derivative), "
                         "must be at least the derivative order 'n' + 1.")
    if order % 2 == 0:
        raise ValueError("'order' (the number of points used to compute the derivative) "
                         "must be odd.")
    # pre-computed for n=1 and 2 and low-order for speed.
    if n == 1:
        if order == 3:
            weights = array([-1,0,1])/2.0
        elif order == 5:
            weights = array([1,-8,0,8,-1])/12.0
        elif order == 7:
            weights = array([-1,9,-45,0,45,-9,1])/60.0
        elif order == 9:
            weights = array([3,-32,168,-672,0,672,-168,32,-3])/840.0
        else:
            weights = central_diff_weights(order,1)
    elif n == 2:
        if order == 3:
            weights = array([1,-2.0,1])
        elif order == 5:
            weights = array([-1,16,-30,16,-1])/12.0
        elif order == 7:
            weights = array([2,-27,270,-490,270,-27,2])/180.0
        elif order == 9:
            weights = array([-9,128,-1008,8064,-14350,8064,-1008,128,-9])/5040.0
        else:
            weights = central_diff_weights(order,2)
    else:
        weights = central_diff_weights(order, n)
    val = 0.0
    ho = order >> 1
    for k in range(order):
        val += weights[k]*func(x0+(k-ho)*dx,*args)
    return val / product((dx,)*n,axis=0)


def pade(an, m):
    """
    Return Pade approximation to a polynomial as the ratio of two polynomials.

    Parameters
    ----------
    an : (N,) array_like
        Taylor series coefficients.
    m : int
        The order of the returned approximating polynomials.

    Returns
    -------
    p, q : Polynomial class
        The pade approximation of the polynomial defined by `an` is
        `p(x)/q(x)`.

    Examples
    --------
    >>> from scipy import misc
    >>> e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0]
    >>> p, q = misc.pade(e_exp, 2)

    >>> e_exp.reverse()
    >>> e_poly = np.poly1d(e_exp)

    Compare ``e_poly(x)`` and the pade approximation ``p(x)/q(x)``

    >>> e_poly(1)
    2.7166666666666668

    >>> p(1)/q(1)
    2.7179487179487181

    """
    from scipy import linalg
    an = asarray(an)
    N = len(an) - 1
    n = N - m
    if n < 0:
        raise ValueError("Order of q <m> must be smaller than len(an)-1.")
    Akj = eye(N+1, n+1)
    Bkj = zeros((N+1, m), 'd')
    for row in range(1, m+1):
        Bkj[row,:row] = -(an[:row])[::-1]
    for row in range(m+1, N+1):
        Bkj[row,:] = -(an[row-m:row])[::-1]
    C = hstack((Akj, Bkj))
    pq = linalg.solve(C, an)
    p = pq[:n+1]
    q = r_[1.0, pq[n+1:]]
    return poly1d(p[::-1]), poly1d(q[::-1])


def lena():
    """
    Function that previously returned an example image

    .. note:: Removed in 0.17

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        This functionality has been removed due to licensing reasons.

    Notes
    -----
    The image previously returned by this function has an incompatible license
    and has been removed from SciPy. Please use `face` or `ascent` instead.

    See Also
    --------
    face, ascent
    """
    raise RuntimeError('lena() is no longer included in SciPy, please use '
                       'ascent() or face() instead')


def ascent():
    """
    Get an 8-bit grayscale bit-depth, 512 x 512 derived image for easy use in demos

    The image is derived from accent-to-the-top.jpg at
    http://www.public-domain-image.com/people-public-domain-images-pictures/

    Parameters
    ----------
    None

    Returns
    -------
    ascent : ndarray
       convenient image to use for testing and demonstration

    Examples
    --------
    >>> import scipy.misc
    >>> ascent = scipy.misc.ascent()
    >>> ascent.shape
    (512, 512)
    >>> ascent.max()
    255

    >>> import matplotlib.pyplot as plt
    >>> plt.gray()
    >>> plt.imshow(ascent)
    >>> plt.show()

    """
    import pickle
    import os
    fname = os.path.join(os.path.dirname(__file__),'ascent.dat')
    with open(fname, 'rb') as f:
        ascent = array(pickle.load(f))
    return ascent


def face(gray=False):
    """
    Get a 1024 x 768, color image of a raccoon face.

    raccoon-procyon-lotor.jpg at http://www.public-domain-image.com

    Parameters
    ----------
    gray : bool, optional
        If True return 8-bit grey-scale image, otherwise return a color image

    Returns
    -------
    face : ndarray
        image of a racoon face

    Examples
    --------
    >>> import scipy.misc
    >>> face = scipy.misc.face()
    >>> face.shape
    (768, 1024, 3)
    >>> face.max()
    255
    >>> face.dtype
    dtype('uint8')

    >>> import matplotlib.pyplot as plt
    >>> plt.gray()
    >>> plt.imshow(face)
    >>> plt.show()

    """
    import bz2
    import os
    with open(os.path.join(os.path.dirname(__file__), 'face.dat'), 'rb') as f:
        rawdata = f.read()
    data = bz2.decompress(rawdata)
    face = fromstring(data, dtype='uint8')
    face.shape = (768, 1024, 3)
    if gray is True:
        face = (0.21 * face[:,:,0] + 0.71 * face[:,:,1] + 0.07 * face[:,:,2]).astype('uint8')
    return face


def inversefunc(f,vmin=None,vmax=None,vminopen=False,vmaxopen=False,accuracy=2):
    r"""Obtain the inverse of a function.
    
    Returns a callable that calculates the numerical inverse of the function
    `f`. In order for the enumerical inverse to exist in an interval, the 
    input function must have monotonic behavior i.e. be a purely decreasing 
    or purely increasing in that interval. By default the interval spans all
    the real numbers, howevere it can be restricted with the `vmin` and `vmax`
    arguments to [`vmin`, `vmax`],  (`vmin`, `vmax`],  [`vmin`, `vmax`) or  
    (`vmin`, `vmax`) depending on `vminopen` and `vmaxopen`.
    
    Parameters
    ----------
    f : callable
        Callable representing the function to be inverted, able to take a
        ndarray or an scalar and return an object of the same kind with the 
        evaluation of the function. The function must not diverge and have a 
        monotonic behavior in the chosen interval.
    vmin : float, optional
        Lower side of the interval. Default None (-Inf).
    vminopen : bool, optional
        Whether the interval has an open end at `vmin`. Default False.  
    vmax : float, optional
        Higher side of the interval. It must be strictly larger than `vmin`. 
        Default None (-Inf).
    vminopen : bool, optional
        Whether the interval has an open end at `vmax`. Default False.   
    
    Returns
    -------
    callable
        Inverse function of `f`. It can take scalars or ndarrays, and return
        objects of the same kind with the calculated inverse values.

    Examples
    --------
    >>> import scipy.misc
    >>> import numpy as np
    >>> cube = lambda x: x**3
    >>> invcube = scipy.misc.inversefunc(cube)
    >>> invcube(27) # Should give 3
    array(3.0000000063797567)
    >>> square = lambda x: x**2
    >>> invsquare = scipy.misc.inversefunc(square, vmin=0, accuracy=4)
    >>> invsquare([4,16,64]) # Should give [2, 4, 8]
    array([ 2.        ,  4.00096317,  8.00028687])
    >>> log = lambda x: np.log10(x)
    >>> invlog = scipy.misc.inversefunc(log, vmin=0, vminopen=True)
    >>> invlog(-2.) # Should give 0.01
    array(0.010001898156620113)
    >>> cos = lambda x: np.cos(x)
    >>> invcos = scipy.misc.inversefunc(cos, vmin=0, vmax=np.pi)
    >>> invcos([1,0,-1]) # Should give [0, pi/2, pi]
    array([  4.31643739e-06,   1.57079633e+00,   3.14158845e+00])
    >>> tan = lambda x: np.tan(x)
    >>> invtan = scipy.misc.inversefunc(tan, 
                                    vmin=-np.pi/2, 
                                    vmax=np.pi/2, 
                                    vminopen=True, 
                                    vmaxopen=True)
    >>> invtan([1,0,-1]) # Should give [pi/4, 0, -pi/4]
    array([ 0.78539955,  0.        , -0.78539955])
    
    """

    if vmin is not None:
        vmin=float(vmin)
    if vmax is not None:
        vmax=float(vmax)
    
    if vmin is not None and vmax is not None:
        if vmin >= vmax:
            raise ValueError("vmin must be less than vmax")
    
    min_kwargs={}
    constraint=None
    
    if vmin is not None and vmax is not None:
        min_kwargs['method']='bounded'
        min_kwargs['bounds']=(vmin,vmax)
        if vminopen:
            constraint=lambda x: (x-vmin>0)
        if vmaxopen:
            constraint=lambda x: (vmax-x>0)
        if vminopen and vminopen:
            constraint=lambda x: (vmax-x>0)*(x-vmin>0)
    elif vmin is not None:
        if vminopen:
            constraint=lambda x: x-vmin>0
            min_kwargs['bracket']=(vmin+1,vmin+2)
        else:
            constraint=lambda x: x-vmin>=0
            min_kwargs['bracket']=(vmin,vmin+2)
        min_kwargs['tol']=10.**(-accuracy-1)
    elif vmax is not None:
        if vmaxopen:
            constraint=lambda x: vmax-x>0
            min_kwargs['bracket']=(vmax-2,vmax-1)
        else:
            constraint=lambda x: vmax-x>=0
            min_kwargs['bracket']=(vmax-2,vmax)
        min_kwargs['tol']=10.**(-accuracy-1)

    def inv(xin):
        xin=np.asarray(xin, dtype=np.float64)
        shapein=xin.shape
        xin=xin.flatten()
        results=xin.copy()*np.nan
        resultsmask=np.zeros(xin.shape,dtype=np.bool)
        
        for j in range(xin.size):
            if constraint:    
                optimizer=lambda x, j=j,f=f: ((((f(x)-xin[j]))**2).sum() 
                                              if constraint(x) else np.inf)
            else:
                optimizer=lambda x, j=j,f=f: (((f(x)-xin[j]))**2).sum()
            result=minimize_scalar(optimizer,**min_kwargs)
            try:
                results[j]=result.x
                resultsmask[j]=result.success
            except:
                resultsmask[j]=False
        if any(resultsmask==False):
            warnings.warn("Trouble calculating inverse for values: "
                          "%s" % str(xin[~resultsmask]),RuntimeWarning)
        
        try:
            np.testing.assert_array_almost_equal(xin,f(results), decimal=accuracy)
        except AssertionError:        
            warnings.warn("Results obtained with less than %g "  
                          "decimal digits of accuracy"
                          %accuracy,RuntimeWarning)
            
        return results.reshape(shapein)
    return inv
