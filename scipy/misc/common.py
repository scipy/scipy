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

__all__ = ['logsumexp', 'central_diff_weights', 'derivative', 'pade', 'lena',
           'ascent', 'face', 'inversefunc']


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


def inversefunc(f,
                domain=None,
                image=None,
                open_domain=None,
                accuracy=2):
    r"""Obtain the inverse of a function.

    Returns a callable that calculates the numerical inverse of the function
    `f`. In order for the numerical inverse to exist in its domain, the
    input function must have, continuous, strictly monotonic behavior i.e. be a
    purely decreasing or purely increasing in that domain. By default the
    domain interval spans all the real numbers, however it can be restricted
    with the `domain` and `open_domain` arguments. The image of the function
    in the interval may be provided, for cases where the function is non
    non continuous right at the end of an open interval.

    Parameters
    ----------
    f : callable
        Callable representing the function to be inverted, able to take a
        ndarray or an scalar and return an object of the same kind with the
        evaluation of the function. The function must not diverge and have a
        monotonic behavior in the chosen interval.
    domain : float, ndarray, optional
        Boundaries of the domain (`domain[0]`, `domain[1]`).
        `domain[1]` must be larger than `domain[0]`.
        None values are assumed to be no boundary in that direction.
        A single scalar value will set it to [`domain`, None].
        Default None (-Inf, Inf).
    open_domain : bool, ndarray, optional
        Whether the domain is an open interval at each of the ends.
        A single scalar boolean will set it to [`open_domain`, `open_domain`].
        Default None [False, False].
    image : float, ndarray, optional
        Image of the function in the domain (`image[0]`, `image[1]`).
        `image[1]` must be larger than `image[0]`.
        None values are assumed to be no boundary in that direction.
        Default None, this is (-Inf, Inf) if domain is None, or the limits
        set by f(domain[0]) and f(domain[1]).
    accuracy : int, optional
        Number of digits for the desired accuracy. It will give a warning
        if the accuracy is worse than this.
        Default 2.

    Returns
    -------
    callable
        Inverse function of `f`. It can take scalars or ndarrays, and return
        objects of the same kind with the calculated inverse values.

    Notes
    -----

    .. versionadded:: 0.19.0

    Examples
    --------
    >>> import scipy.misc
    >>> import numpy as np
    >>> cube = (lambda x: x**3)
    >>> invcube = scipy.misc.inversefunc(cube)
    >>> invcube(27) # Should give 3
    array(3.0000000063797567)
    >>> square = (lambda x: x**2)
    >>> invsquare = scipy.misc.inversefunc(square, domain=0)
    >>> invsquare([4, 16, 64]) # Should give [2, 4, 8]
    array([ 2.,  4.,  8.])
    >>> log = (lambda x: np.log10(x))
    >>> invlog = scipy.misc.inversefunc(log, domain=0, open_domain=True)
    >>> invlog(-2.) # Should give 0.01
    array(0.0099999999882423)
    >>> cos = (lambda x: np.cos(x))
    >>> invcos = scipy.misc.inversefunc(cos, domain=[0, np.pi])
    >>> invcos([1, 0, -1]) # Should give [0, pi / 2, pi]
    array([  5.44203736e-09,   1.57079632e+00,   3.14159265e+00])
    >>> tan = (lambda x: np.tan(x))
    >>> invtan = scipy.misc.inversefunc(tan,
    ...                                 domain=[-np.pi / 2, np.pi / 2],
    ...                                 open_domain=True)
    >>> invtan([1, 0, -1]) # Should give [pi / 4, 0, -pi / 4]
    array([  7.85398163e-01,   1.29246971e-26,  -7.85398163e-01])

    """
    error_domain = ("domain must be a single scalar, or a have two "
                    "elements [xmin, xmax]. Set None, to leave it "
                    "unlimited on one side.")

    # Homogenizing parameters
    if domain is None:
        xmin = None
        xmax = None
    else:
        domain = np.asarray(domain)
        if domain.ndim == 0:
            xmin = float(domain)
            xmax = None
        elif domain.ndim == 1 and domain.size != 2:
            raise ValueError(error_domain)
        elif domain.ndim > 1:
            raise ValueError(error_domain)
        else:
            xmin = (float(domain[0]) if domain[0] is not None else None)
            xmax = (float(domain[1]) if domain[1] is not None else None)

    error_open_domain = ("open_domain must be a single scalar, or a have two "
                         "bool elements [open_xmin, open_xmax].")
    if open_domain is None:
        xmin_open = False
        xmax_open = False
    else:
        open_domain = np.asarray(open_domain)
        if open_domain.ndim == 0:
            xmin_open = bool(open_domain)
            xmax_open = bool(open_domain)
        elif open_domain.ndim == 1 and domain.size != 2:
            raise ValueError(error_open_domain)
        elif open_domain.ndim > 1:
            raise ValueError(error_open_domain)
        else:
            xmin_open = bool(open_domain[0])
            xmax_open = bool(open_domain[1])

    error_image = ("image must be a single scalar, or a have two "
                   "bool elements [ymin, ymax].")
    if image is None:
        ymin = None
        ymax = None
    else:
        image = np.asarray(image)
        if image.ndim != 1 or image.size != 2:
            raise ValueError(error_image)
        else:
            ymin = (float(image[0]) if image[0] is not None else None)
            ymax = (float(image[1]) if image[1] is not None else None)

    if xmin is not None and xmax is not None:
        if xmin >= xmax:
            raise ValueError("domain[0] min must be less than domain[1]")

    if ymin is not None and ymax is not None:
        if ymin >= ymax:
            raise ValueError("image[0] min must be less than image[1]")

    # Calculating if the function is increasing or decreasing, using ref points
    # anywhere in the valid range (Function has to be strictly monotonic)
    if xmin is not None and xmax is not None:
        d = xmax-xmin
        ref1 = xmin + d/4.
        ref2 = xmax - d/4.
    elif xmin is not None:
        ref1 = xmin+1.
        ref2 = xmin+2.
    elif xmax is not None:
        ref1 = xmax-2.
        ref2 = xmax-1.
    else:
        ref1 = 0.
        ref2 = 1.
    trend = np.sign(f(ref2)-f(ref1))

    if trend == 0:
        raise ValueError("Function is not strictly monotonic")

    # Calculating the image by default
    if ymin is None and ((xmin is not None and trend == 1) or
                         (xmax is not None and trend == -1)):
        try:
            with warnings.catch_warnings(record=True):
                ymin = f(xmin) if trend == 1 else f(xmax)
        except:
            raise ValueError("Cannot automatically calculate the lower limit "
                             "of the image please inclue it as a parameter")
    if ymax is None and ((xmax is not None and trend == 1) or
                         (xmin is not None and trend == -1)):
        try:
            with warnings.catch_warnings(record=True):
                ymax = f(xmax) if trend == 1 else f(xmin)
        except:
            raise ValueError("Cannot automatically calculate the upper limit "
                             "of the image please include it as a parameter")

    # Creating bounded function
    def bounded_f(x):
        if xmin is not None and (x < xmin or (x == xmin and xmin_open)):
                val = -1*np.inf*trend
        elif xmax is not None and (x > xmax or (x == xmax and xmax_open)):
                val = np.inf*trend
        else:
            val = f(x)
        return val

    min_kwargs = {}
    min_kwargs['bracket'] = (ref1, ref2)
    min_kwargs['tol'] = 1.48e-08
    min_kwargs['method'] = 'Brent'

    def inv(yin):
        yin = np.asarray(yin, dtype=np.float64)
        shapein = yin.shape
        yin = yin.flatten()

        if ymin is not None:
            if xmin_open:
                mask = yin <= ymin
            else:
                mask = yin < ymin
            if yin[mask].size > 0:
                ValueError("Requested values %s lower than the lower limit"
                           "%g of the image" % (yin[mask], ymin))
        if ymax is not None:
            if xmax_open:
                mask = yin >= ymax
            else:
                mask = yin > ymax
            if yin[mask].size > 0:
                ValueError("Requested values %s higher than the higher limit"
                           "%g of the image" % (yin[mask], ymax))

        results = yin.copy() * np.nan
        resultsmask = np.zeros(yin.shape, dtype=np.bool)

        for j in range(yin.size):
            optimizer = (lambda x, j=j, f=f: (((bounded_f(x) - yin[j]))**2))
            try:
                with warnings.catch_warnings(record=True):
                    result = minimize_scalar(optimizer, **min_kwargs)
                results[j] = result.x
                resultsmask[j] = result.success
            except:
                resultsmask[j] = False
        if any(~resultsmask):
            warnings.warn("Trouble calculating inverse for values: "
                          "%s" % str(yin[~resultsmask]), RuntimeWarning)

        try:
            np.testing.assert_array_almost_equal(yin, f(results),
                                                 decimal=accuracy)
        except AssertionError:
            warnings.warn("Results obtained with less than %g "
                          "decimal digits of accuracy"
                          % accuracy, RuntimeWarning)

        return results.reshape(shapein)
    return inv
