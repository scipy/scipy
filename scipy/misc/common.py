"""
Functions which are common and require SciPy Base and Level 1 SciPy
(special, linalg)
"""

from numpy import exp, log, asarray, arange, newaxis, hstack, product, array, \
                  where, zeros, extract, place, pi, sqrt, eye, poly1d, dot, r_

__all__ = ['logsumexp', 'factorial','factorial2','factorialk','comb',
           'central_diff_weights', 'derivative', 'pade', 'lena']

# XXX: the factorial functions could move to scipy.special, and the others
# to numpy perhaps?

def logsumexp(a):
    """Compute the log of the sum of exponentials of input elements.

    Parameters
    ----------
    a : array_like
        Input array.

    Returns
    -------
    res : ndarray
        The result, ``np.log(np.sum(np.exp(a)))`` calculated in a numerically
        more stable way.

    See Also
    --------
    numpy.logaddexp, numpy.logaddexp2

    Notes
    -----
    Numpy has a logaddexp function which is very similar to `logsumexp`.

    """
    a = asarray(a)
    a_max = a.max()
    return a_max + log((exp(a-a_max)).sum())

def factorial(n,exact=0):
    """
    The factorial function, n! = special.gamma(n+1).

    If exact is 0, then floating point precision is used, otherwise
    exact long integer is computed.

    - Array argument accepted only for exact=0 case.
    - If n<0, the return value is 0.

    Parameters
    ----------
    n : int or array_like of ints
        Calculate ``n!``.  Arrays are only supported with `exact` set
        to False.  If ``n < 0``, the return value is 0.
    exact : bool, optional
        The result can be approximated rapidly using the gamma-formula
        above.  If `exact` is set to True, calculate the
        answer exactly using integer arithmetic. Default is False.

    Returns
    -------
    nf : float or int
        Factorial of `n`, as an integer or a float depending on `exact`.

    Examples
    --------
    >>> arr = np.array([3,4,5])
    >>> sc.factorial(arr, exact=False)
    array([   6.,   24.,  120.])
    >>> sc.factorial(5, exact=True)
    120L

    """
    if exact:
        if n < 0:
            return 0L
        val = 1L
        for k in xrange(1,n+1):
            val *= k
        return val
    else:
        from scipy import special
        n = asarray(n)
        sv = special.errprint(0)
        vals = special.gamma(n+1)
        sv = special.errprint(sv)
        return where(n>=0,vals,0)


def factorial2(n, exact=False):
    """
    Double factorial.

    This is the factorial with every second value skipped, i.e.,
    ``7!! = 7 * 5 * 3 * 1``.  It can be approximated numerically as::

      n!! = special.gamma(n/2+1)*2**((m+1)/2)/sqrt(pi)  n odd
          = 2**(n/2) * (n/2)!                           n even

    Parameters
    ----------
    n : int or array_like
        Calculate ``n!!``.  Arrays are only supported with `exact` set
        to False.  If ``n < 0``, the return value is 0.
    exact : bool, optional
        The result can be approximated rapidly using the gamma-formula
        above (default).  If `exact` is set to True, calculate the
        answer exactly using integer arithmetic.

    Returns
    -------
    nff : float or int
        Double factorial of `n`, as an int or a float depending on
        `exact`.

    Examples
    --------
    >>> factorial2(7, exact=False)
    array(105.00000000000001)
    >>> factorial2(7, exact=True)
    105L

    """
    if exact:
        if n < -1:
            return 0L
        if n <= 0:
            return 1L
        val = 1L
        for k in xrange(n,0,-2):
            val *= k
        return val
    else:
        from scipy import special
        n = asarray(n)
        vals = zeros(n.shape,'d')
        cond1 = (n % 2) & (n >= -1)
        cond2 = (1-(n % 2)) & (n >= -1)
        oddn = extract(cond1,n)
        evenn = extract(cond2,n)
        nd2o = oddn / 2.0
        nd2e = evenn / 2.0
        place(vals,cond1,special.gamma(nd2o+1)/sqrt(pi)*pow(2.0,nd2o+0.5))
        place(vals,cond2,special.gamma(nd2e+1) * pow(2.0,nd2e))
        return vals

def factorialk(n,k,exact=1):
    """
    n(!!...!)  = multifactorial of order k
    k times


    Parameters
    ----------
    n : int, array_like
        Calculate multifactorial. Arrays are only supported with exact
        set to False. If n < 0, the return value is 0.
    exact : bool, optional
        If exact is set to True, calculate the answer exactly using
        integer arithmetic.

    Returns
    -------
    val : int
        Multi factorial of n.

    Raises
    ------
    NotImplementedError
        Raises when exact is False

    Examples
    --------
    >>> sc.factorialk(5, 1, exact=True)
    120L
    >>> sc.factorialk(5, 3, exact=True)
    10L

    """
    if exact:
        if n < 1-k:
            return 0L
        if n<=0:
            return 1L
        val = 1L
        for j in xrange(n,0,-k):
            val = val*j
        return val
    else:
        raise NotImplementedError


def comb(N,k,exact=0):
    """
    The number of combinations of N things taken k at a time.
    This is often expressed as "N choose k".

    Parameters
    ----------
    N : int, array
        Number of things.
    k : int, array
        Number of elements taken.
    exact : int, optional
        If exact is 0, then floating point precision is used, otherwise
        exact long integer is computed.

    Returns
    -------
    val : int, array
        The total number of combinations.

    Notes
    -----
    - Array arguments accepted only for exact=0 case.
    - If k > N, N < 0, or k < 0, then a 0 is returned.

    Examples
    --------
    >>> k = np.array([3, 4])
    >>> n = np.array([10, 10])
    >>> sc.comb(n, k, exact=False)
    array([ 120.,  210.])
    >>> sc.comb(10, 3, exact=True)
    120L

    """
    if exact:
        if (k > N) or (N < 0) or (k < 0):
            return 0L
        val = 1L
        for j in xrange(min(k, N-k)):
            val = (val*(N-j))//(j+1)
        return val
    else:
        from scipy import special
        k,N = asarray(k), asarray(N)
        lgam = special.gammaln
        cond = (k <= N) & (N >= 0) & (k >= 0)
        sv = special.errprint(0)
        vals = exp(lgam(N+1) - lgam(N-k+1) - lgam(k+1))
        sv = special.errprint(sv)
        return where(cond, vals, 0.0)

def central_diff_weights(Np, ndiv=1):
    """
    Return weights for an Np-point central derivative of order ndiv
    assuming equally-spaced function points.

    If weights are in the vector w, then
    derivative is w[0] * f(x-ho*dx) + ... + w[-1] * f(x+h0*dx)

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
    Find the n-th derivative of a function at point x0.

    Given a function, use a central difference formula with spacing `dx` to
    compute the n-th derivative at `x0`.

    Parameters
    ----------
    func : function
        Input function.
    x0 : float
        The point at which nth derivative is found.
    dx : int, optional
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
    >>> def x2(x):
    ...     return x*x
    ...
    >>> derivative(x2, 2)
    4.0

    """
    if order < n + 1:
        raise ValueError("'order' (the number of points used to compute the derivative), "
                         "must be at least the derivative order 'n' + 1.")
    if order % 2 == 0:
        raise ValueError("'order' (the number of points used to compute the derivative) "
                         "must be odd.")
    # pre-computed for n=1 and 2 and low-order for speed.
    if n==1:
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
    elif n==2:
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
    """Given Taylor series coefficients in an, return a Pade approximation to
    the function as the ratio of two polynomials p / q  where the order of q is m.
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
    Get classic image processing example image, Lena, at 8-bit grayscale
    bit-depth, 512 x 512 size.

    Parameters
    ----------
    None

    Returns
    -------
    lena : ndarray
        Lena image

    Examples
    --------
    >>> import scipy.misc
    >>> lena = scipy.misc.lena()
    >>> lena.shape
    (512, 512)
    >>> lena.max()
    245
    >>> lena.dtype
    dtype('int32')

    >>> import matplotlib.pyplot as plt
    >>> plt.gray()
    >>> plt.imshow(lena)
    >>> plt.show()

    """
    import cPickle, os
    fname = os.path.join(os.path.dirname(__file__),'lena.dat')
    f = open(fname,'rb')
    lena = array(cPickle.load(f))
    f.close()
    return lena
