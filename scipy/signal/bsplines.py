from numpy import logical_and, asarray, pi, piecewise, array
from numpy.core.umath import sqrt, exp, less, less_equal, greater_equal

# From splinemodule.c
from .spline import cspline2d, sepfir2d

from scipy.special import comb
from scipy._lib._util import float_factorial

__all__ = ['spline_filter', 'bspline', 'gauss_spline']


def spline_filter(Iin, lmbda=5.0):
    """Smoothing spline (cubic) filtering of a rank-2 array.

    Filter an input data set, `Iin`, using a (cubic) smoothing spline of
    fall-off `lmbda`.

    Parameters
    ----------
    Iin : array_like
        input data set
    lmbda : float, optional
        spline smooghing fall-off value, default is `5.0`.

    Returns
    -------
    res : ndarray
        filterd input data

    Examples
    --------
    We can filter an multi dimentional signal (ex: 2D image) using cubic
    B-spline filter:

    >>> from scipy.signal import spline_filter
    >>> import matplotlib.pyplot as plt
    >>> orig_img = np.eye(20)  # create an image
    >>> orig_img[10, :] = 1.0
    >>> sp_filter = spline_filter(orig_img, lmbda=0.1)
    >>> f, ax = plt.subplots(1, 2, sharex=True)
    >>> for ind, data in enumerate([[orig_img, "original image"],
    ...                             [sp_filter, "spline filter"]]):
    ...     ax[ind].imshow(data[0], cmap='gray_r')
    ...     ax[ind].set_title(data[1])
    >>> plt.tight_layout()
    >>> plt.show()

    """
    intype = Iin.dtype.char
    hcol = array([1.0, 4.0, 1.0], 'f') / 6.0
    if intype in ['F', 'D']:
        Iin = Iin.astype('F')
        ckr = cspline2d(Iin.real, lmbda)
        cki = cspline2d(Iin.imag, lmbda)
        outr = sepfir2d(ckr, hcol, hcol)
        outi = sepfir2d(cki, hcol, hcol)
        out = (outr + 1j * outi).astype(intype)
    elif intype in ['f', 'd']:
        ckr = cspline2d(Iin, lmbda)
        out = sepfir2d(ckr, hcol, hcol)
        out = out.astype(intype)
    else:
        raise TypeError("Invalid data type for Iin")
    return out


_splinefunc_cache = {}


def _bspline_piecefunctions(order):
    """Returns the function defined over the left-side pieces for a bspline of
    a given order.

    The 0th piece is the first one less than 0. The last piece is a function
    identical to 0 (returned as the constant 0). (There are order//2 + 2 total
    pieces).

    Also returns the condition functions that when evaluated return boolean
    arrays for use with `numpy.piecewise`.
    """
    try:
        return _splinefunc_cache[order]
    except KeyError:
        pass

    def condfuncgen(num, val1, val2):
        if num == 0:
            return lambda x: logical_and(less_equal(x, val1),
                                         greater_equal(x, val2))
        elif num == 2:
            return lambda x: less_equal(x, val2)
        else:
            return lambda x: logical_and(less(x, val1),
                                         greater_equal(x, val2))

    last = order // 2 + 2
    if order % 2:
        startbound = -1.0
    else:
        startbound = -0.5
    condfuncs = [condfuncgen(0, 0, startbound)]
    bound = startbound
    for num in range(1, last - 1):
        condfuncs.append(condfuncgen(1, bound, bound - 1))
        bound = bound - 1
    condfuncs.append(condfuncgen(2, 0, -(order + 1) / 2.0))

    # final value of bound is used in piecefuncgen below

    # the functions to evaluate are taken from the left-hand side
    #  in the general expression derived from the central difference
    #  operator (because they involve fewer terms).

    fval = float_factorial(order)

    def piecefuncgen(num):
        Mk = order // 2 - num
        if (Mk < 0):
            return 0  # final function is 0
        coeffs = [(1 - 2 * (k % 2)) * float(comb(order + 1, k, exact=1)) / fval
                  for k in range(Mk + 1)]
        shifts = [-bound - k for k in range(Mk + 1)]

        def thefunc(x):
            res = 0.0
            for k in range(Mk + 1):
                res += coeffs[k] * (x + shifts[k]) ** order
            return res
        return thefunc

    funclist = [piecefuncgen(k) for k in range(last)]

    _splinefunc_cache[order] = (funclist, condfuncs)

    return funclist, condfuncs


def bspline(x, n):
    """B-spline basis function of order n.

    Parameters
    ----------
    x : array_like
        a knot vector
    n : int
        The order of the spline. Must be non-negative, i.e., n >= 0

    Returns
    -------
    res : ndarray
        B-spline basis function values

    See Also
    --------
    cubic : A cubic B-spline.
    quadratic : A quadratic B-spline.

    Notes
    -----
    Uses numpy.piecewise and automatic function-generator.

    Examples
    --------
    We can calculate B-Spline basis function of several orders:

    >>> from scipy.signal import bspline, cubic, quadratic
    >>> bspline(0.0, 1)
    1

    >>> knots = [-1.0, 0.0, -1.0]
    >>> bspline(knots, 2)
    array([0.125, 0.75, 0.125])

    >>> np.array_equal(bspline(knots, 2), quadratic(knots))
    True

    >>> np.array_equal(bspline(knots, 3), cubic(knots))
    True

    """
    ax = -abs(asarray(x))
    # number of pieces on the left-side is (n+1)/2
    funclist, condfuncs = _bspline_piecefunctions(n)
    condlist = [func(ax) for func in condfuncs]
    return piecewise(ax, condlist, funclist)


def gauss_spline(x, n):
    r"""Gaussian approximation to B-spline basis function of order n.

    Parameters
    ----------
    x : array_like
        a knot vector
    n : int
        The order of the spline. Must be non-negative, i.e., n >= 0

    Returns
    -------
    res : ndarray
        B-spline basis function values approximated by a zero-mean Gaussian
        function.

    Notes
    -----
    The B-spline basis function can be approximated well by a zero-mean
    Gaussian function with standard-deviation equal to :math:`\sigma=(n+1)/12`
    for large `n` :

    .. math::  \frac{1}{\sqrt {2\pi\sigma^2}}exp(-\frac{x^2}{2\sigma})

    References
    ----------
    .. [1] Bouma H., Vilanova A., Bescos J.O., ter Haar Romeny B.M., Gerritsen
       F.A. (2007) Fast and Accurate Gaussian Derivatives Based on B-Splines. In:
       Sgallari F., Murli A., Paragios N. (eds) Scale Space and Variational
       Methods in Computer Vision. SSVM 2007. Lecture Notes in Computer
       Science, vol 4485. Springer, Berlin, Heidelberg
    .. [2] http://folk.uio.no/inf3330/scripting/doc/python/SciPy/tutorial/old/node24.html

    Examples
    --------
    We can calculate B-Spline basis functions approximated by a gaussian
    distribution:

    >>> from scipy.signal import gauss_spline, bspline
    >>> knots = np.array([-1.0, 0.0, -1.0])
    >>> gauss_spline(knots, 3)
    array([0.15418033, 0.6909883, 0.15418033])  # may vary

    >>> bspline(knots, 3)
    array([0.16666667, 0.66666667, 0.16666667])  # may vary

    """
    x = asarray(x)
    signsq = (n + 1) / 12.0
    return 1 / sqrt(2 * pi * signsq) * exp(-x ** 2 / 2 / signsq)
