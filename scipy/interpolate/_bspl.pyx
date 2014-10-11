"""
Routines for evaluating and manipulating B-splines.

"""

import numpy as np
cimport numpy as cnp

cimport cython

cdef extern from "src/__fitpack.h":
    void _deBoor_D(double *t, double x, int k, int ell, int m, double *result)

cdef double nan = np.nan

ctypedef double complex double_complex

ctypedef fused double_or_complex:
    double
    double complex


#------------------------------------------------------------------------------
# B-splines
#------------------------------------------------------------------------------

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int find_interval(double[::1] t,
                       int k,
                       double xval, 
                       int prev_l,
                       int extrapolate) nogil:
    """
    Find an interval such that t[interval] <= xval < t[interval+1].

    Uses a linear search with locality, see fitpack's splev.

    Parameters
    ----------
    t : ndarray, shape (nt,)
        Knots
    k : int
        B-spline degree
    xval : double
        value to find the inteval for
    prev_l : int
        interval where the previous value was located.
        if unknown, use any value < k to start the search.
    extrapolate : int
        whether to return the last or the first interval if xval
        is out of bounds.

    Returns
    -------
    interval : int
        Suitable interval or -1 if xval was nan.

    """
    cdef:
        int l
        int n = t.shape[0] - k - 1
        double tb = t[k]
        double te = t[n]

    if xval != xval:
        # nan
        return -1

    if (xval < tb) or (xval > te):
        if extrapolate:
            pass
        else:
            return -1

    l = prev_l if k < prev_l < n else k

    # xval is in support, search for interval s.t. t[interval] <= xval < t[l+1]
    while(xval < t[l] and l != k):
        l -= 1

    l += 1
    while(xval >= t[l] and l != n):
        l += 1

    return l-1


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def evaluate_spline(double[::1] t,
             double_or_complex[:,::1] c,
             int k,
             double[::1] xp,
             int der,
             int extrapolate,
             double_or_complex[:,::1] out):
    """
    Evaluate a spline in the B-spline basis.

    Parameters
    ----------
    t : ndarray, shape (n+k+1)
        knots
    c : ndarray, shape (n, m)
        B-spline coefficients
    xp : ndarray, shape (s,)
        Points to evaluate the spline at.
    der : int
        Order of derivative to evaluate.
    extrapolate : int, optional
        Whether to extrapolate to ouf-of-bounds points, or to return NaNs.
    out : ndarray, shape (s, m)
        Computed values of the spline at each of the input points.
        This argument is modified in-place.

    """

    cdef int ip, jp, n, a
    cdef int i, interval
    cdef double xval

    # shape checks
    if out.shape[0] != xp.shape[0]:
        raise ValueError("out and xp have incompatible shapes")
    if out.shape[1] != c.shape[1]:
        raise ValueError("out and c have incompatible shapes")

    # check derivative order
    if der < 0:
        raise NotImplementedError("Cannot do derivative order %s." % der)

    n = c.shape[0]
    cdef double[::1] work = np.empty(2*k+2, dtype=np.float_)

    # evaluate
    interval = k
    for ip in range(xp.shape[0]):
        xval = xp[ip]

        # Find correct interval
        interval = find_interval(t, k, xval, interval, extrapolate)

        if interval < 0:
            # xval was nan etc
            for jp in range(c.shape[1]):
                out[ip, jp] = nan
            continue

        # Evaluate (k+1) b-splines which are non-zero on the interval.
        # on return, first k+1 elemets of work are B_{m-k},..., B_{m}
        _deBoor_D(&t[0], xval, k, interval, der, &work[0])

        # Form linear combinations
        for jp in range(c.shape[1]):
            out[ip, jp] = 0.
            for a in range(k+1):
                out[ip, jp] = out[ip, jp] + c[interval + a - k, jp] * work[a]


def evaluate_all_bspl(double[::1] t, int k, double xval, int m, int nu=0):
    """Evaluate the ``k+1`` B-splines which are non-zero on interval ``m``.

    Parameters
    ----------
    t : ndarray, shape (nt + k + 1,)
        sorted 1D array of knots
    k : int
        spline order
    xval: float
        argument at which to evaluate the B-splines
    m : int
        index of the left edge of the evaluation interval, ``t[m] <= x < t[m+1]``
    nu : int, optional
        Evaluate derivatives order `nu`. Default is zero.

    Returns
    -------
    ndarray, shape (k+1,)
        The values of B-splines :math:`[B_{m-k}(xval), ..., B_{m}(xval)]` if
        `nu` is zero, otherwise the derivatives of order `nu`.

    Examples
    --------

    A textbook use of this sort of routine is plotting the ``k+1`` polynomial
    pieces which make up a B-spline of order `k`.

    Consider a cubic spline

    >>> k = 3
    >>> t = [0., 2., 2., 3., 4.]   # internal knots
    >>> a, b = t[0], t[-1]    # base interval is [a, b)
    >>> t = [a]*k + t + [b]*k  # add boundary knots

    >>> import matplotlib.pyplot as plt
    >>> xx = np.linspace(a, b, 100)
    >>> plt.plot(xx, BSpline.basis_element(k:-k)(xx),
    ...          'r-', lw=5, alpha=0.5)
    >>> c = ['b', 'g', 'c', 'k']

    Now we use slide an interval ``t[m]..t[m+1]`` along the base interval 
    ``a..b`` and use `evaluate_all_bspl` to compute the restriction of
    the B-spline of interest to this interval:

    >>> for i in range(k+1):
    ...    x1, x2 = t[2*k - i], t[2*k - i + 1]
    ...    xx = np.linspace(x1 - 0.5, x2 + 0.5)
    ...    yy = [evaluate_all_bspl(t, k, x, 2*k - i)[i] for x in xx]
    ...
    ...    plt.plot(xx, yy, c[i] + '--', lw=3, label=str(i))
    ...
    ... plt.grid(True)
    ... plt.legend()
    ... plt.show()

    """
    bbb = np.empty(2*k+2, dtype=np.float_)
    cdef double[::1] work = bbb
    _deBoor_D(&t[0], xval, k, m, nu, &work[0])
    return bbb[:k+1]


@cython.wraparound(False)
@cython.boundscheck(False)
def _colloc(double[::1] x, double[::1] t, int k, int offset=0):
    """Build the B-spline collocation matrix.  

    :math:`B_{j,l} = B_l(x_j)`, so that row ``j`` contains
    all the B-splines which are non-zero at ``x_j``.

    The matrix is constructed in the banded storage, ready for
    consumption by ``scipy.linalg.solve_banded``.
    This routine is not supposed to be called directly, and
    does no error checking.

    Parameters
    x : ndarray, shape (n,)
        sorted 1D array of x values
    t : ndarray, shape (nt + k + 1,)
        sorted 1D array of knots
    k : int
        spline order
    offset : int, optional
        skip this many rows

    Returns
    -------
    Ab : ndarray, shape ((ku+kl+1), nt)
        B-spline collocation matrix in the band storage with
        ``ku`` upper diagonals and ``kl`` lower diagonals.
    (kl, ku) : (int, int)
        the number of lower and upper diagonals
        
    """
    cdef int nt = t.shape[0] - k -1
    cdef int left, j, a, kl, ku, clmn
    cdef double xval

    kl = ku = k
    cdef cnp.ndarray[cnp.float_t, ndim=2] Ab = np.zeros((kl + ku + 1, nt),
            dtype=np.float_)
    cdef cnp.ndarray[cnp.float_t, ndim=1] out = np.empty(2*k+2, dtype=np.float_)
    
    # collocation matrix
    left = k
    for j in range(x.shape[0]):
        xval = x[j]
        # find interval
        left = find_interval(t, k, xval, left, extrapolate=False)

        # fill a row
        _deBoor_D(&t[0], xval, k, left, 0, &out[0])
        # for a full matrix it would be ``A[j + offset, left-k:left+1] = bb``
        # in the banded storage, need to spread the row over
        for a in range(k+1):
            clmn = left - k + a
            Ab[ku + j + offset - clmn, clmn] = out[a]
    return Ab, (kl, ku)


def _handle_lhs_derivatives(double[::1]t, int k, double xval, ab, 
                            kl_ku, deriv_ords, int offset=0):
    """ Fill in the entries of the collocation matrix corresponding to known
    derivatives at xval.
    
    The collocation matrix is in the banded storage, as prepared by
    _colloc and ready for consumption by scipy.linalg.solve_banded.
    No error checking.
    
    Parameters
    ----------
    t : ndarray, shape (nt + k + 1,)
        knots
    k : integer
        B-spline order
    xval : float
        The value at which to evaluate the derivatives at.
    ab : ndarray, shape(kl + ku + 1, nt)
        B-spline collocation matrix.
        This argument is modified *in-place*.
    kl_ku : (integer, integer)
        Number of lower and upper diagonals of ab.
    deriv_ords : 1D ndarray
        Orders of derivatives known at xval
    offset : integer, optional
        Skip this many rows of the matrix ab.
    
    """
    cdef int kl, ku, left, nu, a, clmn, row

    kl, ku = kl_ku
    cdef double[::1] out = np.empty(2*k+2, dtype=np.float_)


    # derivatives @ xval
    left = find_interval(t, k, xval, k, extrapolate=False)
    for row in range(deriv_ords.size):
        nu = deriv_ords[row]
        _deBoor_D(&t[0], xval, k, left, nu, &out[0])
        # if A were a full matrix, it would be just
        # ``A[row + offset, left-k:left+1] = bb``.
        for a in range(k+1):
            clmn = left - k + a
            ab[ku + offset + row - clmn, clmn] = out[a]

