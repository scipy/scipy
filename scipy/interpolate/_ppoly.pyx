"""
Routines for evaluating and manipulating piecewise polynomials in
local power basis.

"""

from .polyint import _Interpolator1D
import numpy as np

cimport cython

cdef double nan = np.nan

ctypedef double complex double_complex

ctypedef fused double_or_complex:
    double
    double complex


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def evaluate(double_or_complex[:,:,::1] c,
             double[::1] x,
             double[::1] xp,
             int dx,
             double_or_complex[:,::1] out):
    """
    Evaluate a piecewise polynomial.

    Parameters
    ----------
    c : ndarray, shape (k, m, n)
        Coefficients local polynomials of order `k-1` in `m` intervals.
        There are `n` polynomials in each interval.
        Coefficient of highest order-term comes first.
    x : ndarray, shape (m+1,)
        Breakpoints of polynomials
    xp : ndarray, shape (r,)
        Points to evaluate the piecewise polynomial at.
    dx : int
        Order of derivative to evaluate.  The derivative is evaluated
        piecewise and may have discontinuities.

    Returns
    -------
    out : ndarray, shape (r, n)
        Value of each polynomial at each of the input points.
        For points outside the span ``x[0] ... x[-1]``,
        ``nan`` is returned.

    """

    cdef int ip, jp, kp, i
    cdef int interval, high, low, mid
    cdef int has_out_of_bounds
    cdef double_or_complex z, res
    cdef double prefactor, xval

    # check derivative order
    if dx < 0:
        raise ValueError("Order of derivative cannot be negative")

    # shape checks
    if out.shape[0] != xp.shape[0]:
        raise ValueError("out and xp have incompatible shapes")
    if out.shape[1] != c.shape[2]:
        raise ValueError("out and c have incompatible shapes")
    if c.shape[1] != x.shape[0] - 1:
        raise ValueError("x and c have incompatible shapes")

    # evaluate
    interval = 0
    has_out_of_bounds = 0

    for ip in range(len(xp)):
        xval = xp[ip]

        # Find correct interval
        i = find_interval(x, xval, interval)
        if i < 0:
            has_out_of_bounds = 1
            for jp in range(c.shape[2]):
                out[ip, jp] = nan
            continue
        interval = i

        # Evaluate the local polynomial(s)
        for jp in range(c.shape[2]):
            out[ip, jp] = evaluate_poly1(xval - x[interval], c, interval, jp, dx)

    return has_out_of_bounds


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def fix_continuity(double_or_complex[:,:,::1] c,
                   double[::1] x,
                   int order):
    """
    Make a piecewise polynomial continuously differentiable to given order.

    Parameters
    ----------
    c : ndarray, shape (k, m, n)
        Coefficients local polynomials of order `k-1` in `m` intervals.
        There are `n` polynomials in each interval.
        Coefficient of highest order-term comes first.

        Coefficients c[-order-1:] are modified in-place.
    x : ndarray, shape (m+1,)
        Breakpoints of polynomials
    order : int
        Order up to which enforce piecewise differentiability.

    """

    cdef int ip, jp, kp, dx
    cdef int interval
    cdef double_or_complex z, res
    cdef double prefactor, xval

    # check derivative order
    if order < 0:
        raise ValueError("Order of derivative cannot be negative")

    # shape checks
    if c.shape[1] != x.shape[0] - 1:
        raise ValueError("x and c have incompatible shapes")
    if order >= c.shape[0] - 1:
        raise ValueError("order too large")
    if order < 0:
        raise ValueError("order negative")

    # evaluate
    for ip in range(1, len(x)-1):
        xval = x[ip]
        interval = ip - 1

        for jp in range(c.shape[2]):
            # ensure continuity for derivatives, starting at the
            # highest one (the lower derivatives depend on the higher
            # ones, but not vice versa)
            for dx in range(order, -1, -1):
                # evaluate dx-th derivative of the polynomial in previous interval
                res = evaluate_poly1(xval - x[interval], c, interval, jp, dx)

                # set dx-th coefficient of polynomial in current
                # interval so that the dx-th derivative is continuous
                for kp in range(dx):
                    res /= kp + 1

                c[c.shape[0] - dx - 1, ip, jp] = res


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def integrate(double_or_complex[:,:,::1] c,
              double[::1] x,
              double a,
              double b,
              double_or_complex[::1] out):
    """
    Compute integral over a piecewise polynomial.

    Parameters
    ----------
    c : ndarray, shape (k, m, n)
        Coefficients local polynomials of order `k-1` in `m` intervals.
    x : ndarray, shape (m+1,)
        Breakpoints of polynomials
    a : double
        Start point of integration.
    b : double
        End point of integration.
    out
        Integral of the piecewise polynomial, assuming the polynomial
        is zero outside the range (x[0], x[-1]).

    """

    cdef int ip, jp, kp, dx
    cdef int start_interval, end_interval, interval, sign
    cdef double_or_complex z, res, va, vb, vtot
    cdef double prefactor, xval

    # shape checks
    if c.shape[1] != x.shape[0] - 1:
        raise ValueError("x and c have incompatible shapes")
    if out.shape[0] != c.shape[2]:
        raise ValueError("x and c have incompatible shapes")

    # fix integration order
    if not (b >= a):
        raise ValueError("Integral bounds not in order")

    # find intervals
    start_interval = find_interval(x, a)
    if start_interval < 0:
        raise ValueError("a out of range")

    end_interval = find_interval(x, b)
    if end_interval < 0:
        raise ValueError("b out of range")

    # evaluate
    for jp in range(c.shape[2]):
        vtot = 0
        for interval in range(start_interval, end_interval+1):
            # local antiderivative, end point
            if interval == end_interval:
                vb = evaluate_poly1(b - x[interval], c, interval, jp, dx=-1)
            else:
                vb = evaluate_poly1(x[interval+1] - x[interval], c, interval, jp, dx=-1)

            # local antiderivative, start point
            if interval == start_interval:
                va = evaluate_poly1(a - x[interval], c, interval, jp, dx=-1)
            else:
                va = evaluate_poly1(0, c, interval, jp, dx=-1)

            # integral
            vtot = vtot + (vb - va)

        out[jp] = vtot


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef int find_interval(double[::1] x,
                       double xval,
                       int prev_interval=0) nogil:
    """
    Find an interval such that x[interval] <= xval < x[interval+1]

    Parameters
    ----------
    x : array of double, shape (m,)
        Piecewise polynomial breakpoints
    xval : double
        Point to find
    prev_interval : int, optional
        Interval where a previous point was found

    Returns
    -------
    interval : int
        Interval where previous point was found, or
        -1 below lower bound, -2 if above upper bound,
        or -3 if nan.

    """
    cdef int interval, high, low, mid
    cdef double a, b

    a = x[0]
    b = x[x.shape[0]-1]

    interval = prev_interval
    if interval < 0 or interval >= x.shape[0]:
        interval = 0

    if not (a <= xval <= b):
        # Out-of-bounds (or nan)
        if xval < a:
            interval = -1 # below
        elif xval > b:
            interval = -2 # above
        else:
            interval = -3 # nan or something
    elif xval == b:
        # Make the interval closed from the right
        interval = x.shape[0] - 2
    else:
        # Find the interval the coordinate is in
        # (binary search with locality)
        if xval >= x[interval]:
            low = interval
            high = x.shape[0]-2
        else:
            low = 0
            high = interval

        if xval < x[low+1]:
            high = low

        while low < high:
            mid = (high + low)//2
            if xval < x[mid]:
                # mid < high
                high = mid
            elif xval >= x[mid + 1]:
                low = mid + 1
            else:
                # x[mid] <= xval < x[mid+1]
                low = mid
                break

        interval = low

    return interval


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef double_or_complex evaluate_poly1(double s, double_or_complex[:,:,::1] c, int ci, int cj, int dx) nogil:
    """
    Evaluate polynomial, derivative, or antiderivative in a single interval.

    Antiderivatives are evaluated assuming zero integration constants.

    Parameters
    ----------
    s : double
        Polynomial x-value
    c : double[:,:,:]
        Polynomial coefficients. c[:,ci,cj] will be used
    ci, cj : int
        Which of the coefs to use
    dx : int
        Order of derivative (> 0) or antiderivative (< 0) to evaluate.

    """
    cdef int kp, k
    cdef double_or_complex res, z
    cdef double prefactor

    res = 0.0
    z = 1.0

    if dx < 0:
        for k in range(-dx):
            z *= s

    for kp in range(c.shape[0]):
        # prefactor of term after differentiation
        if dx == 0:
            prefactor = 1.0
        elif dx > 0:
            # derivative
            if kp < dx:
                continue
            else:
                prefactor = 1.0
                for k in range(kp, kp - dx, -1):
                    prefactor *= k
        else:
            # antiderivative
            prefactor = 1.0
            for k in range(kp, kp - dx):
                prefactor /= k + 1

        res = res + c[c.shape[0] - kp - 1, ci, cj] * z * prefactor

        # compute x**max(k-dx,0)
        if kp < c.shape[0] - 1 and kp >= dx:
            z *= s

    return res
