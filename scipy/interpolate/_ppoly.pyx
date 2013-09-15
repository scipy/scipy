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

    Returns
    -------
    out : ndarray, shape (r, n)
        Value of each polynomial at each of the input points.
        For points outside the span ``x[0] ... x[-1]``,
        ``nan`` is returned.

    """

    cdef int ip, jp, kp
    cdef int interval, high, low, mid
    cdef int has_out_of_bounds
    cdef double_or_complex z, res
    cdef double a, b

    a = x[0]
    b = x[x.shape[0]-1]

    interval = 0
    has_out_of_bounds = 0

    for ip in range(len(xp)):
        # Find correct interval
        if not (a <= xp[ip] <= b):
            # Out-of-bounds (or nan)
            has_out_of_bounds = 1
            for jp in range(c.shape[2]):
                out[ip, jp] = nan
            continue
        elif xp[ip] == b:
            # Make the interval closed from the right
            interval = x.shape[0] - 2
        else:
            # Find the interval the coordinate is in
            # (binary search with locality)
            if xp[ip] >= x[interval]:
                low = interval
                high = x.shape[0]-2
            else:
                low = 0
                high = interval

            if xp[ip] < x[low+1]:
                high = low

            while low < high:
                mid = (high + low)//2
                if xp[ip] < x[mid]:
                    # mid < high
                    high = mid
                elif xp[ip] >= x[mid + 1]:
                    low = mid + 1
                else:
                    # x[mid] <= xp[ip] < x[mid+1]
                    low = mid
                    break

            interval = low

            assert x[interval] <= xp[ip] < x[interval+1]

        # Evaluate the local polynomial(s)
        for jp in range(c.shape[2]):
            out[ip, jp] = 0

        z = 1.0
        for kp in range(c.shape[0]):
            for jp in range(c.shape[2]):
                out[ip, jp] = out[ip, jp] + c[c.shape[0] - kp - 1, interval, jp] * z
            if kp < c.shape[0] - 1:
                z *= xp[ip] - x[interval]

    return has_out_of_bounds
