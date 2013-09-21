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
             int at_endpoint,
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
    at_endpoint : bool
        Evaluate the value of the polynomial at the endpoints of each
        break points, xp[j]=x[j+1]. If True, the parameter `xp` is ignored
        and `out` should be of length `x-1`.

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
    cdef double a, b, prefactor, xval

    # check derivative order
    if dx < 0:
        raise ValueError("Order of derivative cannot be negative")

    # shape checks
    if at_endpoint:
        if out.shape[0] != x.shape[0] - 1:
            raise ValueError("at_endpoint given but out has wrong shape")
    else:
        if out.shape[0] != xp.shape[0]:
            raise ValueError("out and xp have incompatible shapes")
    if out.shape[1] != c.shape[2]:
        raise ValueError("out and c have incompatible shapes")

    # evaluate
    a = x[0]
    b = x[x.shape[0]-1]

    interval = 0
    has_out_of_bounds = 0

    for ip in range(len(xp)):
        if at_endpoint:
            # Evaluate at interval endpoints
            xval = x[ip+1]
            interval = ip
        else:
            xval = xp[ip]

            # Find correct interval
            if not (a <= xval <= b):
                # Out-of-bounds (or nan)
                has_out_of_bounds = 1
                for jp in range(c.shape[2]):
                    out[ip, jp] = nan
                continue
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

                assert x[interval] <= xval < x[interval+1]

        # Evaluate the local polynomial(s)
        for jp in range(c.shape[2]):
            out[ip, jp] = 0

            res = 0
            z = 1.0

            for kp in range(c.shape[0]):
                # prefactor of term after differentiation
                if kp < dx:
                    continue
                else:
                    prefactor = 1.0
                    for i in range(kp, kp - dx, -1):
                        prefactor *= i

                # sum terms
                if prefactor == 1:
                    res = res + c[c.shape[0] - kp - 1, interval, jp] * z
                else:
                    res = res + c[c.shape[0] - kp - 1, interval, jp] * z * prefactor

                # compute x**max(k-dx,0)
                if kp < c.shape[0] - 1 and kp >= dx:
                    z *= xval - x[interval]

            out[ip, jp] = res

    return has_out_of_bounds
