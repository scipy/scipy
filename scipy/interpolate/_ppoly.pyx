"""
Routines for evaluating and manipulating piecewise polynomials in
local power basis.

"""

from .polyint import _Interpolator1D
import numpy as np

cimport cython

cdef double nan = np.nan

cimport libc.stdlib

ctypedef double complex double_complex

ctypedef fused double_or_complex:
    double
    double complex

cdef extern:
    void dgeev_(char *jobvl, char *jobvr, int *n, double *a,
                int *lda, double *wr, double *wi, double *vl, int *ldvl,
                double *vr, int *ldvr, double *work, int *lwork,
                int *info) nogil


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
    out : ndarray, shape (r, n)
        Value of each polynomial at each of the input points.
        For points outside the span ``x[0] ... x[-1]``,
        ``nan`` is returned.
        This argument is modified in-place.

    """

    cdef int ip, jp
    cdef int interval
    cdef int has_out_of_bounds
    cdef double xval

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
    cdef double_or_complex res
    cdef double xval

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
    out : ndarray, shape (n,)
        Integral of the piecewise polynomial, assuming the polynomial
        is zero outside the range (x[0], x[-1]).
        This argument is modified in-place.

    """

    cdef int jp
    cdef int start_interval, end_interval, interval
    cdef double_or_complex va, vb, vtot

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
def real_roots(double[:,:,::1] c, double[::1] x, int report_discont):
    """
    Compute real roots of a real-valued piecewise polynomial function.

    If a section of the piecewise polynomial is identically zero, the
    values (x[begin], nan) are appended to the root list.

    If the piecewise polynomial is not continuous, and the sign
    changes across a breakpoint, the breakpoint is added to the root
    set if `report_discont` is True.

    """
    cdef list roots
    cdef list cur_roots
    cdef int interval, jp, k, i

    cdef double *wr, *wi, last_root, va, vb
    cdef void *workspace

    if c.shape[1] != x.shape[0] - 1:
        raise ValueError("x and c have incompatible shapes")

    if c.shape[0] == 0:
        return np.array([], dtype=float)

    wr = <double*>libc.stdlib.malloc(c.shape[0] * sizeof(double))
    wi = <double*>libc.stdlib.malloc(c.shape[0] * sizeof(double))
    workspace = NULL

    last_root = nan

    roots = []
    try:
        for jp in range(c.shape[2]):
            cur_roots = []
            for interval in range(c.shape[1]):
                # Check for sign change across intervals
                if interval > 0 and report_discont:
                    va = evaluate_poly1(x[interval] - x[interval-1], c, interval-1, jp, 0)
                    vb = evaluate_poly1(0, c, interval, jp, 0)
                    if (va < 0 and vb > 0) or (va > 0 and vb < 0):
                        # sign change between intervals
                        if x[interval] != last_root:
                            last_root = x[interval]
                            cur_roots.append(float(last_root))

                # Compute first the complex roots
                k = croots_poly1(c, interval, jp, wr, wi, &workspace)

                # Check for errors and identically zero values
                if k == -1:
                    # Zero everywhere
                    if x[interval] == x[interval+1]:
                        # Only a point
                        if x[interval] != last_root:
                            last_root = x[interval]
                            cur_roots.append(x[interval])
                    else:
                        # A real interval
                        cur_roots.append(x[interval])
                        cur_roots.append(np.nan)
                        last_root = nan
                    continue
                elif k < -1:
                    # An error occurred
                    raise RuntimeError("Internal error in root finding; "
                                       "please report this bug")
                elif k == 0:
                    # No roots
                    continue

                # Filter real roots
                for i in range(k):
                    # Check real root
                    #
                    # The reality of a root is a decision that can be left to LAPACK,
                    # which has to determine this in any case.
                    if wi[i] != 0:
                        continue

                    # Check interval
                    wr[i] += x[interval]
                    if not (x[interval] <= wr[i] <= x[interval+1]):
                        continue

                    # Add to list
                    if wr[i] != last_root:
                        last_root = wr[i]
                        cur_roots.append(float(last_root))

            # Construct roots
            roots.append(np.array(cur_roots, dtype=float))
    finally:
        if workspace != NULL:
            libc.stdlib.free(workspace)
        libc.stdlib.free(wr)
        libc.stdlib.free(wi)

    return roots


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


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef int croots_poly1(double[:,:,::1] c, int ci, int cj, double* wr, double* wi,
                      void **workspace) nogil:
    """
    Find all complex roots of a local polynomial.

    Parameters
    ----------
    c : ndarray, shape (k, m, n)
         Coefficients of polynomials of order k
    ci, cj : int
         Index of the local polynomial whose coefficients c[:,ci,cj] to use
    wr, wi : double*
         Allocated double arrays of size `k`. The complex roots are stored
         here after call.
    workspace : double**
         Work space pointer. workspace[0] should be NULL on initial
         call.  Multiple subsequent calls with same `k` can share the
         same `workspace`.  If workspace[0] is non-NULL after the
         calls, it must be freed with libc.stdlib.free.

    Returns
    -------
    nroots : int
        How many roots found for the polynomial.
        If `-1`, the polynomial is identically zero.
        If `< -1`, an error occurred.

    Notes
    -----
    Uses LAPACK + the companion matrix method.

    """
    cdef double *a, *work
    cdef int lwork, n, j, order
    cdef int nworkspace, info

    n = c.shape[0]

    # Check actual polynomial order
    for j in range(n):
        if c[j,ci,cj] != 0:
            order = n - 1 - j
            break
    else:
        order = -1

    if order < 0:
        # Zero everywhere
        return -1
    elif order == 0:
        # Nonzero constant polynomial: no roots
        return 0

    # Compute required workspace and allocate it
    lwork = 1 + 8*n

    if workspace[0] == NULL:
        nworkspace = n*n + lwork
        workspace[0] = libc.stdlib.malloc(nworkspace * sizeof(double))

    a = <double*>workspace[0]
    work = a + n*n

    # Initialize the companion matrix, Fortran order
    for j in range(order*order):
        a[j] = 0
    for j in range(order):
        a[j + (order-1)*order] = -c[c.shape[0]-1-j,ci,cj]/c[c.shape[0]-1-order,ci,cj]
        if j + 1 < order:
            a[j+1 + order*j] = 1

    # Compute companion matrix eigenvalues
    info = 0
    dgeev_("N", "N", &order, a, &order, <double*>wr, <double*>wi,
           NULL, &order, NULL, &order, work, &lwork, &info)
    if info != 0:
        # Failure
        return -2

    # Return with roots
    return order
