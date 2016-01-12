"""
Routines for evaluating and manipulating piecewise polynomials in
local power basis.

"""

from .polyint import _Interpolator1D
import numpy as np

cimport cython

cimport libc.stdlib
cimport libc.math

ctypedef double complex double_complex

ctypedef fused double_or_complex:
    double
    double complex

cdef extern from "blas_defs.h":
    void c_dgeev(char *jobvl, char *jobvr, int *n, double *a,
                 int *lda, double *wr, double *wi, double *vl, int *ldvl,
                 double *vr, int *ldvr, double *work, int *lwork,
                 int *info)

cdef extern from "numpy/npy_math.h":
    double nan "NPY_NAN"

#------------------------------------------------------------------------------
# Piecewise power basis polynomials
#------------------------------------------------------------------------------

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def evaluate(double_or_complex[:,:,::1] c,
             double[::1] x,
             double[::1] xp,
             int dx,
             int extrapolate,
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
    extrapolate : int, optional
        Whether to extrapolate to out-of-bounds points based on first
        and last intervals, or to return NaNs.
    out : ndarray, shape (r, n)
        Value of each polynomial at each of the input points.
        This argument is modified in-place.

    """

    cdef int ip, jp
    cdef int interval
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

    for ip in range(len(xp)):
        xval = xp[ip]

        # Find correct interval
        i = find_interval(x, xval, interval, extrapolate)
        if i < 0:
            # xval was nan etc
            for jp in range(c.shape[2]):
                out[ip, jp] = nan
            continue
        else:
            interval = i

        # Evaluate the local polynomial(s)
        for jp in range(c.shape[2]):
            out[ip, jp] = evaluate_poly1(xval - x[interval], c, interval, jp, dx)


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
              int extrapolate,
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
    extrapolate : int, optional
        Whether to extrapolate to out-of-bounds points based on first
        and last intervals, or to return NaNs.
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
    start_interval = find_interval(x, a, 0, extrapolate)
    if start_interval < 0:
        out[:] = nan
        return

    end_interval = find_interval(x, b, 0, extrapolate)
    if end_interval < 0:
        out[:] = nan
        return

    # evaluate
    for jp in range(c.shape[2]):
        vtot = 0
        for interval in range(start_interval, end_interval+1):
            # local antiderivative, end point
            if interval == end_interval:
                vb = evaluate_poly1(b - x[interval], c, interval, jp, -1)
            else:
                vb = evaluate_poly1(x[interval+1] - x[interval], c, interval, jp, -1)

            # local antiderivative, start point
            if interval == start_interval:
                va = evaluate_poly1(a - x[interval], c, interval, jp, -1)
            else:
                va = evaluate_poly1(0, c, interval, jp, -1)

            # integral
            vtot = vtot + (vb - va)

        out[jp] = vtot


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def real_roots(double[:,:,::1] c, double[::1] x, int report_discont,
               int extrapolate):
    """
    Compute real roots of a real-valued piecewise polynomial function.

    If a section of the piecewise polynomial is identically zero, the
    values (x[begin], nan) are appended to the root list.

    If the piecewise polynomial is not continuous, and the sign
    changes across a breakpoint, the breakpoint is added to the root
    set if `report_discont` is True.

    Parameters
    ----------
    c, x
        Polynomial coefficients, as above
    report_discont : int, optional
        Whether to report discontinuities across zero at breakpoints
        as roots
    extrapolate : int, optional
        Whether to consider roots obtained by extrapolating based
        on first and last intervals.

    """
    cdef list roots
    cdef list cur_roots
    cdef int interval, jp, k, i, p

    cdef double *wr
    cdef double *wi
    cdef double last_root, va, vb
    cdef double f, df, dx
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

                    # Refine root by one Newton iteration
                    f = evaluate_poly1(wr[i], c, interval, jp, 0)
                    df = evaluate_poly1(wr[i], c, interval, jp, 1)
                    if df != 0:
                        dx = f/df
                        if abs(dx) < abs(wr[i]):
                            wr[i] = wr[i] - dx

                    # Check interval
                    wr[i] += x[interval]
                    if interval == 0 and extrapolate:
                        # Half-open to the left
                        if not wr[i] <= x[interval+1]:
                            continue
                    elif interval == c.shape[1] - 1 and extrapolate:
                        # Half-open to the right
                        if not wr[i] >= x[interval]:
                            continue
                    else:
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
                       int prev_interval=0,
                       int extrapolate=1) nogil:
    """
    Find an interval such that x[interval] <= xval < x[interval+1]
    or interval == 0 and xval < x[0]
    or interval == n-2 and xval > x[n-1]

    Parameters
    ----------
    x : array of double, shape (m,)
        Piecewise polynomial breakpoints
    xval : double
        Point to find
    prev_interval : int, optional
        Interval where a previous point was found
    extrapolate : int, optional
        Whether to return the last of the first interval if the
        point is out-of-bounds. 

    Returns
    -------
    interval : int
        Suitable interval or -1 if nan.

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
        if xval < a and extrapolate:
            # below
            interval = 0
        elif xval > b and extrapolate:
            # above
            interval = x.shape[0] - 2
        else:
            # nan or no extrapolation
            interval = -1
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
                      void **workspace):
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
         here after call. The roots are sorted in increasing order according
         to the real part.
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
    cdef double *a
    cdef double *work
    cdef double a0, a1, a2, d, br, bi
    cdef int lwork, n, i, j, order
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
    elif order == 1:
        # Low-order polynomial: a0*x + a1
        a0 = c[n-1-order,ci,cj]
        a1 = c[n-1-order+1,ci,cj]
        wr[0] = -a1 / a0
        wi[0] = 0
        return 1
    elif order == 2:
        # Low-order polynomial: a0*x**2 + a1*x + a2
        a0 = c[n-1-order,ci,cj]
        a1 = c[n-1-order+1,ci,cj]
        a2 = c[n-1-order+2,ci,cj]

        d = a1*a1 - 4*a0*a2
        if d < 0:
            # no real roots
            d = libc.math.sqrt(-d)
            wr[0] = -a1/(2*a0)
            wi[0] = -d/(2*a0)
            wr[1] = -a1/(2*a0)
            wi[1] = d/(2*a0)
            return 2

        d = libc.math.sqrt(d)

        # avoid cancellation in subtractions
        if d == 0:
            wr[0] = -a1/(2*a0)
            wi[0] = 0
            wr[1] = -a1/(2*a0)
            wi[1] = 0
        elif a1 < 0:
            wr[0] = (2*a2) / (-a1 + d) # == (-a1 - d)/(2*a0)
            wi[0] = 0
            wr[1] = (-a1 + d) / (2*a0)
            wi[1] = 0
        else:
            wr[0] = (-a1 - d)/(2*a0)
            wi[0] = 0
            wr[1] = (2*a2) / (-a1 - d) # == (-a1 + d)/(2*a0)
            wi[1] = 0

        return 2

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
        a[j + (order-1)*order] = -c[n-1-j,ci,cj]/c[n-1-order,ci,cj]
        if j + 1 < order:
            a[j+1 + order*j] = 1

    # Compute companion matrix eigenvalues
    info = 0
    c_dgeev("N", "N", &order, a, &order, <double*>wr, <double*>wi,
            NULL, &order, NULL, &order, work, &lwork, &info)
    if info != 0:
        # Failure
        return -2

    # Sort roots (insertion sort)
    for i in range(order):
        br = wr[i]
        bi = wi[i]
        for j in range(i - 1, -1, -1):
            if wr[j] > br:
                wr[j+1] = wr[j]
                wi[j+1] = wi[j]
            else:
                wr[j+1] = br
                wi[j+1] = bi
                break
        else:
            wr[0] = br
            wi[0] = bi

    # Return with roots
    return order


def _croots_poly1(double[:,:,::1] c, double_complex[:,:,::1] w):
    """
    Find roots of polynomials.

    This function is for testing croots_poly1

    Parameters
    ----------
    c : ndarray, (k, m, n)
        Coefficients of several order-k polynomials
    w : ndarray, (k, m, n)
        Output argument --- roots of the polynomials.

    """

    cdef double *wr
    cdef double *wi
    cdef void *workspace
    cdef int i, j, k, nroots

    if (c.shape[0] != w.shape[0] or c.shape[1] != w.shape[1]
            or c.shape[2] != w.shape[2]):
        raise ValueError("c and w have incompatible shapes")
    if c.shape[0] <= 0:
        return

    wr = <double*>libc.stdlib.malloc(c.shape[0] * sizeof(double))
    wi = <double*>libc.stdlib.malloc(c.shape[0] * sizeof(double))
    workspace = NULL

    try:
        for i in range(c.shape[1]):
            for j in range(c.shape[2]):
                for k in range(c.shape[0]):
                    w[k,i,j] = nan

                nroots = croots_poly1(c, i, j, wr, wi, &workspace)

                if nroots == -1:
                    continue
                elif nroots < -1 or nroots >= c.shape[0]:
                    raise RuntimeError("root-finding failed")

                for k in range(nroots):
                    w[k,i,j].real = wr[k]
                    w[k,i,j].imag = wi[k]
    finally:
        if workspace != NULL:
            libc.stdlib.free(workspace)
        libc.stdlib.free(wr)
        libc.stdlib.free(wi)


#------------------------------------------------------------------------------
# Piecewise Bernstein basis polynomials
#------------------------------------------------------------------------------

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef double_or_complex evaluate_bpoly1(double_or_complex s,
                                       double_or_complex[:,:,::1] c,
                                       int ci, int cj) nogil:
    """
    Evaluate polynomial in the Bernstein basis in a single interval.

    A Bernstein polynomial is defined as

        .. math:: b_{j, k} = comb(k, j) x^{j} (1-x)^{k-j}

    with ``0 <= x <= 1``.

    Parameters
    ----------
    s : double
        Polynomial x-value
    c : double[:,:,:]
        Polynomial coefficients. c[:,ci,cj] will be used
    ci, cj : int
        Which of the coefs to use

    """
    cdef int k, j
    cdef double_or_complex res, s1, comb

    k = c.shape[0] - 1  # polynomial order
    s1 = 1. - s

    # special-case lowest orders
    if k == 0:
        res = c[0, ci, cj]
    elif k == 1: 
        res = c[0, ci, cj] * s1 + c[1, ci, cj] * s
    elif k == 2:
        res = c[0, ci, cj] * s1*s1 + c[1, ci, cj] * 2.*s1*s + c[2, ci, cj] * s*s
    elif k == 3:
        res = (c[0, ci, cj] * s1*s1*s1 + c[1, ci, cj] * 3.*s1*s1*s +
               c[2, ci, cj] * 3.*s1*s*s + c[3, ci, cj] * s*s*s)
    else:
        # XX: replace with de Casteljau's algorithm if needs be
        res, comb = 0., 1.
        for j in range(k+1):
            res += comb * s**j * s1**(k-j) * c[j, ci, cj]
            comb *= 1. * (k-j) / (j+1.)

    return res


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef double_or_complex evaluate_bpoly1_deriv(double_or_complex s,
                                             double_or_complex[:,:,::1] c,
                                             int ci, int cj,
                                             int nu,
                                             double_or_complex[:,:,::1] wrk) nogil:
    """
    Evaluate the derivative of a polynomial in the Bernstein basis 
    in a single interval.

    A Bernstein polynomial is defined as

        .. math:: b_{j, k} = comb(k, j) x^{j} (1-x)^{k-j}

    with ``0 <= x <= 1``.

    The algorithm is detailed in BPoly._construct_from_derivatives.

    Parameters
    ----------
    s : double
        Polynomial x-value
    c : double[:,:,:]
        Polynomial coefficients. c[:,ci,cj] will be used
    ci, cj : int
        Which of the coefs to use
    nu : int
        Order of the derivative to evaluate. Assumed strictly positive
        (no checks are made).
    wrk : double[:,:,::1]
        A work array, shape (c.shape[0]-nu, 1, 1).

    """
    cdef int k, j, a
    cdef double_or_complex res, term
    cdef double comb, poch

    k = c.shape[0] - 1  # polynomial order

    if nu == 0:
        res = evaluate_bpoly1(s, c, ci, cj)
    else:
        poch = 1.
        for a in range(nu):
            poch *= k - a

        term = 0.
        for a in range(k - nu + 1):
            term, comb = 0., 1.
            for j in range(nu+1):
                term += c[j+a, ci, cj] * (-1)**(j+nu) * comb
                comb *= 1. * (nu-j) / (j+1)
            wrk[a, 0, 0] = term * poch
        res = evaluate_bpoly1(s, wrk, 0, 0)
    return res

#
# Evaluation; only differs from _ppoly by evaluate_poly1 -> evaluate_bpoly1
#
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def evaluate_bernstein(double_or_complex[:,:,::1] c,
             double[::1] x,
             double[::1] xp,
             int nu,
             int extrapolate,
             double_or_complex[:,::1] out):
    """
    Evaluate a piecewise polynomial in the Bernstein basis.

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
    nu : int
        Order of derivative to evaluate.  The derivative is evaluated
        piecewise and may have discontinuities.
    extrapolate : int, optional
        Whether to extrapolate to out-of-bounds points based on first
        and last intervals, or to return NaNs.
    out : ndarray, shape (r, n)
        Value of each polynomial at each of the input points.
        This argument is modified in-place.

    """

    cdef int ip, jp
    cdef int interval
    cdef double xval
    cdef double_or_complex s, ds, ds_nu
    cdef double_or_complex[:,:,::1] wrk

    # check derivative order
    if nu < 0:
        raise NotImplementedError("Cannot do antiderivatives in the B-basis yet.")

    # shape checks
    if out.shape[0] != xp.shape[0]:
        raise ValueError("out and xp have incompatible shapes")
    if out.shape[1] != c.shape[2]:
        raise ValueError("out and c have incompatible shapes")
    if c.shape[1] != x.shape[0] - 1:
        raise ValueError("x and c have incompatible shapes")

    if nu > 0:
        if double_or_complex is double_complex:
            wrk = np.empty((c.shape[0]-nu, 1, 1), dtype=np.complex_)
        else:
            wrk = np.empty((c.shape[0]-nu, 1, 1), dtype=np.float_)
        
    # evaluate
    interval = 0

    for ip in range(len(xp)):
        xval = xp[ip]

        # Find correct interval
        i = find_interval(x, xval, interval, extrapolate)
        if i < 0:
            # xval was nan etc
            for jp in range(c.shape[2]):
                out[ip, jp] = nan
            continue
        else:
            interval = i

        # Evaluate the local polynomial(s)
        ds = x[interval+1] - x[interval]
        ds_nu = ds**nu
        for jp in range(c.shape[2]):
            s = (xval - x[interval]) / ds
            if nu == 0:
                out[ip, jp] = evaluate_bpoly1(s, c, interval, jp)
            else:
                out[ip, jp] = evaluate_bpoly1_deriv(s, c, interval, jp,
                        nu, wrk) / ds_nu
