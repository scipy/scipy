"""
Routines for evaluating and manipulating B-splines.

"""

import numpy as np
cimport numpy as cnp

from numpy cimport npy_intp, npy_int64, npy_int32

cimport cython
from libc.math cimport NAN

cnp.import_array()

cdef extern from "src/__fitpack.h":
    void _deBoor_D(const double *t, double x, int k, int ell, int m, double *result) nogil


ctypedef fused int32_or_int64:
    cnp.npy_int32
    cnp.npy_int64

#------------------------------------------------------------------------------
# B-splines
#------------------------------------------------------------------------------

@cython.wraparound(False)
@cython.boundscheck(False)
cdef inline int find_interval(const double[::1] t,
                       int k,
                       double xval,
                       int prev_l,
                       bint extrapolate) noexcept nogil:
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
        value to find the interval for
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

    if ((xval < tb) or (xval > te)) and not extrapolate:
        return -1

    l = prev_l if k < prev_l < n else k

    # xval is in support, search for interval s.t. t[interval] <= xval < t[l+1]
    while(xval < t[l] and l != k):
        l -= 1

    l += 1
    while(xval >= t[l] and l != n):
        l += 1

    return l-1


# NB: a python wrapper for find_interval. The leading underscore signals
# it's not meant to be user-visible outside of _bsplines.py
@cython.wraparound(False)
@cython.boundscheck(False)
def _find_interval(const double[::1] t,
                   int k,
                   double xval,
                   int prev_l,
                   bint extrapolate):
    return find_interval(t, k, xval, prev_l, extrapolate)


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
def evaluate_spline(const double[::1] t,
             const double[:, ::1] c,
             int k,
             const double[::1] xp,
             int nu,
             bint extrapolate,
             double[:, ::1] out):
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
    nu : int
        Order of derivative to evaluate.
    extrapolate : int, optional
        Whether to extrapolate to ouf-of-bounds points, or to return NaNs.
    out : ndarray, shape (s, m)
        Computed values of the spline at each of the input points.
        This argument is modified in-place.

    """

    cdef int ip, jp, a
    cdef int interval
    cdef double xval

    # shape checks
    if out.shape[0] != xp.shape[0]:
        raise ValueError("out and xp have incompatible shapes")
    if out.shape[1] != c.shape[1]:
        raise ValueError("out and c have incompatible shapes")

    # check derivative order
    if nu < 0:
        raise NotImplementedError("Cannot do derivative order %s." % nu)

    cdef double[::1] work = np.empty(2*k+2, dtype=np.float64)

    # evaluate
    with nogil:
        interval = k
        for ip in range(xp.shape[0]):
            xval = xp[ip]

            # Find correct interval
            interval = find_interval(t, k, xval, interval, extrapolate)

            if interval < 0:
                # xval was nan etc
                for jp in range(c.shape[1]):
                    out[ip, jp] = NAN
                continue

            # Evaluate (k+1) b-splines which are non-zero on the interval.
            # on return, first k+1 elements of work are B_{m-k},..., B_{m}
            _deBoor_D(&t[0], xval, k, interval, nu, &work[0])

            # Form linear combinations
            for jp in range(c.shape[1]):
                out[ip, jp] = 0.
                for a in range(k+1):
                    out[ip, jp] = out[ip, jp] + c[interval + a - k, jp] * work[a]


def evaluate_all_bspl(const double[::1] t, int k, double xval, int m, int nu=0):
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
    >>> t = [0., 1., 2., 3., 4.]   # internal knots
    >>> a, b = t[0], t[-1]    # base interval is [a, b)
    >>> t = np.array([a]*k + t + [b]*k)  # add boundary knots

    >>> import matplotlib.pyplot as plt
    >>> xx = np.linspace(a, b, 100)
    >>> plt.plot(xx, BSpline.basis_element(t[k:-k])(xx),
    ...          lw=3, alpha=0.5, label='basis_element')

    Now we use slide an interval ``t[m]..t[m+1]`` along the base interval
    ``a..b`` and use `evaluate_all_bspl` to compute the restriction of
    the B-spline of interest to this interval:

    >>> for i in range(k+1):
    ...    x1, x2 = t[2*k - i], t[2*k - i + 1]
    ...    xx = np.linspace(x1 - 0.5, x2 + 0.5)
    ...    yy = [evaluate_all_bspl(t, k, x, 2*k - i)[i] for x in xx]
    ...    plt.plot(xx, yy, '--', label=str(i))
    ...
    >>> plt.grid(True)
    >>> plt.legend()
    >>> plt.show()

    """
    bbb = np.empty(2*k+2, dtype=np.float64)
    cdef double[::1] work = bbb
    _deBoor_D(&t[0], xval, k, m, nu, &work[0])
    return bbb[:k+1]


@cython.wraparound(False)
@cython.boundscheck(False)
def _colloc(const double[::1] x, const double[::1] t, int k, double[::1, :] ab,
            int offset=0):
    """Build the B-spline collocation matrix.

    The collocation matrix is defined as :math:`B_{j,l} = B_l(x_j)`,
    so that row ``j`` contains all the B-splines which are non-zero
    at ``x_j``.

    The matrix is constructed in the LAPACK banded storage.
    Basically, for an N-by-N matrix A with ku upper diagonals and
    kl lower diagonals, the shape of the array Ab is (2*kl + ku +1, N),
    where the last kl+ku+1 rows of Ab contain the diagonals of A, and
    the first kl rows of Ab are not referenced.
    For more info see, e.g. the docs for the ``*gbsv`` routine.

    This routine is not supposed to be called directly, and
    does no error checking.

    Parameters
    ----------
    x : ndarray, shape (n,)
        sorted 1D array of x values
    t : ndarray, shape (nt + k + 1,)
        sorted 1D array of knots
    k : int
        spline order
    ab : ndarray, shape (2*kl + ku + 1, nt), F-order
        This parameter is modified in-place.
        On exit: zeroed out.
        On exit: B-spline collocation matrix in the band storage with
        ``ku`` upper diagonals and ``kl`` lower diagonals.
        Here ``kl = ku = k``.
    offset : int, optional
        skip this many rows

    """
    cdef int left, j, a, kl, ku, clmn
    cdef double xval

    kl = ku = k
    cdef double[::1] wrk = np.empty(2*k + 2, dtype=np.float64)

    # collocation matrix
    with nogil:
        left = k
        for j in range(x.shape[0]):
            xval = x[j]
            # find interval
            left = find_interval(t, k, xval, left, extrapolate=False)

            # fill a row
            _deBoor_D(&t[0], xval, k, left, 0, &wrk[0])
            # for a full matrix it would be ``A[j + offset, left-k:left+1] = bb``
            # in the banded storage, need to spread the row over
            for a in range(k+1):
                clmn = left - k + a
                ab[kl + ku + j + offset - clmn, clmn] = wrk[a]


@cython.wraparound(False)
@cython.boundscheck(False)
def _handle_lhs_derivatives(const double[::1]t, int k, double xval,
                            double[::1, :] ab,
                            int kl, int ku,
                            const cnp.npy_long[::1] deriv_ords,
                            int offset=0):
    """ Fill in the entries of the collocation matrix corresponding to known
    derivatives at xval.

    The collocation matrix is in the banded storage, as prepared by _colloc.
    No error checking.

    Parameters
    ----------
    t : ndarray, shape (nt + k + 1,)
        knots
    k : integer
        B-spline order
    xval : float
        The value at which to evaluate the derivatives at.
    ab : ndarray, shape(2*kl + ku + 1, nt), Fortran order
        B-spline collocation matrix.
        This argument is modified *in-place*.
    kl : integer
        Number of lower diagonals of ab.
    ku : integer
        Number of upper diagonals of ab.
    deriv_ords : 1D ndarray
        Orders of derivatives known at xval
    offset : integer, optional
        Skip this many rows of the matrix ab.

    """
    cdef:
        int left, nu, a, clmn, row
        double[::1] wrk = np.empty(2*k+2, dtype=np.float64)

    # derivatives @ xval
    with nogil:
        left = find_interval(t, k, xval, k, extrapolate=False)
        for row in range(deriv_ords.shape[0]):
            nu = deriv_ords[row]
            _deBoor_D(&t[0], xval, k, left, nu, &wrk[0])
            # if A were a full matrix, it would be just
            # ``A[row + offset, left-k:left+1] = bb``.
            for a in range(k+1):
                clmn = left - k + a
                ab[kl + ku + offset + row - clmn, clmn] = wrk[a]


@cython.wraparound(False)
@cython.boundscheck(False)
def _norm_eq_lsq(const double[::1] x,
                 const double[::1] t,
                 int k,
                 const double[:, ::1] y,
                 const double[::1] w,
                 double[::1, :] ab,
                 double[:, ::1] rhs):
    """Construct the normal equations for the B-spline LSQ problem.

    The observation equations are ``A @ c = y``, and the normal equations are
    ``A.T @ A @ c = A.T @ y``. This routine fills in the rhs and lhs for the
    latter.

    The B-spline collocation matrix is defined as :math:`A_{j,l} = B_l(x_j)`,
    so that row ``j`` contains all the B-splines which are non-zero
    at ``x_j``.

    The normal eq matrix has at most `2k+1` bands and is constructed in the
    LAPACK symmetrix banded storage: ``A[i, j] == ab[i-j, j]`` with `i >= j`.
    See the doctsring for `scipy.linalg.cholesky_banded` for more info.

    This routine is not supposed to be called directly, and
    does no error checking.

    Parameters
    ----------
    x : ndarray, shape (n,)
        sorted 1D array of x values
    t : ndarray, shape (nt + k + 1,)
        sorted 1D array of knots
    k : int
        spline order
    y : ndarray, shape (n, s)
        a 2D array of y values. The second dimension contains all trailing
        dimensions of the original array of ordinates.
    w : ndarray, shape(n,)
        Weights.
    ab : ndarray, shape (k+1, n), in Fortran order.
        This parameter is modified in-place.
        On entry: should be zeroed out.
        On exit: LHS of the normal equations.
    rhs : ndarray, shape (n, s), in C order.
        This parameter is modified in-place.
        On entry: should be zeroed out.
        On exit: RHS of the normal equations.

    """
    cdef:
        int j, r, s, row, clmn, left, ci
        double xval, wval
        double[::1] wrk = np.empty(2*k + 2, dtype=np.float64)

    with nogil:
        left = k
        for j in range(x.shape[0]):
            xval = x[j]
            wval = w[j] * w[j]
            # find interval
            left = find_interval(t, k, xval, left, extrapolate=False)

            # non-zero B-splines at xval
            _deBoor_D(&t[0], xval, k, left, 0, &wrk[0])

            # non-zero values of A.T @ A: banded storage w/ lower=True
            # The colloq matrix in full storage would be
            #   A[j, left-k:left+1] = wrk,
            # Here we work out A.T @ A *in the banded storage* w/lower=True
            # see the docstring of `scipy.linalg.cholesky_banded`.
            for r in range(k+1):
                row = left - k + r
                for s in range(r+1):
                    clmn = left - k + s
                    ab[r-s, clmn] += wrk[r] * wrk[s] * wval

                # ... and A.T @ y
                for ci in range(rhs.shape[1]):
                    rhs[row, ci] = rhs[row, ci] + wrk[r] * y[j, ci] * wval

@cython.wraparound(False)
@cython.boundscheck(False)
def _make_design_matrix(const double[::1] x,
                        const double[::1] t,
                        int k,
                        bint extrapolate,
                        int32_or_int64[::1] indices):
    """
    Returns a design matrix in CSR format.

    Note that only indices is passed, but not indptr because indptr is already
    precomputed in the calling Python function design_matrix.

    Parameters
    ----------
    x : array_like, shape (n,)
        Points to evaluate the spline at.
    t : array_like, shape (nt,)
        Sorted 1D array of knots.
    k : int
        B-spline degree.
    extrapolate : bool, optional
        Whether to extrapolate to ouf-of-bounds points.
    indices : ndarray, shape (n * (k + 1),)
        Preallocated indices of the final CSR array.

    Returns
    -------
    data
        The data array of a CSR array of the b-spline design matrix.
        In each row all the basis elements are evaluated at the certain point
        (first row - x[0], ..., last row - x[-1]).

    indices
        The indices array of a CSR array of the b-spline design matrix.
    """
    cdef:
        cnp.npy_intp i, j, m, ind
        cnp.npy_intp n = x.shape[0]
        double[::1] work = np.empty(2*k+2, dtype=float)
        double[::1] data = np.zeros(n * (k + 1), dtype=float)
        double xval
    ind = k
    for i in range(n):
        xval = x[i]

        # Find correct interval. Note that interval >= 0 always as
        # extrapolate=False and out of bound values are already dealt with in
        # design_matrix
        ind = find_interval(t, k, xval, ind, extrapolate)
        _deBoor_D(&t[0], xval, k, ind, 0, &work[0])

        # data[(k + 1) * i : (k + 1) * (i + 1)] = work[:k + 1]
        # indices[(k + 1) * i : (k + 1) * (i + 1)] = np.arange(ind - k, ind + 1)
        for j in range(k + 1):
            m = (k + 1) * i + j
            data[m] = work[j]
            indices[m] = ind - k + j

    return np.asarray(data), np.asarray(indices)



@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
def evaluate_ndbspline(const double[:, ::1] xi,
                       const double[:, ::1] t,
                       const long[::1] len_t,
                       long[::1] k,
                       int[::1] nu,
                       bint extrapolate,
                       const double[::1] c1r,
                       npy_intp num_c_tr,
                       const npy_intp[::1] strides_c1,
                       const npy_intp[:, ::] indices_k1d,
                       double[:, ::1] out,
                      ):
        """Evaluate an N-dim tensor product spline or its derivative.

        Parameters
        ----------
        xi : ndarray, shape(npoints, ndim)
            ``npoints`` values to evaluate the spline at, each value is
            a point in an ``ndim``-dimensional space.
        t : ndarray, shape(ndim, max_len_t)
            Array of knots for each dimension.
            This array packs the tuple of knot arrays per dimension into a single
            2D array. The array is ragged (knot lengths may differ), hence
            the real knots in dimension ``d`` are ``t[d, :len_t[d]]``.
        len_t : ndarray, 1D, shape (ndim,)
            Lengths of the knot arrays, per dimension.
        k : tuple of ints, len(ndim)
            Spline degrees in each dimension.
        nu : ndarray of ints, shape(ndim,)
            Orders of derivatives to compute, per dimension.
        extrapolate : int
            Whether to extrapolate out of bounds or return nans.
        c1r: ndarray, one-dimensional
            Flattened array of coefficients.
            The original N-dimensional coefficient array ``c`` has shape
            ``(n1, ..., nd, ...)`` where each ``ni == len(t[d]) - k[d] - 1``,
            and the second "..." represents trailing dimensions of ``c``.
            In code, given the C-ordered array ``c``, ``c1r`` is
            ``c1 = c.reshape(c.shape[:ndim] + (-1,)); c1r = c1.ravel()``
        num_c_tr : int
            The number of elements of ``c1r``, which correspond to the trailing
            dimensions of ``c``. In code, this is
            ``c1 = c.reshape(c.shape[:ndim] + (-1,)); num_c_tr = c1.shape[-1]``.
        strides_c1 : ndarray, one-dimensional
            Pre-computed strides of the ``c1`` array.
            Note: These are *data* strides, not numpy-style byte strides.
            This array is equivalent to
            ``[stride // s1.dtype.itemsize for stride in s1.strides]``.
        indices_k1d : ndarray, shape((k+1)**ndim, ndim)
            Pre-computed mapping between indices for iterating over a flattened
            array of shape ``[k[d] + 1) for d in range(ndim)`` and
            ndim-dimensional indices of the ``(k+1,)*ndim`` dimensional array.
            This is essentially a transposed version of
            ``np.unravel_index(np.arange((k+1)**ndim), (k+1,)*ndim)``.
        out : ndarray, shape (npoints, num_c_tr)
            Output values of the b-spline at given ``xi`` points.

        Notes
        -----

        This function is essentially equivalent to the following: given an
        N-dimensional vector ``x = (x1, x2, ..., xN)``, iterate over the
        dimensions, form linear combinations of products,
        B(x1) * B(x2) * ... B(xN) of (k+1)**N b-splines which are non-zero
        at ``x``.

        Since b-splines are localized, the sum has (k+1)**N non-zero elements.

        If ``i = (i1, i2, ..., iN)`` is a vector if intervals of the knot
        vectors, ``t[d, id] <= xd < t[d, id+1]``, for ``d=1, 2, ..., N``, then
        the core loop of this function is nothing but

        ```
        result = 0
        iters = [range(i[d] - self.k[d], i[d] + 1) for d in range(ndim)]
        for idx in itertools.product(*iters):
            term = self.c[idx] * np.prod([B(x[d], self.k[d], idx[d], self.t[d])
                                          for d in range(ndim)])
            result += term
        ```

        For efficiency reasons, we iterate over the flattened versions of the
        arrays.

        """
        cdef:
            npy_intp ndim = len(t)

            # 'intervals': indices for a point in xi into the knot arrays t
            npy_intp[::1] i = np.empty(ndim, dtype=np.intp)

            # container for non-zero b-splines at each point in xi
            double[:, ::1] b = np.empty((ndim, max(k) + 1), dtype=float)

            const double[::1] xv     # an ndim-dimensional input point
            double xd               # d-th component of x

            const double[::1] td    # knots in dimension d

            npy_intp kd             # d-th component of k

            npy_intp i_c      # index to loop over range(num_c_tr)
            npy_intp iflat    # index to loop over (k+1)**ndim non-zero terms
            npy_intp volume   # the number of non-zero terms
            const npy_intp[:] idx_b   # ndim-dimensional index corresponding to iflat

            int out_of_bounds
            npy_intp idx_cflat_base, idx
            double factor
            double[::1] wrk = np.empty(2*max(k) + 2, dtype=float)

        if xi.shape[1] != ndim:
            raise ValueError(f"Expacted data points in {ndim}-D space, got"
                             f" {xi.shape[1]}-D points.")

        if out.shape[0] != xi.shape[0]:
            raise ValueError(f"out and xi are inconsistent: expected"
                             f" {xi.shape[0]} output values, got"
                             f" {out.shape[0]}.")
        if out.shape[1] != num_c_tr:
            raise ValueError(f"out and c are inconsistent: num_c={num_c_tr} "
                             f" and out.shape[1] = {out.shape[1]}.")


        with nogil:
            # the number of non-zero terms for each point in ``xi``.
            volume = 1
            for d in range(ndim):
                volume *= k[d] + 1

            ### Iterate over the data points
            for j in range(xi.shape[0]):
                xv = xi[j, :]

                # For each point, iterate over the dimensions
                out_of_bounds = 0
                for d in range(ndim):
                    td = t[d, :len_t[d]]
                    xd = xv[d]
                    kd = k[d]

                    # get the location of x[d] in t[d]
                    i[d] = find_interval(td, kd, xd, kd, extrapolate)

                    if i[d] < 0:
                        out_of_bounds = 1
                        break

                    # compute non-zero b-splines at this value of xd in dimension d
                    _deBoor_D(&td[0], xd, kd, i[d], nu[d], &wrk[0])
                    b[d, :kd+1] = wrk[:kd+1]

                if out_of_bounds:
                    # xd was nan or extrapolate=False: Fill the output array
                    # *for this xv value*, and continue to the next xv in xi.
                    for i_c in range(num_c_tr):
                        out[j, i_c] = NAN
                    continue

                for i_c in range(num_c_tr):
                    out[j, i_c] = 0.0

                # iterate over the direct products of non-zero b-splines
                for iflat in range(volume):
                    idx_b = indices_k1d[iflat, :]
                    # The line above is equivalent to
                    # idx_b = np.unravel_index(iflat, (k+1,)*ndim)

                    # From the indices in ``idx_b``, we prepare to index into
                    # c1.ravel() : for each dimension d, need to shift the index
                    # by ``i[d] - k[d]`` (see the docstring above).
                    #
                    # Since the strides of `c1` are pre-computed, and the array
                    # is already raveled and is guaranteed to be C-ordered, we only
                    # need to compute the base index for iterating over ``num_c_tr``
                    # elements which represent the trailing dimensions of ``c``.
                    #
                    # This all is essentially equivalent to iterating over
                    # idx_cflat = np.ravel_multi_index(tuple(idx_c) + (i_c,),
                    #                                  c1.shape)
                    idx_cflat_base = 0
                    factor = 1.0
                    for d in range(ndim):
                        factor *= b[d, idx_b[d]]
                        idx = idx_b[d] + i[d] - k[d]
                        idx_cflat_base += idx * strides_c1[d]

                    ### collect linear combinations of coef * factor
                    for i_c in range(num_c_tr):
                        out[j, i_c] = out[j, i_c] + c1r[idx_cflat_base + i_c] * factor


@cython.wraparound(False)
@cython.nonecheck(False)
@cython.boundscheck(False)
def _colloc_nd(const double[:, ::1] xvals, tuple t not None, const npy_int32[::1] k):
    """Construct the N-D tensor product collocation matrix as a CSR array.

    In the dense representation, each row of the collocation matrix corresponds
    to a data point and contains non-zero b-spline basis functions which are
    non-zero at this data point.

    Parameters
    ----------
    xvals : ndarray, shape(size, ndim)
        Data points. ``xvals[j, :]`` gives the ``j``-th data point as an
        ``ndim``-dimensional array.
    t : tuple of 1D arrays, length-ndim
        Tuple of knot vectors
    k : ndarray, shape (ndim,)
        Spline degrees

    Returns
    -------
    csr_data, csr_indices, csr_indptr
        The collocation matrix in the CSR array format.

    Notes
    -----
    Algorithm: given `xvals` and the tuple of knots `t`, we construct a tensor
    product spline, i.e. a linear combination of

       B(x1; i1, t1) * B(x2; i2, t2) * ... * B(xN; iN, tN)


    Here ``B(x; i, t)`` is the ``i``-th b-spline defined by the knot vector
    ``t`` evaluated at ``x``.

    Since ``B`` functions are localized, for each point `(x1, ..., xN)` we
    loop over the dimensions, and
    - find the location in the knot array, `t[i] <= x < t[i+1]`,
    - compute all non-zero `B` values
    - place these values into the relevant row

    In the dense representation, the collocation matrix would have had a row per
    data point, and each row has the values of the basis elements (i.e., tensor
    products of B-splines) evaluated at this data point. Since the matrix is very
    sparse (has size = len(x)**ndim, with only (k+1)**ndim non-zero elements per
    row), we construct it in the CSR format.
    """
    cdef:
        npy_intp size = xvals.shape[0]
        npy_intp ndim = xvals.shape[1]

        # 'intervals': indices for a point in xi into the knot arrays t
        npy_intp[::1] i = np.empty(ndim, dtype=np.intp)

        # container for non-zero b-splines at each point in xi
        double[:, ::1] b = np.empty((ndim, max(k) + 1), dtype=float)

        double xd               # d-th component of x
        const double[::1] td    # knots in the dimension d
        npy_intp kd             # d-th component of k

        npy_intp iflat    # index to loop over (k+1)**ndim non-zero terms
        npy_intp volume   # the number of non-zero terms
        npy_intp[:, ::1] _indices_k1d    # tabulated np.unravel_index

        # shifted indices into the data array
        npy_intp[::1] idx_c = np.ones(ndim, dtype=np.intp) * (-101)  # any sentinel would do, really
        npy_intp[::1] cstrides
        npy_intp idx_cflat

        npy_intp[::1] nu = np.zeros(ndim, dtype=np.intp)

        int out_of_bounds
        double factor
        double[::1] wrk = np.empty(2*max(k) + 2, dtype=float)

        # output
        double[::1] csr_data
        npy_int64[::1] csr_indices

        int j, d

    # the number of non-zero b-splines for each data point.
    k1_shape = tuple(kd + 1 for kd in k)
    volume = 1
    for d in range(ndim):
        volume *= k[d] + 1

    # Precompute the shape and strides of the coefficients array.
    # This would have been the NdBSpline coefficients; in the present context
    # this is a helper to compute the indices into the collocation matrix.
    c_shape = tuple(len(t[d]) - k1_shape[d] for d in range(ndim))

    # The computation is equivalent to
    # >>> x = np.empty(c_shape)
    # >>> cstrides = [s // 8 for s in x.strides]
    cs = c_shape[1:] + (1,)
    cstrides = np.cumprod(cs[::-1], dtype=np.intp)[::-1].copy()

    # tabulate flat indices for iterating over the (k+1)**ndim subarray of
    # non-zero b-spline elements
    indices = np.unravel_index(np.arange(volume), k1_shape)
    _indices_k1d = np.asarray(indices, dtype=np.intp).T.copy()

    # Allocate the collocation matrix in the CSR format.
    # If dense, this would have been
    # >>> matr = np.zeros((size, max_row_index), dtype=float)
    csr_indices = np.empty(shape=(size*volume,), dtype=np.int64)
    csr_data = np.empty(shape=(size*volume,), dtype=float)
    csr_indptr = np.arange(0, volume*size + 1, volume, dtype=np.int64)

    # ### Iterate over the data points ###
    for j in range(size):
        xv = xvals[j, :]

        # For each point, iterate over the dimensions
        out_of_bounds = 0
        for d in range(ndim):
            td = t[d]
            xd = xv[d]
            kd = k[d]

            # get the location of x[d] in t[d]
            i[d] = find_interval(td, kd, xd, kd, extrapolate=True)

            if i[d] < 0:
                out_of_bounds = 1
                break

            # compute non-zero b-splines at this value of xd in dimension d
            _deBoor_D(&td[0], xd, kd, i[d], nu[d], &wrk[0])
            b[d, :kd+1] = wrk[:kd+1]

        if out_of_bounds:
            raise ValueError(f"Out of bounds in {d = }, with {xv = }")

        # Iterate over the products of non-zero b-splines and place them
        # into the current row of the design matrix
        for iflat in range(volume):
            # the line below is an unrolled version of
            # idx_b = np.unravel_index(iflat,  tuple(kd+1 for kd in k))
            idx_b = _indices_k1d[iflat, :]

            factor = 1.0
            idx_cflat = 0
            for d in range(ndim):
                factor *= b[d, idx_b[d]]
                idx_c[d] = idx_b[d] + i[d] - k[d]
                idx_cflat += idx_c[d] * cstrides[d]

            # The `idx_cflat` computation above is an unrolled version of
            # idx_cflat = np.ravel_multi_index(tuple(idx_c), c_shape)

            # Fill the row of the collocation matrix in the CSR format.
            # If it were dense, it would have been just
            # >>> matr[j, idx_cflat] = factor

            # Each row of the full matrix has `volume` non-zero elements.
            # Thus the CSR format `indptr` increases in steps of `volume`
            csr_indices[j*volume + iflat] = idx_cflat
            csr_data[j*volume + iflat] = factor

    return np.asarray(csr_data), np.asarray(csr_indices), csr_indptr
