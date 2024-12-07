"""
Routines for evaluating and manipulating B-splines.

"""

import numpy as np
cimport numpy as cnp

from numpy cimport npy_intp, npy_int64, npy_int32

cimport cython
from libc.math cimport NAN

cnp.import_array()

cdef extern from "src/__fitpack.h" namespace "fitpack":
    void _deBoor_D(const double *t, double x, int k, int ell, int m, double *result
    ) noexcept nogil
    npy_int64 _find_interval(const double* tptr, npy_int64 len_t,
                           int k,
                           double xval,
                           npy_int64 prev_l,
                           int extrapolate
    ) noexcept nogil


#------------------------------------------------------------------------------
# B-splines
#------------------------------------------------------------------------------

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.nonecheck(False)
def evaluate_ndbspline(const double[:, ::1] xi,
                       const double[:, ::1] t,
                       const npy_int32[::1] len_t,
                       const npy_int32[::1] k,
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
                    i[d] = _find_interval(&td[0], td.shape[0], kd, xd, kd, extrapolate)

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
def _colloc_nd(const double[:, ::1] xvals,
               const double[:, ::1] _t,
               const npy_int32[::1] len_t,
               const npy_int32[::1] k,
               const npy_intp[:, ::1] _indices_k1d,
               const npy_intp[::1] _cstrides):
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

        # shifted indices into the data array
        npy_intp[::1] idx_c = np.ones(ndim, dtype=np.intp) * (-101)  # any sentinel would do, really
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
    volume = 1
    for d in range(ndim):
        volume *= k[d] + 1

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
            td = _t[d, :len_t[d]]
            xd = xv[d]
            kd = k[d]

            # get the location of x[d] in t[d]
            i[d] = _find_interval(&td[0], td.shape[0], kd, xd, kd, True)

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
                idx_cflat += idx_c[d] * _cstrides[d]

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

