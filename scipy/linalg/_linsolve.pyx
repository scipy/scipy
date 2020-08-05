cimport cython
from libc.stdlib cimport abs

ctypedef fused lapack_t:
    float
    double
    float complex
    double complex
ctypedef fused lapack_cz_t:
    float complex
    double complex
ctypedef fused lapack_sd_t:
    float
    double

ctypedef (int, int) (*f_tri_ptr)(lapack_t[:, ::1])
ctypedef (bint, bint) (*f_sym_ptr)(lapack_t[:, ::1])
ctypedef int (*f_posdiag_ptr)(lapack_t[:, ::1])


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (int, int) _is_tri_c(lapack_t[:, ::1]A) nogil:
    """
    Computes the lower and the upper bandwidth of a square
    matrix. If both bands are zero, it's diagonal, if either
    is zero then lower/upper triangular etc. If all entries
    are full then for nxn matrix the result is (n-1, n-1).

    This is the row-major layout version (C-contiguous).
    """
    cdef size_t n = A.shape[0]
    cdef int r, c, c_nnz, lower_band, upper_band
    lower_band = 0
    upper_band = 0
    for r in xrange(n):
        for c in xrange(n):
            if (A[r, c] != 0) and abs(r - c) > lower_band:
                # Only bother if outside the existing band:
                if (r > c):
                    lower_band = r-c
                else:
                    upper_band = c-r

    return lower_band, upper_band


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (bint, bint) _is_sym_her_c(lapack_t[:, ::1]A) nogil:
    """
    Checks for the symmetry or Hermitianness for a square matrix.
    For reals it doesn't matter however for complex matrices,
    symmetry can still be exploited even though the matrix is not
    Hermitian.
    """
    cdef size_t n = A.shape[0]
    cdef int r, c
    cdef lapack_t x, y
    cdef bint symmetric = 1
    cdef bint hermitian = 1

    # Real branch
    if lapack_t in lapack_sd_t:
        for r in xrange(n):
            for c in xrange(n):
                if r == c:
                    continue
                else:
                    if A[r, c] != A[c, r]:
                        return 0, 0
        return 1, 1
    # Complex branch
    else:
        for r in xrange(n):
            for c in xrange(n):
                # look-up once
                x, y = A[r, c], A[c, r]
                # once caught no need to keep checking continue for sym
                if hermitian:
                    if x != y.conjugate():
                        hermitian = 0

                if x != y:
                    symmetric = 0

                if not (symmetric or hermitian):
                    return 0, 0

        return hermitian, symmetric


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline int _is_posdef_diag_c(lapack_t[:, ::1]A) nogil:
    """
    One of the easy checks for positive/negative definiteness is that the
    diagonal should all be +/-. Here we check the sign of the first element
    and then compare with the rest and early exit if not.
    """
    cdef size_t n = A.shape[0]
    cdef int d
    cdef bint pos

    # Real branch
    if lapack_t in lapack_sd_t:
        if A[0, 0] == 0:
            return 0
        pos = A[0, 0] > 0

        for d in xrange(1, n):
            if (A[d, d] > 0) ^ pos:
                return 0
    else:
        if A[0, 0].real == 0:
            return 0
        pos = A[0, 0].real > 0

        for d in xrange(1, n):
            if (A[d, d].real > 0) ^ pos:
                return 0

    return 1 if pos else 2


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef (int, int) _check_for_structure(lapack_t[:, ::1]A):
    """
    Returns (family, member) structure information given a square C-contiguous
    array.

    Family
    ------
    0 : generic - no structure found
        0: dummy member
    1 : Triangular
        0: diagonal
        1: lower
        2: upper
        3: lower bidiagonal
        4: upper bidiagonal
    2 : Hessenberg
        0: lower
        1: upper
        2: tridiagonal
    3 : Hermitian
        0: non-positive definite
        1: maybe positive definite (survived the simple checks)
        2: maybe negative definite (survived the simple checks)
    4 : complex symmetric
        0: dummy member


    """
    cdef int lb, ub
    cdef f_tri_ptr tri_bands
    cdef f_sym_ptr sym_check
    cdef f_posdiag_ptr posdiag_check

    # Function pointers
    tri_bands = &_is_tri_c[lapack_t]
    sym_check = &_is_sym_her_c[lapack_t]
    posdiag_check = &_is_posdef_diag_c[lapack_t]

    # Triangular/Hess family
    lb, ub = tri_bands(A)
    if lb == 0 and ub == 0:
        return 1, 0  # diagonal
    elif lb == 0:
        if ub == 1:
            return 1, 4  # upper bidiagonal
        else:
            return 1, 1  # upper triangular
    elif ub == 0:
        if lb == 1:
            return 1, 3  # lower bidiagonal
        else:
            return 1, 2  # lower triangular

    if lb == 1 and ub == 1:
        return 2, 2  # tridiagonal
    elif lb == 1:
        return 2, 0  # lower hessenberg
    elif ub == 1:
        return 2, 1  # upper hessenberg

    # Sym/Her/PosDef family
    her, sym = sym_check(A)
    if her:
        return 3, posdiag_check(A)
    if sym:
        return 4, 0

    # No known structure return generic
    return 0, 0

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _check_cholesky_success(lapack_t[:, ::1]A):
    """
    When solving a hermitian linear system, attempting
    to factorize a square matrix often fails quite quickly
    hence the penalty is negligible if the array is not
    pos/neg definite and the benefit is noteworthy if it succeeds.

    Here we try to chol-factorize a hermitian matrix.
    LAPACK ?poXXX routines use one of the hermitian sides
    of the array including the diagonal. Here we use one
    side for this experiment (if succeeds we keep using
    that side). If fails we restore the diagonal and use
    the other side for later usage. Requires extra linear
    memory access and storage for the diagonal.

    We use the lower side for the experiment and upper for
    the failure. Note that the LAPACK functions use Fortran
    layout.
    """
    cdef int info
    cdef char* uplo = 'U'
    cdef int n = <int>a.shape[0]

    if lapack_t is float:
        spotrf(uplo, &n, &a[0, 0], &n, &info)
    elif lapack_t is double:
        dpotrf(uplo, &n, &a[0, 0], &n, &info)
    elif lapack_t is float_complex:
        cpotrf(uplo, &n, &a[0, 0], &n, &info)
    else:
        zpotrf(uplo, &n, &a[0, 0], &n, &info)

    return info


@cython.boundscheck(False)
@cython.wraparound(False)
cdef (lapack_t, lapack_t) _validate_input_arrays(lapack_t[:, :]A,
                                                 lapack_t[:, :]B):
    """
    All LAPACK funcs work with Fortran memory layout, but most use cases have
    C layout since NumPy's defaults. Hence if the data is contiguous in memory
    we work with transposed C arrays and if not we make a copy. Similarly,
    if the dtype is not exactly a member of lapack_t, we have to make a copy
    due to LAPACK specs.

    This provides further surgery options as opposed to letting f2py do the
    memory handling and array copies.
    """


