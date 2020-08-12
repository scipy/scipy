"""
Polyalgorithm functions for scipy.linalg.solve()

Solving linear system equations are central to the numerical linear algebra
routines. And thus, any performance benefit introduced by the solver directly
impacts upstream code elsewhere. The established convention to solve systems
with unknown structure is to perform gradually more expensive checks and branch
off to the structured solvers if possible, finally fallback to generic solvers
if no structure is found. This type of "if-this-use-that" algorithms are called
polyalgorithms. Many numerical algebra tools already use this type of branching
(matlab's backslash operator possibly being the most well-known). It might be
counter-intuitive to waste time on such checks at first, however, it has to be
noted that linear solvers often have O(n**3) complexity. Many checks we
perform are of O(n**2) hence the worst case, "no structure after all checks",
doesn't suffer a huge delay, in turn if structure is found, significant speed
benefits are obtained.

Here, Cythonized versions of these checks and branching logic is provided;
since we use BLAS/LAPACK backend, the arrays have to be prepared strictly to
have only single/double precision dtypes and Fortran-contiguous (column-major)
memory layout to be used in LAPACK calls. Thus the entrance of the solve(A, b)
regularizes the inputs and make copies if not possible. If the datatype is
matched, C-contiguous arrays are transposed and used without any copies.

Next, we start checking the bandwidth of the arrays looping over rows/columns.
If a particular form is detected, then checks stop and relevant solver is used.
If no triangular form is discovered then symmetry/hermitian properties are
checked. If further obvious properties are satisfied (diagonal entries of a
sign-definite matrix should have the same sign and be nonzero) then Cholesky
factorization is tested. Depending on the result relevant solvers are used.
It should be noted that even though Cholesky factorization is comparable to
other solvers in time-complexity, it often fails very quickly in practice and
the worst-case scenarios hardly ever encountered and in fact require special
matrix constructions.

If no structure is detected, standard generic linsolve routines are used.
"""

cimport cython
from libc.stdlib cimport abs
cimport numpy as cnp

ctypedef fused valid_type_t:
    cnp.float32
    cnp.float64
    cnp.complex


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

# The numpy facilities for type casting checks are too slow for small-sized
# arrays and eats away the time-budget for the checkups. Here we provide a
# precomputed version of the numpy.can_cast() table.

# output of np.typecodes['All']
cdef str casting_key = '?bhilqpBHILQPefdgFDGSUVOMm'

#                          f  d  F  D
casting_table = np.array([[1, 1, 1, 1],  # ?
                          [1, 1, 1, 1],  # b
                          [1, 1, 1, 1],  # h
                          [0, 1, 0, 1],  # i
                          [0, 1, 0, 1],  # l
                          [0, 1, 0, 1],  # q
                          [0, 1, 0, 1],  # p
                          [1, 1, 1, 1],  # B
                          [1, 1, 1, 1],  # H
                          [0, 1, 0, 1],  # I
                          [0, 1, 0, 1],  # L
                          [0, 1, 0, 1],  # Q
                          [0, 1, 0, 1],  # P
                          [1, 1, 1, 1],  # e
                          [1, 1, 1, 1],  # f
                          [0, 1, 0, 1],  # d
                          [0, 1, 0, 1],  # g
                          [0, 0, 1, 1],  # F
                          [0, 0, 0, 1],  # D
                          [0, 0, 0, 1],  # G
                          [0, 0, 0, 0],  # S
                          [0, 0, 0, 0],  # U
                          [0, 0, 0, 0],  # V
                          [0, 0, 0, 0],  # O
                          [0, 0, 0, 0],  # M
                          [0, 0, 0, 0]], # m
                         dtype=bool)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef inline int get_lapack_flavor(
    str a,
    str b,
    str K = casting_key,
    const cnp.npy_bool[:, ::1]C = casting_table
    ):
    """
    Given dtype characters, this function finds the smallest LAPACK compatible
    dtype index (counted for "s", "d", "c", "z"). If not found, returns -1.
    """
    cdef size_t ind
    cdef size_t rowa = K.find(a)
    cdef size_t rowb = K.find(b)
    if rowa == -1 or rowb == -1:
        raise TypeError('Unknown data type for linalg.solve()')

    for ind in range(4):
        if C[rowa, ind] and C[rowb, ind]:
            return ind

    return -1


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
cpdef (int, int) _check_for_structure_c(lapack_t[:, ::1]A):
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
cdef int _check_cholesky_success_c(lapack_t[:, ::1]a):
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






