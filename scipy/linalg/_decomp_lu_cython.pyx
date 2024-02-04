# cython: language_level=3

import cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from scipy.linalg.cython_lapack cimport sgetrf, dgetrf, cgetrf, zgetrf
from scipy.linalg._cythonized_array_utils cimport swap_c_and_f_layout

cimport numpy as cnp
cnp.import_array()

ctypedef fused lapack_t:
    cnp.float32_t
    cnp.float64_t
    cnp.complex64_t
    cnp.complex128_t


@cython.nonecheck(False)
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.initializedcheck(False)
cdef void lu_decompose(cnp.ndarray[lapack_t, ndim=2] a,
                       cnp.ndarray[lapack_t, ndim=2] lu,
                       int[::1] perm,
                       bint permute_l) noexcept:
    """LU decomposition and copy operations using ?getrf routines

    This function overwrites inputs. For interfacing LAPACK,
    it creates a memory buffer and copies into with F-order
    then swaps back to C order hence no need for dealing with
    Fortran arrays which are inconvenient.

    After the LU factorization, to minimize the amount of data
    copied, for rectangle arrays, and depending on the size,
    the smaller portion is copied out to U and the rest becomes
    either U or L. The logic is handled by the caller.

            ┌     ┐
            │ \ U │
            │  \  │       ┌                ┐
            │   \ │       │  \             │
         a= │     │  or   │   \      U     │
            │ L   │       │ L  \           │
            │     │       └                ┘
            └     ┘
             tall               wide
         (extract U)         (extract L)

    """
    cdef int m = a.shape[0], n = a.shape[1], mn = min(m, n)
    cdef cnp.npy_intp dims[2]
    cdef int info = 0, ind1, ind2, tmp_int
    cdef lapack_t *aa = <lapack_t *>cnp.PyArray_DATA(a)
    cdef lapack_t *bb
    cdef int *ipiv = <int*>PyMem_Malloc(m * sizeof(int))
    if not ipiv:
        raise MemoryError('scipy.linalg.lu failed to allocate '
                          'required memory.')
    dims[0] = m
    dims[1] = n

    if lapack_t is cnp.float32_t:
        b = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_FLOAT32)
        bb = <cnp.float32_t *>cnp.PyArray_DATA(b)
        swap_c_and_f_layout(aa, bb, m, n)
        sgetrf(&m, &n, bb, &m, ipiv, &info)
    elif lapack_t is cnp.float64_t:
        b = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_FLOAT64)
        bb = <cnp.float64_t *>cnp.PyArray_DATA(b)
        swap_c_and_f_layout(aa, bb, m, n)
        dgetrf(&m, &n, bb, &m, ipiv, &info)
    elif lapack_t is cnp.complex64_t:
        b = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_COMPLEX64)
        bb = <cnp.complex64_t *>cnp.PyArray_DATA(b)
        swap_c_and_f_layout(aa, bb, m, n)
        cgetrf(&m, &n, bb, &m, ipiv, &info)
    else:
        b = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_COMPLEX128)
        bb = <cnp.complex128_t *>cnp.PyArray_DATA(b)
        swap_c_and_f_layout(aa, bb, m, n)
        zgetrf(&m, &n, bb, &m, ipiv, &info)

    if info < 0:
        raise ValueError('scipy.linalg.lu has encountered an internal'
                         ' error in ?getrf routine with invalid value'
                         f' at {-info}-th parameter.')

    # Get the result back to C-contiguous layout and clean-up
    swap_c_and_f_layout(bb, aa, n, m)

    # Convert swaps on A to permutations on L since A = P @ L @ U
    try:
        # Basically we are following the cycles in ipiv
        # and swapping an "np.arange" array for the inverse perm.
        # Initialize perm
        for ind1 in range(m): perm[ind1] = ind1
        for ind1 in range(mn):
            tmp_int = perm[ipiv[ind1]-1]
            perm[ipiv[ind1]-1] = perm[ind1]
            perm[ind1] = tmp_int

        # convert iperm to perm into ipiv and store back into perm
        # as final. Solution without argsort : ipiv[perm] = np.arange(m)
        for ind1 in range(m):
            ipiv[perm[ind1]] = ind1
        for ind1 in range(m):
            perm[ind1] = ipiv[ind1]

    finally:
        PyMem_Free(ipiv)

    # Separation of L and U parts

    if m > n:  # tall, "a" holds bigger L
        # Extract upper right rectangle to lu
        for ind1 in range(mn):  # rows
            lu[ind1, ind1:mn] = a[ind1, ind1:mn]

        for ind1 in range(mn):
            a[ind1, ind1] = 1
            a[ind1, ind1+1:mn] = 0

    else:  # square or fat, "a" holds bigger U

        lu[0, 0] = 1
        for ind1 in range(1, mn):  # rows
            lu[ind1, :ind1] = a[ind1, :ind1]
            lu[ind1, ind1] = 1

        for ind2 in range(mn - 1):  # cols
            for ind1 in range(ind2+1, m):  # rows
                a[ind1, ind2] = 0

    if permute_l:
        # b still exists -> use it as temp array
        # we copy everything to b and pick back
        # rows from b as dictated by perm

        if m > n:
            b[:, :] = a[:, :]
            # memcpy(bb, &a[0, 0], m*mn*sizeof(lapack_t))
            for ind1 in range(m):
                if perm[ind1] == ind1:
                    continue
                else:
                    a[ind1, :] = b[perm[ind1], :]

        else:  # same but for lu array
            b[:mn, :mn] = lu[:, :]
            # memcpy(bb, &lu[0, 0], mn*n*sizeof(lapack_t))
            for ind1 in range(mn):
                if perm[ind1] == ind1:
                    continue
                else:
                    lu[ind1, :] = b[perm[ind1], :mn]


@cython.nonecheck(False)
@cython.initializedcheck(False)
def lu_dispatcher(a, u, piv, permute_l):
    if a.dtype.char == 'f':
        lu_decompose[cnp.float32_t](a, u, piv, permute_l)
    elif a.dtype.char == 'd':
        lu_decompose[cnp.float64_t](a, u, piv, permute_l)
    elif a.dtype.char == 'F':
        lu_decompose[cnp.complex64_t](a, u, piv, permute_l)
    elif a.dtype.char == 'D':
        lu_decompose[cnp.complex128_t](a, u, piv, permute_l)
    else:
        raise TypeError("Unsupported type given to lu_dispatcher")
