"""
Routines for updating QR decompositions

.. versionadded: 0.16.0

"""
# 
# Copyright (C) 2014 Eric Moore
#
# A few references for Updating QR factorizations:
#
# 1. Golub, G. H. & Loan, C. F. van V. Matrix Computations, 3rd Ed.
#    (Johns Hopkins University Press, 1996).
#
# 2. Daniel, J. W., Gragg, W. B., Kaufman, L. & Stewart, G. W.
#    Reorthogonalization and stable algorithms for updating the Gram-Schmidt QR
#    factorization. Math. Comput. 30, 772795 (1976).
#
# 3. Gill, P. E., Golub, G. H., Murray, W. & Saunders, M. A. Methods for
#    modifying matrix factorizations. Math. Comp. 28, 505535 (1974).
#
# 4. Hammarling, S. & Lucas, C. Updating the QR factorization and the least
#    squares problem. 1-73 (The University of Manchester, 2008).
#    at <http://eprints.ma.man.ac.uk/1192/>
# 

cimport cython
cimport libc.stdlib
cimport libc.limits as limits
from libc.string cimport memset
cimport numpy as cnp

cdef int NPY_ANYORDER = 0  
cdef int MEMORY_ERROR = limits.INT_MAX

# These are commented out in the numpy support we cimported above.
# Here I have declared them with as taking void* instead of PyArrayDescr
# and object. In this file, only NULL is passed to these parameters.
cdef extern from *:
    cnp.ndarray PyArray_CheckFromAny(object, void*, int, int, int, void*)
    cnp.ndarray PyArray_FromArray(cnp.ndarray, void*, int)

from scipy.linalg cimport blas_pointers
from scipy.linalg cimport lapack_pointers

import numpy as np

#------------------------------------------------------------------------------
# These are a set of fused type wrappers around the BLAS and LAPACK calls used. 
#------------------------------------------------------------------------------

ctypedef float complex float_complex
ctypedef double complex double_complex
ctypedef fused blas_t:
    float
    double
    float_complex
    double_complex

cdef inline blas_t* index2(blas_t* a, int* as, int i, int j) nogil:
    return a + i*as[0] + j*as[1]

cdef inline blas_t* index1(blas_t* a, int* as, int i) nogil:
    return a + i*as[0]

cdef inline blas_t* row(blas_t* a, int* as, int i) nogil:
    return a + i*as[0]

cdef inline blas_t* col(blas_t* a, int* as, int j) nogil:
    return a + j*as[1]
    
cdef inline void copy(int n, blas_t* x, int incx, blas_t* y, int incy) nogil:
    if blas_t is float:
        blas_pointers.scopy(&n, x, &incx, y, &incy) 
    elif blas_t is double:
        blas_pointers.dcopy(&n, x, &incx, y, &incy)
    elif blas_t is float_complex:
        blas_pointers.ccopy(&n, x, &incx, y, &incy)
    else:
        blas_pointers.zcopy(&n, x, &incx, y, &incy)

cdef inline void swap(int n, blas_t* x, int incx, blas_t* y, int incy) nogil:
    if blas_t is float:
        blas_pointers.sswap(&n, x, &incx, y, &incy) 
    elif blas_t is double:
        blas_pointers.dswap(&n, x, &incx, y, &incy)
    elif blas_t is float_complex:
        blas_pointers.cswap(&n, x, &incx, y, &incy)
    else:
        blas_pointers.zswap(&n, x, &incx, y, &incy)

cdef inline void axpy(int n, blas_t a, blas_t* x, int incx,
                      blas_t* y, int incy) nogil:
    if blas_t is float:
        blas_pointers.saxpy(&n, &a, x, &incx, y, &incy)
    elif blas_t is double:
        blas_pointers.daxpy(&n, &a, x, &incx, y, &incy)
    elif blas_t is float_complex:
        blas_pointers.caxpy(&n, &a, x, &incx, y, &incy)
    else:
        blas_pointers.zaxpy(&n, &a, x, &incx, y, &incy)

cdef inline void lartg(blas_t* a, blas_t* b, blas_t* c, blas_t* s) nogil:
    cdef blas_t g
    if blas_t is float:
        lapack_pointers.slartg(a, b, c, s, &g)
    elif blas_t is double:
        lapack_pointers.dlartg(a, b, c, s, &g)
    elif blas_t is float_complex:
        lapack_pointers.clartg(a, b, <float*>c, s, &g)
    else:
        lapack_pointers.zlartg(a, b, <double*>c, s, &g)
    # make this function more like the BLAS drotg
    a[0] = g
    b[0] = 0

cdef inline void rot(int n, blas_t* x, int incx, blas_t* y, int incy,
                     blas_t c, blas_t s) nogil:
    if blas_t is float:
        blas_pointers.srot(&n, x, &incx, y, &incy, &c, &s) 
    elif blas_t is double:
        blas_pointers.drot(&n, x, &incx, y, &incy, &c, &s)
    elif blas_t is float_complex:
        lapack_pointers.crot(&n, x, &incx, y, &incy, <float*>&c, &s)
    else:
        lapack_pointers.zrot(&n, x, &incx, y, &incy, <double*>&c, &s)

cdef inline void larfg(int n, blas_t* alpha, blas_t* x, int incx,
                       blas_t* tau) nogil:
    if blas_t is float:
        lapack_pointers.slarfg(&n, alpha, x, &incx, tau)
    elif blas_t is double:
        lapack_pointers.dlarfg(&n, alpha, x, &incx, tau)
    elif blas_t is float_complex:
        lapack_pointers.clarfg(&n, alpha, x, &incx, tau)
    else:
        lapack_pointers.zlarfg(&n, alpha, x, &incx, tau)

cdef inline void larf(char* side, int m, int n, blas_t* v, int incv, blas_t tau,
    blas_t* c, int ldc, blas_t* work) nogil:
    if blas_t is float:
        lapack_pointers.slarf(side, &m, &n, v, &incv, &tau, c, &ldc, work)
    elif blas_t is double:
        lapack_pointers.dlarf(side, &m, &n, v, &incv, &tau, c, &ldc, work)
    elif blas_t is float_complex:
        lapack_pointers.clarf(side, &m, &n, v, &incv, &tau, c, &ldc, work)
    else:
        lapack_pointers.zlarf(side, &m, &n, v, &incv, &tau, c, &ldc, work)

#------------------------------------------------------------------------------
# Utility routines
#------------------------------------------------------------------------------

cdef cnp.ndarray PyArray_FromArraySafe(cnp.ndarray arr, void* newtype, int flags):
    """In Numpy 1.5.1, Use of the NPY_F_CONTIGUOUS flag is broken when used
    with a 1D input array.  The work around is to pass NPY_C_CONTIGUOUS since
    for 1D arrays these are equivalent.  This is Numpy's gh-2287, fixed in
    9b8ff38.

    FIXME: Is it worth only applying this for numpy 1.5.1? 
    """
    if arr.ndim == 1 and (flags & cnp.NPY_F_CONTIGUOUS):
        flags = flags & ~cnp.NPY_F_CONTIGUOUS
        flags |= cnp.NPY_C_CONTIGUOUS
    return PyArray_FromArray(arr, newtype, flags)

cdef void blas_t_conj(int n, blas_t* x, int* xs) nogil:
    cdef int j
    if blas_t is float_complex or blas_t is double_complex:
        for j in range(n):
            index1(x, xs, j)[0] = index1(x, xs, j)[0].conjugate()

#------------------------------------------------------------------------------
# QR update routines start here.
#------------------------------------------------------------------------------
cdef void qr_row_delete(int m, int n, blas_t* q, int* qs, blas_t* r, int* rs,
                        int k) nogil:
    cdef int j
    cdef blas_t c, s
    
    if k != 0:
        for j in range(k, 0, -1):
            swap(m, row(q, qs, j), qs[1], row(q, qs, j-1), qs[1])

    blas_t_conj(m, row(q, qs, 0), &qs[1])

    for j in range(m-2, -1, -1):
        lartg(index2(q, qs, 0, j), index2(q, qs, 0, j+1), &c, &s)

        # update the columns of r if there is a nonzero row.
        if j < n:
            rot(n-j, index2(r, rs, j, j), rs[1], index2(r, rs, j+1, j), rs[1],
                    c, s)

        # update the rows of q
        rot(m-1, index2(q, qs, 1, j), qs[0], index2(q, qs, 1, j+1), qs[0],
                c, s.conjugate())

cdef void qr_block_row_delete(int m, int n, blas_t* q, int* qs,
                              blas_t* r, int* rs, int k, int p) nogil:
    cdef int i, j
    cdef blas_t c,s
    cdef blas_t* W
    cdef int* ws

    if k != 0:
        for j in range(k, 0, -1):
            swap(m, row(q, qs, j+p-1), qs[1], row(q, qs, j-1), qs[1])
    
    # W is the block of rows to be removed from q, has shape, (p,m)
    W = q
    ws = qs

    for j in range(p):
        blas_t_conj(m, row(W, ws, j), &ws[1])

    for i in range(p):
        for j in range(m-2, i-1, -1):
            lartg(index2(W, ws, i, j), index2(W, ws, i, j+1), &c, &s)
            
            # update W
            if i+1 < p:
                rot(p-i-1, index2(W, ws, i+1, j), ws[0],
                    index2(W, ws, i+1, j+1), ws[0], c, s)

            # update r if there is a nonzero row.
            if j-i < n:
                rot(n-j+i, index2(r, rs, j, j-i), rs[1],
                    index2(r, rs, j+1, j-i), rs[1], c, s)

            # update q
            rot(m-p, index2(q, qs, p, j), qs[0], index2(q, qs, p, j+1), qs[0],
                c, s.conjugate())

cdef void qr_col_delete(int m, int n, blas_t* q, int* qs, blas_t* r, int* rs,
                        int k) nogil:
    cdef int j

    for j in range(k, n-1):
        copy(m, col(r, rs, j+1), rs[0], col(r, rs, j), rs[0])

    hessenberg_qr(m, n-1, q, qs, r, rs, k)

cdef int qr_block_col_delete(int m, int n, blas_t* q, int* qs,
                              blas_t* r, int* rs, int k, int p) nogil:
    cdef int j
    cdef blas_t* work
    cdef int worksize = max(m, n)

    work = <blas_t*>libc.stdlib.malloc(worksize*sizeof(blas_t))
    if not work:
        return MEMORY_ERROR

    # move the columns to removed to the end
    for j in range(k, n-p):
        copy(m, col(r, rs, j+p), rs[0], col(r, rs, j), rs[0])

    p_subdiag_qr(m, n-p, q, qs, r, rs, k, p, work)

    libc.stdlib.free(work)
    return 0

cdef void hessenberg_qr(int m, int n, blas_t* q, int* qs, blas_t* r, int* rs,
                        int k) nogil:
    """Reduce an upper hessenberg matrix r, to upper triangluar,
    starting in row j.  Apply these transformation to q as well.
    """
    cdef int j
    cdef blas_t c, s
    cdef int limit = min(m-1, n)

    for j in range(k, limit):
        lartg(index2(r, rs, j, j), index2(r, rs, j+1, j), &c, &s)

        # update the rest of r
        if j+1 < m:
            rot(n-j-1, index2(r, rs, j, j+1), rs[1],
                    index2(r, rs, j+1, j+1), rs[1], c, s) 

        # update q
        rot(m, col(q, qs, j), qs[0], col(q, qs, j+1), qs[0], c, s.conjugate())

cdef void p_subdiag_qr(int m, int n, blas_t* q, int* qs, blas_t* r, int* rs,
                       int k, int p, blas_t* work) nogil:
    """ Reduce a matrix r to upper triangular form by eliminating the lower p 
        subdiagionals using reflectors.
        
        q and r must be fortran order here, with work at least max(m,n) long.
    """
    cdef int j
    cdef int last
    cdef blas_t tau
    cdef blas_t rjj 
    cdef int limit = min(m-1, n)
    cdef char* sideR = 'R'
    cdef char* sideL = 'L'

    # R now has p subdiagonal values to be removed starting from col k.
    for j in range(k, limit):
        # length of the reflector
        last = min(p+1, m-j)
        rjj = index2(r, rs, j, j)[0]
        larfg(last, &rjj, index2(r, rs, j+1, j), rs[0], &tau)
        index2(r, rs, j, j)[0] = 1

        # apply the reflector to r if necessary
        if j+1 < n:
            larf(sideL, last, n-j-1, index2(r, rs, j, j), rs[0],
                    tau.conjugate(), index2(r, rs, j, j+1), rs[1], work)

        # apply the reflector to q
        larf(sideR, m, last, index2(r, rs, j, j), rs[0], tau,
                index2(q, qs, 0, j), qs[1], work)

        # rezero the householder vector we no longer need.
        memset(index2(r, rs, j+1, j), 0, (last-1)*sizeof(blas_t))

        # restore the rjj element
        index2(r, rs, j, j)[0] = rjj

cdef validate_array(cnp.ndarray a):
    # here we check that a has positive strides and that its size is small
    # enough to fit in into an int, as BLAS/LAPACK require
    cdef bint copy = False
    cdef bint too_large = False
    cdef int j

    for j in range(a.ndim):
        if a.strides[j] <= 0 or \
                (a.strides[j] / a.descr.itemsize) >= limits.INT_MAX:
            copy = True
        if a.shape[j] >= limits.INT_MAX:
            raise ValueError('Input array to large for use with BLAS')

    if copy:
            return PyArray_FromArraySafe(a, NULL, cnp.NPY_F_CONTIGUOUS)
    return a

cdef validate_qr(object q0, object r0, bint overwrite_q, int q_order,
                 bint overwrite_r, int r_order):
    cdef cnp.ndarray Q
    cdef cnp.ndarray R
    cdef int typecode

    q_order |= cnp.NPY_BEHAVED_NS | cnp.NPY_ELEMENTSTRIDES  
    r_order |= cnp.NPY_BEHAVED_NS | cnp.NPY_ELEMENTSTRIDES

    if not overwrite_q:
        q_order |= cnp.NPY_ENSURECOPY

    if not overwrite_r:
        r_order |= cnp.NPY_ENSURECOPY

    # in the interests of giving better error messages take any number of
    # dimensions here.
    Q = PyArray_CheckFromAny(q0, NULL, 0, 0, q_order, NULL)
    R = PyArray_CheckFromAny(r0, NULL, 0, 0, r_order, NULL)

    if Q.ndim != 2 or R.ndim != 2:
        raise ValueError('Q and R must be 2d')

    typecode = cnp.PyArray_TYPE(Q)

    if typecode != cnp.PyArray_TYPE(R):
        raise ValueError('q and r must have the same type')

    if not (typecode == cnp.NPY_FLOAT or typecode == cnp.NPY_DOUBLE 
            or typecode == cnp.NPY_CFLOAT or typecode == cnp.NPY_CDOUBLE):
        raise ValueError('only floatingcomplex arrays supported')

    if Q.shape[1] != R.shape[0]:
        raise ValueError('Q and R do not have compatible shapes')

    # economic decomposition are not supported.
    if Q.shape[0] != Q.shape[1]:
        raise ValueError('economic mode decompositions are not supported.')

    Q = validate_array(Q)
    R = validate_array(R)

    return Q, R, typecode, Q.shape[0], R.shape[1]
 
cdef void* extract(cnp.ndarray a, int* as):
    if a.ndim == 2:
        as[0] = a.strides[0] / cnp.PyArray_ITEMSIZE(a)
        as[1] = a.strides[1] / cnp.PyArray_ITEMSIZE(a)
    elif a.ndim == 1:
        as[0] = a.strides[0] / cnp.PyArray_ITEMSIZE(a)
        as[1] = 0
    return cnp.PyArray_DATA(a)

@cython.embedsignature(True)
def qr_delete(Q, R, k, p=1, which='row', overwrite_qr=True):
    """QR downdate on row or column deletions

    If ``A = Q R`` is the qr factorization of A, return the qr factorization
    of `A` where `p` rows or columns have been removed starting at row or
    column `k`.

    Parameters
    ----------
    Q : array_like
        Unitary/orthogonal matrix from QR decomposition.
    R : array_like
        Upper trianglar matrix from QR decomposition.
    k : int
        index of the first row or column to delete.
    p : int, optional
        number of rows or columns to delete, defaults to 1.
    which: {'row', 'col'}, optional
        Determines if rows or columns will be deleted, defaults to 'row'
    overwrite_qr : bool, optional
        If True, consume Q and R, overwriting their contents with their
        downdated versions, and returning approriately sized views.  
        Defaults to True.

    Returns
    -------
    Q1 : ndarray
        Updated unitary/orthogonal factor
    R1 : ndarray
        Updated upper triangular factor

    Notes
    -----
    This routine does not guarantee that the diagonal entries of `R1` are
    positive.

    .. versionadded:: 0.16.0

    Examples
    --------
    >>> from scipy import linalg
    >>> a = np.array([[  3.,  -2.,  -2.],
                      [  6.,  -9.,  -3.],
                      [ -3.,  10.,   1.],
                      [  6.,  -7.,   4.],
                      [  7.,   8.,  -6.]])
    >>> q, r = linalg.qr(a)

    Given this q, r decomposition, update q and r when 2 rows are removed.

    >>> q1, r1 = linalg.qr_delete(q, r, 2, 2, 'row', False)
    >>> q1
    array([[ 0.30942637,  0.15347579,  0.93845645],
           [ 0.61885275,  0.71680171, -0.32127338],
           [ 0.72199487, -0.68017681, -0.12681844]])
    >>> r1
    array([[  9.69535971,  -0.4125685 ,  -6.80738023],
           [  0.        , -12.19958144,   1.62370412],
           [  0.        ,   0.        ,  -0.15218213]])

    The update is equivalent, but faster than the following.

    >>> a1 = np.delete(a, slice(2,4), 0)
    >>> a1
    array([[ 3., -2., -2.],
           [ 6., -9., -3.],
           [ 7.,  8., -6.]])
    >>> q_direct, r_direct = linalg.qr(a1)

    Check that we have equivalent results:

    >>> np.dot(q1, r1)
    array([[ 3., -2., -2.],
           [ 6., -9., -3.],
           [ 7.,  8., -6.]])
    >>> np.allclose(np.dot(q1, r1), a1)
    True

    And the updated Q is still unitary:

    >>> np.allclose(np.dot(q1.T, q1), np.eye(3))
    True

    """
    cdef cnp.ndarray q1, r1
    cdef int k1 = k
    cdef int p1 = p
    cdef int typecode, m, n, info
    cdef int qs[2]
    cdef int rs[2]

    if which == 'row':
        q1, r1, typecode, m, n = validate_qr(Q, R, overwrite_qr, NPY_ANYORDER,
                overwrite_qr, NPY_ANYORDER)
        if not (-m <= k1 < m):
            raise ValueError('k is not a valid index')
        if k1 < 0:
            k1 += m
        if k1 + p1 > m or p1 <= 0: 
            raise ValueError('p out of range')
        if p1 == 1:
            if typecode == cnp.NPY_FLOAT:
                qr_row_delete(m, n, <float*>extract(q1, qs), qs,
                    <float*>extract(r1, rs), rs, k1)
            elif typecode == cnp.NPY_DOUBLE:
                qr_row_delete(m, n, <double*>extract(q1, qs), qs,
                    <double*>extract(r1, rs), rs, k1)
            elif typecode == cnp.NPY_CFLOAT:
                qr_row_delete(m, n, <float_complex*>extract(q1, qs), qs,
                    <float_complex*>extract(r1, rs), rs, k1)
            else:  # cnp.NPY_CDOUBLE:
                qr_row_delete(m, n, <double_complex*>extract(q1, qs), qs,
                    <double_complex*>extract(r1, rs), rs, k1)
        else:
            if typecode == cnp.NPY_FLOAT:
                qr_block_row_delete(m, n, <float*>extract(q1, qs), qs,
                    <float*>extract(r1, rs), rs, k1, p1)
            elif typecode == cnp.NPY_DOUBLE:
                qr_block_row_delete(m, n, <double*>extract(q1, qs), qs,
                    <double*>extract(r1, rs), rs, k1, p1)
            elif typecode == cnp.NPY_CFLOAT:
                qr_block_row_delete(m, n, <float_complex*>extract(q1, qs), qs,
                    <float_complex*>extract(r1, rs), rs, k1, p1)
            else:  # cnp.NPY_CDOUBLE:
                qr_block_row_delete(m, n, <double_complex*>extract(q1, qs), qs,
                    <double_complex*>extract(r1, rs), rs, k1, p1)
        return q1[p:, p:], r1[p:, :]
    elif which == 'col':
        if p1 > 1:
            q1, r1, typecode, m, n = validate_qr(Q, R, overwrite_qr,
                    cnp.NPY_F_CONTIGUOUS, overwrite_qr, cnp.NPY_F_CONTIGUOUS)
        else:
            q1, r1, typecode, m, n = validate_qr(Q, R, overwrite_qr,
                    NPY_ANYORDER, overwrite_qr, NPY_ANYORDER)
        if not (-n <= k1 < n):
            raise ValueError('k is not a valid index')
        if k1 < 0:
            k1 += n
        if k1 + p1 > n or p1 <= 0:
            raise ValueError('p out of range')

        if p1 == 1:
            if typecode == cnp.NPY_FLOAT:
                qr_col_delete(m, n, <float*>extract(q1, qs), qs, 
                    <float*>extract(r1, rs), rs, k1)
            elif typecode == cnp.NPY_DOUBLE:
                qr_col_delete(m, n, <double*>extract(q1, qs), qs, 
                    <double*>extract(r1, rs), rs, k1)
            elif typecode == cnp.NPY_CFLOAT:
                qr_col_delete(m, n, <float_complex*>extract(q1, qs), qs, 
                    <float_complex*>extract(r1, rs), rs, k1)
            else:  # cnp.NPY_CDOUBLE:
                qr_col_delete(m, n, <double_complex*>extract(q1, qs), qs, 
                    <double_complex*>extract(r1, rs), rs, k1)
        else:
            if typecode == cnp.NPY_FLOAT:
                info = qr_block_col_delete(m, n, <float*>extract(q1, qs), qs,
                    <float*>extract(r1, rs), rs, k1, p1)
            elif typecode == cnp.NPY_DOUBLE:
                info = qr_block_col_delete(m, n, <double*>extract(q1, qs), qs,
                    <double*>extract(r1, rs), rs, k1, p1)
            elif typecode == cnp.NPY_CFLOAT:
                info = qr_block_col_delete(m, n, <float_complex*>extract(q1, qs), qs,
                    <float_complex*>extract(r1, rs), rs, k1, p1)
            else:  # cnp.NPY_CDOUBLE:
                info = qr_block_col_delete(m, n, <double_complex*>extract(q1, qs), qs,
                    <double_complex*>extract(r1, rs), rs, k1, p1)
            if info == MEMORY_ERROR:
                raise MemoryError('malloc failed')
        return q1, r1[:,:-p]
    else:
        raise ValueError("which must be either 'row' or 'col'")

cnp.import_array()

