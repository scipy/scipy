# cython: boundscheck = False
# cython: wraparound = False
# cython: cdivision = True

from scipy.linalg cimport blas_pointers as blas

cdef inline bint is_contiguous(double[:,:] a, int axis):
    return (a.strides[axis] == sizeof(a[0,0]) or a.shape[axis] == 1)

cpdef float complex cdotc(float complex[:] cx, float complex[:] cy):
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
        int incy = cy.strides[0] // sizeof(cy[0])
    return blas.cdotc(&n, &cx[0], &incx, &cy[0], &incy)

cpdef float complex cdotu(float complex[:] cx, float complex[:] cy):
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
        int incy = cy.strides[0] // sizeof(cy[0])
    return blas.cdotu(&n, &cx[0], &incx, &cy[0], &incy)

cpdef double dasum(double[:] dx):
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
    return blas.dasum(&n, &dx[0], &incx)

cpdef double ddot(double[:] dx, double[:] dy):
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
        int incy = dy.strides[0] // sizeof(dy[0])
    return blas.ddot(&n, &dx[0], &incx, &dy[0], &incy)

cpdef int dgemm(double alpha, double[:,:] a, double[:,:] b, double beta,
                double[:,:] c) except -1:
    cdef:
        char *transa
        char *transb
        int m, n, k, lda, ldb, ldc
        double *a0=&a[0,0]
        double *b0=&b[0,0]
        double *c0=&c[0,0]
    # In the case that c is C contiguous, swap a and b and
    # swap whether or not each of them is transposed.
    # This can be done because a.dot(b) = b.T.dot(a.T).T.
    if is_contiguous(c, 1):
        if is_contiguous(a, 1):
            transb = 'n'
            ldb = (&a[1,0]) - a0 if a.shape[0] > 1 else 1
        elif is_contiguous(a, 0):
            transb = 't'
            ldb = (&a[0,1]) - a0 if a.shape[1] > 1 else 1
        else:
            raise ValueError("Input 'a' is neither C nor Fortran contiguous.")
        if is_contiguous(b, 1):
            transa = 'n'
            lda = (&b[1,0]) - b0 if b.shape[0] > 1 else 1
        elif is_contiguous(b, 0):
            transa = 't'
            lda = (&b[0,1]) - b0 if b.shape[1] > 1 else 1
        else:
            raise ValueError("Input 'b' is neither C nor Fortran contiguous.")
        k = b.shape[0]
        if k != a.shape[1]:
            raise ValueError("Shape mismatch in input arrays.")
        m = b.shape[1]
        n = a.shape[0]
        if n != c.shape[0] or m != c.shape[1]:
            raise ValueError("Output array does not have the correct shape.")
        ldc = (&c[1,0]) - c0 if c.shape[0] > 1 else 1
        blas.dgemm(transa, transb, &m, &n, &k, &alpha, b0, &lda, a0,
                   &ldb, &beta, c0, &ldc)
    elif is_contiguous(c, 0):
        if is_contiguous(a, 1):
            transa = 't'
            lda = (&a[1,0]) - a0 if a.shape[0] > 1 else 1
        elif is_contiguous(a, 0):
            transa = 'n'
            lda = (&a[0,1]) - a0 if a.shape[1] > 1 else 1
        else:
            raise ValueError("Input 'a' is neither C nor Fortran contiguous.")
        if is_contiguous(b, 1):
            transb = 't'
            ldb = (&b[1,0]) - b0 if b.shape[0] > 1 else 1
        elif is_contiguous(b, 0):
            transb = 'n'
            ldb = (&b[0,1]) - b0 if b.shape[1] > 1 else 1
        else:
            raise ValueError("Input 'b' is neither C nor Fortran contiguous.")
        m = a.shape[0]
        k = a.shape[1]
        if k != b.shape[0]:
            raise ValueError("Shape mismatch in input arrays.")
        n = b.shape[1]
        if m != c.shape[0] or n != c.shape[1]:
            raise ValueError("Output array does not have the correct shape.")
        ldc = (&c[0,1]) - c0 if c.shape[1] > 1 else 1
        blas.dgemm(transa, transb, &m, &n, &k, &alpha, a0, &lda, b0,
                   &ldb, &beta, c0, &ldc)
    else:
        raise ValueError("Input 'c' is neither C nor Fortran contiguous.")
    return 0

cpdef double dnrm2(double[:] x):
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return blas.dnrm2(&n, &x[0], &incx)

cpdef double dzasum(double complex[:] zx):
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
    return blas.dzasum(&n, &zx[0], &incx)

cpdef double dznrm2(double complex[:] x):
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return blas.dznrm2(&n, &x[0], &incx)

cpdef int icamax(float complex[:] cx):
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
    return blas.icamax(&n, &cx[0], &incx)

cpdef int idamax(double[:] dx):
    cdef:
        int n = dx.shape[0]
        int incx = dx.strides[0] // sizeof(dx[0])
    return blas.idamax(&n, &dx[0], &incx)

cpdef int isamax(float[:] sx):
    cdef:
        int n = sx.shape[0]
        int incx = sx.strides[0] // sizeof(sx[0])
    return blas.isamax(&n, &sx[0], &incx)

cpdef int izamax(double complex[:] zx):
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
    return blas.izamax(&n, &zx[0], &incx)

cpdef float sasum(float[:] sx):
    cdef:
        int n = sx.shape[0]
        int incx = sx.shape[0] // sizeof(sx[0])
    return blas.sasum(&n, &sx[0], &incx)

cpdef float scasum(float complex[:] cx):
    cdef:
        int n = cx.shape[0]
        int incx = cx.strides[0] // sizeof(cx[0])
    return blas.scasum(&n, &cx[0], &incx)

cpdef float scnrm2(float complex[:] x):
    cdef:
        int n = x.shape[0]
        int incx = x.strides[0] // sizeof(x[0])
    return blas.scnrm2(&n, &x[0], &incx)

cpdef float sdot(float[:] sx, float[:] sy):
    cdef:
        int n = sx.shape[0]
        int incx = sx.strides[0] // sizeof(sx[0])
        int incy = sy.strides[0] // sizeof(sy[0])
    return blas.sdot(&n, &sx[0], &incx, &sy[0], &incy)

cpdef float snrm2(float[:] x):
    cdef:
        int n = x.shape[0]
        int incx = x.shape[0] // sizeof(x[0])
    return blas.snrm2(&n, &x[0], &incx)

cpdef double complex zdotc(double complex[:] zx, double complex[:] zy):
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
        int incy = zy.strides[0] // sizeof(zy[0])
    return blas.zdotc(&n, &zx[0], &incx, &zy[0], &incy)

cpdef double complex zdotu(double complex[:] zx, double complex[:] zy):
    cdef:
        int n = zx.shape[0]
        int incx = zx.strides[0] // sizeof(zx[0])
        int incy = zy.strides[0] // sizeof(zy[0])
    return blas.zdotu(&n, &zx[0], &incx, &zy[0], &incy)
