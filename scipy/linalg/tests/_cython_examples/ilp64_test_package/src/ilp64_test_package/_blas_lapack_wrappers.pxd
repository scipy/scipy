# Internal .pxd declaring wrapper functions with `int` parameters.
# This simulates the pattern used by scikit-learn and statsmodels,
# where the existing internal API uses `int` and must be adapted
# to work with ILP64 builds where BLAS/LAPACK expect 64-bit integers.

cdef double _ddot(int n, const double *x, int incx,
                  const double *y, int incy) noexcept nogil

cdef void _daxpy(int n, double alpha, const double *x, int incx,
                 double *y, int incy) noexcept nogil

cdef double _dnrm2(int n, const double *x, int incx) noexcept nogil

cdef void _dgemm(char *transa, char *transb, int m, int n, int k,
                 double alpha, const double *a, int lda,
                 const double *b, int ldb,
                 double beta, double *c, int ldc) noexcept nogil

cdef int _dgetrf(int m, int n, double *a, int lda,
                 int *ipiv) noexcept nogil

cpdef int get_blas_int_size()
