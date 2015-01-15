ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

ctypedef int caxpy_t(int *n, c *ca, c *cx, int *incx, c *cy, int *incy) nogil
ctypedef int ccopy_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
ctypedef c cdotc_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
ctypedef c cdotu_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
ctypedef int cgemm_t(char *transa, char *transb, int *m, int *n, int *k, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
ctypedef int cgemv_t(char *trans, int *m, int *n, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
ctypedef int cgerc_t(int *m, int *n, c *alpha, c *x, int *incx, c *y, int *incy, c *a, int *lda) nogil
ctypedef int cgeru_t(int *m, int *n, c *alpha, c *x, int *incx, c *y, int *incy, c *a, int *lda) nogil
ctypedef int chemm_t(char *side, char *uplo, int *m, int *n, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
ctypedef int chemv_t(char *uplo, int *n, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
ctypedef int cher_t(char *uplo, int *n, s *alpha, c *x, int *incx, c *a, int *lda) nogil
ctypedef int cher2_t(char *uplo, int *n, c *alpha, c *x, int *incx, c *y, int *incy, c *a, int *lda) nogil
ctypedef int cher2k_t(char *uplo, char *trans, int *n, int *k, c *alpha, c *a, int *lda, c *b, int *ldb, s *beta, c *c, int *ldc) nogil
ctypedef int cherk_t(char *uplo, char *trans, int *n, int *k, s *alpha, c *a, int *lda, s *beta, c *c, int *ldc) nogil
ctypedef int crotg_t(c *ca, c *cb, s *c, c *s) nogil
ctypedef int cscal_t(int *n, c *ca, c *cx, int *incx) nogil
ctypedef int csrot_t(int *n, c *cx, int *incx, c *cy, int *incy, s *c, s *s) nogil
ctypedef int csscal_t(int *n, s *sa, c *cx, int *incx) nogil
ctypedef int cswap_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
ctypedef int csymm_t(char *side, char *uplo, int *m, int *n, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
ctypedef int csyr_t(char *uplo, int *n, c *alpha, c *x, int *incx, c *a, int *lda) nogil
ctypedef int csyr2k_t(char *uplo, char *trans, int *n, int *k, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
ctypedef int csyrk_t(char *uplo, char *trans, int *n, int *k, c *alpha, c *a, int *lda, c *beta, c *c, int *ldc) nogil
ctypedef int ctrmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, c *alpha, c *a, int *lda, c *b, int *ldb) nogil
ctypedef int ctrmv_t(char *uplo, char *trans, char *diag, int *n, c *a, int *lda, c *x, int *incx) nogil
ctypedef d dasum_t(int *n, d *dx, int *incx) nogil
ctypedef int daxpy_t(int *n, d *da, d *dx, int *incx, d *dy, int *incy) nogil
ctypedef int dcopy_t(int *n, d *dx, int *incx, d *dy, int *incy) nogil
ctypedef d ddot_t(int *n, d *dx, int *incx, d *dy, int *incy) nogil
ctypedef int dgemm_t(char *transa, char *transb, int *m, int *n, int *k, d *alpha, d *a, int *lda, d *b, int *ldb, d *beta, d *c, int *ldc) nogil
ctypedef int dgemv_t(char *trans, int *m, int *n, d *alpha, d *a, int *lda, d *x, int *incx, d *beta, d *y, int *incy) nogil
ctypedef int dger_t(int *m, int *n, d *alpha, d *x, int *incx, d *y, int *incy, d *a, int *lda) nogil
ctypedef d dnrm2_t(int *n, d *x, int *incx) nogil
ctypedef int drot_t(int *n, d *dx, int *incx, d *dy, int *incy, d *c, d *s) nogil
ctypedef int drotg_t(d *da, d *db, d *c, d *s) nogil
ctypedef int drotm_t(int *n, d *dx, int *incx, d *dy, int *incy, d *dparam) nogil
ctypedef int drotmg_t(d *dd1, d *dd2, d *dx1, d *dy1, d *dparam) nogil
ctypedef int dscal_t(int *n, d *da, d *dx, int *incx) nogil
ctypedef int dswap_t(int *n, d *dx, int *incx, d *dy, int *incy) nogil
ctypedef int dsymm_t(char *side, char *uplo, int *m, int *n, d *alpha, d *a, int *lda, d *b, int *ldb, d *beta, d *c, int *ldc) nogil
ctypedef int dsymv_t(char *uplo, int *n, d *alpha, d *a, int *lda, d *x, int *incx, d *beta, d *y, int *incy) nogil
ctypedef int dsyr_t(char *uplo, int *n, d *alpha, d *x, int *incx, d *a, int *lda) nogil
ctypedef int dsyr2_t(char *uplo, int *n, d *alpha, d *x, int *incx, d *y, int *incy, d *a, int *lda) nogil
ctypedef int dsyr2k_t(char *uplo, char *trans, int *n, int *k, d *alpha, d *a, int *lda, d *b, int *ldb, d *beta, d *c, int *ldc) nogil
ctypedef int dsyrk_t(char *uplo, char *trans, int *n, int *k, d *alpha, d *a, int *lda, d *beta, d *c, int *ldc) nogil
ctypedef int dtrmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, d *alpha, d *a, int *lda, d *b, int *ldb) nogil
ctypedef int dtrmv_t(char *uplo, char *trans, char *diag, int *n, d *a, int *lda, d *x, int *incx) nogil
ctypedef d dzasum_t(int *n, z *zx, int *incx) nogil
ctypedef d dznrm2_t(int *n, z *x, int *incx) nogil
ctypedef int icamax_t(int *n, c *cx, int *incx) nogil
ctypedef int idamax_t(int *n, d *dx, int *incx) nogil
ctypedef int isamax_t(int *n, s *sx, int *incx) nogil
ctypedef int izamax_t(int *n, z *zx, int *incx) nogil
ctypedef s sasum_t(int *n, s *sx, int *incx) nogil
ctypedef int saxpy_t(int *n, s *sa, s *sx, int *incx, s *sy, int *incy) nogil
ctypedef s scasum_t(int *n, c *cx, int *incx) nogil
ctypedef s scnrm2_t(int *n, c *x, int *incx) nogil
ctypedef int scopy_t(int *n, s *sx, int *incx, s *sy, int *incy) nogil
ctypedef s sdot_t(int *n, s *sx, int *incx, s *sy, int *incy) nogil
ctypedef int sgemm_t(char *transa, char *transb, int *m, int *n, int *k, s *alpha, s *a, int *lda, s *b, int *ldb, s *beta, s *c, int *ldc) nogil
ctypedef int sgemv_t(char *trans, int *m, int *n, s *alpha, s *a, int *lda, s *x, int *incx, s *beta, s *y, int *incy) nogil
ctypedef int sger_t(int *m, int *n, s *alpha, s *x, int *incx, s *y, int *incy, s *a, int *lda) nogil
ctypedef s snrm2_t(int *n, s *x, int *incx) nogil
ctypedef int srot_t(int *n, s *sx, int *incx, s *sy, int *incy, s *c, s *s) nogil
ctypedef int srotg_t(s *sa, s *sb, s *c, s *s) nogil
ctypedef int srotm_t(int *n, s *sx, int *incx, s *sy, int *incy, s *sparam) nogil
ctypedef int srotmg_t(s *sd1, s *sd2, s *sx1, s *sy1, s *sparam) nogil
ctypedef int sscal_t(int *n, s *sa, s *sx, int *incx) nogil
ctypedef int sswap_t(int *n, s *sx, int *incx, s *sy, int *incy) nogil
ctypedef int ssymm_t(char *side, char *uplo, int *m, int *n, s *alpha, s *a, int *lda, s *b, int *ldb, s *beta, s *c, int *ldc) nogil
ctypedef int ssymv_t(char *uplo, int *n, s *alpha, s *a, int *lda, s *x, int *incx, s *beta, s *y, int *incy) nogil
ctypedef int ssyr_t(char *uplo, int *n, s *alpha, s *x, int *incx, s *a, int *lda) nogil
ctypedef int ssyr2_t(char *uplo, int *n, s *alpha, s *x, int *incx, s *y, int *incy, s *a, int *lda) nogil
ctypedef int ssyr2k_t(char *uplo, char *trans, int *n, int *k, s *alpha, s *a, int *lda, s *b, int *ldb, s *beta, s *c, int *ldc) nogil
ctypedef int ssyrk_t(char *uplo, char *trans, int *n, int *k, s *alpha, s *a, int *lda, s *beta, s *c, int *ldc) nogil
ctypedef int strmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, s *alpha, s *a, int *lda, s *b, int *ldb) nogil
ctypedef int strmv_t(char *uplo, char *trans, char *diag, int *n, s *a, int *lda, s *x, int *incx) nogil
ctypedef int zaxpy_t(int *n, z *za, z *zx, int *incx, z *zy, int *incy) nogil
ctypedef int zcopy_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
ctypedef z zdotc_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
ctypedef z zdotu_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
ctypedef int zdrot_t(int *n, z *cx, int *incx, z *cy, int *incy, d *c, d *s) nogil
ctypedef int zdscal_t(int *n, d *da, z *zx, int *incx) nogil
ctypedef int zgemm_t(char *transa, char *transb, int *m, int *n, int *k, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
ctypedef int zgemv_t(char *trans, int *m, int *n, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
ctypedef int zgerc_t(int *m, int *n, z *alpha, z *x, int *incx, z *y, int *incy, z *a, int *lda) nogil
ctypedef int zgeru_t(int *m, int *n, z *alpha, z *x, int *incx, z *y, int *incy, z *a, int *lda) nogil
ctypedef int zhemm_t(char *side, char *uplo, int *m, int *n, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
ctypedef int zhemv_t(char *uplo, int *n, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
ctypedef int zher_t(char *uplo, int *n, d *alpha, z *x, int *incx, z *a, int *lda) nogil
ctypedef int zher2_t(char *uplo, int *n, z *alpha, z *x, int *incx, z *y, int *incy, z *a, int *lda) nogil
ctypedef int zher2k_t(char *uplo, char *trans, int *n, int *k, z *alpha, z *a, int *lda, z *b, int *ldb, d *beta, z *c, int *ldc) nogil
ctypedef int zherk_t(char *uplo, char *trans, int *n, int *k, d *alpha, z *a, int *lda, d *beta, z *c, int *ldc) nogil
ctypedef int zrotg_t(z *ca, z *cb, d *c, z *s) nogil
ctypedef int zscal_t(int *n, z *za, z *zx, int *incx) nogil
ctypedef int zswap_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
ctypedef int zsymm_t(char *side, char *uplo, int *m, int *n, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
ctypedef int zsyr_t(char *uplo, int *n, z *alpha, z *x, int *incx, z *a, int *lda) nogil
ctypedef int zsyr2k_t(char *uplo, char *trans, int *n, int *k, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
ctypedef int zsyrk_t(char *uplo, char *trans, int *n, int *k, z *alpha, z *a, int *lda, z *beta, z *c, int *ldc) nogil
ctypedef int ztrmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, z *alpha, z *a, int *lda, z *b, int *ldb) nogil
ctypedef int ztrmv_t(char *uplo, char *trans, char *diag, int *n, z *a, int *lda, z *x, int *incx) nogil

cdef:
    caxpy_t *caxpy_f
    ccopy_t *ccopy_f
    cdotc_t *cdotc_f
    cdotu_t *cdotu_f
    cgemm_t *cgemm_f
    cgemv_t *cgemv_f
    cgerc_t *cgerc_f
    cgeru_t *cgeru_f
    chemm_t *chemm_f
    chemv_t *chemv_f
    cher_t *cher_f
    cher2_t *cher2_f
    cher2k_t *cher2k_f
    cherk_t *cherk_f
    crotg_t *crotg_f
    cscal_t *cscal_f
    csrot_t *csrot_f
    csscal_t *csscal_f
    cswap_t *cswap_f
    csymm_t *csymm_f
    csyr_t *csyr_f
    csyr2k_t *csyr2k_f
    csyrk_t *csyrk_f
    ctrmm_t *ctrmm_f
    ctrmv_t *ctrmv_f
    dasum_t *dasum_f
    daxpy_t *daxpy_f
    dcopy_t *dcopy_f
    ddot_t *ddot_f
    dgemm_t *dgemm_f
    dgemv_t *dgemv_f
    dger_t *dger_f
    dnrm2_t *dnrm2_f
    drot_t *drot_f
    drotg_t *drotg_f
    drotm_t *drotm_f
    drotmg_t *drotmg_f
    dscal_t *dscal_f
    dswap_t *dswap_f
    dsymm_t *dsymm_f
    dsymv_t *dsymv_f
    dsyr_t *dsyr_f
    dsyr2_t *dsyr2_f
    dsyr2k_t *dsyr2k_f
    dsyrk_t *dsyrk_f
    dtrmm_t *dtrmm_f
    dtrmv_t *dtrmv_f
    dzasum_t *dzasum_f
    dznrm2_t *dznrm2_f
    icamax_t *icamax_f
    idamax_t *idamax_f
    isamax_t *isamax_f
    izamax_t *izamax_f
    sasum_t *sasum_f
    saxpy_t *saxpy_f
    scasum_t *scasum_f
    scnrm2_t *scnrm2_f
    scopy_t *scopy_f
    sdot_t *sdot_f
    sgemm_t *sgemm_f
    sgemv_t *sgemv_f
    sger_t *sger_f
    snrm2_t *snrm2_f
    srot_t *srot_f
    srotg_t *srotg_f
    srotm_t *srotm_f
    srotmg_t *srotmg_f
    sscal_t *sscal_f
    sswap_t *sswap_f
    ssymm_t *ssymm_f
    ssymv_t *ssymv_f
    ssyr_t *ssyr_f
    ssyr2_t *ssyr2_f
    ssyr2k_t *ssyr2k_f
    ssyrk_t *ssyrk_f
    strmm_t *strmm_f
    strmv_t *strmv_f
    zaxpy_t *zaxpy_f
    zcopy_t *zcopy_f
    zdotc_t *zdotc_f
    zdotu_t *zdotu_f
    zdrot_t *zdrot_f
    zdscal_t *zdscal_f
    zgemm_t *zgemm_f
    zgemv_t *zgemv_f
    zgerc_t *zgerc_f
    zgeru_t *zgeru_f
    zhemm_t *zhemm_f
    zhemv_t *zhemv_f
    zher_t *zher_f
    zher2_t *zher2_f
    zher2k_t *zher2k_f
    zherk_t *zherk_f
    zrotg_t *zrotg_f
    zscal_t *zscal_f
    zswap_t *zswap_f
    zsymm_t *zsymm_f
    zsyr_t *zsyr_f
    zsyr2k_t *zsyr2k_f
    zsyrk_t *zsyrk_f
    ztrmm_t *ztrmm_f
    ztrmv_t *ztrmv_f

cpdef float complex cdotc(float complex[:] cx, float complex[:] cy)
cpdef float complex cdotu(float complex[:] cx, float complex[:] cy)
cpdef double dasum(double[:] dx)
cpdef double ddot(double[:] dx, double[:] dy)
cpdef int dgemm(double alpha, double[:,:] a, double[:,:] b, double beta, double[:,:] c) except -1
cpdef double dnrm2(double[:] x)
cpdef double dzasum(double complex[:] zx)
cpdef double dznrm2(double complex[:] x)
cpdef int icamax(float complex[:] cx)
cpdef int idamax(double[:] dx)
cpdef int isamax(float[:] sx)
cpdef int izamax(double complex[:] zx)
cpdef float sasum(float[:] sx)
cpdef float scasum(float complex[:] cx)
cpdef float scnrm2(float complex[:] x)
cpdef float sdot(float[:] sx, float[:] sy)
cpdef float snrm2(float[:] x)
cpdef double complex zdotc(double complex[:] zx, double complex[:] zy)
cpdef double complex zdotu(double complex[:] zx, double complex[:] zy)
