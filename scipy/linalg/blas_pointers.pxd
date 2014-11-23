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
    caxpy_t *caxpy
    ccopy_t *ccopy
    cdotc_t *cdotc
    cdotu_t *cdotu
    cgemm_t *cgemm
    cgemv_t *cgemv
    cgerc_t *cgerc
    cgeru_t *cgeru
    chemm_t *chemm
    chemv_t *chemv
    cher_t *cher
    cher2_t *cher2
    cher2k_t *cher2k
    cherk_t *cherk
    crotg_t *crotg
    cscal_t *cscal
    csrot_t *csrot
    csscal_t *csscal
    cswap_t *cswap
    csymm_t *csymm
    csyr_t *csyr
    csyr2k_t *csyr2k
    csyrk_t *csyrk
    ctrmm_t *ctrmm
    ctrmv_t *ctrmv
    dasum_t *dasum
    daxpy_t *daxpy
    dcopy_t *dcopy
    ddot_t *ddot
    dgemm_t *dgemm
    dgemv_t *dgemv
    dger_t *dger
    dnrm2_t *dnrm2
    drot_t *drot
    drotg_t *drotg
    drotm_t *drotm
    drotmg_t *drotmg
    dscal_t *dscal
    dswap_t *dswap
    dsymm_t *dsymm
    dsymv_t *dsymv
    dsyr_t *dsyr
    dsyr2_t *dsyr2
    dsyr2k_t *dsyr2k
    dsyrk_t *dsyrk
    dtrmm_t *dtrmm
    dtrmv_t *dtrmv
    dzasum_t *dzasum
    dznrm2_t *dznrm2
    icamax_t *icamax
    idamax_t *idamax
    isamax_t *isamax
    izamax_t *izamax
    sasum_t *sasum
    saxpy_t *saxpy
    scasum_t *scasum
    scnrm2_t *scnrm2
    scopy_t *scopy
    sdot_t *sdot
    sgemm_t *sgemm
    sgemv_t *sgemv
    sger_t *sger
    snrm2_t *snrm2
    srot_t *srot
    srotg_t *srotg
    srotm_t *srotm
    srotmg_t *srotmg
    sscal_t *sscal
    sswap_t *sswap
    ssymm_t *ssymm
    ssymv_t *ssymv
    ssyr_t *ssyr
    ssyr2_t *ssyr2
    ssyr2k_t *ssyr2k
    ssyrk_t *ssyrk
    strmm_t *strmm
    strmv_t *strmv
    zaxpy_t *zaxpy
    zcopy_t *zcopy
    zdotc_t *zdotc
    zdotu_t *zdotu
    zdrot_t *zdrot
    zdscal_t *zdscal
    zgemm_t *zgemm
    zgemv_t *zgemv
    zgerc_t *zgerc
    zgeru_t *zgeru
    zhemm_t *zhemm
    zhemv_t *zhemv
    zher_t *zher
    zher2_t *zher2
    zher2k_t *zher2k
    zherk_t *zherk
    zrotg_t *zrotg
    zscal_t *zscal
    zswap_t *zswap
    zsymm_t *zsymm
    zsyr_t *zsyr
    zsyr2k_t *zsyr2k
    zsyrk_t *zsyrk
    ztrmm_t *ztrmm
    ztrmv_t *ztrmv
