ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

ctypedef void caxpy_t(int *n, c *ca, c *cx, int *incx, c *cy, int *incy) nogil
cdef caxpy_t *caxpy_f

ctypedef void ccopy_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
cdef ccopy_t *ccopy_f

ctypedef c cdotc_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
cdef cdotc_t *cdotc_f

ctypedef c cdotu_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
cdef cdotu_t *cdotu_f

ctypedef void cgbmv_t(char *trans, int *m, int *n, int *kl, int *ku, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef cgbmv_t *cgbmv_f

ctypedef void cgemm_t(char *transa, char *transb, int *m, int *n, int *k, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
cdef cgemm_t *cgemm_f

ctypedef void cgemv_t(char *trans, int *m, int *n, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef cgemv_t *cgemv_f

ctypedef void cgerc_t(int *m, int *n, c *alpha, c *x, int *incx, c *y, int *incy, c *a, int *lda) nogil
cdef cgerc_t *cgerc_f

ctypedef void cgeru_t(int *m, int *n, c *alpha, c *x, int *incx, c *y, int *incy, c *a, int *lda) nogil
cdef cgeru_t *cgeru_f

ctypedef void chbmv_t(char *uplo, int *n, int *k, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef chbmv_t *chbmv_f

ctypedef void chemm_t(char *side, char *uplo, int *m, int *n, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
cdef chemm_t *chemm_f

ctypedef void chemv_t(char *uplo, int *n, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef chemv_t *chemv_f

ctypedef void cher_t(char *uplo, int *n, s *alpha, c *x, int *incx, c *a, int *lda) nogil
cdef cher_t *cher_f

ctypedef void cher2_t(char *uplo, int *n, c *alpha, c *x, int *incx, c *y, int *incy, c *a, int *lda) nogil
cdef cher2_t *cher2_f

ctypedef void cher2k_t(char *uplo, char *trans, int *n, int *k, c *alpha, c *a, int *lda, c *b, int *ldb, s *beta, c *c, int *ldc) nogil
cdef cher2k_t *cher2k_f

ctypedef void cherk_t(char *uplo, char *trans, int *n, int *k, s *alpha, c *a, int *lda, s *beta, c *c, int *ldc) nogil
cdef cherk_t *cherk_f

ctypedef void chpmv_t(char *uplo, int *n, c *alpha, c *ap, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef chpmv_t *chpmv_f

ctypedef void chpr_t(char *uplo, int *n, s *alpha, c *x, int *incx, c *ap) nogil
cdef chpr_t *chpr_f

ctypedef void chpr2_t(char *uplo, int *n, c *alpha, c *x, int *incx, c *y, int *incy, c *ap) nogil
cdef chpr2_t *chpr2_f

ctypedef void crotg_t(c *ca, c *cb, s *c, c *s) nogil
cdef crotg_t *crotg_f

ctypedef void cscal_t(int *n, c *ca, c *cx, int *incx) nogil
cdef cscal_t *cscal_f

ctypedef void csrot_t(int *n, c *cx, int *incx, c *cy, int *incy, s *c, s *s) nogil
cdef csrot_t *csrot_f

ctypedef void csscal_t(int *n, s *sa, c *cx, int *incx) nogil
cdef csscal_t *csscal_f

ctypedef void cswap_t(int *n, c *cx, int *incx, c *cy, int *incy) nogil
cdef cswap_t *cswap_f

ctypedef void csymm_t(char *side, char *uplo, int *m, int *n, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
cdef csymm_t *csymm_f

ctypedef void csyr2k_t(char *uplo, char *trans, int *n, int *k, c *alpha, c *a, int *lda, c *b, int *ldb, c *beta, c *c, int *ldc) nogil
cdef csyr2k_t *csyr2k_f

ctypedef void csyrk_t(char *uplo, char *trans, int *n, int *k, c *alpha, c *a, int *lda, c *beta, c *c, int *ldc) nogil
cdef csyrk_t *csyrk_f

ctypedef void ctbmv_t(char *uplo, char *trans, char *diag, int *n, int *k, c *a, int *lda, c *x, int *incx) nogil
cdef ctbmv_t *ctbmv_f

ctypedef void ctbsv_t(char *uplo, char *trans, char *diag, int *n, int *k, c *a, int *lda, c *x, int *incx) nogil
cdef ctbsv_t *ctbsv_f

ctypedef void ctpmv_t(char *uplo, char *trans, char *diag, int *n, c *ap, c *x, int *incx) nogil
cdef ctpmv_t *ctpmv_f

ctypedef void ctpsv_t(char *uplo, char *trans, char *diag, int *n, c *ap, c *x, int *incx) nogil
cdef ctpsv_t *ctpsv_f

ctypedef void ctrmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, c *alpha, c *a, int *lda, c *b, int *ldb) nogil
cdef ctrmm_t *ctrmm_f

ctypedef void ctrmv_t(char *uplo, char *trans, char *diag, int *n, c *a, int *lda, c *x, int *incx) nogil
cdef ctrmv_t *ctrmv_f

ctypedef void ctrsm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, c *alpha, c *a, int *lda, c *b, int *ldb) nogil
cdef ctrsm_t *ctrsm_f

ctypedef void ctrsv_t(char *uplo, char *trans, char *diag, int *n, c *a, int *lda, c *x, int *incx) nogil
cdef ctrsv_t *ctrsv_f

ctypedef d dasum_t(int *n, d *dx, int *incx) nogil
cdef dasum_t *dasum_f

ctypedef void daxpy_t(int *n, d *da, d *dx, int *incx, d *dy, int *incy) nogil
cdef daxpy_t *daxpy_f

ctypedef d dcabs1_t(z *z) nogil
cdef dcabs1_t *dcabs1_f

ctypedef void dcopy_t(int *n, d *dx, int *incx, d *dy, int *incy) nogil
cdef dcopy_t *dcopy_f

ctypedef d ddot_t(int *n, d *dx, int *incx, d *dy, int *incy) nogil
cdef ddot_t *ddot_f

ctypedef void dgbmv_t(char *trans, int *m, int *n, int *kl, int *ku, d *alpha, d *a, int *lda, d *x, int *incx, d *beta, d *y, int *incy) nogil
cdef dgbmv_t *dgbmv_f

ctypedef void dgemm_t(char *transa, char *transb, int *m, int *n, int *k, d *alpha, d *a, int *lda, d *b, int *ldb, d *beta, d *c, int *ldc) nogil
cdef dgemm_t *dgemm_f

ctypedef void dgemv_t(char *trans, int *m, int *n, d *alpha, d *a, int *lda, d *x, int *incx, d *beta, d *y, int *incy) nogil
cdef dgemv_t *dgemv_f

ctypedef void dger_t(int *m, int *n, d *alpha, d *x, int *incx, d *y, int *incy, d *a, int *lda) nogil
cdef dger_t *dger_f

ctypedef d dnrm2_t(int *n, d *x, int *incx) nogil
cdef dnrm2_t *dnrm2_f

ctypedef void drot_t(int *n, d *dx, int *incx, d *dy, int *incy, d *c, d *s) nogil
cdef drot_t *drot_f

ctypedef void drotg_t(d *da, d *db, d *c, d *s) nogil
cdef drotg_t *drotg_f

ctypedef void drotm_t(int *n, d *dx, int *incx, d *dy, int *incy, d *dparam) nogil
cdef drotm_t *drotm_f

ctypedef void drotmg_t(d *dd1, d *dd2, d *dx1, d *dy1, d *dparam) nogil
cdef drotmg_t *drotmg_f

ctypedef void dsbmv_t(char *uplo, int *n, int *k, d *alpha, d *a, int *lda, d *x, int *incx, d *beta, d *y, int *incy) nogil
cdef dsbmv_t *dsbmv_f

ctypedef void dscal_t(int *n, d *da, d *dx, int *incx) nogil
cdef dscal_t *dscal_f

ctypedef d dsdot_t(int *n, s *sx, int *incx, s *sy, int *incy) nogil
cdef dsdot_t *dsdot_f

ctypedef void dspmv_t(char *uplo, int *n, d *alpha, d *ap, d *x, int *incx, d *beta, d *y, int *incy) nogil
cdef dspmv_t *dspmv_f

ctypedef void dspr_t(char *uplo, int *n, d *alpha, d *x, int *incx, d *ap) nogil
cdef dspr_t *dspr_f

ctypedef void dspr2_t(char *uplo, int *n, d *alpha, d *x, int *incx, d *y, int *incy, d *ap) nogil
cdef dspr2_t *dspr2_f

ctypedef void dswap_t(int *n, d *dx, int *incx, d *dy, int *incy) nogil
cdef dswap_t *dswap_f

ctypedef void dsymm_t(char *side, char *uplo, int *m, int *n, d *alpha, d *a, int *lda, d *b, int *ldb, d *beta, d *c, int *ldc) nogil
cdef dsymm_t *dsymm_f

ctypedef void dsymv_t(char *uplo, int *n, d *alpha, d *a, int *lda, d *x, int *incx, d *beta, d *y, int *incy) nogil
cdef dsymv_t *dsymv_f

ctypedef void dsyr_t(char *uplo, int *n, d *alpha, d *x, int *incx, d *a, int *lda) nogil
cdef dsyr_t *dsyr_f

ctypedef void dsyr2_t(char *uplo, int *n, d *alpha, d *x, int *incx, d *y, int *incy, d *a, int *lda) nogil
cdef dsyr2_t *dsyr2_f

ctypedef void dsyr2k_t(char *uplo, char *trans, int *n, int *k, d *alpha, d *a, int *lda, d *b, int *ldb, d *beta, d *c, int *ldc) nogil
cdef dsyr2k_t *dsyr2k_f

ctypedef void dsyrk_t(char *uplo, char *trans, int *n, int *k, d *alpha, d *a, int *lda, d *beta, d *c, int *ldc) nogil
cdef dsyrk_t *dsyrk_f

ctypedef void dtbmv_t(char *uplo, char *trans, char *diag, int *n, int *k, d *a, int *lda, d *x, int *incx) nogil
cdef dtbmv_t *dtbmv_f

ctypedef void dtbsv_t(char *uplo, char *trans, char *diag, int *n, int *k, d *a, int *lda, d *x, int *incx) nogil
cdef dtbsv_t *dtbsv_f

ctypedef void dtpmv_t(char *uplo, char *trans, char *diag, int *n, d *ap, d *x, int *incx) nogil
cdef dtpmv_t *dtpmv_f

ctypedef void dtpsv_t(char *uplo, char *trans, char *diag, int *n, d *ap, d *x, int *incx) nogil
cdef dtpsv_t *dtpsv_f

ctypedef void dtrmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, d *alpha, d *a, int *lda, d *b, int *ldb) nogil
cdef dtrmm_t *dtrmm_f

ctypedef void dtrmv_t(char *uplo, char *trans, char *diag, int *n, d *a, int *lda, d *x, int *incx) nogil
cdef dtrmv_t *dtrmv_f

ctypedef void dtrsm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, d *alpha, d *a, int *lda, d *b, int *ldb) nogil
cdef dtrsm_t *dtrsm_f

ctypedef void dtrsv_t(char *uplo, char *trans, char *diag, int *n, d *a, int *lda, d *x, int *incx) nogil
cdef dtrsv_t *dtrsv_f

ctypedef d dzasum_t(int *n, z *zx, int *incx) nogil
cdef dzasum_t *dzasum_f

ctypedef d dznrm2_t(int *n, z *x, int *incx) nogil
cdef dznrm2_t *dznrm2_f

ctypedef int icamax_t(int *n, c *cx, int *incx) nogil
cdef icamax_t *icamax_f

ctypedef int idamax_t(int *n, d *dx, int *incx) nogil
cdef idamax_t *idamax_f

ctypedef int isamax_t(int *n, s *sx, int *incx) nogil
cdef isamax_t *isamax_f

ctypedef int izamax_t(int *n, z *zx, int *incx) nogil
cdef izamax_t *izamax_f

ctypedef bint lsame_t(char *ca, char *cb) nogil
cdef lsame_t *lsame_f

ctypedef s sasum_t(int *n, s *sx, int *incx) nogil
cdef sasum_t *sasum_f

ctypedef void saxpy_t(int *n, s *sa, s *sx, int *incx, s *sy, int *incy) nogil
cdef saxpy_t *saxpy_f

ctypedef s scasum_t(int *n, c *cx, int *incx) nogil
cdef scasum_t *scasum_f

ctypedef s scnrm2_t(int *n, c *x, int *incx) nogil
cdef scnrm2_t *scnrm2_f

ctypedef void scopy_t(int *n, s *sx, int *incx, s *sy, int *incy) nogil
cdef scopy_t *scopy_f

ctypedef s sdot_t(int *n, s *sx, int *incx, s *sy, int *incy) nogil
cdef sdot_t *sdot_f

ctypedef s sdsdot_t(int *n, s *sb, s *sx, int *incx, s *sy, int *incy) nogil
cdef sdsdot_t *sdsdot_f

ctypedef void sgbmv_t(char *trans, int *m, int *n, int *kl, int *ku, s *alpha, s *a, int *lda, s *x, int *incx, s *beta, s *y, int *incy) nogil
cdef sgbmv_t *sgbmv_f

ctypedef void sgemm_t(char *transa, char *transb, int *m, int *n, int *k, s *alpha, s *a, int *lda, s *b, int *ldb, s *beta, s *c, int *ldc) nogil
cdef sgemm_t *sgemm_f

ctypedef void sgemv_t(char *trans, int *m, int *n, s *alpha, s *a, int *lda, s *x, int *incx, s *beta, s *y, int *incy) nogil
cdef sgemv_t *sgemv_f

ctypedef void sger_t(int *m, int *n, s *alpha, s *x, int *incx, s *y, int *incy, s *a, int *lda) nogil
cdef sger_t *sger_f

ctypedef s snrm2_t(int *n, s *x, int *incx) nogil
cdef snrm2_t *snrm2_f

ctypedef void srot_t(int *n, s *sx, int *incx, s *sy, int *incy, s *c, s *s) nogil
cdef srot_t *srot_f

ctypedef void srotg_t(s *sa, s *sb, s *c, s *s) nogil
cdef srotg_t *srotg_f

ctypedef void srotm_t(int *n, s *sx, int *incx, s *sy, int *incy, s *sparam) nogil
cdef srotm_t *srotm_f

ctypedef void srotmg_t(s *sd1, s *sd2, s *sx1, s *sy1, s *sparam) nogil
cdef srotmg_t *srotmg_f

ctypedef void ssbmv_t(char *uplo, int *n, int *k, s *alpha, s *a, int *lda, s *x, int *incx, s *beta, s *y, int *incy) nogil
cdef ssbmv_t *ssbmv_f

ctypedef void sscal_t(int *n, s *sa, s *sx, int *incx) nogil
cdef sscal_t *sscal_f

ctypedef void sspmv_t(char *uplo, int *n, s *alpha, s *ap, s *x, int *incx, s *beta, s *y, int *incy) nogil
cdef sspmv_t *sspmv_f

ctypedef void sspr_t(char *uplo, int *n, s *alpha, s *x, int *incx, s *ap) nogil
cdef sspr_t *sspr_f

ctypedef void sspr2_t(char *uplo, int *n, s *alpha, s *x, int *incx, s *y, int *incy, s *ap) nogil
cdef sspr2_t *sspr2_f

ctypedef void sswap_t(int *n, s *sx, int *incx, s *sy, int *incy) nogil
cdef sswap_t *sswap_f

ctypedef void ssymm_t(char *side, char *uplo, int *m, int *n, s *alpha, s *a, int *lda, s *b, int *ldb, s *beta, s *c, int *ldc) nogil
cdef ssymm_t *ssymm_f

ctypedef void ssymv_t(char *uplo, int *n, s *alpha, s *a, int *lda, s *x, int *incx, s *beta, s *y, int *incy) nogil
cdef ssymv_t *ssymv_f

ctypedef void ssyr_t(char *uplo, int *n, s *alpha, s *x, int *incx, s *a, int *lda) nogil
cdef ssyr_t *ssyr_f

ctypedef void ssyr2_t(char *uplo, int *n, s *alpha, s *x, int *incx, s *y, int *incy, s *a, int *lda) nogil
cdef ssyr2_t *ssyr2_f

ctypedef void ssyr2k_t(char *uplo, char *trans, int *n, int *k, s *alpha, s *a, int *lda, s *b, int *ldb, s *beta, s *c, int *ldc) nogil
cdef ssyr2k_t *ssyr2k_f

ctypedef void ssyrk_t(char *uplo, char *trans, int *n, int *k, s *alpha, s *a, int *lda, s *beta, s *c, int *ldc) nogil
cdef ssyrk_t *ssyrk_f

ctypedef void stbmv_t(char *uplo, char *trans, char *diag, int *n, int *k, s *a, int *lda, s *x, int *incx) nogil
cdef stbmv_t *stbmv_f

ctypedef void stbsv_t(char *uplo, char *trans, char *diag, int *n, int *k, s *a, int *lda, s *x, int *incx) nogil
cdef stbsv_t *stbsv_f

ctypedef void stpmv_t(char *uplo, char *trans, char *diag, int *n, s *ap, s *x, int *incx) nogil
cdef stpmv_t *stpmv_f

ctypedef void stpsv_t(char *uplo, char *trans, char *diag, int *n, s *ap, s *x, int *incx) nogil
cdef stpsv_t *stpsv_f

ctypedef void strmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, s *alpha, s *a, int *lda, s *b, int *ldb) nogil
cdef strmm_t *strmm_f

ctypedef void strmv_t(char *uplo, char *trans, char *diag, int *n, s *a, int *lda, s *x, int *incx) nogil
cdef strmv_t *strmv_f

ctypedef void strsm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, s *alpha, s *a, int *lda, s *b, int *ldb) nogil
cdef strsm_t *strsm_f

ctypedef void strsv_t(char *uplo, char *trans, char *diag, int *n, s *a, int *lda, s *x, int *incx) nogil
cdef strsv_t *strsv_f

ctypedef void xerbla_t(char *srname, int *info) nogil
cdef xerbla_t *xerbla_f

ctypedef void zaxpy_t(int *n, z *za, z *zx, int *incx, z *zy, int *incy) nogil
cdef zaxpy_t *zaxpy_f

ctypedef void zcopy_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
cdef zcopy_t *zcopy_f

ctypedef z zdotc_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
cdef zdotc_t *zdotc_f

ctypedef z zdotu_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
cdef zdotu_t *zdotu_f

ctypedef void zdrot_t(int *n, z *cx, int *incx, z *cy, int *incy, d *c, d *s) nogil
cdef zdrot_t *zdrot_f

ctypedef void zdscal_t(int *n, d *da, z *zx, int *incx) nogil
cdef zdscal_t *zdscal_f

ctypedef void zgbmv_t(char *trans, int *m, int *n, int *kl, int *ku, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zgbmv_t *zgbmv_f

ctypedef void zgemm_t(char *transa, char *transb, int *m, int *n, int *k, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
cdef zgemm_t *zgemm_f

ctypedef void zgemv_t(char *trans, int *m, int *n, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zgemv_t *zgemv_f

ctypedef void zgerc_t(int *m, int *n, z *alpha, z *x, int *incx, z *y, int *incy, z *a, int *lda) nogil
cdef zgerc_t *zgerc_f

ctypedef void zgeru_t(int *m, int *n, z *alpha, z *x, int *incx, z *y, int *incy, z *a, int *lda) nogil
cdef zgeru_t *zgeru_f

ctypedef void zhbmv_t(char *uplo, int *n, int *k, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zhbmv_t *zhbmv_f

ctypedef void zhemm_t(char *side, char *uplo, int *m, int *n, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
cdef zhemm_t *zhemm_f

ctypedef void zhemv_t(char *uplo, int *n, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zhemv_t *zhemv_f

ctypedef void zher_t(char *uplo, int *n, d *alpha, z *x, int *incx, z *a, int *lda) nogil
cdef zher_t *zher_f

ctypedef void zher2_t(char *uplo, int *n, z *alpha, z *x, int *incx, z *y, int *incy, z *a, int *lda) nogil
cdef zher2_t *zher2_f

ctypedef void zher2k_t(char *uplo, char *trans, int *n, int *k, z *alpha, z *a, int *lda, z *b, int *ldb, d *beta, z *c, int *ldc) nogil
cdef zher2k_t *zher2k_f

ctypedef void zherk_t(char *uplo, char *trans, int *n, int *k, d *alpha, z *a, int *lda, d *beta, z *c, int *ldc) nogil
cdef zherk_t *zherk_f

ctypedef void zhpmv_t(char *uplo, int *n, z *alpha, z *ap, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zhpmv_t *zhpmv_f

ctypedef void zhpr_t(char *uplo, int *n, d *alpha, z *x, int *incx, z *ap) nogil
cdef zhpr_t *zhpr_f

ctypedef void zhpr2_t(char *uplo, int *n, z *alpha, z *x, int *incx, z *y, int *incy, z *ap) nogil
cdef zhpr2_t *zhpr2_f

ctypedef void zrotg_t(z *ca, z *cb, d *c, z *s) nogil
cdef zrotg_t *zrotg_f

ctypedef void zscal_t(int *n, z *za, z *zx, int *incx) nogil
cdef zscal_t *zscal_f

ctypedef void zswap_t(int *n, z *zx, int *incx, z *zy, int *incy) nogil
cdef zswap_t *zswap_f

ctypedef void zsymm_t(char *side, char *uplo, int *m, int *n, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
cdef zsymm_t *zsymm_f

ctypedef void zsyr2k_t(char *uplo, char *trans, int *n, int *k, z *alpha, z *a, int *lda, z *b, int *ldb, z *beta, z *c, int *ldc) nogil
cdef zsyr2k_t *zsyr2k_f

ctypedef void zsyrk_t(char *uplo, char *trans, int *n, int *k, z *alpha, z *a, int *lda, z *beta, z *c, int *ldc) nogil
cdef zsyrk_t *zsyrk_f

ctypedef void ztbmv_t(char *uplo, char *trans, char *diag, int *n, int *k, z *a, int *lda, z *x, int *incx) nogil
cdef ztbmv_t *ztbmv_f

ctypedef void ztbsv_t(char *uplo, char *trans, char *diag, int *n, int *k, z *a, int *lda, z *x, int *incx) nogil
cdef ztbsv_t *ztbsv_f

ctypedef void ztpmv_t(char *uplo, char *trans, char *diag, int *n, z *ap, z *x, int *incx) nogil
cdef ztpmv_t *ztpmv_f

ctypedef void ztpsv_t(char *uplo, char *trans, char *diag, int *n, z *ap, z *x, int *incx) nogil
cdef ztpsv_t *ztpsv_f

ctypedef void ztrmm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, z *alpha, z *a, int *lda, z *b, int *ldb) nogil
cdef ztrmm_t *ztrmm_f

ctypedef void ztrmv_t(char *uplo, char *trans, char *diag, int *n, z *a, int *lda, z *x, int *incx) nogil
cdef ztrmv_t *ztrmv_f

ctypedef void ztrsm_t(char *side, char *uplo, char *transa, char *diag, int *m, int *n, z *alpha, z *a, int *lda, z *b, int *ldb) nogil
cdef ztrsm_t *ztrsm_f

ctypedef void ztrsv_t(char *uplo, char *trans, char *diag, int *n, z *a, int *lda, z *x, int *incx) nogil
cdef ztrsv_t *ztrsv_f
