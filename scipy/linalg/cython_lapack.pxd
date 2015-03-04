ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

# Function pointer type declarations for
# gees and gges families of functions.
ctypedef bint cselect1(c*)
ctypedef bint cselect2(c*, c*)
ctypedef bint dselect2(d*, d*)
ctypedef bint dselect3(d*, d*, d*)
ctypedef bint sselect2(s*, s*)
ctypedef bint sselect3(s*, s*, s*)
ctypedef bint zselect1(z*)
ctypedef bint zselect2(z*, z*)

ctypedef void cbdsqr_t(char *uplo, int *n, int *ncvt, int *nru, int *ncc, s *d, s *e, c *vt, int *ldvt, c *u, int *ldu, c *c, int *ldc, s *rwork, int *info) nogil
cdef cbdsqr_t *cbdsqr_f

ctypedef void cgbbrd_t(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, c *ab, int *ldab, s *d, s *e, c *q, int *ldq, c *pt, int *ldpt, c *c, int *ldc, c *work, s *rwork, int *info) nogil
cdef cgbbrd_t *cgbbrd_f

ctypedef void cgbcon_t(char *norm, int *n, int *kl, int *ku, c *ab, int *ldab, int *ipiv, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cgbcon_t *cgbcon_f

ctypedef void cgbequ_t(int *m, int *n, int *kl, int *ku, c *ab, int *ldab, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef cgbequ_t *cgbequ_f

ctypedef void cgbrfs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgbrfs_t *cgbrfs_f

ctypedef void cgbsv_t(int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgbsv_t *cgbsv_f

ctypedef void cgbsvx_t(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, int *ipiv, char *equed, s *r, s *c, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgbsvx_t *cgbsvx_f

ctypedef void cgbtf2_t(int *m, int *n, int *kl, int *ku, c *ab, int *ldab, int *ipiv, int *info) nogil
cdef cgbtf2_t *cgbtf2_f

ctypedef void cgbtrf_t(int *m, int *n, int *kl, int *ku, c *ab, int *ldab, int *ipiv, int *info) nogil
cdef cgbtrf_t *cgbtrf_f

ctypedef void cgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgbtrs_t *cgbtrs_f

ctypedef void cgebak_t(char *job, char *side, int *n, int *ilo, int *ihi, s *scale, int *m, c *v, int *ldv, int *info) nogil
cdef cgebak_t *cgebak_f

ctypedef void cgebal_t(char *job, int *n, c *a, int *lda, int *ilo, int *ihi, s *scale, int *info) nogil
cdef cgebal_t *cgebal_f

ctypedef void cgebd2_t(int *m, int *n, c *a, int *lda, s *d, s *e, c *tauq, c *taup, c *work, int *info) nogil
cdef cgebd2_t *cgebd2_f

ctypedef void cgebrd_t(int *m, int *n, c *a, int *lda, s *d, s *e, c *tauq, c *taup, c *work, int *lwork, int *info) nogil
cdef cgebrd_t *cgebrd_f

ctypedef void cgecon_t(char *norm, int *n, c *a, int *lda, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cgecon_t *cgecon_f

ctypedef void cgeequ_t(int *m, int *n, c *a, int *lda, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef cgeequ_t *cgeequ_f

ctypedef void cgees_t(char *jobvs, char *sort, cselect1 *select, int *n, c *a, int *lda, int *sdim, c *w, c *vs, int *ldvs, c *work, int *lwork, s *rwork, bint *bwork, int *info) nogil
cdef cgees_t *cgees_f

ctypedef void cgeesx_t(char *jobvs, char *sort, cselect1 *select, char *sense, int *n, c *a, int *lda, int *sdim, c *w, c *vs, int *ldvs, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, bint *bwork, int *info) nogil
cdef cgeesx_t *cgeesx_f

ctypedef void cgeev_t(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *w, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgeev_t *cgeev_f

ctypedef void cgeevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, c *a, int *lda, c *w, c *vl, int *ldvl, c *vr, int *ldvr, int *ilo, int *ihi, s *scale, s *abnrm, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgeevx_t *cgeevx_f

ctypedef void cgegs_t(char *jobvsl, char *jobvsr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vsl, int *ldvsl, c *vsr, int *ldvsr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgegs_t *cgegs_f

ctypedef void cgegv_t(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgegv_t *cgegv_f

ctypedef void cgehd2_t(int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgehd2_t *cgehd2_f

ctypedef void cgehrd_t(int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgehrd_t *cgehrd_f

ctypedef void cgelq2_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgelq2_t *cgelq2_f

ctypedef void cgelqf_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgelqf_t *cgelqf_f

ctypedef void cgels_t(char *trans, int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, c *work, int *lwork, int *info) nogil
cdef cgels_t *cgels_f

ctypedef void cgelsd_t(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, s *s, s *rcond, int *rank, c *work, int *lwork, s *rwork, int *iwork, int *info) nogil
cdef cgelsd_t *cgelsd_f

ctypedef void cgelss_t(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, s *s, s *rcond, int *rank, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgelss_t *cgelss_f

ctypedef void cgelsx_t(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *jpvt, s *rcond, int *rank, c *work, s *rwork, int *info) nogil
cdef cgelsx_t *cgelsx_f

ctypedef void cgelsy_t(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *jpvt, s *rcond, int *rank, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgelsy_t *cgelsy_f

ctypedef void cgeql2_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgeql2_t *cgeql2_f

ctypedef void cgeqlf_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgeqlf_t *cgeqlf_f

ctypedef void cgeqp3_t(int *m, int *n, c *a, int *lda, int *jpvt, c *tau, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgeqp3_t *cgeqp3_f

ctypedef void cgeqpf_t(int *m, int *n, c *a, int *lda, int *jpvt, c *tau, c *work, s *rwork, int *info) nogil
cdef cgeqpf_t *cgeqpf_f

ctypedef void cgeqr2_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgeqr2_t *cgeqr2_f

ctypedef void cgeqrf_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgeqrf_t *cgeqrf_f

ctypedef void cgerfs_t(char *trans, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgerfs_t *cgerfs_f

ctypedef void cgerq2_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgerq2_t *cgerq2_f

ctypedef void cgerqf_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgerqf_t *cgerqf_f

ctypedef void cgesc2_t(int *n, c *a, int *lda, c *rhs, int *ipiv, int *jpiv, s *scale) nogil
cdef cgesc2_t *cgesc2_f

ctypedef void cgesdd_t(char *jobz, int *m, int *n, c *a, int *lda, s *s, c *u, int *ldu, c *vt, int *ldvt, c *work, int *lwork, s *rwork, int *iwork, int *info) nogil
cdef cgesdd_t *cgesdd_f

ctypedef void cgesv_t(int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgesv_t *cgesv_f

ctypedef void cgesvd_t(char *jobu, char *jobvt, int *m, int *n, c *a, int *lda, s *s, c *u, int *ldu, c *vt, int *ldvt, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgesvd_t *cgesvd_f

ctypedef void cgesvx_t(char *fact, char *trans, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, char *equed, s *r, s *c, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgesvx_t *cgesvx_f

ctypedef void cgetc2_t(int *n, c *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef cgetc2_t *cgetc2_f

ctypedef void cgetf2_t(int *m, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef cgetf2_t *cgetf2_f

ctypedef void cgetrf_t(int *m, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef cgetrf_t *cgetrf_f

ctypedef void cgetri_t(int *n, c *a, int *lda, int *ipiv, c *work, int *lwork, int *info) nogil
cdef cgetri_t *cgetri_f

ctypedef void cgetrs_t(char *trans, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgetrs_t *cgetrs_f

ctypedef void cggbak_t(char *job, char *side, int *n, int *ilo, int *ihi, s *lscale, s *rscale, int *m, c *v, int *ldv, int *info) nogil
cdef cggbak_t *cggbak_f

ctypedef void cggbal_t(char *job, int *n, c *a, int *lda, c *b, int *ldb, int *ilo, int *ihi, s *lscale, s *rscale, s *work, int *info) nogil
cdef cggbal_t *cggbal_f

ctypedef void cgges_t(char *jobvsl, char *jobvsr, char *sort, cselect2 *selctg, int *n, c *a, int *lda, c *b, int *ldb, int *sdim, c *alpha, c *beta, c *vsl, int *ldvsl, c *vsr, int *ldvsr, c *work, int *lwork, s *rwork, bint *bwork, int *info) nogil
cdef cgges_t *cgges_f

ctypedef void cggesx_t(char *jobvsl, char *jobvsr, char *sort, cselect2 *selctg, char *sense, int *n, c *a, int *lda, c *b, int *ldb, int *sdim, c *alpha, c *beta, c *vsl, int *ldvsl, c *vsr, int *ldvsr, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef cggesx_t *cggesx_f

ctypedef void cggev_t(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cggev_t *cggev_f

ctypedef void cggevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, int *ilo, int *ihi, s *lscale, s *rscale, s *abnrm, s *bbnrm, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, int *iwork, bint *bwork, int *info) nogil
cdef cggevx_t *cggevx_f

ctypedef void cggglm_t(int *n, int *m, int *p, c *a, int *lda, c *b, int *ldb, c *d, c *x, c *y, c *work, int *lwork, int *info) nogil
cdef cggglm_t *cggglm_f

ctypedef void cgghrd_t(char *compq, char *compz, int *n, int *ilo, int *ihi, c *a, int *lda, c *b, int *ldb, c *q, int *ldq, c *z, int *ldz, int *info) nogil
cdef cgghrd_t *cgghrd_f

ctypedef void cgglse_t(int *m, int *n, int *p, c *a, int *lda, c *b, int *ldb, c *c, c *d, c *x, c *work, int *lwork, int *info) nogil
cdef cgglse_t *cgglse_f

ctypedef void cggqrf_t(int *n, int *m, int *p, c *a, int *lda, c *taua, c *b, int *ldb, c *taub, c *work, int *lwork, int *info) nogil
cdef cggqrf_t *cggqrf_f

ctypedef void cggrqf_t(int *m, int *p, int *n, c *a, int *lda, c *taua, c *b, int *ldb, c *taub, c *work, int *lwork, int *info) nogil
cdef cggrqf_t *cggrqf_f

ctypedef void cggsvd_t(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, c *a, int *lda, c *b, int *ldb, s *alpha, s *beta, c *u, int *ldu, c *v, int *ldv, c *q, int *ldq, c *work, s *rwork, int *iwork, int *info) nogil
cdef cggsvd_t *cggsvd_f

ctypedef void cggsvp_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, c *a, int *lda, c *b, int *ldb, s *tola, s *tolb, int *k, int *l, c *u, int *ldu, c *v, int *ldv, c *q, int *ldq, int *iwork, s *rwork, c *tau, c *work, int *info) nogil
cdef cggsvp_t *cggsvp_f

ctypedef void cgtcon_t(char *norm, int *n, c *dl, c *d, c *du, c *du2, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef cgtcon_t *cgtcon_f

ctypedef void cgtrfs_t(char *trans, int *n, int *nrhs, c *dl, c *d, c *du, c *dlf, c *df, c *duf, c *du2, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgtrfs_t *cgtrfs_f

ctypedef void cgtsv_t(int *n, int *nrhs, c *dl, c *d, c *du, c *b, int *ldb, int *info) nogil
cdef cgtsv_t *cgtsv_f

ctypedef void cgtsvx_t(char *fact, char *trans, int *n, int *nrhs, c *dl, c *d, c *du, c *dlf, c *df, c *duf, c *du2, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgtsvx_t *cgtsvx_f

ctypedef void cgttrf_t(int *n, c *dl, c *d, c *du, c *du2, int *ipiv, int *info) nogil
cdef cgttrf_t *cgttrf_f

ctypedef void cgttrs_t(char *trans, int *n, int *nrhs, c *dl, c *d, c *du, c *du2, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgttrs_t *cgttrs_f

ctypedef void cgtts2_t(int *itrans, int *n, int *nrhs, c *dl, c *d, c *du, c *du2, int *ipiv, c *b, int *ldb) nogil
cdef cgtts2_t *cgtts2_f

ctypedef void chbev_t(char *jobz, char *uplo, int *n, int *kd, c *ab, int *ldab, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chbev_t *chbev_f

ctypedef void chbevd_t(char *jobz, char *uplo, int *n, int *kd, c *ab, int *ldab, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chbevd_t *chbevd_f

ctypedef void chbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, c *ab, int *ldab, c *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chbevx_t *chbevx_f

ctypedef void chbgst_t(char *vect, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, c *x, int *ldx, c *work, s *rwork, int *info) nogil
cdef chbgst_t *chbgst_f

ctypedef void chbgv_t(char *jobz, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chbgv_t *chbgv_f

ctypedef void chbgvd_t(char *jobz, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chbgvd_t *chbgvd_f

ctypedef void chbgvx_t(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, c *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chbgvx_t *chbgvx_f

ctypedef void chbtrd_t(char *vect, char *uplo, int *n, int *kd, c *ab, int *ldab, s *d, s *e, c *q, int *ldq, c *work, int *info) nogil
cdef chbtrd_t *chbtrd_f

ctypedef void checon_t(char *uplo, int *n, c *a, int *lda, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef checon_t *checon_f

ctypedef void cheev_t(char *jobz, char *uplo, int *n, c *a, int *lda, s *w, c *work, int *lwork, s *rwork, int *info) nogil
cdef cheev_t *cheev_f

ctypedef void cheevd_t(char *jobz, char *uplo, int *n, c *a, int *lda, s *w, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef cheevd_t *cheevd_f

ctypedef void cheevr_t(char *jobz, char *range, char *uplo, int *n, c *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, int *isuppz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef cheevr_t *cheevr_f

ctypedef void cheevx_t(char *jobz, char *range, char *uplo, int *n, c *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef cheevx_t *cheevx_f

ctypedef void chegs2_t(int *itype, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef chegs2_t *chegs2_f

ctypedef void chegst_t(int *itype, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef chegst_t *chegst_f

ctypedef void chegv_t(int *itype, char *jobz, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *w, c *work, int *lwork, s *rwork, int *info) nogil
cdef chegv_t *chegv_f

ctypedef void chegvd_t(int *itype, char *jobz, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *w, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chegvd_t *chegvd_f

ctypedef void chegvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chegvx_t *chegvx_f

ctypedef void cherfs_t(char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cherfs_t *cherfs_f

ctypedef void chesv_t(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, c *work, int *lwork, int *info) nogil
cdef chesv_t *chesv_f

ctypedef void chesvx_t(char *fact, char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, int *lwork, s *rwork, int *info) nogil
cdef chesvx_t *chesvx_f

ctypedef void chetd2_t(char *uplo, int *n, c *a, int *lda, s *d, s *e, c *tau, int *info) nogil
cdef chetd2_t *chetd2_f

ctypedef void chetf2_t(char *uplo, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef chetf2_t *chetf2_f

ctypedef void chetrd_t(char *uplo, int *n, c *a, int *lda, s *d, s *e, c *tau, c *work, int *lwork, int *info) nogil
cdef chetrd_t *chetrd_f

ctypedef void chetrf_t(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *lwork, int *info) nogil
cdef chetrf_t *chetrf_f

ctypedef void chetri_t(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *info) nogil
cdef chetri_t *chetri_f

ctypedef void chetrs_t(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef chetrs_t *chetrs_f

ctypedef void chgeqz_t(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, c *h, int *ldh, c *t, int *ldt, c *alpha, c *beta, c *q, int *ldq, c *z, int *ldz, c *work, int *lwork, s *rwork, int *info) nogil
cdef chgeqz_t *chgeqz_f

ctypedef void chpcon_t(char *uplo, int *n, c *ap, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef chpcon_t *chpcon_f

ctypedef void chpev_t(char *jobz, char *uplo, int *n, c *ap, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chpev_t *chpev_f

ctypedef void chpevd_t(char *jobz, char *uplo, int *n, c *ap, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chpevd_t *chpevd_f

ctypedef void chpevx_t(char *jobz, char *range, char *uplo, int *n, c *ap, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chpevx_t *chpevx_f

ctypedef void chpgst_t(int *itype, char *uplo, int *n, c *ap, c *bp, int *info) nogil
cdef chpgst_t *chpgst_f

ctypedef void chpgv_t(int *itype, char *jobz, char *uplo, int *n, c *ap, c *bp, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chpgv_t *chpgv_f

ctypedef void chpgvd_t(int *itype, char *jobz, char *uplo, int *n, c *ap, c *bp, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chpgvd_t *chpgvd_f

ctypedef void chpgvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, c *ap, c *bp, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chpgvx_t *chpgvx_f

ctypedef void chprfs_t(char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef chprfs_t *chprfs_f

ctypedef void chpsv_t(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef chpsv_t *chpsv_f

ctypedef void chpsvx_t(char *fact, char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef chpsvx_t *chpsvx_f

ctypedef void chptrd_t(char *uplo, int *n, c *ap, s *d, s *e, c *tau, int *info) nogil
cdef chptrd_t *chptrd_f

ctypedef void chptrf_t(char *uplo, int *n, c *ap, int *ipiv, int *info) nogil
cdef chptrf_t *chptrf_f

ctypedef void chptri_t(char *uplo, int *n, c *ap, int *ipiv, c *work, int *info) nogil
cdef chptri_t *chptri_f

ctypedef void chptrs_t(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef chptrs_t *chptrs_f

ctypedef void chsein_t(char *side, char *eigsrc, char *initv, bint *select, int *n, c *h, int *ldh, c *w, c *vl, int *ldvl, c *vr, int *ldvr, int *mm, int *m, c *work, s *rwork, int *ifaill, int *ifailr, int *info) nogil
cdef chsein_t *chsein_f

ctypedef void chseqr_t(char *job, char *compz, int *n, int *ilo, int *ihi, c *h, int *ldh, c *w, c *z, int *ldz, c *work, int *lwork, int *info) nogil
cdef chseqr_t *chseqr_f

ctypedef void clacn2_t(int *n, c *v, c *x, s *est, int *kase, int *isave) nogil
cdef clacn2_t *clacn2_f

ctypedef void clacon_t(int *n, c *v, c *x, s *est, int *kase) nogil
cdef clacon_t *clacon_f

ctypedef s clangb_t(char *norm, int *n, int *kl, int *ku, c *ab, int *ldab, s *work) nogil
cdef clangb_t *clangb_f

ctypedef s clange_t(char *norm, int *m, int *n, c *a, int *lda, s *work) nogil
cdef clange_t *clange_f

ctypedef s clangt_t(char *norm, int *n, c *dl, c *d, c *du) nogil
cdef clangt_t *clangt_f

ctypedef s clanhb_t(char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef clanhb_t *clanhb_f

ctypedef s clanhe_t(char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil
cdef clanhe_t *clanhe_f

ctypedef s clanhp_t(char *norm, char *uplo, int *n, c *ap, s *work) nogil
cdef clanhp_t *clanhp_f

ctypedef s clanhs_t(char *norm, int *n, c *a, int *lda, s *work) nogil
cdef clanhs_t *clanhs_f

ctypedef s clanht_t(char *norm, int *n, s *d, c *e) nogil
cdef clanht_t *clanht_f

ctypedef s clansb_t(char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef clansb_t *clansb_f

ctypedef s clansp_t(char *norm, char *uplo, int *n, c *ap, s *work) nogil
cdef clansp_t *clansp_f

ctypedef s clansy_t(char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil
cdef clansy_t *clansy_f

ctypedef s clantb_t(char *norm, char *uplo, char *diag, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef clantb_t *clantb_f

ctypedef s clantp_t(char *norm, char *uplo, char *diag, int *n, c *ap, s *work) nogil
cdef clantp_t *clantp_f

ctypedef s clantr_t(char *norm, char *uplo, char *diag, int *m, int *n, c *a, int *lda, s *work) nogil
cdef clantr_t *clantr_f

ctypedef void clarf_t(char *side, int *m, int *n, c *v, int *incv, c *tau, c *c, int *ldc, c *work) nogil
cdef clarf_t *clarf_f

ctypedef void clarz_t(char *side, int *m, int *n, int *l, c *v, int *incv, c *tau, c *c, int *ldc, c *work) nogil
cdef clarz_t *clarz_f

ctypedef void claswp_t(int *n, c *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef claswp_t *claswp_f

ctypedef void clauum_t(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef clauum_t *clauum_f

ctypedef void cpbcon_t(char *uplo, int *n, int *kd, c *ab, int *ldab, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cpbcon_t *cpbcon_f

ctypedef void cpbequ_t(char *uplo, int *n, int *kd, c *ab, int *ldab, s *s, s *scond, s *amax, int *info) nogil
cdef cpbequ_t *cpbequ_f

ctypedef void cpbrfs_t(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cpbrfs_t *cpbrfs_f

ctypedef void cpbstf_t(char *uplo, int *n, int *kd, c *ab, int *ldab, int *info) nogil
cdef cpbstf_t *cpbstf_f

ctypedef void cpbsv_t(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
cdef cpbsv_t *cpbsv_f

ctypedef void cpbsvx_t(char *fact, char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, char *equed, s *s, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cpbsvx_t *cpbsvx_f

ctypedef void cpbtf2_t(char *uplo, int *n, int *kd, c *ab, int *ldab, int *info) nogil
cdef cpbtf2_t *cpbtf2_f

ctypedef void cpbtrf_t(char *uplo, int *n, int *kd, c *ab, int *ldab, int *info) nogil
cdef cpbtrf_t *cpbtrf_f

ctypedef void cpbtrs_t(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
cdef cpbtrs_t *cpbtrs_f

ctypedef void cpocon_t(char *uplo, int *n, c *a, int *lda, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cpocon_t *cpocon_f

ctypedef void cpoequ_t(int *n, c *a, int *lda, s *s, s *scond, s *amax, int *info) nogil
cdef cpoequ_t *cpoequ_f

ctypedef void cporfs_t(char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cporfs_t *cporfs_f

ctypedef void cposv_t(char *uplo, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef cposv_t *cposv_f

ctypedef void cposvx_t(char *fact, char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, char *equed, s *s, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cposvx_t *cposvx_f

ctypedef void cpotf2_t(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef cpotf2_t *cpotf2_f

ctypedef void cpotrf_t(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef cpotrf_t *cpotrf_f

ctypedef void cpotri_t(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef cpotri_t *cpotri_f

ctypedef void cpotrs_t(char *uplo, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef cpotrs_t *cpotrs_f

ctypedef void cppcon_t(char *uplo, int *n, c *ap, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cppcon_t *cppcon_f

ctypedef void cppequ_t(char *uplo, int *n, c *ap, s *s, s *scond, s *amax, int *info) nogil
cdef cppequ_t *cppequ_f

ctypedef void cpprfs_t(char *uplo, int *n, int *nrhs, c *ap, c *afp, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cpprfs_t *cpprfs_f

ctypedef void cppsv_t(char *uplo, int *n, int *nrhs, c *ap, c *b, int *ldb, int *info) nogil
cdef cppsv_t *cppsv_f

ctypedef void cppsvx_t(char *fact, char *uplo, int *n, int *nrhs, c *ap, c *afp, char *equed, s *s, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cppsvx_t *cppsvx_f

ctypedef void cpptrf_t(char *uplo, int *n, c *ap, int *info) nogil
cdef cpptrf_t *cpptrf_f

ctypedef void cpptri_t(char *uplo, int *n, c *ap, int *info) nogil
cdef cpptri_t *cpptri_f

ctypedef void cpptrs_t(char *uplo, int *n, int *nrhs, c *ap, c *b, int *ldb, int *info) nogil
cdef cpptrs_t *cpptrs_f

ctypedef void cptcon_t(int *n, s *d, c *e, s *anorm, s *rcond, s *rwork, int *info) nogil
cdef cptcon_t *cptcon_f

ctypedef void cpteqr_t(char *compz, int *n, s *d, s *e, c *z, int *ldz, s *work, int *info) nogil
cdef cpteqr_t *cpteqr_f

ctypedef void cptrfs_t(char *uplo, int *n, int *nrhs, s *d, c *e, s *df, c *ef, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cptrfs_t *cptrfs_f

ctypedef void cptsv_t(int *n, int *nrhs, s *d, c *e, c *b, int *ldb, int *info) nogil
cdef cptsv_t *cptsv_f

ctypedef void cptsvx_t(char *fact, int *n, int *nrhs, s *d, c *e, s *df, c *ef, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cptsvx_t *cptsvx_f

ctypedef void cpttrf_t(int *n, s *d, c *e, int *info) nogil
cdef cpttrf_t *cpttrf_f

ctypedef void cpttrs_t(char *uplo, int *n, int *nrhs, s *d, c *e, c *b, int *ldb, int *info) nogil
cdef cpttrs_t *cpttrs_f

ctypedef void cptts2_t(int *iuplo, int *n, int *nrhs, s *d, c *e, c *b, int *ldb) nogil
cdef cptts2_t *cptts2_f

ctypedef void crot_t(int *n, c *cx, int *incx, c *cy, int *incy, s *c, c *s) nogil
cdef crot_t *crot_f

ctypedef void cspcon_t(char *uplo, int *n, c *ap, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef cspcon_t *cspcon_f

ctypedef void cspmv_t(char *uplo, int *n, c *alpha, c *ap, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef cspmv_t *cspmv_f

ctypedef void cspr_t(char *uplo, int *n, c *alpha, c *x, int *incx, c *ap) nogil
cdef cspr_t *cspr_f

ctypedef void csprfs_t(char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef csprfs_t *csprfs_f

ctypedef void cspsv_t(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cspsv_t *cspsv_f

ctypedef void cspsvx_t(char *fact, char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cspsvx_t *cspsvx_f

ctypedef void csptrf_t(char *uplo, int *n, c *ap, int *ipiv, int *info) nogil
cdef csptrf_t *csptrf_f

ctypedef void csptri_t(char *uplo, int *n, c *ap, int *ipiv, c *work, int *info) nogil
cdef csptri_t *csptri_f

ctypedef void csptrs_t(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef csptrs_t *csptrs_f

ctypedef void csrscl_t(int *n, s *sa, c *sx, int *incx) nogil
cdef csrscl_t *csrscl_f

ctypedef void cstedc_t(char *compz, int *n, s *d, s *e, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef cstedc_t *cstedc_f

ctypedef void cstegr_t(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef cstegr_t *cstegr_f

ctypedef void cstein_t(int *n, s *d, s *e, int *m, s *w, int *iblock, int *isplit, c *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef cstein_t *cstein_f

ctypedef void cstemr_t(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, int *m, s *w, c *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef cstemr_t *cstemr_f

ctypedef void csteqr_t(char *compz, int *n, s *d, s *e, c *z, int *ldz, s *work, int *info) nogil
cdef csteqr_t *csteqr_f

ctypedef void csycon_t(char *uplo, int *n, c *a, int *lda, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef csycon_t *csycon_f

ctypedef void csymv_t(char *uplo, int *n, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef csymv_t *csymv_f

ctypedef void csyr_t(char *uplo, int *n, c *alpha, c *x, int *incx, c *a, int *lda) nogil
cdef csyr_t *csyr_f

ctypedef void csyrfs_t(char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef csyrfs_t *csyrfs_f

ctypedef void csysv_t(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, c *work, int *lwork, int *info) nogil
cdef csysv_t *csysv_f

ctypedef void csysvx_t(char *fact, char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, int *lwork, s *rwork, int *info) nogil
cdef csysvx_t *csysvx_f

ctypedef void csytf2_t(char *uplo, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef csytf2_t *csytf2_f

ctypedef void csytrf_t(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *lwork, int *info) nogil
cdef csytrf_t *csytrf_f

ctypedef void csytri_t(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *info) nogil
cdef csytri_t *csytri_f

ctypedef void csytrs_t(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef csytrs_t *csytrs_f

ctypedef void ctbcon_t(char *norm, char *uplo, char *diag, int *n, int *kd, c *ab, int *ldab, s *rcond, c *work, s *rwork, int *info) nogil
cdef ctbcon_t *ctbcon_f

ctypedef void ctbrfs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef ctbrfs_t *ctbrfs_f

ctypedef void ctbtrs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
cdef ctbtrs_t *ctbtrs_f

ctypedef void ctgevc_t(char *side, char *howmny, bint *select, int *n, c *s, int *lds, c *p, int *ldp, c *vl, int *ldvl, c *vr, int *ldvr, int *mm, int *m, c *work, s *rwork, int *info) nogil
cdef ctgevc_t *ctgevc_f

ctypedef void ctgex2_t(bint *wantq, bint *wantz, int *n, c *a, int *lda, c *b, int *ldb, c *q, int *ldq, c *z, int *ldz, int *j1, int *info) nogil
cdef ctgex2_t *ctgex2_f

ctypedef void ctgexc_t(bint *wantq, bint *wantz, int *n, c *a, int *lda, c *b, int *ldb, c *q, int *ldq, c *z, int *ldz, int *ifst, int *ilst, int *info) nogil
cdef ctgexc_t *ctgexc_f

ctypedef void ctgsen_t(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *q, int *ldq, c *z, int *ldz, int *m, s *pl, s *pr, s *dif, c *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ctgsen_t *ctgsen_f

ctypedef void ctgsja_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, c *a, int *lda, c *b, int *ldb, s *tola, s *tolb, s *alpha, s *beta, c *u, int *ldu, c *v, int *ldv, c *q, int *ldq, c *work, int *ncycle, int *info) nogil
cdef ctgsja_t *ctgsja_f

ctypedef void ctgsna_t(char *job, char *howmny, bint *select, int *n, c *a, int *lda, c *b, int *ldb, c *vl, int *ldvl, c *vr, int *ldvr, s *s, s *dif, int *mm, int *m, c *work, int *lwork, int *iwork, int *info) nogil
cdef ctgsna_t *ctgsna_f

ctypedef void ctgsy2_t(char *trans, int *ijob, int *m, int *n, c *a, int *lda, c *b, int *ldb, c *c, int *ldc, c *d, int *ldd, c *e, int *lde, c *f, int *ldf, s *scale, s *rdsum, s *rdscal, int *info) nogil
cdef ctgsy2_t *ctgsy2_f

ctypedef void ctgsyl_t(char *trans, int *ijob, int *m, int *n, c *a, int *lda, c *b, int *ldb, c *c, int *ldc, c *d, int *ldd, c *e, int *lde, c *f, int *ldf, s *scale, s *dif, c *work, int *lwork, int *iwork, int *info) nogil
cdef ctgsyl_t *ctgsyl_f

ctypedef void ctpcon_t(char *norm, char *uplo, char *diag, int *n, c *ap, s *rcond, c *work, s *rwork, int *info) nogil
cdef ctpcon_t *ctpcon_f

ctypedef void ctprfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *ap, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef ctprfs_t *ctprfs_f

ctypedef void ctptri_t(char *uplo, char *diag, int *n, c *ap, int *info) nogil
cdef ctptri_t *ctptri_f

ctypedef void ctptrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *ap, c *b, int *ldb, int *info) nogil
cdef ctptrs_t *ctptrs_f

ctypedef void ctrcon_t(char *norm, char *uplo, char *diag, int *n, c *a, int *lda, s *rcond, c *work, s *rwork, int *info) nogil
cdef ctrcon_t *ctrcon_f

ctypedef void ctrevc_t(char *side, char *howmny, bint *select, int *n, c *t, int *ldt, c *vl, int *ldvl, c *vr, int *ldvr, int *mm, int *m, c *work, s *rwork, int *info) nogil
cdef ctrevc_t *ctrevc_f

ctypedef void ctrexc_t(char *compq, int *n, c *t, int *ldt, c *q, int *ldq, int *ifst, int *ilst, int *info) nogil
cdef ctrexc_t *ctrexc_f

ctypedef void ctrrfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef ctrrfs_t *ctrrfs_f

ctypedef void ctrsen_t(char *job, char *compq, bint *select, int *n, c *t, int *ldt, c *q, int *ldq, c *w, int *m, s *s, s *sep, c *work, int *lwork, int *info) nogil
cdef ctrsen_t *ctrsen_f

ctypedef void ctrsna_t(char *job, char *howmny, bint *select, int *n, c *t, int *ldt, c *vl, int *ldvl, c *vr, int *ldvr, s *s, s *sep, int *mm, int *m, c *work, int *ldwork, s *rwork, int *info) nogil
cdef ctrsna_t *ctrsna_f

ctypedef void ctrsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, c *a, int *lda, c *b, int *ldb, c *c, int *ldc, s *scale, int *info) nogil
cdef ctrsyl_t *ctrsyl_f

ctypedef void ctrti2_t(char *uplo, char *diag, int *n, c *a, int *lda, int *info) nogil
cdef ctrti2_t *ctrti2_f

ctypedef void ctrtri_t(char *uplo, char *diag, int *n, c *a, int *lda, int *info) nogil
cdef ctrtri_t *ctrtri_f

ctypedef void ctrtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef ctrtrs_t *ctrtrs_f

ctypedef void ctzrqf_t(int *m, int *n, c *a, int *lda, c *tau, int *info) nogil
cdef ctzrqf_t *ctzrqf_f

ctypedef void ctzrzf_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef ctzrzf_t *ctzrzf_f

ctypedef void cung2l_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cung2l_t *cung2l_f

ctypedef void cung2r_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cung2r_t *cung2r_f

ctypedef void cungbr_t(char *vect, int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungbr_t *cungbr_f

ctypedef void cunghr_t(int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cunghr_t *cunghr_f

ctypedef void cungl2_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cungl2_t *cungl2_f

ctypedef void cunglq_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cunglq_t *cunglq_f

ctypedef void cungql_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungql_t *cungql_f

ctypedef void cungqr_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungqr_t *cungqr_f

ctypedef void cungr2_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cungr2_t *cungr2_f

ctypedef void cungrq_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungrq_t *cungrq_f

ctypedef void cungtr_t(char *uplo, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungtr_t *cungtr_f

ctypedef void cunm2l_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunm2l_t *cunm2l_f

ctypedef void cunm2r_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunm2r_t *cunm2r_f

ctypedef void cunmbr_t(char *vect, char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmbr_t *cunmbr_f

ctypedef void cunmhr_t(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmhr_t *cunmhr_f

ctypedef void cunml2_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunml2_t *cunml2_f

ctypedef void cunmlq_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmlq_t *cunmlq_f

ctypedef void cunmql_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmql_t *cunmql_f

ctypedef void cunmqr_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmqr_t *cunmqr_f

ctypedef void cunmr2_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunmr2_t *cunmr2_f

ctypedef void cunmr3_t(char *side, char *trans, int *m, int *n, int *k, int *l, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunmr3_t *cunmr3_f

ctypedef void cunmrq_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmrq_t *cunmrq_f

ctypedef void cunmrz_t(char *side, char *trans, int *m, int *n, int *k, int *l, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmrz_t *cunmrz_f

ctypedef void cunmtr_t(char *side, char *uplo, char *trans, int *m, int *n, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmtr_t *cunmtr_f

ctypedef void cupgtr_t(char *uplo, int *n, c *ap, c *tau, c *q, int *ldq, c *work, int *info) nogil
cdef cupgtr_t *cupgtr_f

ctypedef void cupmtr_t(char *side, char *uplo, char *trans, int *m, int *n, c *ap, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cupmtr_t *cupmtr_f

ctypedef void dbdsdc_t(char *uplo, char *compq, int *n, d *d, d *e, d *u, int *ldu, d *vt, int *ldvt, d *q, int *iq, d *work, int *iwork, int *info) nogil
cdef dbdsdc_t *dbdsdc_f

ctypedef void dbdsqr_t(char *uplo, int *n, int *ncvt, int *nru, int *ncc, d *d, d *e, d *vt, int *ldvt, d *u, int *ldu, d *c, int *ldc, d *work, int *info) nogil
cdef dbdsqr_t *dbdsqr_f

ctypedef void ddisna_t(char *job, int *m, int *n, d *d, d *sep, int *info) nogil
cdef ddisna_t *ddisna_f

ctypedef void dgbbrd_t(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, d *ab, int *ldab, d *d, d *e, d *q, int *ldq, d *pt, int *ldpt, d *c, int *ldc, d *work, int *info) nogil
cdef dgbbrd_t *dgbbrd_f

ctypedef void dgbcon_t(char *norm, int *n, int *kl, int *ku, d *ab, int *ldab, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dgbcon_t *dgbcon_f

ctypedef void dgbequ_t(int *m, int *n, int *kl, int *ku, d *ab, int *ldab, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef dgbequ_t *dgbequ_f

ctypedef void dgbrfs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgbrfs_t *dgbrfs_f

ctypedef void dgbsv_t(int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgbsv_t *dgbsv_f

ctypedef void dgbsvx_t(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, int *ipiv, char *equed, d *r, d *c, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgbsvx_t *dgbsvx_f

ctypedef void dgbtf2_t(int *m, int *n, int *kl, int *ku, d *ab, int *ldab, int *ipiv, int *info) nogil
cdef dgbtf2_t *dgbtf2_f

ctypedef void dgbtrf_t(int *m, int *n, int *kl, int *ku, d *ab, int *ldab, int *ipiv, int *info) nogil
cdef dgbtrf_t *dgbtrf_f

ctypedef void dgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgbtrs_t *dgbtrs_f

ctypedef void dgebak_t(char *job, char *side, int *n, int *ilo, int *ihi, d *scale, int *m, d *v, int *ldv, int *info) nogil
cdef dgebak_t *dgebak_f

ctypedef void dgebal_t(char *job, int *n, d *a, int *lda, int *ilo, int *ihi, d *scale, int *info) nogil
cdef dgebal_t *dgebal_f

ctypedef void dgebd2_t(int *m, int *n, d *a, int *lda, d *d, d *e, d *tauq, d *taup, d *work, int *info) nogil
cdef dgebd2_t *dgebd2_f

ctypedef void dgebrd_t(int *m, int *n, d *a, int *lda, d *d, d *e, d *tauq, d *taup, d *work, int *lwork, int *info) nogil
cdef dgebrd_t *dgebrd_f

ctypedef void dgecon_t(char *norm, int *n, d *a, int *lda, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dgecon_t *dgecon_f

ctypedef void dgeequ_t(int *m, int *n, d *a, int *lda, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef dgeequ_t *dgeequ_f

ctypedef void dgees_t(char *jobvs, char *sort, dselect2 *select, int *n, d *a, int *lda, int *sdim, d *wr, d *wi, d *vs, int *ldvs, d *work, int *lwork, bint *bwork, int *info) nogil
cdef dgees_t *dgees_f

ctypedef void dgeesx_t(char *jobvs, char *sort, dselect2 *select, char *sense, int *n, d *a, int *lda, int *sdim, d *wr, d *wi, d *vs, int *ldvs, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef dgeesx_t *dgeesx_f

ctypedef void dgeev_t(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *wr, d *wi, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
cdef dgeev_t *dgeev_f

ctypedef void dgeevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, d *a, int *lda, d *wr, d *wi, d *vl, int *ldvl, d *vr, int *ldvr, int *ilo, int *ihi, d *scale, d *abnrm, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, int *info) nogil
cdef dgeevx_t *dgeevx_f

ctypedef void dgegs_t(char *jobvsl, char *jobvsr, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vsl, int *ldvsl, d *vsr, int *ldvsr, d *work, int *lwork, int *info) nogil
cdef dgegs_t *dgegs_f

ctypedef void dgegv_t(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
cdef dgegv_t *dgegv_f

ctypedef void dgehd2_t(int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgehd2_t *dgehd2_f

ctypedef void dgehrd_t(int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgehrd_t *dgehrd_f

ctypedef void dgelq2_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgelq2_t *dgelq2_f

ctypedef void dgelqf_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgelqf_t *dgelqf_f

ctypedef void dgels_t(char *trans, int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *work, int *lwork, int *info) nogil
cdef dgels_t *dgels_f

ctypedef void dgelsd_t(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *s, d *rcond, int *rank, d *work, int *lwork, int *iwork, int *info) nogil
cdef dgelsd_t *dgelsd_f

ctypedef void dgelss_t(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *s, d *rcond, int *rank, d *work, int *lwork, int *info) nogil
cdef dgelss_t *dgelss_f

ctypedef void dgelsx_t(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *jpvt, d *rcond, int *rank, d *work, int *info) nogil
cdef dgelsx_t *dgelsx_f

ctypedef void dgelsy_t(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *jpvt, d *rcond, int *rank, d *work, int *lwork, int *info) nogil
cdef dgelsy_t *dgelsy_f

ctypedef void dgeql2_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgeql2_t *dgeql2_f

ctypedef void dgeqlf_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgeqlf_t *dgeqlf_f

ctypedef void dgeqp3_t(int *m, int *n, d *a, int *lda, int *jpvt, d *tau, d *work, int *lwork, int *info) nogil
cdef dgeqp3_t *dgeqp3_f

ctypedef void dgeqpf_t(int *m, int *n, d *a, int *lda, int *jpvt, d *tau, d *work, int *info) nogil
cdef dgeqpf_t *dgeqpf_f

ctypedef void dgeqr2_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgeqr2_t *dgeqr2_f

ctypedef void dgeqrf_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgeqrf_t *dgeqrf_f

ctypedef void dgerfs_t(char *trans, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgerfs_t *dgerfs_f

ctypedef void dgerq2_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgerq2_t *dgerq2_f

ctypedef void dgerqf_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgerqf_t *dgerqf_f

ctypedef void dgesc2_t(int *n, d *a, int *lda, d *rhs, int *ipiv, int *jpiv, d *scale) nogil
cdef dgesc2_t *dgesc2_f

ctypedef void dgesdd_t(char *jobz, int *m, int *n, d *a, int *lda, d *s, d *u, int *ldu, d *vt, int *ldvt, d *work, int *lwork, int *iwork, int *info) nogil
cdef dgesdd_t *dgesdd_f

ctypedef void dgesv_t(int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgesv_t *dgesv_f

ctypedef void dgesvd_t(char *jobu, char *jobvt, int *m, int *n, d *a, int *lda, d *s, d *u, int *ldu, d *vt, int *ldvt, d *work, int *lwork, int *info) nogil
cdef dgesvd_t *dgesvd_f

ctypedef void dgesvx_t(char *fact, char *trans, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, char *equed, d *r, d *c, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgesvx_t *dgesvx_f

ctypedef void dgetc2_t(int *n, d *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef dgetc2_t *dgetc2_f

ctypedef void dgetf2_t(int *m, int *n, d *a, int *lda, int *ipiv, int *info) nogil
cdef dgetf2_t *dgetf2_f

ctypedef void dgetrf_t(int *m, int *n, d *a, int *lda, int *ipiv, int *info) nogil
cdef dgetrf_t *dgetrf_f

ctypedef void dgetri_t(int *n, d *a, int *lda, int *ipiv, d *work, int *lwork, int *info) nogil
cdef dgetri_t *dgetri_f

ctypedef void dgetrs_t(char *trans, int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgetrs_t *dgetrs_f

ctypedef void dggbak_t(char *job, char *side, int *n, int *ilo, int *ihi, d *lscale, d *rscale, int *m, d *v, int *ldv, int *info) nogil
cdef dggbak_t *dggbak_f

ctypedef void dggbal_t(char *job, int *n, d *a, int *lda, d *b, int *ldb, int *ilo, int *ihi, d *lscale, d *rscale, d *work, int *info) nogil
cdef dggbal_t *dggbal_f

ctypedef void dgges_t(char *jobvsl, char *jobvsr, char *sort, dselect3 *selctg, int *n, d *a, int *lda, d *b, int *ldb, int *sdim, d *alphar, d *alphai, d *beta, d *vsl, int *ldvsl, d *vsr, int *ldvsr, d *work, int *lwork, bint *bwork, int *info) nogil
cdef dgges_t *dgges_f

ctypedef void dggesx_t(char *jobvsl, char *jobvsr, char *sort, dselect3 *selctg, char *sense, int *n, d *a, int *lda, d *b, int *ldb, int *sdim, d *alphar, d *alphai, d *beta, d *vsl, int *ldvsl, d *vsr, int *ldvsr, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef dggesx_t *dggesx_f

ctypedef void dggev_t(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
cdef dggev_t *dggev_f

ctypedef void dggevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, int *ilo, int *ihi, d *lscale, d *rscale, d *abnrm, d *bbnrm, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, bint *bwork, int *info) nogil
cdef dggevx_t *dggevx_f

ctypedef void dggglm_t(int *n, int *m, int *p, d *a, int *lda, d *b, int *ldb, d *d, d *x, d *y, d *work, int *lwork, int *info) nogil
cdef dggglm_t *dggglm_f

ctypedef void dgghrd_t(char *compq, char *compz, int *n, int *ilo, int *ihi, d *a, int *lda, d *b, int *ldb, d *q, int *ldq, d *z, int *ldz, int *info) nogil
cdef dgghrd_t *dgghrd_f

ctypedef void dgglse_t(int *m, int *n, int *p, d *a, int *lda, d *b, int *ldb, d *c, d *d, d *x, d *work, int *lwork, int *info) nogil
cdef dgglse_t *dgglse_f

ctypedef void dggqrf_t(int *n, int *m, int *p, d *a, int *lda, d *taua, d *b, int *ldb, d *taub, d *work, int *lwork, int *info) nogil
cdef dggqrf_t *dggqrf_f

ctypedef void dggrqf_t(int *m, int *p, int *n, d *a, int *lda, d *taua, d *b, int *ldb, d *taub, d *work, int *lwork, int *info) nogil
cdef dggrqf_t *dggrqf_f

ctypedef void dggsvd_t(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, d *a, int *lda, d *b, int *ldb, d *alpha, d *beta, d *u, int *ldu, d *v, int *ldv, d *q, int *ldq, d *work, int *iwork, int *info) nogil
cdef dggsvd_t *dggsvd_f

ctypedef void dggsvp_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, d *a, int *lda, d *b, int *ldb, d *tola, d *tolb, int *k, int *l, d *u, int *ldu, d *v, int *ldv, d *q, int *ldq, int *iwork, d *tau, d *work, int *info) nogil
cdef dggsvp_t *dggsvp_f

ctypedef void dgtcon_t(char *norm, int *n, d *dl, d *d, d *du, d *du2, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dgtcon_t *dgtcon_f

ctypedef void dgtrfs_t(char *trans, int *n, int *nrhs, d *dl, d *d, d *du, d *dlf, d *df, d *duf, d *du2, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgtrfs_t *dgtrfs_f

ctypedef void dgtsv_t(int *n, int *nrhs, d *dl, d *d, d *du, d *b, int *ldb, int *info) nogil
cdef dgtsv_t *dgtsv_f

ctypedef void dgtsvx_t(char *fact, char *trans, int *n, int *nrhs, d *dl, d *d, d *du, d *dlf, d *df, d *duf, d *du2, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgtsvx_t *dgtsvx_f

ctypedef void dgttrf_t(int *n, d *dl, d *d, d *du, d *du2, int *ipiv, int *info) nogil
cdef dgttrf_t *dgttrf_f

ctypedef void dgttrs_t(char *trans, int *n, int *nrhs, d *dl, d *d, d *du, d *du2, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgttrs_t *dgttrs_f

ctypedef void dgtts2_t(int *itrans, int *n, int *nrhs, d *dl, d *d, d *du, d *du2, int *ipiv, d *b, int *ldb) nogil
cdef dgtts2_t *dgtts2_f

ctypedef void dhgeqz_t(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, d *h, int *ldh, d *t, int *ldt, d *alphar, d *alphai, d *beta, d *q, int *ldq, d *z, int *ldz, d *work, int *lwork, int *info) nogil
cdef dhgeqz_t *dhgeqz_f

ctypedef void dhsein_t(char *side, char *eigsrc, char *initv, bint *select, int *n, d *h, int *ldh, d *wr, d *wi, d *vl, int *ldvl, d *vr, int *ldvr, int *mm, int *m, d *work, int *ifaill, int *ifailr, int *info) nogil
cdef dhsein_t *dhsein_f

ctypedef void dhseqr_t(char *job, char *compz, int *n, int *ilo, int *ihi, d *h, int *ldh, d *wr, d *wi, d *z, int *ldz, d *work, int *lwork, int *info) nogil
cdef dhseqr_t *dhseqr_f

ctypedef bint disnan_t(d *din) nogil
cdef disnan_t *disnan_f

ctypedef void dlacn2_t(int *n, d *v, d *x, int *isgn, d *est, int *kase, int *isave) nogil
cdef dlacn2_t *dlacn2_f

ctypedef void dlacon_t(int *n, d *v, d *x, int *isgn, d *est, int *kase) nogil
cdef dlacon_t *dlacon_f

ctypedef d dlamch_t(char *cmach) nogil
cdef dlamch_t *dlamch_f

ctypedef d dlangb_t(char *norm, int *n, int *kl, int *ku, d *ab, int *ldab, d *work) nogil
cdef dlangb_t *dlangb_f

ctypedef d dlange_t(char *norm, int *m, int *n, d *a, int *lda, d *work) nogil
cdef dlange_t *dlange_f

ctypedef d dlangt_t(char *norm, int *n, d *dl, d *d, d *du) nogil
cdef dlangt_t *dlangt_f

ctypedef d dlanhs_t(char *norm, int *n, d *a, int *lda, d *work) nogil
cdef dlanhs_t *dlanhs_f

ctypedef d dlansb_t(char *norm, char *uplo, int *n, int *k, d *ab, int *ldab, d *work) nogil
cdef dlansb_t *dlansb_f

ctypedef d dlansp_t(char *norm, char *uplo, int *n, d *ap, d *work) nogil
cdef dlansp_t *dlansp_f

ctypedef d dlanst_t(char *norm, int *n, d *d, d *e) nogil
cdef dlanst_t *dlanst_f

ctypedef d dlansy_t(char *norm, char *uplo, int *n, d *a, int *lda, d *work) nogil
cdef dlansy_t *dlansy_f

ctypedef d dlantb_t(char *norm, char *uplo, char *diag, int *n, int *k, d *ab, int *ldab, d *work) nogil
cdef dlantb_t *dlantb_f

ctypedef d dlantp_t(char *norm, char *uplo, char *diag, int *n, d *ap, d *work) nogil
cdef dlantp_t *dlantp_f

ctypedef d dlantr_t(char *norm, char *uplo, char *diag, int *m, int *n, d *a, int *lda, d *work) nogil
cdef dlantr_t *dlantr_f

ctypedef void dlanv2_t(d *a, d *b, d *c, d *d, d *rt1r, d *rt1i, d *rt2r, d *rt2i, d *cs, d *sn) nogil
cdef dlanv2_t *dlanv2_f

ctypedef void dlarf_t(char *side, int *m, int *n, d *v, int *incv, d *tau, d *c, int *ldc, d *work) nogil
cdef dlarf_t *dlarf_f

ctypedef void dlarz_t(char *side, int *m, int *n, int *l, d *v, int *incv, d *tau, d *c, int *ldc, d *work) nogil
cdef dlarz_t *dlarz_f

ctypedef void dlaswp_t(int *n, d *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef dlaswp_t *dlaswp_f

ctypedef void dlauum_t(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dlauum_t *dlauum_f

ctypedef void dopgtr_t(char *uplo, int *n, d *ap, d *tau, d *q, int *ldq, d *work, int *info) nogil
cdef dopgtr_t *dopgtr_f

ctypedef void dopmtr_t(char *side, char *uplo, char *trans, int *m, int *n, d *ap, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dopmtr_t *dopmtr_f

ctypedef void dorg2l_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorg2l_t *dorg2l_f

ctypedef void dorg2r_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorg2r_t *dorg2r_f

ctypedef void dorgbr_t(char *vect, int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgbr_t *dorgbr_f

ctypedef void dorghr_t(int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorghr_t *dorghr_f

ctypedef void dorgl2_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorgl2_t *dorgl2_f

ctypedef void dorglq_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorglq_t *dorglq_f

ctypedef void dorgql_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgql_t *dorgql_f

ctypedef void dorgqr_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgqr_t *dorgqr_f

ctypedef void dorgr2_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorgr2_t *dorgr2_f

ctypedef void dorgrq_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgrq_t *dorgrq_f

ctypedef void dorgtr_t(char *uplo, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgtr_t *dorgtr_f

ctypedef void dorm2l_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dorm2l_t *dorm2l_f

ctypedef void dorm2r_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dorm2r_t *dorm2r_f

ctypedef void dormbr_t(char *vect, char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormbr_t *dormbr_f

ctypedef void dormhr_t(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormhr_t *dormhr_f

ctypedef void dorml2_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dorml2_t *dorml2_f

ctypedef void dormlq_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormlq_t *dormlq_f

ctypedef void dormql_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormql_t *dormql_f

ctypedef void dormqr_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormqr_t *dormqr_f

ctypedef void dormr2_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dormr2_t *dormr2_f

ctypedef void dormr3_t(char *side, char *trans, int *m, int *n, int *k, int *l, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dormr3_t *dormr3_f

ctypedef void dormrq_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormrq_t *dormrq_f

ctypedef void dormrz_t(char *side, char *trans, int *m, int *n, int *k, int *l, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormrz_t *dormrz_f

ctypedef void dormtr_t(char *side, char *uplo, char *trans, int *m, int *n, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormtr_t *dormtr_f

ctypedef void dpbcon_t(char *uplo, int *n, int *kd, d *ab, int *ldab, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dpbcon_t *dpbcon_f

ctypedef void dpbequ_t(char *uplo, int *n, int *kd, d *ab, int *ldab, d *s, d *scond, d *amax, int *info) nogil
cdef dpbequ_t *dpbequ_f

ctypedef void dpbrfs_t(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dpbrfs_t *dpbrfs_f

ctypedef void dpbstf_t(char *uplo, int *n, int *kd, d *ab, int *ldab, int *info) nogil
cdef dpbstf_t *dpbstf_f

ctypedef void dpbsv_t(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
cdef dpbsv_t *dpbsv_f

ctypedef void dpbsvx_t(char *fact, char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, char *equed, d *s, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dpbsvx_t *dpbsvx_f

ctypedef void dpbtf2_t(char *uplo, int *n, int *kd, d *ab, int *ldab, int *info) nogil
cdef dpbtf2_t *dpbtf2_f

ctypedef void dpbtrf_t(char *uplo, int *n, int *kd, d *ab, int *ldab, int *info) nogil
cdef dpbtrf_t *dpbtrf_f

ctypedef void dpbtrs_t(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
cdef dpbtrs_t *dpbtrs_f

ctypedef void dpocon_t(char *uplo, int *n, d *a, int *lda, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dpocon_t *dpocon_f

ctypedef void dpoequ_t(int *n, d *a, int *lda, d *s, d *scond, d *amax, int *info) nogil
cdef dpoequ_t *dpoequ_f

ctypedef void dporfs_t(char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dporfs_t *dporfs_f

ctypedef void dposv_t(char *uplo, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dposv_t *dposv_f

ctypedef void dposvx_t(char *fact, char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, char *equed, d *s, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dposvx_t *dposvx_f

ctypedef void dpotf2_t(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dpotf2_t *dpotf2_f

ctypedef void dpotrf_t(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dpotrf_t *dpotrf_f

ctypedef void dpotri_t(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dpotri_t *dpotri_f

ctypedef void dpotrs_t(char *uplo, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dpotrs_t *dpotrs_f

ctypedef void dppcon_t(char *uplo, int *n, d *ap, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dppcon_t *dppcon_f

ctypedef void dppequ_t(char *uplo, int *n, d *ap, d *s, d *scond, d *amax, int *info) nogil
cdef dppequ_t *dppequ_f

ctypedef void dpprfs_t(char *uplo, int *n, int *nrhs, d *ap, d *afp, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dpprfs_t *dpprfs_f

ctypedef void dppsv_t(char *uplo, int *n, int *nrhs, d *ap, d *b, int *ldb, int *info) nogil
cdef dppsv_t *dppsv_f

ctypedef void dppsvx_t(char *fact, char *uplo, int *n, int *nrhs, d *ap, d *afp, char *equed, d *s, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dppsvx_t *dppsvx_f

ctypedef void dpptrf_t(char *uplo, int *n, d *ap, int *info) nogil
cdef dpptrf_t *dpptrf_f

ctypedef void dpptri_t(char *uplo, int *n, d *ap, int *info) nogil
cdef dpptri_t *dpptri_f

ctypedef void dpptrs_t(char *uplo, int *n, int *nrhs, d *ap, d *b, int *ldb, int *info) nogil
cdef dpptrs_t *dpptrs_f

ctypedef void dptcon_t(int *n, d *d, d *e, d *anorm, d *rcond, d *work, int *info) nogil
cdef dptcon_t *dptcon_f

ctypedef void dpteqr_t(char *compz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *info) nogil
cdef dpteqr_t *dpteqr_f

ctypedef void dptrfs_t(int *n, int *nrhs, d *d, d *e, d *df, d *ef, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *info) nogil
cdef dptrfs_t *dptrfs_f

ctypedef void dptsv_t(int *n, int *nrhs, d *d, d *e, d *b, int *ldb, int *info) nogil
cdef dptsv_t *dptsv_f

ctypedef void dptsvx_t(char *fact, int *n, int *nrhs, d *d, d *e, d *df, d *ef, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *info) nogil
cdef dptsvx_t *dptsvx_f

ctypedef void dpttrf_t(int *n, d *d, d *e, int *info) nogil
cdef dpttrf_t *dpttrf_f

ctypedef void dpttrs_t(int *n, int *nrhs, d *d, d *e, d *b, int *ldb, int *info) nogil
cdef dpttrs_t *dpttrs_f

ctypedef void dptts2_t(int *n, int *nrhs, d *d, d *e, d *b, int *ldb) nogil
cdef dptts2_t *dptts2_f

ctypedef void drscl_t(int *n, d *sa, d *sx, int *incx) nogil
cdef drscl_t *drscl_f

ctypedef void dsbev_t(char *jobz, char *uplo, int *n, int *kd, d *ab, int *ldab, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dsbev_t *dsbev_f

ctypedef void dsbevd_t(char *jobz, char *uplo, int *n, int *kd, d *ab, int *ldab, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsbevd_t *dsbevd_f

ctypedef void dsbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, d *ab, int *ldab, d *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dsbevx_t *dsbevx_f

ctypedef void dsbgst_t(char *vect, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *x, int *ldx, d *work, int *info) nogil
cdef dsbgst_t *dsbgst_f

ctypedef void dsbgv_t(char *jobz, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dsbgv_t *dsbgv_f

ctypedef void dsbgvd_t(char *jobz, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsbgvd_t *dsbgvd_f

ctypedef void dsbgvx_t(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dsbgvx_t *dsbgvx_f

ctypedef void dsbtrd_t(char *vect, char *uplo, int *n, int *kd, d *ab, int *ldab, d *d, d *e, d *q, int *ldq, d *work, int *info) nogil
cdef dsbtrd_t *dsbtrd_f

ctypedef void dsgesv_t(int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *work, s *swork, int *iter, int *info) nogil
cdef dsgesv_t *dsgesv_f

ctypedef void dspcon_t(char *uplo, int *n, d *ap, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dspcon_t *dspcon_f

ctypedef void dspev_t(char *jobz, char *uplo, int *n, d *ap, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dspev_t *dspev_f

ctypedef void dspevd_t(char *jobz, char *uplo, int *n, d *ap, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dspevd_t *dspevd_f

ctypedef void dspevx_t(char *jobz, char *range, char *uplo, int *n, d *ap, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dspevx_t *dspevx_f

ctypedef void dspgst_t(int *itype, char *uplo, int *n, d *ap, d *bp, int *info) nogil
cdef dspgst_t *dspgst_f

ctypedef void dspgv_t(int *itype, char *jobz, char *uplo, int *n, d *ap, d *bp, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dspgv_t *dspgv_f

ctypedef void dspgvd_t(int *itype, char *jobz, char *uplo, int *n, d *ap, d *bp, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dspgvd_t *dspgvd_f

ctypedef void dspgvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, d *ap, d *bp, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dspgvx_t *dspgvx_f

ctypedef void dsprfs_t(char *uplo, int *n, int *nrhs, d *ap, d *afp, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dsprfs_t *dsprfs_f

ctypedef void dspsv_t(char *uplo, int *n, int *nrhs, d *ap, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dspsv_t *dspsv_f

ctypedef void dspsvx_t(char *fact, char *uplo, int *n, int *nrhs, d *ap, d *afp, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dspsvx_t *dspsvx_f

ctypedef void dsptrd_t(char *uplo, int *n, d *ap, d *d, d *e, d *tau, int *info) nogil
cdef dsptrd_t *dsptrd_f

ctypedef void dsptrf_t(char *uplo, int *n, d *ap, int *ipiv, int *info) nogil
cdef dsptrf_t *dsptrf_f

ctypedef void dsptri_t(char *uplo, int *n, d *ap, int *ipiv, d *work, int *info) nogil
cdef dsptri_t *dsptri_f

ctypedef void dsptrs_t(char *uplo, int *n, int *nrhs, d *ap, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dsptrs_t *dsptrs_f

ctypedef void dstebz_t(char *range, char *order, int *n, d *vl, d *vu, int *il, int *iu, d *abstol, d *d, d *e, int *m, int *nsplit, d *w, int *iblock, int *isplit, d *work, int *iwork, int *info) nogil
cdef dstebz_t *dstebz_f

ctypedef void dstedc_t(char *compz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstedc_t *dstedc_f

ctypedef void dstegr_t(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstegr_t *dstegr_f

ctypedef void dstein_t(int *n, d *d, d *e, int *m, d *w, int *iblock, int *isplit, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dstein_t *dstein_f

ctypedef void dstemr_t(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, int *m, d *w, d *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstemr_t *dstemr_f

ctypedef void dsteqr_t(char *compz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *info) nogil
cdef dsteqr_t *dsteqr_f

ctypedef void dsterf_t(int *n, d *d, d *e, int *info) nogil
cdef dsterf_t *dsterf_f

ctypedef void dstev_t(char *jobz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *info) nogil
cdef dstev_t *dstev_f

ctypedef void dstevd_t(char *jobz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstevd_t *dstevd_f

ctypedef void dstevr_t(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstevr_t *dstevr_f

ctypedef void dstevx_t(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dstevx_t *dstevx_f

ctypedef void dsycon_t(char *uplo, int *n, d *a, int *lda, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dsycon_t *dsycon_f

ctypedef void dsyev_t(char *jobz, char *uplo, int *n, d *a, int *lda, d *w, d *work, int *lwork, int *info) nogil
cdef dsyev_t *dsyev_f

ctypedef void dsyevd_t(char *jobz, char *uplo, int *n, d *a, int *lda, d *w, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsyevd_t *dsyevd_f

ctypedef void dsyevr_t(char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsyevr_t *dsyevr_f

ctypedef void dsyevx_t(char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef dsyevx_t *dsyevx_f

ctypedef void dsygs2_t(int *itype, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dsygs2_t *dsygs2_f

ctypedef void dsygst_t(int *itype, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dsygst_t *dsygst_f

ctypedef void dsygv_t(int *itype, char *jobz, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *w, d *work, int *lwork, int *info) nogil
cdef dsygv_t *dsygv_f

ctypedef void dsygvd_t(int *itype, char *jobz, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *w, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsygvd_t *dsygvd_f

ctypedef void dsygvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef dsygvx_t *dsygvx_f

ctypedef void dsyrfs_t(char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dsyrfs_t *dsyrfs_f

ctypedef void dsysv_t(char *uplo, int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, d *work, int *lwork, int *info) nogil
cdef dsysv_t *dsysv_f

ctypedef void dsysvx_t(char *fact, char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *lwork, int *iwork, int *info) nogil
cdef dsysvx_t *dsysvx_f

ctypedef void dsytd2_t(char *uplo, int *n, d *a, int *lda, d *d, d *e, d *tau, int *info) nogil
cdef dsytd2_t *dsytd2_f

ctypedef void dsytf2_t(char *uplo, int *n, d *a, int *lda, int *ipiv, int *info) nogil
cdef dsytf2_t *dsytf2_f

ctypedef void dsytrd_t(char *uplo, int *n, d *a, int *lda, d *d, d *e, d *tau, d *work, int *lwork, int *info) nogil
cdef dsytrd_t *dsytrd_f

ctypedef void dsytrf_t(char *uplo, int *n, d *a, int *lda, int *ipiv, d *work, int *lwork, int *info) nogil
cdef dsytrf_t *dsytrf_f

ctypedef void dsytri_t(char *uplo, int *n, d *a, int *lda, int *ipiv, d *work, int *info) nogil
cdef dsytri_t *dsytri_f

ctypedef void dsytrs_t(char *uplo, int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dsytrs_t *dsytrs_f

ctypedef void dtbcon_t(char *norm, char *uplo, char *diag, int *n, int *kd, d *ab, int *ldab, d *rcond, d *work, int *iwork, int *info) nogil
cdef dtbcon_t *dtbcon_f

ctypedef void dtbrfs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dtbrfs_t *dtbrfs_f

ctypedef void dtbtrs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
cdef dtbtrs_t *dtbtrs_f

ctypedef void dtgevc_t(char *side, char *howmny, bint *select, int *n, d *s, int *lds, d *p, int *ldp, d *vl, int *ldvl, d *vr, int *ldvr, int *mm, int *m, d *work, int *info) nogil
cdef dtgevc_t *dtgevc_f

ctypedef void dtgex2_t(bint *wantq, bint *wantz, int *n, d *a, int *lda, d *b, int *ldb, d *q, int *ldq, d *z, int *ldz, int *j1, int *n1, int *n2, d *work, int *lwork, int *info) nogil
cdef dtgex2_t *dtgex2_f

ctypedef void dtgexc_t(bint *wantq, bint *wantz, int *n, d *a, int *lda, d *b, int *ldb, d *q, int *ldq, d *z, int *ldz, int *ifst, int *ilst, d *work, int *lwork, int *info) nogil
cdef dtgexc_t *dtgexc_f

ctypedef void dtgsen_t(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *q, int *ldq, d *z, int *ldz, int *m, d *pl, d *pr, d *dif, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dtgsen_t *dtgsen_f

ctypedef void dtgsja_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, d *a, int *lda, d *b, int *ldb, d *tola, d *tolb, d *alpha, d *beta, d *u, int *ldu, d *v, int *ldv, d *q, int *ldq, d *work, int *ncycle, int *info) nogil
cdef dtgsja_t *dtgsja_f

ctypedef void dtgsna_t(char *job, char *howmny, bint *select, int *n, d *a, int *lda, d *b, int *ldb, d *vl, int *ldvl, d *vr, int *ldvr, d *s, d *dif, int *mm, int *m, d *work, int *lwork, int *iwork, int *info) nogil
cdef dtgsna_t *dtgsna_f

ctypedef void dtgsy2_t(char *trans, int *ijob, int *m, int *n, d *a, int *lda, d *b, int *ldb, d *c, int *ldc, d *d, int *ldd, d *e, int *lde, d *f, int *ldf, d *scale, d *rdsum, d *rdscal, int *iwork, int *pq, int *info) nogil
cdef dtgsy2_t *dtgsy2_f

ctypedef void dtgsyl_t(char *trans, int *ijob, int *m, int *n, d *a, int *lda, d *b, int *ldb, d *c, int *ldc, d *d, int *ldd, d *e, int *lde, d *f, int *ldf, d *scale, d *dif, d *work, int *lwork, int *iwork, int *info) nogil
cdef dtgsyl_t *dtgsyl_f

ctypedef void dtpcon_t(char *norm, char *uplo, char *diag, int *n, d *ap, d *rcond, d *work, int *iwork, int *info) nogil
cdef dtpcon_t *dtpcon_f

ctypedef void dtprfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *ap, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dtprfs_t *dtprfs_f

ctypedef void dtptri_t(char *uplo, char *diag, int *n, d *ap, int *info) nogil
cdef dtptri_t *dtptri_f

ctypedef void dtptrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *ap, d *b, int *ldb, int *info) nogil
cdef dtptrs_t *dtptrs_f

ctypedef void dtrcon_t(char *norm, char *uplo, char *diag, int *n, d *a, int *lda, d *rcond, d *work, int *iwork, int *info) nogil
cdef dtrcon_t *dtrcon_f

ctypedef void dtrevc_t(char *side, char *howmny, bint *select, int *n, d *t, int *ldt, d *vl, int *ldvl, d *vr, int *ldvr, int *mm, int *m, d *work, int *info) nogil
cdef dtrevc_t *dtrevc_f

ctypedef void dtrexc_t(char *compq, int *n, d *t, int *ldt, d *q, int *ldq, int *ifst, int *ilst, d *work, int *info) nogil
cdef dtrexc_t *dtrexc_f

ctypedef void dtrrfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dtrrfs_t *dtrrfs_f

ctypedef void dtrsen_t(char *job, char *compq, bint *select, int *n, d *t, int *ldt, d *q, int *ldq, d *wr, d *wi, int *m, d *s, d *sep, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dtrsen_t *dtrsen_f

ctypedef void dtrsna_t(char *job, char *howmny, bint *select, int *n, d *t, int *ldt, d *vl, int *ldvl, d *vr, int *ldvr, d *s, d *sep, int *mm, int *m, d *work, int *ldwork, int *iwork, int *info) nogil
cdef dtrsna_t *dtrsna_f

ctypedef void dtrsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, d *a, int *lda, d *b, int *ldb, d *c, int *ldc, d *scale, int *info) nogil
cdef dtrsyl_t *dtrsyl_f

ctypedef void dtrti2_t(char *uplo, char *diag, int *n, d *a, int *lda, int *info) nogil
cdef dtrti2_t *dtrti2_f

ctypedef void dtrtri_t(char *uplo, char *diag, int *n, d *a, int *lda, int *info) nogil
cdef dtrtri_t *dtrtri_f

ctypedef void dtrtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dtrtrs_t *dtrtrs_f

ctypedef void dtzrqf_t(int *m, int *n, d *a, int *lda, d *tau, int *info) nogil
cdef dtzrqf_t *dtzrqf_f

ctypedef void dtzrzf_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dtzrzf_t *dtzrzf_f

ctypedef d dzsum1_t(int *n, z *cx, int *incx) nogil
cdef dzsum1_t *dzsum1_f

ctypedef int icmax1_t(int *n, c *cx, int *incx) nogil
cdef icmax1_t *icmax1_f

ctypedef int ieeeck_t(int *ispec, s *zero, s *one) nogil
cdef ieeeck_t *ieeeck_f

ctypedef int ilaenv_t(int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4) nogil
cdef ilaenv_t *ilaenv_f

ctypedef int iparmq_t(int *ispec, char *name, char *opts, int *n, int *ilo, int *ihi, int *lwork) nogil
cdef iparmq_t *iparmq_f

ctypedef int izmax1_t(int *n, z *cx, int *incx) nogil
cdef izmax1_t *izmax1_f

ctypedef bint lsamen_t(int *n, char *ca, char *cb) nogil
cdef lsamen_t *lsamen_f

ctypedef void sbdsdc_t(char *uplo, char *compq, int *n, s *d, s *e, s *u, int *ldu, s *vt, int *ldvt, s *q, int *iq, s *work, int *iwork, int *info) nogil
cdef sbdsdc_t *sbdsdc_f

ctypedef void sbdsqr_t(char *uplo, int *n, int *ncvt, int *nru, int *ncc, s *d, s *e, s *vt, int *ldvt, s *u, int *ldu, s *c, int *ldc, s *work, int *info) nogil
cdef sbdsqr_t *sbdsqr_f

ctypedef s scsum1_t(int *n, c *cx, int *incx) nogil
cdef scsum1_t *scsum1_f

ctypedef void sdisna_t(char *job, int *m, int *n, s *d, s *sep, int *info) nogil
cdef sdisna_t *sdisna_f

ctypedef void sgbbrd_t(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, s *ab, int *ldab, s *d, s *e, s *q, int *ldq, s *pt, int *ldpt, s *c, int *ldc, s *work, int *info) nogil
cdef sgbbrd_t *sgbbrd_f

ctypedef void sgbcon_t(char *norm, int *n, int *kl, int *ku, s *ab, int *ldab, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sgbcon_t *sgbcon_f

ctypedef void sgbequ_t(int *m, int *n, int *kl, int *ku, s *ab, int *ldab, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef sgbequ_t *sgbequ_f

ctypedef void sgbrfs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgbrfs_t *sgbrfs_f

ctypedef void sgbsv_t(int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgbsv_t *sgbsv_f

ctypedef void sgbsvx_t(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, int *ipiv, char *equed, s *r, s *c, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgbsvx_t *sgbsvx_f

ctypedef void sgbtf2_t(int *m, int *n, int *kl, int *ku, s *ab, int *ldab, int *ipiv, int *info) nogil
cdef sgbtf2_t *sgbtf2_f

ctypedef void sgbtrf_t(int *m, int *n, int *kl, int *ku, s *ab, int *ldab, int *ipiv, int *info) nogil
cdef sgbtrf_t *sgbtrf_f

ctypedef void sgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgbtrs_t *sgbtrs_f

ctypedef void sgebak_t(char *job, char *side, int *n, int *ilo, int *ihi, s *scale, int *m, s *v, int *ldv, int *info) nogil
cdef sgebak_t *sgebak_f

ctypedef void sgebal_t(char *job, int *n, s *a, int *lda, int *ilo, int *ihi, s *scale, int *info) nogil
cdef sgebal_t *sgebal_f

ctypedef void sgebd2_t(int *m, int *n, s *a, int *lda, s *d, s *e, s *tauq, s *taup, s *work, int *info) nogil
cdef sgebd2_t *sgebd2_f

ctypedef void sgebrd_t(int *m, int *n, s *a, int *lda, s *d, s *e, s *tauq, s *taup, s *work, int *lwork, int *info) nogil
cdef sgebrd_t *sgebrd_f

ctypedef void sgecon_t(char *norm, int *n, s *a, int *lda, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sgecon_t *sgecon_f

ctypedef void sgeequ_t(int *m, int *n, s *a, int *lda, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef sgeequ_t *sgeequ_f

ctypedef void sgees_t(char *jobvs, char *sort, sselect2 *select, int *n, s *a, int *lda, int *sdim, s *wr, s *wi, s *vs, int *ldvs, s *work, int *lwork, bint *bwork, int *info) nogil
cdef sgees_t *sgees_f

ctypedef void sgeesx_t(char *jobvs, char *sort, sselect2 *select, char *sense, int *n, s *a, int *lda, int *sdim, s *wr, s *wi, s *vs, int *ldvs, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef sgeesx_t *sgeesx_f

ctypedef void sgeev_t(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *wr, s *wi, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
cdef sgeev_t *sgeev_f

ctypedef void sgeevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, s *a, int *lda, s *wr, s *wi, s *vl, int *ldvl, s *vr, int *ldvr, int *ilo, int *ihi, s *scale, s *abnrm, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, int *info) nogil
cdef sgeevx_t *sgeevx_f

ctypedef void sgegs_t(char *jobvsl, char *jobvsr, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vsl, int *ldvsl, s *vsr, int *ldvsr, s *work, int *lwork, int *info) nogil
cdef sgegs_t *sgegs_f

ctypedef void sgegv_t(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
cdef sgegv_t *sgegv_f

ctypedef void sgehd2_t(int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgehd2_t *sgehd2_f

ctypedef void sgehrd_t(int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgehrd_t *sgehrd_f

ctypedef void sgelq2_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgelq2_t *sgelq2_f

ctypedef void sgelqf_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgelqf_t *sgelqf_f

ctypedef void sgels_t(char *trans, int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *work, int *lwork, int *info) nogil
cdef sgels_t *sgels_f

ctypedef void sgelsd_t(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *s, s *rcond, int *rank, s *work, int *lwork, int *iwork, int *info) nogil
cdef sgelsd_t *sgelsd_f

ctypedef void sgelss_t(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *s, s *rcond, int *rank, s *work, int *lwork, int *info) nogil
cdef sgelss_t *sgelss_f

ctypedef void sgelsx_t(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *jpvt, s *rcond, int *rank, s *work, int *info) nogil
cdef sgelsx_t *sgelsx_f

ctypedef void sgelsy_t(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *jpvt, s *rcond, int *rank, s *work, int *lwork, int *info) nogil
cdef sgelsy_t *sgelsy_f

ctypedef void sgeql2_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgeql2_t *sgeql2_f

ctypedef void sgeqlf_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgeqlf_t *sgeqlf_f

ctypedef void sgeqp3_t(int *m, int *n, s *a, int *lda, int *jpvt, s *tau, s *work, int *lwork, int *info) nogil
cdef sgeqp3_t *sgeqp3_f

ctypedef void sgeqpf_t(int *m, int *n, s *a, int *lda, int *jpvt, s *tau, s *work, int *info) nogil
cdef sgeqpf_t *sgeqpf_f

ctypedef void sgeqr2_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgeqr2_t *sgeqr2_f

ctypedef void sgeqrf_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgeqrf_t *sgeqrf_f

ctypedef void sgerfs_t(char *trans, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgerfs_t *sgerfs_f

ctypedef void sgerq2_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgerq2_t *sgerq2_f

ctypedef void sgerqf_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgerqf_t *sgerqf_f

ctypedef void sgesc2_t(int *n, s *a, int *lda, s *rhs, int *ipiv, int *jpiv, s *scale) nogil
cdef sgesc2_t *sgesc2_f

ctypedef void sgesdd_t(char *jobz, int *m, int *n, s *a, int *lda, s *s, s *u, int *ldu, s *vt, int *ldvt, s *work, int *lwork, int *iwork, int *info) nogil
cdef sgesdd_t *sgesdd_f

ctypedef void sgesv_t(int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgesv_t *sgesv_f

ctypedef void sgesvd_t(char *jobu, char *jobvt, int *m, int *n, s *a, int *lda, s *s, s *u, int *ldu, s *vt, int *ldvt, s *work, int *lwork, int *info) nogil
cdef sgesvd_t *sgesvd_f

ctypedef void sgesvx_t(char *fact, char *trans, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, char *equed, s *r, s *c, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgesvx_t *sgesvx_f

ctypedef void sgetc2_t(int *n, s *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef sgetc2_t *sgetc2_f

ctypedef void sgetf2_t(int *m, int *n, s *a, int *lda, int *ipiv, int *info) nogil
cdef sgetf2_t *sgetf2_f

ctypedef void sgetrf_t(int *m, int *n, s *a, int *lda, int *ipiv, int *info) nogil
cdef sgetrf_t *sgetrf_f

ctypedef void sgetri_t(int *n, s *a, int *lda, int *ipiv, s *work, int *lwork, int *info) nogil
cdef sgetri_t *sgetri_f

ctypedef void sgetrs_t(char *trans, int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgetrs_t *sgetrs_f

ctypedef void sggbak_t(char *job, char *side, int *n, int *ilo, int *ihi, s *lscale, s *rscale, int *m, s *v, int *ldv, int *info) nogil
cdef sggbak_t *sggbak_f

ctypedef void sggbal_t(char *job, int *n, s *a, int *lda, s *b, int *ldb, int *ilo, int *ihi, s *lscale, s *rscale, s *work, int *info) nogil
cdef sggbal_t *sggbal_f

ctypedef void sgges_t(char *jobvsl, char *jobvsr, char *sort, sselect3 *selctg, int *n, s *a, int *lda, s *b, int *ldb, int *sdim, s *alphar, s *alphai, s *beta, s *vsl, int *ldvsl, s *vsr, int *ldvsr, s *work, int *lwork, bint *bwork, int *info) nogil
cdef sgges_t *sgges_f

ctypedef void sggesx_t(char *jobvsl, char *jobvsr, char *sort, sselect3 *selctg, char *sense, int *n, s *a, int *lda, s *b, int *ldb, int *sdim, s *alphar, s *alphai, s *beta, s *vsl, int *ldvsl, s *vsr, int *ldvsr, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef sggesx_t *sggesx_f

ctypedef void sggev_t(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
cdef sggev_t *sggev_f

ctypedef void sggevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, int *ilo, int *ihi, s *lscale, s *rscale, s *abnrm, s *bbnrm, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, bint *bwork, int *info) nogil
cdef sggevx_t *sggevx_f

ctypedef void sggglm_t(int *n, int *m, int *p, s *a, int *lda, s *b, int *ldb, s *d, s *x, s *y, s *work, int *lwork, int *info) nogil
cdef sggglm_t *sggglm_f

ctypedef void sgghrd_t(char *compq, char *compz, int *n, int *ilo, int *ihi, s *a, int *lda, s *b, int *ldb, s *q, int *ldq, s *z, int *ldz, int *info) nogil
cdef sgghrd_t *sgghrd_f

ctypedef void sgglse_t(int *m, int *n, int *p, s *a, int *lda, s *b, int *ldb, s *c, s *d, s *x, s *work, int *lwork, int *info) nogil
cdef sgglse_t *sgglse_f

ctypedef void sggqrf_t(int *n, int *m, int *p, s *a, int *lda, s *taua, s *b, int *ldb, s *taub, s *work, int *lwork, int *info) nogil
cdef sggqrf_t *sggqrf_f

ctypedef void sggrqf_t(int *m, int *p, int *n, s *a, int *lda, s *taua, s *b, int *ldb, s *taub, s *work, int *lwork, int *info) nogil
cdef sggrqf_t *sggrqf_f

ctypedef void sggsvd_t(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, s *a, int *lda, s *b, int *ldb, s *alpha, s *beta, s *u, int *ldu, s *v, int *ldv, s *q, int *ldq, s *work, int *iwork, int *info) nogil
cdef sggsvd_t *sggsvd_f

ctypedef void sggsvp_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, s *a, int *lda, s *b, int *ldb, s *tola, s *tolb, int *k, int *l, s *u, int *ldu, s *v, int *ldv, s *q, int *ldq, int *iwork, s *tau, s *work, int *info) nogil
cdef sggsvp_t *sggsvp_f

ctypedef void sgtcon_t(char *norm, int *n, s *dl, s *d, s *du, s *du2, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sgtcon_t *sgtcon_f

ctypedef void sgtrfs_t(char *trans, int *n, int *nrhs, s *dl, s *d, s *du, s *dlf, s *df, s *duf, s *du2, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgtrfs_t *sgtrfs_f

ctypedef void sgtsv_t(int *n, int *nrhs, s *dl, s *d, s *du, s *b, int *ldb, int *info) nogil
cdef sgtsv_t *sgtsv_f

ctypedef void sgtsvx_t(char *fact, char *trans, int *n, int *nrhs, s *dl, s *d, s *du, s *dlf, s *df, s *duf, s *du2, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgtsvx_t *sgtsvx_f

ctypedef void sgttrf_t(int *n, s *dl, s *d, s *du, s *du2, int *ipiv, int *info) nogil
cdef sgttrf_t *sgttrf_f

ctypedef void sgttrs_t(char *trans, int *n, int *nrhs, s *dl, s *d, s *du, s *du2, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgttrs_t *sgttrs_f

ctypedef void sgtts2_t(int *itrans, int *n, int *nrhs, s *dl, s *d, s *du, s *du2, int *ipiv, s *b, int *ldb) nogil
cdef sgtts2_t *sgtts2_f

ctypedef void shgeqz_t(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, s *h, int *ldh, s *t, int *ldt, s *alphar, s *alphai, s *beta, s *q, int *ldq, s *z, int *ldz, s *work, int *lwork, int *info) nogil
cdef shgeqz_t *shgeqz_f

ctypedef void shsein_t(char *side, char *eigsrc, char *initv, bint *select, int *n, s *h, int *ldh, s *wr, s *wi, s *vl, int *ldvl, s *vr, int *ldvr, int *mm, int *m, s *work, int *ifaill, int *ifailr, int *info) nogil
cdef shsein_t *shsein_f

ctypedef void shseqr_t(char *job, char *compz, int *n, int *ilo, int *ihi, s *h, int *ldh, s *wr, s *wi, s *z, int *ldz, s *work, int *lwork, int *info) nogil
cdef shseqr_t *shseqr_f

ctypedef void slacn2_t(int *n, s *v, s *x, int *isgn, s *est, int *kase, int *isave) nogil
cdef slacn2_t *slacn2_f

ctypedef void slacon_t(int *n, s *v, s *x, int *isgn, s *est, int *kase) nogil
cdef slacon_t *slacon_f

ctypedef s slamch_t(char *cmach) nogil
cdef slamch_t *slamch_f

ctypedef s slangb_t(char *norm, int *n, int *kl, int *ku, s *ab, int *ldab, s *work) nogil
cdef slangb_t *slangb_f

ctypedef s slange_t(char *norm, int *m, int *n, s *a, int *lda, s *work) nogil
cdef slange_t *slange_f

ctypedef s slangt_t(char *norm, int *n, s *dl, s *d, s *du) nogil
cdef slangt_t *slangt_f

ctypedef s slanhs_t(char *norm, int *n, s *a, int *lda, s *work) nogil
cdef slanhs_t *slanhs_f

ctypedef s slansb_t(char *norm, char *uplo, int *n, int *k, s *ab, int *ldab, s *work) nogil
cdef slansb_t *slansb_f

ctypedef s slansp_t(char *norm, char *uplo, int *n, s *ap, s *work) nogil
cdef slansp_t *slansp_f

ctypedef s slanst_t(char *norm, int *n, s *d, s *e) nogil
cdef slanst_t *slanst_f

ctypedef s slansy_t(char *norm, char *uplo, int *n, s *a, int *lda, s *work) nogil
cdef slansy_t *slansy_f

ctypedef s slantb_t(char *norm, char *uplo, char *diag, int *n, int *k, s *ab, int *ldab, s *work) nogil
cdef slantb_t *slantb_f

ctypedef s slantp_t(char *norm, char *uplo, char *diag, int *n, s *ap, s *work) nogil
cdef slantp_t *slantp_f

ctypedef s slantr_t(char *norm, char *uplo, char *diag, int *m, int *n, s *a, int *lda, s *work) nogil
cdef slantr_t *slantr_f

ctypedef void slanv2_t(s *a, s *b, s *c, s *d, s *rt1r, s *rt1i, s *rt2r, s *rt2i, s *cs, s *sn) nogil
cdef slanv2_t *slanv2_f

ctypedef void slarf_t(char *side, int *m, int *n, s *v, int *incv, s *tau, s *c, int *ldc, s *work) nogil
cdef slarf_t *slarf_f

ctypedef void slarz_t(char *side, int *m, int *n, int *l, s *v, int *incv, s *tau, s *c, int *ldc, s *work) nogil
cdef slarz_t *slarz_f

ctypedef void slaswp_t(int *n, s *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef slaswp_t *slaswp_f

ctypedef void slauum_t(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef slauum_t *slauum_f

ctypedef void sopgtr_t(char *uplo, int *n, s *ap, s *tau, s *q, int *ldq, s *work, int *info) nogil
cdef sopgtr_t *sopgtr_f

ctypedef void sopmtr_t(char *side, char *uplo, char *trans, int *m, int *n, s *ap, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sopmtr_t *sopmtr_f

ctypedef void sorg2l_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorg2l_t *sorg2l_f

ctypedef void sorg2r_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorg2r_t *sorg2r_f

ctypedef void sorgbr_t(char *vect, int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgbr_t *sorgbr_f

ctypedef void sorghr_t(int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorghr_t *sorghr_f

ctypedef void sorgl2_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorgl2_t *sorgl2_f

ctypedef void sorglq_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorglq_t *sorglq_f

ctypedef void sorgql_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgql_t *sorgql_f

ctypedef void sorgqr_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgqr_t *sorgqr_f

ctypedef void sorgr2_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorgr2_t *sorgr2_f

ctypedef void sorgrq_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgrq_t *sorgrq_f

ctypedef void sorgtr_t(char *uplo, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgtr_t *sorgtr_f

ctypedef void sorm2l_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sorm2l_t *sorm2l_f

ctypedef void sorm2r_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sorm2r_t *sorm2r_f

ctypedef void sormbr_t(char *vect, char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormbr_t *sormbr_f

ctypedef void sormhr_t(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormhr_t *sormhr_f

ctypedef void sorml2_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sorml2_t *sorml2_f

ctypedef void sormlq_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormlq_t *sormlq_f

ctypedef void sormql_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormql_t *sormql_f

ctypedef void sormqr_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormqr_t *sormqr_f

ctypedef void sormr2_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sormr2_t *sormr2_f

ctypedef void sormr3_t(char *side, char *trans, int *m, int *n, int *k, int *l, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sormr3_t *sormr3_f

ctypedef void sormrq_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormrq_t *sormrq_f

ctypedef void sormrz_t(char *side, char *trans, int *m, int *n, int *k, int *l, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormrz_t *sormrz_f

ctypedef void sormtr_t(char *side, char *uplo, char *trans, int *m, int *n, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormtr_t *sormtr_f

ctypedef void spbcon_t(char *uplo, int *n, int *kd, s *ab, int *ldab, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef spbcon_t *spbcon_f

ctypedef void spbequ_t(char *uplo, int *n, int *kd, s *ab, int *ldab, s *s, s *scond, s *amax, int *info) nogil
cdef spbequ_t *spbequ_f

ctypedef void spbrfs_t(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef spbrfs_t *spbrfs_f

ctypedef void spbstf_t(char *uplo, int *n, int *kd, s *ab, int *ldab, int *info) nogil
cdef spbstf_t *spbstf_f

ctypedef void spbsv_t(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
cdef spbsv_t *spbsv_f

ctypedef void spbsvx_t(char *fact, char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, char *equed, s *s, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef spbsvx_t *spbsvx_f

ctypedef void spbtf2_t(char *uplo, int *n, int *kd, s *ab, int *ldab, int *info) nogil
cdef spbtf2_t *spbtf2_f

ctypedef void spbtrf_t(char *uplo, int *n, int *kd, s *ab, int *ldab, int *info) nogil
cdef spbtrf_t *spbtrf_f

ctypedef void spbtrs_t(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
cdef spbtrs_t *spbtrs_f

ctypedef void spocon_t(char *uplo, int *n, s *a, int *lda, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef spocon_t *spocon_f

ctypedef void spoequ_t(int *n, s *a, int *lda, s *s, s *scond, s *amax, int *info) nogil
cdef spoequ_t *spoequ_f

ctypedef void sporfs_t(char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sporfs_t *sporfs_f

ctypedef void sposv_t(char *uplo, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef sposv_t *sposv_f

ctypedef void sposvx_t(char *fact, char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, char *equed, s *s, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sposvx_t *sposvx_f

ctypedef void spotf2_t(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef spotf2_t *spotf2_f

ctypedef void spotrf_t(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef spotrf_t *spotrf_f

ctypedef void spotri_t(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef spotri_t *spotri_f

ctypedef void spotrs_t(char *uplo, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef spotrs_t *spotrs_f

ctypedef void sppcon_t(char *uplo, int *n, s *ap, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sppcon_t *sppcon_f

ctypedef void sppequ_t(char *uplo, int *n, s *ap, s *s, s *scond, s *amax, int *info) nogil
cdef sppequ_t *sppequ_f

ctypedef void spprfs_t(char *uplo, int *n, int *nrhs, s *ap, s *afp, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef spprfs_t *spprfs_f

ctypedef void sppsv_t(char *uplo, int *n, int *nrhs, s *ap, s *b, int *ldb, int *info) nogil
cdef sppsv_t *sppsv_f

ctypedef void sppsvx_t(char *fact, char *uplo, int *n, int *nrhs, s *ap, s *afp, char *equed, s *s, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sppsvx_t *sppsvx_f

ctypedef void spptrf_t(char *uplo, int *n, s *ap, int *info) nogil
cdef spptrf_t *spptrf_f

ctypedef void spptri_t(char *uplo, int *n, s *ap, int *info) nogil
cdef spptri_t *spptri_f

ctypedef void spptrs_t(char *uplo, int *n, int *nrhs, s *ap, s *b, int *ldb, int *info) nogil
cdef spptrs_t *spptrs_f

ctypedef void sptcon_t(int *n, s *d, s *e, s *anorm, s *rcond, s *work, int *info) nogil
cdef sptcon_t *sptcon_f

ctypedef void spteqr_t(char *compz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *info) nogil
cdef spteqr_t *spteqr_f

ctypedef void sptrfs_t(int *n, int *nrhs, s *d, s *e, s *df, s *ef, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *info) nogil
cdef sptrfs_t *sptrfs_f

ctypedef void sptsv_t(int *n, int *nrhs, s *d, s *e, s *b, int *ldb, int *info) nogil
cdef sptsv_t *sptsv_f

ctypedef void sptsvx_t(char *fact, int *n, int *nrhs, s *d, s *e, s *df, s *ef, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *info) nogil
cdef sptsvx_t *sptsvx_f

ctypedef void spttrf_t(int *n, s *d, s *e, int *info) nogil
cdef spttrf_t *spttrf_f

ctypedef void spttrs_t(int *n, int *nrhs, s *d, s *e, s *b, int *ldb, int *info) nogil
cdef spttrs_t *spttrs_f

ctypedef void sptts2_t(int *n, int *nrhs, s *d, s *e, s *b, int *ldb) nogil
cdef sptts2_t *sptts2_f

ctypedef void srscl_t(int *n, s *sa, s *sx, int *incx) nogil
cdef srscl_t *srscl_f

ctypedef void ssbev_t(char *jobz, char *uplo, int *n, int *kd, s *ab, int *ldab, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef ssbev_t *ssbev_f

ctypedef void ssbevd_t(char *jobz, char *uplo, int *n, int *kd, s *ab, int *ldab, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssbevd_t *ssbevd_f

ctypedef void ssbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, s *ab, int *ldab, s *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef ssbevx_t *ssbevx_f

ctypedef void ssbgst_t(char *vect, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *x, int *ldx, s *work, int *info) nogil
cdef ssbgst_t *ssbgst_f

ctypedef void ssbgv_t(char *jobz, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef ssbgv_t *ssbgv_f

ctypedef void ssbgvd_t(char *jobz, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssbgvd_t *ssbgvd_f

ctypedef void ssbgvx_t(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef ssbgvx_t *ssbgvx_f

ctypedef void ssbtrd_t(char *vect, char *uplo, int *n, int *kd, s *ab, int *ldab, s *d, s *e, s *q, int *ldq, s *work, int *info) nogil
cdef ssbtrd_t *ssbtrd_f

ctypedef void sspcon_t(char *uplo, int *n, s *ap, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sspcon_t *sspcon_f

ctypedef void sspev_t(char *jobz, char *uplo, int *n, s *ap, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef sspev_t *sspev_f

ctypedef void sspevd_t(char *jobz, char *uplo, int *n, s *ap, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sspevd_t *sspevd_f

ctypedef void sspevx_t(char *jobz, char *range, char *uplo, int *n, s *ap, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sspevx_t *sspevx_f

ctypedef void sspgst_t(int *itype, char *uplo, int *n, s *ap, s *bp, int *info) nogil
cdef sspgst_t *sspgst_f

ctypedef void sspgv_t(int *itype, char *jobz, char *uplo, int *n, s *ap, s *bp, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef sspgv_t *sspgv_f

ctypedef void sspgvd_t(int *itype, char *jobz, char *uplo, int *n, s *ap, s *bp, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sspgvd_t *sspgvd_f

ctypedef void sspgvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, s *ap, s *bp, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sspgvx_t *sspgvx_f

ctypedef void ssprfs_t(char *uplo, int *n, int *nrhs, s *ap, s *afp, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef ssprfs_t *ssprfs_f

ctypedef void sspsv_t(char *uplo, int *n, int *nrhs, s *ap, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sspsv_t *sspsv_f

ctypedef void sspsvx_t(char *fact, char *uplo, int *n, int *nrhs, s *ap, s *afp, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sspsvx_t *sspsvx_f

ctypedef void ssptrd_t(char *uplo, int *n, s *ap, s *d, s *e, s *tau, int *info) nogil
cdef ssptrd_t *ssptrd_f

ctypedef void ssptrf_t(char *uplo, int *n, s *ap, int *ipiv, int *info) nogil
cdef ssptrf_t *ssptrf_f

ctypedef void ssptri_t(char *uplo, int *n, s *ap, int *ipiv, s *work, int *info) nogil
cdef ssptri_t *ssptri_f

ctypedef void ssptrs_t(char *uplo, int *n, int *nrhs, s *ap, int *ipiv, s *b, int *ldb, int *info) nogil
cdef ssptrs_t *ssptrs_f

ctypedef void sstebz_t(char *range, char *order, int *n, s *vl, s *vu, int *il, int *iu, s *abstol, s *d, s *e, int *m, int *nsplit, s *w, int *iblock, int *isplit, s *work, int *iwork, int *info) nogil
cdef sstebz_t *sstebz_f

ctypedef void sstedc_t(char *compz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstedc_t *sstedc_f

ctypedef void sstegr_t(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstegr_t *sstegr_f

ctypedef void sstein_t(int *n, s *d, s *e, int *m, s *w, int *iblock, int *isplit, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sstein_t *sstein_f

ctypedef void sstemr_t(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, int *m, s *w, s *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstemr_t *sstemr_f

ctypedef void ssteqr_t(char *compz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *info) nogil
cdef ssteqr_t *ssteqr_f

ctypedef void ssterf_t(int *n, s *d, s *e, int *info) nogil
cdef ssterf_t *ssterf_f

ctypedef void sstev_t(char *jobz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *info) nogil
cdef sstev_t *sstev_f

ctypedef void sstevd_t(char *jobz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstevd_t *sstevd_f

ctypedef void sstevr_t(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstevr_t *sstevr_f

ctypedef void sstevx_t(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sstevx_t *sstevx_f

ctypedef void ssycon_t(char *uplo, int *n, s *a, int *lda, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef ssycon_t *ssycon_f

ctypedef void ssyev_t(char *jobz, char *uplo, int *n, s *a, int *lda, s *w, s *work, int *lwork, int *info) nogil
cdef ssyev_t *ssyev_f

ctypedef void ssyevd_t(char *jobz, char *uplo, int *n, s *a, int *lda, s *w, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssyevd_t *ssyevd_f

ctypedef void ssyevr_t(char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssyevr_t *ssyevr_f

ctypedef void ssyevx_t(char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef ssyevx_t *ssyevx_f

ctypedef void ssygs2_t(int *itype, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef ssygs2_t *ssygs2_f

ctypedef void ssygst_t(int *itype, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef ssygst_t *ssygst_f

ctypedef void ssygv_t(int *itype, char *jobz, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *w, s *work, int *lwork, int *info) nogil
cdef ssygv_t *ssygv_f

ctypedef void ssygvd_t(int *itype, char *jobz, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *w, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssygvd_t *ssygvd_f

ctypedef void ssygvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef ssygvx_t *ssygvx_f

ctypedef void ssyrfs_t(char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef ssyrfs_t *ssyrfs_f

ctypedef void ssysv_t(char *uplo, int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, s *work, int *lwork, int *info) nogil
cdef ssysv_t *ssysv_f

ctypedef void ssysvx_t(char *fact, char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *lwork, int *iwork, int *info) nogil
cdef ssysvx_t *ssysvx_f

ctypedef void ssytd2_t(char *uplo, int *n, s *a, int *lda, s *d, s *e, s *tau, int *info) nogil
cdef ssytd2_t *ssytd2_f

ctypedef void ssytf2_t(char *uplo, int *n, s *a, int *lda, int *ipiv, int *info) nogil
cdef ssytf2_t *ssytf2_f

ctypedef void ssytrd_t(char *uplo, int *n, s *a, int *lda, s *d, s *e, s *tau, s *work, int *lwork, int *info) nogil
cdef ssytrd_t *ssytrd_f

ctypedef void ssytrf_t(char *uplo, int *n, s *a, int *lda, int *ipiv, s *work, int *lwork, int *info) nogil
cdef ssytrf_t *ssytrf_f

ctypedef void ssytri_t(char *uplo, int *n, s *a, int *lda, int *ipiv, s *work, int *info) nogil
cdef ssytri_t *ssytri_f

ctypedef void ssytrs_t(char *uplo, int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
cdef ssytrs_t *ssytrs_f

ctypedef void stbcon_t(char *norm, char *uplo, char *diag, int *n, int *kd, s *ab, int *ldab, s *rcond, s *work, int *iwork, int *info) nogil
cdef stbcon_t *stbcon_f

ctypedef void stbrfs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef stbrfs_t *stbrfs_f

ctypedef void stbtrs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
cdef stbtrs_t *stbtrs_f

ctypedef void stgevc_t(char *side, char *howmny, bint *select, int *n, s *s, int *lds, s *p, int *ldp, s *vl, int *ldvl, s *vr, int *ldvr, int *mm, int *m, s *work, int *info) nogil
cdef stgevc_t *stgevc_f

ctypedef void stgex2_t(bint *wantq, bint *wantz, int *n, s *a, int *lda, s *b, int *ldb, s *q, int *ldq, s *z, int *ldz, int *j1, int *n1, int *n2, s *work, int *lwork, int *info) nogil
cdef stgex2_t *stgex2_f

ctypedef void stgexc_t(bint *wantq, bint *wantz, int *n, s *a, int *lda, s *b, int *ldb, s *q, int *ldq, s *z, int *ldz, int *ifst, int *ilst, s *work, int *lwork, int *info) nogil
cdef stgexc_t *stgexc_f

ctypedef void stgsen_t(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *q, int *ldq, s *z, int *ldz, int *m, s *pl, s *pr, s *dif, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef stgsen_t *stgsen_f

ctypedef void stgsja_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, s *a, int *lda, s *b, int *ldb, s *tola, s *tolb, s *alpha, s *beta, s *u, int *ldu, s *v, int *ldv, s *q, int *ldq, s *work, int *ncycle, int *info) nogil
cdef stgsja_t *stgsja_f

ctypedef void stgsna_t(char *job, char *howmny, bint *select, int *n, s *a, int *lda, s *b, int *ldb, s *vl, int *ldvl, s *vr, int *ldvr, s *s, s *dif, int *mm, int *m, s *work, int *lwork, int *iwork, int *info) nogil
cdef stgsna_t *stgsna_f

ctypedef void stgsy2_t(char *trans, int *ijob, int *m, int *n, s *a, int *lda, s *b, int *ldb, s *c, int *ldc, s *d, int *ldd, s *e, int *lde, s *f, int *ldf, s *scale, s *rdsum, s *rdscal, int *iwork, int *pq, int *info) nogil
cdef stgsy2_t *stgsy2_f

ctypedef void stgsyl_t(char *trans, int *ijob, int *m, int *n, s *a, int *lda, s *b, int *ldb, s *c, int *ldc, s *d, int *ldd, s *e, int *lde, s *f, int *ldf, s *scale, s *dif, s *work, int *lwork, int *iwork, int *info) nogil
cdef stgsyl_t *stgsyl_f

ctypedef void stpcon_t(char *norm, char *uplo, char *diag, int *n, s *ap, s *rcond, s *work, int *iwork, int *info) nogil
cdef stpcon_t *stpcon_f

ctypedef void stprfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *ap, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef stprfs_t *stprfs_f

ctypedef void stptri_t(char *uplo, char *diag, int *n, s *ap, int *info) nogil
cdef stptri_t *stptri_f

ctypedef void stptrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *ap, s *b, int *ldb, int *info) nogil
cdef stptrs_t *stptrs_f

ctypedef void strcon_t(char *norm, char *uplo, char *diag, int *n, s *a, int *lda, s *rcond, s *work, int *iwork, int *info) nogil
cdef strcon_t *strcon_f

ctypedef void strevc_t(char *side, char *howmny, bint *select, int *n, s *t, int *ldt, s *vl, int *ldvl, s *vr, int *ldvr, int *mm, int *m, s *work, int *info) nogil
cdef strevc_t *strevc_f

ctypedef void strexc_t(char *compq, int *n, s *t, int *ldt, s *q, int *ldq, int *ifst, int *ilst, s *work, int *info) nogil
cdef strexc_t *strexc_f

ctypedef void strrfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef strrfs_t *strrfs_f

ctypedef void strsen_t(char *job, char *compq, bint *select, int *n, s *t, int *ldt, s *q, int *ldq, s *wr, s *wi, int *m, s *s, s *sep, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef strsen_t *strsen_f

ctypedef void strsna_t(char *job, char *howmny, bint *select, int *n, s *t, int *ldt, s *vl, int *ldvl, s *vr, int *ldvr, s *s, s *sep, int *mm, int *m, s *work, int *ldwork, int *iwork, int *info) nogil
cdef strsna_t *strsna_f

ctypedef void strsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, s *a, int *lda, s *b, int *ldb, s *c, int *ldc, s *scale, int *info) nogil
cdef strsyl_t *strsyl_f

ctypedef void strti2_t(char *uplo, char *diag, int *n, s *a, int *lda, int *info) nogil
cdef strti2_t *strti2_f

ctypedef void strtri_t(char *uplo, char *diag, int *n, s *a, int *lda, int *info) nogil
cdef strtri_t *strtri_f

ctypedef void strtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef strtrs_t *strtrs_f

ctypedef void stzrqf_t(int *m, int *n, s *a, int *lda, s *tau, int *info) nogil
cdef stzrqf_t *stzrqf_f

ctypedef void stzrzf_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef stzrzf_t *stzrzf_f

ctypedef void xerbla_t(char *srname, int *info) nogil
cdef xerbla_t *xerbla_f

ctypedef void zbdsqr_t(char *uplo, int *n, int *ncvt, int *nru, int *ncc, d *d, d *e, z *vt, int *ldvt, z *u, int *ldu, z *c, int *ldc, d *rwork, int *info) nogil
cdef zbdsqr_t *zbdsqr_f

ctypedef void zcgesv_t(int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, z *x, int *ldx, z *work, c *swork, d *rwork, int *iter, int *info) nogil
cdef zcgesv_t *zcgesv_f

ctypedef void zdrscl_t(int *n, d *sa, z *sx, int *incx) nogil
cdef zdrscl_t *zdrscl_f

ctypedef void zgbbrd_t(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, z *ab, int *ldab, d *d, d *e, z *q, int *ldq, z *pt, int *ldpt, z *c, int *ldc, z *work, d *rwork, int *info) nogil
cdef zgbbrd_t *zgbbrd_f

ctypedef void zgbcon_t(char *norm, int *n, int *kl, int *ku, z *ab, int *ldab, int *ipiv, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zgbcon_t *zgbcon_f

ctypedef void zgbequ_t(int *m, int *n, int *kl, int *ku, z *ab, int *ldab, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef zgbequ_t *zgbequ_f

ctypedef void zgbrfs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgbrfs_t *zgbrfs_f

ctypedef void zgbsv_t(int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgbsv_t *zgbsv_f

ctypedef void zgbsvx_t(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, int *ipiv, char *equed, d *r, d *c, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgbsvx_t *zgbsvx_f

ctypedef void zgbtf2_t(int *m, int *n, int *kl, int *ku, z *ab, int *ldab, int *ipiv, int *info) nogil
cdef zgbtf2_t *zgbtf2_f

ctypedef void zgbtrf_t(int *m, int *n, int *kl, int *ku, z *ab, int *ldab, int *ipiv, int *info) nogil
cdef zgbtrf_t *zgbtrf_f

ctypedef void zgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgbtrs_t *zgbtrs_f

ctypedef void zgebak_t(char *job, char *side, int *n, int *ilo, int *ihi, d *scale, int *m, z *v, int *ldv, int *info) nogil
cdef zgebak_t *zgebak_f

ctypedef void zgebal_t(char *job, int *n, z *a, int *lda, int *ilo, int *ihi, d *scale, int *info) nogil
cdef zgebal_t *zgebal_f

ctypedef void zgebd2_t(int *m, int *n, z *a, int *lda, d *d, d *e, z *tauq, z *taup, z *work, int *info) nogil
cdef zgebd2_t *zgebd2_f

ctypedef void zgebrd_t(int *m, int *n, z *a, int *lda, d *d, d *e, z *tauq, z *taup, z *work, int *lwork, int *info) nogil
cdef zgebrd_t *zgebrd_f

ctypedef void zgecon_t(char *norm, int *n, z *a, int *lda, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zgecon_t *zgecon_f

ctypedef void zgeequ_t(int *m, int *n, z *a, int *lda, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef zgeequ_t *zgeequ_f

ctypedef void zgees_t(char *jobvs, char *sort, zselect1 *select, int *n, z *a, int *lda, int *sdim, z *w, z *vs, int *ldvs, z *work, int *lwork, d *rwork, bint *bwork, int *info) nogil
cdef zgees_t *zgees_f

ctypedef void zgeesx_t(char *jobvs, char *sort, zselect1 *select, char *sense, int *n, z *a, int *lda, int *sdim, z *w, z *vs, int *ldvs, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, bint *bwork, int *info) nogil
cdef zgeesx_t *zgeesx_f

ctypedef void zgeev_t(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *w, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgeev_t *zgeev_f

ctypedef void zgeevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, z *a, int *lda, z *w, z *vl, int *ldvl, z *vr, int *ldvr, int *ilo, int *ihi, d *scale, d *abnrm, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgeevx_t *zgeevx_f

ctypedef void zgegs_t(char *jobvsl, char *jobvsr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vsl, int *ldvsl, z *vsr, int *ldvsr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgegs_t *zgegs_f

ctypedef void zgegv_t(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgegv_t *zgegv_f

ctypedef void zgehd2_t(int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgehd2_t *zgehd2_f

ctypedef void zgehrd_t(int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgehrd_t *zgehrd_f

ctypedef void zgelq2_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgelq2_t *zgelq2_f

ctypedef void zgelqf_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgelqf_t *zgelqf_f

ctypedef void zgels_t(char *trans, int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, z *work, int *lwork, int *info) nogil
cdef zgels_t *zgels_f

ctypedef void zgelsd_t(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, d *s, d *rcond, int *rank, z *work, int *lwork, d *rwork, int *iwork, int *info) nogil
cdef zgelsd_t *zgelsd_f

ctypedef void zgelss_t(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, d *s, d *rcond, int *rank, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgelss_t *zgelss_f

ctypedef void zgelsx_t(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *jpvt, d *rcond, int *rank, z *work, d *rwork, int *info) nogil
cdef zgelsx_t *zgelsx_f

ctypedef void zgelsy_t(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *jpvt, d *rcond, int *rank, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgelsy_t *zgelsy_f

ctypedef void zgeql2_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgeql2_t *zgeql2_f

ctypedef void zgeqlf_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgeqlf_t *zgeqlf_f

ctypedef void zgeqp3_t(int *m, int *n, z *a, int *lda, int *jpvt, z *tau, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgeqp3_t *zgeqp3_f

ctypedef void zgeqpf_t(int *m, int *n, z *a, int *lda, int *jpvt, z *tau, z *work, d *rwork, int *info) nogil
cdef zgeqpf_t *zgeqpf_f

ctypedef void zgeqr2_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgeqr2_t *zgeqr2_f

ctypedef void zgeqrf_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgeqrf_t *zgeqrf_f

ctypedef void zgerfs_t(char *trans, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgerfs_t *zgerfs_f

ctypedef void zgerq2_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgerq2_t *zgerq2_f

ctypedef void zgerqf_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgerqf_t *zgerqf_f

ctypedef void zgesc2_t(int *n, z *a, int *lda, z *rhs, int *ipiv, int *jpiv, d *scale) nogil
cdef zgesc2_t *zgesc2_f

ctypedef void zgesdd_t(char *jobz, int *m, int *n, z *a, int *lda, d *s, z *u, int *ldu, z *vt, int *ldvt, z *work, int *lwork, d *rwork, int *iwork, int *info) nogil
cdef zgesdd_t *zgesdd_f

ctypedef void zgesv_t(int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgesv_t *zgesv_f

ctypedef void zgesvd_t(char *jobu, char *jobvt, int *m, int *n, z *a, int *lda, d *s, z *u, int *ldu, z *vt, int *ldvt, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgesvd_t *zgesvd_f

ctypedef void zgesvx_t(char *fact, char *trans, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, char *equed, d *r, d *c, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgesvx_t *zgesvx_f

ctypedef void zgetc2_t(int *n, z *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef zgetc2_t *zgetc2_f

ctypedef void zgetf2_t(int *m, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zgetf2_t *zgetf2_f

ctypedef void zgetrf_t(int *m, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zgetrf_t *zgetrf_f

ctypedef void zgetri_t(int *n, z *a, int *lda, int *ipiv, z *work, int *lwork, int *info) nogil
cdef zgetri_t *zgetri_f

ctypedef void zgetrs_t(char *trans, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgetrs_t *zgetrs_f

ctypedef void zggbak_t(char *job, char *side, int *n, int *ilo, int *ihi, d *lscale, d *rscale, int *m, z *v, int *ldv, int *info) nogil
cdef zggbak_t *zggbak_f

ctypedef void zggbal_t(char *job, int *n, z *a, int *lda, z *b, int *ldb, int *ilo, int *ihi, d *lscale, d *rscale, d *work, int *info) nogil
cdef zggbal_t *zggbal_f

ctypedef void zgges_t(char *jobvsl, char *jobvsr, char *sort, zselect2 *selctg, int *n, z *a, int *lda, z *b, int *ldb, int *sdim, z *alpha, z *beta, z *vsl, int *ldvsl, z *vsr, int *ldvsr, z *work, int *lwork, d *rwork, bint *bwork, int *info) nogil
cdef zgges_t *zgges_f

ctypedef void zggesx_t(char *jobvsl, char *jobvsr, char *sort, zselect2 *selctg, char *sense, int *n, z *a, int *lda, z *b, int *ldb, int *sdim, z *alpha, z *beta, z *vsl, int *ldvsl, z *vsr, int *ldvsr, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef zggesx_t *zggesx_f

ctypedef void zggev_t(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zggev_t *zggev_f

ctypedef void zggevx_t(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, int *ilo, int *ihi, d *lscale, d *rscale, d *abnrm, d *bbnrm, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, int *iwork, bint *bwork, int *info) nogil
cdef zggevx_t *zggevx_f

ctypedef void zggglm_t(int *n, int *m, int *p, z *a, int *lda, z *b, int *ldb, z *d, z *x, z *y, z *work, int *lwork, int *info) nogil
cdef zggglm_t *zggglm_f

ctypedef void zgghrd_t(char *compq, char *compz, int *n, int *ilo, int *ihi, z *a, int *lda, z *b, int *ldb, z *q, int *ldq, z *z, int *ldz, int *info) nogil
cdef zgghrd_t *zgghrd_f

ctypedef void zgglse_t(int *m, int *n, int *p, z *a, int *lda, z *b, int *ldb, z *c, z *d, z *x, z *work, int *lwork, int *info) nogil
cdef zgglse_t *zgglse_f

ctypedef void zggqrf_t(int *n, int *m, int *p, z *a, int *lda, z *taua, z *b, int *ldb, z *taub, z *work, int *lwork, int *info) nogil
cdef zggqrf_t *zggqrf_f

ctypedef void zggrqf_t(int *m, int *p, int *n, z *a, int *lda, z *taua, z *b, int *ldb, z *taub, z *work, int *lwork, int *info) nogil
cdef zggrqf_t *zggrqf_f

ctypedef void zggsvd_t(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, z *a, int *lda, z *b, int *ldb, d *alpha, d *beta, z *u, int *ldu, z *v, int *ldv, z *q, int *ldq, z *work, d *rwork, int *iwork, int *info) nogil
cdef zggsvd_t *zggsvd_f

ctypedef void zggsvp_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, z *a, int *lda, z *b, int *ldb, d *tola, d *tolb, int *k, int *l, z *u, int *ldu, z *v, int *ldv, z *q, int *ldq, int *iwork, d *rwork, z *tau, z *work, int *info) nogil
cdef zggsvp_t *zggsvp_f

ctypedef void zgtcon_t(char *norm, int *n, z *dl, z *d, z *du, z *du2, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zgtcon_t *zgtcon_f

ctypedef void zgtrfs_t(char *trans, int *n, int *nrhs, z *dl, z *d, z *du, z *dlf, z *df, z *duf, z *du2, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgtrfs_t *zgtrfs_f

ctypedef void zgtsv_t(int *n, int *nrhs, z *dl, z *d, z *du, z *b, int *ldb, int *info) nogil
cdef zgtsv_t *zgtsv_f

ctypedef void zgtsvx_t(char *fact, char *trans, int *n, int *nrhs, z *dl, z *d, z *du, z *dlf, z *df, z *duf, z *du2, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgtsvx_t *zgtsvx_f

ctypedef void zgttrf_t(int *n, z *dl, z *d, z *du, z *du2, int *ipiv, int *info) nogil
cdef zgttrf_t *zgttrf_f

ctypedef void zgttrs_t(char *trans, int *n, int *nrhs, z *dl, z *d, z *du, z *du2, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgttrs_t *zgttrs_f

ctypedef void zgtts2_t(int *itrans, int *n, int *nrhs, z *dl, z *d, z *du, z *du2, int *ipiv, z *b, int *ldb) nogil
cdef zgtts2_t *zgtts2_f

ctypedef void zhbev_t(char *jobz, char *uplo, int *n, int *kd, z *ab, int *ldab, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhbev_t *zhbev_f

ctypedef void zhbevd_t(char *jobz, char *uplo, int *n, int *kd, z *ab, int *ldab, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhbevd_t *zhbevd_f

ctypedef void zhbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, z *ab, int *ldab, z *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhbevx_t *zhbevx_f

ctypedef void zhbgst_t(char *vect, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, z *x, int *ldx, z *work, d *rwork, int *info) nogil
cdef zhbgst_t *zhbgst_f

ctypedef void zhbgv_t(char *jobz, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhbgv_t *zhbgv_f

ctypedef void zhbgvd_t(char *jobz, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhbgvd_t *zhbgvd_f

ctypedef void zhbgvx_t(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, z *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhbgvx_t *zhbgvx_f

ctypedef void zhbtrd_t(char *vect, char *uplo, int *n, int *kd, z *ab, int *ldab, d *d, d *e, z *q, int *ldq, z *work, int *info) nogil
cdef zhbtrd_t *zhbtrd_f

ctypedef void zhecon_t(char *uplo, int *n, z *a, int *lda, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zhecon_t *zhecon_f

ctypedef void zheev_t(char *jobz, char *uplo, int *n, z *a, int *lda, d *w, z *work, int *lwork, d *rwork, int *info) nogil
cdef zheev_t *zheev_f

ctypedef void zheevd_t(char *jobz, char *uplo, int *n, z *a, int *lda, d *w, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zheevd_t *zheevd_f

ctypedef void zheevr_t(char *jobz, char *range, char *uplo, int *n, z *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, int *isuppz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zheevr_t *zheevr_f

ctypedef void zheevx_t(char *jobz, char *range, char *uplo, int *n, z *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zheevx_t *zheevx_f

ctypedef void zhegs2_t(int *itype, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zhegs2_t *zhegs2_f

ctypedef void zhegst_t(int *itype, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zhegst_t *zhegst_f

ctypedef void zhegv_t(int *itype, char *jobz, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *w, z *work, int *lwork, d *rwork, int *info) nogil
cdef zhegv_t *zhegv_f

ctypedef void zhegvd_t(int *itype, char *jobz, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *w, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhegvd_t *zhegvd_f

ctypedef void zhegvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhegvx_t *zhegvx_f

ctypedef void zherfs_t(char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zherfs_t *zherfs_f

ctypedef void zhesv_t(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, z *work, int *lwork, int *info) nogil
cdef zhesv_t *zhesv_f

ctypedef void zhesvx_t(char *fact, char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zhesvx_t *zhesvx_f

ctypedef void zhetd2_t(char *uplo, int *n, z *a, int *lda, d *d, d *e, z *tau, int *info) nogil
cdef zhetd2_t *zhetd2_f

ctypedef void zhetf2_t(char *uplo, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zhetf2_t *zhetf2_f

ctypedef void zhetrd_t(char *uplo, int *n, z *a, int *lda, d *d, d *e, z *tau, z *work, int *lwork, int *info) nogil
cdef zhetrd_t *zhetrd_f

ctypedef void zhetrf_t(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *lwork, int *info) nogil
cdef zhetrf_t *zhetrf_f

ctypedef void zhetri_t(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *info) nogil
cdef zhetri_t *zhetri_f

ctypedef void zhetrs_t(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zhetrs_t *zhetrs_f

ctypedef void zhgeqz_t(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, z *h, int *ldh, z *t, int *ldt, z *alpha, z *beta, z *q, int *ldq, z *z, int *ldz, z *work, int *lwork, d *rwork, int *info) nogil
cdef zhgeqz_t *zhgeqz_f

ctypedef void zhpcon_t(char *uplo, int *n, z *ap, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zhpcon_t *zhpcon_f

ctypedef void zhpev_t(char *jobz, char *uplo, int *n, z *ap, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhpev_t *zhpev_f

ctypedef void zhpevd_t(char *jobz, char *uplo, int *n, z *ap, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhpevd_t *zhpevd_f

ctypedef void zhpevx_t(char *jobz, char *range, char *uplo, int *n, z *ap, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhpevx_t *zhpevx_f

ctypedef void zhpgst_t(int *itype, char *uplo, int *n, z *ap, z *bp, int *info) nogil
cdef zhpgst_t *zhpgst_f

ctypedef void zhpgv_t(int *itype, char *jobz, char *uplo, int *n, z *ap, z *bp, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhpgv_t *zhpgv_f

ctypedef void zhpgvd_t(int *itype, char *jobz, char *uplo, int *n, z *ap, z *bp, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhpgvd_t *zhpgvd_f

ctypedef void zhpgvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, z *ap, z *bp, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhpgvx_t *zhpgvx_f

ctypedef void zhprfs_t(char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zhprfs_t *zhprfs_f

ctypedef void zhpsv_t(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zhpsv_t *zhpsv_f

ctypedef void zhpsvx_t(char *fact, char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zhpsvx_t *zhpsvx_f

ctypedef void zhptrd_t(char *uplo, int *n, z *ap, d *d, d *e, z *tau, int *info) nogil
cdef zhptrd_t *zhptrd_f

ctypedef void zhptrf_t(char *uplo, int *n, z *ap, int *ipiv, int *info) nogil
cdef zhptrf_t *zhptrf_f

ctypedef void zhptri_t(char *uplo, int *n, z *ap, int *ipiv, z *work, int *info) nogil
cdef zhptri_t *zhptri_f

ctypedef void zhptrs_t(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zhptrs_t *zhptrs_f

ctypedef void zhsein_t(char *side, char *eigsrc, char *initv, bint *select, int *n, z *h, int *ldh, z *w, z *vl, int *ldvl, z *vr, int *ldvr, int *mm, int *m, z *work, d *rwork, int *ifaill, int *ifailr, int *info) nogil
cdef zhsein_t *zhsein_f

ctypedef void zhseqr_t(char *job, char *compz, int *n, int *ilo, int *ihi, z *h, int *ldh, z *w, z *z, int *ldz, z *work, int *lwork, int *info) nogil
cdef zhseqr_t *zhseqr_f

ctypedef void zlacn2_t(int *n, z *v, z *x, d *est, int *kase, int *isave) nogil
cdef zlacn2_t *zlacn2_f

ctypedef void zlacon_t(int *n, z *v, z *x, d *est, int *kase) nogil
cdef zlacon_t *zlacon_f

ctypedef d zlangb_t(char *norm, int *n, int *kl, int *ku, z *ab, int *ldab, d *work) nogil
cdef zlangb_t *zlangb_f

ctypedef d zlange_t(char *norm, int *m, int *n, z *a, int *lda, d *work) nogil
cdef zlange_t *zlange_f

ctypedef d zlangt_t(char *norm, int *n, z *dl, z *d, z *du) nogil
cdef zlangt_t *zlangt_f

ctypedef d zlanhb_t(char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef zlanhb_t *zlanhb_f

ctypedef d zlanhe_t(char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil
cdef zlanhe_t *zlanhe_f

ctypedef d zlanhp_t(char *norm, char *uplo, int *n, z *ap, d *work) nogil
cdef zlanhp_t *zlanhp_f

ctypedef d zlanhs_t(char *norm, int *n, z *a, int *lda, d *work) nogil
cdef zlanhs_t *zlanhs_f

ctypedef d zlanht_t(char *norm, int *n, d *d, z *e) nogil
cdef zlanht_t *zlanht_f

ctypedef d zlansb_t(char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef zlansb_t *zlansb_f

ctypedef d zlansp_t(char *norm, char *uplo, int *n, z *ap, d *work) nogil
cdef zlansp_t *zlansp_f

ctypedef d zlansy_t(char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil
cdef zlansy_t *zlansy_f

ctypedef d zlantb_t(char *norm, char *uplo, char *diag, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef zlantb_t *zlantb_f

ctypedef d zlantp_t(char *norm, char *uplo, char *diag, int *n, z *ap, d *work) nogil
cdef zlantp_t *zlantp_f

ctypedef d zlantr_t(char *norm, char *uplo, char *diag, int *m, int *n, z *a, int *lda, d *work) nogil
cdef zlantr_t *zlantr_f

ctypedef void zlarf_t(char *side, int *m, int *n, z *v, int *incv, z *tau, z *c, int *ldc, z *work) nogil
cdef zlarf_t *zlarf_f

ctypedef void zlarz_t(char *side, int *m, int *n, int *l, z *v, int *incv, z *tau, z *c, int *ldc, z *work) nogil
cdef zlarz_t *zlarz_f

ctypedef void zlaswp_t(int *n, z *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef zlaswp_t *zlaswp_f

ctypedef void zlauum_t(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zlauum_t *zlauum_f

ctypedef void zpbcon_t(char *uplo, int *n, int *kd, z *ab, int *ldab, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zpbcon_t *zpbcon_f

ctypedef void zpbequ_t(char *uplo, int *n, int *kd, z *ab, int *ldab, d *s, d *scond, d *amax, int *info) nogil
cdef zpbequ_t *zpbequ_f

ctypedef void zpbrfs_t(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zpbrfs_t *zpbrfs_f

ctypedef void zpbstf_t(char *uplo, int *n, int *kd, z *ab, int *ldab, int *info) nogil
cdef zpbstf_t *zpbstf_f

ctypedef void zpbsv_t(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
cdef zpbsv_t *zpbsv_f

ctypedef void zpbsvx_t(char *fact, char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, char *equed, d *s, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zpbsvx_t *zpbsvx_f

ctypedef void zpbtf2_t(char *uplo, int *n, int *kd, z *ab, int *ldab, int *info) nogil
cdef zpbtf2_t *zpbtf2_f

ctypedef void zpbtrf_t(char *uplo, int *n, int *kd, z *ab, int *ldab, int *info) nogil
cdef zpbtrf_t *zpbtrf_f

ctypedef void zpbtrs_t(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
cdef zpbtrs_t *zpbtrs_f

ctypedef void zpocon_t(char *uplo, int *n, z *a, int *lda, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zpocon_t *zpocon_f

ctypedef void zpoequ_t(int *n, z *a, int *lda, d *s, d *scond, d *amax, int *info) nogil
cdef zpoequ_t *zpoequ_f

ctypedef void zporfs_t(char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zporfs_t *zporfs_f

ctypedef void zposv_t(char *uplo, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zposv_t *zposv_f

ctypedef void zposvx_t(char *fact, char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, char *equed, d *s, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zposvx_t *zposvx_f

ctypedef void zpotf2_t(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zpotf2_t *zpotf2_f

ctypedef void zpotrf_t(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zpotrf_t *zpotrf_f

ctypedef void zpotri_t(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zpotri_t *zpotri_f

ctypedef void zpotrs_t(char *uplo, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zpotrs_t *zpotrs_f

ctypedef void zppcon_t(char *uplo, int *n, z *ap, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zppcon_t *zppcon_f

ctypedef void zppequ_t(char *uplo, int *n, z *ap, d *s, d *scond, d *amax, int *info) nogil
cdef zppequ_t *zppequ_f

ctypedef void zpprfs_t(char *uplo, int *n, int *nrhs, z *ap, z *afp, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zpprfs_t *zpprfs_f

ctypedef void zppsv_t(char *uplo, int *n, int *nrhs, z *ap, z *b, int *ldb, int *info) nogil
cdef zppsv_t *zppsv_f

ctypedef void zppsvx_t(char *fact, char *uplo, int *n, int *nrhs, z *ap, z *afp, char *equed, d *s, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zppsvx_t *zppsvx_f

ctypedef void zpptrf_t(char *uplo, int *n, z *ap, int *info) nogil
cdef zpptrf_t *zpptrf_f

ctypedef void zpptri_t(char *uplo, int *n, z *ap, int *info) nogil
cdef zpptri_t *zpptri_f

ctypedef void zpptrs_t(char *uplo, int *n, int *nrhs, z *ap, z *b, int *ldb, int *info) nogil
cdef zpptrs_t *zpptrs_f

ctypedef void zptcon_t(int *n, d *d, z *e, d *anorm, d *rcond, d *rwork, int *info) nogil
cdef zptcon_t *zptcon_f

ctypedef void zpteqr_t(char *compz, int *n, d *d, d *e, z *z, int *ldz, d *work, int *info) nogil
cdef zpteqr_t *zpteqr_f

ctypedef void zptrfs_t(char *uplo, int *n, int *nrhs, d *d, z *e, d *df, z *ef, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zptrfs_t *zptrfs_f

ctypedef void zptsv_t(int *n, int *nrhs, d *d, z *e, z *b, int *ldb, int *info) nogil
cdef zptsv_t *zptsv_f

ctypedef void zptsvx_t(char *fact, int *n, int *nrhs, d *d, z *e, d *df, z *ef, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zptsvx_t *zptsvx_f

ctypedef void zpttrf_t(int *n, d *d, z *e, int *info) nogil
cdef zpttrf_t *zpttrf_f

ctypedef void zpttrs_t(char *uplo, int *n, int *nrhs, d *d, z *e, z *b, int *ldb, int *info) nogil
cdef zpttrs_t *zpttrs_f

ctypedef void zptts2_t(int *iuplo, int *n, int *nrhs, d *d, z *e, z *b, int *ldb) nogil
cdef zptts2_t *zptts2_f

ctypedef void zrot_t(int *n, z *cx, int *incx, z *cy, int *incy, d *c, z *s) nogil
cdef zrot_t *zrot_f

ctypedef void zspcon_t(char *uplo, int *n, z *ap, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zspcon_t *zspcon_f

ctypedef void zspmv_t(char *uplo, int *n, z *alpha, z *ap, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zspmv_t *zspmv_f

ctypedef void zspr_t(char *uplo, int *n, z *alpha, z *x, int *incx, z *ap) nogil
cdef zspr_t *zspr_f

ctypedef void zsprfs_t(char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zsprfs_t *zsprfs_f

ctypedef void zspsv_t(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zspsv_t *zspsv_f

ctypedef void zspsvx_t(char *fact, char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zspsvx_t *zspsvx_f

ctypedef void zsptrf_t(char *uplo, int *n, z *ap, int *ipiv, int *info) nogil
cdef zsptrf_t *zsptrf_f

ctypedef void zsptri_t(char *uplo, int *n, z *ap, int *ipiv, z *work, int *info) nogil
cdef zsptri_t *zsptri_f

ctypedef void zsptrs_t(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zsptrs_t *zsptrs_f

ctypedef void zstedc_t(char *compz, int *n, d *d, d *e, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zstedc_t *zstedc_f

ctypedef void zstegr_t(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef zstegr_t *zstegr_f

ctypedef void zstein_t(int *n, d *d, d *e, int *m, d *w, int *iblock, int *isplit, z *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef zstein_t *zstein_f

ctypedef void zstemr_t(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, int *m, d *w, z *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef zstemr_t *zstemr_f

ctypedef void zsteqr_t(char *compz, int *n, d *d, d *e, z *z, int *ldz, d *work, int *info) nogil
cdef zsteqr_t *zsteqr_f

ctypedef void zsycon_t(char *uplo, int *n, z *a, int *lda, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zsycon_t *zsycon_f

ctypedef void zsymv_t(char *uplo, int *n, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zsymv_t *zsymv_f

ctypedef void zsyr_t(char *uplo, int *n, z *alpha, z *x, int *incx, z *a, int *lda) nogil
cdef zsyr_t *zsyr_f

ctypedef void zsyrfs_t(char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zsyrfs_t *zsyrfs_f

ctypedef void zsysv_t(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, z *work, int *lwork, int *info) nogil
cdef zsysv_t *zsysv_f

ctypedef void zsysvx_t(char *fact, char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zsysvx_t *zsysvx_f

ctypedef void zsytf2_t(char *uplo, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zsytf2_t *zsytf2_f

ctypedef void zsytrf_t(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *lwork, int *info) nogil
cdef zsytrf_t *zsytrf_f

ctypedef void zsytri_t(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *info) nogil
cdef zsytri_t *zsytri_f

ctypedef void zsytrs_t(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zsytrs_t *zsytrs_f

ctypedef void ztbcon_t(char *norm, char *uplo, char *diag, int *n, int *kd, z *ab, int *ldab, d *rcond, z *work, d *rwork, int *info) nogil
cdef ztbcon_t *ztbcon_f

ctypedef void ztbrfs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef ztbrfs_t *ztbrfs_f

ctypedef void ztbtrs_t(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
cdef ztbtrs_t *ztbtrs_f

ctypedef void ztgevc_t(char *side, char *howmny, bint *select, int *n, z *s, int *lds, z *p, int *ldp, z *vl, int *ldvl, z *vr, int *ldvr, int *mm, int *m, z *work, d *rwork, int *info) nogil
cdef ztgevc_t *ztgevc_f

ctypedef void ztgex2_t(bint *wantq, bint *wantz, int *n, z *a, int *lda, z *b, int *ldb, z *q, int *ldq, z *z, int *ldz, int *j1, int *info) nogil
cdef ztgex2_t *ztgex2_f

ctypedef void ztgexc_t(bint *wantq, bint *wantz, int *n, z *a, int *lda, z *b, int *ldb, z *q, int *ldq, z *z, int *ldz, int *ifst, int *ilst, int *info) nogil
cdef ztgexc_t *ztgexc_f

ctypedef void ztgsen_t(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *q, int *ldq, z *z, int *ldz, int *m, d *pl, d *pr, d *dif, z *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ztgsen_t *ztgsen_f

ctypedef void ztgsja_t(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, z *a, int *lda, z *b, int *ldb, d *tola, d *tolb, d *alpha, d *beta, z *u, int *ldu, z *v, int *ldv, z *q, int *ldq, z *work, int *ncycle, int *info) nogil
cdef ztgsja_t *ztgsja_f

ctypedef void ztgsna_t(char *job, char *howmny, bint *select, int *n, z *a, int *lda, z *b, int *ldb, z *vl, int *ldvl, z *vr, int *ldvr, d *s, d *dif, int *mm, int *m, z *work, int *lwork, int *iwork, int *info) nogil
cdef ztgsna_t *ztgsna_f

ctypedef void ztgsy2_t(char *trans, int *ijob, int *m, int *n, z *a, int *lda, z *b, int *ldb, z *c, int *ldc, z *d, int *ldd, z *e, int *lde, z *f, int *ldf, d *scale, d *rdsum, d *rdscal, int *info) nogil
cdef ztgsy2_t *ztgsy2_f

ctypedef void ztgsyl_t(char *trans, int *ijob, int *m, int *n, z *a, int *lda, z *b, int *ldb, z *c, int *ldc, z *d, int *ldd, z *e, int *lde, z *f, int *ldf, d *scale, d *dif, z *work, int *lwork, int *iwork, int *info) nogil
cdef ztgsyl_t *ztgsyl_f

ctypedef void ztpcon_t(char *norm, char *uplo, char *diag, int *n, z *ap, d *rcond, z *work, d *rwork, int *info) nogil
cdef ztpcon_t *ztpcon_f

ctypedef void ztprfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *ap, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef ztprfs_t *ztprfs_f

ctypedef void ztptri_t(char *uplo, char *diag, int *n, z *ap, int *info) nogil
cdef ztptri_t *ztptri_f

ctypedef void ztptrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *ap, z *b, int *ldb, int *info) nogil
cdef ztptrs_t *ztptrs_f

ctypedef void ztrcon_t(char *norm, char *uplo, char *diag, int *n, z *a, int *lda, d *rcond, z *work, d *rwork, int *info) nogil
cdef ztrcon_t *ztrcon_f

ctypedef void ztrevc_t(char *side, char *howmny, bint *select, int *n, z *t, int *ldt, z *vl, int *ldvl, z *vr, int *ldvr, int *mm, int *m, z *work, d *rwork, int *info) nogil
cdef ztrevc_t *ztrevc_f

ctypedef void ztrexc_t(char *compq, int *n, z *t, int *ldt, z *q, int *ldq, int *ifst, int *ilst, int *info) nogil
cdef ztrexc_t *ztrexc_f

ctypedef void ztrrfs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef ztrrfs_t *ztrrfs_f

ctypedef void ztrsen_t(char *job, char *compq, bint *select, int *n, z *t, int *ldt, z *q, int *ldq, z *w, int *m, d *s, d *sep, z *work, int *lwork, int *info) nogil
cdef ztrsen_t *ztrsen_f

ctypedef void ztrsna_t(char *job, char *howmny, bint *select, int *n, z *t, int *ldt, z *vl, int *ldvl, z *vr, int *ldvr, d *s, d *sep, int *mm, int *m, z *work, int *ldwork, d *rwork, int *info) nogil
cdef ztrsna_t *ztrsna_f

ctypedef void ztrsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, z *a, int *lda, z *b, int *ldb, z *c, int *ldc, d *scale, int *info) nogil
cdef ztrsyl_t *ztrsyl_f

ctypedef void ztrti2_t(char *uplo, char *diag, int *n, z *a, int *lda, int *info) nogil
cdef ztrti2_t *ztrti2_f

ctypedef void ztrtri_t(char *uplo, char *diag, int *n, z *a, int *lda, int *info) nogil
cdef ztrtri_t *ztrtri_f

ctypedef void ztrtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef ztrtrs_t *ztrtrs_f

ctypedef void ztzrqf_t(int *m, int *n, z *a, int *lda, z *tau, int *info) nogil
cdef ztzrqf_t *ztzrqf_f

ctypedef void ztzrzf_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef ztzrzf_t *ztzrzf_f

ctypedef void zung2l_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zung2l_t *zung2l_f

ctypedef void zung2r_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zung2r_t *zung2r_f

ctypedef void zungbr_t(char *vect, int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungbr_t *zungbr_f

ctypedef void zunghr_t(int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zunghr_t *zunghr_f

ctypedef void zungl2_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zungl2_t *zungl2_f

ctypedef void zunglq_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zunglq_t *zunglq_f

ctypedef void zungql_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungql_t *zungql_f

ctypedef void zungqr_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungqr_t *zungqr_f

ctypedef void zungr2_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zungr2_t *zungr2_f

ctypedef void zungrq_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungrq_t *zungrq_f

ctypedef void zungtr_t(char *uplo, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungtr_t *zungtr_f

ctypedef void zunm2l_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunm2l_t *zunm2l_f

ctypedef void zunm2r_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunm2r_t *zunm2r_f

ctypedef void zunmbr_t(char *vect, char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmbr_t *zunmbr_f

ctypedef void zunmhr_t(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmhr_t *zunmhr_f

ctypedef void zunml2_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunml2_t *zunml2_f

ctypedef void zunmlq_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmlq_t *zunmlq_f

ctypedef void zunmql_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmql_t *zunmql_f

ctypedef void zunmqr_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmqr_t *zunmqr_f

ctypedef void zunmr2_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunmr2_t *zunmr2_f

ctypedef void zunmr3_t(char *side, char *trans, int *m, int *n, int *k, int *l, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunmr3_t *zunmr3_f

ctypedef void zunmrq_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmrq_t *zunmrq_f

ctypedef void zunmrz_t(char *side, char *trans, int *m, int *n, int *k, int *l, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmrz_t *zunmrz_f

ctypedef void zunmtr_t(char *side, char *uplo, char *trans, int *m, int *n, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmtr_t *zunmtr_f

ctypedef void zupgtr_t(char *uplo, int *n, z *ap, z *tau, z *q, int *ldq, z *work, int *info) nogil
cdef zupgtr_t *zupgtr_f

ctypedef void zupmtr_t(char *side, char *uplo, char *trans, int *m, int *n, z *ap, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zupmtr_t *zupmtr_f
