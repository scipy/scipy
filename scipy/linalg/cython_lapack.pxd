ctypedef float s
ctypedef double d
ctypedef float complex c
ctypedef double complex z

# Function pointer type declarations for
# gees and gges families of functions.
ctypedef int cselect1(c*)
ctypedef int cselect2(c*, c*)
ctypedef int dselect2(d*, d*)
ctypedef int dselect3(d*, d*, d*)
ctypedef int sselect2(s*, s*)
ctypedef int sselect3(s*, s*, s*)
ctypedef int zselect1(z*)
ctypedef int zselect2(z*, z*)

ctypedef int cgbsv_t(int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, int *ipiv, c *b, int *ldb, int *info) nogil
ctypedef int cgbtrf_t(int *m, int *n, int *kl, int *ku, c *ab, int *ldab, int *ipiv, int *info) nogil
ctypedef int cgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, int *ipiv, c *b, int *ldb, int *info) nogil
ctypedef int cgebal_t(char *job, int *n, c *a, int *lda, int *ilo, int *ihi, s *scale, int *info) nogil
ctypedef int cgees_t(char *jobvs, char *sort, cselect1 *select, int *n, c *a, int *lda, int *sdim, c *w, c *vs, int *ldvs, c *work, int *lwork, s *rwork, int *bwork, int *info) nogil
ctypedef int cgeev_t(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *w, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
ctypedef int cgegv_t(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
ctypedef int cgehrd_t(int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
ctypedef int cgelss_t(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, s *s, s *rcond, int *rank, c *work, int *lwork, s *rwork, int *info) nogil
ctypedef int cgeqp3_t(int *m, int *n, c *a, int *lda, int *jpvt, c *tau, c *work, int *lwork, s *rwork, int *info) nogil
ctypedef int cgeqrf_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
ctypedef int cgerqf_t(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
ctypedef int cgesdd_t(char *jobz, int *m, int *n, c *a, int *lda, s *s, c *u, int *ldu, c *vt, int *ldvt, c *work, int *lwork, s *rwork, int *iwork, int *info) nogil
ctypedef int cgesv_t(int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
ctypedef int cgetrf_t(int *m, int *n, c *a, int *lda, int *ipiv, int *info) nogil
ctypedef int cgetri_t(int *n, c *a, int *lda, int *ipiv, c *work, int *lwork, int *info) nogil
ctypedef int cgetrs_t(char *trans, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
ctypedef int cgges_t(char *jobvsl, char *jobvsr, char *sort, cselect2 *selctg, int *n, c *a, int *lda, c *b, int *ldb, int *sdim, c *alpha, c *beta, c *vsl, int *ldvsl, c *vsr, int *ldvsr, c *work, int *lwork, s *rwork, int *bwork, int *info) nogil
ctypedef int cggev_t(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
ctypedef int chbevd_t(char *jobz, char *uplo, int *n, int *kd, c *ab, int *ldab, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int chbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, c *ab, int *ldab, c *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
ctypedef int cheev_t(char *jobz, char *uplo, int *n, c *a, int *lda, s *w, c *work, int *lwork, s *rwork, int *info) nogil
ctypedef int cheevd_t(char *jobz, char *uplo, int *n, c *a, int *lda, s *w, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int cheevr_t(char *jobz, char *range, char *uplo, int *n, c *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, int *isuppz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int chegv_t(int *itype, char *jobz, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *w, c *work, int *lwork, s *rwork, int *info) nogil
ctypedef int chegvd_t(int *itype, char *jobz, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *w, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int chegvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *iwork, int *ifail, int *info) nogil
ctypedef int claswp_t(int *n, c *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
ctypedef int clauum_t(char *uplo, int *n, c *a, int *lda, int *info) nogil
ctypedef int cpbsv_t(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
ctypedef int cpbtrf_t(char *uplo, int *n, int *kd, c *ab, int *ldab, int *info) nogil
ctypedef int cpbtrs_t(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
ctypedef int cposv_t(char *uplo, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
ctypedef int cpotrf_t(char *uplo, int *n, c *a, int *lda, int *info) nogil
ctypedef int cpotri_t(char *uplo, int *n, c *a, int *lda, int *info) nogil
ctypedef int cpotrs_t(char *uplo, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
ctypedef int ctrsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, c *a, int *lda, c *b, int *ldb, c *c, int *ldc, s *scale, int *info) nogil
ctypedef int ctrtri_t(char *uplo, char *diag, int *n, c *a, int *lda, int *info) nogil
ctypedef int ctrtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
ctypedef int cungqr_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
ctypedef int cungrq_t(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
ctypedef int cunmqr_t(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
ctypedef int dgbsv_t(int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, int *ipiv, d *b, int *ldb, int *info) nogil
ctypedef int dgbtrf_t(int *m, int *n, int *kl, int *ku, d *ab, int *ldab, int *ipiv, int *info) nogil
ctypedef int dgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, int *ipiv, d *b, int *ldb, int *info) nogil
ctypedef int dgebal_t(char *job, int *n, d *a, int *lda, int *ilo, int *ihi, d *scale, int *info) nogil
ctypedef int dgees_t(char *jobvs, char *sort, dselect2 *select, int *n, d *a, int *lda, int *sdim, d *wr, d *wi, d *vs, int *ldvs, d *work, int *lwork, int *bwork, int *info) nogil
ctypedef int dgeev_t(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *wr, d *wi, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
ctypedef int dgegv_t(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *b, int *ldb, d *alpha, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, s *rwork, int *info) nogil
ctypedef int dgehrd_t(int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
ctypedef int dgelss_t(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *s, d *rcond, int *rank, d *work, int *lwork, int *info) nogil
ctypedef int dgeqp3_t(int *m, int *n, d *a, int *lda, int *jpvt, d *tau, d *work, int *lwork, int *info) nogil
ctypedef int dgeqrf_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
ctypedef int dgerqf_t(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
ctypedef int dgesdd_t(char *jobz, int *m, int *n, d *a, int *lda, d *s, d *u, int *ldu, d *vt, int *ldvt, d *work, int *lwork, int *iwork, int *info) nogil
ctypedef int dgesv_t(int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
ctypedef int dgetrf_t(int *m, int *n, d *a, int *lda, int *ipiv, int *info) nogil
ctypedef int dgetri_t(int *n, d *a, int *lda, int *ipiv, d *work, int *lwork, int *info) nogil
ctypedef int dgetrs_t(char *trans, int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
ctypedef int dgges_t(char *jobvsl, char *jobvsr, char *sort, dselect3 *selctg, int *n, d *a, int *lda, d *b, int *ldb, int *sdim, d *alphar, d *alphai, d *beta, d *vsl, int *ldvsl, d *vsr, int *ldvsr, d *work, int *lwork, int *bwork, int *info) nogil
ctypedef int dggev_t(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
ctypedef d dlamch_t(char *cmach) nogil
ctypedef int dlaswp_t(int *n, d *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
ctypedef int dlauum_t(char *uplo, int *n, d *a, int *lda, int *info) nogil
ctypedef int dorgqr_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
ctypedef int dorgrq_t(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
ctypedef int dormqr_t(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
ctypedef int dpbsv_t(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
ctypedef int dpbtrf_t(char *uplo, int *n, int *kd, d *ab, int *ldab, int *info) nogil
ctypedef int dpbtrs_t(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
ctypedef int dposv_t(char *uplo, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
ctypedef int dpotrf_t(char *uplo, int *n, d *a, int *lda, int *info) nogil
ctypedef int dpotri_t(char *uplo, int *n, d *a, int *lda, int *info) nogil
ctypedef int dpotrs_t(char *uplo, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
ctypedef int dsbev_t(char *jobz, char *uplo, int *n, int *kd, d *ab, int *ldab, d *w, d *z, int *ldz, d *work, int *info) nogil
ctypedef int dsbevd_t(char *jobz, char *uplo, int *n, int *kd, d *ab, int *ldab, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int dsbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, d *ab, int *ldab, d *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
ctypedef int dsyev_t(char *jobz, char *uplo, int *n, d *a, int *lda, d *w, d *work, int *lwork, int *info) nogil
ctypedef int dsyevd_t(char *jobz, char *uplo, int *n, d *a, int *lda, d *w, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int dsyevr_t(char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int dsygv_t(int *itype, char *jobz, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *w, d *work, int *lwork, int *info) nogil
ctypedef int dsygvd_t(int *itype, char *jobz, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *w, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int dsygvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *ifail, int *info) nogil
ctypedef int dtrsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, d *a, int *lda, d *b, int *ldb, d *c, int *ldc, d *scale, int *info) nogil
ctypedef int dtrtri_t(char *uplo, char *diag, int *n, d *a, int *lda, int *info) nogil
ctypedef int dtrtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
ctypedef int sgbsv_t(int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, int *ipiv, s *b, int *ldb, int *info) nogil
ctypedef int sgbtrf_t(int *m, int *n, int *kl, int *ku, s *ab, int *ldab, int *ipiv, int *info) nogil
ctypedef int sgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, int *ipiv, s *b, int *ldb, int *info) nogil
ctypedef int sgebal_t(char *job, int *n, s *a, int *lda, int *ilo, int *ihi, s *scale, int *info) nogil
ctypedef int sgees_t(char *jobvs, char *sort, sselect2 *select, int *n, s *a, int *lda, int *sdim, s *wr, s *wi, s *vs, int *ldvs, s *work, int *lwork, int *bwork, int *info) nogil
ctypedef int sgeev_t(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *wr, s *wi, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
ctypedef int sgegv_t(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *b, int *ldb, s *alpha, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, s *rwork, int *info) nogil
ctypedef int sgehrd_t(int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
ctypedef int sgelss_t(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *s, s *rcond, int *rank, s *work, int *lwork, int *info) nogil
ctypedef int sgeqp3_t(int *m, int *n, s *a, int *lda, int *jpvt, s *tau, s *work, int *lwork, int *info) nogil
ctypedef int sgeqrf_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
ctypedef int sgerqf_t(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
ctypedef int sgesdd_t(char *jobz, int *m, int *n, s *a, int *lda, s *s, s *u, int *ldu, s *vt, int *ldvt, s *work, int *lwork, int *iwork, int *info) nogil
ctypedef int sgesv_t(int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
ctypedef int sgetrf_t(int *m, int *n, s *a, int *lda, int *ipiv, int *info) nogil
ctypedef int sgetri_t(int *n, s *a, int *lda, int *ipiv, s *work, int *lwork, int *info) nogil
ctypedef int sgetrs_t(char *trans, int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
ctypedef int sgges_t(char *jobvsl, char *jobvsr, char *sort, sselect3 *selctg, int *n, s *a, int *lda, s *b, int *ldb, int *sdim, s *alphar, s *alphai, s *beta, s *vsl, int *ldvsl, s *vsr, int *ldvsr, s *work, int *lwork, int *bwork, int *info) nogil
ctypedef int sggev_t(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
ctypedef s slamch_t(char *cmach) nogil
ctypedef int slaswp_t(int *n, s *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
ctypedef int slauum_t(char *uplo, int *n, s *a, int *lda, int *info) nogil
ctypedef int sorgqr_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
ctypedef int sorgrq_t(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
ctypedef int sormqr_t(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
ctypedef int spbsv_t(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
ctypedef int spbtrf_t(char *uplo, int *n, int *kd, s *ab, int *ldab, int *info) nogil
ctypedef int spbtrs_t(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
ctypedef int sposv_t(char *uplo, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
ctypedef int spotrf_t(char *uplo, int *n, s *a, int *lda, int *info) nogil
ctypedef int spotri_t(char *uplo, int *n, s *a, int *lda, int *info) nogil
ctypedef int spotrs_t(char *uplo, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
ctypedef int ssbev_t(char *jobz, char *uplo, int *n, int *kd, s *ab, int *ldab, s *w, s *z, int *ldz, s *work, int *info) nogil
ctypedef int ssbevd_t(char *jobz, char *uplo, int *n, int *kd, s *ab, int *ldab, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int ssbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, s *ab, int *ldab, s *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
ctypedef int ssyev_t(char *jobz, char *uplo, int *n, s *a, int *lda, s *w, s *work, int *lwork, int *info) nogil
ctypedef int ssyevd_t(char *jobz, char *uplo, int *n, s *a, int *lda, s *w, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int ssyevr_t(char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int ssygv_t(int *itype, char *jobz, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *w, s *work, int *lwork, int *info) nogil
ctypedef int ssygvd_t(int *itype, char *jobz, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *w, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
ctypedef int ssygvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *ifail, int *info) nogil
ctypedef int strsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, s *a, int *lda, s *b, int *ldb, s *c, int *ldc, s *scale, int *info) nogil
ctypedef int strtri_t(char *uplo, char *diag, int *n, s *a, int *lda, int *info) nogil
ctypedef int strtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
ctypedef int zgbsv_t(int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, int *ipiv, z *b, int *ldb, int *info) nogil
ctypedef int zgbtrf_t(int *m, int *n, int *kl, int *ku, z *ab, int *ldab, int *ipiv, int *info) nogil
ctypedef int zgbtrs_t(char *trans, int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, int *ipiv, z *b, int *ldb, int *info) nogil
ctypedef int zgebal_t(char *job, int *n, z *a, int *lda, int *ilo, int *ihi, d *scale, int *info) nogil
ctypedef int zgees_t(char *jobvs, char *sort, zselect1 *select, int *n, z *a, int *lda, int *sdim, z *w, z *vs, int *ldvs, z *work, int *lwork, d *rwork, int *bwork, int *info) nogil
ctypedef int zgeev_t(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *w, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
ctypedef int zgegv_t(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, s *rwork, int *info) nogil
ctypedef int zgehrd_t(int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
ctypedef int zgelss_t(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, d *s, d *rcond, int *rank, z *work, int *lwork, d *rwork, int *info) nogil
ctypedef int zgeqp3_t(int *m, int *n, z *a, int *lda, int *jpvt, z *tau, z *work, int *lwork, d *rwork, int *info) nogil
ctypedef int zgeqrf_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
ctypedef int zgerqf_t(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
ctypedef int zgesdd_t(char *jobz, int *m, int *n, z *a, int *lda, d *s, z *u, int *ldu, z *vt, int *ldvt, z *work, int *lwork, d *rwork, int *iwork, int *info) nogil
ctypedef int zgesv_t(int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
ctypedef int zgetrf_t(int *m, int *n, z *a, int *lda, int *ipiv, int *info) nogil
ctypedef int zgetri_t(int *n, z *a, int *lda, int *ipiv, z *work, int *lwork, int *info) nogil
ctypedef int zgetrs_t(char *trans, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
ctypedef int zgges_t(char *jobvsl, char *jobvsr, char *sort, zselect2 *selctg, int *n, z *a, int *lda, z *b, int *ldb, int *sdim, z *alpha, z *beta, z *vsl, int *ldvsl, z *vsr, int *ldvsr, z *work, int *lwork, d *rwork, int *bwork, int *info) nogil
ctypedef int zggev_t(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
ctypedef int zhbevd_t(char *jobz, char *uplo, int *n, int *kd, z *ab, int *ldab, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int zhbevx_t(char *jobz, char *range, char *uplo, int *n, int *kd, z *ab, int *ldab, z *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
ctypedef int zheev_t(char *jobz, char *uplo, int *n, z *a, int *lda, d *w, z *work, int *lwork, d *rwork, int *info) nogil
ctypedef int zheevd_t(char *jobz, char *uplo, int *n, z *a, int *lda, d *w, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int zheevr_t(char *jobz, char *range, char *uplo, int *n, z *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, int *isuppz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int zhegv_t(int *itype, char *jobz, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *w, z *work, int *lwork, d *rwork, int *info) nogil
ctypedef int zhegvd_t(int *itype, char *jobz, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *w, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
ctypedef int zhegvx_t(int *itype, char *jobz, char *range, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *iwork, int *ifail, int *info) nogil
ctypedef int zlaswp_t(int *n, z *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
ctypedef int zlauum_t(char *uplo, int *n, z *a, int *lda, int *info) nogil
ctypedef int zpbsv_t(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
ctypedef int zpbtrf_t(char *uplo, int *n, int *kd, z *ab, int *ldab, int *info) nogil
ctypedef int zpbtrs_t(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
ctypedef int zposv_t(char *uplo, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
ctypedef int zpotrf_t(char *uplo, int *n, z *a, int *lda, int *info) nogil
ctypedef int zpotri_t(char *uplo, int *n, z *a, int *lda, int *info) nogil
ctypedef int zpotrs_t(char *uplo, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
ctypedef int ztrsyl_t(char *trana, char *tranb, int *isgn, int *m, int *n, z *a, int *lda, z *b, int *ldb, z *c, int *ldc, d *scale, int *info) nogil
ctypedef int ztrtri_t(char *uplo, char *diag, int *n, z *a, int *lda, int *info) nogil
ctypedef int ztrtrs_t(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
ctypedef int zungqr_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
ctypedef int zungrq_t(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
ctypedef int zunmqr_t(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil

cdef:
    cgbsv_t *cgbsv_f
    cgbtrf_t *cgbtrf_f
    cgbtrs_t *cgbtrs_f
    cgebal_t *cgebal_f
    cgees_t *cgees_f
    cgeev_t *cgeev_f
    cgegv_t *cgegv_f
    cgehrd_t *cgehrd_f
    cgelss_t *cgelss_f
    cgeqp3_t *cgeqp3_f
    cgeqrf_t *cgeqrf_f
    cgerqf_t *cgerqf_f
    cgesdd_t *cgesdd_f
    cgesv_t *cgesv_f
    cgetrf_t *cgetrf_f
    cgetri_t *cgetri_f
    cgetrs_t *cgetrs_f
    cgges_t *cgges_f
    cggev_t *cggev_f
    chbevd_t *chbevd_f
    chbevx_t *chbevx_f
    cheev_t *cheev_f
    cheevd_t *cheevd_f
    cheevr_t *cheevr_f
    chegv_t *chegv_f
    chegvd_t *chegvd_f
    chegvx_t *chegvx_f
    claswp_t *claswp_f
    clauum_t *clauum_f
    cpbsv_t *cpbsv_f
    cpbtrf_t *cpbtrf_f
    cpbtrs_t *cpbtrs_f
    cposv_t *cposv_f
    cpotrf_t *cpotrf_f
    cpotri_t *cpotri_f
    cpotrs_t *cpotrs_f
    ctrsyl_t *ctrsyl_f
    ctrtri_t *ctrtri_f
    ctrtrs_t *ctrtrs_f
    cungqr_t *cungqr_f
    cungrq_t *cungrq_f
    cunmqr_t *cunmqr_f
    dgbsv_t *dgbsv_f
    dgbtrf_t *dgbtrf_f
    dgbtrs_t *dgbtrs_f
    dgebal_t *dgebal_f
    dgees_t *dgees_f
    dgeev_t *dgeev_f
    dgegv_t *dgegv_f
    dgehrd_t *dgehrd_f
    dgelss_t *dgelss_f
    dgeqp3_t *dgeqp3_f
    dgeqrf_t *dgeqrf_f
    dgerqf_t *dgerqf_f
    dgesdd_t *dgesdd_f
    dgesv_t *dgesv_f
    dgetrf_t *dgetrf_f
    dgetri_t *dgetri_f
    dgetrs_t *dgetrs_f
    dgges_t *dgges_f
    dggev_t *dggev_f
    dlamch_t *dlamch_f
    dlaswp_t *dlaswp_f
    dlauum_t *dlauum_f
    dorgqr_t *dorgqr_f
    dorgrq_t *dorgrq_f
    dormqr_t *dormqr_f
    dpbsv_t *dpbsv_f
    dpbtrf_t *dpbtrf_f
    dpbtrs_t *dpbtrs_f
    dposv_t *dposv_f
    dpotrf_t *dpotrf_f
    dpotri_t *dpotri_f
    dpotrs_t *dpotrs_f
    dsbev_t *dsbev_f
    dsbevd_t *dsbevd_f
    dsbevx_t *dsbevx_f
    dsyev_t *dsyev_f
    dsyevd_t *dsyevd_f
    dsyevr_t *dsyevr_f
    dsygv_t *dsygv_f
    dsygvd_t *dsygvd_f
    dsygvx_t *dsygvx_f
    dtrsyl_t *dtrsyl_f
    dtrtri_t *dtrtri_f
    dtrtrs_t *dtrtrs_f
    sgbsv_t *sgbsv_f
    sgbtrf_t *sgbtrf_f
    sgbtrs_t *sgbtrs_f
    sgebal_t *sgebal_f
    sgees_t *sgees_f
    sgeev_t *sgeev_f
    sgegv_t *sgegv_f
    sgehrd_t *sgehrd_f
    sgelss_t *sgelss_f
    sgeqp3_t *sgeqp3_f
    sgeqrf_t *sgeqrf_f
    sgerqf_t *sgerqf_f
    sgesdd_t *sgesdd_f
    sgesv_t *sgesv_f
    sgetrf_t *sgetrf_f
    sgetri_t *sgetri_f
    sgetrs_t *sgetrs_f
    sgges_t *sgges_f
    sggev_t *sggev_f
    slamch_t *slamch_f
    slaswp_t *slaswp_f
    slauum_t *slauum_f
    sorgqr_t *sorgqr_f
    sorgrq_t *sorgrq_f
    sormqr_t *sormqr_f
    spbsv_t *spbsv_f
    spbtrf_t *spbtrf_f
    spbtrs_t *spbtrs_f
    sposv_t *sposv_f
    spotrf_t *spotrf_f
    spotri_t *spotri_f
    spotrs_t *spotrs_f
    ssbev_t *ssbev_f
    ssbevd_t *ssbevd_f
    ssbevx_t *ssbevx_f
    ssyev_t *ssyev_f
    ssyevd_t *ssyevd_f
    ssyevr_t *ssyevr_f
    ssygv_t *ssygv_f
    ssygvd_t *ssygvd_f
    ssygvx_t *ssygvx_f
    strsyl_t *strsyl_f
    strtri_t *strtri_f
    strtrs_t *strtrs_f
    zgbsv_t *zgbsv_f
    zgbtrf_t *zgbtrf_f
    zgbtrs_t *zgbtrs_f
    zgebal_t *zgebal_f
    zgees_t *zgees_f
    zgeev_t *zgeev_f
    zgegv_t *zgegv_f
    zgehrd_t *zgehrd_f
    zgelss_t *zgelss_f
    zgeqp3_t *zgeqp3_f
    zgeqrf_t *zgeqrf_f
    zgerqf_t *zgerqf_f
    zgesdd_t *zgesdd_f
    zgesv_t *zgesv_f
    zgetrf_t *zgetrf_f
    zgetri_t *zgetri_f
    zgetrs_t *zgetrs_f
    zgges_t *zgges_f
    zggev_t *zggev_f
    zhbevd_t *zhbevd_f
    zhbevx_t *zhbevx_f
    zheev_t *zheev_f
    zheevd_t *zheevd_f
    zheevr_t *zheevr_f
    zhegv_t *zhegv_f
    zhegvd_t *zhegvd_f
    zhegvx_t *zhegvx_f
    zlaswp_t *zlaswp_f
    zlauum_t *zlauum_f
    zpbsv_t *zpbsv_f
    zpbtrf_t *zpbtrf_f
    zpbtrs_t *zpbtrs_f
    zposv_t *zposv_f
    zpotrf_t *zpotrf_f
    zpotri_t *zpotri_f
    zpotrs_t *zpotrs_f
    ztrsyl_t *ztrsyl_f
    ztrtri_t *ztrtri_f
    ztrtrs_t *ztrtrs_f
    zungqr_t *zungqr_f
    zungrq_t *zungrq_f
    zunmqr_t *zunmqr_f

cpdef double dlamch(char cmach)
cpdef float slamch(char cmach)
