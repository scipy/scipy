#ifndef LAPACK_H
#define LAPACK_H

#include "fortran.h"

/* typedefs for the FORTRAN COMPLEX*8 and COMPLEX*16 data types */

/* COMPLEX*8  type */

#ifndef __COMPLEX8_T
#define __COMPLEX8_T 1
    typedef struct {
	float  re, im;
    } complex8_t;
#endif

/* COMPLEX*16 type */

#ifndef __COMPLEX16_T
#define __COMPLEX16_T 1
    typedef struct {
	double re, im;
    } complex16_t;
#endif

extern void F77(cgbcon)(char* norm, int* n, int* kl, int* ku,
        complex8_t ab[], int* ldab, int ipiv[], float* anorm,
        float* rcond, complex8_t work[], float rwork[], int* info,
        int len_norm);

extern void F77(cgbsv)(int* n, int* kl, int* ku, int* nrhs, complex8_t ab[],
        int* ldab, int ipiv[], complex8_t b[], int* ldb, int* info);

extern void F77(cgbsvx)(char* fact, char* trans, int* n, int* kl, int* ku,
        int* nrhs, complex8_t ab[], int* ldab, complex8_t afb[],
        int* ldafb, int ipiv[], char* equed, float r[], float c[],
        complex8_t b[], int* ldb, complex8_t x[], int* ldx,
        float* rcond, float ferr[], float berr[], complex8_t work[],
        float rwork[], int* info, int len_fact, int len_trans,
        int len_equed);

extern void F77(cgbtrf)(int* m, int* n, int* kl, int* ku, complex8_t ab[],
        int* ldab, int ipiv[], int* info);

extern void F77(cgbtrs)(char* trans, int* n, int* kl, int* ku, int* nrhs,
        complex8_t ab[], int* ldab, int ipiv[], complex8_t b[],
        int* ldb, int* info, int len_trans);

extern void F77(cgecon)(char* norm, int* n, complex8_t a[], int* lda,
        float* anorm, float* rcond, complex8_t work[], float rwork[],
        int* info, int len_norm);

extern void F77(cgees)(char* jobvs, char* sort, int (*select)(), int* n,
        complex8_t a[], int* lda, int* sdim, complex8_t w[],
        complex8_t vs[], int* ldvs, complex8_t work[], int* lwork,
        float rwork[], int bwork[], int* info, int len_jobvs,
        int len_sort);

extern void F77(cgeesx)(char* jobvs, char* sort, int (*select)(),
        char* sense, int* n, complex8_t a[], int* lda, int* sdim,
        complex8_t w[], complex8_t vs[], int* ldvs, float* rconde,
        float* rcondv, complex8_t work[], int* lwork, float rwork[],
        int bwork[], int* info, int len_jobvs, int len_sort,
        int len_sense);

extern void F77(cgeev)(char* jobvl, char* jobvr, int* n, complex8_t a[],
        int* lda, complex8_t w[], complex8_t vl[], int* ldvl,
        complex8_t vr[], int* ldvr, complex8_t work[], int* lwork,
        float rwork[], int* info, int len_jobvl, int len_jobvr);

extern void F77(cgeevx)(char* balanc, char* jobvl, char* jobvr, char* sense,
        int* n, complex8_t a[], int* lda, complex8_t w[],
        complex8_t vl[], int* ldvl, complex8_t vr[], int* ldvr,
        int* ilo, int* ihi, float scale[], float* abnrm, float rconde[],
        float rcondv[], complex8_t work[], int* lwork, float rwork[],
        int* info, int len_balanc, int len_jobvl, int len_jobvr,
        int len_sense);

extern void F77(cgegs)(char* jobvsl, char* jobvsr, int* n, complex8_t a[],
        int* lda, complex8_t b[], int* ldb, complex8_t alpha[],
        complex8_t beta[], complex8_t vsl[], int* ldvsl,
        complex8_t vsr[], int* ldvsr, complex8_t work[], int* lwork,
        float rwork[], int* info, int len_jobvsl, int len_jobvsr);

extern void F77(cgegv)(char* jobvl, char* jobvr, int* n, complex8_t a[],
        int* lda, complex8_t b[], int* ldb, complex8_t alpha[],
        complex8_t beta[], complex8_t vl[], int* ldvl, complex8_t vr[],
        int* ldvr, complex8_t work[], int* lwork, float rwork[],
        int* info, int len_jobvl, int len_jobvr);

extern void F77(cgelqf)(int* m, int* n, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cgels)(char* trans, int* m, int* n, int* nrhs,
        complex8_t a[], int* lda, complex8_t b[], int* ldb,
        complex8_t work[], int* lwork, int* info, int len_trans);

extern void F77(cgelss)(int* m, int* n, int* nrhs, complex8_t a[], int* lda,
        complex8_t b[], int* ldb, float s[], float* rcond, int* rank,
        complex8_t work[], int* lwork, float rwork[], int* info);

extern void F77(cgelsx)(int* m, int* n, int* nrhs, complex8_t a[], int* lda,
        complex8_t b[], int* ldb, int jpvt[], float* rcond, int* rank,
        complex8_t work[], float rwork[], int* info);

extern void F77(cgeqlf)(int* m, int* n, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cgeqpf)(int* m, int* n, complex8_t a[], int* lda,
        int jpvt[], complex8_t tau[], complex8_t work[], float rwork[],
        int* info);

extern void F77(cgeqrf)(int* m, int* n, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cgerqf)(int* m, int* n, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cgesv)(int* n, int* nrhs, complex8_t a[], int* lda,
        int ipiv[], complex8_t b[], int* ldb, int* info);

extern void F77(cgesvd)(char* jobu, char* jobvt, int* m, int* n,
        complex8_t a[], int* lda, float s[], complex8_t u[], int* ldu,
        complex8_t vt[], int* ldvt, complex8_t work[], int* lwork,
        float rwork[], int* info, int len_jobu, int len_jobvt);

extern void F77(cgesvx)(char* fact, char* trans, int* n, int* nrhs,
        complex8_t a[], int* lda, complex8_t af[], int* ldaf,
        int ipiv[], char* equed, float r[], float c[], complex8_t b[],
        int* ldb, complex8_t x[], int* ldx, float* rcond, float ferr[],
        float berr[], complex8_t work[], float rwork[], int* info,
        int len_fact, int len_trans, int len_equed);

extern void F77(cgetrf)(int* m, int* n, complex8_t a[], int* lda,
        int ipiv[], int* info);

extern void F77(cgetri)(int* n, complex8_t a[], int* lda, int ipiv[],
        complex8_t work[], int* lwork, int* info);

extern void F77(cgetrs)(char* trans, int* n, int* nrhs, complex8_t a[],
        int* lda, int ipiv[], complex8_t b[], int* ldb, int* info,
        int len_trans);

extern void F77(cggglm)(int* n, int* m, int* p, complex8_t a[], int* lda,
        complex8_t b[], int* ldb, complex8_t d[], complex8_t x[],
        complex8_t y[], complex8_t work[], int* lwork, int* info);

extern void F77(cgglse)(int* m, int* n, int* p, complex8_t a[], int* lda,
        complex8_t b[], int* ldb, complex8_t c[], complex8_t d[],
        complex8_t x[], complex8_t work[], int* lwork, int* info);

extern void F77(cggqrf)(int* m, int* n, int* p, complex8_t a[], int* lda,
        complex8_t taua[], complex8_t b[], int* ldb, complex8_t taub[],
        complex8_t work[], int* lwork, int* info);

extern void F77(cggrqf)(int* m, int* n, int* p, complex8_t a[], int* lda,
        complex8_t taua[], complex8_t b[], int* ldb, complex8_t taub[],
        complex8_t work[], int* lwork, int* info);

extern void F77(cggsvd)(char* jobu, char* jobv, char* jobq, int* m, int* n,
        int* p, int* k, int* l, complex8_t a[], int* lda,
        complex8_t b[], int* ldb, float alpha[], float beta[],
        complex8_t u[], int* ldu, complex8_t v[], int* ldv,
        complex8_t q[], int* ldq, complex8_t work[], float rwork[],
        int iwork[], int* info, int len_jobu, int len_jobv,
        int len_jobq);

extern void F77(cgtcon)(char* norm, int* n, complex8_t dl[], complex8_t d[],
        complex8_t du[], complex8_t du2[], int ipiv[], float* anorm,
        float* rcond, complex8_t work[], int* info, int len_norm);

extern void F77(cgtsv)(int* n, int* nrhs, complex8_t dl[], complex8_t d[],
        complex8_t du[], complex8_t b[], int* ldb, int* info);

extern void F77(cgtsvx)(char* fact, char* trans, int* n, int* nrhs,
        complex8_t dl[], complex8_t d[], complex8_t du[],
        complex8_t dlf[], complex8_t df[], complex8_t duf[],
        complex8_t du2[], int ipiv[], complex8_t b[], int* ldb,
        complex8_t x[], int* ldx, float* rcond, float ferr[],
        float berr[], complex8_t work[], float rwork[], int* info,
        int len_fact, int len_trans);

extern void F77(cgttrf)(int* n, complex8_t dl[], complex8_t d[],
        complex8_t du[], complex8_t du2[], int ipiv[], int* info);

extern void F77(cgttrs)(char* trans, int* n, int* nrhs, complex8_t dl[],
        complex8_t d[], complex8_t du[], complex8_t du2[], int ipiv[],
        complex8_t b[], int* ldb, int* info, int len_trans);

extern void F77(chbev)(char* jobz, char* uplo, int* n, int* kd,
        complex8_t ab[], int* ldab, float w[], complex8_t z[], int* ldz,
        complex8_t work[], float rwork[], int* info, int len_jobz,
        int len_uplo);

extern void F77(chbevd)(char* jobz, char* uplo, int* n, int* kd,
        complex8_t ab[], int* ldab, float w[], complex8_t z[], int* ldz,
        complex8_t work[], int* lwork, float rwork[], int* lrwork,
        int iwork[], int* liwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(chbevx)(char* jobz, char* range, char* uplo, int* n,
        int* kd, complex8_t ab[], int* ldab, complex8_t q[], int* ldq,
        float* vl, float* vu, int* il, int* iu, float* abstol, int* m,
        float w[], complex8_t z[], int* ldz, complex8_t work[],
        float rwork[], int iwork[], int ifail[], int* info,
        int len_jobz, int len_range, int len_uplo);

extern void F77(chbgv)(char* jobz, char* uplo, int* n, int* ka, int* kb,
        complex8_t ab[], int* ldab, complex8_t bb[], int* ldbb,
        float w[], complex8_t z[], int* ldz, complex8_t work[],
        float rwork[], int* info, int len_jobz, int len_uplo);

extern void F77(checon)(char* uplo, int* n, complex8_t a[], int* lda,
        int ipiv[], float* anorm, float* rcond, complex8_t work[],
        int* info, int len_uplo);

extern void F77(cheev)(char* jobz, char* uplo, int* n, complex8_t a[],
        int* lda, float w[], complex8_t work[], int* lwork,
        float rwork[], int* info, int len_jobz, int len_uplo);

extern void F77(cheevd)(char* jobz, char* uplo, int* n, complex8_t a[],
        int* lda, float w[], complex8_t work[], int* lwork,
        float rwork[], int* lrwork, int iwork[], int* liwork, int* info,
        int len_jobz, int len_uplo);

extern void F77(cheevx)(char* jobz, char* range, char* uplo, int* n,
        complex8_t a[], int* lda, float* vl, float* vu, int* il,
        int* iu, float* abstol, int* m, float w[], complex8_t z[],
        int* ldz, complex8_t work[], int* lwork, float rwork[],
        int iwork[], int ifail[], int* info, int len_jobz,
        int len_range, int len_uplo);

extern void F77(chegv)(int* itype, char* jobz, char* uplo, int* n,
        complex8_t a[], int* lda, complex8_t b[], int* ldb, float w[],
        complex8_t work[], int* lwork, float rwork[], int* info,
        int len_jobz, int len_uplo);

extern void F77(chesv)(char* uplo, int* n, int* nrhs, complex8_t a[],
        int* lda, int ipiv[], complex8_t b[], int* ldb,
        complex8_t work[], int* lwork, int* info, int len_uplo);

extern void F77(chesvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex8_t a[], int* lda, complex8_t af[], int* ldaf,
        int ipiv[], complex8_t b[], int* ldb, complex8_t x[], int* ldx,
        float* rcond, float ferr[], float berr[], complex8_t work[],
        int* lwork, float rwork[], int* info, int len_fact,
        int len_uplo);

extern void F77(chetrf)(char* uplo, int* n, complex8_t a[], int* lda,
        int ipiv[], complex8_t work[], int* lwork, int* info,
        int len_uplo);

extern void F77(chetri)(char* uplo, int* n, complex8_t a[], int* lda,
        int ipiv[], complex8_t work[], int* info, int len_uplo);

extern void F77(chetrs)(char* uplo, int* n, int* nrhs, complex8_t a[],
        int* lda, int ipiv[], complex8_t b[], int* ldb, int* info,
        int len_uplo);

extern void F77(chpcon)(char* uplo, int* n, complex8_t ap[], int ipiv[],
        float* anorm, float* rcond, complex8_t work[], int* info,
        int len_uplo);

extern void F77(chpev)(char* jobz, char* uplo, int* n, complex8_t ap[],
        float w[], complex8_t z[], int* ldz, complex8_t work[],
        float rwork[], int* info, int len_jobz, int len_uplo);

extern void F77(chpevd)(char* jobz, char* uplo, int* n, complex8_t ap[],
        float w[], complex8_t z[], int* ldz, complex8_t work[],
        int* lwork, float rwork[], int* lrwork, int iwork[],
        int* liwork, int* info, int len_jobz, int len_uplo);

extern void F77(chpevx)(char* jobz, char* range, char* uplo, int* n,
        complex8_t ap[], float* vl, float* vu, int* il, int* iu,
        float* abstol, int* m, float w[], complex8_t z[], int* ldz,
        complex8_t work[], float rwork[], int iwork[], int ifail[],
        int* info, int len_jobz, int len_range, int len_uplo);

extern void F77(chpgv)(int* itype, char* jobz, char* uplo, int* n,
        complex8_t ap[], complex8_t bp[], float w[], complex8_t z[],
        int* ldz, complex8_t work[], float rwork[], int* info,
        int len_jobz, int len_uplo);

extern void F77(chpsv)(char* uplo, int* n, int* nrhs, complex8_t ap[],
        int ipiv[], complex8_t b[], int* ldb, int* info, int len_uplo);

extern void F77(chpsvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex8_t ap[], complex8_t afp[], int ipiv[], complex8_t b[],
        int* ldb, complex8_t x[], int* ldx, float* rcond, float ferr[],
        float berr[], complex8_t work[], float rwork[], int* info,
        int len_fact, int len_uplo);

extern void F77(chptrf)(char* uplo, int* n, complex8_t ap[], int ipiv[],
        int* info, int len_uplo);

extern void F77(chptri)(char* uplo, int* n, complex8_t ap[], int ipiv[],
        complex8_t work[], int* info, int len_uplo);

extern void F77(chptrs)(char* uplo, int* n, int* nrhs, complex8_t ap[],
        int ipiv[], complex8_t b[], int* ldb, int* info, int len_uplo);

extern float F77(clangb)(char* norm, int* n, int* kl, int* ku,
        complex8_t ab[], int* ldab, float rwork[], int len_norm);

extern float F77(clange)(char* norm, int* m, int* n, complex8_t a[],
        int* lda, float rwork[], int len_norm);

extern float F77(clangt)(char* norm, int* n, complex8_t dl[],
        complex8_t d[], complex8_t du[], int len_norm);

extern float F77(clanhb)(char* norm, char* uplo, int* n, int* kd,
        complex8_t ab[], int* ldab, float rwork[], int len_norm,
        int len_uplo);

extern float F77(clanhe)(char* norm, char* uplo, int* n, complex8_t a[],
        int* lda, float rwork[], int len_norm, int len_uplo);

extern float F77(clanhp)(char* norm, char* uplo, int* n, complex8_t ap[],
        float rwork[], int len_norm, int len_uplo);

extern float F77(clanht)(char* norm, int* n, float d[], complex8_t e[],
        int len_norm);

extern float F77(clansb)(char* norm, char* uplo, int* n, int* kd,
        complex8_t ab[], int* ldab, float rwork[], int len_norm,
        int len_uplo);

extern float F77(clansp)(char* norm, char* uplo, int* n, complex8_t ap[],
        float rwork[], int len_norm, int len_uplo);

extern float F77(clansy)(char* norm, char* uplo, int* n, complex8_t a[],
        int* lda, float rwork[], int len_norm, int len_uplo);

extern void F77(cpbcon)(char* uplo, int* n, int* kd, complex8_t ab[],
        int* ldab, float* anorm, float* rcond, complex8_t work[],
        float rwork[], int* info, int len_uplo);

extern void F77(cpbsv)(char* uplo, int* n, int* kd, int* nrhs,
        complex8_t ab[], int* ldab, complex8_t b[], int* ldb, int* info,
        int len_uplo);

extern void F77(cpbsvx)(char* fact, char* uplo, int* n, int* kd, int* nrhs,
        complex8_t ab[], int* ldab, complex8_t afb[], int* ldafb,
        char* equed, float s[], complex8_t b[], int* ldb,
        complex8_t x[], int* ldx, float* rcond, float ferr[],
        float berr[], complex8_t work[], float rwork[], int* info,
        int len_fact, int len_uplo, int len_equed);

extern void F77(cpbtrf)(char* uplo, int* n, int* kd, complex8_t ab[],
        int* ldab, int* info, int len_uplo);

extern void F77(cpbtrs)(char* uplo, int* n, int* kd, int* nrhs,
        complex8_t ab[], int* ldab, complex8_t b[], int* ldb, int* info,
        int len_uplo);

extern void F77(cpocon)(char* uplo, int* n, complex8_t a[], int* lda,
        float* anorm, float* rcond, complex8_t work[], float rwork[],
        int* info, int len_uplo);

extern void F77(cposv)(char* uplo, int* n, int* nrhs, complex8_t a[],
        int* lda, complex8_t b[], int* ldb, int* info, int len_uplo);

extern void F77(cposvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex8_t a[], int* lda, complex8_t af[], int* ldaf,
        char* equed, float s[], complex8_t b[], int* ldb,
        complex8_t x[], int* ldx, float* rcond, float ferr[],
        float berr[], complex8_t work[], float rwork[], int* info,
        int len_fact, int len_uplo, int len_equed);

extern void F77(cpotrf)(char* uplo, int* n, complex8_t a[], int* lda,
        int* info, int len_uplo);

extern void F77(cpotri)(char* uplo, int* n, complex8_t a[], int* lda,
        int* info, int len_uplo);

extern void F77(cpotrs)(char* uplo, int* n, int* nrhs, complex8_t a[],
        int* lda, complex8_t b[], int* ldb, int* info, int len_uplo);

extern void F77(cppcon)(char* uplo, int* n, complex8_t ap[], float* anorm,
        float* rcond, complex8_t work[], float rwork[], int* info,
        int len_uplo);

extern void F77(cppsv)(char* uplo, int* n, int* nrhs, complex8_t ap[],
        complex8_t b[], int* ldb, int* info, int len_uplo);

extern void F77(cppsvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex8_t ap[], complex8_t afp[], char* equed, float s[],
        complex8_t b[], int* ldb, complex8_t x[], int* ldx,
        float* rcond, float ferr[], float berr[], complex8_t work[],
        float rwork[], int* info, int len_fact, int len_uplo,
        int len_equed);

extern void F77(cpptrf)(char* uplo, int* n, complex8_t ap[], int* info,
        int len_uplo);

extern void F77(cpptri)(char* uplo, int* n, complex8_t ap[], int* info,
        int len_uplo);

extern void F77(cpptrs)(char* uplo, int* n, int* nrhs, complex8_t ap[],
        complex8_t b[], int* ldb, int* info, int len_uplo);

extern void F77(cptcon)(int* n, float d[], complex8_t e[], float* anorm,
        float* rcond, float rwork[], int* info);

extern void F77(cptsv)(int* n, int* nrhs, float d[], complex8_t e[],
        complex8_t b[], int* ldb, int* info);

extern void F77(cptsvx)(char* fact, int* n, int* nrhs, float d[],
        complex8_t e[], float df[], complex8_t ef[], complex8_t b[],
        int* ldb, complex8_t x[], int* ldx, float* rcond, float ferr[],
        float berr[], complex8_t work[], float rwork[], int* info,
        int len_fact);

extern void F77(cpttrf)(int* n, float d[], complex8_t e[], int* info);

extern void F77(cpttrs)(char* uplo, int* n, int* nrhs, float d[],
        complex8_t e[], complex8_t b[], int* ldb, int* info,
        int len_uplo);

extern float F77(cputime)(float* tzero);

extern void F77(cspcon)(char* uplo, int* n, complex8_t ap[], int ipiv[],
        float* anorm, float* rcond, complex8_t work[], int* info,
        int len_uplo);

extern void F77(cspsv)(char* uplo, int* n, int* nrhs, complex8_t ap[],
        int ipiv[], complex8_t b[], int* ldb, int* info, int len_uplo);

extern void F77(cspsvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex8_t ap[], complex8_t afp[], int ipiv[], complex8_t b[],
        int* ldb, complex8_t x[], int* ldx, float* rcond, float ferr[],
        float berr[], complex8_t work[], float rwork[], int* info,
        int len_fact, int len_uplo);

extern void F77(csptrf)(char* uplo, int* n, complex8_t ap[], int ipiv[],
        int* info, int len_uplo);

extern void F77(csptri)(char* uplo, int* n, complex8_t ap[], int ipiv[],
        complex8_t work[], int* info, int len_uplo);

extern void F77(csptrs)(char* uplo, int* n, int* nrhs, complex8_t ap[],
        int ipiv[], complex8_t b[], int* ldb, int* info, int len_uplo);

extern void F77(cstedc)(char* jobz, int* n, float d[], float e[],
        complex8_t z[], int* ldz, complex8_t work[], int* lwork,
        float rwork[], int* lrwork, int iwork[], int* liwork, int* info,
        int len_jobz);

extern void F77(csteqr)(char* jobz, int* n, float d[], float e[],
        complex8_t z[], int* ldz, float work[], int* info,
        int len_jobz);

extern void F77(csycon)(char* uplo, int* n, complex8_t a[], int* lda,
        int ipiv[], float* anorm, float* rcond, complex8_t work[],
        int* info, int len_uplo);

extern void F77(csysv)(char* uplo, int* n, int* nrhs, complex8_t a[],
        int* lda, int ipiv[], complex8_t b[], int* ldb,
        complex8_t work[], int* lwork, int* info, int len_uplo);

extern void F77(csysvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex8_t a[], int* lda, complex8_t af[], int* ldaf,
        int ipiv[], complex8_t b[], int* ldb, complex8_t x[], int* ldx,
        float* rcond, float ferr[], float berr[], complex8_t work[],
        int* lwork, float rwork[], int* info, int len_fact,
        int len_uplo);

extern void F77(csytrf)(char* uplo, int* n, complex8_t a[], int* lda,
        int ipiv[], complex8_t work[], int* lwork, int* info,
        int len_uplo);

extern void F77(csytri)(char* uplo, int* n, complex8_t a[], int* lda,
        int ipiv[], complex8_t work[], int* info, int len_uplo);

extern void F77(csytrs)(char* uplo, int* n, int* nrhs, complex8_t a[],
        int* lda, int ipiv[], complex8_t b[], int* ldb, int* info,
        int len_uplo);

extern void F77(ctbcon)(char* norm, char* uplo, char* diag, int* n, int* kd,
        complex8_t ab[], int* ldab, float* rcond, complex8_t work[],
        float rwork[], int* info, int len_norm, int len_uplo,
        int len_diag);

extern void F77(ctbtrs)(char* uplo, char* trans, char* diag, int* n,
        int* kd, int* nrhs, complex8_t ab[], int* ldab, complex8_t b[],
        int* ldb, int* info, int len_uplo, int len_trans, int len_diag);

extern void F77(ctpcon)(char* norm, char* uplo, char* diag, int* n,
        complex8_t ap[], float* rcond, complex8_t work[], float rwork[],
        int* info, int len_norm, int len_uplo, int len_diag);

extern void F77(ctptri)(char* uplo, char* diag, int* n, complex8_t ap[],
        int* info, int len_uplo, int len_diag);

extern void F77(ctptrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, complex8_t ap[], complex8_t b[], int* ldb, int* info,
        int len_uplo, int len_trans, int len_diag);

extern void F77(ctrcon)(char* norm, char* uplo, char* diag, int* n,
        complex8_t a[], int* lda, float* rcond, complex8_t work[],
        float rwork[], int* info, int len_norm, int len_uplo,
        int len_diag);

extern void F77(ctrtri)(char* uplo, char* diag, int* n, complex8_t a[],
        int* lda, int* info, int len_uplo, int len_diag);

extern void F77(ctrtrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, complex8_t a[], int* lda, complex8_t b[], int* ldb,
        int* info, int len_uplo, int len_trans, int len_diag);

extern void F77(ctzrqf)(int* m, int* n, complex8_t a[], int* lda,
        complex8_t tau[], int* info);

extern void F77(cunglq)(int* m, int* n, int* k, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cungql)(int* m, int* n, int* k, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cungqr)(int* m, int* n, int* k, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cungrq)(int* m, int* n, int* k, complex8_t a[], int* lda,
        complex8_t tau[], complex8_t work[], int* lwork, int* info);

extern void F77(cunmlq)(char* side, char* trans, int* m, int* n, int* k,
        complex8_t a[], int* lda, complex8_t tau[], complex8_t c[],
        int* ldc, complex8_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(cunmql)(char* side, char* trans, int* m, int* n, int* k,
        complex8_t a[], int* lda, complex8_t tau[], complex8_t c[],
        int* ldc, complex8_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(cunmqr)(char* side, char* trans, int* m, int* n, int* k,
        complex8_t a[], int* lda, complex8_t tau[], complex8_t c[],
        int* ldc, complex8_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(cunmrq)(char* side, char* trans, int* m, int* n, int* k,
        complex8_t a[], int* lda, complex8_t tau[], complex8_t c[],
        int* ldc, complex8_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(dgbcon)(char* norm, int* n, int* kl, int* ku, double ab[],
        int* ldab, int ipiv[], double* anorm, double* rcond,
        double work[], int iwork[], int* info, int len_norm);

extern void F77(dgbsv)(int* n, int* kl, int* ku, int* nrhs, double ab[],
        int* ldab, int ipiv[], double b[], int* ldb, int* info);

extern void F77(dgbsvx)(char* fact, char* trans, int* n, int* kl, int* ku,
        int* nrhs, double ab[], int* ldab, double afb[], int* ldafb,
        int ipiv[], char* equed, double r[], double c[], double b[],
        int* ldb, double x[], int* ldx, double* rcond, double ferr[],
        double berr[], double work[], int iwork[], int* info,
        int len_fact, int len_trans, int len_equed);

extern void F77(dgbtrf)(int* m, int* n, int* kl, int* ku, double ab[],
        int* ldab, int ipiv[], int* info);

extern void F77(dgbtrs)(char* trans, int* n, int* kl, int* ku, int* nrhs,
        double ab[], int* ldab, int ipiv[], double b[], int* ldb,
        int* info, int len_trans);

extern void F77(dgecon)(char* norm, int* n, double a[], int* lda,
        double* anorm, double* rcond, double work[], int iwork[],
        int* info, int len_norm);

extern void F77(dgees)(char* jobvs, char* sort, int (*select)(), int* n,
        double a[], int* lda, int* sdim, double wr[], double wi[],
        double vs[], int* ldvs, double work[], int* lwork, int bwork[],
        int* info, int len_jobvs, int len_sort);

extern void F77(dgeesx)(char* jobvs, char* sort, int (*select)(),
        char* sense, int* n, double a[], int* lda, int* sdim,
        double wr[], double wi[], double vs[], int* ldvs,
        double* rconde, double* rcondv, double work[], int* lwork,
        int iwork[], int* liwork, int bwork[], int* info, int len_jobvs,
        int len_sort, int len_sense);

extern void F77(dgeev)(char* jobvl, char* jobvr, int* n, double a[],
        int* lda, double wr[], double wi[], double vl[], int* ldvl,
        double vr[], int* ldvr, double work[], int* lwork, int* info,
        int len_jobvl, int len_jobvr);

extern void F77(dgeevx)(char* balanc, char* jobvl, char* jobvr, char* sense,
        int* n, double a[], int* lda, double wr[], double wi[],
        double vl[], int* ldvl, double vr[], int* ldvr, int* ilo,
        int* ihi, double scale[], double* abnrm, double rconde[],
        double rcondv[], double work[], int* lwork, int iwork[],
        int* info, int len_balanc, int len_jobvl, int len_jobvr,
        int len_sense);

extern void F77(dgegs)(char* jobvsl, char* jobvsr, int* n, double a[],
        int* lda, double b[], int* ldb, double alphar[],
        double alphai[], double beta[], double vsl[], int* ldvsl,
        double vsr[], int* ldvsr, double work[], int* lwork, int* info,
        int len_jobvsl, int len_jobvsr);

extern void F77(dgegv)(char* jobvl, char* jobvr, int* n, double a[],
        int* lda, double b[], int* ldb, double alphar[],
        double alphai[], double beta[], double vl[], int* ldvl,
        double vr[], int* ldvr, double work[], int* lwork, int* info,
        int len_jobvl, int len_jobvr);

extern void F77(dgelqf)(int* m, int* n, double a[], int* lda, double tau[],
        double work[], int* lwork, int* info);

extern void F77(dgels)(char* trans, int* m, int* n, int* nrhs, double a[],
        int* lda, double b[], int* ldb, double work[], int* lwork,
        int* info, int len_trans);

extern void F77(dgelss)(int* m, int* n, int* nrhs, double a[], int* lda,
        double b[], int* ldb, double s[], double* rcond, int* rank,
        double work[], int* lwork, int* info);

extern void F77(dgelsx)(int* m, int* n, int* nrhs, double a[], int* lda,
        double b[], int* ldb, int jpvt[], double* rcond, int* rank,
        double work[], int* info);

extern void F77(dgeqlf)(int* m, int* n, double a[], int* lda, double tau[],
        double work[], int* lwork, int* info);

extern void F77(dgeqpf)(int* m, int* n, double a[], int* lda, int jpvt[],
        double tau[], double work[], int* info);

extern void F77(dgeqrf)(int* m, int* n, double a[], int* lda, double tau[],
        double work[], int* lwork, int* info);

extern void F77(dgerqf)(int* m, int* n, double a[], int* lda, double tau[],
        double work[], int* lwork, int* info);

extern void F77(dgesv)(int* n, int* nrhs, double a[], int* lda, int ipiv[],
        double b[], int* ldb, int* info);

extern void F77(dgesvd)(char* jobu, char* jobvt, int* m, int* n, double a[],
        int* lda, double s[], double u[], int* ldu, double vt[],
        int* ldvt, double work[], int* lwork, int* info, int len_jobu,
        int len_jobvt);

extern void F77(dgesvx)(char* fact, char* trans, int* n, int* nrhs,
        double a[], int* lda, double af[], int* ldaf, int ipiv[],
        char* equed, double r[], double c[], double b[], int* ldb,
        double x[], int* ldx, double* rcond, double ferr[],
        double berr[], double work[], int iwork[], int* info,
        int len_fact, int len_trans, int len_equed);

extern void F77(dgetrf)(int* m, int* n, double a[], int* lda, int ipiv[],
        int* info);

extern void F77(dgetri)(int* n, double a[], int* lda, int ipiv[],
        double work[], int* lwork, int* info);

extern void F77(dgetrs)(char* trans, int* n, int* nrhs, double a[],
        int* lda, int ipiv[], double b[], int* ldb, int* info,
        int len_trans);

extern void F77(dggglm)(int* n, int* m, int* p, double a[], int* lda,
        double b[], int* ldb, double d[], double x[], double y[],
        double work[], int* lwork, int* info);

extern void F77(dgglse)(int* m, int* n, int* p, double a[], int* lda,
        double b[], int* ldb, double c[], double d[], double x[],
        double work[], int* lwork, int* info);

extern void F77(dggqrf)(int* m, int* n, int* p, double a[], int* lda,
        double taua[], double b[], int* ldb, double taub[],
        double work[], int* lwork, int* info);

extern void F77(dggrqf)(int* m, int* n, int* p, double a[], int* lda,
        double taua[], double b[], int* ldb, double taub[],
        double work[], int* lwork, int* info);

extern void F77(dggsvd)(char* jobu, char* jobv, char* jobq, int* m, int* n,
        int* p, int* k, int* l, double a[], int* lda, double b[],
        int* ldb, double alpha[], double beta[], double u[], int* ldu,
        double v[], int* ldv, double q[], int* ldq, double work[],
        int iwork[], int* info, int len_jobu, int len_jobv,
        int len_jobq);

extern void F77(dgtcon)(char* norm, int* n, double dl[], double d[],
        double du[], double du2[], int ipiv[], double* anorm,
        double* rcond, double work[], int iwork[], int* info,
        int len_norm);

extern void F77(dgtsv)(int* n, int* nrhs, double dl[], double d[],
        double du[], double b[], int* ldb, int* info);

extern void F77(dgtsvx)(char* fact, char* trans, int* n, int* nrhs,
        double dl[], double d[], double du[], double dlf[], double df[],
        double duf[], double du2[], int ipiv[], double b[], int* ldb,
        double x[], int* ldx, double* rcond, double ferr[],
        double berr[], double work[], int iwork[], int* info,
        int len_fact, int len_trans);

extern void F77(dgttrf)(int* n, double dl[], double d[], double du[],
        double du2[], int ipiv[], int* info);

extern void F77(dgttrs)(char* trans, int* n, int* nrhs, double dl[],
        double d[], double du[], double du2[], int ipiv[], double b[],
        int* ldb, int* info, int len_trans);

extern double F77(dlamch)(char* name, int len_name);

extern double F77(dlangb)(char* norm, int* n, int* kl, int* ku, double ab[],
        int* ldab, double work[], int len_norm);

extern double F77(dlange)(char* norm, int* m, int* n, double a[], int* lda,
        double work[], int len_norm);

extern double F77(dlangt)(char* norm, int* n, double dl[], double d[],
        double du[], int len_norm);

extern double F77(dlansb)(char* norm, char* uplo, int* n, int* kd,
        double ab[], int* ldab, double work[], int len_norm,
        int len_uplo);

extern double F77(dlansp)(char* norm, char* uplo, int* n, double ap[],
        double work[], int len_norm, int len_uplo);

extern double F77(dlanst)(char* norm, int* n, double d[], double e[],
        int len_norm);

extern double F77(dlansy)(char* norm, char* uplo, int* n, double a[],
        int* lda, double work[], int len_norm, int len_uplo);

extern void F77(dorglq)(int* m, int* n, int* k, double a[], int* lda,
        double tau[], double work[], int* lwork, int* info);

extern void F77(dorgql)(int* m, int* n, int* k, double a[], int* lda,
        double tau[], double work[], int* lwork, int* info);

extern void F77(dorgqr)(int* m, int* n, int* k, double a[], int* lda,
        double tau[], double work[], int* lwork, int* info);

extern void F77(dorgrq)(int* m, int* n, int* k, double a[], int* lda,
        double tau[], double work[], int* lwork, int* info);

extern void F77(dormlq)(char* side, char* trans, int* m, int* n, int* k,
        double a[], int* lda, double tau[], double c[], int* ldc,
        double work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(dormql)(char* side, char* trans, int* m, int* n, int* k,
        double a[], int* lda, double tau[], double c[], int* ldc,
        double work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(dormqr)(char* side, char* trans, int* m, int* n, int* k,
        double a[], int* lda, double tau[], double c[], int* ldc,
        double work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(dormrq)(char* side, char* trans, int* m, int* n, int* k,
        double a[], int* lda, double tau[], double c[], int* ldc,
        double work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(dpbcon)(char* uplo, int* n, int* kd, double ab[], int* ldab,
        double* anorm, double* rcond, double work[], int iwork[],
        int* info, int len_uplo);

extern void F77(dpbsv)(char* uplo, int* n, int* kd, int* nrhs, double ab[],
        int* ldab, double b[], int* ldb, int* info, int len_uplo);

extern void F77(dpbsvx)(char* fact, char* uplo, int* n, int* kd, int* nrhs,
        double ab[], int* ldab, double afb[], int* ldafb, char* equed,
        double s[], double b[], int* ldb, double x[], int* ldx,
        double* rcond, double ferr[], double berr[], double work[],
        int iwork[], int* info, int len_fact, int len_uplo,
        int len_equed);

extern void F77(dpbtrf)(char* uplo, int* n, int* kd, double ab[], int* ldab,
        int* info, int len_uplo);

extern void F77(dpbtrs)(char* uplo, int* n, int* kd, int* nrhs, double ab[],
        int* ldab, double b[], int* ldb, int* info, int len_uplo);

extern void F77(dpocon)(char* uplo, int* n, double a[], int* lda,
        double* anorm, double* rcond, double work[], int iwork[],
        int* info, int len_uplo);

extern void F77(dposv)(char* uplo, int* n, int* nrhs, double a[], int* lda,
        double b[], int* ldb, int* info, int len_uplo);

extern void F77(dposvx)(char* fact, char* uplo, int* n, int* nrhs,
        double a[], int* lda, double af[], int* ldaf, char* equed,
        double s[], double b[], int* ldb, double x[], int* ldx,
        double* rcond, double ferr[], double berr[], double work[],
        int iwork[], int* info, int len_fact, int len_uplo,
        int len_equed);

extern void F77(dpotrf)(char* uplo, int* n, double a[], int* lda, int* info,
        int len_uplo);

extern void F77(dpotri)(char* uplo, int* n, double a[], int* lda, int* info,
        int len_uplo);

extern void F77(dpotrs)(char* uplo, int* n, int* nrhs, double a[], int* lda,
        double b[], int* ldb, int* info, int len_uplo);

extern void F77(dppcon)(char* uplo, int* n, double ap[], double* anorm,
        double* rcond, double work[], int iwork[], int* info,
        int len_uplo);

extern void F77(dppsv)(char* uplo, int* n, int* nrhs, double ap[],
        double b[], int* ldb, int* info, int len_uplo);

extern void F77(dppsvx)(char* fact, char* uplo, int* n, int* nrhs,
        double ap[], double afp[], char* equed, double s[], double b[],
        int* ldb, double x[], int* ldx, double* rcond, double ferr[],
        double berr[], double work[], int iwork[], int* info,
        int len_fact, int len_uplo, int len_equed);

extern void F77(dpptrf)(char* uplo, int* n, double ap[], int* info,
        int len_uplo);

extern void F77(dpptri)(char* uplo, int* n, double ap[], int* info,
        int len_uplo);

extern void F77(dpptrs)(char* uplo, int* n, int* nrhs, double ap[],
        double b[], int* ldb, int* info, int len_uplo);

extern void F77(dptcon)(int* n, double d[], double e[], double* anorm,
        double* rcond, double work[], int* info);

extern void F77(dptsv)(int* n, int* nrhs, double d[], double e[],
        double b[], int* ldb, int* info);

extern void F77(dptsvx)(char* fact, int* n, int* nrhs, double d[],
        double e[], double df[], double ef[], double b[], int* ldb,
        double x[], int* ldx, double* rcond, double ferr[],
        double berr[], double work[], int* info, int len_fact);

extern void F77(dpttrf)(int* n, double d[], double e[], int* info);

extern void F77(dpttrs)(int* n, int* nrhs, double d[], double e[],
        double b[], int* ldb, int* info);

extern void F77(dsbev)(char* jobz, char* uplo, int* n, int* kd, double ab[],
        int* ldab, double w[], double z[], int* ldz, double work[],
        int* info, int len_jobz, int len_uplo);

extern void F77(dsbevd)(char* jobz, char* uplo, int* n, int* kd,
        double ab[], int* ldab, double w[], double z[], int* ldz,
        double work[], int* lwork, int iwork[], int* liwork, int* info,
        int len_jobz, int len_uplo);

extern void F77(dsbevx)(char* jobz, char* range, char* uplo, int* n,
        int* kd, double ab[], int* ldab, double q[], int* ldq,
        double* vl, double* vu, int* il, int* iu, double* abstol,
        int* m, double w[], double z[], int* ldz, double work[],
        int iwork[], int ifail[], int* info, int len_jobz,
        int len_range, int len_uplo);

extern void F77(dsbgv)(char* jobz, char* uplo, int* n, int* ka, int* kb,
        double ab[], int* ldab, double bb[], int* ldbb, double w[],
        double z[], int* ldz, double work[], int* info, int len_jobz,
        int len_uplo);

extern void F77(dspcon)(char* uplo, int* n, double ap[], int ipiv[],
        double* anorm, double* rcond, double work[], int iwork[],
        int* info, int len_uplo);

extern void F77(dspev)(char* jobz, char* uplo, int* n, double ap[],
        double w[], double z[], int* ldz, double work[], int* info,
        int len_jobz, int len_uplo);

extern void F77(dspevd)(char* jobz, char* uplo, int* n, double ap[],
        double w[], double z[], int* ldz, double work[], int* lwork,
        int iwork[], int* liwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(dspevx)(char* jobz, char* range, char* uplo, int* n,
        double ap[], double* vl, double* vu, int* il, int* iu,
        double* abstol, int* m, double w[], double z[], int* ldz,
        double work[], int iwork[], int ifail[], int* info,
        int len_jobz, int len_range, int len_uplo);

extern void F77(dspgv)(int* itype, char* jobz, char* uplo, int* n,
        double ap[], double bp[], double w[], double z[], int* ldz,
        double work[], int* info, int len_jobz, int len_uplo);

extern void F77(dspsv)(char* uplo, int* n, int* nrhs, double ap[],
        int ipiv[], double b[], int* ldb, int* info, int len_uplo);

extern void F77(dspsvx)(char* fact, char* uplo, int* n, int* nrhs,
        double ap[], double afp[], int ipiv[], double b[], int* ldb,
        double x[], int* ldx, double* rcond, double ferr[],
        double berr[], double work[], int iwork[], int* info,
        int len_fact, int len_uplo);

extern void F77(dsptrf)(char* uplo, int* n, double ap[], int ipiv[],
        int* info, int len_uplo);

extern void F77(dsptri)(char* uplo, int* n, double ap[], int ipiv[],
        double work[], int* info, int len_uplo);

extern void F77(dsptrs)(char* uplo, int* n, int* nrhs, double ap[],
        int ipiv[], double b[], int* ldb, int* info, int len_uplo);

extern void F77(dstedc)(char* jobz, int* n, double d[], double e[],
        double z[], int* ldz, double work[], int* lwork, int iwork[],
        int* liwork, int* info, int len_jobz);

extern void F77(dsteqr)(char* jobz, int* n, double d[], double e[],
        double z[], int* ldz, double work[], int* info, int len_jobz);

extern void F77(dstev)(char* jobz, int* n, double d[], double e[],
        double z[], int* ldz, double work[], int* info, int len_jobz);

extern void F77(dstevd)(char* jobz, int* n, double d[], double e[],
        double z[], int* ldz, double work[], int* lwork, int iwork[],
        int* liwork, int* info, int len_jobz);

extern void F77(dstevx)(char* jobz, char* range, int* n, double d[],
        double e[], double* vl, double* vu, int* il, int* iu,
        double* abstol, int* m, double w[], double z[], int* ldz,
        double work[], int iwork[], int ifail[], int* info,
        int len_jobz, int len_range);

extern void F77(dsycon)(char* uplo, int* n, double a[], int* lda,
        int ipiv[], double* anorm, double* rcond, double work[],
        int iwork[], int* info, int len_uplo);

extern void F77(dsyev)(char* jobz, char* uplo, int* n, double a[], int* lda,
        double w[], double work[], int* lwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(dsyevd)(char* jobz, char* uplo, int* n, double a[],
        int* lda, double w[], double work[], int* lwork, int iwork[],
        int* liwork, int* info, int len_jobz, int len_uplo);

extern void F77(dsyevx)(char* jobz, char* range, char* uplo, int* n,
        double a[], int* lda, double* vl, double* vu, int* il, int* iu,
        double* abstol, int* m, double w[], double z[], int* ldz,
        double work[], int* lwork, int iwork[], int ifail[], int* info,
        int len_jobz, int len_range, int len_uplo);

extern void F77(dsygv)(int* itype, char* jobz, char* uplo, int* n,
        double a[], int* lda, double b[], int* ldb, double w[],
        double work[], int* lwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(dsysv)(char* uplo, int* n, int* nrhs, double a[], int* lda,
        int ipiv[], double b[], int* ldb, double work[], int* lwork,
        int* info, int len_uplo);

extern void F77(dsysvx)(char* fact, char* uplo, int* n, int* nrhs,
        double a[], int* lda, double af[], int* ldaf, int ipiv[],
        double b[], int* ldb, double x[], int* ldx, double* rcond,
        double ferr[], double berr[], double work[], int* lwork,
        int iwork[], int* info, int len_fact, int len_uplo);

extern void F77(dsytrf)(char* uplo, int* n, double a[], int* lda,
        int ipiv[], double work[], int* lwork, int* info, int len_uplo);

extern void F77(dsytri)(char* uplo, int* n, double a[], int* lda,
        int ipiv[], double work[], int* info, int len_uplo);

extern void F77(dsytrs)(char* uplo, int* n, int* nrhs, double a[], int* lda,
        int ipiv[], double b[], int* ldb, int* info, int len_uplo);

extern void F77(dtbcon)(char* norm, char* uplo, char* diag, int* n, int* kd,
        double ab[], int* ldab, double* rcond, double work[],
        int iwork[], int* info, int len_norm, int len_uplo,
        int len_diag);

extern void F77(dtbtrs)(char* uplo, char* trans, char* diag, int* n,
        int* kd, int* nrhs, double ab[], int* ldab, double b[],
        int* ldb, int* info, int len_uplo, int len_trans, int len_diag);

extern void F77(dtpcon)(char* norm, char* uplo, char* diag, int* n,
        double ap[], double* rcond, double work[], int iwork[],
        int* info, int len_norm, int len_uplo, int len_diag);

extern void F77(dtptri)(char* uplo, char* diag, int* n, double ap[],
        int* info, int len_uplo, int len_diag);

extern void F77(dtptrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, double ap[], double b[], int* ldb, int* info,
        int len_uplo, int len_trans, int len_diag);

extern void F77(dtrcon)(char* norm, char* uplo, char* diag, int* n,
        double a[], int* lda, double* rcond, double work[], int iwork[],
        int* info, int len_norm, int len_uplo, int len_diag);

extern void F77(dtrtri)(char* uplo, char* diag, int* n, double a[],
        int* lda, int* info, int len_uplo, int len_diag);

extern void F77(dtrtrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, double a[], int* lda, double b[], int* ldb,
        int* info, int len_uplo, int len_trans, int len_diag);

extern void F77(dtzrqf)(int* m, int* n, double a[], int* lda, double tau[],
        int* info);

extern int F77(ilaenv)(int* ispec, char* name, char* opts, int* n1, int* n2,
        int* n3, int* n4, int len_name, int len_opts);

extern void F77(sgbcon)(char* norm, int* n, int* kl, int* ku, float ab[],
        int* ldab, int ipiv[], float* anorm, float* rcond, float work[],
        int iwork[], int* info, int len_norm);

extern void F77(sgbsv)(int* n, int* kl, int* ku, int* nrhs, float ab[],
        int* ldab, int ipiv[], float b[], int* ldb, int* info);

extern void F77(sgbsvx)(char* fact, char* trans, int* n, int* kl, int* ku,
        int* nrhs, float ab[], int* ldab, float afb[], int* ldafb,
        int ipiv[], char* equed, float r[], float c[], float b[],
        int* ldb, float x[], int* ldx, float* rcond, float ferr[],
        float berr[], float work[], int iwork[], int* info,
        int len_fact, int len_trans, int len_equed);

extern void F77(sgbtrf)(int* m, int* n, int* kl, int* ku, float ab[],
        int* ldab, int ipiv[], int* info);

extern void F77(sgbtrs)(char* trans, int* n, int* kl, int* ku, int* nrhs,
        float ab[], int* ldab, int ipiv[], float b[], int* ldb,
        int* info, int len_trans);

extern void F77(sgecon)(char* norm, int* n, float a[], int* lda,
        float* anorm, float* rcond, float work[], int iwork[],
        int* info, int len_norm);

extern void F77(sgees)(char* jobvs, char* sort, int (*select)(), int* n,
        float a[], int* lda, int* sdim, float wr[], float wi[],
        float vs[], int* ldvs, float work[], int* lwork, int bwork[],
        int* info, int len_jobvs, int len_sort);

extern void F77(sgeesx)(char* jobvs, char* sort, int (*select)(),
        char* sense, int* n, float a[], int* lda, int* sdim, float wr[],
        float wi[], float vs[], int* ldvs, float* rconde, float* rcondv,
        float work[], int* lwork, int iwork[], int* liwork, int bwork[],
        int* info, int len_jobvs, int len_sort, int len_sense);

extern void F77(sgeev)(char* jobvl, char* jobvr, int* n, float a[],
        int* lda, float wr[], float wi[], float vl[], int* ldvl,
        float vr[], int* ldvr, float work[], int* lwork, int* info,
        int len_jobvl, int len_jobvr);

extern void F77(sgeevx)(char* balanc, char* jobvl, char* jobvr, char* sense,
        int* n, float a[], int* lda, float wr[], float wi[], float vl[],
        int* ldvl, float vr[], int* ldvr, int* ilo, int* ihi,
        float scale[], float* abnrm, float rconde[], float rcondv[],
        float work[], int* lwork, int iwork[], int* info,
        int len_balanc, int len_jobvl, int len_jobvr, int len_sense);

extern void F77(sgegs)(char* jobvsl, char* jobvsr, int* n, float a[],
        int* lda, float b[], int* ldb, float alphar[], float alphai[],
        float beta[], float vsl[], int* ldvsl, float vsr[], int* ldvsr,
        float work[], int* lwork, int* info, int len_jobvsl,
        int len_jobvsr);

extern void F77(sgegv)(char* jobvl, char* jobvr, int* n, float a[],
        int* lda, float b[], int* ldb, float alphar[], float alphai[],
        float beta[], float vl[], int* ldvl, float vr[], int* ldvr,
        float work[], int* lwork, int* info, int len_jobvl,
        int len_jobvr);

extern void F77(sgelqf)(int* m, int* n, float a[], int* lda, float tau[],
        float work[], int* lwork, int* info);

extern void F77(sgels)(char* trans, int* m, int* n, int* nrhs, float a[],
        int* lda, float b[], int* ldb, float work[], int* lwork,
        int* info, int len_trans);

extern void F77(sgelss)(int* m, int* n, int* nrhs, float a[], int* lda,
        float b[], int* ldb, float s[], float* rcond, int* rank,
        float work[], int* lwork, int* info);

extern void F77(sgelsx)(int* m, int* n, int* nrhs, float a[], int* lda,
        float b[], int* ldb, int jpvt[], float* rcond, int* rank,
        float work[], int* info);

extern void F77(sgeqlf)(int* m, int* n, float a[], int* lda, float tau[],
        float work[], int* lwork, int* info);

extern void F77(sgeqpf)(int* m, int* n, float a[], int* lda, int jpvt[],
        float tau[], float work[], int* info);

extern void F77(sgeqrf)(int* m, int* n, float a[], int* lda, float tau[],
        float work[], int* lwork, int* info);

extern void F77(sgerqf)(int* m, int* n, float a[], int* lda, float tau[],
        float work[], int* lwork, int* info);

extern void F77(sgesv)(int* n, int* nrhs, float a[], int* lda, int ipiv[],
        float b[], int* ldb, int* info);

extern void F77(sgesvd)(char* jobu, char* jobvt, int* m, int* n, float a[],
        int* lda, float s[], float u[], int* ldu, float vt[], int* ldvt,
        float work[], int* lwork, int* info, int len_jobu,
        int len_jobvt);

extern void F77(sgesvx)(char* fact, char* trans, int* n, int* nrhs,
        float a[], int* lda, float af[], int* ldaf, int ipiv[],
        char* equed, float r[], float c[], float b[], int* ldb,
        float x[], int* ldx, float* rcond, float ferr[], float berr[],
        float work[], int iwork[], int* info, int len_fact,
        int len_trans, int len_equed);

extern void F77(sgetrf)(int* m, int* n, float a[], int* lda, int ipiv[],
        int* info);

extern void F77(sgetri)(int* n, float a[], int* lda, int ipiv[],
        float work[], int* lwork, int* info);

extern void F77(sgetrs)(char* trans, int* n, int* nrhs, float a[], int* lda,
        int ipiv[], float b[], int* ldb, int* info, int len_trans);

extern void F77(sggglm)(int* n, int* m, int* p, float a[], int* lda,
        float b[], int* ldb, float d[], float x[], float y[],
        float work[], int* lwork, int* info);

extern void F77(sgglse)(int* m, int* n, int* p, float a[], int* lda,
        float b[], int* ldb, float c[], float d[], float x[],
        float work[], int* lwork, int* info);

extern void F77(sggqrf)(int* m, int* n, int* p, float a[], int* lda,
        float taua[], float b[], int* ldb, float taub[], float work[],
        int* lwork, int* info);

extern void F77(sggrqf)(int* m, int* n, int* p, float a[], int* lda,
        float taua[], float b[], int* ldb, float taub[], float work[],
        int* lwork, int* info);

extern void F77(sggsvd)(char* jobu, char* jobv, char* jobq, int* m, int* n,
        int* p, int* k, int* l, float a[], int* lda, float b[],
        int* ldb, float alpha[], float beta[], float u[], int* ldu,
        float v[], int* ldv, float q[], int* ldq, float work[],
        int iwork[], int* info, int len_jobu, int len_jobv,
        int len_jobq);

extern void F77(sgtcon)(char* norm, int* n, float dl[], float d[],
        float du[], float du2[], int ipiv[], float* anorm, float* rcond,
        float work[], int iwork[], int* info, int len_norm);

extern void F77(sgtsv)(int* n, int* nrhs, float dl[], float d[], float du[],
        float b[], int* ldb, int* info);

extern void F77(sgtsvx)(char* fact, char* trans, int* n, int* nrhs,
        float dl[], float d[], float du[], float dlf[], float df[],
        float duf[], float du2[], int ipiv[], float b[], int* ldb,
        float x[], int* ldx, float* rcond, float ferr[], float berr[],
        float work[], int iwork[], int* info, int len_fact,
        int len_trans);

extern void F77(sgttrf)(int* n, float dl[], float d[], float du[],
        float du2[], int ipiv[], int* info);

extern void F77(sgttrs)(char* trans, int* n, int* nrhs, float dl[],
        float d[], float du[], float du2[], int ipiv[], float b[],
        int* ldb, int* info, int len_trans);

extern float F77(slamch)(char* name, int len_name);

extern float F77(slangb)(char* norm, int* n, int* kl, int* ku, float ab[],
        int* ldab, float work[], int len_norm);

extern float F77(slange)(char* norm, int* m, int* n, float a[], int* lda,
        float work[], int len_norm);

extern float F77(slangt)(char* norm, int* n, float dl[], float d[],
        float du[], int len_norm);

extern float F77(slansb)(char* norm, char* uplo, int* n, int* kd,
        float ab[], int* ldab, float work[], int len_norm,
        int len_uplo);

extern float F77(slansp)(char* norm, char* uplo, int* n, float ap[],
        float work[], int len_norm, int len_uplo);

extern float F77(slanst)(char* norm, int* n, float d[], float e[],
        int len_norm);

extern float F77(slansy)(char* norm, char* uplo, int* n, float a[],
        int* lda, float work[], int len_norm, int len_uplo);

extern void F77(sorglq)(int* m, int* n, int* k, float a[], int* lda,
        float tau[], float work[], int* lwork, int* info);

extern void F77(sorgql)(int* m, int* n, int* k, float a[], int* lda,
        float tau[], float work[], int* lwork, int* info);

extern void F77(sorgqr)(int* m, int* n, int* k, float a[], int* lda,
        float tau[], float work[], int* lwork, int* info);

extern void F77(sorgrq)(int* m, int* n, int* k, float a[], int* lda,
        float tau[], float work[], int* lwork, int* info);

extern void F77(sormlq)(char* side, char* trans, int* m, int* n, int* k,
        float a[], int* lda, float tau[], float c[], int* ldc,
        float work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(sormql)(char* side, char* trans, int* m, int* n, int* k,
        float a[], int* lda, float tau[], float c[], int* ldc,
        float work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(sormqr)(char* side, char* trans, int* m, int* n, int* k,
        float a[], int* lda, float tau[], float c[], int* ldc,
        float work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(sormrq)(char* side, char* trans, int* m, int* n, int* k,
        float a[], int* lda, float tau[], float c[], int* ldc,
        float work[], int* lwork, int* info, int len_side,
        int len_trans);

extern void F77(spbcon)(char* uplo, int* n, int* kd, float ab[], int* ldab,
        float* anorm, float* rcond, float work[], int iwork[],
        int* info, int len_uplo);

extern void F77(spbsv)(char* uplo, int* n, int* kd, int* nrhs, float ab[],
        int* ldab, float b[], int* ldb, int* info, int len_uplo);

extern void F77(spbsvx)(char* fact, char* uplo, int* n, int* kd, int* nrhs,
        float ab[], int* ldab, float afb[], int* ldafb, char* equed,
        float s[], float b[], int* ldb, float x[], int* ldx,
        float* rcond, float ferr[], float berr[], float work[],
        int iwork[], int* info, int len_fact, int len_uplo,
        int len_equed);

extern void F77(spbtrf)(char* uplo, int* n, int* kd, float ab[], int* ldab,
        int* info, int len_uplo);

extern void F77(spbtrs)(char* uplo, int* n, int* kd, int* nrhs, float ab[],
        int* ldab, float b[], int* ldb, int* info, int len_uplo);

extern void F77(spocon)(char* uplo, int* n, float a[], int* lda,
        float* anorm, float* rcond, float work[], int iwork[],
        int* info, int len_uplo);

extern void F77(sposv)(char* uplo, int* n, int* nrhs, float a[], int* lda,
        float b[], int* ldb, int* info, int len_uplo);

extern void F77(sposvx)(char* fact, char* uplo, int* n, int* nrhs,
        float a[], int* lda, float af[], int* ldaf, char* equed,
        float s[], float b[], int* ldb, float x[], int* ldx,
        float* rcond, float ferr[], float berr[], float work[],
        int iwork[], int* info, int len_fact, int len_uplo,
        int len_equed);

extern void F77(spotrf)(char* uplo, int* n, float a[], int* lda, int* info,
        int len_uplo);

extern void F77(spotri)(char* uplo, int* n, float a[], int* lda, int* info,
        int len_uplo);

extern void F77(spotrs)(char* uplo, int* n, int* nrhs, float a[], int* lda,
        float b[], int* ldb, int* info, int len_uplo);

extern void F77(sppcon)(char* uplo, int* n, float ap[], float* anorm,
        float* rcond, float work[], int iwork[], int* info,
        int len_uplo);

extern void F77(sppsv)(char* uplo, int* n, int* nrhs, float ap[], float b[],
        int* ldb, int* info, int len_uplo);

extern void F77(sppsvx)(char* fact, char* uplo, int* n, int* nrhs,
        float ap[], float afp[], char* equed, float s[], float b[],
        int* ldb, float x[], int* ldx, float* rcond, float ferr[],
        float berr[], float work[], int iwork[], int* info,
        int len_fact, int len_uplo, int len_equed);

extern void F77(spptrf)(char* uplo, int* n, float ap[], int* info,
        int len_uplo);

extern void F77(spptri)(char* uplo, int* n, float ap[], int* info,
        int len_uplo);

extern void F77(spptrs)(char* uplo, int* n, int* nrhs, float ap[],
        float b[], int* ldb, int* info, int len_uplo);

extern void F77(sptcon)(int* n, float d[], float e[], float* anorm,
        float* rcond, float work[], int* info);

extern void F77(sptsv)(int* n, int* nrhs, float d[], float e[], float b[],
        int* ldb, int* info);

extern void F77(sptsvx)(char* fact, int* n, int* nrhs, float d[], float e[],
        float df[], float ef[], float b[], int* ldb, float x[],
        int* ldx, float* rcond, float ferr[], float berr[],
        float work[], int* info, int len_fact);

extern void F77(spttrf)(int* n, float d[], float e[], int* info);

extern void F77(spttrs)(int* n, int* nrhs, float d[], float e[], float b[],
        int* ldb, int* info);

extern void F77(ssbev)(char* jobz, char* uplo, int* n, int* kd, float ab[],
        int* ldab, float w[], float z[], int* ldz, float work[],
        int* info, int len_jobz, int len_uplo);

extern void F77(ssbevd)(char* jobz, char* uplo, int* n, int* kd, float ab[],
        int* ldab, float w[], float z[], int* ldz, float work[],
        int* lwork, int iwork[], int* liwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(ssbevx)(char* jobz, char* range, char* uplo, int* n,
        int* kd, float ab[], int* ldab, float q[], int* ldq, float* vl,
        float* vu, int* il, int* iu, float* abstol, int* m, float w[],
        float z[], int* ldz, float work[], int iwork[], int ifail[],
        int* info, int len_jobz, int len_range, int len_uplo);

extern void F77(ssbgv)(char* jobz, char* uplo, int* n, int* ka, int* kb,
        float ab[], int* ldab, float bb[], int* ldbb, float w[],
        float z[], int* ldz, float work[], int* info, int len_jobz,
        int len_uplo);

extern void F77(sspcon)(char* uplo, int* n, float ap[], int ipiv[],
        float* anorm, float* rcond, float work[], int iwork[],
        int* info, int len_uplo);

extern void F77(sspev)(char* jobz, char* uplo, int* n, float ap[],
        float w[], float z[], int* ldz, float work[], int* info,
        int len_jobz, int len_uplo);

extern void F77(sspevd)(char* jobz, char* uplo, int* n, float ap[],
        float w[], float z[], int* ldz, float work[], int* lwork,
        int iwork[], int* liwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(sspevx)(char* jobz, char* range, char* uplo, int* n,
        float ap[], float* vl, float* vu, int* il, int* iu,
        float* abstol, int* m, float w[], float z[], int* ldz,
        float work[], int iwork[], int ifail[], int* info, int len_jobz,
        int len_range, int len_uplo);

extern void F77(sspgv)(int* itype, char* jobz, char* uplo, int* n,
        float ap[], float bp[], float w[], float z[], int* ldz,
        float work[], int* info, int len_jobz, int len_uplo);

extern void F77(sspsv)(char* uplo, int* n, int* nrhs, float ap[],
        int ipiv[], float b[], int* ldb, int* info, int len_uplo);

extern void F77(sspsvx)(char* fact, char* uplo, int* n, int* nrhs,
        float ap[], float afp[], int ipiv[], float b[], int* ldb,
        float x[], int* ldx, float* rcond, float ferr[], float berr[],
        float work[], int iwork[], int* info, int len_fact,
        int len_uplo);

extern void F77(ssptrf)(char* uplo, int* n, float ap[], int ipiv[],
        int* info, int len_uplo);

extern void F77(ssptri)(char* uplo, int* n, float ap[], int ipiv[],
        float work[], int* info, int len_uplo);

extern void F77(ssptrs)(char* uplo, int* n, int* nrhs, float ap[],
        int ipiv[], float b[], int* ldb, int* info, int len_uplo);

extern void F77(sstedc)(char* jobz, int* n, float d[], float e[], float z[],
        int* ldz, float work[], int* lwork, int iwork[], int* liwork,
        int* info, int len_jobz);

extern void F77(ssteqr)(char* jobz, int* n, float d[], float e[], float z[],
        int* ldz, float work[], int* info, int len_jobz);

extern void F77(sstev)(char* jobz, int* n, float d[], float e[], float z[],
        int* ldz, float work[], int* info, int len_jobz);

extern void F77(sstevd)(char* jobz, int* n, float d[], float e[], float z[],
        int* ldz, float work[], int* lwork, int iwork[], int* liwork,
        int* info, int len_jobz);

extern void F77(sstevx)(char* jobz, char* range, int* n, float d[],
        float e[], float* vl, float* vu, int* il, int* iu,
        float* abstol, int* m, float w[], float z[], int* ldz,
        float work[], int iwork[], int ifail[], int* info, int len_jobz,
        int len_range);

extern void F77(ssycon)(char* uplo, int* n, float a[], int* lda, int ipiv[],
        float* anorm, float* rcond, float work[], int iwork[],
        int* info, int len_uplo);

extern void F77(ssyev)(char* jobz, char* uplo, int* n, float a[], int* lda,
        float w[], float work[], int* lwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(ssyevd)(char* jobz, char* uplo, int* n, float a[], int* lda,
        float w[], float work[], int* lwork, int iwork[], int* liwork,
        int* info, int len_jobz, int len_uplo);

extern void F77(ssyevx)(char* jobz, char* range, char* uplo, int* n,
        float a[], int* lda, float* vl, float* vu, int* il, int* iu,
        float* abstol, int* m, float w[], float z[], int* ldz,
        float work[], int* lwork, int iwork[], int ifail[], int* info,
        int len_jobz, int len_range, int len_uplo);

extern void F77(ssygv)(int* itype, char* jobz, char* uplo, int* n,
        float a[], int* lda, float b[], int* ldb, float w[],
        float work[], int* lwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(ssysv)(char* uplo, int* n, int* nrhs, float a[], int* lda,
        int ipiv[], float b[], int* ldb, float work[], int* lwork,
        int* info, int len_uplo);

extern void F77(ssysvx)(char* fact, char* uplo, int* n, int* nrhs,
        float a[], int* lda, float af[], int* ldaf, int ipiv[],
        float b[], int* ldb, float x[], int* ldx, float* rcond,
        float ferr[], float berr[], float work[], int* lwork,
        int iwork[], int* info, int len_fact, int len_uplo);

extern void F77(ssytrf)(char* uplo, int* n, float a[], int* lda, int ipiv[],
        float work[], int* lwork, int* info, int len_uplo);

extern void F77(ssytri)(char* uplo, int* n, float a[], int* lda, int ipiv[],
        float work[], int* info, int len_uplo);

extern void F77(ssytrs)(char* uplo, int* n, int* nrhs, float a[], int* lda,
        int ipiv[], float b[], int* ldb, int* info, int len_uplo);

extern void F77(stbcon)(char* norm, char* uplo, char* diag, int* n, int* kd,
        float ab[], int* ldab, float* rcond, float work[], int iwork[],
        int* info, int len_norm, int len_uplo, int len_diag);

extern void F77(stbtrs)(char* uplo, char* trans, char* diag, int* n,
        int* kd, int* nrhs, float ab[], int* ldab, float b[], int* ldb,
        int* info, int len_uplo, int len_trans, int len_diag);

extern void F77(stpcon)(char* norm, char* uplo, char* diag, int* n,
        float ap[], float* rcond, float work[], int iwork[], int* info,
        int len_norm, int len_uplo, int len_diag);

extern void F77(stptri)(char* uplo, char* diag, int* n, float ap[],
        int* info, int len_uplo, int len_diag);

extern void F77(stptrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, float ap[], float b[], int* ldb, int* info,
        int len_uplo, int len_trans, int len_diag);

extern void F77(strcon)(char* norm, char* uplo, char* diag, int* n,
        float a[], int* lda, float* rcond, float work[], int iwork[],
        int* info, int len_norm, int len_uplo, int len_diag);

extern void F77(strtri)(char* uplo, char* diag, int* n, float a[], int* lda,
        int* info, int len_uplo, int len_diag);

extern void F77(strtrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, float a[], int* lda, float b[], int* ldb, int* info,
        int len_uplo, int len_trans, int len_diag);

extern void F77(stzrqf)(int* m, int* n, float a[], int* lda, float tau[],
        int* info);

extern float F77(walltime)(float* tzero);

extern void F77(xerbla)(char* name, int* iarg, int len_name);

extern void F77(zgbcon)(char* norm, int* n, int* kl, int* ku,
        complex16_t ab[], int* ldab, int ipiv[], double* anorm,
        double* rcond, complex16_t work[], double rwork[], int* info,
        int len_norm);

extern void F77(zgbsv)(int* n, int* kl, int* ku, int* nrhs,
        complex16_t ab[], int* ldab, int ipiv[], complex16_t b[],
        int* ldb, int* info);

extern void F77(zgbsvx)(char* fact, char* trans, int* n, int* kl, int* ku,
        int* nrhs, complex16_t ab[], int* ldab, complex16_t afb[],
        int* ldafb, int ipiv[], char* equed, double r[], double c[],
        complex16_t b[], int* ldb, complex16_t x[], int* ldx,
        double* rcond, double ferr[], double berr[], complex16_t work[],
        double rwork[], int* info, int len_fact, int len_trans,
        int len_equed);

extern void F77(zgbtrf)(int* m, int* n, int* kl, int* ku, complex16_t ab[],
        int* ldab, int ipiv[], int* info);

extern void F77(zgbtrs)(char* trans, int* n, int* kl, int* ku, int* nrhs,
        complex16_t ab[], int* ldab, int ipiv[], complex16_t b[],
        int* ldb, int* info, int len_trans);

extern void F77(zgecon)(char* norm, int* n, complex16_t a[], int* lda,
        double* anorm, double* rcond, complex16_t work[],
        double rwork[], int* info, int len_norm);

extern void F77(zgees)(char* jobvs, char* sort, int (*select)(), int* n,
        complex16_t a[], int* lda, int* sdim, complex16_t w[],
        complex16_t vs[], int* ldvs, complex16_t work[], int* lwork,
        double rwork[], int bwork[], int* info, int len_jobvs,
        int len_sort);

extern void F77(zgeesx)(char* jobvs, char* sort, int (*select)(),
        char* sense, int* n, complex16_t a[], int* lda, int* sdim,
        complex16_t w[], complex16_t vs[], int* ldvs, double* rconde,
        double* rcondv, complex16_t work[], int* lwork, double rwork[],
        int bwork[], int* info, int len_jobvs, int len_sort,
        int len_sense);

extern void F77(zgeev)(char* jobvl, char* jobvr, int* n, complex16_t a[],
        int* lda, complex16_t w[], complex16_t vl[], int* ldvl,
        complex16_t vr[], int* ldvr, complex16_t work[], int* lwork,
        double rwork[], int* info, int len_jobvl, int len_jobvr);

extern void F77(zgeevx)(char* balanc, char* jobvl, char* jobvr, char* sense,
        int* n, complex16_t a[], int* lda, complex16_t w[],
        complex16_t vl[], int* ldvl, complex16_t vr[], int* ldvr,
        int* ilo, int* ihi, double scale[], double* abnrm,
        double rconde[], double rcondv[], complex16_t work[],
        int* lwork, double rwork[], int* info, int len_balanc,
        int len_jobvl, int len_jobvr, int len_sense);

extern void F77(zgegs)(char* jobvsl, char* jobvsr, int* n, complex16_t a[],
        int* lda, complex16_t b[], int* ldb, complex16_t alpha[],
        complex16_t beta[], complex16_t vsl[], int* ldvsl,
        complex16_t vsr[], int* ldvsr, complex16_t work[], int* lwork,
        double rwork[], int* info, int len_jobvsl, int len_jobvsr);

extern void F77(zgegv)(char* jobvl, char* jobvr, int* n, complex16_t a[],
        int* lda, complex16_t b[], int* ldb, complex16_t alpha[],
        complex16_t beta[], complex16_t vl[], int* ldvl,
        complex16_t vr[], int* ldvr, complex16_t work[], int* lwork,
        double rwork[], int* info, int len_jobvl, int len_jobvr);

extern void F77(zgelqf)(int* m, int* n, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zgels)(char* trans, int* m, int* n, int* nrhs,
        complex16_t a[], int* lda, complex16_t b[], int* ldb,
        complex16_t work[], int* lwork, int* info, int len_trans);

extern void F77(zgelss)(int* m, int* n, int* nrhs, complex16_t a[],
        int* lda, complex16_t b[], int* ldb, double s[], double* rcond,
        int* rank, complex16_t work[], int* lwork, double rwork[],
        int* info);

extern void F77(zgelsx)(int* m, int* n, int* nrhs, complex16_t a[],
        int* lda, complex16_t b[], int* ldb, int jpvt[], double* rcond,
        int* rank, complex16_t work[], double rwork[], int* info);

extern void F77(zgeqlf)(int* m, int* n, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zgeqpf)(int* m, int* n, complex16_t a[], int* lda,
        int jpvt[], complex16_t tau[], complex16_t work[],
        double rwork[], int* info);

extern void F77(zgeqrf)(int* m, int* n, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zgerqf)(int* m, int* n, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zgesv)(int* n, int* nrhs, complex16_t a[], int* lda,
        int ipiv[], complex16_t b[], int* ldb, int* info);

extern void F77(zgesvd)(char* jobu, char* jobvt, int* m, int* n,
        complex16_t a[], int* lda, double s[], complex16_t u[],
        int* ldu, complex16_t vt[], int* ldvt, complex16_t work[],
        int* lwork, double rwork[], int* info, int len_jobu,
        int len_jobvt);

extern void F77(zgesvx)(char* fact, char* trans, int* n, int* nrhs,
        complex16_t a[], int* lda, complex16_t af[], int* ldaf,
        int ipiv[], char* equed, double r[], double c[],
        complex16_t b[], int* ldb, complex16_t x[], int* ldx,
        double* rcond, double ferr[], double berr[], complex16_t work[],
        double rwork[], int* info, int len_fact, int len_trans,
        int len_equed);

extern void F77(zgetrf)(int* m, int* n, complex16_t a[], int* lda,
        int ipiv[], int* info);

extern void F77(zgetri)(int* n, complex16_t a[], int* lda, int ipiv[],
        complex16_t work[], int* lwork, int* info);

extern void F77(zgetrs)(char* trans, int* n, int* nrhs, complex16_t a[],
        int* lda, int ipiv[], complex16_t b[], int* ldb, int* info,
        int len_trans);

extern void F77(zggglm)(int* n, int* m, int* p, complex16_t a[], int* lda,
        complex16_t b[], int* ldb, complex16_t d[], complex16_t x[],
        complex16_t y[], complex16_t work[], int* lwork, int* info);

extern void F77(zgglse)(int* m, int* n, int* p, complex16_t a[], int* lda,
        complex16_t b[], int* ldb, complex16_t c[], complex16_t d[],
        complex16_t x[], complex16_t work[], int* lwork, int* info);

extern void F77(zggqrf)(int* m, int* n, int* p, complex16_t a[], int* lda,
        complex16_t taua[], complex16_t b[], int* ldb,
        complex16_t taub[], complex16_t work[], int* lwork, int* info);

extern void F77(zggrqf)(int* m, int* n, int* p, complex16_t a[], int* lda,
        complex16_t taua[], complex16_t b[], int* ldb,
        complex16_t taub[], complex16_t work[], int* lwork, int* info);

extern void F77(zggsvd)(char* jobu, char* jobv, char* jobq, int* m, int* n,
        int* p, int* k, int* l, complex16_t a[], int* lda,
        complex16_t b[], int* ldb, double alpha[], double beta[],
        complex16_t u[], int* ldu, complex16_t v[], int* ldv,
        complex16_t q[], int* ldq, complex16_t work[], double rwork[],
        int iwork[], int* info, int len_jobu, int len_jobv,
        int len_jobq);

extern void F77(zgtcon)(char* norm, int* n, complex16_t dl[],
        complex16_t d[], complex16_t du[], complex16_t du2[],
        int ipiv[], double* anorm, double* rcond, complex16_t work[],
        int* info, int len_norm);

extern void F77(zgtsv)(int* n, int* nrhs, complex16_t dl[], complex16_t d[],
        complex16_t du[], complex16_t b[], int* ldb, int* info);

extern void F77(zgtsvx)(char* fact, char* trans, int* n, int* nrhs,
        complex16_t dl[], complex16_t d[], complex16_t du[],
        complex16_t dlf[], complex16_t df[], complex16_t duf[],
        complex16_t du2[], int ipiv[], complex16_t b[], int* ldb,
        complex16_t x[], int* ldx, double* rcond, double ferr[],
        double berr[], complex16_t work[], double rwork[], int* info,
        int len_fact, int len_trans);

extern void F77(zgttrf)(int* n, complex16_t dl[], complex16_t d[],
        complex16_t du[], complex16_t du2[], int ipiv[], int* info);

extern void F77(zgttrs)(char* trans, int* n, int* nrhs, complex16_t dl[],
        complex16_t d[], complex16_t du[], complex16_t du2[],
        int ipiv[], complex16_t b[], int* ldb, int* info,
        int len_trans);

extern void F77(zhbev)(char* jobz, char* uplo, int* n, int* kd,
        complex16_t ab[], int* ldab, double w[], complex16_t z[],
        int* ldz, complex16_t work[], double rwork[], int* info,
        int len_jobz, int len_uplo);

extern void F77(zhbevd)(char* jobz, char* uplo, int* n, int* kd,
        complex16_t ab[], int* ldab, double w[], complex16_t z[],
        int* ldz, complex16_t work[], int* lwork, double rwork[],
        int* lrwork, int iwork[], int* liwork, int* info, int len_jobz,
        int len_uplo);

extern void F77(zhbevx)(char* jobz, char* range, char* uplo, int* n,
        int* kd, complex16_t ab[], int* ldab, complex16_t q[], int* ldq,
        double* vl, double* vu, int* il, int* iu, double* abstol,
        int* m, double w[], complex16_t z[], int* ldz,
        complex16_t work[], double rwork[], int iwork[], int ifail[],
        int* info, int len_jobz, int len_range, int len_uplo);

extern void F77(zhbgv)(char* jobz, char* uplo, int* n, int* ka, int* kb,
        complex16_t ab[], int* ldab, complex16_t bb[], int* ldbb,
        double w[], complex16_t z[], int* ldz, complex16_t work[],
        double rwork[], int* info, int len_jobz, int len_uplo);

extern void F77(zhecon)(char* uplo, int* n, complex16_t a[], int* lda,
        int ipiv[], double* anorm, double* rcond, complex16_t work[],
        int* info, int len_uplo);

extern void F77(zheev)(char* jobz, char* uplo, int* n, complex16_t a[],
        int* lda, double w[], complex16_t work[], int* lwork,
        double rwork[], int* info, int len_jobz, int len_uplo);

extern void F77(zheevd)(char* jobz, char* uplo, int* n, complex16_t a[],
        int* lda, double w[], complex16_t work[], int* lwork,
        double rwork[], int* lrwork, int iwork[], int* liwork,
        int* info, int len_jobz, int len_uplo);

extern void F77(zheevx)(char* jobz, char* range, char* uplo, int* n,
        complex16_t a[], int* lda, double* vl, double* vu, int* il,
        int* iu, double* abstol, int* m, double w[], complex16_t z[],
        int* ldz, complex16_t work[], int* lwork, double rwork[],
        int iwork[], int ifail[], int* info, int len_jobz,
        int len_range, int len_uplo);

extern void F77(zhegv)(int* itype, char* jobz, char* uplo, int* n,
        complex16_t a[], int* lda, complex16_t b[], int* ldb,
        double w[], complex16_t work[], int* lwork, double rwork[],
        int* info, int len_jobz, int len_uplo);

extern void F77(zhesv)(char* uplo, int* n, int* nrhs, complex16_t a[],
        int* lda, int ipiv[], complex16_t b[], int* ldb,
        complex16_t work[], int* lwork, int* info, int len_uplo);

extern void F77(zhesvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex16_t a[], int* lda, complex16_t af[], int* ldaf,
        int ipiv[], complex16_t b[], int* ldb, complex16_t x[],
        int* ldx, double* rcond, double ferr[], double berr[],
        complex16_t work[], int* lwork, double rwork[], int* info,
        int len_fact, int len_uplo);

extern void F77(zhetrf)(char* uplo, int* n, complex16_t a[], int* lda,
        int ipiv[], complex16_t work[], int* lwork, int* info,
        int len_uplo);

extern void F77(zhetri)(char* uplo, int* n, complex16_t a[], int* lda,
        int ipiv[], complex16_t work[], int* info, int len_uplo);

extern void F77(zhetrs)(char* uplo, int* n, int* nrhs, complex16_t a[],
        int* lda, int ipiv[], complex16_t b[], int* ldb, int* info,
        int len_uplo);

extern void F77(zhpcon)(char* uplo, int* n, complex16_t ap[], int ipiv[],
        double* anorm, double* rcond, complex16_t work[], int* info,
        int len_uplo);

extern void F77(zhpev)(char* jobz, char* uplo, int* n, complex16_t ap[],
        double w[], complex16_t z[], int* ldz, complex16_t work[],
        double rwork[], int* info, int len_jobz, int len_uplo);

extern void F77(zhpevd)(char* jobz, char* uplo, int* n, complex16_t ap[],
        double w[], complex16_t z[], int* ldz, complex16_t work[],
        int* lwork, double rwork[], int* lrwork, int iwork[],
        int* liwork, int* info, int len_jobz, int len_uplo);

extern void F77(zhpevx)(char* jobz, char* range, char* uplo, int* n,
        complex16_t ap[], double* vl, double* vu, int* il, int* iu,
        double* abstol, int* m, double w[], complex16_t z[], int* ldz,
        complex16_t work[], double rwork[], int iwork[], int ifail[],
        int* info, int len_jobz, int len_range, int len_uplo);

extern void F77(zhpgv)(int* itype, char* jobz, char* uplo, int* n,
        complex16_t ap[], complex16_t bp[], double w[], complex16_t z[],
        int* ldz, complex16_t work[], double rwork[], int* info,
        int len_jobz, int len_uplo);

extern void F77(zhpsv)(char* uplo, int* n, int* nrhs, complex16_t ap[],
        int ipiv[], complex16_t b[], int* ldb, int* info, int len_uplo);

extern void F77(zhpsvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex16_t ap[], complex16_t afp[], int ipiv[],
        complex16_t b[], int* ldb, complex16_t x[], int* ldx,
        double* rcond, double ferr[], double berr[], complex16_t work[],
        double rwork[], int* info, int len_fact, int len_uplo);

extern void F77(zhptrf)(char* uplo, int* n, complex16_t ap[], int ipiv[],
        int* info, int len_uplo);

extern void F77(zhptri)(char* uplo, int* n, complex16_t ap[], int ipiv[],
        complex16_t work[], int* info, int len_uplo);

extern void F77(zhptrs)(char* uplo, int* n, int* nrhs, complex16_t ap[],
        int ipiv[], complex16_t b[], int* ldb, int* info, int len_uplo);

extern double F77(zlangb)(char* norm, int* n, int* kl, int* ku,
        complex16_t ab[], int* ldab, double rwork[], int len_norm);

extern double F77(zlange)(char* norm, int* m, int* n, complex16_t a[],
        int* lda, double rwork[], int len_norm);

extern double F77(zlangt)(char* norm, int* n, complex16_t dl[],
        complex16_t d[], complex16_t du[], int len_norm);

extern double F77(zlanhb)(char* norm, char* uplo, int* n, int* kd,
        complex16_t ab[], int* ldab, double rwork[], int len_norm,
        int len_uplo);

extern double F77(zlanhe)(char* norm, char* uplo, int* n, complex16_t a[],
        int* lda, double rwork[], int len_norm, int len_uplo);

extern double F77(zlanhp)(char* norm, char* uplo, int* n, complex16_t ap[],
        double rwork[], int len_norm, int len_uplo);

extern double F77(zlanht)(char* norm, int* n, double d[], complex16_t e[],
        int len_norm);

extern double F77(zlansb)(char* norm, char* uplo, int* n, int* kd,
        complex16_t ab[], int* ldab, double rwork[], int len_norm,
        int len_uplo);

extern double F77(zlansp)(char* norm, char* uplo, int* n, complex16_t ap[],
        double rwork[], int len_norm, int len_uplo);

extern double F77(zlansy)(char* norm, char* uplo, int* n, complex16_t a[],
        int* lda, double rwork[], int len_norm, int len_uplo);

extern void F77(zpbcon)(char* uplo, int* n, int* kd, complex16_t ab[],
        int* ldab, double* anorm, double* rcond, complex16_t work[],
        double rwork[], int* info, int len_uplo);

extern void F77(zpbsv)(char* uplo, int* n, int* kd, int* nrhs,
        complex16_t ab[], int* ldab, complex16_t b[], int* ldb,
        int* info, int len_uplo);

extern void F77(zpbsvx)(char* fact, char* uplo, int* n, int* kd, int* nrhs,
        complex16_t ab[], int* ldab, complex16_t afb[], int* ldafb,
        char* equed, double s[], complex16_t b[], int* ldb,
        complex16_t x[], int* ldx, double* rcond, double ferr[],
        double berr[], complex16_t work[], double rwork[], int* info,
        int len_fact, int len_uplo, int len_equed);

extern void F77(zpbtrf)(char* uplo, int* n, int* kd, complex16_t ab[],
        int* ldab, int* info, int len_uplo);

extern void F77(zpbtrs)(char* uplo, int* n, int* kd, int* nrhs,
        complex16_t ab[], int* ldab, complex16_t b[], int* ldb,
        int* info, int len_uplo);

extern void F77(zpocon)(char* uplo, int* n, complex16_t a[], int* lda,
        double* anorm, double* rcond, complex16_t work[],
        double rwork[], int* info, int len_uplo);

extern void F77(zposv)(char* uplo, int* n, int* nrhs, complex16_t a[],
        int* lda, complex16_t b[], int* ldb, int* info, int len_uplo);

extern void F77(zposvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex16_t a[], int* lda, complex16_t af[], int* ldaf,
        char* equed, double s[], complex16_t b[], int* ldb,
        complex16_t x[], int* ldx, double* rcond, double ferr[],
        double berr[], complex16_t work[], double rwork[], int* info,
        int len_fact, int len_uplo, int len_equed);

extern void F77(zpotrf)(char* uplo, int* n, complex16_t a[], int* lda,
        int* info, int len_uplo);

extern void F77(zpotri)(char* uplo, int* n, complex16_t a[], int* lda,
        int* info, int len_uplo);

extern void F77(zpotrs)(char* uplo, int* n, int* nrhs, complex16_t a[],
        int* lda, complex16_t b[], int* ldb, int* info, int len_uplo);

extern void F77(zppcon)(char* uplo, int* n, complex16_t ap[], double* anorm,
        double* rcond, complex16_t work[], double rwork[], int* info,
        int len_uplo);

extern void F77(zppsv)(char* uplo, int* n, int* nrhs, complex16_t ap[],
        complex16_t b[], int* ldb, int* info, int len_uplo);

extern void F77(zppsvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex16_t ap[], complex16_t afp[], char* equed, double s[],
        complex16_t b[], int* ldb, complex16_t x[], int* ldx,
        double* rcond, double ferr[], double berr[], complex16_t work[],
        double rwork[], int* info, int len_fact, int len_uplo,
        int len_equed);

extern void F77(zpptrf)(char* uplo, int* n, complex16_t ap[], int* info,
        int len_uplo);

extern void F77(zpptri)(char* uplo, int* n, complex16_t ap[], int* info,
        int len_uplo);

extern void F77(zpptrs)(char* uplo, int* n, int* nrhs, complex16_t ap[],
        complex16_t b[], int* ldb, int* info, int len_uplo);

extern void F77(zptcon)(int* n, double d[], complex16_t e[], double* anorm,
        double* rcond, double rwork[], int* info);

extern void F77(zptsv)(int* n, int* nrhs, double d[], complex16_t e[],
        complex16_t b[], int* ldb, int* info);

extern void F77(zptsvx)(char* fact, int* n, int* nrhs, double d[],
        complex16_t e[], double df[], complex16_t ef[], complex16_t b[],
        int* ldb, complex16_t x[], int* ldx, double* rcond,
        double ferr[], double berr[], complex16_t work[],
        double rwork[], int* info, int len_fact);

extern void F77(zpttrf)(int* n, double d[], complex16_t e[], int* info);

extern void F77(zpttrs)(char* uplo, int* n, int* nrhs, double d[],
        complex16_t e[], complex16_t b[], int* ldb, int* info,
        int len_uplo);

extern void F77(zspcon)(char* uplo, int* n, complex16_t ap[], int ipiv[],
        double* anorm, double* rcond, complex16_t work[], int* info,
        int len_uplo);

extern void F77(zspsv)(char* uplo, int* n, int* nrhs, complex16_t ap[],
        int ipiv[], complex16_t b[], int* ldb, int* info, int len_uplo);

extern void F77(zspsvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex16_t ap[], complex16_t afp[], int ipiv[],
        complex16_t b[], int* ldb, complex16_t x[], int* ldx,
        double* rcond, double ferr[], double berr[], complex16_t work[],
        double rwork[], int* info, int len_fact, int len_uplo);

extern void F77(zsptrf)(char* uplo, int* n, complex16_t ap[], int ipiv[],
        int* info, int len_uplo);

extern void F77(zsptri)(char* uplo, int* n, complex16_t ap[], int ipiv[],
        complex16_t work[], int* info, int len_uplo);

extern void F77(zsptrs)(char* uplo, int* n, int* nrhs, complex16_t ap[],
        int ipiv[], complex16_t b[], int* ldb, int* info, int len_uplo);

extern void F77(zstedc)(char* jobz, int* n, double d[], double e[],
        complex16_t z[], int* ldz, complex16_t work[], int* lwork,
        double rwork[], int* lrwork, int iwork[], int* liwork,
        int* info, int len_jobz);

extern void F77(zsteqr)(char* jobz, int* n, double d[], double e[],
        complex16_t z[], int* ldz, double work[], int* info,
        int len_jobz);

extern void F77(zsycon)(char* uplo, int* n, complex16_t a[], int* lda,
        int ipiv[], double* anorm, double* rcond, complex16_t work[],
        int* info, int len_uplo);

extern void F77(zsysv)(char* uplo, int* n, int* nrhs, complex16_t a[],
        int* lda, int ipiv[], complex16_t b[], int* ldb,
        complex16_t work[], int* lwork, int* info, int len_uplo);

extern void F77(zsysvx)(char* fact, char* uplo, int* n, int* nrhs,
        complex16_t a[], int* lda, complex16_t af[], int* ldaf,
        int ipiv[], complex16_t b[], int* ldb, complex16_t x[],
        int* ldx, double* rcond, double ferr[], double berr[],
        complex16_t work[], int* lwork, double rwork[], int* info,
        int len_fact, int len_uplo);

extern void F77(zsytrf)(char* uplo, int* n, complex16_t a[], int* lda,
        int ipiv[], complex16_t work[], int* lwork, int* info,
        int len_uplo);

extern void F77(zsytri)(char* uplo, int* n, complex16_t a[], int* lda,
        int ipiv[], complex16_t work[], int* info, int len_uplo);

extern void F77(zsytrs)(char* uplo, int* n, int* nrhs, complex16_t a[],
        int* lda, int ipiv[], complex16_t b[], int* ldb, int* info,
        int len_uplo);

extern void F77(ztbcon)(char* norm, char* uplo, char* diag, int* n, int* kd,
        complex16_t ab[], int* ldab, double* rcond, complex16_t work[],
        double rwork[], int* info, int len_norm, int len_uplo,
        int len_diag);

extern void F77(ztbtrs)(char* uplo, char* trans, char* diag, int* n,
        int* kd, int* nrhs, complex16_t ab[], int* ldab,
        complex16_t b[], int* ldb, int* info, int len_uplo,
        int len_trans, int len_diag);

extern void F77(ztpcon)(char* norm, char* uplo, char* diag, int* n,
        complex16_t ap[], double* rcond, complex16_t work[],
        double rwork[], int* info, int len_norm, int len_uplo,
        int len_diag);

extern void F77(ztptri)(char* uplo, char* diag, int* n, complex16_t ap[],
        int* info, int len_uplo, int len_diag);

extern void F77(ztptrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, complex16_t ap[], complex16_t b[], int* ldb,
        int* info, int len_uplo, int len_trans, int len_diag);

extern void F77(ztrcon)(char* norm, char* uplo, char* diag, int* n,
        complex16_t a[], int* lda, double* rcond, complex16_t work[],
        double rwork[], int* info, int len_norm, int len_uplo,
        int len_diag);

extern void F77(ztrtri)(char* uplo, char* diag, int* n, complex16_t a[],
        int* lda, int* info, int len_uplo, int len_diag);

extern void F77(ztrtrs)(char* uplo, char* trans, char* diag, int* n,
        int* nrhs, complex16_t a[], int* lda, complex16_t b[], int* ldb,
        int* info, int len_uplo, int len_trans, int len_diag);

extern void F77(ztzrqf)(int* m, int* n, complex16_t a[], int* lda,
        complex16_t tau[], int* info);

extern void F77(zunglq)(int* m, int* n, int* k, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zungql)(int* m, int* n, int* k, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zungqr)(int* m, int* n, int* k, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zungrq)(int* m, int* n, int* k, complex16_t a[], int* lda,
        complex16_t tau[], complex16_t work[], int* lwork, int* info);

extern void F77(zunmlq)(char* side, char* trans, int* m, int* n, int* k,
        complex16_t a[], int* lda, complex16_t tau[], complex16_t c[],
        int* ldc, complex16_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(zunmql)(char* side, char* trans, int* m, int* n, int* k,
        complex16_t a[], int* lda, complex16_t tau[], complex16_t c[],
        int* ldc, complex16_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(zunmqr)(char* side, char* trans, int* m, int* n, int* k,
        complex16_t a[], int* lda, complex16_t tau[], complex16_t c[],
        int* ldc, complex16_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(zunmrq)(char* side, char* trans, int* m, int* n, int* k,
        complex16_t a[], int* lda, complex16_t tau[], complex16_t c[],
        int* ldc, complex16_t work[], int* lwork, int* info,
        int len_side, int len_trans);

extern void F77(dlacpy)(char *UPLO, int *M, int *N, double *A, int *LDA, 
                        double *B, int *LDB, int len_uplo);

extern void F77(dlaset)(char *UPLO, int *M, int *N, double *ALPHA, 
                        double *BETA, double *A, int *LDA, int len_uplo );

extern void F77(dlarnv)(int *IDIST, int *ISEED, int *N, double *X);

#endif
