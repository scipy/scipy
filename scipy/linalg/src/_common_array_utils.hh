/*
 * LAPACK declarations.
 */
#pragma once

namespace sp_linalg {

using namespace lapack;
using wrapper::f32;
using wrapper::f64;
using wrapper::c64;
using wrapper::c128;
using wrapper::real_of_t;
using wrapper::complex_of_t;
using wrapper::is_complex_v;


/*
 * Generate type overloads, to map from C array types (f32, f64, c64, c128)
 * to LAPACK prefixes, "sdcz".
 */
#define GEN_GETRF(PREFIX, TYPE) \
inline void \
call_getrf(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getrf)(m, n, a, lda, ipiv, info); \
};

GEN_GETRF(s,f32)
GEN_GETRF(d,f64)
GEN_GETRF(c,c64)
GEN_GETRF(z,c128)


#define GEN_GETRS(PREFIX, TYPE) \
inline void \
call_getrs(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getrs)(trans, n, nrhs, a, lda, ipiv, b, ldb, info); \
};

GEN_GETRS(s,f32)
GEN_GETRS(d,f64)
GEN_GETRS(c,c64)
GEN_GETRS(z,c128)


#define GEN_GETRI(PREFIX, TYPE) \
inline void \
call_getri(CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getri)(n, a, lda, ipiv, work, lwork, info); \
};

GEN_GETRI(s,f32)
GEN_GETRI(d,f64)
GEN_GETRI(c,c64)
GEN_GETRI(z,c128)


// NB: iwork for real arrays or rwork for complex arrays
#define GEN_GECON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_gecon(char* norm, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## gecon)(norm, n, a, lda, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_GECON(s, f32, f32, CBLAS_INT)
GEN_GECON(d, f64, f64, CBLAS_INT)
GEN_GECON(c, c64, f32, f32)
GEN_GECON(z, c128, f64, f64)


#define GEN_TRTRI(PREFIX, TYPE) \
inline void \
call_trtri(char* uplo, char *diag, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## trtri)(uplo, diag, n, a, lda, info); \
};

GEN_TRTRI(s, f32)
GEN_TRTRI(d, f64)
GEN_TRTRI(c, c64)
GEN_TRTRI(z, c128)


#define GEN_TRCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_trcon(char* norm, char *uplo, char *diag, CBLAS_INT *n, CTYPE *a, CBLAS_INT *lda, RTYPE *rcond, CTYPE *work, void *irwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## trcon)(norm, uplo, diag, n, a, lda, rcond, work, (WTYPE *)irwork, info); \
};

GEN_TRCON(s, f32, f32, CBLAS_INT)
GEN_TRCON(d, f64, f64, CBLAS_INT)
GEN_TRCON(c, c64, f32, f32)
GEN_TRCON(z, c128, f64, f64)


#define GEN_TRTRS(PREFIX, TYPE) \
inline void \
call_trtrs(char* uplo, char *trans, char *diag, CBLAS_INT* n, CBLAS_INT* nrhs, TYPE* a, CBLAS_INT* lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## trtrs)(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info); \
};

GEN_TRTRS(s, f32)
GEN_TRTRS(d, f64)
GEN_TRTRS(c, c64)
GEN_TRTRS(z, c128)


#define GEN_POTRF(PREFIX, TYPE) \
inline void \
call_potrf(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potrf)(uplo, n, a, lda, info); \
};

GEN_POTRF(s, f32)
GEN_POTRF(d, f64)
GEN_POTRF(c, c64)
GEN_POTRF(z, c128)


#define GEN_POTRI(PREFIX, TYPE) \
inline void \
call_potri(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potri)(uplo, n, a, lda, info); \
};

GEN_POTRI(s, f32)
GEN_POTRI(d, f64)
GEN_POTRI(c, c64)
GEN_POTRI(z, c128)


// NB: iwork for real arrays or rwork for complex arrays
#define GEN_POCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_pocon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## pocon)(uplo, n, a, lda, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_POCON(s, f32, f32, CBLAS_INT)
GEN_POCON(d, f64, f64, CBLAS_INT)
GEN_POCON(c, c64, f32, f32)
GEN_POCON(z, c128, f64, f64)


#define GEN_POTRS(PREFIX, TYPE) \
inline void \
call_potrs(char* uplo, CBLAS_INT* n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potrs)(uplo, n, nrhs, a, lda, b, ldb, info); \
};

GEN_POTRS(s, f32)
GEN_POTRS(d, f64)
GEN_POTRS(c, c64)
GEN_POTRS(z, c128)


#define GEN_SYTRF(PREFIX, TYPE) \
inline void \
call_sytrf(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sytrf)(uplo, n, a, lda, ipiv, work, lwork, info); \
};

GEN_SYTRF(s, f32)
GEN_SYTRF(d, f64)
GEN_SYTRF(c, c64)
GEN_SYTRF(z, c128)


// dispatch to sSYtrf for "float hermitian"
#define GEN_HETRF(PREFIX, L_PREFIX, TYPE) \
inline void \
call_hetrf(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## L_PREFIX ## trf)(uplo, n, a, lda, ipiv, work, lwork, info); \
};

GEN_HETRF(s, sy, f32)
GEN_HETRF(d, sy, f64)
GEN_HETRF(c, he, c64)
GEN_HETRF(z, he, c128)


#define GEN_SYTRI(PREFIX, TYPE) \
inline void \
call_sytri(char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## sytri)(uplo, n, a, lda, ipiv, work, info); \
};

GEN_SYTRI(s, f32)
GEN_SYTRI(d, f64)
GEN_SYTRI(c, c64)
GEN_SYTRI(z, c128)


// dispatch to sSYtri for "float hermitian"
#define GEN_HETRI(PREFIX, L_PREFIX, TYPE) \
inline void \
call_hetri(char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## L_PREFIX ## tri)(uplo, n, a, lda, ipiv, work, info); \
};

GEN_HETRI(s, sy, f32)
GEN_HETRI(d, sy, f64)
GEN_HETRI(c, he, c64)
GEN_HETRI(z, he, c128)


// NB: iwork for real arrays only, no rwork for complex routines (10 arguments for s- d- variants; 9 arguments for c- and z- variants)
#define GEN_SYCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_sycon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sycon)(uplo, n, a, lda, ipiv, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_SYCON(s, f32, f32, CBLAS_INT)
GEN_SYCON(d, f64, f64, CBLAS_INT)

#define GEN_SYCON_CZ(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_sycon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sycon)(uplo, n, a, lda, ipiv, anorm, rcond, work, info); \
};

GEN_SYCON_CZ(c, c64, f32, f32)
GEN_SYCON_CZ(z, c128, f64, f64)


// dispatch to sSYcon for "float hermitian"
#define GEN_HECON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_hecon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sycon)(uplo, n, a, lda, ipiv, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_HECON(s, f32, f32, CBLAS_INT)
GEN_HECON(d, f64, f64, CBLAS_INT)

#define GEN_HECON_CZ(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
call_hecon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, CBLAS_INT *ipiv, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## hecon)(uplo, n, a, lda, ipiv, anorm, rcond, work, info); \
};

GEN_HECON_CZ(c, c64, f32, f32)
GEN_HECON_CZ(z, c128, f64, f64)


#define GEN_SYTRS(PREFIX, TYPE) \
inline void \
call_sytrs(char* uplo, CBLAS_INT* n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## sytrs)(uplo, n, nrhs, a, lda, ipiv, b, ldb, info); \
};

GEN_SYTRS(s, f32)
GEN_SYTRS(d, f64)
GEN_SYTRS(c, c64)
GEN_SYTRS(z, c128)


// dispatch to sSYtrs for "float hermitian"
#define GEN_HETRS(PREFIX, L_PREFIX, TYPE) \
inline void \
call_hetrs(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## L_PREFIX ## trs)(uplo, n, nrhs, a, lda, ipiv, b, ldb, info); \
};

GEN_HETRS(s, sy, f32)
GEN_HETRS(d, sy, f64)
GEN_HETRS(c, he, c64)
GEN_HETRS(z, he, c128)


#define GEN_GTTRF(PREFIX, TYPE) \
inline void \
call_gttrf(CBLAS_INT *n, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gttrf)(n, dl, d, du, du2, ipiv, info); \
};

GEN_GTTRF(s, f32)
GEN_GTTRF(d, f64)
GEN_GTTRF(c, c64)
GEN_GTTRF(z, c128)


#define GEN_GTTRS(PREFIX, TYPE) \
inline void \
call_gttrs(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gttrs)(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info); \
};

GEN_GTTRS(s, f32)
GEN_GTTRS(d, f64)
GEN_GTTRS(c, c64)
GEN_GTTRS(z, c128)


#define GEN_GTCON(PREFIX, TYPE) \
inline void \
call_gtcon(char *norm, CBLAS_INT *n, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, TYPE *anorm, TYPE *rcond, TYPE *work, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gtcon)(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info); \
};

GEN_GTCON(s, f32)
GEN_GTCON(d, f64)


// NB: `iwork` is not used for c- and z- variants of ?gtcon
#define GEN_GTCON_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gtcon(char *norm, CBLAS_INT *n, TYPE *dl, TYPE *d, TYPE *du, TYPE *du2, CBLAS_INT *ipiv, RTYPE *anorm, RTYPE *rcond, TYPE *work, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gtcon)(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info); \
};

GEN_GTCON_CZ(c, c64, f32)
GEN_GTCON_CZ(z, c128, f64)


#define GEN_GBTRF(PREFIX, TYPE) \
inline void \
call_gbtrf(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *kl, CBLAS_INT *ku, TYPE *ab, CBLAS_INT *ldab, CBLAS_INT *ipiv, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gbtrf)(m, n, kl, ku, ab, ldab, ipiv, info); \
};

GEN_GBTRF(s, f32)
GEN_GBTRF(d, f64)
GEN_GBTRF(c, c64)
GEN_GBTRF(z, c128)


#define GEN_GBTRS(PREFIX, TYPE) \
inline void \
call_gbtrs(char *trans, CBLAS_INT *n, CBLAS_INT *kl, CBLAS_INT *ku, CBLAS_INT *nrhs, TYPE *ab, CBLAS_INT *ldab, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gbtrs)(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info); \
};

GEN_GBTRS(s, f32)
GEN_GBTRS(d, f64)
GEN_GBTRS(c, c64)
GEN_GBTRS(z, c128)


// s- and d- versions of `gbcon` need integer iwork.
#define GEN_GBCON(PREFIX, TYPE) \
inline void \
call_gbcon(char *norm, CBLAS_INT *n, CBLAS_INT *kl, CBLAS_INT *ku, TYPE *ab, CBLAS_INT *ldab, CBLAS_INT *ipiv, TYPE *anorm, TYPE *rcond, TYPE *work, void *irwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gbcon)(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, (CBLAS_INT *)irwork, info); \
};

GEN_GBCON(s, f32)
GEN_GBCON(d, f64)


// c- and z- variants need floating type rwork instead of iwork.
#define GEN_GBCON_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gbcon(char *norm, CBLAS_INT *n, CBLAS_INT *kl, CBLAS_INT *ku, TYPE *ab, CBLAS_INT *ldab, CBLAS_INT *ipiv, RTYPE *anorm, RTYPE *rcond, TYPE *work, void *irwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gbcon)(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, (RTYPE *)irwork, info); \
};

GEN_GBCON_CZ(c, c64, f32)
GEN_GBCON_CZ(z, c128, f64)


#define GEN_GEQRF(PREFIX, TYPE) \
inline void \
call_geqrf(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *tau, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## geqrf)(m, n, a, lda, tau, work, lwork, info); \
};

GEN_GEQRF(s, f32);
GEN_GEQRF(d, f64);
GEN_GEQRF(c, c64);
GEN_GEQRF(z, c128);


// N.B. `rwork` is not used for `s` and `d` variants, so swallowed prior to calling LAPACK
#define GEN_GEQP3(PREFIX, TYPE) \
inline void \
call_geqp3(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *jpvt, TYPE *tau, TYPE *work, CBLAS_INT *lwork, void *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## geqp3)(m, n, a, lda, jpvt, tau, work, lwork, info); \
};

GEN_GEQP3(s, f32);
GEN_GEQP3(d, f64);


#define GEN_GEQP3_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_geqp3(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *jpvt, TYPE *tau, TYPE *work, CBLAS_INT *lwork, void *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## geqp3)(m, n, a, lda, jpvt, tau, work, lwork, (RTYPE *)rwork, info); \
};

GEN_GEQP3_CZ(c, c64, f32);
GEN_GEQP3_CZ(z, c128, f64);


// NB: wrap {s-,d-}orgqr for reals and {c-,z-}ungqr for complex
#define GEN_OR_UN_GQR(PREFIX, TYPE) \
inline void \
call_or_un_gqr(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *k, TYPE *a, CBLAS_INT *lda, TYPE *tau, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gqr)(m, n, k, a, lda, tau, work, lwork, info); \
};

GEN_OR_UN_GQR(sor, f32)
GEN_OR_UN_GQR(dor, f64)
GEN_OR_UN_GQR(cun, c64)
GEN_OR_UN_GQR(zun, c128)


/*
 * ?GESVD wrappers.
 *
 * We need to wrap over:
 *   - four type variants, s-, d-, c-, and zgesvd;
 *   - complex variants, c- and z-, receive the `rwork` argument, while s- and d- variants do not.
 * Thus,
 *   - `call_gesvd` has four overloads;
 *   - all variants receive the `rwork` argument; c- and z- variants forward it to LAPACK,
 *     and s- and d- variants swallow it.
 */
inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, f32 *a, CBLAS_INT *lda,
    f32 *s, f32 *u, CBLAS_INT *ldu, f32 *vt, CBLAS_INT *ldvt, f32 *work, CBLAS_INT *lwork,
    f32 *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(sgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
};

inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, f64 *a, CBLAS_INT *lda,
    f64 *s, f64 *u, CBLAS_INT *ldu, f64 *vt, CBLAS_INT *ldvt, f64 *work, CBLAS_INT *lwork,
    f64 *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(dgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
};

inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, c64 *a, CBLAS_INT *lda,
    f32 *s, c64 *u, CBLAS_INT *ldu, c64 *vt, CBLAS_INT *ldvt,
    c64 *work, CBLAS_INT *lwork, f32 *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(cgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
};

inline void call_gesvd(
    char *jobu, char *jobvt, CBLAS_INT *m, CBLAS_INT *n, c128 *a, CBLAS_INT *lda,
    f64 *s, c128 *u, CBLAS_INT *ldu, c128 *vt, CBLAS_INT *ldvt,
    c128 *work, CBLAS_INT *lwork, f64 *rwork, CBLAS_INT *info)
{
    BLAS_FUNC(zgesvd)(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
};



/*
 * ?GESDD wrappers.
 *
 * The logic is similar to ?gesdd:
 *   - we overload for four type variants, s-, d-, c-, and z-;
 *   - we forward `rwork` to c- and z- LAPACK functions, and swallow it for s- and d-;
 *
 */

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, f32 *a, CBLAS_INT *lda, f32 *s, f32 *u, CBLAS_INT *ldu,
    f32 *vt, CBLAS_INT *ldvt, f32 *work, CBLAS_INT *lwork, f32 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(sgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
};

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, f64 *a, CBLAS_INT *lda, f64 *s, f64 *u, CBLAS_INT *ldu,
    f64 *vt, CBLAS_INT *ldvt, f64 *work, CBLAS_INT *lwork, f64 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(dgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
};

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, c64 *a, CBLAS_INT *lda, f32 *s, c64 *u, CBLAS_INT *ldu,
    c64 *vt, CBLAS_INT *ldvt, c64 *work, CBLAS_INT *lwork, f32 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(cgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
};

inline void call_gesdd(
    char *jobz, CBLAS_INT *m, CBLAS_INT *n, c128 *a, CBLAS_INT *lda, f64 *s, c128 *u, CBLAS_INT *ldu,
    c128 *vt, CBLAS_INT *ldvt, c128 *work, CBLAS_INT *lwork, f64 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
{
    BLAS_FUNC(zgesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
};


// NB: s- and d- variants ignore the rwork argument (because LAPACK routines do not have it
#define GEN_GELSS_SD(PREFIX, TYPE, RTYPE) \
inline void \
call_gelss(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelss)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info); \
};

GEN_GELSS_SD(s, f32, f32)
GEN_GELSS_SD(d, f64, f64)


#define GEN_GELSS_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gelss(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelss)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info); \
};

GEN_GELSS_CZ(c, c64, f32)
GEN_GELSS_CZ(z, c128, f64)


// NB: s- and d- variants ignore the rwork argument (because LAPACK routines do not have it
#define GEN_GELSD_SD(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsd(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsd)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info); \
};

GEN_GELSD_SD(s, f32, f32)
GEN_GELSD_SD(d, f64, f64)


#define GEN_GELSD_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsd(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *s, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsd)(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info); \
};

GEN_GELSD_CZ(c, c64, f32)
GEN_GELSD_CZ(z, c128, f64)


// NB: s- and d- variants ignore the rwork argument (because LAPACK routines do not have it
#define GEN_GELSY_SD(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsy(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsy)(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info); \
};

GEN_GELSY_SD(s, f32, f32)
GEN_GELSY_SD(d, f64, f64)


#define GEN_GELSY_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_gelsy(CBLAS_INT *m, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *jpvt, RTYPE *rcond, CBLAS_INT *rank, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## gelsy)(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info); \
};

GEN_GELSY_CZ(c, c64, f32)
GEN_GELSY_CZ(z, c128, f64)


/*
 * ?GEEV wrappers.
 *
 * We need to wrap over:
 *   - four type variants, s-, d-, c-, and zgeev;
 *   - complex variants, c- and z-, receive the `rwork` argument, while s- and d- variants do not.
 *   - s- and d- variants return real and imaginary parts of eigenvalues separately, in *wr and *wi arrays
 *     c- and z- variants return a single complex array, *w, instead
 * Thus,
 *   - `call_geev` has four overloads;
 *   - all variants receive the `rwork` argument; c- and z- variants forward it to LAPACK,
 *     and s- and d- variants swallow it.
 *   - all variants have *wr and *wi arguments, both of the same type as *a
 *     (real for real *a, complex for complex *a);
 *     real-valued overloads, s- and d-, only fill *wr and ignore the *wi argument.
 */
#define GEN_GEEV_SD(PREFIX, TYPE) \
inline void \
call_geev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *wr, TYPE *wi, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork,  TYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## geev)(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info); \
};

GEN_GEEV_SD(s, f32)
GEN_GEEV_SD(d, f64)

#define GEN_GEEV_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_geev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *wr, TYPE *wi, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    /* ignore wi */ \
    BLAS_FUNC(PREFIX ## geev)(jobvl, jobvr, n, a, lda, wr, vl, ldvl, vr, ldvr, work, lwork, rwork, info); \
};

GEN_GEEV_CZ(c, c64, f32)
GEN_GEEV_CZ(z, c128, f64)


/*
 * Wrappers for ?GGEV
 *
 * The design is similar to that of ?GEEV wrappers: all overloads receive *rwork and *alphar, *alphai,
 */
#define GEN_GGEV_SD(PREFIX, TYPE) \
inline void \
call_ggev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, TYPE *alphar, TYPE *alphai, TYPE *beta, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork, TYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## ggev)(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info); \
};

GEN_GGEV_SD(s, f32)
GEN_GGEV_SD(d, f64)


#define GEN_GGEV_CZ(PREFIX, TYPE, RTYPE) \
inline void \
call_ggev(char *jobvl, char *jobvr, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, TYPE *alphar, TYPE *alphai, TYPE *beta, TYPE *vl, CBLAS_INT *ldvl, TYPE *vr, CBLAS_INT *ldvr, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## ggev)(jobvl, jobvr, n, a, lda, b, ldb, alphar, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info); \
};

GEN_GGEV_CZ(c, c64, f32)
GEN_GGEV_CZ(z, c128, f64)


/*
 * Wrappers for ?SY/HEEVR
 *
 * Discriminate between real and complex cases; the latter receive `rwork`, the former gobble that input.
 */
#define GEN_SYEVR(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_evr(char *jobz, char *range, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *vl, RTYPE *vu, CBLAS_INT *il, CBLAS_INT *iu, RTYPE *abstol, CBLAS_INT *m, RTYPE *w, TYPE *z, CBLAS_INT *ldz, CBLAS_INT *isuppz, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *lrwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## syevr)(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info); \
};

GEN_SYEVR(s, f32, f32);
GEN_SYEVR(d, f64, f64);


#define GEN_HEEVR(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_evr(char *jobz, char *range, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *vl, RTYPE *vu, CBLAS_INT *il, CBLAS_INT *iu, RTYPE *abstol, CBLAS_INT *m, RTYPE *w, TYPE *z, CBLAS_INT *ldz, CBLAS_INT *isuppz, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *lrwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## heevr)(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork, info); \
};

GEN_HEEVR(c, c64, f32);
GEN_HEEVR(z, c128, f64);


/*
 * Wrappers for ?SY/HEEV
 *
 * Discrimate between real and complex cases; the latter receive `rwork`, the former gobble that input.
 */
#define GEN_SYEV(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_ev(char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## syev)(jobz, uplo, n, a, lda, w, work, lwork, info); \
};

GEN_SYEV(s, f32, f32);
GEN_SYEV(d, f64, f64);


#define GEN_HEEV(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_ev(char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## heev)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info); \
};

GEN_HEEV(c, c64, f32);
GEN_HEEV(z, c128, f64);


/*
 * Wrappers for ?SY/HEEVD
 *
 * Discriminate between real and complex cases; the latter receive `rwork`, the former gobble that input
 */
#define GEN_SYEVD(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_evd(char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *lrwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## syevd)(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info); \
};

GEN_SYEVD(s, f32, f32);
GEN_SYEVD(d, f64, f64);


#define GEN_HEEVD(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_evd(char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *lrwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## heevd)(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info); \
};

GEN_HEEVD(c, c64, f32);
GEN_HEEVD(z, c128, f64);


/*
 * Wrappers for ?SY/HEEVX
 *
 * Discriminate between real and complex cases; the latter receive `rwork`, the former gobble that input
 */
#define GEN_SYEVX(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_evx(char *jobz, char *range, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *vl, RTYPE *vu, CBLAS_INT *il, CBLAS_INT *iu, RTYPE *abstol, CBLAS_INT *m, RTYPE *w, TYPE *z, CBLAS_INT *ldz, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## syevx)(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info); \
};

GEN_SYEVX(s, f32, f32);
GEN_SYEVX(d, f64, f64);


#define GEN_HEEVX(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_evx(char *jobz, char *range, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, RTYPE *vl, RTYPE *vu, CBLAS_INT *il, CBLAS_INT *iu, RTYPE *abstol, CBLAS_INT *m, RTYPE *w, TYPE *z, CBLAS_INT *ldz, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## heevx)(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info); \
};

GEN_HEEVX(c, c64, f32);
GEN_HEEVX(z, c128, f64);


/*
 * Wrappers for ?SY/HEGV
 *
 * Discriminate between real and complex cases due to `rwork`
 */
#define GEN_SYGV(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_gv(CBLAS_INT *itype, char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## sygv)(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info); \
};

GEN_SYGV(s, f32, f32);
GEN_SYGV(d, f64, f64);


#define GEN_HEGV(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_gv(CBLAS_INT *itype, char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## hegv)(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info); \
};

GEN_HEGV(c, c64, f32);
GEN_HEGV(z, c128, f64);


/*
 * Wrappers for ?SY/HEGVD
 *
 * Discriminate between real and complex cases due to `rwork`
 */
#define GEN_SYGVD(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_gvd(CBLAS_INT *itype, char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *lrwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## sygvd)(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info); \
};

GEN_SYGVD(s, f32, f32);
GEN_SYGVD(d, f64, f64);


#define GEN_HEGVD(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_gvd(CBLAS_INT *itype, char *jobz, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *w, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *lrwork, CBLAS_INT *iwork, CBLAS_INT *liwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## hegvd)(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info); \
};

GEN_HEGVD(c, c64, f32);
GEN_HEGVD(z, c128, f64);


/*
 * Wrappers for ?SY/HEGVX
 *
 * Disriminate between real and complex cases due to `rwork`
 */
#define GEN_SYGVX(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_gvx(CBLAS_INT *itype, char *jobz, char *range, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *vl, RTYPE *vu, CBLAS_INT *il, CBLAS_INT *iu, RTYPE *abstol, CBLAS_INT *m, RTYPE *w, TYPE *z, CBLAS_INT *ldz, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## sygvx)(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info); \
};

GEN_SYGVX(s, f32, f32);
GEN_SYGVX(d, f64, f64)


#define GEN_HEGVX(PREFIX, TYPE, RTYPE) \
inline void \
call_sy_he_gvx(CBLAS_INT *itype, char *jobz, char *range, char *uplo, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, RTYPE *vl, RTYPE *vu, CBLAS_INT *il, CBLAS_INT *iu, RTYPE *abstol, CBLAS_INT *m, RTYPE *w, TYPE *z, CBLAS_INT *ldz, TYPE *work, CBLAS_INT *lwork, RTYPE *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## hegvx)(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info); \
};

GEN_HEGVX(c, c64, f32);
GEN_HEGVX(z, c128, f64);


// Structure tags; python side maps assume_a strings to these values
enum St : Py_ssize_t
{
    NONE = -1,
    GENERAL = 0,
    DIAGONAL = 11,
    TRIDIAGONAL = 31,
    BANDED = 41,
    UPPER_TRIANGULAR = 21,
    LOWER_TRIANGULAR = 22,
    POS_DEF = 101,
    SYM = 201,
    HER = 211
};

// QR mode tags; python side maps mode strings to these values
enum QR_mode : Py_ssize_t
{
    FULL = 1,
    R = 11,
    RAW_MODE = 21,
    ECONOMIC = 31
};

// Eigh driver tags; pythn side maps driver strings to these values
enum Eigh_driver : Py_ssize_t
{
    EV = 1,
    EVD = 2,
    EVR = 3,
    EVX = 4,
    GV = 10,
    GVD = 11,
    GVX = 12
};


/*
 * Rich return object
 */
struct SliceStatus {
    Py_ssize_t slice_num;
    Py_ssize_t structure;
    int is_singular;
    int is_ill_conditioned;
    double rcond;
    Py_ssize_t lapack_info;
};


void init_status(SliceStatus& slice_status, npy_intp idx, St slice_structure) {
    slice_status.slice_num = idx;
    slice_status.structure = (Py_ssize_t)slice_structure;
    slice_status.is_singular = 0;
    slice_status.is_ill_conditioned = 0;
    slice_status.rcond = 0;
    slice_status.lapack_info = 0;
}


typedef std::vector<SliceStatus> SliceStatusVec;


/*
 * When looping over slices, each `solve_slice_XYZ` return the `status` struct,
 * which records one of three possible outcomes:
 *
 *  - the slice is singular or LAPACK returned non-zero `info`
 *  - the slice is detected to be ill-conditioned
 *  - all is well
 *
 * In the first two cases, we record the `status` in `vec_status`.
 * For non-recoverable errors (singularity of non-zero info), we terminate the loop over the slices.
 */
int
_detect_problems(const SliceStatus& slice_status, SliceStatusVec& vec_status) {
    if ((slice_status.lapack_info < 0) || (slice_status.is_singular)) {
        vec_status.push_back(slice_status);
        return 1;
    }
    else if (slice_status.is_ill_conditioned) {
        vec_status.push_back(slice_status);
    }
    return 0;
}


/*
 * lwork defensive handler :
 *  cf https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/lapack.py#L1004
 *
 *  Round floating-point lwork returned by lapack to integer.
 *
 *  Several LAPACK routines compute optimal values for LWORK, which
 *  they return in a floating-point variable. However, for large
 *  values of LWORK, single-precision floating point is not sufficient
 *  to hold the exact value --- some LAPACK versions (<= 3.5.0 at
 *  least) truncate the returned integer to single precision and in
 *  some cases this can be smaller than the required value.
 *
 *  The fudge_factor comes from
 *  https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/_basic.py#L1154
 *
 *  A fudge factor of 1% was added in commit https://github.com/scipy/scipy/commit/dfb543c147c
 *  to avoid a "curious segfault with 500x500 matrices and OpenBLAS".
 */
template<typename T>
CBLAS_INT _calc_lwork(T _lwrk, f64 fudge_factor=1.0) {
    using real_type = real_of_t<T>;

    real_type value = std::real(_lwrk) * fudge_factor;
    if((std::is_same<real_type, f32>::value) ||
       (std::is_same<real_type, c64>::value)
    ) {
        // Single-precision routine -- take next fp value to work
        // around possible truncation in LAPACK code
        value = std::nextafter(value, std::numeric_limits<real_type>::infinity());
    }

    CBLAS_INT lwork;
    if ((value < 0) || !(value <= (real_type)std::numeric_limits<CBLAS_INT>::max())) {
        // Too large lwork required - Computation cannot be performed with standard LAPACK
        lwork = -1;
    }
    else {
        lwork = value > 0 ? (CBLAS_INT)value : 1;
    }
    return lwork;
}


/*
 * Given an array of ndim >= 2, compute a pointer to the start of the 2D "slice" idx
 */
template<typename T>
T* compute_slice_ptr(npy_intp idx, T *Am_data, npy_intp ndim, npy_intp *shape, npy_intp *strides) {
    npy_intp offset = 0;
    npy_intp temp_idx = idx;
    for (int i = ndim - 3; i >= 0; i--) {
        offset += (temp_idx % shape[i]) * strides[i];
        temp_idx /= shape[i];
    }
    T* slice_ptr = (T *)(Am_data + (offset/sizeof(T)));
    return slice_ptr;
}


/*
 * Copy n-by-m slice from slice_ptr to dst.
 */
template<typename T>
void copy_slice(T* dst, const T* slice_ptr, const npy_intp n, const npy_intp m, const npy_intp s2, const npy_intp s1) {

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i * m + j] = *(slice_ptr + (i*s2/sizeof(T)) + (j*s1/sizeof(T)));
        }
    }
}


/*
 * Copy n-by-m C-order slice from slice_ptr to dst in F-order.
 *
 * `src` is n-by-m, strided
 * `dst` is ldb-by-m, F-ordered.
 *
 * The default is to have src and dst of the same size (ldb=-1 means ldb=n).
 */
template<typename T>
void copy_slice_F(T* dst, const T* slice_ptr, const npy_intp n, const npy_intp m, const npy_intp s2, const npy_intp s1, npy_intp ldb=-1) {

    if (ldb == -1) {ldb = n;}

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i + j*ldb] = *(slice_ptr + (i*s2/sizeof(T)) + (j*s1/sizeof(T)));  // == src[i*m + j]
        }
    }
}


/*
 * Copy n-by-m F-ordered `src` to C-ordered `dst`.
 *
 * `src` is ldb-by-m, F-ordered
 * `dst` is n-by-m, C-ordered
 *
 * The default is to have src and dst of the same size (ldb=-1 means ldb=n)
 */
template<typename T>
void copy_slice_F_to_C(T* dst, const T* src, const npy_intp n, const npy_intp m, npy_intp ldb=-1) {

    if (ldb == -1) {ldb = n;}

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i*m + j] = src[i + j*ldb];
        }
    }
}


/*
 * Copy one triangle of a strided m-by-n `src` to C-ordered `dst`.
 * Note that by specifying m and/or n smaller than the input array `src`
 * it is possible to only copy a subset of the triangle.
 *
 * Only elements in the `uplo` triangle are copied;
 * `dst` is assumed zero-initialized for the other triangle.
 *
 * The function does not use the symmetry of `src` — it reads exactly the
 * triangle specified by `uplo`. The caller is responsible for ensuring
 * that the correct triangle is populated in `src`.
 *
 * Examples (matrix dimension m x n):
 *   C-contiguous src (s0=n, s1=1): src[i*n + j] -> dst[i*n + j]
 *      both sequential in the inner loop - effectively a partial memcpy.
 *   F-contiguous src (s0=1, s2=m): src[i + j*m] -> dst[i*n + j],
 *      column-sequential reads, row-sequential writes.
 */
template<typename T>
void copy_triangle_to_C(T *dst, const T *src, const npy_intp m, const npy_intp n, const npy_intp s0, const npy_intp s1, const char uplo) {
    if (uplo == 'L') {
        for (npy_intp i = 0; i < m; i++) {
            npy_intp stop = std::min(i + 1, n);
            for (npy_intp j = 0; j < stop; j++) {
                dst[i * n + j] = *(src + i * s0 + j * s1);
            }
        }
    } else {
        for (npy_intp i = 0; i < m; i++) {
            for (npy_intp j = i; j < n; j++) {
                dst[i * n + j] = *(src + i * s0 + j * s1);
            }
        }
    }
}


/*
 * 1-norm of a matrix
 */

template<typename T>
real_of_t<T>
norm1_(T* A, const npy_intp n)
{
    using real_type = real_of_t<T>;

    real_type norm = 0.0;
    for (CBLAS_INT i = 0; i < n; i++) {
        real_type tmp = 0.0;
        for (CBLAS_INT j = 0; j < n; j++) {
            tmp += std::abs(A[i * n + j]);
        }

        if (tmp > norm) { norm = tmp; }
    }

    return norm;
}


template<typename T>
real_of_t<T>
norm1_sym_herm_upper(T* A, T* work, const npy_intp n)
{
    using real_type = real_of_t<T>;

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { rwork[i] = std::abs(A[i]);
     }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) {
        // only loop over the upper triangle
        rwork[i] += std::abs(A[i*n + i]);
        for (j = i+1; j < n; j++) {
            temp = std::abs(A[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
real_of_t<T>
norm1_sym_herm_lower(T* A, T* work, const npy_intp n)
{
    using real_type = real_of_t<T>;

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    for (i = 0; i < n; i++) { rwork[i] = 0.0; }

    for (i=0; i < n; i++) {
        rwork[i] += std::abs(A[i*n + i]);
        for (j=0; j < i; j++) {
            temp = std::abs(A[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }

    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
real_of_t<T>
norm1_sym_herm(char uplo, T *A, T *work, const npy_intp n) {
    // NB: transpose for the F order
    if (uplo == 'U') {return norm1_sym_herm_lower(A, work, n);}
    else if (uplo == 'L') {return norm1_sym_herm_upper(A, work, n);}
    else {throw std::runtime_error("uplo at norms");}
}


template<typename T>
real_of_t<T>
norm1_tridiag(T* dl, T *d, T *du, T *work, const npy_intp n) {
    using real_type = real_of_t<T>;

    real_type *rwork = (real_type *)work;

    npy_intp i;
    for (i=0; i<n; i++) {
        rwork[i] = std::abs(d[i]);
    }
    for (i=0; i<n-1; i++) {
        rwork[i] += std::abs(dl[i]);
    }
    for (i=1; i<n-1; i++) {
        rwork[i] += std::abs(du[i-1]);
    }

    real_type temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}

/*
 * Compute the 1 norm of a matrix `A`, but assume it is already in its banded
 * form `ab` as constructed by `to_banded`. It is assumed that the size of `ab`
 * is always such that its number of rows is `2 * kl + ku + 1`.
 */
template <typename T>
real_of_t<T>
norm1_banded(T* ab, const npy_intp kl, const npy_intp ku, T* work, const npy_intp n) {
    using real_type = real_of_t<T>;

    real_type *rwork = (real_type *)work;

    npy_intp i, j;
    npy_intp ldab = 2 * kl + ku + 1;

    for (i = 0; i < n; i++) {
        rwork[i] = std::abs(ab[i * ldab + kl + ku]);
    }

    for (i = 0; i < kl; i++) { // run over lower bands
        for (j = 0; j < n - i - 1; j++) {
            rwork[j] += std::abs(ab[j * ldab + kl + ku + i + 1]);
        }
    }


    for (i = 0; i < ku; i++) { // run over upper bands
        for (j = i + 1; j < n; j++) {
            rwork[j] += std::abs(ab[j * ldab + kl + ku - i - 1]);
        }
    }

    real_type temp = 0.0;
    for (i = 0; i < n; i++) {if (rwork[i] > temp) {temp = rwork[i];} }
    return temp;
}


/***************************
 ***  Structure detection
 ***************************/

template<typename T>
void
bandwidth(T* data, npy_intp n, npy_intp m, npy_intp* lower_band, npy_intp* upper_band)
{
    T zero = T(0.);

    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t c = 0; c < m-1; c++)
    {
        for (Py_ssize_t r = n-1; r > c + lb; r--)
        {
            if (data[c*n + r] != zero) { lb = r - c; break; }
        }
        if (c + lb + 1 > m) { break; }
    }
    for (Py_ssize_t c = m-1; c > 0; c--)
    {
        for (Py_ssize_t r = 0; r < c - ub; r++)
        {
            if (data[c*n + r] != zero) { ub = c - r; break; }

        }
        if (c <= ub) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


/*
 * Overload of the original `bandwidth` function that allows to take into
 * account the strides of the matrix to avoid having to explicitly set a
 * flag regarding the ordering of the matrix.
 *
 * The addressing is done using `npy_intp` instead of `Py_ssize_t` for
 * consistency.
 */
template<typename T>
void
bandwidth_strided(T* data, npy_intp n, npy_intp m, npy_intp s1, npy_intp s2, npy_intp *lower_band, npy_intp *upper_band)
{
    T zero = T(0.);

    s1 = s1 / sizeof(T);
    s2 = s2 / sizeof(T);
    npy_intp lb = 0, ub = 0;
    for (npy_intp c = 0; c < m-1; c++) {
        for (npy_intp r = n-1; r > c + lb; r--) {
            if (data[c * s2 + r * s1] != zero) { lb = r - c; break; }
        }
        if (c + lb + 1 > m) { break; }
    }
    for (npy_intp c = m-1; c > 0; c--) {
        for (npy_intp r = 0; r < c - ub; r++) {
            if (data[c * s2 + r * s1] != zero) { ub = c - r; break; }
        }
        if (c <= ub) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


template<typename T>
void
detect_bandwidths(T* data, npy_intp ndim, npy_intp outer_size, npy_intp *shape, npy_intp *strides, npy_intp *kl, npy_intp *ku, npy_intp *ldab_max) {
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        T* slice_ptr = compute_slice_ptr(idx, data, ndim, shape, strides);

        bandwidth_strided(slice_ptr, shape[ndim-2], shape[ndim-1], strides[ndim-2], strides[ndim-1], &kl[idx], &ku[idx]);
        if (2 * kl[idx] + ku[idx] + 1 > *ldab_max) {
            *ldab_max = 2 * kl[idx] + ku[idx] + 1;
        }
    }
}


template<typename T>
std::tuple<bool, bool>
is_sym_or_herm(const T *data, npy_intp n) {
    // Return a pair of (is_symmetric, is_hermitian)
    bool all_sym = true, all_herm = true;

    for (npy_intp i=0; i < n; i++) {
        for (npy_intp j=0; j < n; j++) {
            T elem1 = data[i*n + j];
            T elem2 = data[i + j*n];
            all_sym = all_sym && (elem1 == elem2);
            all_herm = all_herm && (elem1 == std::conj(elem2));
            if(!(all_sym || all_herm)) {
                // short-circuit : it's neither symmetric not hermitian
                return std::make_tuple(false, false);
            }
        }
    }
    return std::make_tuple(all_sym, all_herm);
}


template<typename T>
inline void
swap_cf(T* src, T* dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    T *bb = dst;
    T *aa = src;
    if ((r < 16) && (c < 16)) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf(src, dst, r2, c, n);
            swap_cf(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf(src, dst, r, c2, n);
            swap_cf(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
}


/*
 * Common matrices
 */

// fill np.triu(a) from np.tril(a) or np.tril(a) from np.triu(a)
template<typename T>
inline void
fill_other_triangle(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                if constexpr(is_complex_v<T>) {
                    data[j + i*n] = std::conj(data[i + j*n]);
                } else {
                    data[j + i*n] = data[i + j*n];
                }
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i+1; j++){
                if constexpr(is_complex_v<T>) {
                    data[j + i*n] = std::conj(data[i + j*n]);
                } else {
                    data[j + i*n] = data[i + j*n];
                }
            }
        }
    }
}

// XXX deduplicate (conj or noconj)
template<typename T>
inline void
fill_other_triangle_noconj(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                data[j + i*n] = data[i + j*n];
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i+1; j++){
                data[j + i*n] = data[i + j*n];
            }
        }
    }
}


/*
 * Helper for converting an NxN matrix to tridiagonal form:
 * extract the diagonals of `data` (a full NxN matrix) into du, d, dl
 *
 * the upper and lower subdiagonals, `dl` and `du`, are have length N-1
 * the main diagonal, `d`, has length N
 *
 */
template<typename T>
inline void
to_tridiag(const T *data, npy_intp N, T *du, T *d, T *dl) {
    for (npy_intp i=0; i<N; i++) {
        d[i] = data[i + i*N];
    }
    for (npy_intp i=0; i<N-1; i++) {
        dl[i] = data[i*N + i + 1];
    }
    for (npy_intp i=0; i<N-1; i++) {
        du[i] = data[(i+1)*N + i];
    }
}


/*
 * Helper function for reshuffling a banded slice into the appropriate
 * structure for ?gbcon and ?gbtrf. `s1` and `s2` contain the strides in
 * the column and row direction (`ndim` - 2 and `ndim` - 1, respectively).
 * The result is stored in `ab` in Fortran order.
 *
 * It is assumed that `ab` provides at least `ldab` x `n` memory elements,
 * where ldab >= 2 * `kl` + `ku` + 1
 *
 * Reference: https://www.netlib.org/lapack/explore-html/df/dd6/group__gbtrf_ga682f53142f0398f83f5461c277d23ba2.html#ga682f53142f0398f83f5461c277d23ba2
 */
template<typename T>
inline void
to_banded(const T *data, npy_intp n, npy_intp kl, npy_intp ku, npy_intp ldab, T *ab, npy_intp s1, npy_intp s2) {
    s1 = s1 / sizeof(T);
    s2 = s2 / sizeof(T);
    npy_intp i, j;

    // main diagonal
    for (i = 0; i < n; i++) {
        ab[(i + 1) * ldab - kl - 1] = data[i * (s1 + s2)];
    }

    // lower bands
    for (i = 0; i < kl; i++) {
        for (j = 0; j < n - i - 1; j++) {
            ab[(j + 1) * ldab - kl + i] = data[(i + 1) * s1 + j * (s1 + s2)];
        }
    }

    // upper bands
    for (i = 0; i < ku; i++) {
        for (j = i + 1; j < n; j++) {
            ab[(j + 1) * ldab - kl - i - 2] = data[(i + 1) * s2 + (j - i - 1) * (s1 + s2)];
        }
    }
}


template<typename T>
inline void
zero_other_triangle(char uplo, T *data, const npy_intp m, npy_intp n = -1, npy_intp lda = -1) {
    if (n == -1) { n = m; }
    if (lda == -1) { lda = m; }

    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<m; j++){
                data[j + i*lda] = 0.0;
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            npy_intp stop = std::min(i, m);
            for (npy_intp j=0; j < stop; j++){
                data[j + i*lda] = 0.0;
            }
        }
    }
}

} // namespace sp_linalg
