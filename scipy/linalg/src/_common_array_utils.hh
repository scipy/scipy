/*
 * LAPACK declarations.
 */
#ifndef _SCIPY_COMMON_ARRAY_UTILS_H
#define _SCIPY_COMMON_ARRAY_UTILS_H
#include "Python.h"
#include "numpy/npy_math.h"
#include "npy_cblas.h"

/*
 * declare LAPACK prototypes
 */

extern "C" {

/* ?GETRF */
CBLAS_INT
BLAS_FUNC(sgetrf)(CBLAS_INT *m, CBLAS_INT *n, float a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(dgetrf)(CBLAS_INT *m, CBLAS_INT *n, double a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(cgetrf)(CBLAS_INT *m, CBLAS_INT *n, npy_complex64 a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(zgetrf)(CBLAS_INT *m, CBLAS_INT *n, npy_complex128 a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);


/* ?GETRI */
CBLAS_INT
BLAS_FUNC(sgetri)(CBLAS_INT *n, float a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  float work[], CBLAS_INT *lwork, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(dgetri)(CBLAS_INT *n, double a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  double work[], CBLAS_INT *lwork, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(cgetri)(CBLAS_INT *n, npy_complex64 a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  npy_complex64 work[], CBLAS_INT *lwork, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(zgetri)(CBLAS_INT *n, npy_complex128 a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  npy_complex128 work[], CBLAS_INT *lwork, CBLAS_INT *info
);


/* ?GETRS */
void BLAS_FUNC(sgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, CBLAS_INT *ipiv, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, CBLAS_INT *ipiv, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(cgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(zgetrs)(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *ipiv, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);


/* ?GECON */
void BLAS_FUNC(sgecon)(char* norm, CBLAS_INT* n, float* a,       CBLAS_INT* lda, float* anorm,  float* rcond,  float* work,       CBLAS_INT* iwork, CBLAS_INT* info);
void BLAS_FUNC(dgecon)(char* norm, CBLAS_INT* n, double* a,      CBLAS_INT* lda, double* anorm, double* rcond, double* work,      CBLAS_INT* iwork, CBLAS_INT* info);
void BLAS_FUNC(cgecon)(char* norm, CBLAS_INT* n, npy_complex64* a,  CBLAS_INT* lda, float* anorm,  float* rcond,  npy_complex64* work,  float* rwork,     CBLAS_INT* info);
void BLAS_FUNC(zgecon)(char* norm, CBLAS_INT* n, npy_complex128* a, CBLAS_INT* lda, double* anorm, double* rcond, npy_complex128* work, double* rwork,    CBLAS_INT* info);


/* ?TRTRI */
void BLAS_FUNC(strtri)(char *uplo, char *diag, CBLAS_INT *n, float* a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(dtrtri)(char *uplo, char *diag, CBLAS_INT *n, double* a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(ctrtri)(char *uplo, char *diag, CBLAS_INT *n, npy_complex64* a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(ztrtri)(char *uplo, char *diag, CBLAS_INT *n, npy_complex128* a, CBLAS_INT *lda, CBLAS_INT *info);

/* ?TRCON */
void BLAS_FUNC(strcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, float *a, CBLAS_INT *lda, float *rcond, float *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(dtrcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, double *a, CBLAS_INT *lda, double *rcond, double *work, CBLAS_INT *iwork, CBLAS_INT *info);
void BLAS_FUNC(ctrcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, float *rcond, npy_complex64 *work, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(ztrcon)(char *norm, char *uplo, char *diag, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, double *rcond, npy_complex128 *work, double *rwork, CBLAS_INT *info);

/* ?TRTRS */
void BLAS_FUNC(strtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dtrtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(ctrtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(ztrtrs)(char *uplo, char *trans, char *diag, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);

/* ?POTRF */
void BLAS_FUNC(spotrf)(char *uplo, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(dpotrf)(char *uplo, CBLAS_INT *n, double *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(cpotrf)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(zpotrf)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *info);

/* ?POTRI */
void BLAS_FUNC(spotri)(char *uplo, CBLAS_INT *n, float *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(dpotri)(char *uplo, CBLAS_INT *n, double *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(cpotri)(char *uplo, CBLAS_INT *n, npy_complex64 *a, CBLAS_INT *lda, CBLAS_INT *info);
void BLAS_FUNC(zpotri)(char *uplo, CBLAS_INT *n, npy_complex128 *a, CBLAS_INT *lda, CBLAS_INT *info);

/* ?POCON */
void BLAS_FUNC(spocon)(char *uplo, CBLAS_INT *n, float* a, CBLAS_INT *lda, float *anorm, float *rcond, float* work, CBLAS_INT* iwork, CBLAS_INT *info);
void BLAS_FUNC(dpocon)(char *uplo, CBLAS_INT *n, double* a, CBLAS_INT *lda, double *anorm, double *rcond, double* work, CBLAS_INT* iwork, CBLAS_INT *info);
void BLAS_FUNC(cpocon)(char *uplo, CBLAS_INT *n, npy_complex64* a, CBLAS_INT *lda, float *anorm, float *rcond, npy_complex64* work, float *rwork, CBLAS_INT *info);
void BLAS_FUNC(zpocon)(char *uplo, CBLAS_INT *n, npy_complex128* a, CBLAS_INT *lda, double *anorm, double *rcond, npy_complex128* work, double *rwork, CBLAS_INT *info);

/* ?POTRS*/
void BLAS_FUNC(spotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, float *a, CBLAS_INT *lda, float *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(dpotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, double *a, CBLAS_INT *lda, double *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(cpotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex64 *a, CBLAS_INT *lda, npy_complex64 *b, CBLAS_INT *ldb, CBLAS_INT *info);
void BLAS_FUNC(zpotrs)(char *uplo, CBLAS_INT *n, CBLAS_INT *nrhs, npy_complex128 *a, CBLAS_INT *lda, npy_complex128 *b, CBLAS_INT *ldb, CBLAS_INT *info);


} // extern "C"


/*
 * Generate type overloads, to map from C array types (float, double, npy_complex64, npy_complex128)
 * to LAPACK prefixes, "sdcz".
 */
#define GEN_GETRF(PREFIX, TYPE) \
inline void \
getrf(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getrf)(m, n, a, lda, ipiv, info); \
};

GEN_GETRF(s,float)
GEN_GETRF(d,double)
GEN_GETRF(c,npy_complex64)
GEN_GETRF(z,npy_complex128)


#define GEN_GETRS(PREFIX, TYPE) \
inline void \
getrs(char *trans, CBLAS_INT *n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getrs)(trans, n, nrhs, a, lda, ipiv, b, ldb, info); \
};

GEN_GETRS(s,float)
GEN_GETRS(d,double)
GEN_GETRS(c,npy_complex64)
GEN_GETRS(z,npy_complex128)


#define GEN_GETRI(PREFIX, TYPE) \
inline void \
getri(CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getri)(n, a, lda, ipiv, work, lwork, info); \
};

GEN_GETRI(s,float)
GEN_GETRI(d,double)
GEN_GETRI(c,npy_complex64)
GEN_GETRI(z,npy_complex128)


// NB: iwork for real arrays or rwork for complex arrays
#define GEN_GECON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
gecon(char* norm, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## gecon)(norm, n, a, lda, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_GECON(s, float, float, CBLAS_INT)
GEN_GECON(d, double, double, CBLAS_INT)
GEN_GECON(c, npy_complex64, float, float)
GEN_GECON(z, npy_complex128, double, double)


#define GEN_TRTRI(PREFIX, TYPE) \
inline void \
trtri(char* uplo, char *diag, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## trtri)(uplo, diag, n, a, lda, info); \
};

GEN_TRTRI(s, float)
GEN_TRTRI(d, double)
GEN_TRTRI(c, npy_complex64)
GEN_TRTRI(z, npy_complex128)


#define GEN_TRCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
trcon(char* norm, char *uplo, char *diag, CBLAS_INT *n, CTYPE *a, CBLAS_INT *lda, RTYPE *rcond, CTYPE *work, void *irwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## trcon)(norm, uplo, diag, n, a, lda, rcond, work, (WTYPE *)irwork, info); \
};

GEN_TRCON(s, float, float, CBLAS_INT)
GEN_TRCON(d, double, double, CBLAS_INT)
GEN_TRCON(c, npy_complex64, float, float)
GEN_TRCON(z, npy_complex128, double, double)


#define GEN_TRTRS(PREFIX, TYPE) \
inline void \
trtrs(char* uplo, char *trans, char *diag, CBLAS_INT* n, CBLAS_INT* nrhs, TYPE* a, CBLAS_INT* lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## trtrs)(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info); \
};

GEN_TRTRS(s, float)
GEN_TRTRS(d, double)
GEN_TRTRS(c, npy_complex64)
GEN_TRTRS(z, npy_complex128)


#define GEN_POTRF(PREFIX, TYPE) \
inline void \
potrf(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potrf)(uplo, n, a, lda, info); \
};

GEN_POTRF(s, float)
GEN_POTRF(d, double)
GEN_POTRF(c, npy_complex64)
GEN_POTRF(z, npy_complex128)


#define GEN_POTRI(PREFIX, TYPE) \
inline void \
potri(char* uplo, CBLAS_INT* n, TYPE* a, CBLAS_INT* lda, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potri)(uplo, n, a, lda, info); \
};

GEN_POTRI(s, float)
GEN_POTRI(d, double)
GEN_POTRI(c, npy_complex64)
GEN_POTRI(z, npy_complex128)


// NB: iwork for real arrays or rwork for complex arrays
#define GEN_POCON(PREFIX, CTYPE, RTYPE, WTYPE) \
inline void \
pocon(char* uplo, CBLAS_INT* n, CTYPE* a, CBLAS_INT* lda, RTYPE* anorm, RTYPE* rcond, CTYPE* work, void *irwork, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## pocon)(uplo, n, a, lda, anorm, rcond, work, (WTYPE *)irwork, info); \
};

GEN_POCON(s, float, float, CBLAS_INT)
GEN_POCON(d, double, double, CBLAS_INT)
GEN_POCON(c, npy_complex64, float, float)
GEN_POCON(z, npy_complex128, double, double)


#define GEN_POTRS(PREFIX, TYPE) \
inline void \
potrs(char* uplo, CBLAS_INT* n, CBLAS_INT *nrhs, TYPE *a, CBLAS_INT *lda, TYPE *b, CBLAS_INT *ldb, CBLAS_INT* info) \
{ \
    BLAS_FUNC(PREFIX ## potrs)(uplo, n, nrhs, a, lda, b, ldb, info); \
};

GEN_POTRS(s, float)
GEN_POTRS(d, double)
GEN_POTRS(c, npy_complex64)
GEN_POTRS(z, npy_complex128)



// Structure tags; python side maps assume_a strings to these values
enum St : Py_ssize_t
{
    NONE = -1,
    GENERAL = 0,
    UPPER_TRIANGULAR = 21,
    LOWER_TRIANGULAR = 22,
    POS_DEF = 101,
    POS_DEF_UPPER = 111,
    POS_DEF_LOWER = 112,
};


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
 */
template<typename T>
void copy_slice_F(T* dst, const T* slice_ptr, const npy_intp n, const npy_intp m, const npy_intp s2, const npy_intp s1) {

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i + j*n] = *(slice_ptr + (i*s2/sizeof(T)) + (j*s1/sizeof(T)));  // == src[i*m + j]
        }
    }
}

/*
 * Copy n-by-m F-ordered `src` to C-ordered `dst`.
 */
template<typename T>
void copy_slice_F_to_C(T* dst, const T* src, const npy_intp n, const npy_intp m) {

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i*m + j] = src[i + j*n];
        }
    }
}





/*
 * 1-norm of a matrix
 */

template<typename T>
typename type_traits<T>::real_type
norm1_(T* A, T* work, const npy_intp n)
{
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;
    value_type *pA = reinterpret_cast<value_type *>(A);

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { rwork[i] = std::abs(pA[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { rwork[j] += std::abs(pA[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
typename type_traits<T>::real_type
norm1_sym_herm_upper(T* A, T* work, const npy_intp n)
{
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;
    value_type *pA = reinterpret_cast<value_type *>(A);

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { rwork[i] = std::abs(pA[i]);
     }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) {
        // only loop over the upper triangle
        rwork[i] += std::abs(pA[i*n + i]);
        for (j = i+1; j < n; j++) {
            temp = std::abs(pA[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
typename type_traits<T>::real_type
norm1_sym_herm_lower(T* A, T* work, const npy_intp n)
{
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;
    value_type *pA = reinterpret_cast<value_type *>(A);

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    for (i = 0; i < n; i++) { rwork[i] = 0.0; }

    for (i=0; i < n; i++) {
        rwork[i] += std::abs(pA[i*n + i]);
        for (j=0; j < i; j++) {
            temp = std::abs(pA[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }

    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
typename type_traits<T>::real_type
norm1_sym_herm(char uplo, T *A, T *work, const npy_intp n) {
    // NB: transpose for the F order
    if (uplo == 'U') {return norm1_sym_herm_lower(A, work, n);}
    else if (uplo == 'L') {return norm1_sym_herm_upper(A, work, n);}
    else {throw std::runtime_error("uplo at norms");}
}


/***************************
 ***  Structure detection
 ***************************/

template<typename T>
void
bandwidth(T* data, npy_intp n, npy_intp m, npy_intp* lower_band, npy_intp* upper_band)
{
    using value_type = typename type_traits<T>::value_type;
    value_type *p_data = reinterpret_cast<value_type *>(data);
    value_type zero = value_type(0.);

    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t c = 0; c < m-1; c++)
    {
        for (Py_ssize_t r = n-1; r > c + lb; r--)
        {
            if (p_data[c*n + r] != zero) { lb = r - c; break; }
        }
        if (c + lb + 1 > m) { break; }
    }
    for (Py_ssize_t c = m-1; c > 0; c--)
    {
        for (Py_ssize_t r = 0; r < c - ub; r++)
        {
            if (p_data[c*n + r] != zero) { ub = c - r; break; }

        }
        if (c <= ub) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}




template<typename T>
bool
is_sym_herm(const T *data, npy_intp n) {
    using value_type = typename type_traits<T>::value_type;
    const value_type *p_data = reinterpret_cast<const value_type *>(data);

    for (npy_intp i=0; i < n; i++) {
        for (npy_intp j=0; j < n; j++) {
            if(p_data[i*n + j] != std::conj(p_data[i + j*n])) { return false; }
        }
    }
    return true;
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
                data[j + i*n] = conj(data[i + j*n]);
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i+1; j++){
                data[j + i*n] = conj(data[i + j*n]);
            }
        }
    }
}



template<typename T>
inline void
zero_other_triangle(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                data[j + i*n] = numeric_limits<T>::zero;
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i; j++){
                data[j + i*n] = numeric_limits<T>::zero;
            }
        }
    }
}


template<typename T>
inline void
nan_matrix(T * data, npy_intp n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            data[i * n + j] = numeric_limits<T>::nan;
        }
    }
}
#endif




