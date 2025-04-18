/*
 * LAPACK declarations.
 */
#pragma once
#include "Python.h"
#include "numpy/npy_math.h"
#include "npy_cblas.h"

/*
 * declare LAPACK prototypes
 */


/* ?GESV */
extern "C" {
CBLAS_INT
BLAS_FUNC(sgesv)(CBLAS_INT *n, CBLAS_INT *nrhs, float a[], CBLAS_INT *lda,
                 CBLAS_INT ipiv[], float b[], CBLAS_INT *ldb, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(dgesv)(CBLAS_INT *n, CBLAS_INT *nrhs, double a[], CBLAS_INT *lda,
                 CBLAS_INT ipiv[], double b[], CBLAS_INT *ldb, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(cgesv)(CBLAS_INT *n, CBLAS_INT *nrhs, npy_cfloat a[], CBLAS_INT *lda,
                 CBLAS_INT ipiv[], npy_cfloat b[], CBLAS_INT *ldb, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(zgesv)(CBLAS_INT *n, CBLAS_INT *nrhs, npy_cdouble a[], CBLAS_INT *lda,
                 CBLAS_INT ipiv[], npy_cdouble b[], CBLAS_INT *ldb, CBLAS_INT *info
);

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
BLAS_FUNC(cgetrf)(CBLAS_INT *m, CBLAS_INT *n, npy_cfloat a[], CBLAS_INT *lda,
                  CBLAS_INT ipiv[], CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(zgetrf)(CBLAS_INT *m, CBLAS_INT *n, npy_cdouble a[], CBLAS_INT *lda,
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
BLAS_FUNC(cgetri)(CBLAS_INT *n, npy_cfloat a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  npy_cfloat work[], CBLAS_INT *lwork, CBLAS_INT *info
);
CBLAS_INT
BLAS_FUNC(zgetri)(CBLAS_INT *n, npy_cdouble a[], CBLAS_INT *lda, CBLAS_INT ipiv[],
                  npy_cdouble work[], CBLAS_INT *lwork, CBLAS_INT *info
);


} // extern "C"


/*
 * Generate type overloads, to map from C array types (float, double, npy_cfloat, npy_cdouble)
 * to LAPACK prefixes, "sdcz".
 */
#define GEN_GETRF(PREFIX, TYPE) \
inline void \
GETRF(CBLAS_INT *m, CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getrf)(m, n, a, lda, ipiv, info); \
};


#define GEN_GETRI(PREFIX, TYPE) \
inline void \
GETRI(CBLAS_INT *n, TYPE *a, CBLAS_INT *lda, CBLAS_INT *ipiv, TYPE *work, CBLAS_INT *lwork, CBLAS_INT *info) \
{ \
    BLAS_FUNC(PREFIX ## getri)(n, a, lda, ipiv, work, lwork, info); \
};


GEN_GETRF(s,float)
GEN_GETRF(d,double)
GEN_GETRF(c,npy_cfloat)
GEN_GETRF(z,npy_cdouble)

GEN_GETRI(s,float)
GEN_GETRI(d,double)
GEN_GETRI(c,npy_cfloat)
GEN_GETRI(z,npy_cdouble)

