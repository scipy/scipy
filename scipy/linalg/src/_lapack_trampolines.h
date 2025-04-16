/*
 * LAPACK declarations and call trampolines.
 *
 * For each LAPACK function, `?FUNC`,
 *     1. declare the LAPACK prototypes
 *     2. declare the struct to hold the arguments. In a constructor, give arguments
 *        "default" values. E.g. set LDA to N; pointer to arrays, set to NULL.
 *        For a LAPACK function ?FUNC, the struct is `func_data_t`.
 *     3. declare/define the `call_func` overloads to map from C array types
 *        (float, double, npy_cfloat, npy_cdouble) to LAPACK prefixes, "sdcz".
 */
#pragma once
#include "Python.h"
#include "numpy/npy_math.h"
#include "npy_cblas.h"

typedef CBLAS_INT fortran_int;

/*
 * declare LAPACK prototypes
 */


/* ?GESV */
extern "C" {
fortran_int
BLAS_FUNC(sgesv)(fortran_int *n, fortran_int *nrhs, float a[], fortran_int *lda,
                 fortran_int ipiv[], float b[], fortran_int *ldb, fortran_int *info
);
fortran_int
BLAS_FUNC(dgesv)(fortran_int *n, fortran_int *nrhs, double a[], fortran_int *lda,
                 fortran_int ipiv[], double b[], fortran_int *ldb, fortran_int *info
);
fortran_int
BLAS_FUNC(cgesv)(fortran_int *n, fortran_int *nrhs, npy_cfloat a[], fortran_int *lda,
                 fortran_int ipiv[], npy_cfloat b[], fortran_int *ldb, fortran_int *info
);
fortran_int
BLAS_FUNC(zgesv)(fortran_int *n, fortran_int *nrhs, npy_cdouble a[], fortran_int *lda,
                 fortran_int ipiv[], npy_cdouble b[], fortran_int *ldb, fortran_int *info
);

/* ?GETRF */
fortran_int
BLAS_FUNC(sgetrf)(fortran_int *m, fortran_int *n, float a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);
fortran_int
BLAS_FUNC(dgetrf)(fortran_int *m, fortran_int *n, double a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);
fortran_int
BLAS_FUNC(cgetrf)(fortran_int *m, fortran_int *n, npy_cfloat a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);
fortran_int
BLAS_FUNC(zgetrf)(fortran_int *m, fortran_int *n, npy_cdouble a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);


/* ?GETRI */
fortran_int
BLAS_FUNC(sgetri)(fortran_int *n, float a[], fortran_int *lda, fortran_int ipiv[],
                  float work[], fortran_int *lwork, fortran_int *info
);
fortran_int
BLAS_FUNC(dgetri)(fortran_int *n, double a[], fortran_int *lda, fortran_int ipiv[],
                  double work[], fortran_int *lwork, fortran_int *info
);
fortran_int
BLAS_FUNC(cgetri)(fortran_int *n, npy_cfloat a[], fortran_int *lda, fortran_int ipiv[],
                  npy_cfloat work[], fortran_int *lwork, fortran_int *info
);
fortran_int
BLAS_FUNC(zgetri)(fortran_int *n, npy_cdouble a[], fortran_int *lda, fortran_int ipiv[],
                  npy_cdouble work[], fortran_int *lwork, fortran_int *info
);


} // extern "C"


/*
 * Hold the GESV related variables, handle allocation/deallocation.
 */
template<typename T>
struct gesv_data_t {
    fortran_int n;
    fortran_int nrhs;
    T *a;
    fortran_int lda;
    fortran_int *ipiv;
    fortran_int ldb;
    T *b;
    fortran_int info;

    gesv_data_t(fortran_int n_) :
        n(n_), nrhs(n_), lda(n_), ldb(n_), info(0)
    {
        a = (T *)malloc(n*n*sizeof(T));
        b = (T *)malloc(n*n*sizeof(T));
        ipiv = (fortran_int *)malloc(n*sizeof(fortran_int));
        if ((a==NULL) || (b == NULL) || (ipiv == NULL)) {
            PyErr_NoMemory();
            info = -1;
        }
    };

    ~gesv_data_t() {
        free(a);
        free(b);
        free(ipiv);
    };
};


/*
 * Trampoline from a C type to the BLAS prefix (sdcz)
 */
inline void call_gesv(gesv_data_t<float>& gesv_data) {
    BLAS_FUNC(sgesv)(&gesv_data.n, &gesv_data.nrhs, gesv_data.a, &gesv_data.lda,
                     gesv_data.ipiv, gesv_data.b, &gesv_data.ldb, &gesv_data.info
    );
}

inline void call_gesv(gesv_data_t<double>& gesv_data) {
    BLAS_FUNC(dgesv)(&gesv_data.n, &gesv_data.nrhs, gesv_data.a, &gesv_data.lda,
                     gesv_data.ipiv, gesv_data.b, &gesv_data.ldb, &gesv_data.info
    );
}

inline void call_gesv(gesv_data_t<npy_cfloat>& gesv_data) {
    BLAS_FUNC(cgesv)(&gesv_data.n, &gesv_data.nrhs, gesv_data.a, &gesv_data.lda,
                     gesv_data.ipiv, gesv_data.b, &gesv_data.ldb, &gesv_data.info
    );
}

inline void call_gesv(gesv_data_t<npy_cdouble>& gesv_data) {
    BLAS_FUNC(zgesv)(&gesv_data.n, &gesv_data.nrhs, gesv_data.a, &gesv_data.lda,
                     gesv_data.ipiv, gesv_data.b, &gesv_data.ldb, &gesv_data.info
    );
}


/*
 * Hold the GETRF related variables.
 *
 * No allocations/deallocations, initialize all array pointers to NULL and all other
 * out variables to sentinels (-101 etc).
 */
template<typename T>
struct getrf_data_t {
    fortran_int m;
    fortran_int n;
    T *a;
    fortran_int lda;
    fortran_int *ipiv;
    fortran_int info;

    getrf_data_t(fortran_int m_, fortran_int n_) :
        m(m_), n(n_), a(NULL), lda(n_), ipiv(NULL), info(-101)
    {};
};
inline void
call_getrf(getrf_data_t<float>& data) {
    BLAS_FUNC(sgetrf)(&data.m, &data.n, data.a, &data.lda, data.ipiv, &data.info);
}
inline void
call_getrf(getrf_data_t<double>& data) {
    BLAS_FUNC(dgetrf)(&data.m, &data.n, data.a, &data.lda, data.ipiv, &data.info);
}
inline void
call_getrf(getrf_data_t<npy_cfloat>& data) {
    BLAS_FUNC(cgetrf)(&data.m, &data.n, data.a, &data.lda, data.ipiv, &data.info);
}
inline void
call_getrf(getrf_data_t<npy_cdouble>& data) {
    BLAS_FUNC(zgetrf)(&data.m, &data.n, data.a, &data.lda, data.ipiv, &data.info);
}




/*
 * Hold the GETRI related variables.
 * Perform the workspace query, allocate the work array.
 * Is only usable after a call to ?getrf.
 */
template<typename T>
struct getri_data_t {
    fortran_int n;
    T *a;
    fortran_int lda;
    fortran_int *ipiv;
    T *work;
    fortran_int lwork;
    fortran_int info;

    getri_data_t(getrf_data_t<T>& data) :
        n(data.n), a(data.a), lda(data.lda), ipiv(data.ipiv), work(NULL), lwork(-1), info(-101)
    {};
};


inline void
call_getri(getri_data_t<float>& data) {
    BLAS_FUNC(sgetri)(&data.n, data.a, &data.lda, data.ipiv, data.work, &data.lwork, &data.info);
}
inline void
call_getri(getri_data_t<double>& data) {
    BLAS_FUNC(dgetri)(&data.n, data.a, &data.lda, data.ipiv, data.work, &data.lwork, &data.info);
}
inline void
call_getri(getri_data_t<npy_cfloat>& data) {
    BLAS_FUNC(cgetri)(&data.n, data.a, &data.lda, data.ipiv, data.work, &data.lwork, &data.info);
}
inline void
call_getri(getri_data_t<npy_cdouble>& data) {
    BLAS_FUNC(zgetri)(&data.n, data.a, &data.lda, data.ipiv, data.work, &data.lwork, &data.info);
}
