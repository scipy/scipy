/*
 * LAPACK declarations and call trampolines.
 */
#pragma once
#include "Python.h"
#include <complex>
#include "numpy/npy_math.h"
#include "npy_cblas.h"

typedef CBLAS_INT fortran_int;


/* ?GESV */
extern "C" fortran_int
BLAS_FUNC(sgesv)(fortran_int *n, fortran_int *nrhs,
                 float a[], fortran_int *lda,
                 fortran_int ipiv[],
                 float b[], fortran_int *ldb,
                 fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(dgesv)(fortran_int *n, fortran_int *nrhs,
                 double a[], fortran_int *lda,
                 fortran_int ipiv[],
                 double b[], fortran_int *ldb,
                 fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(cgesv)(fortran_int *n, fortran_int *nrhs,
                 npy_cfloat a[], fortran_int *lda,
                 fortran_int ipiv[],
                 npy_cfloat b[], fortran_int *ldb,
                 fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(zgesv)(fortran_int *n, fortran_int *nrhs,
                 npy_cdouble a[], fortran_int *lda,
                 fortran_int ipiv[],
                 npy_cdouble b[], fortran_int *ldb,
                 fortran_int *info
);

/* ?GETRF */
extern "C" fortran_int
BLAS_FUNC(sgetrf)(fortran_int *m, fortran_int *n, float a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(dgetrf)(fortran_int *m, fortran_int *n, double a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(cgetrf)(fortran_int *m, fortran_int *n, npy_cfloat a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(zgetrf)(fortran_int *m, fortran_int *n, npy_cdouble a[], fortran_int *lda,
                  fortran_int ipiv[], fortran_int *info
);


/* ?GETRI */
extern "C" fortran_int
BLAS_FUNC(sgetri)(fortran_int *n, float a[], fortran_int *lda,
                  fortran_int ipiv[], float work[], fortran_int *lwork, fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(dgetri)(fortran_int *n, double a[], fortran_int *lda,
                  fortran_int ipiv[], double work[], fortran_int *lwork, fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(cgetri)(fortran_int *n, npy_cfloat a[], fortran_int *lda,
                  fortran_int ipiv[], npy_cfloat work[], fortran_int *lwork, fortran_int *info
);
extern "C" fortran_int
BLAS_FUNC(zgetri)(fortran_int *n, npy_cdouble a[], fortran_int *lda,
                  fortran_int ipiv[], npy_cdouble work[], fortran_int *lwork, fortran_int *info
);



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
    BLAS_FUNC(sgesv)(
        &gesv_data.n,
        &gesv_data.nrhs,
        gesv_data.a,
        &gesv_data.lda,
        gesv_data.ipiv,
        gesv_data.b,
        &gesv_data.ldb,
        &gesv_data.info
    );
}

inline void call_gesv(gesv_data_t<double>& gesv_data) {
    BLAS_FUNC(dgesv)(
        &gesv_data.n,
        &gesv_data.nrhs,
        gesv_data.a,
        &gesv_data.lda,
        gesv_data.ipiv,
        gesv_data.b,
        &gesv_data.ldb,
        &gesv_data.info
    );
}

inline void call_gesv(gesv_data_t<npy_cfloat>& gesv_data) {
    BLAS_FUNC(cgesv)(
        &gesv_data.n,
        &gesv_data.nrhs,
        gesv_data.a,
        &gesv_data.lda,
        gesv_data.ipiv,
        gesv_data.b,
        &gesv_data.ldb,
        &gesv_data.info
    );
}

inline void call_gesv(gesv_data_t<npy_cdouble>& gesv_data) {
    BLAS_FUNC(zgesv)(
        &gesv_data.n,
        &gesv_data.nrhs,
        gesv_data.a,
        &gesv_data.lda,
        gesv_data.ipiv,
        gesv_data.b,
        &gesv_data.ldb,
        &gesv_data.info
    );
}


/*
 * Hold the GETRF related variables, handle allocation/deallocation.
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
        m(m_), n(n_), lda(n_), info(-101)
    {
        a = (T *)malloc(m*n*sizeof(T));
        ipiv = (fortran_int *)malloc(n*sizeof(fortran_int));
        if ((a == NULL) || (ipiv == NULL)) {
            PyErr_NoMemory();
            info = -1;
        }
        info = 0;
    };

    ~getrf_data_t() {
        free(a);
        free(ipiv);
    };
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
 * Grab a real part of a possibly complex array.
 * This is for the work queries.
 * There must be a better way, I'm sure.
 * XXX: move together with numeric_limits etc
 */
inline float real_part(float value){ return value; }
inline double real_part(double value){ return value; }
inline float real_part(npy_cfloat value){ return npy_crealf(value); }
inline double real_part(npy_cdouble value){return npy_creal(value); }


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

    getri_data_t(getrf_data_t<T>& data) {
        assert(data.m == data.n);
        n = data.n;
        a = data.a;
        lda = data.lda;
        ipiv = data.ipiv;

        /*
         * Workspace query.
         */
        lwork = -1;
        work = (T *)malloc(10*sizeof(T));
        if (work == NULL) {
            PyErr_NoMemory();
        }

        call_getri(*this);
        if (info != 0) {
            PyErr_SetString(PyExc_ValueError, "?getri: lwork allocation failed.");
        }

        // Grab the value of `lwork`. The factor of 1.01 is from
        // https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/_basic.py#L1154
        lwork = (fortran_int)(1.01 * real_part(work[0]));
        free(work);

        // Finally, allocate
        work = (T *)malloc(lwork*sizeof(T));
        if (work == NULL) {
            PyErr_NoMemory();
            info = -101;
        }
    }

    ~getri_data_t() {
        free(work);
    };
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
