/* -*- c -*- */

/*
 *****************************************************************************
 **                            INCLUDES                                     **
 *****************************************************************************
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "../_build_utils/src/npy_cblas.h"

#include <cstddef>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <utility>

#include <iostream>


/*
 *****************************************************************************
 *                    BLAS/LAPACK calling macros                             *
 *****************************************************************************
 */

typedef CBLAS_INT         fortran_int;

typedef struct { float r, i; } f2c_complex;
typedef struct { double r, i; } f2c_doublecomplex;
/* typedef long int (*L_fp)(); */

typedef float             fortran_real;
typedef double            fortran_doublereal;
typedef f2c_complex       fortran_complex;
typedef f2c_doublecomplex fortran_doublecomplex;


extern "C" fortran_int
BLAS_FUNC(sgesv)(fortran_int *n, fortran_int *nrhs,
             float a[], fortran_int *lda,
             fortran_int ipiv[],
             float b[], fortran_int *ldb,
             fortran_int *info);
extern "C" fortran_int
BLAS_FUNC(dgesv)(fortran_int *n, fortran_int *nrhs,
             double a[], fortran_int *lda,
             fortran_int ipiv[],
             double b[], fortran_int *ldb,
             fortran_int *info);
extern "C" fortran_int
BLAS_FUNC(cgesv)(fortran_int *n, fortran_int *nrhs,
             f2c_complex a[], fortran_int *lda,
             fortran_int ipiv[],
             f2c_complex b[], fortran_int *ldb,
             fortran_int *info);
extern "C" fortran_int
BLAS_FUNC(zgesv)(fortran_int *n, fortran_int *nrhs,
             f2c_doublecomplex a[], fortran_int *lda,
             fortran_int ipiv[],
             f2c_doublecomplex b[], fortran_int *ldb,
             fortran_int *info);


/*
 *****************************************************************************
 **                      Some handy constants                               **
 *****************************************************************************
 */

template<typename T>
struct numeric_limits;

template<>
struct numeric_limits<float> {
static constexpr float one = 1.0f;
};
constexpr float numeric_limits<float>::one;

template<>
struct numeric_limits<double> {
static constexpr double one = 1.0;
};
constexpr double numeric_limits<double>::one;

template<>
struct numeric_limits<npy_cfloat> {
static constexpr npy_cfloat one = {1.0f};
};
constexpr npy_cfloat numeric_limits<npy_cfloat>::one;

template<>
struct numeric_limits<f2c_complex> {
static constexpr f2c_complex one = {1.0f, 0.0f};
};

template<>
struct numeric_limits<npy_cdouble> {
static constexpr npy_cdouble one = {1.0};
};
constexpr npy_cdouble numeric_limits<npy_cdouble>::one;

template<>
struct numeric_limits<f2c_doublecomplex> {
static constexpr f2c_doublecomplex one = {1.0};
};
constexpr f2c_doublecomplex numeric_limits<f2c_doublecomplex>::one;


/*
 *****************************************************************************
 **             rearranging of 2D matrices using blas                       **
 *****************************************************************************
 */
               /* identity square matrix generation */
template<typename typ>
static inline void
identity_matrix(typ *matrix, size_t n)
{
    size_t i;
    /* in IEEE floating point, zeroes are represented as bitwise 0 */
    memset((void *)matrix, 0, n*n*sizeof(typ));

    for (i = 0; i < n; ++i)
    {
        *matrix = numeric_limits<typ>::one;
        matrix += n+1;
    }
}


/* -------------------------------------------------------------------------- */
                  /* Solve family (includes inv) */

template<typename typ>
struct GESV_PARAMS_t
{
    typ *A; /* A is (N, N) of base type */
    typ *B; /* B is (N, NRHS) of base type */
    fortran_int * IPIV; /* IPIV is (N) */

    fortran_int N;
    fortran_int NRHS;
    fortran_int LDA;
    fortran_int LDB;
};


static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_real> *params)
{
    fortran_int rv;
    BLAS_FUNC(sgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}

static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_doublereal> *params)
{
    fortran_int rv;
    BLAS_FUNC(dgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}

static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_complex> *params)
{
    fortran_int rv;
    BLAS_FUNC(cgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}

static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_doublecomplex> *params)
{
    fortran_int rv;
    BLAS_FUNC(zgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}


/*
 * Initialize the parameters to use in for the lapack function _heev
 * Handles buffer allocation
 */
template<typename ftyp>
static inline int
init_gesv(GESV_PARAMS_t<ftyp> *params, fortran_int N, fortran_int NRHS)
{
    //npy_uint8 *mem_buff = NULL;
    //npy_uint8 *a, *b, *ipiv;
    size_t safe_N = N;
    size_t safe_NRHS = NRHS;
    fortran_int ld = std::max(N, 1);
    //mem_buff = (npy_uint8 *)malloc(safe_N * safe_N * sizeof(ftyp) +
    //                  safe_N * safe_NRHS*sizeof(ftyp) +
    //                 safe_N * sizeof(fortran_int));

    params->A = (ftyp *)malloc(safe_N * safe_N * sizeof(ftyp));
    params->B = (ftyp *)malloc(safe_N * safe_NRHS * sizeof(ftyp));
    params->IPIV = (fortran_int *)malloc(safe_N * sizeof(fortran_int));

    if (!params->A || !params->B || !params->IPIV) {
        goto error;
    }
  //  a = mem_buff;
  //  b = a + safe_N * safe_N * sizeof(ftyp);
  //  ipiv = b + safe_N * safe_NRHS * sizeof(ftyp);

  //  params->A = (ftyp*)a;
  //  params->B = (ftyp*)b;
  //  params->IPIV = (fortran_int*)ipiv;
    params->N = N;
    params->NRHS = NRHS;
    params->LDA = ld;
    params->LDB = ld;

    return 1;
 error:
    free(params->A);
    free(params->B);
    free(params->IPIV);
    memset(params, 0, sizeof(*params));

    return 0;
}


template<typename ftyp>
static inline void
release_gesv(GESV_PARAMS_t<ftyp> *params)
{
    free(params->A);
    //free(params->B); // batched_inv returns B
    free(params->IPIV);
    //memset(params, 0, sizeof(*params));
}



/*
 *****************************************************************************
 *                    Python side                                            *
 *****************************************************************************
 */

static PyObject*
py_batched_inv(PyObject* self, PyObject *args) {

    PyObject *py_arr = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_arr)) {
        return NULL;
    }

    if (!PyArray_CheckExact(py_arr)) {
        // python side has done asarray already
        return NULL;
    }

    // now we know it's an array
    PyArrayObject *arr = (PyArrayObject *)py_arr;

    // Allocate the output buffer
    npy_intp *arr_strides = PyArray_STRIDES(arr);
    npy_intp *arr_shape = PyArray_SHAPE(arr);
    int ndim = PyArray_NDIM(arr);
    PyArrayObject *ainv = (PyArrayObject *)PyArray_SimpleNew(ndim, arr_shape, NPY_DOUBLE);

    assert(ndims == 3);     // assumes the batch dimensions have been reshaped away on the python side

    // loop over the subarrays holding square matrices  
    npy_intp n_batch = arr_shape[0];  // the batch dimension
    for (npy_intp j=0; j < n_batch; j++) {

        GESV_PARAMS_t<double> params; // dgesv only
        fortran_int n = (fortran_int)arr_shape[ndim-1]; // arr.shape = (..., n, n)
        if (!init_gesv(&params, n, n)){
            // memory error
            release_gesv(&params);
            return PyErr_NoMemory();
        }

        // copy input arr : FIXME: C/F-ordering, strides
        memcpy(params.A,
               (char *)PyArray_DATA(arr) + j*arr_strides[0],
               n * n * sizeof(double));

        // copy the identity matrix for b
        identity_matrix((double*)params.B, n);

        fortran_int info;
        BLAS_FUNC(dgesv)(&params.N, &params.NRHS,
              params.A, &params.LDA,
              params.IPIV,
              params.B,
              &params.LDB,
              &info);  

        if (info != 0) {
            PyErr_SetString(PyExc_ValueError, "info != 0");
            return NULL;
        }

        // copy B into the Ainv buffer: B is 'behaved', can memcpy
        memcpy((char *)PyArray_DATA(ainv) + j*n*n*sizeof(double),
               params.B,
               n*n*sizeof(double));

        // clean up
        release_gesv(&params);

    }
    return (PyObject *)ainv;
}

/////////////////////////////////////

static PyMethodDef _batched_lapack_methods[] = {
    //...
    {"inv", py_batched_inv, METH_VARARGS, 
     "batched inv"},
    //...
    {NULL, NULL, 0, NULL}        /* Sentinel */
};



static struct PyModuleDef _lapackmodule = {
    PyModuleDef_HEAD_INIT,
    "_batched_lapack",   /* name of module */
    NULL, //spam_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    _batched_lapack_methods
};


PyMODINIT_FUNC
PyInit__batched_lapack(void)
{
    import_array();
    return PyModule_Create(&_lapackmodule);
}



