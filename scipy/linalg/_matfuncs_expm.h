#ifndef _MATFUNCS_EXPM_H
#define _MATFUNCS_EXPM_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <math.h>
#include "numpy/arrayobject.h"

#if defined(_MSC_VER)
    // MSVC nonsense
    #include <complex.h>
    #define EXPM_Z _Dcomplex
    #define EXPM_C _Fcomplex
#else
    // C99-compliant compilers
    #include <complex.h>
    #define EXPM_Z double complex
    #define EXPM_C float complex
#endif


// BLAS and LAPACK functions used
void saxpy_(int* n, float* sa, float* sx, int* incx, float* sy, int* incy);
void scopy_(int* n, float* dx, int* incx, float* dy, int* incy);
void sgemm_(char* transa, char* transb, int* m, int* n, int* k, float* alpha, float* a, int* lda, float* b, int* ldb, float* beta, float* c, int* ldc);
void sgemv_(char* trans, int* m, int* n, float* alpha, float* a, int* lda, float* x, int* incx, float* beta, float* y, int* incy);
void sgetrf_(int* m, int* n, float* a, int* lda, int* ipiv, int* info);
void sgetrs_(char* trans, int* n, int* nrhs, float* a, int* lda, int* ipiv, float* b, int* ldb, int* info);
void slacn2_(int* n, float* v, float* x, int* isgn, float* est, int* kase, int* isave);
void sscal_(int* n, float* sa, float* sx, int* incx);

void daxpy_(int* n, double* sa, double* sx, int* incx, double* sy, int* incy);
void dcopy_(int* n, double* dx, int* incx, double* dy, int* incy);
void dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* a, int* lda, double* b, int* ldb, double* beta, double* c, int* ldc);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
void dgetrs_(char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
void dlacn2_(int* n, double* v, double* x, int* isgn, double* est, int* kase, int* isave);
void dscal_(int* n, double* sa, double* sx, int* incx);

void caxpy_(int* n, EXPM_C* sa, EXPM_C* sx, int* incx, EXPM_C* sy, int* incy);
void ccopy_(int* n, EXPM_C* dx, int* incx, EXPM_C* dy, int* incy);
void cgemm_(char* transa, char* transb, int* m, int* n, int* k, EXPM_C* alpha, EXPM_C* a, int* lda, EXPM_C* b, int* ldb, EXPM_C* beta, EXPM_C* c, int* ldc);
void cgemv_(char* trans, int* m, int* n, EXPM_C* alpha, EXPM_C* a, int* lda, EXPM_C* x, int* incx, EXPM_C* beta, EXPM_C* y, int* incy);
void cgetrf_(int* m, int* n, EXPM_C* a, int* lda, int* ipiv, int* info);
void cgetrs_(char* trans, int* n, int* nrhs, EXPM_C* a, int* lda, int* ipiv, EXPM_C* b, int* ldb, int* info);
void clacn2_(int* n, EXPM_C* v, EXPM_C* x, float* est, int* kase, int* isave);
void csscal_(int* n, float* sa, EXPM_C* sx, int* incx);

void zaxpy_(int* n, EXPM_Z* sa, EXPM_Z* sx, int* incx, EXPM_Z* sy, int* incy);
void zcopy_(int* n, EXPM_Z* dx, int* incx, EXPM_Z* dy, int* incy);
void zgemm_(char* transa, char* transb, int* m, int* n, int* k, EXPM_Z* alpha, EXPM_Z* a, int* lda, EXPM_Z* b, int* ldb, EXPM_Z* beta, EXPM_Z* c, int* ldc);
void zgemv_(char* trans, int* m, int* n, EXPM_Z* alpha, EXPM_Z* a, int* lda, EXPM_Z* x, int* incx, EXPM_Z* beta, EXPM_Z* y, int* incy);
void zgetrf_(int* m, int* n, EXPM_Z* a, int* lda, int* ipiv, int* info);
void zgetrs_(char* trans, int* n, int* nrhs, EXPM_Z* a, int* lda, int* ipiv, EXPM_Z* b, int* ldb, int* info);
void zlacn2_(int* n, EXPM_Z* v, EXPM_Z* x, double* est, int* kase, int* isave);
void zdscal_(int* n, double* sa, EXPM_Z* sx, int* incx);


void pick_pade_structure_s(float* Am, const Py_ssize_t n, int* m, int* s);
void pick_pade_structure_d(double* Am, const Py_ssize_t n, int* m, int* s);
void pick_pade_structure_c(EXPM_C* Am, const Py_ssize_t n, int* m, int* s);
void pick_pade_structure_z(EXPM_Z* Am, const Py_ssize_t n, int* m, int* s);

void pade_UV_calc_s(float* Am, const Py_ssize_t n, const int m, int* info);
void pade_UV_calc_d(double* Am, const Py_ssize_t n, const int m, int* info);
void pade_UV_calc_c(EXPM_C* Am, const Py_ssize_t n, const int m, int* info);
void pade_UV_calc_z(EXPM_Z* Am, const Py_ssize_t n, const int m, int* info);


#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}

static PyObject* expm_error;

static char doc_pps[] = "[m, s] = pick_pade_structure(Am)";

static PyObject*
pick_pade_structure(PyObject *dummy, PyObject *args) {
    Py_ssize_t n, n1;
    int m = 0, s = 0;
    PyArrayObject *ap_Am=NULL;

    if (!PyArg_ParseTuple(args, ("O!"), &PyArray_Type, (PyObject **)&ap_Am)) return NULL;

    if (!(PyArray_IS_C_CONTIGUOUS(ap_Am)) ||
        ((PyArray_TYPE(ap_Am) != NPY_FLOAT64) &&
         (PyArray_TYPE(ap_Am) != NPY_FLOAT32) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX64) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX128)) ||
        (PyArray_NDIM(ap_Am) != 3)
        )
    {
        PYERR(expm_error, "Input must be a 3D C-contiguous array with size"
                          " (5, n, n) that is of type float32, float64,"
                          " complex64, or complex128.");
    }

    n = PyArray_DIMS(ap_Am)[2];
    n1 = PyArray_DIMS(ap_Am)[1];
    if (n != n1)
    {
        PYERR(expm_error, "Last two dimensions of the input must be the same.");
    }
    // Make the call based on dtype
    switch (PyArray_TYPE(ap_Am))
    {
        case (NPY_FLOAT32):
        {
            float* Ams = (float*)PyArray_DATA(ap_Am);
            pick_pade_structure_s(Ams, n, &m, &s);
            break;
        }
        case (NPY_FLOAT64):
        {
            double* Amd = (double*)PyArray_DATA(ap_Am);
            pick_pade_structure_d(Amd, n, &m, &s);
            break;
        }
        case (NPY_COMPLEX64):
        {
            EXPM_C* Amc = (EXPM_C*)PyArray_DATA(ap_Am);
            pick_pade_structure_c(Amc, n, &m, &s);
            break;
        }
        case (NPY_COMPLEX128):
        {
            EXPM_Z* Amz = (EXPM_Z*)PyArray_DATA(ap_Am);
            pick_pade_structure_z(Amz, n, &m, &s);
            break;
        }
    }

    return Py_BuildValue("ii", m, s);
}


static char doc_puv[] = "info = pade_UV_calc(Am, m)";

static PyObject*
pade_UV_calc(PyObject *dummy, PyObject *args) {
    Py_ssize_t n, n1;
    int m = 0, info = 0;
    PyArrayObject *ap_Am=NULL;

    if (!PyArg_ParseTuple(args, ("O!i"), &PyArray_Type, (PyObject **)&ap_Am, &m)) return NULL;
    if (!(PyArray_IS_C_CONTIGUOUS(ap_Am)) ||
        ((PyArray_TYPE(ap_Am) != NPY_FLOAT64) &&
         (PyArray_TYPE(ap_Am) != NPY_FLOAT32) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX64) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX128)) ||
        (PyArray_NDIM(ap_Am) != 3)
        )
    {
        PYERR(expm_error, "Input must be a 3D C-contiguous array with size"
                          " (5, n, n) that is of type float32, float64,"
                          " complex64, or complex128.");
    }

    n = PyArray_DIMS(ap_Am)[2];
    n1 = PyArray_DIMS(ap_Am)[1];
    if (n != n1)
    {
        PYERR(expm_error, "Last two dimensions of the input must be the same.");
    }
    // Make the call based on dtype
    switch (PyArray_TYPE(ap_Am))
    {
        case (NPY_FLOAT32):
        {
            float* Ams = (float*)PyArray_DATA(ap_Am);
            pade_UV_calc_s(Ams, n, m, &info);
            break;
        }
        case (NPY_FLOAT64):
        {
            double* Amd = (double*)PyArray_DATA(ap_Am);
            pade_UV_calc_d(Amd, n, m, &info);
            break;
        }
        case (NPY_COMPLEX64):
        {
            EXPM_C* Amc = (EXPM_C*)PyArray_DATA(ap_Am);
            pade_UV_calc_c(Amc, n, m, &info);
            break;
        }
        case (NPY_COMPLEX128):
        {
            EXPM_Z* Amz = (EXPM_Z*)PyArray_DATA(ap_Am);
            pade_UV_calc_z(Amz, n, m, &info);
            break;
        }
    }

    return Py_BuildValue("i", info);
}

static struct PyMethodDef expm_module_methods[] = {
  {"pick_pade_structure", pick_pade_structure, METH_VARARGS, doc_pps},
  {"pade_UV_calc"       , pade_UV_calc       , METH_VARARGS, doc_puv},
  {NULL                 , NULL               , 0           , NULL   }
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_matfuncs_expm",
    NULL,
    -1,
    expm_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__matfuncs_expm(void)
{
    PyObject *module, *mdict;

    import_array();

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }

    mdict = PyModule_GetDict(module);
    if (mdict == NULL) {
        return NULL;
    }
    expm_error = PyErr_NewException("_matfuncs_expm.error", NULL, NULL);
    if (expm_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", expm_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}

#undef EXPM_C
#undef EXPM_Z

#endif // _MATFUNCS_EXPM_H
