#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "../_build_utils/src/npy_cblas.h"
#include "../_build_utils/src/fortran_defs.h"
#include <math.h>

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}
#define DAXPY BLAS_FUNC(daxpy)
#define DCOPY BLAS_FUNC(dcopy)
#define DDOT BLAS_FUNC(ddot)
#define DNRM2 BLAS_FUNC(dnrm2)
#define DSCAL BLAS_FUNC(dscal)
#define DPOTRF BLAS_FUNC(dpotrf)
#define DTRTRS BLAS_FUNC(dtrtrs)


static void DAXPY(const int* n, const double* da, const double* dx,
                  const int* incx, double* dy, const int* incy);
static void DCOPY(const int* n, const double* dx, const int *incx,
                  double *dy, const int *incy);
static double DDOT(const int *n, const double *dx, const int *incx,
                   const double *dy, const int *incy);
static double DNRM2(const int *n, const double *x, const int *incx);
static void DSCAL(const int *n, const double *da, double *dx,
                  const int *incx);
static void DPORTF(const char* uplo, const int *n, double *a, const int *lda,
                   int *info);
static void DTRTRS(const char* uplo, const char *trans, const char *diag,
                   const int *n, const int *nrhs,
                   const double *a, const int *lda, double *b, const int *ldb,
                   int *info);


static PyObject* lbfgsb_error;

static void setulb(
    const int n, const int m, double* x, double* l, double* u, int* nbd,
    double* f, double* g, const double factr, const double pgtol,
    double* wa, int* iwa, int* task, int* task_msg, const int iprint,
    int* lsave, int* isave, double* dsave, const int maxls
);

static char doc_setulb[] = "setulb";

static PyObject*
lbfgsb_setulb(PyObject *self, PyObject *args)
{
    int m, n, task, task_msg, iprint, maxls;
    double f, factr, pgtol;
    int *nbd, *iwa, *isave, *lsave;
    double *x, *l, *u, *g, *wa, *dsave;

    PyArrayObject *ap_x = NULL, *ap_l = NULL, *ap_u = NULL, *ap_g = NULL;
    PyArrayObject *ap_nbd = NULL, *ap_wa = NULL, *ap_iwa = NULL;
    PyArrayObject *ap_lsave = NULL, *ap_isave = NULL, *ap_dsave = NULL;

    // Check argument types
    // m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave,maxls
    if (!(PyArg_ParseTuple(args,
                           "iO!O!O!O!dO!ddO!O!iiiO!O!O!i",
                           &m,                                      // i
                           &PyArray_Type, (PyObject **)&ap_x,       // O!
                           &PyArray_Type, (PyObject **)&ap_l,       // O!
                           &PyArray_Type, (PyObject **)&ap_u,       // O!
                           &PyArray_Type, (PyObject **)&ap_nbd,     // O!
                           &f,                                      // d
                           &PyArray_Type, (PyObject **)&ap_g,       // O!
                           &factr,                                  // d
                           &pgtol,                                  // d
                           &PyArray_Type, (PyObject **)&ap_wa,      // O!
                           &PyArray_Type, (PyObject **)&ap_iwa,     // O!
                           &task,                                   // i
                           &task_msg,                               // i
                           &iprint,                                 // i
                           &PyArray_Type, (PyObject **)&ap_lsave,   // O!
                           &PyArray_Type, (PyObject **)&ap_isave,   // O!
                           &PyArray_Type, (PyObject **)&ap_dsave,   // O!
                           &maxls                                   // i
                           ))) { return NULL; }

    // Check if arrays are contiguous and with the right dtype.
    // All arrays are 1D hence both F and C contiguous flags are set to True.
    if (!(PyArray_IS_C_CONTIGUOUS(ap_x)) || (PyArray_TYPE(ap_x) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (x) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_l)) || (PyArray_TYPE(ap_l) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (l) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_u)) || (PyArray_TYPE(ap_u) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (u) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_wa)) || (PyArray_TYPE(ap_wa) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (wa) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_iwa)) || (PyArray_TYPE(ap_iwa) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (iwa) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_dsave)) || (PyArray_TYPE(ap_dsave) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (dsave) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_nbd)) || (PyArray_TYPE(ap_nbd) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (nbd) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_lsave)) || (PyArray_TYPE(ap_lsave) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (lsave) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_isave)) || (PyArray_TYPE(ap_isave) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (isave) must be a contiguous array of type int32.");
    }
    n = PyArray_DIMS(ap_x)[0];

    x = (double *)PyArray_DATA(ap_x);
    l = (double *)PyArray_DATA(ap_l);
    u = (double *)PyArray_DATA(ap_u);
    nbd = (int *)PyArray_DATA(ap_nbd);
    g = (double *)PyArray_DATA(ap_g);
    wa = (double *)PyArray_DATA(ap_wa);
    iwa = (int *)PyArray_DATA(ap_iwa);
    lsave = (int *)PyArray_DATA(ap_lsave);
    isave = (int *)PyArray_DATA(ap_isave);
    dsave = (double *)PyArray_DATA(ap_dsave);

    // Make the function call
    setulb(n, m, x, l, u, nbd, &f, g, factr, pgtol, wa, iwa, &task, &task_msg,
           iprint, lsave, isave, dsave, maxls);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef lbfgsb_module_methods[] = {
  {"setulb", lbfgsb_setulb, METH_VARARGS, doc_setulb},
  {NULL,     NULL,          0,            NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_lbfgsb",
    NULL,
    -1,
    lbfgsb_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__lbfgsb(void)
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
    lbfgsb_error = PyErr_NewException ("_lbfgsb.error", NULL, NULL);
    if (lbfgsb_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", lbfgsb_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}


#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
