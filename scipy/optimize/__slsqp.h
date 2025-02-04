#ifndef __SLSQPLIB_H
#define __SLSQPLIB_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}
static PyObject* slsqp_error;

#include <math.h>
#include "__nnls.h"

void daxpy_(int* n, double* sa, double* sx, int* incx, double* sy, int* incy);
double ddot_(int* n, double* dx, int* incx, double* dy, int* incy);
void dgemv_(char* trans, int* m, int* n, double* alpha, double* a, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
void dgeqr2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);
void dgeqrf_(int* m, int* n, double* a, int* lda, double* tau, double* work, double* lwork, int* info);
void dgerq2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);
void dlarf_(char* side, int* m, int* n, double* v, int* incv, double* tau, double* c, int* ldc, double* work);
void dlarfgp_(int* n, double* alpha, double* x, int* incx, double* tau);
void dlartgp_(double* f, double* g, double* cs, double* sn, double* r);
double dnrm2_(int* n, double* x, int* incx);
void dorm2r_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* info);
void dormr2_(char* side, char* trans, int* m, int* n, int* k, double* a, int* lda, double* tau, double* c, int* ldc, double* work, int* info);
void dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);
void dtrsv_(char* uplo, char* trans, char* diag, int* n, double* a, int* lda, double* x, int* incx);

static void ldp(int m, int n, double* g, double* h, double* x, double* buffer, int* indices, double* xnorm, int* mode);
static void lsi(int ma, int mg, int n, double* a, double* b, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);


static PyObject*
nnls(PyObject *dummy, PyObject *args) {

    int maxiter, info = 0;
    PyArrayObject* ap_A=NULL;
    PyArrayObject* ap_b=NULL;
    double* buffer;
    double rnorm;

    // Get the input array
    if (!PyArg_ParseTuple(args,
                         ("O!O!i"),
                         &PyArray_Type, (PyObject **)&ap_A,
                         &PyArray_Type, (PyObject **)&ap_b,
                         &maxiter)
        )
    {
        return NULL;
    }

    // Check for dtype compatibility
    if ((PyArray_TYPE(ap_A) != NPY_FLOAT64) || (PyArray_TYPE(ap_b) != NPY_FLOAT64))
    {
        PYERR(slsqp_error, "Inputs to nnls must be of type numpy.float64.");
    }

    int ndim = PyArray_NDIM(ap_A);              // Number of dimensions
    if (ndim != 2)
    {
        PYERR(slsqp_error, "Input array A must be 2D.");
    }
    npy_intp* shape = PyArray_SHAPE(ap_A);       // Array shape
    npy_intp m = shape[0];                       // Number of rows
    npy_intp n = shape[1];                       // Number of columns

    int ndim_b = PyArray_NDIM(ap_b);             // Number of dimensions
    npy_intp* shape_b = PyArray_SHAPE(ap_b);     // Array shape
    if (ndim_b == 1)
    {
        if (shape_b[0] != m)
        {
            PYERR(slsqp_error, "Input array b must have the same number of rows as A.");
        }
    } else if (ndim_b == 2) {
        if (shape_b[0] != m)
        {
            PYERR(slsqp_error, "Input array b must have the same number of rows as A.");
        }
        if (shape_b[1] != 1)
        {
            PYERR(slsqp_error, "Input array b must have only one column.");
        }
    } else {
        PYERR(slsqp_error, "Input array b must be 1D or 2D with one column.");
    }

    // Allocate memory for the algorithm,
    // A is m x n, b is m, x is n, w is n, zz is m, work is n
    // total m*(n+2) + 3*n
    //indices is n
    buffer = malloc((m*(n+2) + 3*n)*sizeof(double));
    if (buffer == NULL)
    {
        PYERR(slsqp_error, "Memory allocation failed.");
    }
    int *indices = malloc(n*sizeof(int));
    if (indices == NULL)
    {
        free(buffer);
        PYERR(slsqp_error, "Memory allocation failed.");
    }

    double* x = &buffer[0];                 // Solution vector x (n)
    double* a = &buffer[n];                 // Matrix A (m x n)
    double* b = &buffer[n*m + n];           // Vector b (m)
    double* w = &buffer[(n+1)*m + n];       // Vector w (n)
    double* zz = &buffer[(n+1)*m + 2*n];    // Vector zz (m)
    double* work = &buffer[(n+2)*m + 2*n];  // Work vector (n)

    npy_intp* restrict strides = PyArray_STRIDES(ap_A);
    double* restrict data_A = (double *)PyArray_DATA(ap_A);
    npy_intp* restrict stride_b = PyArray_STRIDES(ap_b);
    // If b is 2D then pick the stride of the first dimension
    npy_intp rc_stride = (ndim_b == 1 ? stride_b[0] : stride_b[1]);
    double* restrict data_b = (double *)PyArray_DATA(ap_b);

    // Copy the data from the numpy array
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            a[i + j*m] = data_A[(j*strides[1] + i*strides[0])/sizeof(double)];
        }
    }
    for (int i = 0; i < m; i++)
    {
        b[i] = data_b[(i * rc_stride)/sizeof(double)];
    }

    // Call nnls
    nonnegative_lsq_imp((int)m, (int)n, a, b, x, w, zz, work, indices, maxiter, &rnorm, &info);
    // x is the first n elements of buffer, shrink buffer to n elements
    free(indices);
    double* mem_ret = realloc(buffer, n*sizeof(double));
    // Very unlikely, but just in case
    if (mem_ret == NULL)
    {
        free(buffer);
        PYERR(slsqp_error, "Memory reallocation failed.");
    }
    npy_intp shape_ret[1] = {n};
    PyArrayObject* ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(1, shape_ret, NPY_FLOAT64, mem_ret);
    // Return the result
    return Py_BuildValue("Ndi",PyArray_Return(ap_ret), rnorm, info);

}


static PyObject*
ldp_wrapper(PyObject *dummy, PyObject *args)
{
    npy_intp m, n;
    int mode;
    double xnorm;
    PyArrayObject *ap_g=NULL, *ap_h=NULL;


    if (!PyArg_ParseTuple(args, "O!O!",
                          &PyArray_Type, (PyObject **)&ap_g,
                          &PyArray_Type, (PyObject **)&ap_h))
    {
        return NULL;
    }

    if ((PyArray_TYPE(ap_g) != NPY_FLOAT64) || (PyArray_TYPE(ap_h) != NPY_FLOAT64))
    {
        PYERR(slsqp_error, "Inputs to ldp must be of type numpy.float64.");
    }

    int ndim_g = PyArray_NDIM(ap_g);
    if (ndim_g != 2)
    {
        PYERR(slsqp_error, "Input array g must be 2D.");
    }
    npy_intp* shape_g = PyArray_SHAPE(ap_g);
    m = shape_g[0];
    n = shape_g[1];

    int ndim_h = PyArray_NDIM(ap_h);
    if (ndim_h != 1)
    {
        PYERR(slsqp_error, "Input array h must be 1D.");
    }
    npy_intp* shape_h = PyArray_SHAPE(ap_h);
    if (shape_h[0] != m)
    {
        PYERR(slsqp_error, "Input array h must have the same shape as (m,).");
    }
    // Over allocate indices a bit for transposing problems and other uses
    int* indices = malloc((m+n)*sizeof(int));
    if (indices == NULL) { PYERR(slsqp_error, "Memory allocation for indices failed."); }

    // Allocate memory for the problem data and the algorithm.
    double* mem_ret = malloc(((m*n + m) + (m+2)*(n+1) + 3*m + n)*sizeof(double));
    if (mem_ret == NULL)
    {
        free(indices);
        PYERR(slsqp_error, "Memory allocation for buffer failed.");
    }
    double* restrict x = &mem_ret[0];
    double* restrict g = &mem_ret[n];
    double* restrict h = &mem_ret[n + m*n];
    double* restrict buffer = &mem_ret[(m*n + m + n)];

    // Copy the data from the numpy arrays
    double* data_g = (double*)PyArray_DATA(ap_g);
    double* data_h = (double*)PyArray_DATA(ap_h);
    npy_intp* restrict strides_g = PyArray_STRIDES(ap_g);
    npy_intp* restrict stride_h = PyArray_STRIDES(ap_h);
    npy_intp row_stride = (strides_g[0]/sizeof(double));
    npy_intp col_stride = (strides_g[1]/sizeof(double));

    // Copy the data from the numpy arrays in Fortran order
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            g[i + j*m] = data_g[i*row_stride + j*col_stride];
        }
    }
    row_stride = (stride_h[0]/sizeof(double));
    for (int i = 0; i < m; i++) { h[i] = data_h[i*row_stride]; }

    ldp((int)m, (int)n, g, h, x, buffer, indices, &xnorm, &mode);
    free(indices);
    // Truncate it to the size of the solution vector x.
    double* mem_x = realloc(mem_ret, n*sizeof(double));
    if (mem_x == NULL) { PYERR(slsqp_error, "Memory reallocation failed."); }
    npy_intp shape_x[1] = {n};
    PyArrayObject* ap_x = (PyArrayObject*)PyArray_SimpleNewFromData(1, shape_x, NPY_FLOAT64, mem_x);

    return Py_BuildValue("N", PyArray_Return(ap_x));
}


static PyObject*
lsi_wrapper(PyObject *dummy, PyObject *args)
{
    npy_intp ma, mg, n;
    int mode;
    double xnorm;
    PyArrayObject *ap_a=NULL, *ap_b=NULL, *ap_g=NULL, *ap_h=NULL;

    if (!PyArg_ParseTuple(args, "O!O!O!O!",
                          &PyArray_Type, (PyObject **)&ap_a,
                          &PyArray_Type, (PyObject **)&ap_b,
                          &PyArray_Type, (PyObject **)&ap_g,
                          &PyArray_Type, (PyObject **)&ap_h))
    {
        return NULL;
    }

    if ((PyArray_TYPE(ap_a) != NPY_FLOAT64) || (PyArray_TYPE(ap_b) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_g) != NPY_FLOAT64) || (PyArray_TYPE(ap_h) != NPY_FLOAT64))
    {
        PYERR(slsqp_error, "Inputs to lsi must be of type numpy.float64.");
    }

    int ndim_a = PyArray_NDIM(ap_a);
    if (ndim_a != 2)
    {
        PYERR(slsqp_error, "Input array a must be 2D.");
    }
    npy_intp* shape_a = PyArray_SHAPE(ap_a);
    ma = shape_a[0];
    n = shape_a[1];

    int ndim_b = PyArray_NDIM(ap_b);
    if (ndim_b != 1)
    {
        PYERR(slsqp_error, "Input array b must be 1D.");
    }
    npy_intp* shape_b = PyArray_SHAPE(ap_b);
    if (shape_b[0] != ma)
    {
        PYERR(slsqp_error, "Input array b must have the same shape as (ma,).");
    }

    int ndim_g = PyArray_NDIM(ap_g);
    if (ndim_g != 2)
    {
        PYERR(slsqp_error, "Input array g must be 2D.");
    }
    npy_intp* shape_g = PyArray_SHAPE(ap_g);
    mg = shape_g[0];
    if (shape_g[1] != n)
    {
        PYERR(slsqp_error, "Input array g must have the same number of columns as a.");
    }

    int ndim_h = PyArray_NDIM(ap_h);
    if (ndim_h != 1)
    {
        PYERR(slsqp_error, "Input array h must be 1D.");
    }
    npy_intp* shape_h = PyArray_SHAPE(ap_h);
    if (shape_h[0] != mg)
    {
        PYERR(slsqp_error, "Input array h must have the same shape as (mg,).");
    }

    // Over allocate indices a bit for transposing problems and other uses
    int* indices = malloc((ma+n)*sizeof(int));
    if (indices == NULL) { PYERR(slsqp_error, "Memory allocation for indices failed."); }

    // Allocate memory for the problem data and the algorithm.
    double* mem_ret = malloc(((ma*n + ma) + (mg*n + mg) + (mg+2)*(n+1) + 3*mg + n)*sizeof(double));
    if (mem_ret == NULL)
    {
        free(indices);
        PYERR(slsqp_error, "Memory allocation for buffer failed.");
    }
    double* restrict x = &mem_ret[0];
    double* restrict a = &mem_ret[n];
    double* restrict b = &mem_ret[n + ma*n];
    double* restrict g = &mem_ret[n + ma*n + ma];
    double* restrict h = &mem_ret[n + ma*n + ma + mg*n];
    double* restrict buffer = &mem_ret[n + ma*n + ma + mg*n + mg];

    // Copy the data from the numpy arrays
    double* data_a = (double*)PyArray_DATA(ap_a);
    double* data_b = (double*)PyArray_DATA(ap_b);
    double* data_g = (double*)PyArray_DATA(ap_g);
    double* data_h = (double*)PyArray_DATA(ap_h);
    npy_intp* restrict strides_a = PyArray_STRIDES(ap_a);
    npy_intp* restrict stride_b = PyArray_STRIDES(ap_b);
    npy_intp* restrict strides_g = PyArray_STRIDES(ap_g);
    npy_intp* restrict stride_h = PyArray_STRIDES(ap_h);
    npy_intp row_stride = (strides_a[0]/sizeof(double));
    npy_intp col_stride = (strides_a[1]/sizeof(double));

    // Copy the data from the numpy arrays in Fortran order
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < ma; i++) {
            a[i + j*ma] = data_a[i*row_stride + j*col_stride];
        }
    }
    row_stride = (stride_b[0]/sizeof(double));
    for (int i = 0; i < ma; i++) { b[i] = data_b[i*row_stride]; }

    row_stride = (strides_g[0]/sizeof(double));
    col_stride = (strides_g[1]/sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < mg; i++) {
            g[i + j*mg] = data_g[i*row_stride + j*col_stride];
        }
    }
    row_stride = (stride_h[0]/sizeof(double));
    for (int i = 0; i < mg; i++) { h[i] = data_h[i*row_stride]; }

    lsi((int)ma, (int)mg, (int)n, a, b, g, h, x, buffer, indices, &xnorm, &mode);
    printf("mode: %d\n", mode);
    printf("xnorm: %f\n", xnorm);
    free(indices);
    // Truncate it to the size of the solution vector x.
    double* mem_x = realloc(mem_ret, n*sizeof(double));
    if (mem_x == NULL) { PYERR(slsqp_error, "Memory reallocation failed."); }
    npy_intp shape_x[1] = {n};
    PyArrayObject* ap_x = (PyArrayObject*)PyArray_SimpleNewFromData(1, shape_x, NPY_FLOAT64, mem_x);

    return Py_BuildValue("N", PyArray_Return(ap_x));
}


static char doc_nnls[] = ("Compute the nonnegative least squares solution.\n\n"
                           "    x, info = nnls(A)\n\n");
static char doc_lsi_wrapper[] = ("Compute the least squares solution subject to "
                                 "inequality constraints.\n\n"
                                 "    x = lsi_wrapper(A, b, G, h)\n\n");

static char doc_ldp_wrapper[] = ("Compute the least distance solution.\n\n"
                                 "    x = ldp_wrapper(G, h)\n\n");

static struct PyMethodDef slsqplib_module_methods[] = {
  {"nnls", nnls, METH_VARARGS, doc_nnls},
  {"ldp_wrapper", ldp_wrapper, METH_VARARGS, doc_ldp_wrapper},
  {"lsi_wrapper", lsi_wrapper, METH_VARARGS, doc_lsi_wrapper},
  {NULL, NULL, 0, NULL}
};

struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_slsqplib",
    NULL,
    -1,
    slsqplib_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__slsqplib(void)
{
    import_array();

    PyObject* module = PyModule_Create(&moduledef);
    if (module == NULL) { return NULL; }
    PyObject* mdict = PyModule_GetDict(module);
    if (mdict == NULL) { return NULL; }
    slsqp_error = PyErr_NewException("_slsqplib.error", NULL, NULL);
    if (slsqp_error == NULL) { return NULL; }
    if (PyDict_SetItemString(mdict, "error", slsqp_error)) { return NULL; }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}


#endif // __SLSQPLIB_H