/*
 *  ALGORITHM 733, COLLECTED ALGORITHMS FROM ACM.
 *  TRANSACTIONS ON MATHEMATICAL SOFTWARE,
 *  VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 262-281.
 *  https://doi.org/10.1145/192115.192124
 *
 *
 *  https://web.archive.org/web/20170106155705/http://permalink.gmane.org/gmane.comp.python.scientific.devel/6725
 *  ------
 *  From: Deborah Cotton <cotton@hq.acm.org>
 *  Date: Fri, 14 Sep 2007 12:35:55 -0500
 *  Subject: RE: Algorithm License requested
 *  To: Alan Isaac
 *
 *  Prof. Issac,
 *
 *  In that case, then because the author consents to [the ACM] releasing
 *  the code currently archived at http://www.netlib.org/toms/733 under the
 *  BSD license, the ACM hereby releases this code under the BSD license.
 *
 *  Regards,
 *
 *  Deborah Cotton, Copyright & Permissions
 *  ACM Publications
 *  2 Penn Plaza, Suite 701**
 *  New York, NY 10121-0701
 *  permissions@acm.org
 *  212.869.7440 ext. 652
 *  Fax. 212.869.0481
 *  ------
*/

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
void dgelsy_(int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, int* jpvt, double* rcond, int* rank, double* work, int* lwork, int* info);
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
void dtpmv_(char* uplo, char* trans, char* diag, int* n, double* ap, double* x, int* incx);
void dtpsv_(char* uplo, char* trans, char* diag, int* n, double* ap, double* x, int* incx);
void dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);
void dtrsv_(char* uplo, char* trans, char* diag, int* n, double* a, int* lda, double* x, int* incx);

// static void ldp(int m, int n, double* g, double* h, double* x, double* buffer, int* indices, double* xnorm, int* mode);
// static void lsi(int ma, int mg, int n, double* a, double* b, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);
// static void lsei(int ma, int me, int mg, int n, double* a, double* b, double* e, double* f, double* g, double* h, double* x, double* buffer, int* jw, double* xnorm, int* mode);
// static void lsq(int m, int meq, int n, int nl, double* S, double* t, double* C, double* d, double* xl, double* xu, double* x, double* y, double* buffer, int* jw, int* mode);
static void __slsqp_body(struct SLSQP_static_vars S, double* C, double* d, double* sol, double* mult, double* xl, double* xu, double* funx, double* gradx, double* buffer, int* indices);

struct SLSQP_static_vars {
    double acc;
    double alpha;
    double f0;
    double gs;
    double h1;
    double h2;
    double h3;
    double h4;
    double t;
    double t0;
    double tol;
    int exact;
    int inconsistent;
    int reset;
    int iter;
    int itermax;
    int line;
    int m;
    int meq;
    int mode;
    int n;
};


// Some helper x macros to pack and unpack the SLSQP_static_vars struct and
// the Python dictionary.

#define STRUCT_DOUBLE_FIELD_NAMES X(acc) X(alpha) X(f0) X(gs) X(h1) X(h2) X(h3) X(h4) X(t) X(t0) X(tol)
#define STRUCT_INT_FIELD_NAMES X(exact) X(inconsistent) X(reset) X(iter) X(itermax) X(line) X(meq) X(mode)
#define STRUCT_FIELD_NAMES STRUCT_INT_FIELD_NAMES STRUCT_DOUBLE_FIELD_NAMES


static PyObject*
nnls(PyObject* Py_UNUSED(dummy), PyObject* args) {

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
    // A is m x n, b is m, x is n, w is n, zz is m
    // total m*(n+2) + 2*n
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
    __nnls((int)m, (int)n, a, b, x, w, zz, indices, maxiter, &rnorm, &info);
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
ldp_wrapper(PyObject* Py_UNUSED(dummy), PyObject* args)
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
    double* mem_ret = malloc(((m*n + m) + (m+2)*(n+1) + 2*m + n)*sizeof(double));
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
lsi_wrapper(PyObject* Py_UNUSED(dummy), PyObject* args)
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
    double* mem_ret = calloc(((ma*n + ma) + (mg*n + mg) + (mg+2)*(n+1) + 2*mg + n), sizeof(double));
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

    free(indices);
    // Truncate it to the size of the solution vector x.
    double* mem_x = realloc(mem_ret, n*sizeof(double));
    if (mem_x == NULL) { PYERR(slsqp_error, "Memory reallocation failed."); }
    npy_intp shape_x[1] = {n};
    PyArrayObject* ap_x = (PyArrayObject*)PyArray_SimpleNewFromData(1, shape_x, NPY_FLOAT64, mem_x);

    return Py_BuildValue("Ndi", PyArray_Return(ap_x), xnorm, mode);
}


static PyObject*
lsei_wrapper(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    npy_intp ma, me, mg, n;
    int mode = 0;
    double xnorm =0.0;
    PyArrayObject *ap_a=NULL, *ap_b=NULL, *ap_e=NULL, *ap_f=NULL, *ap_g=NULL, *ap_h=NULL;

    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!",
                          &PyArray_Type, (PyObject **)&ap_a,
                          &PyArray_Type, (PyObject **)&ap_b,
                          &PyArray_Type, (PyObject **)&ap_e,
                          &PyArray_Type, (PyObject **)&ap_f,
                          &PyArray_Type, (PyObject **)&ap_g,
                          &PyArray_Type, (PyObject **)&ap_h))
    {
        return NULL;
    }

    if ((PyArray_TYPE(ap_a) != NPY_FLOAT64) || (PyArray_TYPE(ap_b) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_e) != NPY_FLOAT64) || (PyArray_TYPE(ap_f) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_g) != NPY_FLOAT64) || (PyArray_TYPE(ap_h) != NPY_FLOAT64))
    {
        PYERR(slsqp_error, "Inputs to lsei must be of type numpy.float64.");
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

    int ndim_e = PyArray_NDIM(ap_e);
    if (ndim_e != 2)
    {
        PYERR(slsqp_error, "Input array e must be 2D.");
    }
    npy_intp* shape_e = PyArray_SHAPE(ap_e);
    me = shape_e[0];
    if (shape_e[1] != n)
    {
        PYERR(slsqp_error, "Input array e must have the same number of columns as a.");
    }

    int ndim_f = PyArray_NDIM(ap_f);
    if (ndim_f != 1)
    {
        PYERR(slsqp_error, "Input array f must be 1D.");
    }
    npy_intp* shape_f = PyArray_SHAPE(ap_f);
    if (shape_f[0] != me)
    {
        PYERR(slsqp_error, "Input array f must have the same shape as (me,).");
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
    // A : ma*n, b : ma
    // E : me*n, f : me
    // G : mg*n, h : mg
    // x : n
    // for LSEI buffer for multipliers, residuals, and others : [mg + 2*me + ma]
    // for subarrays A2, G2 : [(ma + mg)*(n - me)]
    // for the later call to LSI and LDP: ((mg+2)*((n - me)+1) + 2*mg + (n - me))
    npy_intp total_size = ((ma + mg + me)*(n + 2) + me + n + (ma + mg)*(n - me) +
                            (mg + 2)*((n - me) + 1) + 2*mg + (n - me));
    double* mem_ret = malloc(total_size*sizeof(double));
    if (mem_ret == NULL)
    {
        free(indices);
        PYERR(slsqp_error, "Memory allocation for buffer failed.");
    }
    double* restrict x = &mem_ret[0];
    double* restrict a = &mem_ret[n];
    double* restrict b = &mem_ret[n*ma + n];
    double* restrict e = &mem_ret[n*ma + n + ma];
    double* restrict f = &mem_ret[n*ma + n + ma + n*me];
    double* restrict g = &mem_ret[n*ma + n + ma + n*me + me];
    double* restrict h = &mem_ret[n*ma + n + ma + n*me + me + n*mg];
    double* restrict buffer = &mem_ret[n*ma + n + ma + n*me + me + n*mg + mg];

    // Copy the data from the numpy arrays
    double* data_a = (double*)PyArray_DATA(ap_a);
    double* data_b = (double*)PyArray_DATA(ap_b);
    double* data_e = (double*)PyArray_DATA(ap_e);
    double* data_f = (double*)PyArray_DATA(ap_f);
    double* data_g = (double*)PyArray_DATA(ap_g);
    double* data_h = (double*)PyArray_DATA(ap_h);
    npy_intp* restrict strides_a = PyArray_STRIDES(ap_a);
    npy_intp* restrict stride_b = PyArray_STRIDES(ap_b);
    npy_intp* restrict strides_e = PyArray_STRIDES(ap_e);
    npy_intp* restrict stride_f = PyArray_STRIDES(ap_f);
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
    row_stride = (strides_e[0]/sizeof(double));
    col_stride = (strides_e[1]/sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < me; i++) {
            e[i + j*me] = data_e[i*row_stride + j*col_stride];
        }
    }
    row_stride = (strides_g[0]/sizeof(double));
    col_stride = (strides_g[1]/sizeof(double));
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < mg; i++) {
            g[i + j*mg] = data_g[i*row_stride + j*col_stride];
        }
    }

    row_stride = (stride_b[0]/sizeof(double));
    for (int i = 0; i < ma; i++) { b[i] = data_b[i*row_stride]; }
    row_stride = (stride_f[0]/sizeof(double));
    for (int i = 0; i < me; i++) { f[i] = data_f[i*row_stride]; }
    row_stride = (stride_h[0]/sizeof(double));
    for (int i = 0; i < mg; i++) { h[i] = data_h[i*row_stride]; }

    lsei((int)ma, (int)me, (int)mg, (int)n, a, b, e, f, g, h, x, buffer, indices, &xnorm, &mode);
    free(indices);
    // Truncate it to the size of the solution vector x.
    double* mem_x = realloc(mem_ret, n*sizeof(double));
    if (mem_x == NULL) { free(mem_ret);PYERR(slsqp_error, "Memory reallocation failed."); }
    npy_intp shape_x[1] = {n};
    PyArrayObject* ap_x = (PyArrayObject*)PyArray_SimpleNewFromData(1, shape_x, NPY_FLOAT64, mem_x);

    return Py_BuildValue("Ndi", PyArray_Return(ap_x), xnorm, mode);

}


static PyObject*
lsq_wrapper(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    npy_intp m, meq, n, augmented;
    int mode = 0;
    PyArrayObject *ap_Lf=NULL, *ap_gradx=NULL, *ap_C=NULL, *ap_d=NULL, *ap_xl=NULL, *ap_xu=NULL;

    if (!PyArg_ParseTuple(args, "nnnnO!O!O!O!O!O!",
                          &m, &meq, &n, &augmented,
                          &PyArray_Type, (PyObject **)&ap_Lf,
                          &PyArray_Type, (PyObject **)&ap_gradx,
                          &PyArray_Type, (PyObject **)&ap_C,
                          &PyArray_Type, (PyObject **)&ap_d,
                          &PyArray_Type, (PyObject **)&ap_xl,
                          &PyArray_Type, (PyObject **)&ap_xu))
    {
        return NULL;
    }

    if ((PyArray_TYPE(ap_Lf) != NPY_FLOAT64) || (PyArray_TYPE(ap_gradx) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_C) != NPY_FLOAT64) || (PyArray_TYPE(ap_d) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_xl) != NPY_FLOAT64) || (PyArray_TYPE(ap_xu) != NPY_FLOAT64))
    {
        PYERR(slsqp_error, "Inputs to lsq must be of type numpy.float64.");
    }

    int ndim_Lf = PyArray_NDIM(ap_Lf);
    if (ndim_Lf != 1) { PYERR(slsqp_error, "Input array S must be 1D."); }
    npy_intp* shape_Lf = PyArray_SHAPE(ap_Lf);
    int ndim_gradx = PyArray_NDIM(ap_gradx);
    if (ndim_gradx != 1) { PYERR(slsqp_error, "Input array t must be 1D."); }
    npy_intp* shape_gradx = PyArray_SHAPE(ap_gradx);

    if (shape_Lf[0] != (n*(n+1)/2)) { PYERR(slsqp_error, "Input array S must have the same number of elements as n*(n+1)/2."); }
    if (shape_gradx[0] != n) { PYERR(slsqp_error, "Input array gradx must have the same number of elements as n."); }

    // Over allocate indices a bit for transposing problems and other uses
    int* indices = malloc((m+2*n)*sizeof(int));
    if (indices == NULL) { PYERR(slsqp_error, "Memory allocation for indices failed."); }

    // Allocate memory for the problem data and the algorithm.
    size_t lsei_size = ((n+1) + (m-meq) + meq)*(n+3) + ((n+2) + (m-meq))*(n-meq)
                       + ((m-meq) + 2)*((n+1-meq)+1) + 2*(m-meq) + (n+1-meq)
                       + meq + n;
    size_t lsq_size = (n+1)*(n+2) + (n+2)*meq + meq + (m - meq + 2*n)*(n+1) + (m - meq) + 2*(n+1) + m*n;
    size_t total_size = lsei_size + lsq_size;

    double* mem_ret = calloc(total_size, sizeof(double));
    if (mem_ret == NULL)
    {
        free(indices);
        PYERR(slsqp_error, "Memory allocation for buffer failed.");
    }

    double* restrict Lf =     &mem_ret[0];
    double* restrict gradx =  &mem_ret[n*(n+1)/2];
    double* restrict C =      &mem_ret[n*(n+1)/2 + n];
    double* restrict d =      &mem_ret[n*(n+1)/2 + n + n*m];
    double* restrict xl =     &mem_ret[n*(n+1)/2 + n + n*m + m];
    double* restrict xu =     &mem_ret[n*(n+1)/2 + n + n*m + m + n];
    double* restrict x =      &mem_ret[n*(n+1)/2 + n + n*m + m + n + n];
    double* restrict y =      &mem_ret[n*(n+1)/2 + n + n*m + m + n + n + n];
    double* restrict buffer = &mem_ret[n*(n+1)/2 + n + n*m + m + n + n + n + m + n + n];

    // Copy the data from the numpy arrays
    double* data_Lf = (double*)PyArray_DATA(ap_Lf);
    double* data_gradx = (double*)PyArray_DATA(ap_gradx);
    double* data_C = (double*)PyArray_DATA(ap_C);
    double* data_d = (double*)PyArray_DATA(ap_d);
    double* data_xl = (double*)PyArray_DATA(ap_xl);
    double* data_xu = (double*)PyArray_DATA(ap_xu);
    npy_intp* restrict strides_Lf = PyArray_STRIDES(ap_Lf);
    npy_intp* restrict strides_gradx = PyArray_STRIDES(ap_gradx);
    npy_intp* restrict strides_C = PyArray_STRIDES(ap_C);
    npy_intp* restrict strides_d = PyArray_STRIDES(ap_d);
    npy_intp* restrict strides_xl = PyArray_STRIDES(ap_xl);
    npy_intp* restrict strides_xu = PyArray_STRIDES(ap_xu);
    npy_intp row_stride = (strides_Lf[0]/sizeof(double));
    npy_intp row_stride_gradx = (strides_gradx[0]/sizeof(double));
    npy_intp row_stride_C = (strides_C[0]/sizeof(double));
    npy_intp col_stride_C = (strides_C[1]/sizeof(double));
    npy_intp row_stride_d = (strides_d[0]/sizeof(double));
    npy_intp row_stride_xl = (strides_xl[0]/sizeof(double));
    npy_intp row_stride_xu = (strides_xu[0]/sizeof(double));

    // Copy the data from the numpy arrays in Fortran order
    for (int i = 0; i < n*(n+1)/2; i++) { Lf[i] = data_Lf[i*row_stride]; }

    // Only C matrix is 2D use both strides to store in Fortran order
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            C[i + j*m] = data_C[i*row_stride_C + j*col_stride_C];
        }
    }

    for (int i = 0; i < m; i++) { d[i] = data_d[i*row_stride_d]; }

    for (int i = 0; i < n; i++) {
        gradx[i] = data_gradx[i*row_stride_gradx];
        xl[i] = data_xl[i*row_stride_xl];
        xu[i] = data_xu[i*row_stride_xu];
    }

    // Call lsq
    lsq((int)m, (int)meq, (int)n, (int)augmented, 0, Lf, gradx, C, d, xl, xu, x, y, buffer, indices, &mode);
    free(indices);
    // Carry the solution and multipliers to the top of the buffer.
    for (int i = 0; i < m + 3*n; i++) { mem_ret[i] = mem_ret[n*(n+1)/2 + n + n*m + m + n + n + i];}
    npy_intp newshape_y[1] = {m + 2*n};
    PyArrayObject* ap_y = (PyArrayObject*)PyArray_SimpleNew(1, newshape_y, NPY_FLOAT64);
    double* data_y = (double*)PyArray_DATA(ap_y);
    for (int i = n; i < m + 3*n; i++) { data_y[i-n] = mem_ret[i]; }

    // We are only interested in the solution vector x and mode
    // Truncate it to the size of the solution vector x.
    double* mem_x = realloc(mem_ret, n*sizeof(double));
    if (mem_x == NULL) { free(mem_ret);PYERR(slsqp_error, "Memory reallocation failed."); }
    npy_intp shape_x[1] = {n};
    PyArrayObject* ap_x = (PyArrayObject*)PyArray_SimpleNewFromData(1, shape_x, NPY_FLOAT64, mem_x);

    return Py_BuildValue("NNi", PyArray_Return(ap_x), PyArray_Return(ap_y), mode);
}


static PyObject*
slsqp(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    // TODO:
    // Exact calloc size (needed?)

    PyArrayObject *ap_gradx=NULL, *ap_C=NULL, *ap_d=NULL, *ap_mult=NULL;
    PyArrayObject *ap_sol =NULL, *ap_xl=NULL, *ap_xu=NULL, *ap_buffer=NULL;
    PyArrayObject* ap_indices=NULL;
    PyObject* input_dict = NULL;
    double funx;
    struct SLSQP_static_vars Vars;

    // The Python input should provide with a dictionary that maps to the struct
    // SLSQP_static_vars. Necessary fields that would make the algorithm change
    // behavior are m, meq, n, acc, maxiter, and mode. The rest can be left as zero.
    // Changing values mid run is not recommended as they hold the internal state
    // of the algorithm.
    // The reason why they are returned is to make the algorithm stateless.

    // The required arrays C, d, x, xl, xu, gradx, sol are passed as numpy arrays.
    // The remaining arrays are going to be allocated in the buffer.

    if (!PyArg_ParseTuple(args, "O!dO!O!O!O!O!O!",
                          &PyDict_Type, (PyObject **)&input_dict,
                          &funx,
                          &PyArray_Type, (PyObject **)&ap_C,
                          &PyArray_Type, (PyObject **)&ap_d,
                          &PyArray_Type, (PyObject **)&ap_gradx,
                          &PyArray_Type, (PyObject **)&ap_sol,
                          &PyArray_Type, (PyObject **)&ap_mult,
                          &PyArray_Type, (PyObject **)&ap_xl,
                          &PyArray_Type, (PyObject **)&ap_xu,
                          &PyArray_Type, (PyObject **)&ap_buffer,
                          &PyArray_Type, (PyObject **)&ap_indices))
    {
        return NULL;
    }

    // Initialize the struct that will be populated from dict with zeros
    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for doubles and ints.
    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (name##_obj == NULL) { PYERR(slsqp_error, #name " not found in the dictionary."); } \
        Vars.name = PyFloat_AsDouble(name##_obj); \
    STRUCT_DOUBLE_FIELD_NAMES
    #undef X

    #define X(name) name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (name##_obj == NULL) { PYERR(slsqp_error, #name " not found in the dictionary."); } \
        Vars.name = PyLong_AsLong(name##_obj); \
    STRUCT_INT_FIELD_NAMES
    #undef X

    if ((PyArray_TYPE(ap_C) != NPY_FLOAT64) || (PyArray_TYPE(ap_d) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_gradx) != NPY_FLOAT64) || (PyArray_TYPE(ap_sol) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_xl) != NPY_FLOAT64) || (PyArray_TYPE(ap_xu) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_buffer) != NPY_FLOAT64) || (PyArray_TYPE(ap_indices) != NPY_INT32))
    {
        PYERR(slsqp_error, "All inputs to slsqp must be of type numpy.float64, "
                           "except \"indices\" which must be of numpy.int32.");
    }

    // Buffer is 1D hence both F and C contiguous, test with either of them.
    if (!PyArray_IS_C_CONTIGUOUS(ap_buffer))
    {
        PYERR(slsqp_error, "Input array buffer must be 1d contiguous.");
    }

    // Derive the number of variables from the solution vector length.
    int ndim_sol = PyArray_NDIM(ap_sol);
    if (ndim_sol != 1) { PYERR(slsqp_error, "Input array sol must be 1D."); }
    npy_intp* shape_sol = PyArray_SHAPE(ap_sol);
    if ((int)shape_sol[0] != Vars.n) {
        PYERR(slsqp_error, "Input array \"sol\" must have at least n elements.");
    }
    int ndim_mult = PyArray_NDIM(ap_mult);
    if (ndim_mult != 1) { PYERR(slsqp_error, "Input array \"mult\" must be 1D."); }
    npy_intp* shape_mult = PyArray_SHAPE(ap_mult);
    if ((int)shape_mult[0] != 2*Vars.n + Vars.m + 2) {
        PYERR(slsqp_error, "Input array \"mult\" must have m + 2*n + 2 elements.");
    }

    // Derive the number of constraints from the row number of A
    int ndim_C = PyArray_NDIM(ap_C);
    if (ndim_C != 2) { PYERR(slsqp_error, "Input array \"C\" must be 2D."); }
    npy_intp* shape_C = PyArray_SHAPE(ap_C);
    if (Vars.m = (int)shape_C[0]) { PYERR(slsqp_error, "Input array \"C\" must have  \"m\" rows."); }
    if (Vars.n = (int)shape_C[1]) { PYERR(slsqp_error, "Input array \"C\" must have  \"n\" columns."); }

    int ndim_d = PyArray_NDIM(ap_d);
    if (ndim_d != 1) { PYERR(slsqp_error, "Input array C must be 1D."); }
    npy_intp* shape_d = PyArray_SHAPE(ap_d);
    if (Vars.m != (int)shape_d[0]) { PYERR(slsqp_error, "Input array d must have the same number of rows as C."); }

    int ndim_gradx = PyArray_NDIM(ap_gradx);
    if (ndim_gradx != 1) { PYERR(slsqp_error, "Input array gradx must be 1D."); }
    npy_intp* shape_gradx = PyArray_SHAPE(ap_gradx);
    int ndim_xl = PyArray_NDIM(ap_xl);
    if (ndim_xl != 1) { PYERR(slsqp_error, "Input array xl must be 1D."); }
    npy_intp* shape_xl = PyArray_SHAPE(ap_xl);
    int ndim_xu = PyArray_NDIM(ap_xu);
    if (ndim_xu != 1) { PYERR(slsqp_error, "Input array xu must be 1D."); }
    npy_intp* shape_xu = PyArray_SHAPE(ap_xu);

    __slsqp_body(Vars, &funx, ap_gradx, ap_C, ap_d, ap_sol, ap_mult, ap_xl, ap_xu, ap_buffer, ap_indices);

    // Map struct variables back to dictionary.
    #define X(name) \
        PyObject* name##_obj = PyFloat_FromDouble(Vars.name); \
        if (PyDict_SetItemString(input_dict, #name, name##_obj)) { PYERR(slsqp_error, "Setting " #name " failed."); }
    STRUCT_DOUBLE_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyLong_FromLong(Vars.name); \
        if (PyDict_SetItemString(input_dict, #name, name##_obj)) { PYERR(slsqp_error, "Setting " #name " failed."); }
    STRUCT_INT_FIELD_NAMES
    #undef X

    return;

};


// The commented wrappers are for debugging purposes and they will be optimized out
static char doc_nnls[] = ("Compute the nonnegative least squares solution.\n\n"
                           "    x, info = nnls(A)\n\n");

static char doc_slsqp[] = (
    "Sequential Least Squares Programming (SLSQP) optimizer.\n\n"
    "    x, info = slsqp(S: dict, funx: np.float64, "
    "gradx: NDArray, C: NDarray, d: NDArray, "
    "sol: NDArray, xl: NDArray, xu: NDArray, buffer: NDArray, indices: NDArray)"
    "\n\n");
/*
static char doc_lsq_wrapper[] = ("Convert the given compact representation of the"
                                 " linearized problem to an LSEI problem\n\n");

static char doc_lsei_wrapper[] = ("Compute the least squares solution subject to "
                                 "equality and inequality constraints.\n\n"
                                 "    x, xnorm, mode = lsei_wrapper(A, b, E, f, G, h)\n\n");

static char doc_lsi_wrapper[] = ("Compute the least squares solution subject to "
                                 "inequality constraints.\n\n"
                                 "    x, xnorm, mode = lsi_wrapper(A, b, G, h)\n\n");

static char doc_ldp_wrapper[] = ("Compute the least distance solution.\n\n"
                                 "    x = ldp_wrapper(G, h)\n\n");
*/

// Sentinel terminated method list.
static struct PyMethodDef slsqplib_module_methods[] = {
  {"nnls", nnls, METH_VARARGS, doc_nnls},
  {"slsqp", slsqp, METH_VARARGS, doc_slsqp},
  // {"ldp_wrapper", ldp_wrapper, METH_VARARGS, doc_ldp_wrapper},
  // {"lsi_wrapper", lsi_wrapper, METH_VARARGS, doc_lsi_wrapper},
  // {"lsei_wrapper", lsei_wrapper, METH_VARARGS, doc_lsei_wrapper},
  // {"lsq_wrapper", lsq_wrapper, METH_VARARGS, doc_lsq_wrapper},
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