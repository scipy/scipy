#ifndef __NNLS_H
#define __NNLS_H

#include "__slsqp.h"

double dnrm2_(int* n, double* x, int* incx);
void dlarf_(char* side, int* m, int* n, double* v, int* incv, double* tau, double* c, int* ldc, double* work);
void dlarfgp_(int* n, double* alpha, double* x, int* incx, double* tau);
void dlartgp_(double* f, double* g, double* cs, double* sn, double* r);


void
nonnegative_lsq_imp(const npy_intp m, const npy_intp n,
                    double* restrict a, double* restrict b,
                    double* restrict x, double* restrict w,
                    double* restrict zz, double* restrict work,
                    int* restrict indices, const int maxiter,
                    double* rnorm, int* info);


static PyObject*
nnls(PyObject *dummy, PyObject *args) {

    int maxiter, info = 0;
    PyArrayObject *ap_A=NULL, *ap_b=NULL;
    double *buffer;
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
    nonnegative_lsq_imp(m, n, a, b, x, w, zz, work, indices, maxiter, &rnorm, &info);
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

#endif