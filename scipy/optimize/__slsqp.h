/*
 * This file and the accompanying __slsqp.c file are the C translations of the
 * Fortran77 code of the SLSQP algorithm for the SciPy project and hence inherits
 * SciPy license. The original Fortran code is available at
 * http://www.netlib.org/toms/733 written by Dieter Kraft, see:
 *
 *  ALGORITHM 733, COLLECTED ALGORITHMS FROM ACM.
 *  TRANSACTIONS ON MATHEMATICAL SOFTWARE,
 *  VOL. 20, NO. 3, SEPTEMBER, 1994, PP. 262-281.
 *  https://doi.org/10.1145/192115.192124
 *
 *
 * The original Fortran code is released for use under BSD license, with the
 * following statement from the original license holder ACM publications:
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

// BLAS/LAPACK function prototypes used in SLSQP
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
void dscal_(int* n, double* da, double* dx, int* incx);
void dtpmv_(char* uplo, char* trans, char* diag, int* n, double* ap, double* x, int* incx);
void dtpsv_(char* uplo, char* trans, char* diag, int* n, double* ap, double* x, int* incx);
void dtrsm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* a, int* lda, double* b, int* ldb);
void dtrsv_(char* uplo, char* trans, char* diag, int* n, double* a, int* lda, double* x, int* incx);


// The SLSQP_vars struct holds the state of the algorithm and passed to Python
// and back such that it is thread-safe.
struct SLSQP_vars {
    double acc, alpha, f0, gs, h1, h2, h3, h4, t, t0, tol;
    int exact, inconsistent, reset, iter, itermax, line, m, meq, mode, n;
};


void __slsqp_body(struct SLSQP_vars* S, double* funx, double* gradx, double* C, double* d, double* sol, double* mult, double* xl, double* xu, double* buffer, int* indices);


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
slsqp(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyArrayObject *ap_gradx=NULL, *ap_C=NULL, *ap_d=NULL, *ap_mult=NULL;
    PyArrayObject *ap_sol =NULL, *ap_xl=NULL, *ap_xu=NULL, *ap_buffer=NULL;
    PyArrayObject* ap_indices=NULL;
    PyObject* input_dict = NULL;
    double funx;
    struct SLSQP_vars Vars;

    // The Python input should provide with a dictionary that maps to the struct
    // SLSQP_vars. Necessary fields that would make the algorithm change
    // behavior are m, meq, n, acc, maxiter, and mode. The rest can be left as zero.
    // Changing values mid run is not recommended as they hold the internal state
    // of the algorithm.
    // The reason why they are returned is to make the algorithm stateless.

    // The required arrays C, d, x, xl, xu, gradx, sol are passed as numpy arrays.
    // The remaining arrays are going to be allocated in the buffer.

    if (!PyArg_ParseTuple(args, "O!dO!O!O!O!O!O!O!O!O!",
                          &PyDict_Type, (PyObject **)&input_dict,
                          &funx,
                          &PyArray_Type, (PyObject **)&ap_gradx,
                          &PyArray_Type, (PyObject **)&ap_C,
                          &PyArray_Type, (PyObject **)&ap_d,
                          &PyArray_Type, (PyObject **)&ap_sol,
                          &PyArray_Type, (PyObject **)&ap_mult,
                          &PyArray_Type, (PyObject **)&ap_xl,
                          &PyArray_Type, (PyObject **)&ap_xu,
                          &PyArray_Type, (PyObject **)&ap_buffer,
                          &PyArray_Type, (PyObject **)&ap_indices))
    {
        return NULL;
    }

    // Some helper x macros to pack and unpack the SLSQP_vars struct and
    // the Python dictionary.

    #define STRUCT_DOUBLE_FIELD_NAMES X(acc) X(alpha) X(f0) X(gs) X(h1) X(h2) X(h3) X(h4) X(t) X(t0) X(tol)
    #define STRUCT_INT_FIELD_NAMES X(exact) X(inconsistent) X(reset) X(iter) X(itermax) X(line) X(m) X(meq) X(mode) X(n)
    #define STRUCT_FIELD_NAMES STRUCT_INT_FIELD_NAMES STRUCT_DOUBLE_FIELD_NAMES

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for doubles and ints.
    // Initialize the struct that will be populated from dict with zeros
    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    // PyDict_GetItemString returns a borrowed reference.
    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(slsqp_error, #name " not found in the dictionary."); } \
        Vars.name = PyFloat_AsDouble(name##_obj);
    STRUCT_DOUBLE_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(slsqp_error, #name " not found in the dictionary."); } \
        Vars.name = (int)PyLong_AsLong(name##_obj);
        STRUCT_INT_FIELD_NAMES
    #undef X

    // Basic error checks for the numpy arrays.
    if ((PyArray_TYPE(ap_C) != NPY_FLOAT64) || (PyArray_TYPE(ap_d) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_gradx) != NPY_FLOAT64) || (PyArray_TYPE(ap_sol) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_xl) != NPY_FLOAT64) || (PyArray_TYPE(ap_xu) != NPY_FLOAT64) ||
        (PyArray_TYPE(ap_buffer) != NPY_FLOAT64) || (PyArray_TYPE(ap_indices) != NPY_INT32))
    {
        PYERR(slsqp_error, "All inputs to slsqp must be of type numpy.float64, "
                           "except \"indices\" which must be of numpy.int32.");
    }

    // Buffer is 1D hence both F and C contiguous, test with either of them.
    if (!PyArray_IS_C_CONTIGUOUS(ap_buffer)) { PYERR(slsqp_error, "Input array buffer must be 1d contiguous."); }

    // Derive the number of variables from the solution vector length.
    int ndim_sol = PyArray_NDIM(ap_sol);
    npy_intp* shape_sol = PyArray_SHAPE(ap_sol);
    int ndim_mult = PyArray_NDIM(ap_mult);
    npy_intp* shape_mult = PyArray_SHAPE(ap_mult);
    int ndim_C = PyArray_NDIM(ap_C);
    int ndim_d = PyArray_NDIM(ap_d);
    int ndim_gradx = PyArray_NDIM(ap_gradx);
    int ndim_xl = PyArray_NDIM(ap_xl);
    int ndim_xu = PyArray_NDIM(ap_xu);

    if (ndim_sol != 1) { PYERR(slsqp_error, "Input array sol must be 1D."); }
    if ((int)shape_sol[0] != Vars.n) { PYERR(slsqp_error, "Input array \"sol\" must have at least n elements."); }
    if (ndim_mult != 1) { PYERR(slsqp_error, "Input array \"mult\" must be 1D."); }
    if ((int)shape_mult[0] != 2*Vars.n + Vars.m + 2) { PYERR(slsqp_error, "Input array \"mult\" must have m + 2*n + 2 elements."); }
    if (ndim_C != 2) { PYERR(slsqp_error, "Input array \"C\" must be 2D."); }
    if (ndim_d != 1) { PYERR(slsqp_error, "Input array d must be 1D."); }
    if (ndim_gradx != 1) { PYERR(slsqp_error, "Input array gradx must be 1D."); }
    if (ndim_xl != 1) { PYERR(slsqp_error, "Input array xl must be 1D."); }
    if (ndim_xu != 1) { PYERR(slsqp_error, "Input array xu must be 1D."); }

    double* gradx_data = (double*)PyArray_DATA(ap_gradx);
    double* C_data = (double*)PyArray_DATA(ap_C);
    double* d_data = (double*)PyArray_DATA(ap_d);
    double* restrict sol_data = (double*)PyArray_DATA(ap_sol);
    double* mult_data = (double*)PyArray_DATA(ap_mult);
    double* restrict xl_data = (double*)PyArray_DATA(ap_xl);
    double* restrict xu_data = (double*)PyArray_DATA(ap_xu);
    double* buffer_data = (double*)PyArray_DATA(ap_buffer);
    int* indices_data = (int*)PyArray_DATA(ap_indices);

    __slsqp_body(&Vars, &funx, gradx_data, C_data, d_data, sol_data, mult_data, xl_data, xu_data, buffer_data, indices_data);

    // During the intermediate steps, there can be a few ULPs of bound violations,
    // hence we clamp the solution if given, to the finite bound values when mode
    // is 1 or -1.
    if ((Vars.mode == 1) || (Vars.mode == -1))
    {
        int n = Vars.n;
        for (int i = 0; i < n; i++)
        {
            if ((!isnan(xl_data[i])) && (sol_data[i] < xl_data[i])) { sol_data[i] = xl_data[i]; }
            else if ((!isnan(xu_data[i])) && (sol_data[i] > xu_data[i])) { sol_data[i] = xu_data[i]; }
        }
    }

    // Map struct variables back to dictionary.
    // Py_XXX_FromXXX returns a new reference, hence needs to be decremented.

    #define X(name) do { \
            PyObject* tmp_##name = PyFloat_FromDouble(Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
            Py_XDECREF(tmp_##name); \
            PYERR(slsqp_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_DOUBLE_FIELD_NAMES
    #undef X

    #define X(name) do { \
            PyObject* tmp_##name = PyLong_FromLong((long)Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
                Py_XDECREF(tmp_##name); \
                PYERR(slsqp_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_INT_FIELD_NAMES
    #undef X
    #undef STRUCT_FIELD_NAMES
    #undef STRUCT_INT_FIELD_NAMES
    #undef STRUCT_DOUBLE_FIELD_NAMES

    Py_RETURN_NONE;

};


static char doc_nnls[] = ("Compute the nonnegative least squares solution.\n\n"
                           "    x, info = nnls(A)\n\n");


static char doc_slsqp[] = (
    "Sequential Least Squares Programming (SLSQP) optimizer.\n\n"
    "    x, info = slsqp(S: dict, funx: np.float64, "
    "gradx: NDArray, C: NDarray, d: NDArray, "
    "sol: NDArray, xl: NDArray, xu: NDArray, buffer: NDArray, indices: NDArray)"
    "\n\n");


// Sentinel terminated method list.
static struct PyMethodDef slsqplib_module_methods[] = {
  {"nnls", nnls, METH_VARARGS, doc_nnls},
  {"slsqp", slsqp, METH_VARARGS, doc_slsqp},
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
