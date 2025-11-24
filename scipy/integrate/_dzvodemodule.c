/**
 * @file _dzvodemodule.c
 * @brief Python extension module for DVODE and ZVODE C implementations
 *
 * This module provides Python bindings for the thread-safe C implementations
 * of DVODE (real ODE solver) and ZVODE (complex ODE solver). It replaces the
 * f2py-based Fortran wrappers with a direct C extension that offers:
 *   - Thread safety (no global state)
 *   - Better error handling
 *   - Efficient memory management
 *   - Support for both real and complex ODEs
 *
 * The interface closely matches the original Fortran VODE/ZVODE API to maintain
 * compatibility with existing SciPy code.
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "ccallback.h"
#include <string.h>
#include <complex.h>

#include "src/vode.h"
#include "src/zvode.h"

// Error handling macros
#define PYERR(errobj, message) { \
    PyErr_SetString(errobj, message); \
    goto fail; \
}

#define PYERR2(errobj, message) { \
    PyErr_Print(); \
    PyErr_SetString(errobj, message); \
    goto fail; \
}

// Module-level error object
static PyObject *dvode_error;

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Copy array from NumPy array to Fortran column-major layout with padding
 *
 * Copies matrix data from a NumPy array (with arbitrary memory layout) to Fortran
 * column-major format with potential padding. Uses NumPy stride information to
 * correctly handle any input layout (C-order, F-order, or other).
 *
 * @param f Destination array (Fortran column-major)
 * @param ldf Leading dimension of destination array
 * @param nrows Number of rows to copy
 * @param ncols Number of columns to copy
 * @param src_array Source NumPy array (2D)
 */
static void
copy_array_to_fortran(double* restrict f, const int ldf, const int nrows,
                     const int ncols, PyArrayObject* src_array)
{
    int i, j;
    char *src_data = PyArray_DATA(src_array);
    int ndim = PyArray_NDIM(src_array);

    if (ndim == 0) {
        // 0D array: single scalar
        f[0] = *(double*)src_data;
    } else if (ndim == 1) {
        // 1D array: row vector
        npy_intp stride0 = PyArray_STRIDE(src_array, 0);
        for (j = 0; j < ncols; ++j) {
            f[j] = *(double*)(src_data + stride0*j);
        }
    } else {
        // 2D array: use full stride information
        npy_intp stride0 = PyArray_STRIDE(src_array, 0);
        npy_intp stride1 = PyArray_STRIDE(src_array, 1);

        for (j = 0; j < ncols; ++j) {
            for (i = 0; i < nrows; ++i) {
                // f[i,j] = src[i,j] using stride information
                f[ldf*j + i] = *(double*)(src_data + stride0*i + stride1*j);
            }
        }
    }
}

/**
 * @brief Copy complex array from NumPy array to Fortran column-major layout with padding
 *
 * Copies complex matrix data from a NumPy array (with arbitrary memory layout) to Fortran
 * column-major format with potential padding. Uses NumPy stride information to
 * correctly handle any input layout (C-order, F-order, or other).
 *
 * @param f Destination array (Fortran column-major)
 * @param ldf Leading dimension of destination array
 * @param nrows Number of rows to copy
 * @param ncols Number of columns to copy
 * @param src_array Source NumPy array (2D, complex)
 */
static void
copy_complex_array_to_fortran(ZVODE_CPLX_TYPE* restrict f, const int ldf,
                              const int nrows, const int ncols,
                              PyArrayObject* src_array)
{
    int i, j;
    char *src_data = PyArray_DATA(src_array);
    int ndim = PyArray_NDIM(src_array);

    if (ndim == 0) {
        // 0D array: single scalar
        f[0] = *(ZVODE_CPLX_TYPE*)src_data;
    } else if (ndim == 1) {
        // 1D array: row vector
        npy_intp stride0 = PyArray_STRIDE(src_array, 0);
        for (j = 0; j < ncols; ++j) {
            f[j] = *(ZVODE_CPLX_TYPE*)(src_data + stride0*j);
        }
    } else {
        // 2D array: use full stride information
        npy_intp stride0 = PyArray_STRIDE(src_array, 0);
        npy_intp stride1 = PyArray_STRIDE(src_array, 1);

        for (j = 0; j < ncols; ++j) {
            for (i = 0; i < nrows; ++i) {
                // f[i,j] = src[i,j] using stride information
                f[ldf*j + i] = *(ZVODE_CPLX_TYPE*)(src_data + stride0*i + stride1*j);
            }
        }
    }
}

// ============================================================================
// Callback Infrastructure
// ============================================================================

/**
 * @brief Callback structure for DVODE/ZVODE user functions
 */
typedef struct {
    PyObject *ode_function;    // User's ODE function f(t, y) or f(y, t)
    PyObject *jac_function;    // User's Jacobian function (optional)
    PyObject *func_args;       // Extra arguments tuple
    int jac_type;              // Jacobian type: 0=full, 3=banded
} dvode_callback_t;

// Thread-local storage for callbacks
static SCIPY_TLS dvode_callback_t* current_dvode_callback = NULL;
static SCIPY_TLS dvode_callback_t* current_zvode_callback = NULL;

// ============================================================================
// DVODE Callback Thunks (Real ODEs)
// ============================================================================

/**
 * @brief DVODE function thunk - interfaces Python callback with C code
 */
static void
dvode_function_thunk(int neq, double t, double* y, double* ydot, double* rpar, int* ipar)
{
    (void)rpar;  // Unused
    (void)ipar;  // Unused

    if (!current_dvode_callback || !current_dvode_callback->ode_function) {
        PyErr_SetString(PyExc_RuntimeError, "DVODE callback not initialized");
        return;
    }

    npy_intp dims[1] = {neq};

    // Create capsule for y data (no copy)
    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) { return; }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*)y);
    if (!py_y) {
        Py_DECREF(capsule);
        return;
    }

    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        Py_DECREF(py_y);
        return;
    }

    // Build argument tuple: (t, y, *extra_args)
    PyObject *args_tuple;
    PyObject *py_t = PyFloat_FromDouble(t);
    if (!py_t) {
        Py_DECREF(py_y);
        return;
    }

    if (current_dvode_callback->func_args) {
        Py_ssize_t nargs = PyTuple_Size(current_dvode_callback->func_args);
        args_tuple = PyTuple_New(2 + nargs);
        PyTuple_SET_ITEM(args_tuple, 0, py_t);
        PyTuple_SET_ITEM(args_tuple, 1, py_y);
        for (Py_ssize_t i = 0; i < nargs; i++) {
            PyObject *arg = PyTuple_GET_ITEM(current_dvode_callback->func_args, i);
            Py_INCREF(arg);
            PyTuple_SET_ITEM(args_tuple, 2 + i, arg);
        }
    } else {
        args_tuple = PyTuple_Pack(2, py_t, py_y);
    }

    // Call Python function
    PyObject *result = PyObject_CallObject(current_dvode_callback->ode_function, args_tuple);
    Py_DECREF(args_tuple);
    if (!result) { return; }

    // Convert result to array
    PyArrayObject *result_array = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (!result_array) {
        Py_DECREF(result);
        return;
    }

    // Validate dimensions
    if (PyArray_NDIM(result_array) > 1 || PyArray_Size((PyObject*)result_array) != neq) {
        PyErr_SetString(PyExc_ValueError, "DVODE function must return array of length NEQ");
        Py_DECREF(result_array);
        Py_DECREF(result);
        return;
    }

    // Copy result to output
    double *result_data = (double*)PyArray_DATA(result_array);
    memcpy(ydot, result_data, neq * sizeof(double));

    Py_DECREF(result_array);
    Py_DECREF(result);
}

/**
 * @brief DVODE Jacobian thunk - interfaces Python Jacobian callback with C code
 */
static void
dvode_jacobian_thunk(int neq, double t, double* y, int ml, int mu,
                     double* pd, int nrowpd, double* rpar, int* ipar)
{
    (void)rpar;  // Unused
    (void)ipar;  // Unused

    if (!current_dvode_callback || !current_dvode_callback->jac_function) {
        PyErr_SetString(PyExc_RuntimeError, "DVODE Jacobian callback not initialized");
        return;
    }

    npy_intp dims[1] = {neq};

    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) { return; }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT64, (void*)y);
    if (!py_y) {
        Py_DECREF(capsule);
        return;
    }

    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        Py_DECREF(py_y);
        return;
    }

    // Build argument tuple: (t, y, *extra_args)
    PyObject *jac_args_tuple;
    PyObject *py_t = PyFloat_FromDouble(t);
    if (!py_t) {
        Py_DECREF(py_y);
        return;
    }

    if (current_dvode_callback->func_args) {
        Py_ssize_t nargs = PyTuple_Size(current_dvode_callback->func_args);
        jac_args_tuple = PyTuple_New(2 + nargs);
        PyTuple_SET_ITEM(jac_args_tuple, 0, py_t);
        PyTuple_SET_ITEM(jac_args_tuple, 1, py_y);
        for (Py_ssize_t i = 0; i < nargs; i++) {
            PyObject *arg = PyTuple_GET_ITEM(current_dvode_callback->func_args, i);
            Py_INCREF(arg);
            PyTuple_SET_ITEM(jac_args_tuple, 2 + i, arg);
        }
    } else {
        jac_args_tuple = PyTuple_Pack(2, py_t, py_y);
    }

    // Call Python Jacobian function
    PyObject *result = PyObject_CallObject(current_dvode_callback->jac_function, jac_args_tuple);
    Py_DECREF(jac_args_tuple);
    if (!result) { return; }

    PyArrayObject *result_array = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (!result_array) {
        Py_DECREF(result);
        return;
    }

    // Determine expected dimensions based on Jacobian type
    npy_intp expected_nrows, expected_ncols;
    expected_ncols = neq;
    if (current_dvode_callback->jac_type == 3) {
        // Banded Jacobian: user provides compressed format (ml + mu + 1, neq)
        expected_nrows = ml + mu + 1;
    } else {
        // Full Jacobian
        expected_nrows = neq;
    }

    // Validate result dimensions
    npy_intp ndim = PyArray_NDIM(result_array);
    if (ndim > 2) {
        PyErr_Format(PyExc_RuntimeError,
            "The Jacobian array must be at most two dimensional, but got ndim=%ld.", ndim);
        Py_DECREF(result_array);
        Py_DECREF(result);
        return;
    }

    npy_intp *dims_result = PyArray_DIMS(result_array);
    int dim_error = 0;
    if (ndim == 0) {
        if ((expected_nrows != 1) || (expected_ncols != 1)) {
            dim_error = 1;
        }
    } else if (ndim == 1) {
        if ((expected_nrows != 1) || (dims_result[0] != expected_ncols)) {
            dim_error = 1;
        }
    } else if (ndim == 2) {
        if ((dims_result[0] != expected_nrows) || (dims_result[1] != expected_ncols)) {
            dim_error = 1;
        }
    }

    if (dim_error) {
        const char *jac_type_str = (current_dvode_callback->jac_type == 3) ? "banded " : "";
        PyErr_Format(PyExc_RuntimeError,
            "Expected a %sJacobian array with shape (%ld, %ld), but got (%ld, %ld)",
            jac_type_str, expected_nrows, expected_ncols,
            (ndim >= 1) ? dims_result[0] : 1,
            (ndim >= 2) ? dims_result[1] : ((ndim == 1) ? dims_result[0] : 1));
        Py_DECREF(result_array);
        Py_DECREF(result);
        return;
    }

    // Copy result to output array with proper layout handling

    if ((current_dvode_callback->jac_type == 0) && PyArray_IS_C_CONTIGUOUS(result_array)) {
        // Full Jacobian in C-contiguous format - can use memcpy
        double *src_data = (double*)PyArray_DATA(result_array);
        memcpy(pd, src_data, neq * nrowpd * sizeof(double));
    } else {
        // Use stride-aware copy for any other layout (F-contiguous, banded, etc.)
        npy_intp m;

        if (current_dvode_callback->jac_type == 3) {
            // Banded Jacobian: user provides compressed format (ml + mu + 1 rows)
            // DVODE expects it in work array with leading dimension nrowpd
            m = ml + mu + 1;
        } else {
            // Full Jacobian
            m = neq;
        }

        copy_array_to_fortran(pd, nrowpd, m, neq, result_array);
    }

    Py_DECREF(result_array);
    Py_DECREF(result);
}

// ============================================================================
// ZVODE Callback Thunks (Complex ODEs)
// ============================================================================

/**
 * @brief ZVODE function thunk - interfaces Python callback with C code
 */
static void
zvode_function_thunk(int neq, double t, ZVODE_CPLX_TYPE* y, ZVODE_CPLX_TYPE* ydot,
                     ZVODE_CPLX_TYPE* rpar, int* ipar)
{
    (void)rpar;  // Unused
    (void)ipar;  // Unused

    if (!current_zvode_callback || !current_zvode_callback->ode_function) {
        PyErr_SetString(PyExc_RuntimeError, "ZVODE callback not initialized");
        return;
    }

    npy_intp dims[1] = {neq};

    // Create capsule for y data (no copy)
    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) { return; }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_COMPLEX128, (void*)y);
    if (!py_y) {
        Py_DECREF(capsule);
        return;
    }

    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        Py_DECREF(py_y);
        return;
    }

    // Build argument tuple: (t, y, *extra_args)
    PyObject *args_tuple;
    PyObject *py_t = PyFloat_FromDouble(t);
    if (!py_t) {
        Py_DECREF(py_y);
        return;
    }

    if (current_zvode_callback->func_args) {
        Py_ssize_t nargs = PyTuple_Size(current_zvode_callback->func_args);
        args_tuple = PyTuple_New(2 + nargs);
        PyTuple_SET_ITEM(args_tuple, 0, py_t);
        PyTuple_SET_ITEM(args_tuple, 1, py_y);
        for (Py_ssize_t i = 0; i < nargs; i++) {
            PyObject *arg = PyTuple_GET_ITEM(current_zvode_callback->func_args, i);
            Py_INCREF(arg);
            PyTuple_SET_ITEM(args_tuple, 2 + i, arg);
        }
    } else {
        args_tuple = PyTuple_Pack(2, py_t, py_y);
    }

    // Call Python function
    PyObject *result = PyObject_CallObject(current_zvode_callback->ode_function, args_tuple);
    Py_DECREF(args_tuple);
    if (!result) { return; }

    // Convert result to complex array
    PyArrayObject *result_array = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_COMPLEX128, NPY_ARRAY_IN_ARRAY);
    if (!result_array) {
        Py_DECREF(result);
        return;
    }

    // Validate dimensions
    if (PyArray_NDIM(result_array) > 1 || PyArray_Size((PyObject*)result_array) != neq) {
        PyErr_SetString(PyExc_ValueError, "ZVODE function must return complex array of length NEQ");
        Py_DECREF(result_array);
        Py_DECREF(result);
        return;
    }

    // Copy result to output
    ZVODE_CPLX_TYPE *result_data = (ZVODE_CPLX_TYPE*)PyArray_DATA(result_array);
    memcpy(ydot, result_data, neq * sizeof(ZVODE_CPLX_TYPE));

    Py_DECREF(result_array);
    Py_DECREF(result);
}

/**
 * @brief ZVODE Jacobian thunk - interfaces Python Jacobian callback with C code
 */
static void
zvode_jacobian_thunk(int neq, double t, ZVODE_CPLX_TYPE* y, int ml, int mu,
                     ZVODE_CPLX_TYPE* pd, int nrowpd, ZVODE_CPLX_TYPE* rpar, int* ipar)
{
    (void)rpar;  // Unused
    (void)ipar;  // Unused

    if (!current_zvode_callback || !current_zvode_callback->jac_function) {
        PyErr_SetString(PyExc_RuntimeError, "ZVODE Jacobian callback not initialized");
        return;
    }

    npy_intp dims[1] = {neq};

    PyObject *capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!capsule) { return; }

    PyObject* py_y = PyArray_SimpleNewFromData(1, dims, NPY_COMPLEX128, (void*)y);
    if (!py_y) {
        Py_DECREF(capsule);
        return;
    }

    if (PyArray_SetBaseObject((PyArrayObject*)py_y, capsule) < 0) {
        Py_DECREF(py_y);
        return;
    }

    // Build argument tuple: (t, y, *extra_args)
    PyObject *jac_args_tuple;
    PyObject *py_t = PyFloat_FromDouble(t);
    if (!py_t) {
        Py_DECREF(py_y);
        return;
    }

    if (current_zvode_callback->func_args) {
        Py_ssize_t nargs = PyTuple_Size(current_zvode_callback->func_args);
        jac_args_tuple = PyTuple_New(2 + nargs);
        PyTuple_SET_ITEM(jac_args_tuple, 0, py_t);
        PyTuple_SET_ITEM(jac_args_tuple, 1, py_y);
        for (Py_ssize_t i = 0; i < nargs; i++) {
            PyObject *arg = PyTuple_GET_ITEM(current_zvode_callback->func_args, i);
            Py_INCREF(arg);
            PyTuple_SET_ITEM(jac_args_tuple, 2 + i, arg);
        }
    } else {
        jac_args_tuple = PyTuple_Pack(2, py_t, py_y);
    }

    // Call Python Jacobian function
    PyObject *result = PyObject_CallObject(current_zvode_callback->jac_function, jac_args_tuple);
    Py_DECREF(jac_args_tuple);
    if (!result) { return; }

    PyArrayObject *result_array = (PyArrayObject*)PyArray_FROM_OTF(result, NPY_COMPLEX128, NPY_ARRAY_IN_ARRAY);
    if (!result_array) {
        Py_DECREF(result);
        return;
    }

    // Determine expected dimensions based on Jacobian type
    npy_intp expected_nrows, expected_ncols;
    expected_ncols = neq;
    if (current_zvode_callback->jac_type == 3) {
        // Banded Jacobian: user provides compressed format (ml + mu + 1, neq)
        expected_nrows = ml + mu + 1;
    } else {
        // Full Jacobian
        expected_nrows = neq;
    }

    // Validate result dimensions
    npy_intp ndim = PyArray_NDIM(result_array);
    if (ndim > 2) {
        PyErr_Format(PyExc_RuntimeError,
            "The Jacobian array must be at most two dimensional, but got ndim=%ld.", ndim);
        Py_DECREF(result_array);
        Py_DECREF(result);
        return;
    }

    npy_intp *dims_result = PyArray_DIMS(result_array);
    int dim_error = 0;
    if (ndim == 0) {
        if ((expected_nrows != 1) || (expected_ncols != 1)) {
            dim_error = 1;
        }
    } else if (ndim == 1) {
        if ((expected_nrows != 1) || (dims_result[0] != expected_ncols)) {
            dim_error = 1;
        }
    } else if (ndim == 2) {
        if ((dims_result[0] != expected_nrows) || (dims_result[1] != expected_ncols)) {
            dim_error = 1;
        }
    }

    if (dim_error) {
        const char *jac_type_str = (current_zvode_callback->jac_type == 3) ? "banded " : "";
        PyErr_Format(PyExc_RuntimeError,
            "Expected a %sJacobian array with shape (%ld, %ld), but got (%ld, %ld)",
            jac_type_str, expected_nrows, expected_ncols,
            (ndim >= 1) ? dims_result[0] : 1,
            (ndim >= 2) ? dims_result[1] : ((ndim == 1) ? dims_result[0] : 1));
        Py_DECREF(result_array);
        Py_DECREF(result);
        return;
    }

    // Copy result to output array with proper layout handling
    npy_intp m;

    if (current_zvode_callback->jac_type == 3) {
        // Banded Jacobian: user provides compressed format (ml + mu + 1 rows)
        // ZVODE expects it in work array with leading dimension nrowpd
        m = ml + mu + 1;
    } else {
        // Full Jacobian
        m = neq;
    }

    copy_complex_array_to_fortran(pd, nrowpd, m, neq, result_array);

    Py_DECREF(result_array);
    Py_DECREF(result);
}

// ============================================================================
// Callback Management
// ============================================================================

static int
setup_dvode_callback(dvode_callback_t *callback, PyObject *ode_func,
                     PyObject *jac_func, PyObject *extra_args,
                     int jac_type)
{
    callback->ode_function = ode_func;
    callback->jac_function = jac_func;
    callback->jac_type = jac_type;
    callback->func_args = NULL;

    if (ode_func) { Py_INCREF(ode_func); }
    if (jac_func && jac_func != Py_None) { Py_INCREF(jac_func); }

    if (extra_args && PyTuple_Check(extra_args)) {
        Py_INCREF(extra_args);
        callback->func_args = extra_args;
    } else if (!extra_args || extra_args == Py_None) {
        callback->func_args = NULL;
    } else {
        PyErr_SetString(PyExc_TypeError, "Extra arguments must be a tuple");
        return -1;
    }

    return 0;
}

static void
cleanup_dvode_callback(dvode_callback_t *callback)
{
    if (callback->func_args) {
        Py_DECREF(callback->func_args);
        callback->func_args = NULL;
    }
    if (callback->ode_function) {
        Py_DECREF(callback->ode_function);
        callback->ode_function = NULL;
    }
    if (callback->jac_function && callback->jac_function != Py_None) {
        Py_DECREF(callback->jac_function);
        callback->jac_function = NULL;
    }
}


// ============================================================================
// DVODE Python Interface
// ============================================================================

static char doc_dvode[] =
    "y, t, istate = dvode(f, jac, y0, t0, t1, rtol, atol, itask, istate, "
    "rwork, iwork, mf, [args])\n\n"
    "Solve real ODE system dy/dt = f(t,y) using variable-coefficient ODE solver.\n\n"
    "Parameters\n"
    "----------\n"
    "f : callable\n"
    "    Right-hand side function f(t, y, *args) returning dy/dt.\n"
    "jac : callable or None\n"
    "    Jacobian function jac(t, y, *args) returning df/dy (optional).\n"
    "y0 : array_like\n"
    "    Initial condition vector.\n"
    "t0 : float\n"
    "    Initial time.\n"
    "t1 : float\n"
    "    Output time.\n"
    "rtol : array_like\n"
    "    Relative tolerance(s).\n"
    "atol : array_like\n"
    "    Absolute tolerance(s).\n"
    "itask : int\n"
    "    Task indicator (1-5).\n"
    "istate : int\n"
    "    State indicator (1=first call, 2=continue, 3=continue with new params).\n"
    "rwork : array_like\n"
    "    Real workspace array (preserves solver state between calls).\n"
    "iwork : array_like\n"
    "    Integer workspace array (preserves solver state between calls).\n"
    "mf : int\n"
    "    Method flag.\n"
    "args : tuple, optional\n"
    "    Extra arguments passed to f and jac.\n\n"
    "Returns\n"
    "-------\n"
    "y : ndarray\n"
    "    Solution at t1.\n"
    "t : float\n"
    "    Actual time reached (may differ from t1 if error occurred).\n"
    "istate : int\n"
    "    State indicator (2=success, negative=error).\n";

static PyObject*
dvode_wrapper(PyObject* Py_UNUSED(dummy), PyObject* args, PyObject* kwdict)
{
    // Python inputs
    PyObject *fcn, *jac, *y0_obj, *rtol_obj, *atol_obj;
    PyObject *f_extra_args = NULL, *jac_extra_args = NULL;
    PyObject *state_doubles_obj = NULL, *state_ints_obj = NULL;
    double t0, t1;
    int itask, istate_in, mf;
    PyObject *rwork_obj, *iwork_obj;

    // Arrays
    PyArrayObject *ap_y = NULL, *ap_rtol = NULL, *ap_atol = NULL;
    PyArrayObject *ap_rwork = NULL, *ap_iwork = NULL;
    PyArrayObject *ap_state_doubles = NULL, *ap_state_ints = NULL;

    // Work variables
    int neq, itol, iopt = 1, lrw, liw;
    double *y, t, *rtol, *atol, *rwork;
    double *state_doubles = NULL;
    int *iwork, *state_ints = NULL;
    dvode_callback_t callback = {0};
    vode_common_struct_t solver_state = {0};

    static char *kwlist[] = {"f", "jac", "y0", "t0", "t1", "rtol", "atol",
                             "itask", "istate", "rwork", "iwork", "mf",
                             "f_params", "jac_params", "state_doubles", "state_ints", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwdict, "OOOddOOiiOOi|OOOO", kwlist,
                                     &fcn, &jac, &y0_obj, &t0, &t1,
                                     &rtol_obj, &atol_obj, &itask, &istate_in,
                                     &rwork_obj, &iwork_obj, &mf,
                                     &f_extra_args, &jac_extra_args,
                                     &state_doubles_obj, &state_ints_obj))
    {
        return NULL;
    }

    // Validate callable
    if (!PyCallable_Check(fcn)) {
        PYERR(PyExc_TypeError, "f must be callable");
    }
    if (jac != Py_None && !PyCallable_Check(jac)) {
        PYERR(PyExc_TypeError, "jac must be callable or None");
    }

    // y0: create contiguous writable copy for output
    ap_y = (PyArrayObject*)PyArray_FROMANY(y0_obj, NPY_DOUBLE, 1, 1,
                                           NPY_ARRAY_ENSURECOPY | NPY_ARRAY_CARRAY);
    if (!ap_y) { PYERR(PyExc_ValueError, "y0 must be convertible to 1-D double array"); }
    neq = (int)PyArray_SIZE(ap_y);
    if (neq <= 0) { PYERR(PyExc_ValueError, "y0 must have positive length"); }
    y = (double*)PyArray_DATA(ap_y);

    // Tolerances
    ap_rtol = (PyArrayObject*)PyArray_ContiguousFromObject(rtol_obj, NPY_DOUBLE, 0, 1);
    if (!ap_rtol) { PYERR(PyExc_ValueError, "rtol must be a scalar or 1-D array"); }
    ap_atol = (PyArrayObject*)PyArray_ContiguousFromObject(atol_obj, NPY_DOUBLE, 0, 1);
    if (!ap_atol) { PYERR(PyExc_ValueError, "atol must be a scalar or 1-D array"); }

    rtol = (double*)PyArray_DATA(ap_rtol);
    atol = (double*)PyArray_DATA(ap_atol);

    // Determine itol
    npy_intp rtol_size = PyArray_SIZE(ap_rtol);
    npy_intp atol_size = PyArray_SIZE(ap_atol);
    if (rtol_size <= 1 && atol_size <= 1) {
        itol = 1;
    } else if (rtol_size <= 1) {
        itol = 2;
    } else if (atol_size <= 1) {
        itol = 3;
    } else {
        itol = 4;
    }

    // Work arrays: ensure contiguous and correct dtype
    // Python maintains these arrays between calls for continuation (istate=2)
    ap_rwork = (PyArrayObject*)PyArray_ContiguousFromObject(rwork_obj, NPY_DOUBLE, 1, 1);
    if (!ap_rwork) { PYERR(PyExc_ValueError, "rwork must be a 1-D double array"); }
    rwork = (double*)PyArray_DATA(ap_rwork);
    lrw = (int)PyArray_SIZE(ap_rwork);

    ap_iwork = (PyArrayObject*)PyArray_ContiguousFromObject(iwork_obj, NPY_INT, 1, 1);
    if (!ap_iwork) { PYERR(PyExc_ValueError, "iwork must be a 1-D int array"); }
    iwork = (int*)PyArray_DATA(ap_iwork);
    liw = (int)PyArray_SIZE(ap_iwork);

    // Setup callback for DVODE
    // mf = 10*meth + miter, where miter=4,5 means banded Jacobian
    // So mf=14,15 (Adams banded) or mf=24,25 (BDF banded)
    int miter = mf % 10;
    int jac_type = (miter == 4 || miter == 5) ? 3 : 0;  // 3=banded, 0=full
    // For now, use f_extra_args for both f and jac callbacks (matching LSODA behavior)
    // The jac_extra_args parameter is accepted but currently unused
    if (setup_dvode_callback(&callback, fcn, jac, f_extra_args, jac_type) < 0) {
        goto fail;
    }

    // State persistence arrays: validate and get data pointers (no copy)
    if (state_doubles_obj == NULL || state_doubles_obj == Py_None ||
        state_ints_obj == NULL || state_ints_obj == Py_None) {
        PYERR(PyExc_ValueError, "state_doubles and state_ints are required");
    }

    if (!PyArray_Check(state_doubles_obj) || !PyArray_Check(state_ints_obj)) {
        PYERR(PyExc_TypeError, "state_doubles and state_ints must be NumPy arrays");
    }

    ap_state_doubles = (PyArrayObject*)state_doubles_obj;
    ap_state_ints = (PyArrayObject*)state_ints_obj;

    if (PyArray_TYPE(ap_state_doubles) != NPY_DOUBLE || !PyArray_ISCARRAY(ap_state_doubles) ||
        PyArray_NDIM(ap_state_doubles) != 1 || PyArray_DIM(ap_state_doubles, 0) != VODE_STATE_DOUBLE_SIZE) {
        PyErr_Format(PyExc_ValueError,
            "state_doubles must be a C-contiguous 1-D double array of size %d", VODE_STATE_DOUBLE_SIZE);
        goto fail;
    }

    if (PyArray_TYPE(ap_state_ints) != NPY_INT32 || !PyArray_ISCARRAY(ap_state_ints) ||
        PyArray_NDIM(ap_state_ints) != 1 || PyArray_DIM(ap_state_ints, 0) != VODE_STATE_INT_SIZE) {
        PyErr_Format(PyExc_ValueError,
            "state_ints must be a C-contiguous 1-D int32 array of size %d", VODE_STATE_INT_SIZE);
        goto fail;
    }

    state_doubles = (double*)PyArray_DATA(ap_state_doubles);
    state_ints = (int*)PyArray_DATA(ap_state_ints);

    // Set istate and t
    int istate = istate_in;
    t = t0;

    // Restore solver state from arrays for continuation calls (istate=2)
    if (istate == 2) {
        unpack_vode_state(state_doubles, state_ints, &solver_state);
    }

    // Activate callback and call DVODE
    current_dvode_callback = &callback;

    dvode(&solver_state, dvode_function_thunk, neq, y, &t, &t1, itol, rtol, atol, &itask, &istate, &iopt,
        rwork, lrw, iwork, liw, (jac != Py_None) ? dvode_jacobian_thunk : NULL, mf, NULL, NULL);

    current_dvode_callback = NULL;

    // Check for Python errors during callback
    if (PyErr_Occurred()) {
        cleanup_dvode_callback(&callback);
        goto fail;
    }

    // Save solver state to arrays for next call
    pack_vode_state(&solver_state, state_doubles, state_ints);

    // Cleanup callback
    cleanup_dvode_callback(&callback);

    // Build return tuple: (y, t, istate)
    PyObject *py_t = PyFloat_FromDouble(t);
    PyObject *py_istate = PyLong_FromLong(istate);
    PyObject *result = PyTuple_Pack(3, (PyObject*)ap_y, py_t, py_istate);

    Py_DECREF(py_t);
    Py_DECREF(py_istate);
    Py_DECREF(ap_rtol);
    Py_DECREF(ap_atol);
    Py_DECREF(ap_rwork);
    Py_DECREF(ap_iwork);
    Py_DECREF(ap_y);

    return result;

fail:
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_rtol);
    Py_XDECREF(ap_atol);
    Py_XDECREF(ap_rwork);
    Py_XDECREF(ap_iwork);
    return NULL;
}

// ============================================================================
// ZVODE Python Interface
// ============================================================================

static char doc_zvode[] =
    "y, t, istate = zvode(f, jac, y0, t0, t1, rtol, atol, itask, istate, "
    "zwork, rwork, iwork, mf, [args])\n\n"
    "Solve complex ODE system dy/dt = f(t,y) using variable-coefficient ODE solver.\n\n"
    "Parameters\n"
    "----------\n"
    "f : callable\n"
    "    Right-hand side function f(t, y, *args) returning dy/dt (complex).\n"
    "jac : callable or None\n"
    "    Jacobian function jac(t, y, *args) returning df/dy (optional).\n"
    "y0 : array_like, complex\n"
    "    Initial condition vector.\n"
    "t0 : float\n"
    "    Initial time (real).\n"
    "t1 : float\n"
    "    Output time (real).\n"
    "rtol : array_like, real\n"
    "    Relative tolerance(s).\n"
    "atol : array_like, real\n"
    "    Absolute tolerance(s).\n"
    "itask : int\n"
    "    Task indicator (1-5).\n"
    "istate : int\n"
    "    State indicator (1=first call, 2=continue, 3=continue with new params).\n"
    "zwork : array_like, complex\n"
    "    Complex workspace array.\n"
    "rwork : array_like, real\n"
    "    Real workspace array.\n"
    "iwork : array_like, int\n"
    "    Integer workspace array.\n"
    "mf : int\n"
    "    Method flag.\n"
    "args : tuple, optional\n"
    "    Extra arguments passed to f and jac.\n\n"
    "Returns\n"
    "-------\n"
    "y : ndarray, complex\n"
    "    Solution at t1.\n"
    "t : float\n"
    "    Actual time reached (may differ from t1 if error occurred).\n"
    "istate : int\n"
    "    State indicator (2=success, negative=error).\n";

static PyObject*
zvode_wrapper(PyObject* Py_UNUSED(dummy), PyObject* args, PyObject* kwdict)
{

    // Python inputs
    PyObject *fcn, *jac, *y0_obj, *rtol_obj, *atol_obj;
    PyObject *f_extra_args = NULL, *jac_extra_args = NULL;
    PyObject *state_doubles_obj = NULL, *state_ints_obj = NULL;
    double t0, t1;
    int itask, istate_in, mf;
    PyObject *zwork_obj, *rwork_obj, *iwork_obj;

    // Arrays
    PyArrayObject *ap_y = NULL, *ap_rtol = NULL, *ap_atol = NULL;
    PyArrayObject *ap_zwork = NULL, *ap_rwork = NULL, *ap_iwork = NULL;
    PyArrayObject *ap_state_doubles = NULL, *ap_state_ints = NULL;

    // Work variables
    int neq, itol, iopt = 1, lzw, lrw, liw;
    ZVODE_CPLX_TYPE *y, *zwork;
    double t, *rtol, *atol, *rwork;
    double *state_doubles = NULL;
    int *iwork, *state_ints = NULL;
    dvode_callback_t callback = {0};
    zvode_common_struct_t solver_state = {0};

    static char *kwlist[] = {"f", "jac", "y0", "t0", "t1", "rtol", "atol",
                             "itask", "istate", "zwork", "rwork", "iwork", "mf",
                             "f_params", "jac_params", "state_doubles", "state_ints", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwdict, "OOOddOOiiOOOi|OOOO", kwlist,
                                     &fcn, &jac, &y0_obj, &t0, &t1,
                                     &rtol_obj, &atol_obj, &itask, &istate_in,
                                     &zwork_obj, &rwork_obj, &iwork_obj, &mf,
                                     &f_extra_args, &jac_extra_args,
                                     &state_doubles_obj, &state_ints_obj))
    {
        return NULL;
    }

    // Validate callable
    if (!PyCallable_Check(fcn)) {
        PYERR(PyExc_TypeError, "f must be callable");
    }
    if (jac != Py_None && !PyCallable_Check(jac)) {
        PYERR(PyExc_TypeError, "jac must be callable or None");
    }

    // y0: create contiguous writable copy for output (complex)
    ap_y = (PyArrayObject*)PyArray_FROMANY(y0_obj, NPY_COMPLEX128, 1, 1,
                                           NPY_ARRAY_ENSURECOPY | NPY_ARRAY_CARRAY);
    if (!ap_y) { PYERR(PyExc_ValueError, "y0 must be convertible to 1-D complex array"); }
    neq = (int)PyArray_SIZE(ap_y);
    if (neq <= 0) { PYERR(PyExc_ValueError, "y0 must have positive length"); }
    y = (ZVODE_CPLX_TYPE*)PyArray_DATA(ap_y);

    // Convert tolerances (real arrays, read-only)
    ap_rtol = (PyArrayObject*)PyArray_ContiguousFromObject(rtol_obj, NPY_DOUBLE, 0, 1);
    if (!ap_rtol) { PYERR(PyExc_ValueError, "rtol must be a scalar or 1-D array"); }

    ap_atol = (PyArrayObject*)PyArray_ContiguousFromObject(atol_obj, NPY_DOUBLE, 0, 1);
    if (!ap_atol) { PYERR(PyExc_ValueError, "atol must be a scalar or 1-D array"); }

    rtol = (double*)PyArray_DATA(ap_rtol);
    atol = (double*)PyArray_DATA(ap_atol);

    // Determine itol
    npy_intp rtol_size = PyArray_SIZE(ap_rtol);
    npy_intp atol_size = PyArray_SIZE(ap_atol);
    if (rtol_size <= 1 && atol_size <= 1) {
        itol = 1;
    } else if (rtol_size <= 1) {
        itol = 2;
    } else if (atol_size <= 1) {
        itol = 3;
    } else {
        itol = 4;
    }

    // Work arrays: ensure contiguous and correct dtype
    // Python maintains these arrays between calls for continuation (istate=2)
    ap_zwork = (PyArrayObject*)PyArray_ContiguousFromObject(zwork_obj, NPY_COMPLEX128, 1, 1);
    if (!ap_zwork) { PYERR(PyExc_ValueError, "zwork must be a 1-D complex array"); }
    zwork = (ZVODE_CPLX_TYPE*)PyArray_DATA(ap_zwork);
    lzw = (int)PyArray_SIZE(ap_zwork);

    ap_rwork = (PyArrayObject*)PyArray_ContiguousFromObject(rwork_obj, NPY_DOUBLE, 1, 1);
    if (!ap_rwork) { PYERR(PyExc_ValueError, "rwork must be a 1-D double array"); }
    rwork = (double*)PyArray_DATA(ap_rwork);
    lrw = (int)PyArray_SIZE(ap_rwork);

    ap_iwork = (PyArrayObject*)PyArray_ContiguousFromObject(iwork_obj, NPY_INT, 1, 1);
    if (!ap_iwork) { PYERR(PyExc_ValueError, "iwork must be a 1-D int array"); }
    iwork = (int*)PyArray_DATA(ap_iwork);
    liw = (int)PyArray_SIZE(ap_iwork);

    // Setup callback for ZVODE
    // mf = 10*meth + miter, where miter=4,5 means banded Jacobian
    // So mf=14,15 (Adams banded) or mf=24,25 (BDF banded)
    int miter = mf % 10;
    int jac_type = (miter == 4 || miter == 5) ? 3 : 0;  // 3=banded, 0=full

    // For now, use f_extra_args for both f and jac callbacks (matching LSODA behavior)
    // The jac_extra_args parameter is accepted but currently unused
    if (setup_dvode_callback(&callback, fcn, jac, f_extra_args, jac_type) < 0) {
        goto fail;
    }

    // State persistence arrays: validate and get data pointers (no copy)
    if (state_doubles_obj == NULL || state_doubles_obj == Py_None ||
        state_ints_obj == NULL || state_ints_obj == Py_None) {
        PYERR(PyExc_ValueError, "state_doubles and state_ints are required");
    }

    if (!PyArray_Check(state_doubles_obj) || !PyArray_Check(state_ints_obj)) {
        PYERR(PyExc_TypeError, "state_doubles and state_ints must be NumPy arrays");
    }

    ap_state_doubles = (PyArrayObject*)state_doubles_obj;
    ap_state_ints = (PyArrayObject*)state_ints_obj;

    if (PyArray_TYPE(ap_state_doubles) != NPY_DOUBLE || !PyArray_ISCARRAY(ap_state_doubles) ||
        PyArray_NDIM(ap_state_doubles) != 1 || PyArray_DIM(ap_state_doubles, 0) != ZVODE_STATE_DOUBLE_SIZE) {
        PyErr_Format(PyExc_ValueError,
            "state_doubles must be a C-contiguous 1-D double array of size %d", ZVODE_STATE_DOUBLE_SIZE);
        goto fail;
    }

    if (PyArray_TYPE(ap_state_ints) != NPY_INT32 || !PyArray_ISCARRAY(ap_state_ints) ||
        PyArray_NDIM(ap_state_ints) != 1 || PyArray_DIM(ap_state_ints, 0) != ZVODE_STATE_INT_SIZE) {
        PyErr_Format(PyExc_ValueError,
            "state_ints must be a C-contiguous 1-D int32 array of size %d", ZVODE_STATE_INT_SIZE);
        goto fail;
    }

    state_doubles = (double*)PyArray_DATA(ap_state_doubles);
    state_ints = (int*)PyArray_DATA(ap_state_ints);

    // Set istate and t
    int istate = istate_in;
    t = t0;

    // Restore solver state from arrays for continuation calls (istate=2)
    if (istate == 2) {
        unpack_zvode_state(state_doubles, state_ints, &solver_state);
    }

    // Activate callback and call ZVODE
    current_zvode_callback = &callback;

    zvode(&solver_state, zvode_function_thunk, neq, y, &t, &t1, itol, rtol, atol, &itask, &istate, &iopt,
          zwork, lzw, rwork, lrw, iwork, liw, (jac != Py_None) ? zvode_jacobian_thunk : NULL, mf, NULL, NULL);

    current_zvode_callback = NULL;

    // Check for Python errors during callback
    if (PyErr_Occurred()) {
        cleanup_dvode_callback(&callback);
        goto fail;
    }

    // Save solver state to arrays for next call
    pack_zvode_state(&solver_state, state_doubles, state_ints);

    // Cleanup callback
    cleanup_dvode_callback(&callback);

    // Build return tuple: (y, t, istate)
    PyObject *py_t = PyFloat_FromDouble(t);
    PyObject *py_istate = PyLong_FromLong(istate);
    PyObject *result = PyTuple_Pack(3, (PyObject*)ap_y, py_t, py_istate);

    Py_DECREF(py_t);
    Py_DECREF(py_istate);
    Py_DECREF(ap_rtol);
    Py_DECREF(ap_atol);
    Py_DECREF(ap_zwork);
    Py_DECREF(ap_rwork);
    Py_DECREF(ap_iwork);
    Py_DECREF(ap_y);

    return result;

fail:
    Py_XDECREF(ap_y);
    Py_XDECREF(ap_rtol);
    Py_XDECREF(ap_atol);
    Py_XDECREF(ap_zwork);
    Py_XDECREF(ap_rwork);
    Py_XDECREF(ap_iwork);
    return NULL;
}

static struct PyMethodDef dzvode_module_methods[] = {
    {"dvode", (PyCFunction)(void(*)(void))dvode_wrapper, METH_VARARGS | METH_KEYWORDS, doc_dvode},
    {"zvode", (PyCFunction)(void(*)(void))zvode_wrapper, METH_VARARGS | METH_KEYWORDS, doc_zvode},
    {NULL, NULL, 0, NULL}
};

static int
dzvode_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }

    dvode_error = PyErr_NewException("scipy.integrate._vode.error", NULL, NULL);
    if (dvode_error == NULL) {
        return -1;
    }
    Py_INCREF(dvode_error);
    if (PyModule_AddObject(module, "error", dvode_error) < 0) {
        Py_DECREF(dvode_error);
        return -1;
    }

    return 0;
}

static struct PyModuleDef_Slot dzvode_module_slots[] = {
    {Py_mod_exec, dzvode_module_exec},
#if PY_VERSION_HEX >= 0x030c00f0  // Python 3.12+
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};

static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_vode",
    .m_doc = "C-based DVODE and ZVODE ODE solvers from ODEPACK collection",
    .m_size = 0,
    .m_methods = dzvode_module_methods,
    .m_slots = dzvode_module_slots
};

PyMODINIT_FUNC
PyInit__vode(void)
{
    return PyModuleDef_Init(&moduledef);
}
