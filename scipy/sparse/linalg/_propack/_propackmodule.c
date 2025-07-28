/*
 * Python C extension module for PROPACK SVD routines
 *
 * This module provides a thread-safe interface to the PROPACK C library for computing
 * partial SVDs of large matrices.
 *
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "PROPACK/include/propack/propack.h"
#include <stdint.h>

/* Module-level error object */
static PyObject* PropackError;

/*
 * Callback context structure for user-supplied (r)matvec functions.
 *
 * This structure holds the information needed for Python callbacks.
 * Following the pattern is similar to SciPy callback.h strategy, we store this
 * in a thread-local variable that can be accessed by the callback function.
 */
typedef struct {
    PyObject* py_callback;          // Python callback function (_AProd instance)
    PyObject* py_dparm;             // Python dparm array (for compatibility)
    PyObject* py_iparm;             // Python iparm array (for compatibility)
    int m, n;                       // Matrix dimensions for validation
    int error_occurred;             // Flag to signal Python errors
    PyObject* error_type;           // Stored exception type if error occurred
    PyObject* error_value;          // Stored exception value if error occurred
    PyObject* error_traceback;      // Stored exception traceback if error occurred
} propack_callback_info_t;


// Thread-local storage for callback context
static propack_callback_info_t* _active_callback_info = NULL;

/* ===== CALLBACK INFRA ===== */

// Error handling macros
#define HANDLE_ERROR_AND_CLEANUP(info, x_view, y_view, gstate, msg) do { \
    (info)->error_occurred = 1; \
    PyErr_SetString(PyExc_RuntimeError, (msg)); \
    PyErr_Fetch(&(info)->error_type, &(info)->error_value, &(info)->error_traceback); \
    Py_XDECREF((PyObject*)(x_view)); \
    Py_XDECREF((PyObject*)(y_view)); \
    PyGILState_Release(gstate); \
} while(0)

#define STORE_PYTHON_ERROR(info) do { \
    (info)->error_occurred = 1; \
    PyErr_Fetch(&(info)->error_type, &(info)->error_value, &(info)->error_traceback); \
} while(0)

#define CLEANUP_AND_RETURN(x_view, y_view, gstate) do { \
    Py_XDECREF((PyObject*)(x_view)); \
    Py_XDECREF((PyObject*)(y_view)); \
    PyGILState_Release(gstate); \
} while(0)

/*
 * Helper function to create a zero-copy NumPy array view
 *
 * Creating NumPy array views with non-owning data is discouraged by NumPy C API
 * guidelines. Instead the array is created by setting up the base object using
 * a PyCapsule. See PyArray_SetBaseObject documentation for details.
 *
 * Parameters:
 * - data: pointer to external memory (a pointer to the PROPACK workspace)
 * - size: number of elements in the array
 * - typenum: NumPy data type (NPY_FLOAT64, NPY_FLOAT32, etc. and not C type)
 *
 * Returns:
 * - PyArrayObject* on success, NULL on failure (also setting the exception)
 */
static PyArrayObject*
create_propack_array_view(const void* data, npy_intp size, int typenum)
{
    npy_intp dims[1] = {size};

    // Create NumPy array view on external memory
    PyArrayObject* array = (PyArrayObject*)PyArray_SimpleNewFromData(1, dims, typenum, (void*)data);
    // PyArray_SimpleNewFromData already set exception
    if (array == NULL) { return NULL; }

    // Create PyCapsule to tell NumPy that someone else owns the memory
    PyObject* capsule = PyCapsule_New((void*)data, NULL, NULL);
    if (capsule == NULL) { Py_DECREF((PyObject*)array); return NULL; }

    // PyArray_SetBaseObject steals a reference to the capsule
    if (PyArray_SetBaseObject(array, capsule) < 0) {
        /* capsule reference was stolen even on failure, so don't DECREF it */
        Py_DECREF((PyObject*)array);
        return NULL;
    }
    return array;
}

/*
 * Macro to generate PROPACK callback functions for different data types
 *
 * This macro generates the complete callback function with:
 * - Proper GIL handling
 * - Zero-copy array view creation using the helper function
 * - Consistent error handling and cleanup
 * - Python callback invocation
 *
 * Parameters:
 * - suffix: function name suffix (d, s, c, z)
 * - ctype: C data type (double, float, npy_cfloat, npy_cdouble)
 * - numpy_type: NumPy type constant (NPY_FLOAT64, NPY_FLOAT32, etc.)
 */
#define DEFINE_PROPACK_CALLBACK(suffix, ctype, numpy_type) \
static void \
propack_callback_##suffix(int transa, int m, int n, const ctype* x, ctype* y, ctype* dparm, int* iparm) \
{ \
    /* Retrieve callback context from module-level variable */ \
    propack_callback_info_t* info = _active_callback_info; \
    \
    /* This should never happen in normal operation */ \
    if (info == NULL) { return; } \
    \
    /* Early return if previous error occurred */ \
    if (info->error_occurred) { return; } \
    \
    /* Acquire Python GIL for callback execution */ \
    PyGILState_STATE gstate = PyGILState_Ensure(); \
    \
    /* Create zero-copy array views using helper function */ \
    PyArrayObject* x_view = create_propack_array_view(x, n, numpy_type); \
    PyArrayObject* y_view = create_propack_array_view(y, m, numpy_type); \
    \
    if (x_view == NULL || y_view == NULL) { \
        HANDLE_ERROR_AND_CLEANUP(info, x_view, y_view, gstate, "Failed to create NumPy array views"); \
        return; \
    } \
    \
    /* \
     * Call Python callback: callback(transa, m, n, x_view, y_view, dparm, iparm) \
     * The Python callback fills y_view based on x_view input. \
     * Since these are views, changes are directly reflected in PROPACK workspace. \
     */ \
    PyObject* result = PyObject_CallFunction(info->py_callback, "iiiOOOO", \
          transa, m, n, x_view, y_view,  info->py_dparm, info->py_iparm ); \
    \
    /* Check for Python exceptions */ \
    if (result == NULL) { STORE_PYTHON_ERROR(info); } \
    Py_DECREF(result); \
    \
    /* Clean up array views (but not the underlying data!) */ \
    CLEANUP_AND_RETURN(x_view, y_view, gstate); \
}

DEFINE_PROPACK_CALLBACK(s, float, NPY_FLOAT32)
DEFINE_PROPACK_CALLBACK(d, double, NPY_FLOAT64)
DEFINE_PROPACK_CALLBACK(c, PROPACK_CPLXF_TYPE, NPY_COMPLEX64)
DEFINE_PROPACK_CALLBACK(z, PROPACK_CPLX_TYPE, NPY_COMPLEX128)

/*
 * Cleanup callback context
 *
 * Releases all Python object references and resets the context structure.
 */
static void cleanup_propack_context(propack_callback_info_t* info) {
    if (info == NULL) return;

    Py_XDECREF(info->py_callback);
    Py_XDECREF(info->py_dparm);
    Py_XDECREF(info->py_iparm);

    /* Clean up stored exception information */
    Py_XDECREF(info->error_type);
    Py_XDECREF(info->error_value);
    Py_XDECREF(info->error_traceback);

    /* Reset structure to avoid leftovers */
    info->py_callback = NULL;
    info->py_dparm = NULL;
    info->py_iparm = NULL;
    info->error_type = NULL;
    info->error_value = NULL;
    info->error_traceback = NULL;
    info->m = 0;
    info->n = 0;
    info->error_occurred = 0;
}


/* =================================================== */
/* ============= PROPACK WRAPPERS ==================== */
/* =================================================== */


// jobu, jobv, m, n, k, aprod, u, v, tol, *works, doption, ioption, dparm, iparm

/*
    slansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_s aprod,
            float* U, int ldu, float* sigma, float* bnd, float* V, int ldv,
            float tolin, float* work, int lwork, int* iwork,
            float* doption, int* ioption, int* info, float* dparm, int* iparm,
            uint64_t* rng_state)
*/


// Single precision IRL SVD wrapper function slansvd
static PyObject*
propack_slansvd(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int jobu, jobv, k, kmax, m, n, propack_info;
    float tol;
    PyObject* py_callback;
    PyArrayObject *U, *V, *work, *iwork, *doption, *ioption, *sparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiifOO!O!O!O!O!O!O!O!O!",
                          &jobu, &jobv, &m, &n, &k, &kmax,                       // iiiiii
                          &tol,                                                  // f
                          &py_callback,                                          // O
                          &PyArray_Type, &U,                                     // O!
                          &PyArray_Type, &V,                                     // O!
                          &PyArray_Type, &work,                                  // O!
                          &PyArray_Type, &iwork,                                 // O!
                          &PyArray_Type, &doption,                               // O! 
                          &PyArray_Type, &ioption,                               // O!
                          &PyArray_Type, &sparm,                                 // O!
                          &PyArray_Type, &iparm,                                 // O!    
                          &PyArray_Type, &ap_rng_state                           // O!
                        ))
    {
        return NULL;
    }

    if (!PyCallable_Check(py_callback)) { PyErr_SetString(PyExc_TypeError, "Callback must be callable"); return NULL; }

    // 0-Initialize callback context structure
    propack_callback_info_t info = {0};

    // Store Python objects in context and increment refcount
    info.py_callback = py_callback;
    Py_INCREF(py_callback);

    info.py_dparm = (PyObject*)sparm;
    Py_INCREF((PyObject*)sparm);

    info.py_iparm = (PyObject*)iparm;
    Py_INCREF((PyObject*)iparm);

    // Store additional context fields
    info.m = m;
    info.n = n;
    info.error_occurred = 0;

    // Set module-level callback context
    _active_callback_info = &info;

    // Create result arrays for singular values and bounds
    npy_intp result_size = neig;
    PyArrayObject* sigma = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_FLOAT32);
    PyArrayObject* bnd = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_FLOAT32);

    if (!sigma || !bnd) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate result arrays");
        Py_XDECREF((PyObject*)sigma);
        Py_XDECREF((PyObject*)bnd);
        _active_callback_info = NULL;
        cleanup_propack_context(&info);
        return NULL;
    }

    // Call PROPACK slansvd_irl function
    slansvd(
        jobu, jobv,                                       // whether to compute U, V
        m, n, k, kmax,                                    // matrix and algorithm dimensions
        propack_callback_s,                               // Py callback function
        (float*)PyArray_DATA(U), PyArray_DIM(U, 1),       // U matrix and leading dimension
        (float*)PyArray_DATA(sigma),                      // singular values output
        (float*)PyArray_DATA(bnd),                        // error bounds output
        (float*)PyArray_DATA(V), PyArray_DIM(V, 1),       // V matrix and leading dimension
        tol,                                              // tolerance
        (float*)PyArray_DATA(work), PyArray_SIZE(work),   // main workspace
        (int*)PyArray_DATA(iwork),                        // integer workspace
        (float*)PyArray_DATA(doption),                    // float options array
        (int*)PyArray_DATA(ioption),                      // integer options array
        &propack_info,                                    // return code
        (float*)PyArray_DATA(sparm),                      // float parameter array (unused)
        (int*)PyArray_DATA(iparm),                        // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)             // random number state
    );

    // Clear module-level callback context
    _active_callback_info = NULL;

    // Check for Python callback error
    if (info.error_occurred) {
        // Restore exception from the callback
        if (info.error_type) {
            PyErr_Restore(info.error_type, info.error_value, info.error_traceback);
            info.error_type = NULL;
            info.error_value = NULL;
            info.error_traceback = NULL;
        } else {
            PyErr_SetString(PropackError, "Error occurred in Python callback");
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    if (propack_info != 0)
    {
        if (propack_info > 0) {
            PyErr_Format(PropackError, "PROPACK found invariant subspace of dimension %d", propack_info);
        } else {
            PyErr_Format(PropackError, "PROPACK failed to converge (info=%d)", propack_info);
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    // Clean up callback context
    cleanup_propack_context(&info);
    return Py_BuildValue("(NNNNI)", PyArray_Return(U), PyArray_Return(sigma), PyArray_Return(V), PyArray_Return(bnd), propack_info);
}





// Single precision IRL SVD wrapper function slansvd_irl
static PyObject*
propack_slansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    float tol;
    PyObject* py_callback;
    PyArrayObject *U, *V, *work, *iwork, *doption, *ioption, *sparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiifOO!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_callback,
                         &PyArray_Type, &U,
                         &PyArray_Type, &V,
                         &PyArray_Type, &work,
                         &PyArray_Type, &iwork,
                         &PyArray_Type, &doption,
                         &PyArray_Type, &ioption,
                         &PyArray_Type, &sparm,
                         &PyArray_Type, &iparm,
                         &PyArray_Type, &ap_rng_state)) {
        return NULL;
    }

    if (!PyCallable_Check(py_callback)) {PyErr_SetString(PyExc_TypeError, "Callback must be callable"); return NULL; }

    // 0-Initialize callback context structure
    propack_callback_info_t info = {0};

    // Store Python objects in context and increment refcount
    info.py_callback = py_callback;
    Py_INCREF(py_callback);

    info.py_dparm = (PyObject*)sparm;
    Py_INCREF((PyObject*)sparm);

    info.py_iparm = (PyObject*)iparm;
    Py_INCREF((PyObject*)iparm);

    // Store additional context fields
    info.m = m;
    info.n = n;
    info.error_occurred = 0;

    // Set module-level callback context
    _active_callback_info = &info;

    int dim = shifts + neig;

    // Create result arrays for singular values and bounds
    npy_intp result_size = neig;
    PyArrayObject* sigma = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_FLOAT32);
    PyArrayObject* bnd = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_FLOAT32);

    if (!sigma || !bnd) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate result arrays");
        Py_XDECREF((PyObject*)sigma);
        Py_XDECREF((PyObject*)bnd);
        _active_callback_info = NULL;
        cleanup_propack_context(&info);
        return NULL;
    }

    // Call PROPACK slansvd_irl function
    slansvd_irl(
        which,                                            // which singular values to compute
        jobu, jobv,                                       // whether to compute U, V
        m, n, dim, shifts,                                // matrix and algorithm dimensions
        &neig,                                            // number of converged values (input/output)
        maxiter,                                          // maximum iterations
        propack_callback_s,                               // Py callback function
        (float*)PyArray_DATA(U), PyArray_DIM(U, 1),       // U matrix and leading dimension
        (float*)PyArray_DATA(sigma),                      // singular values output
        (float*)PyArray_DATA(bnd),                        // error bounds output
        (float*)PyArray_DATA(V), PyArray_DIM(V, 1),       // V matrix and leading dimension
        tol,                                              // tolerance
        (float*)PyArray_DATA(work), PyArray_SIZE(work),   // main workspace
        (int*)PyArray_DATA(iwork),                        // integer workspace
        (float*)PyArray_DATA(doption),                    // float options array
        (int*)PyArray_DATA(ioption),                      // integer options array
        &propack_info,                                    // return code
        (float*)PyArray_DATA(sparm),                      // float parameter array (unused)
        (int*)PyArray_DATA(iparm),                        // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)             // random number state
    );

    // Clear module-level callback context
    _active_callback_info = NULL;

    // Check for Python callback error
    if (info.error_occurred) {
        // Restore exception from the callback
        if (info.error_type) {
            PyErr_Restore(info.error_type, info.error_value, info.error_traceback);
            info.error_type = NULL;
            info.error_value = NULL;
            info.error_traceback = NULL;
        } else {
            PyErr_SetString(PropackError, "Error occurred in Python callback");
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    if (propack_info != 0)
    {
        if (propack_info > 0) {
            PyErr_Format(PropackError, "PROPACK found invariant subspace of dimension %d", propack_info);
        } else {
            PyErr_Format(PropackError, "PROPACK failed to converge (info=%d)", propack_info);
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    // Clean up callback context
    cleanup_propack_context(&info);
    return Py_BuildValue("(NNNNI)", PyArray_Return(U), PyArray_Return(sigma), PyArray_Return(V), PyArray_Return(bnd), propack_info);
}


// Double precision IRL SVD wrapper function dlansvd_irl
static PyObject*
propack_dlansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    double tol;
    PyObject* py_callback;
    PyArrayObject *U, *V, *work, *iwork, *doption, *ioption, *dparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiidOO!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_callback,
                         &PyArray_Type, &U,
                         &PyArray_Type, &V,
                         &PyArray_Type, &work,
                         &PyArray_Type, &iwork,
                         &PyArray_Type, &doption,
                         &PyArray_Type, &ioption,
                         &PyArray_Type, &dparm,
                         &PyArray_Type, &iparm,
                         &PyArray_Type, &ap_rng_state)) {
        return NULL;
    }

    if (!PyCallable_Check(py_callback)) {PyErr_SetString(PyExc_TypeError, "Callback must be callable"); return NULL; }

    // 0-Initialize callback context structure
    propack_callback_info_t info = {0};

    // Store Python objects in context and increment refcount
    info.py_callback = py_callback;
    Py_INCREF(py_callback);

    info.py_dparm = (PyObject*)dparm;
    Py_INCREF((PyObject*)dparm);

    info.py_iparm = (PyObject*)iparm;
    Py_INCREF((PyObject*)iparm);

    // Store additional context fields
    info.m = m;
    info.n = n;
    info.error_occurred = 0;

    // Set module-level callback context
    _active_callback_info = &info;

    int dim = shifts + neig;

    // Create result arrays for singular values and bounds
    npy_intp result_size = neig;
    PyArrayObject* sigma = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_DOUBLE);
    PyArrayObject* bnd = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_DOUBLE);

    if (!sigma || !bnd) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate result arrays");
        Py_XDECREF((PyObject*)sigma);
        Py_XDECREF((PyObject*)bnd);
        _active_callback_info = NULL;
        cleanup_propack_context(&info);
        return NULL;
    }

    // Call PROPACK dlansvd_irl function
    dlansvd_irl(
        which,                                            // which singular values to compute
        jobu, jobv,                                       // whether to compute U, V
        m, n, dim, shifts,                                // matrix and algorithm dimensions
        &neig,                                            // number of converged values (input/output)
        maxiter,                                          // maximum iterations
        propack_callback_d,                               // Py callback function
        (double*)PyArray_DATA(U), PyArray_DIM(U, 1),      // U matrix and leading dimension
        (double*)PyArray_DATA(sigma),                     // singular values output
        (double*)PyArray_DATA(bnd),                       // error bounds output
        (double*)PyArray_DATA(V), PyArray_DIM(V, 1),      // V matrix and leading dimension
        tol,                                              // tolerance
        (double*)PyArray_DATA(work), PyArray_SIZE(work),  // main workspace
        (int*)PyArray_DATA(iwork),                        // integer workspace
        (double*)PyArray_DATA(doption),                   // double options array
        (int*)PyArray_DATA(ioption),                      // integer options array
        &propack_info,                                    // return code
        (double*)PyArray_DATA(dparm),                     // double parameter array (unused)
        (int*)PyArray_DATA(iparm),                        // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)             // random number state
    );

    // Clear module-level callback context
    _active_callback_info = NULL;

    // Check for Python callback error
    if (info.error_occurred) {
        // Restore exception from the callback
        if (info.error_type) {
            PyErr_Restore(info.error_type, info.error_value, info.error_traceback);
            info.error_type = NULL;
            info.error_value = NULL;
            info.error_traceback = NULL;
        } else {
            PyErr_SetString(PropackError, "Error occurred in Python callback");
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    if (propack_info != 0)
    {
        if (propack_info > 0) {
            PyErr_Format(PropackError, "PROPACK found invariant subspace of dimension %d", propack_info);
        } else {
            PyErr_Format(PropackError, "PROPACK failed to converge (info=%d)", propack_info);
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    // Clean up callback context
    cleanup_propack_context(&info);
    return Py_BuildValue("(NNNNI)", PyArray_Return(U), PyArray_Return(sigma), PyArray_Return(V), PyArray_Return(bnd), propack_info);
}


// Single precision complex IRL SVD wrapper function clansvd_irl
static PyObject*
propack_clansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args) {

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    float tol;
    PyObject* py_callback;
    PyArrayObject *U, *V, *work, *cwork, *iwork, *doption, *ioption, *cparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiifOO!O!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_callback,
                         &PyArray_Type, &U,
                         &PyArray_Type, &V,
                         &PyArray_Type, &work,
                         &PyArray_Type, &cwork,
                         &PyArray_Type, &iwork,
                         &PyArray_Type, &doption,
                         &PyArray_Type, &ioption,
                         &PyArray_Type, &cparm,
                         &PyArray_Type, &iparm,
                         &PyArray_Type, &ap_rng_state)) {
        return NULL;
    }

    if (!PyCallable_Check(py_callback)) {PyErr_SetString(PyExc_TypeError, "Callback must be callable"); return NULL; }

    // 0-Initialize callback context structure
    propack_callback_info_t info = {0};

    // Store Python objects in context and increment refcount
    info.py_callback = py_callback;
    Py_INCREF(py_callback);

    info.py_dparm = (PyObject*)cparm;
    Py_INCREF((PyObject*)cparm);

    info.py_iparm = (PyObject*)iparm;
    Py_INCREF((PyObject*)iparm);

    // Store additional context fields
    info.m = m;
    info.n = n;
    info.error_occurred = 0;

    // Set module-level callback context
    _active_callback_info = &info;

    int dim = shifts + neig;

    // Create result arrays for singular values and bounds
    npy_intp result_size = neig;
    PyArrayObject* sigma = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_FLOAT32);
    PyArrayObject* bnd = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_FLOAT32);

    if (!sigma || !bnd) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate result arrays");
        Py_XDECREF((PyObject*)sigma);
        Py_XDECREF((PyObject*)bnd);
        _active_callback_info = NULL;
        cleanup_propack_context(&info);
        return NULL;
    }

    // Call PROPACK clansvd_irl function
    clansvd_irl(
        which,                                                         // which singular values to compute
        jobu, jobv,                                                    // whether to compute U, V
        m, n, dim, shifts,                                             // matrix and algorithm dimensions
        &neig,                                                         // number of converged values (input/output)
        maxiter,                                                       // maximum iterations
        propack_callback_c,                                            // Py callback function
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(U), PyArray_DIM(U, 1),       // U matrix and leading dimension
        (float*)PyArray_DATA(sigma),                                   // singular values output
        (float*)PyArray_DATA(bnd),                                     // error bounds output
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(V), PyArray_DIM(V, 1),       // V matrix and leading dimension
        tol,                                                           // tolerance
        (float*)PyArray_DATA(work), PyArray_SIZE(work),                // float workspace
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(cwork), PyArray_SIZE(cwork), // float complex workspace
        (int*)PyArray_DATA(iwork),                                     // integer workspace
        (float*)PyArray_DATA(doption),                                 // float options array
        (int*)PyArray_DATA(ioption),                                   // integer options array
        &propack_info,                                                 // return code
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(cparm),                      // float parameter array (unused)
        (int*)PyArray_DATA(iparm),                                     // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)                          // random number state
    );

    // Clear module-level callback context
    _active_callback_info = NULL;

    // Check for Python callback error
    if (info.error_occurred) {
        // Restore exception from the callback
        if (info.error_type) {
            PyErr_Restore(info.error_type, info.error_value, info.error_traceback);
            info.error_type = NULL;
            info.error_value = NULL;
            info.error_traceback = NULL;
        } else {
            PyErr_SetString(PropackError, "Error occurred in Python callback");
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    if (propack_info != 0)
    {
        if (propack_info > 0) {
            PyErr_Format(PropackError, "PROPACK found invariant subspace of dimension %d", propack_info);
        } else {
            PyErr_Format(PropackError, "PROPACK failed to converge (info=%d)", propack_info);
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    // Clean up callback context
    cleanup_propack_context(&info);
    return Py_BuildValue("(NNNNI)", PyArray_Return(U), PyArray_Return(sigma), PyArray_Return(V), PyArray_Return(bnd), propack_info);
}


// Double precision complex IRL SVD wrapper function zlansvd_irl
static PyObject*
propack_zlansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args) {

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    double tol;
    PyObject* py_callback;
    PyArrayObject *U, *V, *work, *cwork, *iwork, *doption, *ioption, *zparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiidOO!O!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_callback,
                         &PyArray_Type, &U,
                         &PyArray_Type, &V,
                         &PyArray_Type, &work,
                         &PyArray_Type, &cwork,
                         &PyArray_Type, &iwork,
                         &PyArray_Type, &doption,
                         &PyArray_Type, &ioption,
                         &PyArray_Type, &zparm,
                         &PyArray_Type, &iparm,
                         &PyArray_Type, &ap_rng_state)) {
        return NULL;
    }

    if (!PyCallable_Check(py_callback)) {
        PyErr_SetString(PyExc_TypeError, "Callback must be callable");
        return NULL;
    }

    // 0-Initialize callback context structure
    propack_callback_info_t info = {0};

    // Store Python objects in context and increment refcount
    info.py_callback = py_callback;
    Py_INCREF(py_callback);

    info.py_dparm = (PyObject*)zparm;
    Py_INCREF((PyObject*)zparm);

    info.py_iparm = (PyObject*)iparm;
    Py_INCREF((PyObject*)iparm);

    // Store additional context fields
    info.m = m;
    info.n = n;
    info.error_occurred = 0;

    // Set module-level callback context
    _active_callback_info = &info;

    int dim = shifts + neig;

    // Create result arrays for singular values and bounds
    npy_intp result_size = neig;
    PyArrayObject* sigma = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_DOUBLE);
    PyArrayObject* bnd = (PyArrayObject*)PyArray_SimpleNew(1, &result_size, NPY_DOUBLE);

    if (!sigma || !bnd) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate result arrays");
        Py_XDECREF((PyObject*)sigma);
        Py_XDECREF((PyObject*)bnd);
        _active_callback_info = NULL;
        cleanup_propack_context(&info);
        return NULL;
    }

    // Call PROPACK zlansvd_irl function
    zlansvd_irl(
        which,                                                         // which singular values to compute
        jobu, jobv,                                                    // whether to compute U, V
        m, n, dim, shifts,                                             // matrix and algorithm dimensions
        &neig,                                                         // number of converged values (input/output)
        maxiter,                                                       // maximum iterations
        propack_callback_z,                                            // our callback function
        (PROPACK_CPLX_TYPE*)PyArray_DATA(U), PyArray_DIM(U, 1),        // U matrix and leading dimension
        (double*)PyArray_DATA(sigma),                                  // singular values output
        (double*)PyArray_DATA(bnd),                                    // error bounds output
        (PROPACK_CPLX_TYPE*)PyArray_DATA(V), PyArray_DIM(V, 1),        // V matrix and leading dimension
        tol,                                                           // tolerance
        (double*)PyArray_DATA(work), PyArray_SIZE(work),               // double workspace
        (PROPACK_CPLX_TYPE*)PyArray_DATA(cwork), PyArray_SIZE(cwork),  // double complex workspace
        (int*)PyArray_DATA(iwork),                                     // integer workspace
        (double*)PyArray_DATA(doption),                                // double options array
        (int*)PyArray_DATA(ioption),                                   // integer options array
        &propack_info,                                                 // return code
        (PROPACK_CPLX_TYPE*)PyArray_DATA(zparm),                       // double parameter array (unused)
        (int*)PyArray_DATA(iparm),                                     // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)                          // random number state
    );

    // Clear module-level callback context
    _active_callback_info = NULL;

    // Check for Python callback error
    if (info.error_occurred) {
        // Restore exception from the callback
        if (info.error_type) {
            PyErr_Restore(info.error_type, info.error_value, info.error_traceback);
            info.error_type = NULL;
            info.error_value = NULL;
            info.error_traceback = NULL;
        } else {
            PyErr_SetString(PropackError, "Error occurred in Python callback");
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    if (propack_info != 0) {
        if (propack_info > 0) {
            PyErr_Format(PropackError, "PROPACK found invariant subspace of dimension %d", propack_info);
        } else {
            PyErr_Format(PropackError, "PROPACK failed to converge (info=%d)", propack_info);
        }
        Py_DECREF((PyObject*)sigma);
        Py_DECREF((PyObject*)bnd);
        cleanup_propack_context(&info);
        return NULL;
    }

    // Clean up callback context
    cleanup_propack_context(&info);
    return Py_BuildValue("(NNNNI)", PyArray_Return(U), PyArray_Return(sigma), PyArray_Return(V), PyArray_Return(bnd), propack_info);
}


/*
 * Module method table
 *
 * Currently contains only dlansvd_irl for demonstration.
 * Additional functions (slansvd_irl, clansvd_irl, zlansvd_irl, and
 * non-IRL versions) would be added here following the same pattern.
 */
static PyMethodDef propack_methods[] = {
    {"slansvd_irl", propack_slansvd_irl, METH_VARARGS, "Single precision implicitly restarted Lanczos SVD"},
    {"dlansvd_irl", propack_dlansvd_irl, METH_VARARGS, "Double precision implicitly restarted Lanczos SVD"},
    {"clansvd_irl", propack_clansvd_irl, METH_VARARGS, "Single precision complex implicitly restarted Lanczos SVD"},
    {"zlansvd_irl", propack_clansvd_irl, METH_VARARGS, "Double precision complex implicitly restarted Lanczos SVD"},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};


static int propack_exec(PyObject* module) {
    // Initialize NumPy C API; import_array()
    if (PyArray_ImportNumPyAPI() < 0) { return -1; }

    if (PyErr_Occurred()) { return -1; }

    PropackError = PyErr_NewException("_propack.PropackError", NULL, NULL);
    if (PropackError == NULL) { return -1; }

    // Add exception to module
    if (PyModule_AddObject(module, "PropackError", PropackError) < 0) {
        Py_DECREF(PropackError);
        return -1;
    }
    return 0;
}

// Module slots for multi-phase initialization
static PyModuleDef_Slot propack_slots[] = {
    {Py_mod_exec, propack_exec},
    {0, NULL}
};


// Module definition structure using designated initializers
static struct PyModuleDef propack_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_propack",
    .m_doc = "PROPACK SVD routines",
    .m_size = 0,  /* Stateless module */
    .m_methods = propack_methods,
    .m_slots = propack_slots,
    .m_traverse = NULL,
    .m_clear = NULL,
    .m_free = NULL
};


PyMODINIT_FUNC PyInit__propack(void) {
    return PyModuleDef_Init(&propack_module);
}
