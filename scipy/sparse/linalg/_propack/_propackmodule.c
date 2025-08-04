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
#include "ccallback.h"
#include <stdint.h>

/* Module-level error object */
static PyObject* PropackError;


// Define the callback signature for PROPACK
static ccallback_signature_t propack_signatures[] = {
    {"void (int, int, int, float *, float *, float *, int *)", 0},
    {NULL, 0}  // Sentinel
};


typedef struct {
    ccallback_t callback;   // Callback struct
    PyObject *py_dparm;     // Unused, but kept for compatibility
    PyObject *py_iparm;     // Unused, but kept for compatibility
    PyObject *args_tuple;   // Pre-allocated args tuple for efficiency
    int m, n;               // Matrix dimensions
} propack_ccallback_info_t;


// Thread-local storage for callback context
static SCIPY_TLS propack_ccallback_info_t* current_propack_callback = NULL;


static void
propack_callback_s_thunk(int transa, int m, int n, float* x, float* y, float* dparm, int* iparm)
{
    // Silence unused parameter warnings for dparm/iparm
    (void)dparm; (void)iparm;

    if ((!current_propack_callback) || (!current_propack_callback->callback.py_function)) { return; }

    npy_intp x_dims[1] = {n};
    npy_intp y_dims[1] = {m};

    // Create PyCapsules for memory management
    PyObject *x_capsule = PyCapsule_New((void*)x, NULL, NULL);
    PyObject *y_capsule = PyCapsule_New((void*)y, NULL, NULL);
    if (!x_capsule || !y_capsule) {
        Py_XDECREF(x_capsule);
        Py_XDECREF(y_capsule);
        return;
    }

    // Create array views
    PyObject* py_x = PyArray_SimpleNewFromData(1, x_dims, NPY_FLOAT32, (void*)x);
    PyObject* py_y = PyArray_SimpleNewFromData(1, y_dims, NPY_FLOAT32, (void*)y);

    if (!py_x || !py_y) {
        Py_XDECREF(x_capsule);
        Py_XDECREF(y_capsule);
        Py_XDECREF(py_x);
        Py_XDECREF(py_y);
        return;
    }

    // Set base objects
    if (PyArray_SetBaseObject((PyArrayObject*)py_x, x_capsule) < 0 ||
        PyArray_SetBaseObject((PyArrayObject*)py_y, y_capsule) < 0) {
        // capsule references stolen even on failure
        Py_XDECREF(py_x);
        Py_XDECREF(py_y);
        return;
    }

    // Update the pre-allocated tuple with changing values
    // Only transa, py_x, py_y change between calls
    PyTuple_SetItem(current_propack_callback->args_tuple, 0, PyLong_FromLong(transa));
    PyTuple_SetItem(current_propack_callback->args_tuple, 3, py_x);
    PyTuple_SetItem(current_propack_callback->args_tuple, 4, py_y);

    // Call Python function
    PyObject *result = PyObject_CallObject(current_propack_callback->callback.py_function,
                                           current_propack_callback->args_tuple);

    Py_XDECREF(result);
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

    // 0-Initialize callback context structure
    propack_callback_info_t info = {0};

    // Prepare ccallback with signatures and PARSE flag for legacy support
    if (ccallback_prepare(&info.callback, propack_signatures, py_callback, CCALLBACK_PARSE) != 0) { return NULL; }



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
