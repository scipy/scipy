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


typedef struct {
    PyObject *py_func;
    PyObject *args_tuple;
} propack_callback_t;


// Thread-local storage for callback context
static SCIPY_TLS propack_callback_t* current_propack_callback = NULL;


static void
propack_callback_s_thunk(int transa, int m, int n, float* x, float* y, float* dparm, int* iparm)
{
    // Silence unused parameter warnings for dparm/iparm
    (void)dparm; (void)iparm;

    if ((!current_propack_callback) || (!current_propack_callback->py_func)) { return; }

    npy_intp x_dims[1] = {transa ? m : n};
    npy_intp y_dims[1] = {transa ? n : m};

    // Create PyCapsules, no need for destructors since we use them only for passing data
    PyObject *x_capsule = PyCapsule_New((void*)x, NULL, NULL);
    PyObject *y_capsule = PyCapsule_New((void*)y, NULL, NULL);
    if ((!x_capsule) || (!y_capsule)) {
        Py_XDECREF(x_capsule);
        Py_XDECREF(y_capsule);
        return;
    }

    // Create array views
    PyObject* py_x = PyArray_SimpleNewFromData(1, x_dims, NPY_FLOAT32, (void*)x);
    PyObject* py_y = PyArray_SimpleNewFromData(1, y_dims, NPY_FLOAT32, (void*)y);

    if ((!py_x) || (!py_y)) {
        Py_XDECREF(x_capsule);Py_XDECREF(y_capsule);
        Py_XDECREF(py_x);Py_XDECREF(py_y);
        return;
    }

    // Set base objects
    if ((PyArray_SetBaseObject((PyArrayObject*)py_x, x_capsule) < 0) ||
        (PyArray_SetBaseObject((PyArrayObject*)py_y, y_capsule) < 0)) {
        Py_XDECREF(py_x);
        Py_XDECREF(py_y);
        return;
    }

    PyTuple_SetItem(current_propack_callback->args_tuple, 0, PyLong_FromLong(transa));
    PyTuple_SetItem(current_propack_callback->args_tuple, 1, PyLong_FromLong(m));
    PyTuple_SetItem(current_propack_callback->args_tuple, 2, PyLong_FromLong(n));
    PyTuple_SetItem(current_propack_callback->args_tuple, 3, py_x);
    PyTuple_SetItem(current_propack_callback->args_tuple, 4, py_y);

    // Remaining arguments dparm, iparm are not used in SciPy.

    // Call Python function, all operations are done in-place hence result is discarded.
    PyObject* result = PyObject_CallObject(current_propack_callback->py_func,
                                           current_propack_callback->args_tuple);

    Py_XDECREF(result);
}


static void
propack_callback_d_thunk(int transa, int m, int n, double* x, double* y, double* dparm, int* iparm)
{
    // Silence unused parameter warnings for dparm/iparm
    (void)dparm; (void)iparm;

    if ((!current_propack_callback) || (!current_propack_callback->py_func)) { return; }

    npy_intp x_dims[1] = {transa ? m : n};
    npy_intp y_dims[1] = {transa ? n : m};

    // Create PyCapsules, no need for destructors since we use them only for passing data
    PyObject *x_capsule = PyCapsule_New((void*)x, NULL, NULL);
    PyObject *y_capsule = PyCapsule_New((void*)y, NULL, NULL);
    if ((!x_capsule) || (!y_capsule)) {
        Py_XDECREF(x_capsule);
        Py_XDECREF(y_capsule);
        return;
    }

    // Create array views
    PyObject* py_x = PyArray_SimpleNewFromData(1, x_dims, NPY_FLOAT64, (void*)x);
    PyObject* py_y = PyArray_SimpleNewFromData(1, y_dims, NPY_FLOAT64, (void*)y);

    if ((!py_x) || (!py_y)) {
        Py_XDECREF(x_capsule);Py_XDECREF(y_capsule);
        Py_XDECREF(py_x);Py_XDECREF(py_y);
        return;
    }

    // Set base objects
    if ((PyArray_SetBaseObject((PyArrayObject*)py_x, x_capsule) < 0) ||
        (PyArray_SetBaseObject((PyArrayObject*)py_y, y_capsule) < 0)) {
        // capsule references stolen even on failure
        Py_XDECREF(py_x);
        Py_XDECREF(py_y);
        return;
    }

    PyTuple_SetItem(current_propack_callback->args_tuple, 0, PyLong_FromLong(transa));
    PyTuple_SetItem(current_propack_callback->args_tuple, 1, PyLong_FromLong(m));
    PyTuple_SetItem(current_propack_callback->args_tuple, 2, PyLong_FromLong(n));
    PyTuple_SetItem(current_propack_callback->args_tuple, 3, py_x);
    PyTuple_SetItem(current_propack_callback->args_tuple, 4, py_y);

    // Remaining arguments dparm, iparm are not used in SciPy.

    // Call Python function, all operations are done in-place hence result is discarded.
    PyObject *result = PyObject_CallObject(current_propack_callback->py_func,
                                           current_propack_callback->args_tuple);

    Py_XDECREF(result);
}


static void
propack_callback_c_thunk(int transa, int m, int n, PROPACK_CPLXF_TYPE* x, PROPACK_CPLXF_TYPE* y, PROPACK_CPLXF_TYPE* dparm, int* iparm)
{
    // Silence unused parameter warnings for dparm/iparm
    (void)dparm; (void)iparm;

    if ((!current_propack_callback) || (!current_propack_callback->py_func)) { return; }

    npy_intp x_dims[1] = {transa ? m : n};
    npy_intp y_dims[1] = {transa ? n : m};

    // Create PyCapsules, no need for destructors since we use them only for passing data
    PyObject *x_capsule = PyCapsule_New((void*)x, NULL, NULL);
    PyObject *y_capsule = PyCapsule_New((void*)y, NULL, NULL);
    if ((!x_capsule) || (!y_capsule)) {
        Py_XDECREF(x_capsule);
        Py_XDECREF(y_capsule);
        return;
    }

    // Create array views
    PyObject* py_x = PyArray_SimpleNewFromData(1, x_dims, NPY_COMPLEX64, (void*)x);
    PyObject* py_y = PyArray_SimpleNewFromData(1, y_dims, NPY_COMPLEX64, (void*)y);

    if ((!py_x) || (!py_y)) {
        Py_XDECREF(x_capsule);Py_XDECREF(y_capsule);
        Py_XDECREF(py_x);Py_XDECREF(py_y);
        return;
    }

    // Set base objects
    if ((PyArray_SetBaseObject((PyArrayObject*)py_x, x_capsule) < 0) ||
        (PyArray_SetBaseObject((PyArrayObject*)py_y, y_capsule) < 0)) {
        // capsule references stolen even on failure
        Py_XDECREF(py_x);
        Py_XDECREF(py_y);
        return;
    }

    PyTuple_SetItem(current_propack_callback->args_tuple, 0, PyLong_FromLong(transa));
    PyTuple_SetItem(current_propack_callback->args_tuple, 1, PyLong_FromLong(m));
    PyTuple_SetItem(current_propack_callback->args_tuple, 2, PyLong_FromLong(n));
    PyTuple_SetItem(current_propack_callback->args_tuple, 3, py_x);
    PyTuple_SetItem(current_propack_callback->args_tuple, 4, py_y);

    // Remaining arguments dparm, iparm are not used in SciPy.

    // Call Python function, all operations are done in-place hence result is discarded.
    PyObject *result = PyObject_CallObject(current_propack_callback->py_func,
                                           current_propack_callback->args_tuple);

    Py_XDECREF(result);
}


static void
propack_callback_z_thunk(int transa, int m, int n, PROPACK_CPLX_TYPE* x, PROPACK_CPLX_TYPE* y, PROPACK_CPLX_TYPE* dparm, int* iparm)
{
    // Silence unused parameter warnings for dparm/iparm
    (void)dparm; (void)iparm;

    if ((!current_propack_callback) || (!current_propack_callback->py_func)) { return; }

    npy_intp x_dims[1] = {transa ? m : n};
    npy_intp y_dims[1] = {transa ? n : m};

    // Create PyCapsules, no need for destructors since we use them only for passing data
    PyObject *x_capsule = PyCapsule_New((void*)x, NULL, NULL);
    PyObject *y_capsule = PyCapsule_New((void*)y, NULL, NULL);
    if ((!x_capsule) || (!y_capsule)) {
        Py_XDECREF(x_capsule);
        Py_XDECREF(y_capsule);
        return;
    }

    // Create array views
    PyObject* py_x = PyArray_SimpleNewFromData(1, x_dims, NPY_COMPLEX128, (void*)x);
    PyObject* py_y = PyArray_SimpleNewFromData(1, y_dims, NPY_COMPLEX128, (void*)y);

    if ((!py_x) || (!py_y)) {
        Py_XDECREF(x_capsule);Py_XDECREF(y_capsule);
        Py_XDECREF(py_x);Py_XDECREF(py_y);
        return;
    }

    // Set base objects
    if ((PyArray_SetBaseObject((PyArrayObject*)py_x, x_capsule) < 0) ||
        (PyArray_SetBaseObject((PyArrayObject*)py_y, y_capsule) < 0)) {
        // capsule references stolen even on failure
        Py_XDECREF(py_x);
        Py_XDECREF(py_y);
        return;
    }

    PyTuple_SetItem(current_propack_callback->args_tuple, 0, PyLong_FromLong(transa));
    PyTuple_SetItem(current_propack_callback->args_tuple, 1, PyLong_FromLong(m));
    PyTuple_SetItem(current_propack_callback->args_tuple, 2, PyLong_FromLong(n));
    PyTuple_SetItem(current_propack_callback->args_tuple, 3, py_x);
    PyTuple_SetItem(current_propack_callback->args_tuple, 4, py_y);

    // Remaining arguments dparm, iparm are not used in SciPy.

    // Call Python function, all operations are done in-place hence result is discarded.
    PyObject *result = PyObject_CallObject(current_propack_callback->py_func,
                                           current_propack_callback->args_tuple);

    Py_XDECREF(result);
}



/* =================================================== */
/* ============= PROPACK WRAPPERS ==================== */
/* =================================================== */

static PyObject*
propack_slansvd(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int jobu, jobv, k, kmax, m, n, propack_info;
    float tol;
    PyObject* py_aprod;
    PyArrayObject *U, *V, *sigma, *bnd, *work, *iwork, *doption, *ioption, *dparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiifOO!O!O!O!O!O!O!O!O!O!O!",
                          &jobu, &jobv, &m, &n, &k, &kmax,                       // iiiiii
                          &tol,                                                  // f
                          &py_aprod,                                             // O
                          &PyArray_Type, &U,                                     // O!
                          &PyArray_Type, &sigma,                                 // O!
                          &PyArray_Type, &bnd,                                   // O!
                          &PyArray_Type, &V,                                     // O!
                          &PyArray_Type, &work,                                  // O!
                          &PyArray_Type, &iwork,                                 // O!
                          &PyArray_Type, &doption,                               // O!
                          &PyArray_Type, &ioption,                               // O!
                          &PyArray_Type, &dparm,                                 // O! (ignored)
                          &PyArray_Type, &iparm,                                 // O! (ignored)
                          &PyArray_Type, &ap_rng_state                           // O!
                        ))
    {
        return NULL;
    }

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }

    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    // Call PROPACK slansvd_irl function
    slansvd(
        jobu, jobv,                                       // whether to compute U, V
        m, n, k, kmax,                                    // matrix and algorithm dimensions
        (PROPACK_aprod_s)propack_callback_s_thunk,        // Py callback function
        (float*)PyArray_DATA(U), PyArray_DIM(U, 0),       // U matrix and leading dimension
        (float*)PyArray_DATA(sigma),                      // singular values output
        (float*)PyArray_DATA(bnd),                        // error bounds output
        (float*)PyArray_DATA(V), PyArray_DIM(V, 0),       // V matrix and leading dimension
        tol,                                              // tolerance
        (float*)PyArray_DATA(work), PyArray_SIZE(work),   // main workspace
        (int*)PyArray_DATA(iwork),                        // integer workspace
        (float*)PyArray_DATA(doption),                    // float options array
        (int*)PyArray_DATA(ioption),                      // integer options array
        &propack_info,                                    // return code
        (float*)PyArray_DATA(dparm),                      // float parameter array (unused)
        (int*)PyArray_DATA(iparm),                        // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)             // random number state
    );

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);

}


static PyObject*
propack_dlansvd(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int jobu, jobv, k, kmax, m, n, propack_info;
    double tol;
    PyObject* py_aprod;
    PyArrayObject *U, *V, *sigma, *bnd, *work, *iwork, *doption, *ioption, *dparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiidOO!O!O!O!O!O!O!O!O!O!O!",
                          &jobu, &jobv, &m, &n, &k, &kmax,                       // iiiiii
                          &tol,                                                  // d
                          &py_aprod,                                             // O
                          &PyArray_Type, &U,                                     // O!
                          &PyArray_Type, &sigma,                                 // O!
                          &PyArray_Type, &bnd,                                   // O!
                          &PyArray_Type, &V,                                     // O!
                          &PyArray_Type, &work,                                  // O!
                          &PyArray_Type, &iwork,                                 // O!
                          &PyArray_Type, &doption,                               // O!
                          &PyArray_Type, &ioption,                               // O!
                          &PyArray_Type, &dparm,                                 // O! (ignored)
                          &PyArray_Type, &iparm,                                 // O! (ignored)
                          &PyArray_Type, &ap_rng_state                           // O!
                        ))
    {
        return NULL;
    }

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }
    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    dlansvd(
        jobu, jobv,                                       // whether to compute U, V
        m, n, k, kmax,                                    // matrix and algorithm dimensions
        (PROPACK_aprod_d)propack_callback_d_thunk,        // Py callback function
        (double*)PyArray_DATA(U), PyArray_DIM(U, 0),      // U matrix and leading dimension
        (double*)PyArray_DATA(sigma),                     // singular values output
        (double*)PyArray_DATA(bnd),                       // error bounds output
        (double*)PyArray_DATA(V), PyArray_DIM(V, 0),      // V matrix and leading dimension
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

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);

}


static PyObject*
propack_clansvd(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int jobu, jobv, k, kmax, m, n, propack_info;
    float tol;
    PyObject* py_aprod;
    PyArrayObject *U, *V, *sigma, *bnd, *work, *cwork, *iwork, *doption, *ioption, *dparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiifOO!O!O!O!O!O!O!O!O!O!O!O!",
                          &jobu, &jobv, &m, &n, &k, &kmax,                       // iiiiii
                          &tol,                                                  // f
                          &py_aprod,                                             // O
                          &PyArray_Type, &U,                                     // O!
                          &PyArray_Type, &sigma,                                 // O!
                          &PyArray_Type, &bnd,                                   // O!
                          &PyArray_Type, &V,                                     // O!
                          &PyArray_Type, &work,                                  // O!
                          &PyArray_Type, &cwork,                                 // O!
                          &PyArray_Type, &iwork,                                 // O!
                          &PyArray_Type, &doption,                               // O!
                          &PyArray_Type, &ioption,                               // O!
                          &PyArray_Type, &dparm,                                 // O! (ignored)
                          &PyArray_Type, &iparm,                                 // O! (ignored)
                          &PyArray_Type, &ap_rng_state                           // O!
                        ))
    {
        return NULL;
    }

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }
    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    clansvd(
        jobu, jobv,                                                    // whether to compute U, V
        m, n, k, kmax,                                                 // matrix and algorithm dimensions
        (PROPACK_aprod_c)propack_callback_c_thunk,                     // Py callback function
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(U), PyArray_DIM(U, 0),       // U matrix and leading dimension
        (float*)PyArray_DATA(sigma),                                   // singular values output
        (float*)PyArray_DATA(bnd),                                     // error bounds output
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(V), PyArray_DIM(V, 0),       // V matrix and leading dimension
        tol,                                                           // tolerance
        (float*)PyArray_DATA(work), PyArray_SIZE(work),                // main workspace
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(cwork), PyArray_SIZE(cwork), // main workspace
        (int*)PyArray_DATA(iwork),                                     // integer workspace
        (float*)PyArray_DATA(doption),                                 // float options array
        (int*)PyArray_DATA(ioption),                                   // integer options array
        &propack_info,                                                 // return code
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(dparm),                      // double parameter array (unused)
        (int*)PyArray_DATA(iparm),                                     // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)                          // random number state
    );

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);

}


static PyObject*
propack_zlansvd(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int jobu, jobv, k, kmax, m, n, propack_info;
    double tol;
    PyObject* py_aprod;
    PyArrayObject *U, *V, *sigma, *bnd, *work, *cwork, *iwork, *doption, *ioption, *dparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiidOO!O!O!O!O!O!O!O!O!O!O!O!",
                          &jobu, &jobv, &m, &n, &k, &kmax,                       // iiiiii
                          &tol,                                                  // d
                          &py_aprod,                                             // O
                          &PyArray_Type, &U,                                     // O!
                          &PyArray_Type, &sigma,                                 // O!
                          &PyArray_Type, &bnd,                                   // O!
                          &PyArray_Type, &V,                                     // O!
                          &PyArray_Type, &work,                                  // O!
                          &PyArray_Type, &cwork,                                 // O!
                          &PyArray_Type, &iwork,                                 // O!
                          &PyArray_Type, &doption,                               // O!
                          &PyArray_Type, &ioption,                               // O!
                          &PyArray_Type, &dparm,                                 // O! (ignored)
                          &PyArray_Type, &iparm,                                 // O! (ignored)
                          &PyArray_Type, &ap_rng_state                           // O!
                        ))
    {
        return NULL;
    }

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }
    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    zlansvd(
        jobu, jobv,                                                    // whether to compute U, V
        m, n, k, kmax,                                                 // matrix and algorithm dimensions
        (PROPACK_aprod_z)propack_callback_z_thunk,                     // Py callback function
        (PROPACK_CPLX_TYPE*)PyArray_DATA(U), PyArray_DIM(U, 0),        // U matrix and leading dimension
        (double*)PyArray_DATA(sigma),                                  // singular values output
        (double*)PyArray_DATA(bnd),                                    // error bounds output
        (PROPACK_CPLX_TYPE*)PyArray_DATA(V), PyArray_DIM(V, 0),        // V matrix and leading dimension
        tol,                                                           // tolerance
        (double*)PyArray_DATA(work), PyArray_SIZE(work),               // main workspace
        (PROPACK_CPLX_TYPE*)PyArray_DATA(cwork), PyArray_SIZE(cwork),  // main workspace
        (int*)PyArray_DATA(iwork),                                     // integer workspace
        (double*)PyArray_DATA(doption),                                // double options array
        (int*)PyArray_DATA(ioption),                                   // integer options array
        &propack_info,                                                 // return code
        (PROPACK_CPLX_TYPE*)PyArray_DATA(dparm),                       // double parameter array (unused)
        (int*)PyArray_DATA(iparm),                                     // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)                          // random number state
    );

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);

}


// Single precision IRL SVD wrapper function slansvd_irl
static PyObject*
propack_slansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    float tol;
    PyObject* py_aprod;
    PyArrayObject *U, *sigma, *bnd, *V, *work, *iwork, *doption, *ioption, *sparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiifOO!O!O!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_aprod,
                         &PyArray_Type, &U,
                         &PyArray_Type, &sigma,
                         &PyArray_Type, &bnd,
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

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }
    int dim = shifts + neig;

    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    slansvd_irl(
        which,                                            // which singular values to compute
        jobu, jobv,                                       // whether to compute U, V
        m, n, dim, shifts,                                // matrix and algorithm dimensions
        &neig,                                            // number of converged values (input/output)
        maxiter,                                          // maximum iterations
        (PROPACK_aprod_s)propack_callback_s_thunk,        // Py callback function
        (float*)PyArray_DATA(U), PyArray_DIM(U, 0),       // U matrix and leading dimension
        (float*)PyArray_DATA(sigma),                      // singular values output
        (float*)PyArray_DATA(bnd),                        // error bounds output
        (float*)PyArray_DATA(V), PyArray_DIM(V, 0),       // V matrix and leading dimension
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

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);

}

// Double precision IRL SVD wrapper function dlansvd_irl
static PyObject*
propack_dlansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    double tol;
    PyObject* py_aprod;
    PyArrayObject *U, *sigma, *bnd, *V, *work, *iwork, *doption, *ioption, *sparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiidOO!O!O!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_aprod,
                         &PyArray_Type, &U,
                         &PyArray_Type, &sigma,
                         &PyArray_Type, &bnd,
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

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }
    int dim = shifts + neig;
    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    dlansvd_irl(
        which,                                             // which singular values to compute
        jobu, jobv,                                        // whether to compute U, V
        m, n, dim, shifts,                                 // matrix and algorithm dimensions
        &neig,                                             // number of converged values (input/output)
        maxiter,                                           // maximum iterations
        (PROPACK_aprod_d)propack_callback_d_thunk,         // Py callback function
        (double*)PyArray_DATA(U), PyArray_DIM(U, 0),       // U matrix and leading dimension
        (double*)PyArray_DATA(sigma),                      // singular values output
        (double*)PyArray_DATA(bnd),                        // error bounds output
        (double*)PyArray_DATA(V), PyArray_DIM(V, 0),       // V matrix and leading dimension
        tol,                                               // tolerance
        (double*)PyArray_DATA(work), PyArray_SIZE(work),   // main workspace
        (int*)PyArray_DATA(iwork),                         // integer workspace
        (double*)PyArray_DATA(doption),                    // double options array
        (int*)PyArray_DATA(ioption),                       // integer options array
        &propack_info,                                     // return code
        (double*)PyArray_DATA(sparm),                      // double parameter array (unused)
        (int*)PyArray_DATA(iparm),                         // integer parameter array (unused)
        (uint64_t*)PyArray_DATA(ap_rng_state)              // random number state
    );

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);

}


// Single precision complex IRL SVD wrapper function clansvd_irl
static PyObject*
propack_clansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args) {

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    float tol;
    PyObject* py_aprod;
    PyArrayObject *U, *sigma, *bnd, *V, *work, *cwork, *iwork, *doption, *ioption, *cparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiifOO!O!O!O!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_aprod,
                         &PyArray_Type, &U,
                         &PyArray_Type, &sigma,
                         &PyArray_Type, &bnd,
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

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }
    int dim = shifts + neig;
    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    // Call PROPACK clansvd_irl function
    clansvd_irl(
        which,                                                         // which singular values to compute
        jobu, jobv,                                                    // whether to compute U, V
        m, n, dim, shifts,                                             // matrix and algorithm dimensions
        &neig,                                                         // number of converged values (input/output)
        maxiter,                                                       // maximum iterations
        (PROPACK_aprod_c)propack_callback_c_thunk,                     // Py callback function
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(U), PyArray_DIM(U, 0),       // U matrix and leading dimension
        (float*)PyArray_DATA(sigma),                                   // singular values output
        (float*)PyArray_DATA(bnd),                                     // error bounds output
        (PROPACK_CPLXF_TYPE*)PyArray_DATA(V), PyArray_DIM(V, 0),       // V matrix and leading dimension
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

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);
}


// Double precision complex IRL SVD wrapper function zlansvd_irl
static PyObject*
propack_zlansvd_irl(PyObject* Py_UNUSED(dummy), PyObject* args) {

    int which, jobu, jobv, m, n, shifts, neig, maxiter, propack_info;
    double tol;
    PyObject* py_aprod;
    PyArrayObject *U, *sigma, *bnd, *V, *work, *cwork, *iwork, *doption, *ioption, *zparm, *iparm, *ap_rng_state;

    if (!PyArg_ParseTuple(args, "iiiiiiiidOO!O!O!O!O!O!O!O!O!O!O!O!",
                         &which, &jobu, &jobv, &m, &n, &shifts, &neig, &maxiter, &tol,
                         &py_aprod,
                         &PyArray_Type, &U,
                         &PyArray_Type, &sigma,
                         &PyArray_Type, &bnd,
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

    if (!PyCallable_Check(py_aprod)) {PyErr_SetString(PyExc_TypeError, "scipy.sparse.linalg: Callback argument must be a callable."); return NULL; }
    int dim = shifts + neig;
    propack_callback_t aprod_callback;
    aprod_callback.py_func = py_aprod;
    Py_XINCREF(py_aprod);
    PyObject *args_tuple = PyTuple_New(5);
    // Set all entries to None. PyTuple_SetItem steals references, so we need to INCREF None
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 0, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 1, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 2, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 3, Py_None);
    Py_INCREF(Py_None); PyTuple_SetItem(args_tuple, 4, Py_None);
    aprod_callback.args_tuple = args_tuple;
    current_propack_callback = &aprod_callback;

    // Call PROPACK zlansvd_irl function
    zlansvd_irl(
        which,                                                         // which singular values to compute
        jobu, jobv,                                                    // whether to compute U, V
        m, n, dim, shifts,                                             // matrix and algorithm dimensions
        &neig,                                                         // number of converged values (input/output)
        maxiter,                                                       // maximum iterations
        (PROPACK_aprod_z)propack_callback_z_thunk,                     // our callback function
        (PROPACK_CPLX_TYPE*)PyArray_DATA(U), PyArray_DIM(U, 0),        // U matrix and leading dimension
        (double*)PyArray_DATA(sigma),                                  // singular values output
        (double*)PyArray_DATA(bnd),                                    // error bounds output
        (PROPACK_CPLX_TYPE*)PyArray_DATA(V), PyArray_DIM(V, 0),        // V matrix and leading dimension
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

    current_propack_callback = NULL;
    Py_XDECREF(py_aprod);
    Py_XDECREF(args_tuple);
    return Py_BuildValue("i", propack_info);
}


static PyMethodDef propack_methods[] = {
    {"slansvd", propack_slansvd, METH_VARARGS, "Single precision SVD"},
    {"dlansvd", propack_dlansvd, METH_VARARGS, "Double precision SVD"},
    {"clansvd", propack_clansvd, METH_VARARGS, "Single precision complex SVD"},
    {"zlansvd", propack_zlansvd, METH_VARARGS, "Double precision complex SVD"},
    {"slansvd_irl", propack_slansvd_irl, METH_VARARGS, "Single precision implicitly restarted Lanczos SVD"},
    {"dlansvd_irl", propack_dlansvd_irl, METH_VARARGS, "Double precision implicitly restarted Lanczos SVD"},
    {"clansvd_irl", propack_clansvd_irl, METH_VARARGS, "Single precision complex implicitly restarted Lanczos SVD"},
    {"zlansvd_irl", propack_zlansvd_irl, METH_VARARGS, "Double precision complex implicitly restarted Lanczos SVD"},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};


static int propack_exec(PyObject* module) {
    // Initialize NumPy C API; import_array()
    if (_import_array() < 0) { return -1; }

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
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
    {0, NULL}
};


// Module definition structure using designated initializers
static struct PyModuleDef propack_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_propack",
    .m_doc = "PROPACK SVD routines",
    .m_size = 0,
    .m_methods = propack_methods,
    .m_slots = propack_slots,
};


PyMODINIT_FUNC PyInit__propack(void) {
    return PyModuleDef_Init(&propack_module);
}
