/*
 * Test code for ccallback.h
 *
 * This also is an internal "best-practices" code example on how to write
 * low-level callback code. (In the examples below, it is assumed the semantics
 * and signatures of the callbacks in test_call_* are fixed by some 3rd party
 * library e.g., implemented in FORTRAN, and they are not necessarily the optimal
 * way.)
 *
 * The general structure of callbacks is the following:
 *
 * - entry point function, callable from Python, calls 3rd party library code
 *   (test_call_*)
 * - 3rd party library code, calls the callback (included as a part of test_call_*)
 * - callback thunk, has the signature expected by the 3rd party library, and
 *   translates the callback call to a Python function or user-provided
 *   low-level function call (thunk_*).
 *
 * The *thunk_simple* function shows how to write a callback thunk that
 * dispatches to user-provided code (written in Python, C, Cython etc.).
 *
 * The *call_simple* function shows how to setup and teardown the ccallback_t
 * data structure in the entry point.
 *
 * The *call_nodata* and *thunk_nodata* show what to do when you need a
 * callback function for some code where it's not possible to pass along a
 * custom data pointer.
 *
 * The *call_nonlocal* and *thunk_nonlocal* show how to use setjmp/longjmp
 * to obtain a nonlocal return on error conditions, in cases where there's no
 * mechanism to interrupt computation. Note that this is the last-resort option,
 * and only safe if there is no memory allocation between setjmp/longjmp (or you
 * need to add additional cleanup yourself).
 *
 */

#include <Python.h>
#include <setjmp.h>

#include "ccallback.h"


#define ERROR_VALUE 2


/*
 * Example 3rd party library code, to be interfaced with
 */

static double library_call_simple(double value, int *error_flag, double (*callback)(double, int*, void*),
                                  void *data)
{
    *error_flag = 0;
    return callback(value, error_flag, data);
}


static double library_call_nodata(double value, int *error_flag, double (*callback)(double, int*))
{
    *error_flag = 0;
    return callback(value, error_flag);
}


static double library_call_nonlocal(double value, double (*callback)(double))
{
    return callback(value);
}



/*
 * Callback thunks for the different cases
 */

static double test_thunk_simple(double a, int *error_flag, void *data)
{
    ccallback_t *callback = (ccallback_t*)data;
    double result = 0;
    int error = 0;

    if (callback->py_function) {
        PyGILState_STATE state = PyGILState_Ensure();
        PyObject *res, *res2;

        res = PyObject_CallFunction(callback->py_function, "d", a);

        if (res == NULL) {
            error = 1;
        }
        else {
            res2 = PyNumber_Float(res);
            if (res2 == NULL) {
                error = 1;
            }
            else {
                result = PyFloat_AsDouble(res2);
                if (PyErr_Occurred()) {
                    error = 1;
                }
                Py_DECREF(res2);
            }
            Py_DECREF(res);
        }

        PyGILState_Release(state);
    }
    else {
        if (callback->signature->value == 0) {
            result = ((double(*)(double, int *, void *))callback->c_function)(
                a, &error, callback->user_data);
        }
        else {
            result = ((double(*)(double, double, int *, void *))callback->c_function)(
                a, 0.0, &error, callback->user_data);
        }
    }

    if (error) {
        *error_flag = 1;
    }

    return result;
}


static double test_thunk_nodata(double a, int *error_flag)
{
    ccallback_t *callback = ccallback_obtain();
    return test_thunk_simple(a, error_flag, (void *)callback);
}


static double test_thunk_nonlocal(double a)
{
    ccallback_t *callback = ccallback_obtain();
    double result;
    int error_flag = 0;

    result = test_thunk_simple(a, &error_flag, (void *)callback);

    if (error_flag) {
        longjmp(callback->error_buf, 1);
    }

    return result;
}


/*
 * Caller entry point functions
 */

static ccallback_signature_t signatures[] = {
    {"double (double, int *, void *)", 0},
    {"double (double, double, int *, void *)", 1},
#if NPY_SIZEOF_INT == NPY_SIZEOF_SHORT
    {"double (double, short *, void *)", 0},
    {"double (double, double, short *, void *)", 1},
#endif
#if NPY_SIZEOF_INT == NPY_SIZEOF_LONG
    {"double (double, long *, void *)", 0},
    {"double (double, double, long *, void *)", 1},
#endif
    {NULL}
};

static PyObject *test_call_simple(PyObject *obj, PyObject *args)
{
    PyObject *callback_obj;
    double value, result;
    ccallback_t callback;
    int error_flag;
    int ret;

    if (!PyArg_ParseTuple(args, "Od", &callback_obj, &value)) {
        return NULL;
    }

    ret = ccallback_prepare(&callback, signatures, callback_obj, CCALLBACK_DEFAULTS);
    if (ret != 0) {
        return NULL;
    }

    /* Call 3rd party library code */
    Py_BEGIN_ALLOW_THREADS
    result = library_call_simple(value, &error_flag, test_thunk_simple, (void *)&callback);
    Py_END_ALLOW_THREADS

    ccallback_release(&callback);

    if (error_flag) {
        return NULL;
    }
    else {
        return PyFloat_FromDouble(result);
    }
}


static PyObject *test_call_nodata(PyObject *obj, PyObject *args)
{
    PyObject *callback_obj;
    double value, result;
    ccallback_t callback;
    int ret;
    int error_flag;

    if (!PyArg_ParseTuple(args, "Od", &callback_obj, &value)) {
        return NULL;
    }

    ret = ccallback_prepare(&callback, signatures, callback_obj, CCALLBACK_OBTAIN);
    if (ret != 0) {
        return NULL;
    }

    /* Call 3rd party library code */
    Py_BEGIN_ALLOW_THREADS
    result = library_call_nodata(value, &error_flag, test_thunk_nodata);
    Py_END_ALLOW_THREADS

    ccallback_release(&callback);

    if (error_flag) {
        return NULL;
    }
    else {
        return PyFloat_FromDouble(result);
    }
}


static PyObject *test_call_nonlocal(PyObject *obj, PyObject *args)
{
    PyObject *callback_obj;
    double value, result;
    int ret;
    ccallback_t callback;
    PyThreadState *_save = NULL;

    if (!PyArg_ParseTuple(args, "Od", &callback_obj, &value)) {
        return NULL;
    }

    ret = ccallback_prepare(&callback, signatures, callback_obj, CCALLBACK_OBTAIN);
    if (ret != 0) {
        /* Immediate error return */
        return NULL;
    }

    /* Nonlocal return */
    _save = PyEval_SaveThread();

    if (setjmp(callback.error_buf) != 0) {
        /* Nonlocal error return */
        PyEval_RestoreThread(_save);
        ccallback_release(&callback);
        return NULL;
    }

    /* Call 3rd party library code */
    result = library_call_nonlocal(value, test_thunk_nonlocal);

    PyEval_RestoreThread(_save);

    ccallback_release(&callback);

    return PyFloat_FromDouble(result);
}


/*
 * Functions for testing the PyCapsule interface
 */

static char *test_plus1_signature = "double (double, int *, void *)";


static double test_plus1_callback(double a, int *error_flag, void *user_data)
{
    if (a == ERROR_VALUE) {
        PyGILState_STATE state = PyGILState_Ensure();
        *error_flag = 1;
        PyErr_SetString(PyExc_ValueError, "ERROR_VALUE encountered!");
        PyGILState_Release(state);
        return 0;
    }

    if (user_data == NULL) {
        return a + 1;
    }
    else {
        return a + *(double *)user_data;
    }
}


static PyObject *test_get_plus1_capsule(PyObject *obj, PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyCapsule_New((void *)test_plus1_callback, test_plus1_signature, NULL);
}


static char *test_plus1b_signature = "double (double, double, int *, void *)";


static double test_plus1b_callback(double a, double b, int *error_flag, void *user_data)
{
    return test_plus1_callback(a, error_flag, user_data) + b;
}


static PyObject *test_get_plus1b_capsule(PyObject *obj, PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyCapsule_New((void *)test_plus1b_callback, test_plus1b_signature, NULL);
}


static char *test_plus1bc_signature = "double (double, double, double, int *, void *)";


static double test_plus1bc_callback(double a, double b, double c, int *error_flag, void *user_data)
{
    return test_plus1_callback(a, error_flag, user_data) + b + c;
}


static PyObject *test_get_plus1bc_capsule(PyObject *obj, PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyCapsule_New((void *)test_plus1bc_callback, test_plus1bc_signature, NULL);
}


static void data_capsule_destructor(PyObject *capsule)
{
    void *data;
    data = PyCapsule_GetPointer(capsule, NULL);
    free(data);
}

static PyObject *test_get_data_capsule(PyObject *obj, PyObject *args)
{
    double *data;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    data = (double *)malloc(sizeof(double));
    if (data == NULL) {
        return PyErr_NoMemory();
    }

    *data = 2.0;

    return PyCapsule_New((void *)data, NULL, data_capsule_destructor);
}


/*
 * Initialize the module
 */

static PyMethodDef test_ccallback_methods[] = {
    {"test_call_simple", (PyCFunction)test_call_simple, METH_VARARGS, ""},
    {"test_call_nodata", (PyCFunction)test_call_nodata, METH_VARARGS, ""},
    {"test_call_nonlocal", (PyCFunction)test_call_nonlocal, METH_VARARGS, ""},
    {"test_get_plus1_capsule", (PyCFunction)test_get_plus1_capsule, METH_VARARGS, ""},
    {"test_get_plus1b_capsule", (PyCFunction)test_get_plus1b_capsule, METH_VARARGS, ""},
    {"test_get_plus1bc_capsule", (PyCFunction)test_get_plus1bc_capsule, METH_VARARGS, ""},
    {"test_get_data_capsule", (PyCFunction)test_get_data_capsule, METH_VARARGS, ""},
    {"test_get_data_capsule", (PyCFunction)test_get_data_capsule, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef test_ccallback_module = {
    PyModuleDef_HEAD_INIT,
    "_test_ccallback",
    NULL,
    -1,
    test_ccallback_methods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC
PyInit__test_ccallback(void)
{
    PyObject *module;

    module = PyModule_Create(&test_ccallback_module);
    if (module == NULL) {
        return module;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;

}
