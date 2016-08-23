/*
 * Test code for ccallback.h
 *
 * This also is an internal "best-practices" code example on how to write
 * low-level callback code.
 *
 * The *trampoline_simple* function shows how to write the wrapper callback
 * function that dispatches to user-provided code (written in Python, C, Cython etc.).
 *
 * The *call_simple* function shows how to setup and teardown the ccallback_t
 * data structure.
 *
 * The *call_nodata* and *trampoline_nodata* show what to do when you need a
 * callback function for some code where it's not possible to pass along a
 * custom data pointer.
 *
 * The *call_nonlocal* and *trampoline_nonlocal* show how to use setjmp/longjmp
 * to obtain a nonlocal return on error conditions, in cases where there's no
 * mechanism to interrupt computation. Note that this is the last-resort option,
 * and only safe if there is no memory allocation between setjmp/longjmp (or you
 * need to add additonal cleanup yourself).
 *
 */

#include <setjmp.h>
#include <Python.h>

#include "ccallback.h"


#define ERROR_VALUE 2


/*
 * Trampolines for the different cases considered
 */

static double trampoline_simple(double a, int *error_flag, void *data)
{
    ccallback_t *callback = (ccallback_t*)data;
    double result = 0;
    int error = 0;

    if (callback->py_function) {
        PyGILState_STATE state = PyGILState_Ensure();
        PyObject *res, *res2;

        if (callback->user_data == NULL) {
            res = PyObject_CallFunction(callback->py_function, "d", a);
        }
        else {
            res = PyObject_CallFunction(callback->py_function, "dO", a,
                                        callback->user_data);
        }
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
        result = ((double(*)(double, int *, void *))callback->c_function)(
            a, &error, callback->user_data);
    }

    if (error) {
        *error_flag = 1;
    }

    return result;
}


static double trampoline_nodata(double a, int *error_flag)
{
    ccallback_t *callback = ccallback_obtain();
    return trampoline_simple(a, error_flag, (void *)callback);
}


static double trampoline_nonlocal(double a)
{
    ccallback_t *callback = ccallback_obtain();
    double result;
    int error_flag = 0;

    result = trampoline_simple(a, &error_flag, (void *)callback);

    if (error_flag) {
        longjmp(callback->error_buf, 1);
    }

    return result;
}


/*
 * Caller entry point functions
 */

static PyObject *call_simple(PyObject *obj, PyObject *args)
{
    PyObject *callback_obj;
    double value, result;
    ccallback_t callback;
    int error_flag = 0;
    int ret;

    if (!PyArg_ParseTuple(args, "Od", &callback_obj, &value)) {
        return NULL;
    }

    ret = ccallback_prepare(&callback, "double (double, int *, void *)", callback_obj,
                            CCALLBACK_DEFAULTS);
    if (ret != 0) {
        return NULL;
    }

    /* Call */
    result = trampoline_simple(value, &error_flag, (void *)&callback);

    ccallback_release(&callback);

    if (error_flag) {
        return NULL;
    }
    else {
        return PyFloat_FromDouble(result);
    }
}


static PyObject *call_nodata(PyObject *obj, PyObject *args)
{
    PyObject *callback_obj;
    double value, result;
    ccallback_t callback;
    int ret;
    int error_flag = 0;

    if (!PyArg_ParseTuple(args, "Od", &callback_obj, &value)) {
        return NULL;
    }

    ret = ccallback_prepare(&callback, "double (double, int *, void *)", callback_obj,
                            CCALLBACK_OBTAIN);
    if (ret != 0) {
        return NULL;
    }

    /* Call */
    result = trampoline_nodata(value, &error_flag);

    ccallback_release(&callback);

    if (error_flag) {
        return NULL;
    }
    else {
        return PyFloat_FromDouble(result);
    }
}


static PyObject *call_nonlocal(PyObject *obj, PyObject *args)
{
    PyObject *callback_obj;
    double value, result;
    ccallback_t callback;
    int ret;

    if (!PyArg_ParseTuple(args, "Od", &callback_obj, &value)) {
        return NULL;
    }

    ret = ccallback_prepare(&callback, "double (double, int *, void *)", callback_obj,
                            CCALLBACK_OBTAIN);
    if (ret != 0) {
        /* Immediate error return */
        return NULL;
    }

    /* Nonlocal return */
    if (setjmp(callback.error_buf) != 0) {
        /* Nonlocal error return */
        ccallback_release(&callback);
        return NULL;
    }

    /* Call */
    result = trampoline_nonlocal(value);

    ccallback_release(&callback);

    return PyFloat_FromDouble(result);
}


/*
 * Test functions
 */

static char *plus1_signature = "double (double, int *, void *)";


static double plus1_callback(double a, int *error_flag, void *user_data)
{
    if (a == ERROR_VALUE) {
        *error_flag = 1;
        PyErr_SetString(PyExc_ValueError, "ERROR_VALUE encountered!");
        return 0;
    }

    if (user_data == NULL) {
        return a + 1;
    }
    else {
        return a + PyFloat_AsDouble((PyObject *)user_data);
    }
}


static PyObject *get_plus1_capsule(PyObject *obj, PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyCapsule_New((void *)plus1_callback, plus1_signature, NULL);
}


/*
 * Initialize the module
 */


static PyMethodDef test_ccallback_methods[] = {
    {"call_simple", (PyCFunction)call_simple, METH_VARARGS, ""},
    {"call_nodata", (PyCFunction)call_nodata, METH_VARARGS, ""},
    {"call_nonlocal", (PyCFunction)call_nonlocal, METH_VARARGS, ""},
    {"get_plus1_capsule", (PyCFunction)get_plus1_capsule, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};


#if PY_VERSION_HEX >= 0x03000000

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
    return PyModule_Create(&test_ccallback_module);
}


#else

PyMODINIT_FUNC
init_test_ccallback(void)
{
    Py_InitModule("_test_ccallback", test_ccallback_methods);
}

#endif
