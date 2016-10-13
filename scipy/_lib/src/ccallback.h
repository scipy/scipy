/*
 * ccallback
 *
 * Callback function interface, supporting
 *
 * (1) pure Python functions
 * (2) plain C functions wrapped in PyCapsules (with Cython CAPI style signatures)
 * (3) ctypes function pointers
 * (4) cffi function pointers
 *
 * This is done avoiding magic or code generation, so you need to write some
 * boilerplate code manually.
 */

/* Boilerplate code that you need to write:

   ----------------------------------------------------------------------------
   static char *my_signatures[] = {
       "double (double, double)",
       "double (double, double, void *)",
       NULL
   };

   double my_thunk(double a, double b) {
       ccallback_t *callback = ccallback_obtain();
       double result;
       int error = 0;

       if (callback->py_function) {
           // You need to deal with GIL management yourself -- if you released
           // it, you need to also reobtain it.
           PyGILState_STATE state = PyGILState_Ensure();
           PyObject *res, *res2;

           res = PyObject_CallFunction(callback->py_function, "dd", a, b);

           if (res == NULL) {
               error = 1;
           }
           else {
               result = PyFloat_AsDouble(res);
               if (PyErr_Occurred()) {
                   error = 1;
               }
           }

           PyGILState_Release(state);
       }
       else {
           switch (callback->signature_idx) {
           case 0:
               result = ((double(*)(double, double))callback->c_function)(a, b);
               break;
           case 1:
               result = ((double(*)(double, double, void *))callback->c_function)(a, b,
                       callback->user_data);
               break;
           }
       }

       if (error) {
           // Bail out via longjmp. Note that this is not always safe. If your
           // library supports a different way of bailing out, use that instead.
           longjmp(callback.error_buf, 1);
       }

       return result;
   }

   void my_entry_point(PyObject *callback_obj) {
       ccallback_t callback;
       int callback_res;
       ...

       callback_res = ccallback_prepare(&callback, my_signatures,
               callback_obj, CCALLBACK_DEFAULTS | CCALLBACK_OBTAIN);
       if (callback_res != 0) {
           // Callback preparation failed, Python error is already set, so bail out
           return NULL;
       }

       if (setjmp(callback.error_buf) != 0) {
           // Python error during callback --- we arrive here from longjmp.
           // Only needed if you use the longjmp nonlocal return method.
           // Be sure to be aware of the many caveats associated with setjmp.
           callback_release(&callback);
           return NULL;
       }

       ...

       call_some_library_that_needs_callback(..., &my_thunk);

       ...

       callback_release(&callback);

       return;
   }
   ----------------------------------------------------------------------------

   If your library cannot pass along the extra data pointer:

   ----------------------------------------------------------------------------
   double my_thunk(double a, double b) {
       ccallback_t *callback = ccallback_obtain();
       ...
   }
   ----------------------------------------------------------------------------

   In addition, add CCALLBACK_OBTAIN to flags passed to
   ccallback_prepare. This will have a performance impact.

   See _test_ccallback.c for full examples.
 */

#ifndef CCALLBACK_H_
#define CCALLBACK_H_


#include <Python.h>
#include <setjmp.h>

#define CCALLBACK_DEFAULTS 0x0
#define CCALLBACK_OBTAIN   0x1
#define CCALLBACK_PARSE    0x2

typedef struct ccallback ccallback_t;

struct ccallback {
    void *c_function;
    PyObject *py_function;
    void *user_data;
    jmp_buf error_buf;
    ccallback_t *prev_callback;
    int signature_index;

    /* Unused variables that can be used by the thunk etc. code for any purpose */
    long info;
    void *info_p;
};


/*
 * Thread-local storage
 */

#if defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ >= 4)))

static __thread ccallback_t *_active_ccallback = NULL;

static void *ccallback__get_thread_local(void)
{
    return (void *)_active_ccallback;
}

static int ccallback__set_thread_local(void *value)
{
    _active_ccallback = value;
}

static ccallback_t *ccallback_obtain(void)
{
    return (ccallback_t *)ccallback__get_thread_local();
}

#else

static void *ccallback__get_thread_local(void)
{
    PyObject *local_dict, *capsule;
    void *callback_ptr;

    local_dict = PyThreadState_GetDict();
    if (local_dict == NULL) {
        Py_FatalError("scipy/ccallback: failed to get local thread state");
    }

    capsule = PyDict_GetItemString(local_dict, "__scipy_ccallback");
    if (capsule == NULL) {
        return NULL;
    }

    callback_ptr = PyCapsule_GetPointer(capsule, NULL);
    if (callback_ptr == NULL) {
        Py_FatalError("scipy/ccallback: invalid callback state");
    }

    return callback_ptr;
}

static int ccallback__set_thread_local(void *value)
{
    PyObject *local_dict;

    local_dict = PyThreadState_GetDict();
    if (local_dict == NULL) {
        Py_FatalError("scipy/ccallback: failed to get local thread state");
    }

    if (value == NULL) {
        return PyDict_DelItemString(local_dict, "__scipy_ccallback");
    }
    else {
        PyObject *capsule;
        int ret;

        capsule = PyCapsule_New(value, NULL, NULL);
        if (capsule == NULL) {
            return -1;
        }
        ret = PyDict_SetItemString(local_dict, "__scipy_ccallback", capsule);
        Py_DECREF(capsule);
        return ret;
    }
}

static ccallback_t *ccallback_obtain(void)
{
    PyGILState_STATE state;
    ccallback_t *callback_ptr;

    state = PyGILState_Ensure();

    callback_ptr = (ccallback_t *)ccallback__get_thread_local();
    if (callback_ptr == NULL) {
        Py_FatalError("scipy/ccallback: failed to get thread local state");
    }

    PyGILState_Release(state);

    return callback_ptr;
}

#endif


/* Set up callback. */
static int ccallback_prepare(ccallback_t *callback, char **signatures, PyObject *callback_obj, int flags)
{
    static PyTypeObject *lowlevelcallable_type = NULL;
    PyObject *callback_obj2 = NULL;
    PyObject *capsule = NULL;

    if (lowlevelcallable_type == NULL) {
        PyObject *module;

        module = PyImport_ImportModule("scipy._lib._ccallback");
        if (module == NULL) {
            goto error;
        }

        lowlevelcallable_type = (PyTypeObject *)PyObject_GetAttrString(module, "LowLevelCallable");
        if (lowlevelcallable_type == NULL) {
            goto error;
        }
    }

    if ((flags & CCALLBACK_PARSE) && !PyObject_TypeCheck(callback_obj, lowlevelcallable_type)) {
        /* Parse callback */
        callback_obj2 = PyObject_CallMethod(lowlevelcallable_type, "_parse_callback",
                                            "O", callback_obj);
        if (callback_obj2 == NULL) {
            goto error;
        }

        callback_obj = callback_obj2;

        if (PyCapsule_CheckExact(callback_obj)) {
            capsule = callback_obj;
        }
    }

    if (PyCallable_Check(callback_obj)) {
        /* Python callable */
        callback->py_function = callback_obj;
        callback->c_function = NULL;
        callback->user_data = NULL;
        callback->signature_index = -1;
    }
    else if (PyObject_TypeCheck(callback_obj, lowlevelcallable_type) &&
             PyTuple_Check(callback_obj) &&
             PyCallable_Check(PyTuple_GET_ITEM(callback_obj, 0))) {
        /* Python callable in LowLevelCallable */
        callback->py_function = PyTuple_GET_ITEM(callback_obj, 0);
        callback->c_function = NULL;
        callback->user_data = NULL;
        callback->signature_index = -1;
    }
    else if (capsule != NULL ||
             PyObject_TypeCheck(callback_obj, lowlevelcallable_type) &&
             PyTuple_Check(callback_obj) &&
             PyCapsule_CheckExact(PyTuple_GET_ITEM(callback_obj, 0))) {
        /* PyCapsule in LowLevelCallable (or parse result from above) */
        void *ptr, *user_data;
        char **sig;
        const char *name;

        if (capsule == NULL) {
            capsule = PyTuple_GET_ITEM(callback_obj, 0);
        }

        name = PyCapsule_GetName(capsule);
        if (PyErr_Occurred()) {
            goto error;
        }
        
        callback->signature_index = 0;
        for (sig = signatures; *sig != NULL; ++sig, ++callback->signature_index) {
            if (name && strcmp(name, *sig) == 0) {
                break;
            }
        }

        ptr = PyCapsule_GetPointer(capsule, *sig);
        if (ptr == NULL) {
            PyErr_SetString(PyExc_ValueError, "Invalid function signature in PyCapsule");
            goto error;
        }

        user_data = PyCapsule_GetContext(capsule);
        if (PyErr_Occurred()) {
            goto error;
        }

        callback->py_function = NULL;
        callback->c_function = ptr;
        callback->user_data = user_data;
    }
    else {
        PyErr_SetString(PyExc_ValueError, "invalid callable given");
        goto error;
    }

    if (flags & CCALLBACK_OBTAIN) {
        callback->prev_callback = ccallback__get_thread_local();
        ccallback__set_thread_local((void *)callback);
    }
    else {
        callback->prev_callback = NULL;
    }

    Py_XDECREF(callback_obj2);
    return 0;

error:
    Py_XDECREF(callback_obj2);
    return -1;
}


/* Tear down callback. */
static void ccallback_release(ccallback_t *callback)
{
    if (callback->prev_callback != NULL) {
        ccallback__set_thread_local(callback->prev_callback);
    }
    callback->prev_callback = NULL;
    callback->c_function = NULL;
    callback->py_function = NULL;
}

#endif /* CCALLBACK_H_ */
