/*
 * ccallback
 *
 * Callback function interface, supporting
 *
 * (1) pure Python functions
 * (2) Cython functions
 * (3) plain C functions wrapped in PyCapsules
 *
 * This is done avoiding magic or code generation, so you need to write some
 * boilerplate code manually.
 */

/* Boilerplate code that you need to write:

   double my_trampoline(double a, double b, void *data) {
       ccallback_t *callback = (ccallback_t*)data;
       double result;
       int error = 0;

       if (callback->py_function) {
           // You need to deal with GIL management yourself -- if you released
           // it, you need to also reobtain it.
           PyGILState_STATE state = PyGILState_Ensure();
           PyObject *res, *res2;

           if (callback->user_data == NULL) {
               res = PyObject_CallFunction(callback->py_function, "dd", a, b);
           }
           else {
               res = PyObject_CallFunction(callback->py_function, "ddO", a, b, callback->user_data);
           }
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
           result = ((double(*)(double, double, void *))callback->c_function)(a, b,
                   callback->user_data);
       }

       if (error) {
           // Bail out via longjmp. Note that this is not always safe. If your
           // libary supports a different way of bailing out, use that instead.
           longjmp(callback.error_buf, 1);
       }

       return result;
   }

   void my_entry_point(PyObject *callback_obj) {
       ccallback_t callback;
       int callback_res;
       ...

       callback_res = ccallback_prepare("double (double *, double *, void *)",
               callback_obj, CCALLBACK_DEFAULTS);
       if (callback_res != 0) {
           // Callback preparation failed, Python error is already set, so bail out
           return NULL;
       }

       if (setjmp(callback.error_buf) != 0) {
           // Python error during callback --- we arrive here from longjmp.
           // Only needed if you use the longjmp nonlocal return method.
           callback_release(&callback);
           return NULL;
       }

       ...

       call_some_library_that_needs_callback(..., &my_trampoline, (void*)&callback);

       ...

       callback_release(&callback);

       return;
   }

   If your library cannot pass along the extra data pointer:

   double my_trampoline(double a, double b) {
       ccallback_t *callback = ccallback_obtain();
       ...
   }

   In addition, add CCALLBACK_NEED_OBTAIN to flags passed to
   ccallback_prepare. This will have a performance impact.

   See _test_ccallback.c for full examples.
 */

#ifndef CCALLBACK_H_
#define CCALLBACK_H_


#include <Python.h>


#define CCALLBACK_DEFAULTS 0x0
#define CCALLBACK_OBTAIN   0x1


typedef struct ccallback ccallback_t;

struct ccallback {
    void *c_function;
    PyObject *py_function;
    void *user_data;
    jmp_buf error_buf;
    ccallback_t *prev_callback;
};


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


void *ccallback__parse_data_pointer(PyObject *data_obj, int *success)
{
    void *data_ptr;

    if (PyCapsule_CheckExact(data_obj)) {
        data_ptr = PyCapsule_GetPointer(data_obj, NULL);
        if (data_ptr == NULL) {
            *success = 0;
            return NULL;
        }
        *success = 1;
        return data_ptr;
    }

    *success = 1;
    return (void*)data_obj;
}


PyObject *ccallback__ctypes_get_type(PyObject *ctypes_module, char *name, Py_ssize_t name_length,
                                     int pointer_level)
{
    char name_buffer[256];
    PyObject *ctypes_type = NULL, *ctypes_ptr = NULL;

    /* Void pointer has a different name */
    if (pointer_level > 0 && strncmp(name, "void", name_length) == 0) {
        name = "void_p";
        name_length = 6;
        --pointer_level;
    }

    /* Format ctypes type name */
    name_buffer[0] = 'c';
    name_buffer[1] = '_';
    if (name_length > sizeof(name_buffer) - 3) {
        name_length = sizeof(name_buffer) - 3;
    }
    strncpy(name_buffer + 2, name, name_length);
    name_buffer[name_length + 2] = '\0';

    /* Get type */
    ctypes_type = PyObject_GetAttrString(ctypes_module, name_buffer);
    if (ctypes_type == NULL) {
        return NULL;
    }

    while (pointer_level > 0) {
        ctypes_ptr = PyObject_CallMethod(ctypes_module, "POINTER", "O", ctypes_type);
        Py_DECREF(ctypes_type);
        if (ctypes_ptr == NULL) {
            return NULL;
        }
        ctypes_type = ctypes_ptr;
        --pointer_level;
    }

    return ctypes_type;
}


void ccallback__parse_carg(char *p, char **end, Py_ssize_t *name_length, int *pointer_level)
{
    char *p2, *p3;

    p2 = p;
    while (*p2 != '\0' && *p2 != '(' && *p2 != ')' && *p2 != ',') {
        ++p2;
    }

    *pointer_level = 0;

    if (p2 > p) {
        p3 = p2 - 1;
        while (*p3 == '*') {
            /* Pointer argument */
            ++*pointer_level;
            --p3;
            while (*p3 == ' ') {
                --p3;
            }
        }
        while (*p3 == ' ') {
            --p3;
        }
        *name_length = p3 - p + 1;
    }
    else {
        name_length = 0;
    }

    *end = p2;
}


typedef struct {
    PyObject_HEAD
    char *b_ptr;
} ccallback__cfuncptr_object;


static void *ccallback__get_ctypes_function_pointer(PyObject * obj)
{
    return (*((void **) (((ccallback__cfuncptr_object *) (obj))->b_ptr)));
}


static int ccallback__parse_ctypes(ccallback_t *callback, PyObject *callback_obj, char *signature)
{
    PyObject *ctypes_module = NULL, *cfuncptr_type = NULL, *expected_type = NULL,
        *argtypes = NULL, *restype = NULL, *arg = NULL;
    PyObject *cfunc;
    void *user_data;
    int return_value = 0;

    ctypes_module = PyImport_ImportModule("ctypes");
    if (ctypes_module == NULL) {
        /* No ctypes, skip detection */
        PyErr_Clear();
        return_value = 0;
        goto done;
    }

    cfuncptr_type = PyObject_GetAttrString(ctypes_module, "_CFuncPtr");
    if (cfuncptr_type == NULL) {
        return_value = -1;
        goto done;
    }

    if (PyObject_TypeCheck(callback_obj, (PyTypeObject *)cfuncptr_type)) {
        cfunc = callback_obj;
        user_data = NULL;
    }
    else if (PyTuple_CheckExact(callback_obj) && PyTuple_GET_SIZE(callback_obj) == 2 &&
             PyObject_TypeCheck(PyTuple_GET_ITEM(callback_obj, 0), (PyTypeObject *)cfuncptr_type))
    {
        int success;
        cfunc = PyTuple_GET_ITEM(callback_obj, 0);
        user_data = ccallback__parse_data_pointer(PyTuple_GET_ITEM(callback_obj, 1), &success);
        if (!success) {
            return_value = -1;
            goto done;
        }
    }
    else {
        return_value = 0;
        goto done;
    }

    restype = PyObject_GetAttrString(cfunc, "restype");
    if (restype == NULL) {
        return_value = -1;
        goto done;
    }

    argtypes = PyObject_GetAttrString(cfunc, "argtypes");
    if (argtypes == NULL) {
        return_value = -1;
        goto done;
    }

    /* Parse and check the signature string */
    {
        char *p = signature;
        char *p2, *p3;
        Py_ssize_t num_args, cur_arg;
        Py_ssize_t name_length;
        int pointer_level;

        /* Return type */
        ccallback__parse_carg(p, &p2, &name_length, &pointer_level);
        expected_type = ccallback__ctypes_get_type(ctypes_module, p, name_length, pointer_level);
        if (expected_type == NULL) {
            return_value = -1;
            goto done;
        }

        if (restype != expected_type) {
            PyObject *expected_type_str;
            return_value = -1;
            expected_type_str = PyObject_Str(expected_type);
            if (expected_type_str == NULL) {
                goto done;
            }
            PyErr_Format(PyExc_ValueError, "ctypes function restype is not '%s'",
                         PyString_AsString(expected_type_str));
            Py_DECREF(expected_type_str);
            goto done;
        }

        Py_DECREF(expected_type);
        expected_type = NULL;

        /* Check arguments */
        p = p2;
        if (*p == '(') {
            ++p;
        }

        cur_arg = 0;
        num_args = PySequence_Size(argtypes);

        while (*p != '\0' && *p != ')') {
            if (num_args <= cur_arg) {
                PyErr_SetString(PyExc_ValueError, "ctypes function takes too few arguments");
                return_value = -1;
                goto done;
            }

            /* Obtain argument type name */
            ccallback__parse_carg(p, &p2, &name_length, &pointer_level);
            expected_type = ccallback__ctypes_get_type(ctypes_module, p, name_length, pointer_level);
            if (expected_type == NULL) {
                return_value = -1;
                goto done;
            }

            arg = PySequence_GetItem(argtypes, cur_arg);
            if (arg == NULL) {
                return_value = -1;
                goto done;
            }

            if (arg != expected_type) {
                PyObject *expected_type_str, *s2;
                return_value = -1;
                expected_type_str = PyObject_Str(expected_type);
                if (expected_type_str == NULL) {
                    goto done;
                }
                PyErr_Format(PyExc_ValueError, "ctypes function argtypes[%d] != %s",
                             cur_arg, PyString_AsString(expected_type_str));
                Py_DECREF(expected_type_str);
                goto done;
            }

            Py_DECREF(arg);
            Py_DECREF(expected_type);
            arg = NULL;
            expected_type = NULL;

            p = p2;
            if (*p == ',') {
                ++p;
                while (*p == ' ') {
                    ++p;
                }
            }

            cur_arg++;
        }
    }

    /* All checks passed */
    return_value = 1;
    callback->c_function = ccallback__get_ctypes_function_pointer(cfunc);
    callback->py_function = NULL;
    callback->user_data = user_data;

done:
    Py_XDECREF(ctypes_module);
    Py_XDECREF(cfuncptr_type);
    Py_XDECREF(expected_type);
    Py_XDECREF(argtypes);
    Py_XDECREF(restype);
    Py_XDECREF(arg);
    return return_value;
}


/*
 * Public functions
 */

/* Find the callback function active for the current thread. */
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


/* Set up callback. */
static int ccallback_prepare(ccallback_t *callback, char *signature, PyObject *callback_obj, int flags)
{
    /* Parse callback_obj */
    if (PyCapsule_CheckExact(callback_obj)) {
        /* C function in a PyCapsule + user data as context */
        void *ptr, *user_data;

        ptr = PyCapsule_GetPointer(callback_obj, signature);
        if (ptr == NULL) {
            PyErr_SetString(PyExc_ValueError, "Invalid function signature in PyCapsule");
            return -1;
        }

        user_data = PyCapsule_GetContext(callback_obj);
        if (PyErr_Occurred()) {
            return -1;
        }

        callback->py_function = NULL;
        callback->c_function = ptr;
        callback->user_data = user_data;
    }
    else if (PyTuple_CheckExact(callback_obj) && PyTuple_GET_SIZE(callback_obj) == 2 &&
             PyCapsule_CheckExact(PyTuple_GET_ITEM(callback_obj, 0))) {
        /* C function in a PyCapsule + explicit user data (possibly in capsule) */
        void *ptr, *data_ptr;
        int success;

        ptr = PyCapsule_GetPointer(PyTuple_GET_ITEM(callback_obj, 0), signature);
        if (ptr == NULL) {
            PyErr_SetString(PyExc_ValueError, "Invalid function signature in PyCapsule");
            return -1;
        }

        data_ptr = ccallback__parse_data_pointer(PyTuple_GET_ITEM(callback_obj, 1), &success);
        if (!success) {
            return -1;
        }

        callback->py_function = NULL;
        callback->c_function = ptr;
        callback->user_data = data_ptr;
    }
    else if (PyTuple_CheckExact(callback_obj) &&
             (PyTuple_GET_SIZE(callback_obj) == 2 || PyTuple_GET_SIZE(callback_obj) == 3) &&
             PyModule_CheckExact(PyTuple_GET_ITEM(callback_obj, 0)) &&
             PyString_CheckExact(PyTuple_GET_ITEM(callback_obj, 1))) {
        /* Module + function name + optional user data --- Cython callable */
        void *ptr, *data_ptr;
        PyObject *pyx_api, *module_dict, *capsule;
        int success;

        module_dict = PyModule_GetDict(PyTuple_GET_ITEM(callback_obj, 0));
        if (module_dict == NULL) {
            return -1;
        }

        pyx_api = PyDict_GetItemString(module_dict, "__pyx_capi__");
        if (pyx_api == NULL) {
            PyErr_SetString(PyExc_KeyError, "Module is not a Cython module containing __pyx_capi__");
            return -1;
        }

        capsule = PyDict_GetItem(pyx_api, PyTuple_GET_ITEM(callback_obj, 1));
        if (capsule == NULL) {
            PyErr_SetString(PyExc_KeyError, "Function not found in __pyx_capi__");
            return -1;
        }

        ptr = PyCapsule_GetPointer(capsule, signature);
        if (ptr == NULL) {
            PyErr_SetString(PyExc_ValueError, "Invalid function signature in PyCapsule");
            return -1;
        }

        if (PyTuple_GET_SIZE(callback_obj) == 3) {
            data_ptr = ccallback__parse_data_pointer(PyTuple_GET_ITEM(callback_obj, 2), &success);
            if (!success) {
                return -1;
            }
        }
        else {
            data_ptr = NULL;
        }

        callback->py_function = NULL;
        callback->c_function = ptr;
        callback->user_data = data_ptr;
    }
    else {
        int ret;

        /* Try ctypes */
        ret = ccallback__parse_ctypes(callback, callback_obj, signature);
        if (ret == -1) {
            return -1;
        }
        else if (ret == 1) {
            /* ok */
        }
        else if (PyCallable_Check(callback_obj)) {
            /* Python callable */
            callback->py_function = callback_obj;
            callback->c_function = NULL;
            callback->user_data = NULL;
        }
        else if (PyTuple_CheckExact(callback_obj) &&
                 PyTuple_GET_SIZE(callback_obj) == 2 &&
                 PyCallable_Check(PyTuple_GET_ITEM(callback_obj, 0))) {
            /* Python callable + user data */
            callback->py_function = PyTuple_GET_ITEM(callback_obj, 0);
            callback->c_function = NULL;
            callback->user_data = PyTuple_GET_ITEM(callback_obj, 1);
        }
        else {
            PyErr_SetString(PyExc_ValueError, "invalid callable given");
            return -1;
        }
    }

    if (flags & CCALLBACK_OBTAIN) {
        callback->prev_callback = ccallback__get_thread_local();
        ccallback__set_thread_local((void *)callback);
    }
    else {
        callback->prev_callback = NULL;
    }

    return 0;
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
