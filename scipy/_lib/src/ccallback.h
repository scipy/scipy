/*
 * ccallback
 *
 * Callback function interface, supporting
 *
 * (1) pure Python functions
 * (2) Cython functions
 * (3) plain C functions wrapped in PyCapsules
 * (4) ctypes function pointers
 * (5) cffi function pointers
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


#define CCALLBACK_MAX_ARGS 128


#if PY_VERSION_HEX >= 0x03000000

#define PY3K 1

#define PyUString_Type PyUnicode_Type
#define PyUString_Check PyUnicode_Check
#define PyUString_CheckExact PyUnicode_CheckExact
#define PyUStringObject PyUnicodeObject
#define PyUString_FromString PyUnicode_FromString
#define PyUString_FromStringAndSize PyUnicode_FromStringAndSize
#define PyUString_FromFormat PyUnicode_FromFormat
#define PyUString_Concat PyUnicode_Concat2
#define PyUString_ConcatAndDel PyUnicode_ConcatAndDel
#define PyUString_GET_SIZE PyUnicode_GET_SIZE
#define PyUString_Size PyUnicode_Size
#define PyUString_InternFromString PyUnicode_InternFromString
#define PyUString_Format PyUnicode_Format
#define PyUString_AsUTF8String PyUnicode_AsUTF8String

static int PyUString_isequal(PyObject *obj, char *s, Py_ssize_t length)
{
    PyObject *bytes;
    char *s2;
    int ret;

    bytes = PyUnicode_AsUTF8String(obj);
    if (bytes == NULL) {
        return -1;
    }

    s2 = PyBytes_AsString(bytes);
    if (s2 == NULL) {
        Py_DECREF(bytes);
        return -1;
    }

    ret = strncmp(s2, s, length);
    if (ret == 0) {
        ret = !(s2[length] == '\0');
    }
    Py_DECREF(bytes);

    return (ret == 0);
}

#else

#undef PY3K

#define PyUString_Type PyString_Type
#define PyUString_Check PyString_Check
#define PyUString_CheckExact PyString_CheckExact
#define PyUStringObject PyStringObject
#define PyUString_FromString PyString_FromString
#define PyUString_FromStringAndSize PyString_FromStringAndSize
#define PyUString_FromFormat PyString_FromFormat
#define PyUString_Concat PyString_Concat
#define PyUString_ConcatAndDel PyString_ConcatAndDel
#define PyUString_GET_SIZE PyString_GET_SIZE
#define PyUString_Size PyString_Size
#define PyUString_InternFromString PyString_InternFromString
#define PyUString_Format PyString_Format

static int PyUString_isequal(PyObject *obj, char *s, Py_ssize_t length)
{
    char *s2;
    int ret;

    s2 = PyBytes_AsString(obj);
    if (s2 == NULL) {
        return -1;
    }

    ret = strncmp(s2, s, length);
    if (ret == 0) {
        ret = !(s2[length] == '\0');
    }
    return (ret == 0);
}

#endif


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


/*
 * Internal: thread-local storage
 */

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


/*
 * Internal: unpack user_data from PyCapsule
 */
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


/*
 * Internal: parse function signatures
 */

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


int ccallback__parse_signature(char *signature, char **types, Py_ssize_t *lengths, int *pointer_levels,
                               int max_items)
{
    char *p, *p2;
    int cur_item;

    if (max_items <= 0) {
        return -1;
    }

    p = signature;

    /* Parse return value */
    types[0] = p;
    ccallback__parse_carg(p, &p2, lengths, pointer_levels);

    /* Parse arguments */
    p = p2;
    if (*p == '(') {
        ++p;
    }

    cur_item = 1;
    while (*p != '\0' && *p != ')') {
        if (max_items <= cur_item) {
            return -1;
        }

        /* Parse one argument */
        types[cur_item] = p;
        ccallback__parse_carg(p, &p2, lengths + cur_item, pointer_levels + cur_item);

        /* Next argument */
        p = p2;
        if (*p == ',') {
            ++p;
            while (*p == ' ') {
                ++p;
            }
        }

        cur_item++;
    }

    return cur_item - 1;
}


/*
 * Internal: obtain function pointers from ctypes
 */

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
    if (name_length > (Py_ssize_t)sizeof(name_buffer) - 3) {
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
    static PyObject *ctypes_module = NULL, *cfuncptr_type = NULL;
    static int import_failed = 0;

    PyObject *expected_type = NULL,
        *argtypes = NULL, *restype = NULL, *arg = NULL;
    PyObject *cfunc;
    void *user_data;
    int return_value = 0;

    if (import_failed) {
        return 0;
    }

    if (ctypes_module == NULL) {
        ctypes_module = PyImport_ImportModule("ctypes");
        if (ctypes_module == NULL) {
            /* No ctypes, skip detection */
            PyErr_Clear();
            import_failed = 1;
            return_value = 0;
            goto done;
        }
    }

    if (cfuncptr_type == NULL) {
        cfuncptr_type = PyObject_GetAttrString(ctypes_module, "_CFuncPtr");
        if (cfuncptr_type == NULL) {
            goto error;
        }
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
            goto error;
        }
    }
    else {
        return_value = 0;
        goto done;
    }

    restype = PyObject_GetAttrString(cfunc, "restype");
    if (restype == NULL) {
        goto error;
    }

    argtypes = PyObject_GetAttrString(cfunc, "argtypes");
    if (argtypes == NULL) {
        goto error;
    }

    /* Parse and check the signature string */
    {
        int num_args, cur_arg;
        Py_ssize_t size;

        char *types[CCALLBACK_MAX_ARGS];
        int pointer_levels[CCALLBACK_MAX_ARGS];
        Py_ssize_t lengths[CCALLBACK_MAX_ARGS];

        num_args = ccallback__parse_signature(signature, types, lengths, pointer_levels, CCALLBACK_MAX_ARGS);
        if (num_args < 0) {
            PyErr_Format(PyExc_RuntimeError, "scipy/ccallback: failed to parse signature '%s'",
                         signature);
            goto error;
        }

        /* Return type */
        expected_type = ccallback__ctypes_get_type(ctypes_module, types[0], lengths[0], pointer_levels[0]);
        if (expected_type == NULL) {
            goto error;
        }

        if (restype != expected_type) {
#if PY3K
            PyErr_Format(PyExc_ValueError, "ctypes function restype is not '%S'",
                         expected_type);
#else
            PyObject *expected_type_str;
            expected_type_str = PyObject_Str(expected_type);
            if (expected_type_str == NULL) {
                goto error;
            }
            PyErr_Format(PyExc_ValueError, "ctypes function restype is not '%s'",
                         PyBytes_AsString(expected_type_str));
            Py_DECREF(expected_type_str);
#endif
            goto error;
        }

        Py_DECREF(expected_type);
        expected_type = NULL;

        /* Check arguments */
        size = PySequence_Size(argtypes);
        if (size == -1) {
            goto error;
        }
        else if (size != num_args) {
            PyErr_SetString(PyExc_ValueError, "ctypes function takes wrong number of arguments");
            goto error;
        }

        for (cur_arg = 0; cur_arg < num_args; ++cur_arg) {
            /* Obtain argument type name */
            expected_type = ccallback__ctypes_get_type(ctypes_module, types[cur_arg+1], lengths[cur_arg+1],
                                                       pointer_levels[cur_arg+1]);
            if (expected_type == NULL) {
                goto error;
            }

            arg = PySequence_GetItem(argtypes, cur_arg);
            if (arg == NULL) {
                goto error;
            }

            if (arg != expected_type) {
#if PY3K
                PyErr_Format(PyExc_ValueError, "ctypes function argtypes[%d] != %S",
                             cur_arg, expected_type);
#else
                PyObject *expected_type_str, *s2;
                expected_type_str = PyObject_Str(expected_type);
                if (expected_type_str == NULL) {
                    goto error;
                }
                PyErr_Format(PyExc_ValueError, "ctypes function argtypes[%d] != %s",
                             cur_arg, PyBytes_AsString(expected_type_str));
                Py_DECREF(expected_type_str);
#endif
                goto error;
            }

            Py_DECREF(arg);
            Py_DECREF(expected_type);
            arg = NULL;
            expected_type = NULL;
        }
    }

    /* All checks passed */
    return_value = 1;
    callback->c_function = ccallback__get_ctypes_function_pointer(cfunc);
    callback->py_function = NULL;
    callback->user_data = user_data;
    goto done;

error:
    return_value = -1;

done:
    Py_XDECREF(expected_type);
    Py_XDECREF(argtypes);
    Py_XDECREF(restype);
    Py_XDECREF(arg);
    return return_value;
}


/*
 * Internal: obtain function pointers from cffi
 */

static int ccallback__check_cffi_type(PyObject *typeobj, char *name, Py_ssize_t length, int pointer_level)
{
    int return_value = 0;

    PyObject *cname = NULL;
    PyObject *kind = NULL;
    PyObject *arg = NULL;

    int ret;
    int ptr_level = 0;

    /* Count pointer level */
    arg = typeobj;
    Py_INCREF(arg);
    while (1) {
        PyObject *item;

        kind = PyObject_GetAttrString(arg, "kind");
        if (kind == NULL) {
            goto error;
        }

        ret = PyUString_isequal(kind, "pointer", 7);
        if (ret == -1) {
            goto error;
        }
        else if (!ret) {
            break;
        }

        Py_DECREF(kind);
        kind = NULL;

        ++ptr_level;

        item = PyObject_GetAttrString(arg, "item");
        if (item == NULL) {
            goto error;
        }
        Py_DECREF(arg);
        arg = item;
    }

    /* Check pointer level */
    if (ptr_level != pointer_level) {
        return_value = 0;
        goto done;
    }

    /* Check object cname */
    cname = PyObject_GetAttrString(arg, "cname");
    if (cname == NULL) {
        goto error;
    }

    ret = PyUString_isequal(cname, name, length);
    if (ret == -1) {
        goto error;
    }
    else if (!ret) {
        return_value = 0;
        goto done;
    }

    /* Match! */
    return_value = 1;
    goto done;

error:
    return_value = -1;

done:
    Py_XDECREF(cname);
    Py_XDECREF(kind);
    Py_XDECREF(arg);

    return return_value;
}


static int ccallback__parse_cffi(ccallback_t *callback, PyObject *callback_obj, char *signature)
{
    static PyObject *ffi = NULL, *cdata_type = NULL;
    static int import_failed = 0;

    int return_value = 0;
    PyObject *cffi_type = NULL;
    PyObject *attr = NULL;
    PyObject *result = NULL;
    PyObject *args = NULL;
    PyObject *arg = NULL;

    PyObject *cfunc;
    int ret;
    Py_ssize_t size;
    void *user_data;
    void *funcptr;

    int cur_arg, num_args;
    char *types[CCALLBACK_MAX_ARGS];
    int pointer_levels[CCALLBACK_MAX_ARGS];
    Py_ssize_t lengths[CCALLBACK_MAX_ARGS];

    if (import_failed) {
        return 0;
    }

    /* ffi = cffi.FFI() */
    if (ffi == NULL) {
        PyObject *cffi_module;
        cffi_module = PyImport_ImportModule("cffi");
        if (cffi_module == NULL) {
            /* cffi not available -- skip */
            PyErr_Clear();
            import_failed = 1;
            return 0;
        }

        ffi = PyObject_CallMethod(cffi_module, "FFI", "");
        Py_DECREF(cffi_module);
        if (ffi == NULL) {
            goto error;
        }
    }

    if (cdata_type == NULL) {
        cdata_type = PyObject_GetAttrString(ffi, "CData");
        if (cdata_type == NULL) {
            goto error;
        }
    }

    /* Deal with tuple form */
    if (PyObject_TypeCheck(callback_obj, (PyTypeObject *)cdata_type)) {
        cfunc = callback_obj;
        user_data = NULL;
    }
    else if (PyTuple_CheckExact(callback_obj) && PyTuple_GET_SIZE(callback_obj) == 2 &&
             PyObject_TypeCheck(PyTuple_GET_ITEM(callback_obj, 0), (PyTypeObject *)cdata_type))
    {
        int success;
        cfunc = PyTuple_GET_ITEM(callback_obj, 0);
        user_data = ccallback__parse_data_pointer(PyTuple_GET_ITEM(callback_obj, 1), &success);
        if (!success) {
            goto error;
        }
    }
    else {
        return_value = 0;
        goto done;
    }

    /* Get cffi object type */
    cffi_type = PyObject_CallMethod(ffi, "typeof", "O", cfunc);
    if (cffi_type == NULL) {
        goto error;
    }

    /* Check object type */
    attr = PyObject_GetAttrString(cffi_type, "kind");
    if (attr == NULL) {
        goto error;
    }

    ret = PyUString_isequal(attr, "function", 8);
    if (ret != 1) {
        goto error;
    }

    Py_DECREF(attr);
    attr = NULL;

    /* Parse expected signature */
    num_args = ccallback__parse_signature(signature, types, lengths, pointer_levels, CCALLBACK_MAX_ARGS);
    if (num_args < 0) {
        PyErr_Format(PyExc_RuntimeError, "scipy/ccallback: failed to parse signature '%s'",
                     signature);
        goto error;
    }

    /* Check return value */
    result = PyObject_GetAttrString(cffi_type, "result");
    if (result == NULL) {
        goto error;
    }

    ret = ccallback__check_cffi_type(result, types[0], lengths[0], pointer_levels[0]);
    if (ret == -1) {
        goto error;
    }
    else if (!ret) {
        return_value = -1;
        PyErr_SetString(PyExc_ValueError, "cffi object has mismatching return value");
        goto done;
    }

    Py_DECREF(result);
    result = NULL;

    /* Check arguments */
    args = PyObject_GetAttrString(cffi_type, "args");
    if (args == NULL) {
        goto error;
    }

    size = PySequence_Size(args);
    if (size == -1) {
        goto error;
    }
    else if (size != num_args) {
        PyErr_SetString(PyExc_ValueError, "cffi function takes wrong number of arguments");
        goto error;
    }

    for (cur_arg = 0; cur_arg < num_args; ++cur_arg) {
        arg = PySequence_GetItem(args, cur_arg);
        if (arg == NULL) {
            goto error;
        }

        ret = ccallback__check_cffi_type(arg, types[cur_arg+1], lengths[cur_arg+1],
                                         pointer_levels[cur_arg+1]);
        if (ret == -1) {
            goto error;
        }
        else if (!ret) {
            PyObject *expected_str;

            expected_str = PyUString_FromStringAndSize(types[cur_arg+1], lengths[cur_arg+1]);
            if (expected_str == NULL) {
                goto error;
            }
#if PY3K
            PyErr_Format(PyExc_ValueError, "cffi function args[%d] != %U", cur_arg, expected_str);
#else
            PyErr_Format(PyExc_ValueError, "ctypes function argtypes[%d] != %s",
                         cur_arg, PyBytes_AsString(expected_str));
#endif
            Py_DECREF(expected_str);
            goto error;
        }

        Py_DECREF(arg);
        arg = NULL;
    }

    /* Get function pointer */
    arg = PyObject_CallMethod(ffi, "cast", "sO", "unsigned long long", cfunc);
    if (arg == NULL) {
        goto error;
    }

    result = PyNumber_Long(arg);
    if (result == NULL) {
        goto error;
    }

    funcptr = PyLong_AsVoidPtr(result);
    if (PyErr_Occurred()) {
        goto error;
    }

    /* Finished */
    return_value = 1;
    callback->c_function = funcptr;
    callback->py_function = NULL;
    callback->user_data = user_data;
    goto done;

error:
    return_value = -1;

done:
    Py_XDECREF(cffi_type);
    Py_XDECREF(attr);
    Py_XDECREF(result);
    Py_XDECREF(args);
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
             PyUString_CheckExact(PyTuple_GET_ITEM(callback_obj, 1))
        ) {
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
        else {
            /* Try cffi */
            ret = ccallback__parse_cffi(callback, callback_obj, signature);
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
