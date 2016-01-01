/*
 * sparsetools.cxx
 *
 * Python module wrapping the sparsetools C++ routines.
 *
 * Each C++ routine is templated vs. an integer (I) and a data (T) parameter.
 * The `generate_sparsetools.py` script generates `*_impl.h` headers
 * that contain thunk functions with a datatype-based switch statement calling
 * each templated instantiation.
 *
 * `generate_sparsetools.py` also generates a PyMethodDef list of Python
 * routines and the corresponding functions call the thunk functions via
 * `call_thunk`.
 *
 * The `call_thunk` function below determines the templated I and T data types
 * based on the Python arguments. It then allocates arrays with pointers to
 * the raw data, with appropriate types, and calls the thunk function after
 * that.
 *
 * The types of arguments are specified by a "spec". This is given in a format
 * where one character represents one argument. The one-character values are
 * listed below in the call_spec function.
 */

#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include <Python.h>

#include <string>
#include <stdexcept>
#include <vector>
#include <cstdlib>

#include "numpy/ndarrayobject.h"

#include "sparsetools.h"

#define MAX_ARGS 16

#if PY_VERSION_HEX >= 0x03000000
#define PyInt_AsSsize_t PyLong_AsSsize_t
#define PyInt_FromSsize_t PyLong_FromSsize_t
#endif

static const int supported_I_typenums[] = {NPY_INT32, NPY_INT64};
static const int n_supported_I_typenums = sizeof(supported_I_typenums) / sizeof(int);

static const int supported_T_typenums[] = {NPY_BOOL,
                                           NPY_BYTE, NPY_UBYTE,
                                           NPY_SHORT, NPY_USHORT,
                                           NPY_INT, NPY_UINT,
                                           NPY_LONG, NPY_ULONG,
                                           NPY_LONGLONG, NPY_ULONGLONG,
                                           NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE,
                                           NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE};
static const int n_supported_T_typenums = sizeof(supported_T_typenums) / sizeof(int);

static PyObject *array_from_std_vector_and_free(int typenum, void *p);
static void *allocate_std_vector_typenum(int typenum);
static void free_std_vector_typenum(int typenum, void *p);
static PyObject *c_array_from_object(PyObject *obj, int typenum, int is_output);


/*
 * Call a thunk function, dealing with input and output arrays.
 *
 * Resolves the templated <integer> and <data> dtypes from the `args` argument
 * list.
 *
 * Parameters
 * ----------
 * ret_spec : {'i', 'v'}
 *     Return value spec. 'i' for integer, 'v' for void.
 * spec
 *     String whose each character specifies a types of an
 *     argument:
 *
 *     'i': <integer> scalar
 *     'I': <integer> array
 *     'T': <data> array
 *     'V': std::vector<integer>
 *     'W': std::vector<data>
 *     'B': npy_bool array
 *     '*': indicates that the next argument is an output argument
 * thunk : Py_ssize_t thunk(int I_typenum, int T_typenum, void **)
 *     Thunk function to call. It is passed a void** array of pointers to
 *     arguments, constructed according to `spec`. The types of data pointed
 *     to by each element agree with I_typenum and T_typenum, or are bools.
 * args
 *     Python tuple containing unprocessed arguments.
 *
 * Returns
 * -------
 * return_value
 *     The Python return value
 *
 */
NPY_VISIBILITY_HIDDEN PyObject *
call_thunk(char ret_spec, const char *spec, thunk_t *thunk, PyObject *args)
{
    void *arg_list[MAX_ARGS];
    PyObject *arg_arrays[MAX_ARGS];
    int is_output[MAX_ARGS];
    PyObject *return_value = NULL;
    int I_typenum = NPY_INT32;
    int T_typenum = -1;
    int VW_count = 0;
    int I_in_arglist = 0;
    int T_in_arglist = 0;
    int next_is_output = 0;
    int j, k, arg_j;
    const char *p;
    Py_ssize_t ret;
    Py_ssize_t max_array_size = 0;
    NPY_BEGIN_THREADS_DEF;

    if (!PyTuple_Check(args)) {
        PyErr_SetString(PyExc_ValueError, "args is not a tuple");
        return NULL;
    }

    for (j = 0; j < MAX_ARGS; ++j) {
        arg_list[j] = NULL;
        arg_arrays[j] = NULL;
        is_output[j] = 0;
    }


    /*
     * Detect data types in the signature
     */
    arg_j = 0;
    j = 0;
    for (p = spec; *p != '\0'; ++p, ++j, ++arg_j) {
        const int *supported_typenums;
        int n_supported_typenums;
        int cur_typenum;
        PyObject *arg;
        PyArray_Descr *dtype;

        if (j >= MAX_ARGS) {
            PyErr_SetString(PyExc_ValueError,
                            "internal error: too many arguments in spec");
            goto fail;
        }

        is_output[j] = next_is_output;
        next_is_output = 0;

        switch (*p) {
        case '*':
            next_is_output = 1;
            --j;
            --arg_j;
            continue;
        case 'i':
            /* Integer scalars */
            arg = PyTuple_GetItem(args, arg_j);
            if (arg == NULL) {
                goto fail;
            }
            Py_INCREF(arg);
            arg_arrays[j] = arg;
            continue;
        case 'I':
            /* Integer arrays */
            supported_typenums = supported_I_typenums;
            n_supported_typenums = n_supported_I_typenums;
            cur_typenum = I_typenum;
            I_in_arglist = 1;
            break;
        case 'T':
            /* Data arrays */
            supported_typenums = supported_T_typenums;
            n_supported_typenums = n_supported_T_typenums;
            cur_typenum = T_typenum;
            T_in_arglist = 1;
            break;
        case 'B':
            /* Boolean arrays */
            arg = PyTuple_GetItem(args, arg_j);
            if (arg == NULL) {
                goto fail;
            }
            arg_arrays[j] = c_array_from_object(arg, NPY_BOOL, is_output[j]);
            if (arg_arrays[j] == NULL) {
                goto fail;
            }
            continue;
        case 'V':
            /* std::vector integer output array */
            I_in_arglist = 1;
            --arg_j;
            VW_count += 1;
            continue;
        case 'W':
            /* std::vector data output array */
            T_in_arglist = 1;
            --arg_j;
            VW_count += 1;
            continue;
        default:
            PyErr_SetString(PyExc_ValueError, "unknown character in spec");
            goto fail;
        }

        arg = PyTuple_GetItem(args, arg_j);
        if (arg == NULL) {
            goto fail;
        }
        arg_arrays[j] = c_array_from_object(arg, -1, is_output[j]);
        if (arg_arrays[j] == NULL) {
            goto fail;
        }

        /* Find a compatible supported data type */
        dtype = PyArray_DESCR(arg_arrays[j]);
        for (k = 0; k < n_supported_typenums; ++k) {
            if (PyArray_CanCastSafely(dtype->type_num, supported_typenums[k]) &&
                (cur_typenum == -1 || PyArray_CanCastSafely(cur_typenum, supported_typenums[k])))
            {
                cur_typenum = supported_typenums[k];
                break;
            }
        }
        if (k == n_supported_typenums) {
            PyErr_SetString(PyExc_ValueError,
                            "unsupported data types in input");
            goto fail;
        }

        if (*p == 'I') {
            I_typenum = cur_typenum;
        }
        else {
            T_typenum = cur_typenum;
        }
    }

    if (arg_j != PyTuple_Size(args)) {
        PyErr_SetString(PyExc_ValueError, "too many arguments");
        goto fail;
    }

    if ((I_in_arglist && I_typenum == -1) ||
        (T_in_arglist && T_typenum == -1)) {
        PyErr_SetString(PyExc_ValueError,
                        "unsupported data types in input");
        goto fail;
    }


    /*
     * Cast and extract argument arrays
     */
    j = 0;
    for (p = spec; *p != '\0'; ++p, ++j) {
        PyObject *arg;
        int cur_typenum;

        if (*p == '*') {
            --j;
            continue;
        }
        else if (*p == 'i') {
            /* Integer scalars */
            Py_ssize_t value;

            value = PyInt_AsSsize_t(arg_arrays[j]);
            if (PyErr_Occurred()) {
                goto fail;
            }

            if (PyArray_EquivTypenums(I_typenum, NPY_INT64)
                    && value == (npy_int64)value) {
                arg_list[j] = std::malloc(sizeof(npy_int64));
                *(npy_int64*)arg_list[j] = (npy_int64)value;
            }
            else if (PyArray_EquivTypenums(I_typenum, NPY_INT32)
                     && value == (npy_int32)value) {
                arg_list[j] = std::malloc(sizeof(npy_int32));
                *(npy_int32*)arg_list[j] = (npy_int32)value;
            }
            else {
                PyErr_SetString(PyExc_ValueError,
                                "could not convert integer scalar");
                goto fail;
            }
            continue;
        }
        else if (*p == 'B') {
            /* Boolean arrays already cast */
        }
        else if (*p == 'V') {
            arg_list[j] = allocate_std_vector_typenum(I_typenum);
            if (arg_list[j] == NULL) {
                goto fail;
            }
            continue;
        }
        else if (*p == 'W') {
            arg_list[j] = allocate_std_vector_typenum(T_typenum);
            if (arg_list[j] == NULL) {
                goto fail;
            }
            continue;
        }
        else {
            cur_typenum = (*p == 'I' || *p == 'i') ? I_typenum : T_typenum;

            /* Cast if necessary */
            arg = arg_arrays[j];
            if (PyArray_EquivTypenums(PyArray_DESCR(arg)->type_num, cur_typenum)) {
                /* No cast needed. */
            }
            else if (!is_output[j] || PyArray_CanCastSafely(cur_typenum, PyArray_DESCR(arg)->type_num)) {
                /* Cast needed. Output arrays require safe cast back. */
                arg_arrays[j] = c_array_from_object(arg, cur_typenum, is_output[j]);
                Py_DECREF(arg);
                if (arg_arrays[j] == NULL) {
                    goto fail;
                }
            }
            else {
                /* Cast back into output array was not safe. */
                PyErr_SetString(PyExc_ValueError,
                                "Output dtype not compatible with inputs.");
                goto fail;
            }
        }

        /* Grab value */
        arg_list[j] = PyArray_DATA(arg_arrays[j]);

        /* Find maximum array size */
        if (PyArray_SIZE(arg_arrays[j]) > max_array_size) {
            max_array_size = PyArray_SIZE(arg_arrays[j]);
        }
    }


    /*
     * Call thunk
     */
    if (max_array_size > 100) {
        /* Threshold GIL release: it's not a free operation */
        NPY_BEGIN_THREADS;
    }
    try {
        ret = thunk(I_typenum, T_typenum, arg_list);
        NPY_END_THREADS;
    } catch (const std::bad_alloc &e) {
        NPY_END_THREADS;
        PyErr_SetString(PyExc_MemoryError, e.what());
        goto fail;
    } catch (const std::exception &e) {
        NPY_END_THREADS;
        PyErr_SetString(PyExc_RuntimeError, e.what());
        goto fail;
    }

    /*
     * Generate return value;
     */

    switch (ret_spec) {
    case 'i':
        return_value = PyInt_FromSsize_t(ret);
        break;
    case 'v':
        Py_INCREF(Py_None);
        return_value = Py_None;
        break;
    default:
        PyErr_SetString(PyExc_ValueError,
                        "internal error: invalid return value spec");
    }

    /*
     * Convert any std::vector output arrays to arrays
     */
    if (VW_count > 0) {
        PyObject *new_ret;
        PyObject *old_ret = return_value;
        int pos;

        return_value = NULL;

        new_ret = PyTuple_New(VW_count + (old_ret == Py_None ? 0 : 1));
        if (new_ret == NULL) {
            goto fail;
        }
        if (old_ret != Py_None) {
            PyTuple_SET_ITEM(new_ret, 0, old_ret);
            pos = 1;
        }
        else {
            Py_DECREF(old_ret);
            pos = 0;
        }

        j = 0;
        for (p = spec; *p != '\0'; ++p, ++j) {
            if (*p == '*') {
                --j;
                continue;
            }
            else if (*p == 'V' || *p == 'W') {
                PyObject *arg;
                if (*p == 'V') {
                    arg = array_from_std_vector_and_free(I_typenum, arg_list[j]);
                } else {
                    arg = array_from_std_vector_and_free(T_typenum, arg_list[j]);
                }
                arg_list[j] = NULL;
                if (arg == NULL) {
                    Py_XDECREF(new_ret);
                    goto fail;
                }
                PyTuple_SET_ITEM(new_ret, pos, arg);
                ++pos;
            }
        }

        return_value = new_ret;
    }


fail:
    /*
     * Cleanup
     */
    for (j = 0; j < MAX_ARGS; ++j) {
        Py_XDECREF(arg_arrays[j]);
        if (spec[j] == 'i' && arg_list[j] != NULL) {
            std::free(arg_list[j]);
        }
        else if (spec[j] == 'V' && arg_list[j] != NULL) {
            free_std_vector_typenum(I_typenum, arg_list[j]);
        }
        else if (spec[j] == 'W' && arg_list[j] != NULL) {
            free_std_vector_typenum(T_typenum, arg_list[j]);
        }
    }
    return return_value;
}


/*
 * Helper functions for dealing with std::vector templated instantiation.
 */

static void *allocate_std_vector_typenum(int typenum)
{
#define PROCESS(ntype, ctype)                                   \
    if (PyArray_EquivTypenums(typenum, ntype)) {                \
        return (void*)(new std::vector<ctype>());               \
    }

    try {
        PROCESS(NPY_BOOL, npy_bool_wrapper);
        PROCESS(NPY_BYTE, npy_byte);
        PROCESS(NPY_UBYTE, npy_ubyte);
        PROCESS(NPY_SHORT, npy_short);
        PROCESS(NPY_USHORT, npy_ushort);
        PROCESS(NPY_INT, npy_int);
        PROCESS(NPY_UINT, npy_uint);
        PROCESS(NPY_LONG, npy_long);
        PROCESS(NPY_ULONG, npy_ulong);
        PROCESS(NPY_LONGLONG, npy_longlong);
        PROCESS(NPY_ULONGLONG, npy_ulonglong);
        PROCESS(NPY_FLOAT, npy_float);
        PROCESS(NPY_DOUBLE, npy_double);
        PROCESS(NPY_LONGDOUBLE, npy_longdouble);
        PROCESS(NPY_CFLOAT, npy_cfloat_wrapper);
        PROCESS(NPY_CDOUBLE, npy_cdouble_wrapper);
        PROCESS(NPY_CLONGDOUBLE, npy_clongdouble_wrapper);
    } catch (std::exception &e) {
        /* failed */
    }

#undef PROCESS

    PyErr_SetString(PyExc_RuntimeError,
                    "failed to allocate std::vector");
    return NULL;
}

static void free_std_vector_typenum(int typenum, void *p)
{
#define PROCESS(ntype, ctype)                                   \
    if (PyArray_EquivTypenums(typenum, ntype)) {                \
        delete ((std::vector<ctype>*)p);                        \
    }

    PROCESS(NPY_BOOL, npy_bool_wrapper);
    PROCESS(NPY_BYTE, npy_byte);
    PROCESS(NPY_UBYTE, npy_ubyte);
    PROCESS(NPY_SHORT, npy_short);
    PROCESS(NPY_USHORT, npy_ushort);
    PROCESS(NPY_INT, npy_int);
    PROCESS(NPY_UINT, npy_uint);
    PROCESS(NPY_LONG, npy_long);
    PROCESS(NPY_ULONG, npy_ulong);
    PROCESS(NPY_LONGLONG, npy_longlong);
    PROCESS(NPY_ULONGLONG, npy_ulonglong);
    PROCESS(NPY_FLOAT, npy_float);
    PROCESS(NPY_DOUBLE, npy_double);
    PROCESS(NPY_LONGDOUBLE, npy_longdouble);
    PROCESS(NPY_CFLOAT, npy_cfloat_wrapper);
    PROCESS(NPY_CDOUBLE, npy_cdouble_wrapper);
    PROCESS(NPY_CLONGDOUBLE, npy_clongdouble_wrapper);

#undef PROCESS
}

static PyObject *array_from_std_vector_and_free(int typenum, void *p)
{
#define PROCESS(ntype, ctype)                                   \
    if (PyArray_EquivTypenums(typenum, ntype)) {                \
        std::vector<ctype> *v = (std::vector<ctype>*)p;         \
        npy_intp length = v->size();                            \
        PyObject *obj = PyArray_SimpleNew(1, &length, typenum); \
        if (length > 0) {                                       \
            memcpy(PyArray_DATA(obj), &((*v)[0]),               \
                   sizeof(ctype)*length);                       \
        }                                                       \
        delete v;                                               \
        return obj;                                             \
    }

    PROCESS(NPY_BOOL, npy_bool_wrapper);
    PROCESS(NPY_BYTE, npy_byte);
    PROCESS(NPY_UBYTE, npy_ubyte);
    PROCESS(NPY_SHORT, npy_short);
    PROCESS(NPY_USHORT, npy_ushort);
    PROCESS(NPY_INT, npy_int);
    PROCESS(NPY_UINT, npy_uint);
    PROCESS(NPY_LONG, npy_long);
    PROCESS(NPY_ULONG, npy_ulong);
    PROCESS(NPY_LONGLONG, npy_longlong);
    PROCESS(NPY_ULONGLONG, npy_ulonglong);
    PROCESS(NPY_FLOAT, npy_float);
    PROCESS(NPY_DOUBLE, npy_double);
    PROCESS(NPY_LONGDOUBLE, npy_longdouble);
    PROCESS(NPY_CFLOAT, npy_cfloat_wrapper);
    PROCESS(NPY_CDOUBLE, npy_cdouble_wrapper);
    PROCESS(NPY_CLONGDOUBLE, npy_clongdouble_wrapper);

#undef PROCESS

    PyErr_SetString(PyExc_RuntimeError,
                    "failed to convert std::vector output array");
    return NULL;
}

static PyObject *c_array_from_object(PyObject *obj, int typenum, int is_output)
{
    if (!is_output) {
        if (typenum == -1) {
            return PyArray_FROM_OF(obj, NPY_C_CONTIGUOUS);
        }
        else {
            return PyArray_FROM_OTF(obj, typenum, NPY_C_CONTIGUOUS);
        }
    }
    else {
        if (typenum == -1) {
            return PyArray_FROM_OF(obj, NPY_C_CONTIGUOUS|NPY_WRITEABLE|NPY_UPDATEIFCOPY);
        }
        else {
            return PyArray_FROM_OTF(obj, typenum, NPY_C_CONTIGUOUS|NPY_WRITEABLE|NPY_UPDATEIFCOPY);
        }
    }
}


/*
 * Python module initialization
 */

extern "C" {

#include "sparsetools_impl.h"

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_sparsetools",
    NULL,
    -1,
    sparsetools_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit__sparsetools(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    import_array();
    return m;
}
#else
PyMODINIT_FUNC init_sparsetools(void) {
    PyObject *m;
    m = Py_InitModule("_sparsetools", sparsetools_methods);
    import_array();
    if (m == NULL) {
        Py_FatalError("can't initialize module _sparsetools");
    }
}
#endif

} /* extern "C" */
