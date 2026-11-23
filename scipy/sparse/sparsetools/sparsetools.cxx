/*
 * sparsetools.cxx
 *
 * Python module wrapping the sparsetools C++ routines.
 *
 * Each C++ routine is templated vs. an integer (I) and a data (T) parameter.
 * Each compilation unit (`bsr.cxx`, `csr.cxx`, `csc.cxx`, `other.cxx`) builds
 * a `constexpr` table of `SparsetoolsRoutine` descriptors whose argument specs
 * and per-(I, T) kernel function pointers are deduced at compile time from
 * C++ function signatures (see `sparsetools.h`).
 *
 * `call_routine` below determines the templated I and T data types based on
 * the Python arguments, prepares contiguous NumPy array buffers and scalars,
 * and dispatches to the matching compiled kernel.
 */

#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include <pybind11/pybind11.h>

#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include "numpy/ndarrayobject.h"

#include "sparsetools.h"
#include "util.h"

namespace py = pybind11;

#define MAX_ARGS 16

static const int supported_I_typenums[] = {NPY_INT32, NPY_INT64};
static const int n_supported_I_typenums = sizeof(supported_I_typenums) / sizeof(int);

static const int supported_T_typenums[] = {
#define PROCESS(ntype, ctype) ntype,
    SPTOOLS_FOR_EACH_DATA_TYPE_CODE(PROCESS)
#undef PROCESS
};
static const int n_supported_T_typenums = sizeof(supported_T_typenums) / sizeof(int);

static PyObject *array_from_std_vector_and_free(int typenum, void *p);
static void *allocate_std_vector_typenum(int typenum);
static void free_std_vector_typenum(int typenum, void *p);
static PyObject *c_array_from_object(PyObject *obj, int typenum, int is_output);

static inline int get_I_index(int I_typenum) {
    if (PyArray_EquivTypenums(I_typenum, NPY_INT32)) {
        return sparsetools_detail::kIndexTypeIndex<npy_int32>;
    }
    if (PyArray_EquivTypenums(I_typenum, NPY_INT64)) {
        return sparsetools_detail::kIndexTypeIndex<npy_int64>;
    }
    return -1;
}

static inline int get_T_index(int T_typenum) {
    switch (T_typenum) {
#define PROCESS(ntype, ctype) \
    case ntype: return sparsetools_detail::kDataTypeIndex<ctype>;
    SPTOOLS_FOR_EACH_DATA_TYPE_CODE(PROCESS)
#undef PROCESS
    default:
        return -1;
    }
}

/*
 * Call a SparsetoolsRoutine, dealing with input and output arrays.
 *
 * Resolves the templated <integer> and <data> dtypes from the `args` argument
 * list via the routine's arg spec:
 *     'i': <integer> scalar
 *     'l': 64-bit integer scalar
 *     'I': <integer> array
 *     'T': <data> array
 *     'V': std::vector<integer>
 *     'W': std::vector<data>
 *     'B': npy_bool array
 *     '*': indicates that the next argument is an output argument
 */
static PyObject *
call_routine(const SparsetoolsRoutine& routine, PyObject *args)
{
    const char *spec = routine.arg_spec;
    const char ret_spec = routine.ret_spec;
    void *arg_list[MAX_ARGS];
    npy_int64 scalar_vals[MAX_ARGS];
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
    PY_LONG_LONG ret;
    Py_ssize_t max_array_size = 0;
    bool should_release_gil = false;
    SparsetoolsKernelFn kernel = NULL;
    NPY_BEGIN_THREADS_DEF;

    if (!PyTuple_Check(args)) {
        PyErr_SetString(PyExc_ValueError, "args is not a tuple");
        return NULL;
    }

    for (j = 0; j < MAX_ARGS; ++j) {
        arg_list[j] = NULL;
        scalar_vals[j] = 0;
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
        case 'l':
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

        for (k = 0; k < n_supported_typenums; ++k) {
            if (PyArray_CanCastSafely(PyArray_TYPE((PyArrayObject *)arg_arrays[j]), supported_typenums[k]) &&
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
        else if (*p == 'i' || *p == 'l') {
            /* Integer scalars */
            PY_LONG_LONG value = PyLong_AsLongLong(arg_arrays[j]);
            if (PyErr_Occurred()) {
                goto fail;
            }

            if ((*p == 'l' || PyArray_EquivTypenums(I_typenum, NPY_INT64))
                    && value == (npy_int64)value) {
                scalar_vals[j] = (npy_int64)value;
                arg_list[j] = &scalar_vals[j];
            }
            else if (*p == 'i' && PyArray_EquivTypenums(I_typenum, NPY_INT32)
                     && value == (npy_int32)value) {
                scalar_vals[j] = (npy_int64)value;
                arg_list[j] = &scalar_vals[j];
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
            if (PyArray_EquivTypenums(PyArray_TYPE((PyArrayObject *)arg), cur_typenum)) {
                /* No cast needed. */
            }
            else if (!is_output[j] || PyArray_CanCastSafely(cur_typenum, PyArray_TYPE((PyArrayObject *)arg))) {
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
        arg_list[j] = PyArray_DATA((PyArrayObject *)arg_arrays[j]);

        /* Find maximum array size */
        if (PyArray_SIZE((PyArrayObject *)arg_arrays[j]) > max_array_size) {
            max_array_size = PyArray_SIZE((PyArrayObject *)arg_arrays[j]);
        }
    }


    /*
     * Select the kernel to invoke
     */
    switch (routine.dispatch_kind) {
    case SparsetoolsDispatchKind::kNone:
        kernel = routine.kernels[0];
        break;
    case SparsetoolsDispatchKind::kIndex: {
        int i_idx = get_I_index(I_typenum);
        if (i_idx >= 0) {
            kernel = routine.kernels[i_idx];
        }
        break;
    }
    case SparsetoolsDispatchKind::kIndexData: {
        int i_idx = get_I_index(I_typenum);
        int t_idx = get_T_index(T_typenum);
        if (i_idx >= 0 && t_idx >= 0) {
            kernel = routine.kernels[
                static_cast<std::size_t>(i_idx) *
                    sparsetools_detail::kNumDataTypes +
                static_cast<std::size_t>(t_idx)];
        }
        break;
    }
    }
    if (kernel == NULL) {
        PyErr_SetString(PyExc_RuntimeError,
                        "internal error: invalid argument typenums");
        goto fail;
    }

    /*
     * Invoke the kenel
     */
    should_release_gil = max_array_size > 100;
    if (should_release_gil) {
        /* Threshold GIL release: it's not a free operation */
        NPY_BEGIN_THREADS;
    }
    try {
        ret = kernel(arg_list);
        if (should_release_gil) NPY_END_THREADS;
    } catch (const std::bad_alloc &e) {
        if (should_release_gil) NPY_END_THREADS;
        PyErr_SetString(PyExc_MemoryError, e.what());
        goto fail;
    } catch (const std::exception &e) {
        if (should_release_gil) NPY_END_THREADS;
        PyErr_SetString(PyExc_RuntimeError, e.what());
        goto fail;
    }

    /*
     * Generate return value;
     */

    switch (ret_spec) {
    case 'i':
    case 'l':
        return_value = PyLong_FromLongLong(ret);
        if (return_value == NULL) goto fail;
        break;
    case 'v':
        Py_INCREF(Py_None);
        return_value = Py_None;
        break;
    default:
        PyErr_SetString(PyExc_ValueError,
                        "internal error: invalid return value spec");
        goto fail;
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
            Py_DECREF(old_ret);
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
                if (arg == NULL) {
                    Py_XDECREF(new_ret);
                    goto fail;
                }
                PyTuple_SET_ITEM(new_ret, pos, arg);
                arg_list[j] = NULL;
                ++pos;
            }
        }

        return_value = new_ret;
    }


fail:
    /*
     * Cleanup
     */
    for (j = 0, p = spec; *p != '\0'; ++p, ++j) {
        if (*p == '*') {
            --j;
            continue;
        }
        if (is_output[j] && arg_arrays[j] != NULL && PyArray_Check(arg_arrays[j])) {
            PyArray_ResolveWritebackIfCopy((PyArrayObject *)arg_arrays[j]);
        }
        Py_XDECREF(arg_arrays[j]);
        if (*p == 'V' && arg_list[j] != NULL) {
            free_std_vector_typenum(I_typenum, arg_list[j]);
        }
        else if (*p == 'W' && arg_list[j] != NULL) {
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
    try {
        switch (typenum) {
#define PROCESS(ntype, ctype)                 \
        case ntype:                           \
            return static_cast<void *>(new std::vector<ctype>());
        SPTOOLS_FOR_EACH_DATA_TYPE_CODE(PROCESS)
#undef PROCESS
        default: break;
        }
    } catch (const std::exception &) {
        /* failed */
    }

    PyErr_SetString(PyExc_RuntimeError,
                    "failed to allocate std::vector");
    return NULL;
}

static void free_std_vector_typenum(int typenum, void *p)
{
    switch (typenum) {
#define PROCESS(ntype, ctype)                       \
    case ntype:                                     \
        delete static_cast<std::vector<ctype> *>(p); \
        return;
    SPTOOLS_FOR_EACH_DATA_TYPE_CODE(PROCESS)
#undef PROCESS
    default: return;
    }
}

static PyObject *array_from_std_vector_and_free(int typenum, void *p)
{
    switch (typenum) {
#define PROCESS(ntype, ctype)                                                 \
    case ntype: {                                                             \
        auto *v = static_cast<std::vector<ctype> *>(p);                       \
        npy_intp length = static_cast<npy_intp>(v->size());                   \
        PyObject *obj = PyArray_SimpleNew(1, &length, typenum);               \
        if (obj == NULL) return NULL;                                         \
        if (length > 0) {                                                     \
            std::memcpy(PyArray_DATA((PyArrayObject *)obj), v->data(),        \
                        sizeof(ctype) * static_cast<std::size_t>(length));    \
        }                                                                     \
        delete v;                                                             \
        return obj;                                                           \
     }
    SPTOOLS_FOR_EACH_DATA_TYPE_CODE(PROCESS)
#undef PROCESS
    default: break;
    }

    PyErr_SetString(PyExc_RuntimeError,
                    "failed to convert std::vector output array");
    return NULL;
}

static PyObject *c_array_from_object(PyObject *obj, int typenum, int is_output)
{
    if (!is_output) {
        if (typenum == -1) {
            return PyArray_FROM_OF(obj, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_NOTSWAPPED);
        }
        else {
            return PyArray_FROM_OTF(obj, typenum, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_NOTSWAPPED);
        }
    }
    else {
        int flags = NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_WRITEABLE|NPY_ARRAY_WRITEBACKIFCOPY|NPY_ARRAY_NOTSWAPPED;
        if (typenum == -1) {
            return PyArray_FROM_OF(obj, flags);
        }
        else {
            return PyArray_FROM_OTF(obj, typenum, flags);
        }
    }
}


/*
 * Python module initialization
 */
#if PYBIND11_VERSION_MAJOR >= 3
PYBIND11_MODULE(_sparsetools, m,
                py::mod_gil_not_used(),
                py::multiple_interpreters::per_interpreter_gil())
#else
PYBIND11_MODULE(_sparsetools, m, py::mod_gil_not_used())
#endif
{
    if (_import_array() != 0) {
      throw py::error_already_set();
    }
    for (const SparsetoolsRoutineSpan unit : {
             get_bsr_routines(),
             get_csr_routines(),
             get_csc_routines(),
             get_other_routines(),
         }) {
        for (std::size_t i = 0; i < unit.size; ++i) {
            const SparsetoolsRoutine &r = unit.data[i];
            m.def(r.name, [&r](py::args args) {
                PyObject *result = call_routine(r, args.ptr());
                if (!result) {
                    throw py::error_already_set();
                }
                return py::reinterpret_steal<py::object>(result);
            });
        }
    }
}

