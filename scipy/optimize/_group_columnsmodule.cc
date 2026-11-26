/*
 * Column grouping for the sparse finite difference Jacobian estimation of
 * `scipy.optimize._numdiff`; see `_numdiff.group_columns`.
 *
 * Two columns can share a group if no row has a nonzero entry in both. The columns are taken
 * greedily, in order: an ungrouped column opens a new group, which every later ungrouped column
 * joins that does not intersect the union of the columns in the group so far.
 *
 * The loops run without the GIL.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <algorithm>
#include <cstdint>
#include <memory>

#include "numpy/arrayobject.h"


namespace {

struct ArrayDeleter {
    void operator()(PyArrayObject *p) const noexcept { Py_XDECREF(p); }
};


/*
 * Owning reference to a NumPy array; `release()` hands the reference over.
 */
using OwnedArray = std::unique_ptr<PyArrayObject, ArrayDeleter>;


// Releases the GIL for its lifetime.
class GILRelease {
  public:
    GILRelease() noexcept : state_(PyEval_SaveThread()) {}
    ~GILRelease() { PyEval_RestoreThread(state_); }
    GILRelease(const GILRelease &) = delete;
    GILRelease &operator=(const GILRelease &) = delete;
  private:
    PyThreadState *state_;
};


/*
 * `obj` as an aligned, native-endian array of `typenum` with `ndim` dimensions, and contiguous
 * if `flags` asks for it. Copies only when needed; casts only when the cast is safe.
 */
OwnedArray as_array(PyObject *obj, int typenum, int ndim, int flags, const char *name)
{
    OwnedArray arr(reinterpret_cast<PyArrayObject *>(PyArray_FROMANY(obj, typenum, 0, 0, flags)));
    if (arr && PyArray_NDIM(arr.get()) != ndim) {
        PyErr_Format(PyExc_ValueError, "`%s` must be a %d-D array, got %d-D.",
                     name, ndim, PyArray_NDIM(arr.get()));
        arr.reset();
    }
    return arr;
}


// A new (n,) intp array filled with -1, for the group of each column.
OwnedArray new_groups(npy_intp n)
{
    OwnedArray groups(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(1, &n, NPY_INTP, 0)));
    if (groups) {
        npy_intp *g = static_cast<npy_intp *>(PyArray_DATA(groups.get()));
        std::fill(g, g + n, npy_intp(-1));
    }
    return groups;
}


/*
 * Greedy grouping of the columns of the dense m-by-n structure at `A`, with byte strides `s0`
 * between rows and `s1` between columns; an entry > 0 is a nonzero. `u` is a workspace of m
 * elements.
 */
void group_dense_impl(const char *A, const npy_intp s0, const npy_intp s1, const npy_intp m,
                      const npy_intp n, npy_intp *groups, unsigned char *u) noexcept
{
    auto nonzero = [&](npy_intp k, npy_intp j) -> unsigned char {
        return *reinterpret_cast<const int32_t *>(A + k * s0 + j * s1) > 0;
    };
    npy_intp current_group = 0;
    for (npy_intp i = 0; i < n; i++) {
        if (groups[i] >= 0) { continue; }  // A group was already assigned.
        groups[i] = current_group;
        bool all_grouped = true;

        // Here we store the union of grouped columns.
        for (npy_intp k = 0; k < m; k++) { u[k] = nonzero(k, i); }

        for (npy_intp j = i + 1; j < n; j++) {
            if (groups[j] >= 0) { continue; }
            all_grouped = false;

            // Add the j-th column to the group if it does not intersect the union.
            bool intersect = false;
            for (npy_intp k = 0; k < m; k++) {
                if (u[k] && nonzero(k, j)) { intersect = true; break; }
            }
            if (!intersect) {
                for (npy_intp k = 0; k < m; k++) { u[k] |= nonzero(k, j); }
                groups[j] = current_group;
            }
        }
        if (all_grouped) { break; }
        current_group++;
    }
}


/*
 * Greedy grouping of the columns of an m-by-n CSC structure with valid `indices`, `indptr`.
 * `u` is a workspace of m elements.
 */
template<typename I>
void group_sparse_impl(const I *indices, const I *indptr, const npy_intp m, const npy_intp n,
                       npy_intp *groups, unsigned char *u) noexcept
{
    npy_intp current_group = 0;
    for (npy_intp i = 0; i < n; i++) {
        if (groups[i] >= 0) { continue; }
        groups[i] = current_group;
        bool all_grouped = true;

        std::fill(u, u + m, static_cast<unsigned char>(0));
        for (I k = indptr[i]; k < indptr[i + 1]; k++) { u[indices[k]] = 1; }

        for (npy_intp j = i + 1; j < n; j++) {
            if (groups[j] >= 0) { continue; }
            all_grouped = false;

            bool intersect = false;
            for (I k = indptr[j]; k < indptr[j + 1]; k++) {
                if (u[indices[k]]) { intersect = true; break; }
            }
            if (!intersect) {
                for (I k = indptr[j]; k < indptr[j + 1]; k++) { u[indices[k]] = 1; }
                groups[j] = current_group;
            }
        }
        if (all_grouped) { break; }
        current_group++;
    }
}


/*
 * Whether `indptr` (n + 1 entries) is nondecreasing within [0, nnz], and the `indices` it
 * refers to are within [0, m). The loops index the workspace with them.
 */
template<typename I>
bool is_valid_csc(const I *indices, const npy_intp nnz, const I *indptr, const npy_intp m,
                  const npy_intp n) noexcept
{
    if (indptr[0] < 0 || static_cast<npy_intp>(indptr[n]) > nnz) { return false; }
    for (npy_intp j = 0; j < n; j++) {
        if (indptr[j] > indptr[j + 1]) { return false; }
    }
    for (npy_intp k = 0; k < static_cast<npy_intp>(indptr[n]); k++) {
        if (indices[k] < 0 || static_cast<npy_intp>(indices[k]) >= m) { return false; }
    }
    return true;
}


bool is_int32_array(PyObject *obj)
{
    if (!PyArray_Check(obj)) { return false; }
    PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(obj);
    return PyArray_ISSIGNED(arr) && PyArray_ITEMSIZE(arr) == 4;
}


PyDoc_STRVAR(doc_group_dense,
"group_dense(m, n, A)\n"
"\n"
"Group the columns of a dense structure for sparse finite differencing.\n"
"\n"
"Parameters\n"
"----------\n"
"m, n : int\n"
"    Shape of `A`.\n"
"A : (m, n) array_like of int32\n"
"    The structure; an entry > 0 is a nonzero.\n"
"\n"
"Returns\n"
"-------\n"
"groups : (n,) ndarray of intp\n"
"    The group of each column, from 0 to the number of groups minus 1.\n");


PyObject *py_group_dense(PyObject *Py_UNUSED(self), PyObject *args)
{
    Py_ssize_t m, n;
    PyObject *A_obj;
    if (!PyArg_ParseTuple(args, "nnO", &m, &n, &A_obj)) { return NULL; }
    if (m < 0 || n < 0) {
        PyErr_SetString(PyExc_ValueError, "`m` and `n` must be nonnegative.");
        return NULL;
    }
    // Any strides, without a copy; `group_columns` passes A[:, order], which is F-contiguous,
    // so the column scans are contiguous.
    OwnedArray A = as_array(A_obj, NPY_INT32, 2, NPY_ARRAY_ALIGNED, "A");
    if (!A) { return NULL; }
    if (PyArray_DIM(A.get(), 0) != m || PyArray_DIM(A.get(), 1) != n) {
        PyErr_SetString(PyExc_ValueError, "`A` must have the shape (m, n).");
        return NULL;
    }
    OwnedArray groups = new_groups(n);
    if (!groups) { return NULL; }

    unsigned char *u = static_cast<unsigned char *>(PyMem_RawMalloc(m > 0 ? m : 1));
    if (u == NULL) { return PyErr_NoMemory(); }
    {
        GILRelease nogil;
        group_dense_impl(PyArray_BYTES(A.get()), PyArray_STRIDE(A.get(), 0),
                         PyArray_STRIDE(A.get(), 1), m, n,
                         static_cast<npy_intp *>(PyArray_DATA(groups.get())), u);
    }
    PyMem_RawFree(u);
    return reinterpret_cast<PyObject *>(groups.release());
}


PyDoc_STRVAR(doc_group_sparse,
"group_sparse(m, n, indices, indptr)\n"
"\n"
"Group the columns of a CSC structure for sparse finite differencing.\n"
"\n"
"Parameters\n"
"----------\n"
"m, n : int\n"
"    Shape of the structure.\n"
"indices, indptr : array_like of int32 or int64\n"
"    The CSC index arrays; `indptr` has n + 1 entries.\n"
"\n"
"Returns\n"
"-------\n"
"groups : (n,) ndarray of intp\n"
"    The group of each column, from 0 to the number of groups minus 1.\n");


PyObject *py_group_sparse(PyObject *Py_UNUSED(self), PyObject *args)
{
    Py_ssize_t m, n;
    PyObject *indices_obj, *indptr_obj;
    if (!PyArg_ParseTuple(args, "nnOO", &m, &n, &indices_obj, &indptr_obj)) { return NULL; }
    if (m < 0 || n < 0) {
        PyErr_SetString(PyExc_ValueError, "`m` and `n` must be nonnegative.");
        return NULL;
    }
    // int32 index arrays as they are, anything else as int64
    const int typenum = (is_int32_array(indices_obj) && is_int32_array(indptr_obj))
                        ? NPY_INT32 : NPY_INT64;
    OwnedArray indices = as_array(indices_obj, typenum, 1, NPY_ARRAY_IN_ARRAY, "indices");
    if (!indices) { return NULL; }
    OwnedArray indptr = as_array(indptr_obj, typenum, 1, NPY_ARRAY_IN_ARRAY, "indptr");
    if (!indptr) { return NULL; }
    if (PyArray_DIM(indptr.get(), 0) != n + 1) {
        PyErr_SetString(PyExc_ValueError, "`indptr` must have n + 1 entries.");
        return NULL;
    }
    const npy_intp nnz = PyArray_DIM(indices.get(), 0);

    auto run = [&](auto tag) -> PyObject * {
        using I = decltype(tag);
        const I *ind = static_cast<const I *>(PyArray_DATA(indices.get()));
        const I *ptr = static_cast<const I *>(PyArray_DATA(indptr.get()));
        if (!is_valid_csc(ind, nnz, ptr, m, n)) {
            PyErr_SetString(PyExc_ValueError,
                            "Invalid CSC structure: `indptr` must be nondecreasing within "
                            "[0, len(indices)], and `indices` within [0, m).");
            return NULL;
        }
        OwnedArray groups = new_groups(n);
        if (!groups) { return NULL; }
        unsigned char *u = static_cast<unsigned char *>(PyMem_RawMalloc(m > 0 ? m : 1));
        if (u == NULL) { return PyErr_NoMemory(); }
        {
            GILRelease nogil;
            group_sparse_impl(ind, ptr, m, n, static_cast<npy_intp *>(PyArray_DATA(groups.get())), u);
        }
        PyMem_RawFree(u);
        return reinterpret_cast<PyObject *>(groups.release());
    };
    return (typenum == NPY_INT32) ? run(int32_t{}) : run(int64_t{});
}


PyMethodDef group_columns_methods[] = {
    {"group_dense", py_group_dense, METH_VARARGS, doc_group_dense},
    {"group_sparse", py_group_sparse, METH_VARARGS, doc_group_sparse},
    {NULL, NULL, 0, NULL}
};


int group_columns_module_exec(PyObject *Py_UNUSED(module))
{
    if (_import_array() < 0) { return -1; }
    return 0;
}


PyModuleDef_Slot group_columns_module_slots[] = {
    {Py_mod_exec, (void *)group_columns_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};


// Designated initializers require C++20; MSVC /std:c++17 rejects them.
PyModuleDef group_columns_module = {
    PyModuleDef_HEAD_INIT,                /* m_base */
    "_group_columns",                     /* m_name */
    NULL,                                 /* m_doc */
    0,                                    /* m_size */
    group_columns_methods,                /* m_methods */
    group_columns_module_slots,           /* m_slots */
    NULL,                                 /* m_traverse */
    NULL,                                 /* m_clear */
    NULL                                  /* m_free */
};

}  // namespace


PyMODINIT_FUNC
PyInit__group_columns(void)
{
    return PyModuleDef_Init(&group_columns_module);
}
