/*
 * The shift register loop of `scipy.signal.max_len_seq`.
 *
 * The loop runs without the GIL.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

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


// Whether the loop can write into `a`: a 1-D, aligned, C-contiguous, writeable int8 array.
bool is_int8_output(PyArrayObject *a)
{
    return PyArray_NDIM(a) == 1 && PyArray_TYPE(a) == NPY_INT8 && PyArray_ISCARRAY(a);
}


/*
 * Run the shift register for `length` steps. `state` is a ring buffer of `nbits` bits, read from
 * position idx = 0 on: each step emits state[idx] into `seq`, replaces it by its XOR with the
 * bits `taps` positions further on, and advances idx. The taps are within [0, nbits].
 *
 * @return The final idx, the position of the first bit of the next run.
 */
npy_intp shift_register(const npy_intp *taps, const npy_intp n_taps, int8_t *state,
                        const npy_intp nbits, const npy_intp length, int8_t *seq) noexcept
{
    npy_intp idx = 0;
    for (npy_intp i = 0; i < length; i++) {
        int8_t feedback = state[idx];
        seq[i] = feedback;
        for (npy_intp t = 0; t < n_taps; t++) {
            // (taps[t] + idx) % nbits, for a sum below 2*nbits
            npy_intp k = taps[t] + idx;
            if (k >= nbits) { k -= nbits; }
            feedback ^= state[k];
        }
        state[idx] = feedback;
        if (++idx == nbits) { idx = 0; }
    }
    return idx;
}


PyDoc_STRVAR(doc_max_len_seq_inner,
"_max_len_seq_inner(taps, state, nbits, length, seq)\n"
"\n"
"Fast inner loop of `max_len_seq`.\n"
"\n"
"Parameters\n"
"----------\n"
"taps : (T,) array_like of int\n"
"    The taps of the shift register, within [0, nbits].\n"
"state : (nbits,) ndarray of int8\n"
"    The initial state; overwritten.\n"
"nbits : int\n"
"    Number of bits of the shift register.\n"
"length : int\n"
"    Number of bits to generate.\n"
"seq : (N,) ndarray of int8, N >= length\n"
"    Receives the generated bits.\n"
"\n"
"Returns\n"
"-------\n"
"state : (nbits,) ndarray of int8\n"
"    The final state, rolled so that a further call continues the sequence.\n");


PyObject *py_max_len_seq_inner(PyObject *Py_UNUSED(self), PyObject *args)
{
    PyObject *taps_obj;
    PyArrayObject *state, *seq;
    Py_ssize_t nbits, length;
    if (!PyArg_ParseTuple(args, "OO!nnO!", &taps_obj, &PyArray_Type, (PyObject **)&state,
                          &nbits, &length, &PyArray_Type, (PyObject **)&seq)) {
        return NULL;
    }
    if (nbits < 1 || length < 0) {
        PyErr_SetString(PyExc_ValueError, "`nbits` must be positive and `length` nonnegative.");
        return NULL;
    }
    if (!is_int8_output(state) || PyArray_DIM(state, 0) != nbits) {
        PyErr_SetString(PyExc_ValueError, "`state` must be a writeable, C-contiguous int8 "
                        "array of size nbits.");
        return NULL;
    }
    if (!is_int8_output(seq) || PyArray_DIM(seq, 0) < length) {
        PyErr_SetString(PyExc_ValueError, "`seq` must be a writeable, C-contiguous int8 "
                        "array of at least `length` elements.");
        return NULL;
    }
    OwnedArray taps(reinterpret_cast<PyArrayObject *>(
        PyArray_FROMANY(taps_obj, NPY_INTP, 1, 1, NPY_ARRAY_IN_ARRAY)));
    if (!taps) { return NULL; }
    const npy_intp n_taps = PyArray_DIM(taps.get(), 0);
    const npy_intp *tp = static_cast<const npy_intp *>(PyArray_DATA(taps.get()));
    // The loop indexes `state` with the taps.
    for (npy_intp t = 0; t < n_taps; t++) {
        if (tp[t] < 0 || tp[t] > nbits) {
            PyErr_SetString(PyExc_ValueError, "`taps` must be within [0, nbits].");
            return NULL;
        }
    }

    npy_intp dims[1] = {nbits};
    OwnedArray rolled(reinterpret_cast<PyArrayObject *>(PyArray_EMPTY(1, dims, NPY_INT8, 0)));
    if (!rolled) { return NULL; }

    int8_t *st = static_cast<int8_t *>(PyArray_DATA(state));
    int8_t *out = static_cast<int8_t *>(PyArray_DATA(rolled.get()));
    {
        GILRelease nogil;
        const npy_intp idx = shift_register(tp, n_taps, st, nbits, length,
                                            static_cast<int8_t *>(PyArray_DATA(seq)));
        // Roll the state so that, with idx = 0, the next run is in the right place.
        for (npy_intp k = 0; k < nbits; k++) {
            const npy_intp j = idx + k;
            out[k] = st[j < nbits ? j : j - nbits];
        }
    }
    return reinterpret_cast<PyObject *>(rolled.release());
}


PyMethodDef max_len_seq_inner_methods[] = {
    {"_max_len_seq_inner", py_max_len_seq_inner, METH_VARARGS, doc_max_len_seq_inner},
    {NULL, NULL, 0, NULL}
};


int max_len_seq_inner_module_exec(PyObject *Py_UNUSED(module))
{
    if (_import_array() < 0) { return -1; }
    return 0;
}


PyModuleDef_Slot max_len_seq_inner_module_slots[] = {
    {Py_mod_exec, (void *)max_len_seq_inner_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};


// Designated initializers require C++20; MSVC /std:c++17 rejects them.
PyModuleDef max_len_seq_inner_module = {
    PyModuleDef_HEAD_INIT,                /* m_base */
    "_max_len_seq_inner",                 /* m_name */
    NULL,                                 /* m_doc */
    0,                                    /* m_size */
    max_len_seq_inner_methods,            /* m_methods */
    max_len_seq_inner_module_slots,       /* m_slots */
    NULL,                                 /* m_traverse */
    NULL,                                 /* m_clear */
    NULL                                  /* m_free */
};

}  // namespace


PyMODINIT_FUNC
PyInit__max_len_seq_inner(void)
{
    return PyModuleDef_Init(&max_len_seq_inner_module);
}
