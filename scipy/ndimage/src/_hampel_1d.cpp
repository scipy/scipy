// Copyright (c) 2026 Gideon Genadi Kogan
// Distributed under the SciPy BSD-3-Clause license; see LICENSE.txt.
//
// 1-D Hampel filter (MAD-based outlier detection) for scipy.ndimage.
//
// The local median is supplied by the caller (typically computed by
// scipy.ndimage.median_filter, which itself uses the O(log W) 1-D rank
// filter), so this module only contributes the MAD computation and outlier
// replacement, in O(log W) per sample, using two indexed heaps with an
// equilibrium loop.

#include "Python.h"
#include "numpy/arrayobject.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <vector>

namespace {

typedef enum {
    NEAREST = 0,
    WRAP = 1,
    REFLECT = 2,
    MIRROR = 3,
    CONSTANT = 4,
} Mode;

// Indexed binary heap with O(log W) eager deletion by external index.
template <typename T>
class IndexedHeap {
public:
    struct Node { T value; int w_idx; };

    std::vector<Node> heap;
    std::vector<int>& pos_map;
    bool is_max_heap;

    IndexedHeap(std::vector<int>& p_map, bool max_h)
        : pos_map(p_map), is_max_heap(max_h) {}

    inline bool compare(T a, T b) const {
        return is_max_heap ? (a > b) : (a < b);
    }

    void swap_nodes(int i, int j) {
        std::swap(heap[i], heap[j]);
        pos_map[heap[i].w_idx] = i;
        pos_map[heap[j].w_idx] = j;
    }

    void shift_up(int i) {
        while (i > 0) {
            int p = (i - 1) / 2;
            if (compare(heap[i].value, heap[p].value)) {
                swap_nodes(i, p);
                i = p;
            } else {
                break;
            }
        }
    }

    void shift_down(int i) {
        int n = (int)heap.size();
        while (2 * i + 1 < n) {
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            int best = left;
            if (right < n && compare(heap[right].value, heap[left].value)) {
                best = right;
            }
            if (compare(heap[best].value, heap[i].value)) {
                swap_nodes(i, best);
                i = best;
            } else {
                break;
            }
        }
    }

    void push(T val, int w_idx) {
        heap.push_back({val, w_idx});
        int idx = (int)heap.size() - 1;
        pos_map[w_idx] = idx;
        shift_up(idx);
    }

    void erase(int w_idx) {
        if (w_idx < 0 || w_idx >= (int)pos_map.size() || pos_map[w_idx] == -1) {
            return;
        }
        int idx = pos_map[w_idx];
        int n = (int)heap.size();
        if (idx == n - 1) {
            heap.pop_back();
            pos_map[w_idx] = -1;
            return;
        }
        swap_nodes(idx, n - 1);
        heap.pop_back();
        pos_map[w_idx] = -1;
        if (idx < (int)heap.size()) {
            shift_up(idx);
            shift_down(idx);
        }
    }

    T peek() const {
        if (!heap.empty()) return heap[0].value;
        return is_max_heap ? -std::numeric_limits<T>::infinity()
                           :  std::numeric_limits<T>::infinity();
    }

    int size() const { return (int)heap.size(); }
};

// Tracks the kth-smallest value of a sliding set using two heaps, with
// O(log W) target-rank changes.
template <typename T>
class RankTracker {
public:
    std::vector<int> pos_max;
    std::vector<int> pos_min;
    IndexedHeap<T> max_heap;
    IndexedHeap<T> min_heap;
    int target_k;

    RankTracker(int universe_size, int target_k_init)
        : pos_max(universe_size, -1), pos_min(universe_size, -1),
          max_heap(pos_max, true), min_heap(pos_min, false),
          target_k(target_k_init) {}

    void insert(T val, int w_idx) {
        if (max_heap.size() == 0 || val <= max_heap.peek()) {
            max_heap.push(val, w_idx);
        } else {
            min_heap.push(val, w_idx);
        }
        balance();
    }

    void remove(int w_idx) {
        if (w_idx >= 0 && w_idx < (int)pos_max.size() && pos_max[w_idx] != -1) {
            max_heap.erase(w_idx);
        } else if (w_idx >= 0 && w_idx < (int)pos_min.size() && pos_min[w_idx] != -1) {
            min_heap.erase(w_idx);
        }
        balance();
    }

    void set_rank(int new_k) {
        target_k = new_k;
        balance();
    }

    void balance() {
        while (max_heap.size() > target_k && max_heap.size() > 0) {
            auto node = max_heap.heap[0];
            max_heap.erase(node.w_idx);
            min_heap.push(node.value, node.w_idx);
        }
        while (max_heap.size() < target_k && min_heap.size() > 0) {
            auto node = min_heap.heap[0];
            min_heap.erase(node.w_idx);
            max_heap.push(node.value, node.w_idx);
        }
    }

    T get_kth_value() const { return max_heap.peek(); }
    T get_kplus1_value() const { return min_heap.peek(); }
};

// Build a padded copy of the signal of length N + 2*Z so that the centred
// window of width W = 2*Z (+1) fits at every output position. Boundary
// extension follows scipy.ndimage's mode codes, matching what the caller's
// median_filter() used to produce the input median signal.
template <typename T>
void pad_signal(const T *signal, npy_intp N, int Z, T *padded,
                int mode, T cval)
{
    // Centre copy.
    for (npy_intp i = 0; i < N; ++i) {
        padded[Z + i] = signal[i];
    }
    if (N == 0) {
        for (int i = 0; i < 2 * Z; ++i) {
            padded[i] = (mode == CONSTANT) ? cval : (T)0;
        }
        return;
    }
    for (int i = 0; i < Z; ++i) {
        npy_intp li, ri;
        switch (mode) {
        case REFLECT:
            // ... d c b a | a b c d | d c b a ...
            li = (i < N) ? i : (N - 1);
            ri = (N - 1 - i >= 0) ? (N - 1 - i) : 0;
            padded[Z - 1 - i] = signal[li];
            padded[Z + N + i] = signal[ri];
            break;
        case MIRROR:
            // ... d c b | a b c d | c b a ...
            li = (i + 1 < N) ? (i + 1) : (N - 1);
            ri = (N - 2 - i >= 0) ? (N - 2 - i) : 0;
            padded[Z - 1 - i] = signal[li];
            padded[Z + N + i] = signal[ri];
            break;
        case NEAREST:
            padded[Z - 1 - i] = signal[0];
            padded[Z + N + i] = signal[N - 1];
            break;
        case WRAP: {
            // periodic
            npy_intp lp = ((-(i + 1)) % N + N) % N;
            npy_intp rp = i % N;
            padded[Z - 1 - i] = signal[lp];
            padded[Z + N + i] = signal[rp];
            break;
        }
        case CONSTANT:
            padded[Z - 1 - i] = cval;
            padded[Z + N + i] = cval;
            break;
        default:
            padded[Z - 1 - i] = signal[0];
            padded[Z + N + i] = signal[N - 1];
            break;
        }
    }
}

// Core Hampel kernel. Returns 0 on success, -1 on memory failure.
//
// For every position i in [0, N) the local median M = median_in[i] is
// already known. We slide a window of width W centred on i (over a
// boundary-padded copy of the signal) and use the equilibrium loop from the
// attached algorithm to extract MAD(i) in O(log W). If
// |signal[i] - M| > threshold * MAD(i) the sample is flagged as an outlier
// and replaced by M in `filtered`; otherwise the original value is kept.
template <typename T>
int _hampel_1d(const T *signal, const T *median_in, npy_intp N, int win_len,
               double threshold, int mode, T cval,
               T *filtered, npy_bool *changed)
{
    if (win_len <= 0) return -1;
    int W = win_len;
    int Z = W / 2;

    // Trivial pass-through for degenerate windows.
    if (W <= 1) {
        for (npy_intp i = 0; i < N; ++i) {
            filtered[i] = signal[i];
            changed[i] = 0;
        }
        return 0;
    }

    npy_intp padded_N = N + 2 * Z;
    std::unique_ptr<T[]> padded;
    try {
        padded = std::unique_ptr<T[]>(new T[padded_N]);
    } catch (std::bad_alloc&) {
        return -1;
    }
    pad_signal(signal, N, Z, padded.get(), mode, cval);

    // Equilibrium-loop parameters following the attached design.
    // Internal "median rank" in 1-based form: M_rank = Z + 1.
    int M_rank = Z + 1;
    int k_L = Z / 2;
    int k_U = Z - k_L;

    RankTracker<T> upper_tracker((int)padded_N, M_rank + k_U);
    RankTracker<T> lower_tracker((int)padded_N, M_rank - k_L - 1);

    try {
        // Initial window: padded[0 .. W-1].
        for (int j = 0; j < W; ++j) {
            upper_tracker.insert(padded[j], j);
            lower_tracker.insert(padded[j], j);
        }
    } catch (std::bad_alloc&) {
        return -1;
    }

    const T INF = std::numeric_limits<T>::infinity();

    // Window i covers padded[i .. i+W-1]; its centre maps to signal[i].
    for (npy_intp i = 0; i < N; ++i) {
        if (i > 0) {
            int old_idx = (int)(i - 1);
            int new_idx = (int)(i + W - 1);
            T new_val = padded[new_idx];
            upper_tracker.remove(old_idx);
            lower_tracker.remove(old_idx);
            try {
                upper_tracker.insert(new_val, new_idx);
                lower_tracker.insert(new_val, new_idx);
            } catch (std::bad_alloc&) {
                return -1;
            }
        }

        T M = median_in[i];

        // Re-balance k_U / k_L until the (k_U)th value above M and the
        // (k_L)th value below M jointly bracket the same rank distance.
        while (true) {
            T X_upper_k = upper_tracker.get_kth_value();
            T X_upper_next = upper_tracker.get_kplus1_value();

            T X_lower_k = lower_tracker.get_kplus1_value();
            T X_lower_prev = lower_tracker.get_kth_value();

            T A_km1 = (k_U > 0) ? (X_upper_k - M) : -INF;
            T A_k   = ((M_rank + k_U) < W) ? (X_upper_next - M) : INF;
            T B_km1 = (k_L > 0) ? (M - X_lower_k) : -INF;
            T B_k   = ((M_rank - k_L - 1) > 0) ? (M - X_lower_prev) : INF;

            if (k_U > 0 && A_km1 > B_k) {
                k_U--;
                k_L++;
                upper_tracker.set_rank(M_rank + k_U);
                lower_tracker.set_rank(M_rank - k_L - 1);
            } else if (k_L > 0 && B_km1 > A_k) {
                k_L--;
                k_U++;
                upper_tracker.set_rank(M_rank + k_U);
                lower_tracker.set_rank(M_rank - k_L - 1);
            } else {
                break;
            }
        }

        T cand_A = (k_U > 0) ? (upper_tracker.get_kth_value() - M) : (T)0;
        T cand_B = (k_L > 0) ? (M - lower_tracker.get_kplus1_value()) : (T)0;
        T mad = std::max({cand_A, cand_B, (T)0});

        T abs_dev = std::fabs(signal[i] - M);
        if (abs_dev > (T)threshold * mad) {
            filtered[i] = M;
            changed[i] = 1;
        } else {
            filtered[i] = signal[i];
            changed[i] = 0;
        }
    }
    return 0;
}

// Python entry point:
//     hampel(signal, median, win_len, threshold, mode, cval,
//            filtered_out, changed_out)
// - signal, median, filtered_out: same length, same dtype (float32/float64).
// - changed_out: bool array, same length.
static PyObject *hampel(PyObject *self, PyObject *args)
{
    PyObject *sig_obj, *med_obj, *cval_obj, *out_obj, *chg_obj;
    int win_len, mode;
    double threshold;
    int status = 0;

    if (!PyArg_ParseTuple(args, "OOidiOOO",
                          &sig_obj, &med_obj, &win_len, &threshold,
                          &mode, &cval_obj, &out_obj, &chg_obj)) {
        return NULL;
    }

    PyArrayObject *signal = (PyArrayObject *)PyArray_FROM_OTF(
        sig_obj, NPY_NOTYPE, NPY_ARRAY_IN_ARRAY);
    if (signal == NULL) return NULL;
    PyArrayObject *median = (PyArrayObject *)PyArray_FROM_OTF(
        med_obj, NPY_NOTYPE, NPY_ARRAY_IN_ARRAY);
    if (median == NULL) { Py_DECREF(signal); return NULL; }
    PyArrayObject *filt = (PyArrayObject *)PyArray_FROM_OTF(
        out_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);
    if (filt == NULL) { Py_DECREF(signal); Py_DECREF(median); return NULL; }
    PyArrayObject *chg = (PyArrayObject *)PyArray_FROM_OTF(
        chg_obj, NPY_BOOL, NPY_ARRAY_INOUT_ARRAY2);
    if (chg == NULL) {
        Py_DECREF(signal); Py_DECREF(median); Py_DECREF(filt);
        return NULL;
    }

    if (PyArray_TYPE(signal) != PyArray_TYPE(median) ||
        PyArray_TYPE(signal) != PyArray_TYPE(filt)) {
        Py_DECREF(signal); Py_DECREF(median);
        Py_DECREF(filt); Py_DECREF(chg);
        PyErr_SetString(PyExc_TypeError,
                        "signal, median and filtered arrays must share dtype");
        return NULL;
    }

    npy_intp N = PyArray_SIZE(signal);
    if (PyArray_SIZE(median) != N || PyArray_SIZE(filt) != N ||
        PyArray_SIZE(chg) != N) {
        Py_DECREF(signal); Py_DECREF(median);
        Py_DECREF(filt); Py_DECREF(chg);
        PyErr_SetString(PyExc_ValueError,
                        "input arrays must all have the same size");
        return NULL;
    }
    int type = PyArray_TYPE(signal);
    npy_bool *c_chg = (npy_bool *)PyArray_DATA(chg);

    switch (type) {
    case NPY_FLOAT: {
        float *s = (float *)PyArray_DATA(signal);
        float *m = (float *)PyArray_DATA(median);
        float *f = (float *)PyArray_DATA(filt);
        float cval = (float)PyFloat_AsDouble(cval_obj);
        status = _hampel_1d<float>(s, m, N, win_len, threshold,
                                   mode, cval, f, c_chg);
        break;
    }
    case NPY_DOUBLE: {
        double *s = (double *)PyArray_DATA(signal);
        double *m = (double *)PyArray_DATA(median);
        double *f = (double *)PyArray_DATA(filt);
        double cval = PyFloat_AsDouble(cval_obj);
        status = _hampel_1d<double>(s, m, N, win_len, threshold,
                                    mode, cval, f, c_chg);
        break;
    }
    default:
        PyErr_SetString(PyExc_TypeError,
                        "hampel filter only supports float32 and float64");
        break;
    }

    if (status == -1 && !PyErr_Occurred()) {
        PyErr_SetString(PyExc_MemoryError,
                        "failed to allocate memory for hampel filter");
    }

    Py_DECREF(signal); Py_DECREF(median);
    Py_DECREF(filt); Py_DECREF(chg);
    if (PyErr_Occurred()) return NULL;
    Py_RETURN_NONE;
}

static PyMethodDef hampel_methods[] = {
    {"hampel", hampel, METH_VARARGS,
     "1D Hampel filter (MAD-based outlier detection)."},
    {NULL, NULL, 0, NULL}
};

static int
_hampel_1d_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }
    if (PyModule_AddStringConstant(module, "__author__",
                                   "Gideon Genadi Kogan") < 0) {
        return -1;
    }
    return 0;
}

static PyModuleDef_Slot _hampel_1d_slots[] = {
    {Py_mod_exec, (void*)_hampel_1d_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};

static struct PyModuleDef _hampel_1d_module = {
    /* m_base     */ PyModuleDef_HEAD_INIT,
    /* m_name     */ "_hampel_1d",
    /* m_doc      */ "1-D Hampel filter (MAD-based outlier detection).\n"
                     "Author: Gideon Genadi Kogan, 2026.",
    /* m_size     */ 0,
    /* m_methods  */ hampel_methods,
    /* m_slots    */ _hampel_1d_slots,
    /* m_traverse */ NULL,
    /* m_clear    */ NULL,
    /* m_free     */ NULL
};

} // anonymous namespace

PyMODINIT_FUNC
PyInit__hampel_1d(void)
{
    return PyModuleDef_Init(&_hampel_1d_module);
}
