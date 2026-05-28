// Copyright (c) 2026 Gideon Genadi Kogan
// Distributed under the SciPy BSD-3-Clause license; see LICENSE.txt.
//
// 1-D Hampel filter (MAD-based outlier detection) for scipy.ndimage.
//
// The local median is supplied by the caller (typically computed by
// scipy.ndimage.median_filter, which itself uses the O(log W) 1-D rank
// filter), so this module only contributes the MAD computation and outlier
// replacement in O(log W) per sample, using two indexed heaps with an
// equilibrium loop.
//
// Memory footprint per kernel invocation is O(W) (not O(N + W)): the heaps
// and per-slot index maps are W-sized circular structures keyed by
// window-slot in [0, W), and the boundary-extended signal is read on the
// fly rather than materialised.

#include "Python.h"
#include "numpy/arrayobject.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

// Per-dtype sentinel values used as +/- "infinity" in the equilibrium loop
// and as the value reported when peeking into an empty heap. For floats we
// use IEEE infinities; for int64 we fall back to numeric_limits::max/lowest,
// which are safe because the loop only compares them, never arithmetics on
// them.
template <typename T>
struct Sentinel {
    static inline T pos() { return std::numeric_limits<T>::max(); }
    static inline T neg() { return std::numeric_limits<T>::lowest(); }
};
template <>
struct Sentinel<float> {
    static inline float pos() { return  std::numeric_limits<float>::infinity(); }
    static inline float neg() { return -std::numeric_limits<float>::infinity(); }
};
template <>
struct Sentinel<double> {
    static inline double pos() { return  std::numeric_limits<double>::infinity(); }
    static inline double neg() { return -std::numeric_limits<double>::infinity(); }
};

// |a - b| as a non-negative double. For int64 a direct subtraction may
// overflow; promoting to double first avoids UB and preserves enough
// precision for the subsequent threshold comparison.
template <typename T>
inline double abs_diff(T a, T b) {
    return std::fabs((double)a - (double)b);
}

// Signed distance (a - b) computed in double to avoid int64 overflow.
template <typename T>
inline double signed_dist(T a, T b) {
    return (double)a - (double)b;
}

// Read padded signal element at absolute padded position j in [0, N + 2Z).
// Position j corresponds to logical signal index (j - Z). Boundary handling
// matches scipy.ndimage's mode codes exactly (same conventions used by the
// caller's median_filter() so that the supplied median trace lines up).
template <typename T>
static inline T read_padded(const T *signal, npy_intp N, npy_intp j, int Z,
                            int mode, T cval) {
    npy_intp idx = j - Z;
    if (idx >= 0 && idx < N) return signal[idx];
    if (idx < 0) {
        npy_intp k;
        switch (mode) {
        case REFLECT:
            // ... d c b a | a b c d | d c b a ...
            k = -idx - 1;
            return signal[(k < N) ? k : (N - 1)];
        case MIRROR:
            // ... d c b | a b c d | c b a ...
            k = -idx;
            return signal[(k < N) ? k : (N - 1)];
        case NEAREST:
            return signal[0];
        case WRAP:
            return signal[((idx % N) + N) % N];
        case CONSTANT:
            return cval;
        default:
            return signal[0];
        }
    } else {
        npy_intp r = idx - N;
        npy_intp k;
        switch (mode) {
        case REFLECT:
            k = N - 1 - r;
            return signal[(k >= 0) ? k : 0];
        case MIRROR:
            k = N - 2 - r;
            return signal[(k >= 0) ? k : 0];
        case NEAREST:
            return signal[N - 1];
        case WRAP:
            return signal[idx % N];
        case CONSTANT:
            return cval;
        default:
            return signal[N - 1];
        }
    }
}

// Indexed binary heap with O(log W) eager deletion and O(log W) in-place
// update by external window-slot id in [0, W).
template <typename T>
class IndexedHeap {
public:
    struct Node { T value; int slot; };

    std::vector<Node> heap;
    std::vector<int>& pos_map;   // pos_map[slot] = heap index, or -1
    bool is_max_heap;

    IndexedHeap(std::vector<int>& p_map, bool max_h)
        : pos_map(p_map), is_max_heap(max_h) {}

    inline bool compare(T a, T b) const {
        return is_max_heap ? (a > b) : (a < b);
    }

    inline void swap_nodes(int i, int j) {
        std::swap(heap[i], heap[j]);
        pos_map[heap[i].slot] = i;
        pos_map[heap[j].slot] = j;
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

    void push(T val, int slot) {
        heap.push_back({val, slot});
        int idx = (int)heap.size() - 1;
        pos_map[slot] = idx;
        shift_up(idx);
    }

    void erase(int slot) {
        int idx = pos_map[slot];
        if (idx < 0) return;
        int n = (int)heap.size();
        if (idx == n - 1) {
            heap.pop_back();
            pos_map[slot] = -1;
            return;
        }
        swap_nodes(idx, n - 1);
        heap.pop_back();
        pos_map[slot] = -1;
        if (idx < (int)heap.size()) {
            shift_up(idx);
            shift_down(idx);
        }
    }

    // Overwrite the value held at `slot` and re-establish heap order. The
    // caller guarantees `slot` is currently in this heap. The value moves
    // either up or down (not both): if the new value is "better" for this
    // heap (larger for a max-heap, smaller for a min-heap) it sifts up,
    // otherwise it sifts down.
    void update(int slot, T new_val) {
        int idx = pos_map[slot];
        T old_val = heap[idx].value;
        heap[idx].value = new_val;
        if (compare(new_val, old_val)) {
            shift_up(idx);
        } else {
            shift_down(idx);
        }
    }

    inline T peek() const {
        if (!heap.empty()) return heap[0].value;
        return is_max_heap ? Sentinel<T>::neg() : Sentinel<T>::pos();
    }

    inline int size() const { return (int)heap.size(); }
};

// Tracks the kth-smallest value of a sliding W-element multiset using two
// W-slot indexed heaps. Supports O(log W) target-rank changes and O(log W)
// in-place value updates at a known slot.
template <typename T>
class RankTracker {
public:
    std::vector<int> pos_max;   // sized W
    std::vector<int> pos_min;   // sized W
    IndexedHeap<T> max_heap;
    IndexedHeap<T> min_heap;
    int target_k;

    RankTracker(int W, int target_k_init)
        : pos_max(W, -1), pos_min(W, -1),
          max_heap(pos_max, true), min_heap(pos_min, false),
          target_k(target_k_init) {}

    void insert(T val, int slot) {
        if (max_heap.size() == 0 || val <= max_heap.peek()) {
            max_heap.push(val, slot);
        } else {
            min_heap.push(val, slot);
        }
        balance();
    }

    // Replace the value at `slot` (currently in one of the two heaps) with
    // `new_val`. Sizes are unchanged; only values shift. A single
    // boundary-swap of the heap roots restores the cross-heap invariant in
    // the worst case (only one slot's value changed, so at most one
    // element can cross the boundary).
    void replace(int slot, T new_val) {
        if (pos_max[slot] >= 0) {
            max_heap.update(slot, new_val);
        } else {
            min_heap.update(slot, new_val);
        }
        if (max_heap.size() > 0 && min_heap.size() > 0 &&
            max_heap.peek() > min_heap.peek()) {
            swap_tops();
        }
    }

    void set_rank(int new_k) {
        target_k = new_k;
        balance();
    }

    void balance() {
        while (max_heap.size() > target_k && max_heap.size() > 0) {
            auto node = max_heap.heap[0];
            max_heap.erase(node.slot);
            min_heap.push(node.value, node.slot);
        }
        while (max_heap.size() < target_k && min_heap.size() > 0) {
            auto node = min_heap.heap[0];
            min_heap.erase(node.slot);
            max_heap.push(node.value, node.slot);
        }
    }

    inline T get_kth_value() const { return max_heap.peek(); }
    inline T get_kplus1_value() const { return min_heap.peek(); }

private:
    // Swap the roots of max_heap and min_heap (one element migrates each
    // way) and re-heap from the new tops. O(log W).
    void swap_tops() {
        int s_max = max_heap.heap[0].slot;
        int s_min = min_heap.heap[0].slot;
        T v_max = max_heap.heap[0].value;
        T v_min = min_heap.heap[0].value;
        max_heap.heap[0] = {v_min, s_min};
        min_heap.heap[0] = {v_max, s_max};
        pos_max[s_max] = -1;
        pos_max[s_min] = 0;
        pos_min[s_max] = 0;
        pos_min[s_min] = -1;
        max_heap.shift_down(0);
        min_heap.shift_down(0);
    }
};

// Core Hampel kernel. Returns 0 on success, -1 on memory failure.
//
// For every position i in [0, N) the local median M = median_in[i] is
// already known. We slide a window of width W centred on i (over the
// boundary-extended signal, read on the fly) and use the equilibrium loop
// to extract MAD(i) in O(log W). If |signal[i] - M| > threshold * MAD(i)
// the sample is flagged as an outlier and replaced by M in `filtered`;
// otherwise the original value is kept.
template <typename T>
int _hampel_1d(const T *signal, const T *median_in, npy_intp N, int win_len,
               double threshold, int mode, T cval,
               T *filtered, npy_bool *changed)
{
    if (win_len <= 0) return -1;
    int W = win_len;
    int Z = W / 2;

    // Trivial pass-through for degenerate windows or empty input.
    if (W <= 1 || N == 0) {
        for (npy_intp i = 0; i < N; ++i) {
            filtered[i] = signal[i];
            changed[i] = 0;
        }
        return 0;
    }

    // Equilibrium-loop parameters: 1-based internal "median rank".
    int M_rank = Z + 1;
    int k_L = Z / 2;
    int k_U = Z - k_L;

    try {
        RankTracker<T> upper_tracker(W, M_rank + k_U);
        RankTracker<T> lower_tracker(W, M_rank - k_L - 1);

        // Initial window covers padded positions [0, W). Each one lives at
        // slot (j % W) == j.
        for (int j = 0; j < W; ++j) {
            T val = read_padded(signal, N, (npy_intp)j, Z, mode, cval);
            upper_tracker.insert(val, j);
            lower_tracker.insert(val, j);
        }

        const double D_NEG_INF = -std::numeric_limits<double>::infinity();
        const double D_POS_INF =  std::numeric_limits<double>::infinity();

        // Window i covers padded[i .. i+W-1]; its centre maps to signal[i].
        // Per slide (i >= 1) one value is replaced in-place at slot
        // ((i - 1) % W) == ((i + W - 1) % W): padded[i-1] leaves,
        // padded[i+W-1] enters.
        for (npy_intp i = 0; i < N; ++i) {
            if (i > 0) {
                int slot = (int)((i - 1) % W);
                T new_val = read_padded(signal, N, i + W - 1, Z, mode, cval);
                upper_tracker.replace(slot, new_val);
                lower_tracker.replace(slot, new_val);
            }

            T M = median_in[i];

            // Re-balance k_U / k_L until the (k_U)th value above M and the
            // (k_L)th value below M jointly bracket the same rank distance.
            // All distances are computed in double to avoid int64 overflow.
            while (true) {
                T X_upper_k = upper_tracker.get_kth_value();
                T X_upper_next = upper_tracker.get_kplus1_value();

                T X_lower_k = lower_tracker.get_kplus1_value();
                T X_lower_prev = lower_tracker.get_kth_value();

                double A_km1 = (k_U > 0) ? signed_dist(X_upper_k, M) : D_NEG_INF;
                double A_k   = ((M_rank + k_U) < W) ? signed_dist(X_upper_next, M) : D_POS_INF;
                double B_km1 = (k_L > 0) ? signed_dist(M, X_lower_k) : D_NEG_INF;
                double B_k   = ((M_rank - k_L - 1) > 0) ? signed_dist(M, X_lower_prev) : D_POS_INF;

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

            double cand_A = (k_U > 0) ? signed_dist(upper_tracker.get_kth_value(), M) : 0.0;
            double cand_B = (k_L > 0) ? signed_dist(M, lower_tracker.get_kplus1_value()) : 0.0;
            double mad = std::max({cand_A, cand_B, 0.0});

            double abs_dev = abs_diff(signal[i], M);
            // All operands are double, avoiding int64 overflow on
            // threshold * mad and giving a dtype-agnostic comparison.
            if (abs_dev > threshold * mad) {
                filtered[i] = M;
                changed[i] = 1;
            } else {
                filtered[i] = signal[i];
                changed[i] = 0;
            }
        }
    } catch (std::bad_alloc&) {
        return -1;
    }
    return 0;
}

// Python entry point:
//     hampel(signal, median, win_len, threshold, mode, cval,
//            filtered_out, changed_out)
// - signal, median, filtered_out: same length, same dtype
//   (float32 / float64 / int64).
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
        PyArray_DiscardWritebackIfCopy(filt);
        PyArray_DiscardWritebackIfCopy(chg);
        Py_DECREF(signal); Py_DECREF(median);
        Py_DECREF(filt); Py_DECREF(chg);
        PyErr_SetString(PyExc_TypeError,
                        "signal, median and filtered arrays must share dtype");
        return NULL;
    }

    npy_intp N = PyArray_SIZE(signal);
    if (PyArray_SIZE(median) != N || PyArray_SIZE(filt) != N ||
        PyArray_SIZE(chg) != N) {
        PyArray_DiscardWritebackIfCopy(filt);
        PyArray_DiscardWritebackIfCopy(chg);
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
        double cval_d = PyFloat_AsDouble(cval_obj);
        if (cval_d == -1.0 && PyErr_Occurred()) {
            PyArray_DiscardWritebackIfCopy(filt);
            PyArray_DiscardWritebackIfCopy(chg);
            Py_DECREF(signal); Py_DECREF(median);
            Py_DECREF(filt); Py_DECREF(chg);
            return NULL;
        }
        float cval = (float)cval_d;
        Py_BEGIN_ALLOW_THREADS
        status = _hampel_1d<float>(s, m, N, win_len, threshold,
                                   mode, cval, f, c_chg);
        Py_END_ALLOW_THREADS
        break;
    }
    case NPY_DOUBLE: {
        double *s = (double *)PyArray_DATA(signal);
        double *m = (double *)PyArray_DATA(median);
        double *f = (double *)PyArray_DATA(filt);
        double cval = PyFloat_AsDouble(cval_obj);
        if (cval == -1.0 && PyErr_Occurred()) {
            PyArray_DiscardWritebackIfCopy(filt);
            PyArray_DiscardWritebackIfCopy(chg);
            Py_DECREF(signal); Py_DECREF(median);
            Py_DECREF(filt); Py_DECREF(chg);
            return NULL;
        }
        Py_BEGIN_ALLOW_THREADS
        status = _hampel_1d<double>(s, m, N, win_len, threshold,
                                    mode, cval, f, c_chg);
        Py_END_ALLOW_THREADS
        break;
    }
    case NPY_INT64: {
        npy_int64 *s = (npy_int64 *)PyArray_DATA(signal);
        npy_int64 *m = (npy_int64 *)PyArray_DATA(median);
        npy_int64 *f = (npy_int64 *)PyArray_DATA(filt);
        npy_int64 cval = (npy_int64)PyLong_AsLongLong(cval_obj);
        if (cval == -1 && PyErr_Occurred()) {
            PyArray_DiscardWritebackIfCopy(filt);
            PyArray_DiscardWritebackIfCopy(chg);
            Py_DECREF(signal); Py_DECREF(median);
            Py_DECREF(filt); Py_DECREF(chg);
            return NULL;
        }
        Py_BEGIN_ALLOW_THREADS
        status = _hampel_1d<npy_int64>(s, m, N, win_len, threshold,
                                       mode, cval, f, c_chg);
        Py_END_ALLOW_THREADS
        break;
    }
    default:
        PyErr_SetString(PyExc_TypeError,
                        "hampel filter only supports float32, float64, "
                        "and int64");
        break;
    }

    if (status == -1 && !PyErr_Occurred()) {
        PyErr_SetString(PyExc_MemoryError,
                        "failed to allocate memory for hampel filter");
    }

    if (PyErr_Occurred()) {
        PyArray_DiscardWritebackIfCopy(filt);
        PyArray_DiscardWritebackIfCopy(chg);
    } else {
        PyArray_ResolveWritebackIfCopy(filt);
        PyArray_ResolveWritebackIfCopy(chg);
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
