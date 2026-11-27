/*
 * The exact two-sided p-value of the two-sample Kolmogorov-Smirnov test for samples of
 * different sizes; `_compute_outer_prob_inside_method` of `scipy.stats._stats_py`.
 *
 * The recursion runs without the GIL.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <new>
#include <utility>
#include <vector>


namespace {

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
 * The recursion of `_compute_outer_prob_inside_method`, for m >= n, with mg = m/g, ng = n/g and
 * the window [0, maxj) of the first column. `A` holds lenA ones on entry.
 *
 * @return false if a window would leave `A`, which the choice of lenA rules out.
 */
bool outer_prob_inside(const int64_t m, const int64_t n, const int64_t mg, const int64_t ng,
                       const int64_t h, int64_t maxj, double *A, const int64_t lenA,
                       double *result) noexcept
{
    int64_t minj = 0;
    int64_t curlen = maxj - minj;
    // Initialize the first column
    std::fill(A + minj, A + maxj, 0.0);
    for (int64_t i = 1; i <= m; i++) {
        // Generate the next column. First calculate the sliding window
        const int64_t lastminj = minj, lastlen = curlen;
        minj = std::max<int64_t>(
            static_cast<int64_t>(std::floor(static_cast<double>(ng * i - h) / mg)) + 1, 0);
        minj = std::min(minj, n);
        maxj = std::min<int64_t>(
            static_cast<int64_t>(std::ceil(static_cast<double>(ng * i + h) / mg)), n + 1);
        if (maxj <= minj) {
            *result = 1.0;
            return true;
        }
        // The loop below reads A up to index maxj - 1 - lastminj.
        if (maxj - lastminj > lenA) { return false; }

        // Now fill in the values. We cannot use cumsum, unfortunately.
        double val = (minj == 0) ? 0.0 : 1.0;
        for (int64_t jj = 0; jj < maxj - minj; jj++) {
            const int64_t j = jj + minj;
            val = (A[jj + minj - lastminj] * i + val * j) / (i + j);
            A[jj] = val;
        }
        curlen = maxj - minj;
        if (lastlen > curlen) {
            // Set some carried-over elements to 1
            std::fill(A + curlen, A + std::min(lastlen, lenA), 1.0);
        }
    }
    *result = A[maxj - minj - 1];
    return true;
}


PyDoc_STRVAR(doc_compute_outer_prob_inside_method,
"_compute_outer_prob_inside_method(m, n, g, h)\n"
"\n"
"Count the proportion of paths that do not stay strictly inside two\n"
"diagonal lines.\n"
"\n"
"Parameters\n"
"----------\n"
"m : int\n"
"    m > 0\n"
"n : int\n"
"    n > 0\n"
"g : int\n"
"    g is greatest common divisor of m and n\n"
"h : int\n"
"    0 <= h <= lcm(m,n)\n"
"\n"
"Returns\n"
"-------\n"
"p : float\n"
"    The proportion of paths that do not stay inside the two lines.\n"
"\n"
"The classical algorithm counts the integer lattice paths from (0, 0)\n"
"to (m, n) which satisfy |x/m - y/n| < h / lcm(m, n).\n"
"The paths make steps of size +1 in either positive x or positive y\n"
"directions.\n"
"We are, however, interested in 1 - proportion to computes p-values,\n"
"so we change the recursion to compute 1 - p directly while staying\n"
"within the \"inside method\" a described by Hodges.\n"
"\n"
"We generally follow Hodges' treatment of Drion/Gnedenko/Korolyuk.\n"
"Hodges, J.L. Jr.,\n"
"\"The Significance Probability of the Smirnov Two-Sample Test,\"\n"
"Arkiv fiur Matematik, 3, No. 43 (1958), 469-86.\n"
"\n"
"For the recursion for 1-p see\n"
"Viehmann, T.: \"Numerically more stable computation of the p-values\n"
"for the two-sample Kolmogorov-Smirnov test,\" arXiv: 2102.08037\n");


PyObject *py_compute_outer_prob_inside_method(PyObject *Py_UNUSED(self), PyObject *args)
{
    long long m_in, n_in, g_in, h_in;
    if (!PyArg_ParseTuple(args, "LLLL", &m_in, &n_in, &g_in, &h_in)) { return NULL; }
    int64_t m = m_in, n = n_in;
    const int64_t g = g_in, h = h_in;
    if (m <= 0 || n <= 0 || g <= 0 || h < 0 || m % g != 0 || n % g != 0) {
        PyErr_SetString(PyExc_ValueError, "Expected m > 0, n > 0, g > 0 a common divisor of "
                        "m and n, and h >= 0.");
        return NULL;
    }

    // Probability is symmetrical in m, n. Computation below uses m >= n.
    if (m < n) { std::swap(m, n); }
    const int64_t mg = m / g, ng = n / g;

    // Count the integer lattice paths from (0, 0) to (m, n) which satisfy
    // |nx/g - my/g| < h.
    // Compute matrix A such that:
    //  A(x, 0) = A(0, y) = 1
    //  A(x, y) = A(x, y-1) + A(x-1, y), for x,y>=1, except that
    //  A(x, y) = 0 if |x/m - y/n|>= h
    // Probability is A(m, n)/binom(m+n, n)
    // Only need to preserve a single column of A, and only a sliding window of it.
    // This is an integer calculation, but the entries are essentially binomial
    // coefficients, hence grow quickly. Instead of the counts, the recursion carries
    // the proportions, which stay within [0, 1].
    const int64_t maxj = std::min<int64_t>(
        static_cast<int64_t>(std::ceil(static_cast<double>(h) / mg)), n + 1);
    // A vector long enough to hold the maximum window needed
    const int64_t lenA = std::min<int64_t>(2 * maxj + 2, n + 1);
    std::vector<double> A;
    try {
        A.assign(static_cast<size_t>(lenA), 1.0);
    } catch (const std::bad_alloc &) {
        return PyErr_NoMemory();
    }

    double result = 0.0;
    bool ok;
    {
        GILRelease nogil;
        ok = outer_prob_inside(m, n, mg, ng, h, maxj, A.data(), lenA, &result);
    }
    if (!ok) {
        PyErr_SetString(PyExc_RuntimeError, "_compute_outer_prob_inside_method: internal "
                        "error, a window of the recursion left its buffer.");
        return NULL;
    }
    return PyFloat_FromDouble(result);
}


PyMethodDef ks_2samp_exact_methods[] = {
    {"_compute_outer_prob_inside_method", py_compute_outer_prob_inside_method, METH_VARARGS,
     doc_compute_outer_prob_inside_method},
    {NULL, NULL, 0, NULL}
};


PyModuleDef_Slot ks_2samp_exact_module_slots[] = {
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};


// Designated initializers require C++20; MSVC /std:c++17 rejects them.
PyModuleDef ks_2samp_exact_module = {
    PyModuleDef_HEAD_INIT,                /* m_base */
    "_ks_2samp_exact",                    /* m_name */
    NULL,                                 /* m_doc */
    0,                                    /* m_size */
    ks_2samp_exact_methods,               /* m_methods */
    ks_2samp_exact_module_slots,          /* m_slots */
    NULL,                                 /* m_traverse */
    NULL,                                 /* m_clear */
    NULL                                  /* m_free */
};

}  // namespace


PyMODINIT_FUNC
PyInit__ks_2samp_exact(void)
{
    return PyModuleDef_Init(&ks_2samp_exact_module);
}
