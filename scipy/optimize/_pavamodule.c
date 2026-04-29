#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdint.h>


/**
 * @brief Pool Adjacent Violators Algorithm (PAVA) for isotonic regression.
 *
 * Implements Algorithm 1 of Busing, F. M. T. A. (2022). "Monotone Regression:
 * A Simple and Fast O(n) PAVA Implementation." Journal of Statistical Software,
 * Code Snippets, 102(1), 1-25. https://doi.org/10.18637/jss.v102.c01
 *
 * The response variable @p x must be pre-sorted according to some covariate,
 * e.g. @c x = x[np.argsort(z)]. All three arrays are modified in-place.
 *
 * @param[in,out] x  Response variable array of length @p n. On return, contains
 *                    the isotonic regression solution.
 * @param[in,out] w  Case weights array of length @p n. On return, contains the
 *                    sum of weights for each block.
 * @param[in,out] r  Block index array of length @p n + 1. On return,
 *                    @c x[r[i]:r[i+1]] gives the i-th constant block.
 * @param[in]     n  Length of @p x and @p w.
 *
 * @return Number of blocks in the isotonic regression solution.
 *
 * @note Deviations from the published algorithm:
 *   - Translated to 0-based indices.
 *   - @c xb, @c wb, @c sb instead of @c x, @c w and @c S (name collisions).
 *   - @c xb_prev and @c wb_prev instead of @c x_hat and @c w_hat.
 *   - ERROR CORRECTED: Lines 9 and 10 use index @c i instead of @c b.
 *   - MODIFIED: Lines 11 and 22 use @c >= instead of @c > to get correct block
 *     indices. Otherwise equal values can land in separate blocks, e.g.
 *     @c x=[2,2] would produce @c r=[0,1,2] instead of @c r=[0,2].
 */
static int64_t
pava_algorithm(double* restrict x, double* restrict w, int64_t* restrict r, const int64_t n)
{
    int64_t b, i, f, k, t;
    double xb, wb, sb, xb_prev, wb_prev, xk;

    // procedure monotone(n, x, w)                 //  1: x in expected order
    r[0] = 0;                                      //  2: initialize index 0
    r[1] = 1;                                      //  3: initialize index 1
    b = 0;                                         //  4: block index
    xb_prev = x[b];                                //  5: prev block value
    wb_prev = w[b];                                //  6: prev block weight
    for (i = 1; i < n; ++i) {                      //  7: loop over elements
        b++;                                       //  8: increase blocks
        xb = x[i];                                 //  9: current value (i)
        wb = w[i];                                 // 10: current weight (i)
        sb = 0;
        if (xb_prev >= xb) {                       // 11: down violation (>=)
            b--;                                   // 12: decrease blocks
            sb = wb_prev * xb_prev + wb * xb;      // 13: weighted block sum
            wb += wb_prev;                         // 14: new block weight
            xb = sb / wb;                          // 15: new block value
            while (i < n - 1 && xb >= x[i + 1]) {  // 16: repair up
                i++;
                sb += w[i] * x[i];                 // 18: new weighted sum
                wb += w[i];
                xb = sb / wb;
            }
            while (b > 0 && x[b - 1] >= xb) {      // 22: repair down (>=)
                b--;
                sb += w[b] * x[b];
                wb += w[b];
                xb = sb / wb;                      // 26: new block value
            }
        }
        x[b] = xb_prev = xb;                       // 29: save block value
        w[b] = wb_prev = wb;                       // 30: save block weight
        r[b + 1] = i + 1;                          // 31: save block index
    }

    f = n - 1;                                     // 33: initialize "from"
    for (k = b; k >= 0; --k) {                     // 34: loop over blocks
        t = r[k];                                  // 35: set "to" index
        xk = x[k];
        for (i = f; i >= t; --i) {                 // 37: from downto to
            x[i] = xk;                             // 38: set to block value
        }
        f = t - 1;                                 // 40: new "from"
    }
    return b + 1;                                  // no. of blocks
}


static PyObject*
pava(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyArrayObject *ap_x = NULL, *ap_w = NULL, *ap_r = NULL;
    double* x, *w;
    int64_t* r;
    int64_t n, nb;

    if (!PyArg_ParseTuple(args, "O!O!O!",
                          &PyArray_Type, (PyObject**)&ap_x,
                          &PyArray_Type, (PyObject**)&ap_w,
                          &PyArray_Type, (PyObject**)&ap_r))
    {
        return NULL;
    }

    if (!(PyArray_IS_C_CONTIGUOUS(ap_x)) || (PyArray_TYPE(ap_x) != NPY_FLOAT64)) {
        PyErr_SetString(PyExc_TypeError, "Argument (x) must be a contiguous array of type float64.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_w)) || (PyArray_TYPE(ap_w) != NPY_FLOAT64)) {
        PyErr_SetString(PyExc_TypeError, "Argument (w) must be a contiguous array of type float64.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_r)) || (PyArray_TYPE(ap_r) != NPY_INT64)) {
        PyErr_SetString(PyExc_TypeError, "Argument (r) must be a contiguous array of type int64.");
        return NULL;
    }
    if (PyArray_NDIM(ap_x) != 1) {
        PyErr_Format(PyExc_ValueError,
                     "array has incorrect number of dimensions: %d; expected 1",
                     PyArray_NDIM(ap_x));
        return NULL;
    }

    n = (int64_t)PyArray_DIMS(ap_x)[0];
    x = (double *)PyArray_DATA(ap_x);
    w = (double *)PyArray_DATA(ap_w);
    r = (int64_t *)PyArray_DATA(ap_r);

    nb = pava_algorithm(x, w, r, n);
    PyObject* nb_py = PyLong_FromLong((long)nb);
    if (nb_py == NULL) { return NULL; }
    return Py_BuildValue("OOON", (PyObject*)ap_x, (PyObject*)ap_w, (PyObject*)ap_r, nb_py);
}


static char doc_pava[] = (
    "pava(x, w, r)\n\n"
    "Pool adjacent violators algorithm (PAVA) for isotonic regression.\n\n"
    "The routine modifies the input arguments x, w and r inplace.\n\n"
    "Parameters\n"
    "----------\n"
    "x : ndarray, shape (n,)\n"
    "    Response variable, float64. Modified in-place to contain the solution.\n"
    "w : ndarray, shape (n,)\n"
    "    Case weights, float64. Modified in-place.\n"
    "r : ndarray, shape (n+1,)\n"
    "    Block index array, intp. Modified in-place.\n\n"
    "Returns\n"
    "-------\n"
    "x : ndarray\n"
    "    The isotonic regression solution.\n"
    "w : ndarray\n"
    "    The array of weights for each block.\n"
    "r : ndarray\n"
    "    The array of indices for each block, such that x[r[i]:r[i+1]]\n"
    "    is the i-th block with all elements having the same value.\n"
    "b : int\n"
    "    Number of blocks.\n"
);


static struct PyMethodDef pava_module_methods[] = {
    {"pava", pava, METH_VARARGS, doc_pava},
    {NULL, NULL, 0, NULL}
};


static int module_exec(PyObject *module) {

    if (_import_array() < 0) { return -1; }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif
    return 0;
}


static struct PyModuleDef_Slot pava_slots[] = {
    {Py_mod_exec, module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_pava",
    .m_size = 0,
    .m_methods = pava_module_methods,
    .m_slots = pava_slots
};


PyMODINIT_FUNC
PyInit__pava(void)
{
    return PyModuleDef_Init(&moduledef);
}
