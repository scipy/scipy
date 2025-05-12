#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <string>
#include <vector>
#include "numpy/arrayobject.h"
#include "__fitpack.h"


/*
 * Check obj to be an ndim C contiguous array of specified typenum
 *
 *  sort of mimic checks cython does for `double[::1] obj`
 */
static int
check_array(PyObject *obj, npy_intp ndim, int typenum) {

    int cond = (PyArray_Check(obj) &&
               (PyArray_TYPE((PyArrayObject*)obj) == typenum) &&
               (PyArray_NDIM((PyArrayObject*)obj) == ndim) &&
                PyArray_CHKFLAGS((PyArrayObject*)obj, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS)
    );

    if(!cond) {
        std::string msg = ("Expected a " + std::to_string(ndim) + "-dim C contiguous array " +
                           " of dtype = " + std::to_string(typenum) + "\n");
        // XXX: name the dtype from typenum? Also type of arg if not array
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        return 0;
    }
    return 1;
}


/*
def _fpknot(const double[::1] x,
            const double[::1] t,
            int k,
            const double[::1] residuals):
*/
static PyObject*
py_fpknot(PyObject* self, PyObject *args)
{
    PyObject *py_x = NULL, *py_t = NULL, *py_residuals = NULL;
    int k;

    if(!PyArg_ParseTuple(args, "OOiO", &py_x, &py_t, &k, &py_residuals)) {
        return NULL;
    }

    if (!(check_array(py_x, 1, NPY_DOUBLE) &&
          check_array(py_t, 1, NPY_DOUBLE) &&
          check_array(py_residuals, 1, NPY_DOUBLE))
    ) {
        return NULL;
    }

    PyArrayObject *a_x = (PyArrayObject *)py_x;
    PyArrayObject *a_t = (PyArrayObject *)py_t;
    PyArrayObject *a_residuals = (PyArrayObject *)py_residuals;

    npy_intp len_x = PyArray_DIM(a_x, 0);
    npy_intp len_r = PyArray_DIM(a_residuals, 0);

    if (len_x != len_r) {
        std::string msg = ("len(x) = " + std::to_string(len_x) + " != " +
                          std::to_string(len_r) + " = len(residuals)");
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        return NULL;
    }

    // heavy lifting happens here
    try {
        double new_knot = fitpack::fpknot(
            static_cast<const double *>(PyArray_DATA(a_x)), PyArray_DIM(a_x, 0),
            static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
            k,
            static_cast<const double *>(PyArray_DATA(a_residuals))
        );
        return PyFloat_FromDouble(new_knot);
    }
    catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;        
    };
}


/*
 * def _fpback(const double[:, ::1] R, ssize_t nc,  # (R, offset, nc) triangular => offset is range(nc)
 *             const double[:, ::1] y
 */
static PyObject*
py_fpback(PyObject* self, PyObject *args)
{
    PyObject *py_R = NULL, *py_y = NULL;
    Py_ssize_t nc;

    if(!PyArg_ParseTuple(args, "OnO", &py_R, &nc, &py_y)) {
        return NULL;
    }

    if (!(check_array(py_R, 2, NPY_DOUBLE) && check_array(py_y, 2, NPY_DOUBLE))) {
        return NULL;
    }

    PyArrayObject *a_R = (PyArrayObject *)py_R;
    PyArrayObject *a_y = (PyArrayObject *)py_y;

    // check consistency of array sizes
    Py_ssize_t m = PyArray_DIM(a_R, 0);
    Py_ssize_t nz = PyArray_DIM(a_R, 1);

    if (PyArray_DIM(a_y, 0) != m) {
        std::string msg = ("len(y) = " + std::to_string(PyArray_DIM(a_y, 0)) + " != " +
                  std::to_string(m) + " = m");
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        return NULL;
    }
    if (nc > m) {
        std::string msg = "nc = " + std::to_string(nc) + " > m = " + std::to_string(m);
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        return NULL;        
    }

    // allocate the output buffer
    npy_intp dims[2] = {nc, PyArray_DIM(a_y, 1)};
    PyArrayObject *a_c = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (a_c == NULL) {
        PyErr_NoMemory();
        return NULL;     
    }

    try {
        // heavy lifting happens here
        fitpack::fpback(static_cast<const double *>(PyArray_DATA(a_R)), m, nz,
                        nc,
                        static_cast<const double *>(PyArray_DATA(a_y)), PyArray_DIM(a_y, 1),
                        static_cast<double *>(PyArray_DATA(a_c))
        );
    }
    catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }


    return (PyObject *)a_c;
}


/*
 * def _qr_reduce(double[:, ::1] a, ssize_t[::1] offset, ssize_t nc,   # A packed
 *              double[:, ::1] y,
 *              ssize_t startrow=1
 * ):
 */
static PyObject*
py_qr_reduce(PyObject* self, PyObject *args, PyObject *kwargs)
{
    PyObject *py_a = NULL, *py_offs = NULL, *py_y = NULL;
    Py_ssize_t nc;
    Py_ssize_t startrow=1; // The default is one (start from the 2nd row)

    // XXX: if the overhead is large, flip back to positional only arguments
    const char *kwlist[] = {"a", "offset", "nc", "y", "startrow", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OOnO|n", const_cast<char **>(kwlist),
                                    &py_a, &py_offs, &nc, &py_y, &startrow)) {
        return NULL;
    }

    if (!(check_array(py_a, 2, NPY_DOUBLE) &&
          check_array(py_offs, 1, NPY_INT64) &&
          check_array(py_y, 2, NPY_DOUBLE))) {
        return NULL;
    }

    PyArrayObject *a_a = (PyArrayObject *)py_a;
    PyArrayObject *a_y = (PyArrayObject *)py_y;
    PyArrayObject *a_offs = (PyArrayObject *)py_offs;

    try {
        // heavy lifting happens here, *in-place*
        fitpack::qr_reduce(
            // a(m, nz), packed
            static_cast<double *>(PyArray_DATA(a_a)), PyArray_DIM(a_a, 0), PyArray_DIM(a_a, 1),
            // offset(m)
            static_cast<int64_t *>(PyArray_DATA(a_offs)),
            // if a were dense, it would have been a(m, nc)
            nc,
            // y(m, ydim2)
            static_cast<double *>(PyArray_DATA(a_y)), PyArray_DIM(a_y, 1),
            startrow
        );
    }
    catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}



/*
 * def _data_matrix(const double[::1] x,
 *                  const double[::1] t,
 *                  int k,
 *                  const double[::1] w,
 *                  bint extrapolate=False):
 */
static PyObject*
py_data_matrix(PyObject *self, PyObject *args)
{
    PyObject *py_x = NULL, *py_t = NULL, *py_w = NULL;
    int64_t nc;
    int k;
    int extrapolate = 0;   // default is False

    if(!PyArg_ParseTuple(args, "OOiO|p", &py_x, &py_t, &k, &py_w, &extrapolate)) {
        return NULL;
    }

    if (!(check_array(py_x, 1, NPY_DOUBLE) &&
          check_array(py_t, 1, NPY_DOUBLE) &&
          check_array(py_w, 1, NPY_DOUBLE))) {
        return NULL;
    }

    PyArrayObject *a_x = (PyArrayObject *)py_x;
    PyArrayObject *a_t = (PyArrayObject *)py_t;
    PyArrayObject *a_w = (PyArrayObject *)py_w;

    // sanity check sizes
    if (PyArray_DIM(a_w, 0) != PyArray_DIM(a_x, 0)) {
        std::string msg = ("len(w) = " + std::to_string(PyArray_DIM(a_w, 0)) + " != " +
                           "len(x) = " + std::to_string(PyArray_DIM(a_x, 0)));
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        return NULL;
    }

    // allocate temp and output arrays
    npy_intp m = PyArray_DIM(a_x, 0);
    npy_intp dims[2] = {m, k+1};
    PyArrayObject *a_A = (PyArrayObject*)PyArray_EMPTY(2, dims, NPY_DOUBLE, 0);
    // np.zeros(m, dtype=np.intp)
    PyArrayObject *a_offs = (PyArrayObject*)PyArray_ZEROS(1, dims, NPY_INT64, 0);
    std::vector<double> wrk(2*k + 2);

    if ((a_A == NULL) || (a_offs == NULL) || (wrk.data() == NULL)) {
        PyErr_NoMemory();
        Py_XDECREF(a_A);
        Py_XDECREF(a_offs);
        return NULL;
    }

    try {
        // heavy lifting happens here
        fitpack::data_matrix(
            static_cast<const double *>(PyArray_DATA(a_x)), m,
            static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
            k,
            static_cast<const double *>(PyArray_DATA(a_w)),
            extrapolate,
            static_cast<double *>(PyArray_DATA(a_A)),     // output: (A, offset, nc)
            static_cast<int64_t*>(PyArray_DATA(a_offs)),
            &nc,
            wrk.data()
        );

        // np.asarray(A), np.asarray(offset), int(nc)
        PyObject *py_nc = PyLong_FromSsize_t(static_cast<Py_ssize_t>(nc));
        return Py_BuildValue("(NNN)", PyArray_Return(a_A), PyArray_Return(a_offs), py_nc);
    }
    catch (const std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}


static char doc_coloc[] =
    "Build the B-spline colocation matrix. \n"
    "\n"
    "The colocation matrix is defined as :math:`B_{j,l} = B_l(x_j)`, \n"
    "so that row ``j`` contains all the B-splines which are non-zero \n"
    "at ``x_j``. \n"
    "\n"
    "The matrix is constructed in the LAPACK banded storage. \n"
    "Basically, for an N-by-N matrix A with ku upper diagonals and \n"
    "kl lower diagonals, the shape of the array Ab is (2*kl + ku +1, N), \n"
    "where the last kl+ku+1 rows of Ab contain the diagonals of A, and \n"
    "the first kl rows of Ab are not referenced. \n"
    "For more info see, e.g. the docs for the ``*gbsv`` routine. \n"
    "\n"
    "This routine is not supposed to be called directly, and \n"
    "does no error checking. \n"
    "\n"
    "Parameters \n"
    "---------- \n"
    "x : ndarray, shape (n,) \n"
    "    sorted 1D array of x values \n"
    "t : ndarray, shape (nt + k + 1,) \n"
    "    sorted 1D array of knots \n"
    "k : int \n"
    "    spline order \n"
    "ab : ndarray, shape (2*kl + ku + 1, nt), F-order \n"
    "    This parameter is modified in-place. \n"
    "    On exit: B-spline colocation matrix in the band storage with \n"
    "    ``ku`` upper diagonals and ``kl`` lower diagonals. \n"
    "    Here ``kl = ku = k``. \n"
    "offset : int, optional \n"
    "    skip this many rows \n";
/*
 * def _coloc(const double[::1] x,
 *             const double[::1] t,
 *             int k,
 *             double[::1, :] ab,  // NOTE: receive the transposed ab
 *             int offset=0):
 *
 * Similar to data_matrix:
 *   - data_matrix constructs the array in the "packed" storage;
 *   - _coloc constructs the array in the LAPACK banded storage;
 */
static PyObject*
py_coloc(PyObject *self, PyObject *args)
{
    PyObject *py_x = NULL, *py_t = NULL, *py_abT = NULL;
    int k, offset = 0;

    if(!PyArg_ParseTuple(args, "OOiOi", &py_x, &py_t, &k, &py_abT, &offset)) {
        return NULL;
    }

    if (!(check_array(py_x, 1, NPY_DOUBLE) &&
          check_array(py_t, 1, NPY_DOUBLE) &&
          check_array(py_abT, 2, NPY_DOUBLE))) {
        return NULL;
    }
    PyArrayObject *a_x = (PyArrayObject *)py_x;
    PyArrayObject *a_t = (PyArrayObject *)py_t;
    PyArrayObject *a_abT = (PyArrayObject *)py_abT;

    // allocate the temp storage
    std::vector<double> wrk(2*k + 2);

    // heavy lifting happens here
    try {
        fitpack::_coloc_matrix(
            static_cast<const double *>(PyArray_DATA(a_x)), PyArray_DIM(a_x, 0),
            static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
            k,
            // abT.shape[1] is nbands because ab.shape == (nbands, nt) and abT is ab.T
            static_cast<double *>(PyArray_DATA(a_abT)), PyArray_DIM(a_abT, 1),
            offset,
            wrk.data()
        );

        Py_RETURN_NONE;
    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}

static char doc_norm_eq_lsq[] =
    "Construct the normal equations for the B-spline LSQ problem. \n"
    " \n"
    "The observation equations are ``A @ c = y``, and the normal equations are \n"
    "``A.T @ A @ c = A.T @ y``. This routine fills in the rhs and lhs for the \n"
    "latter. \n"
    " \n"
    "The B-spline collocation matrix is defined as :math:`A_{j,l} = B_l(x_j)`, \n"
    "so that row ``j`` contains all the B-splines which are non-zero \n"
    "at ``x_j``. \n"
    " \n"
    "The normal eq matrix has at most `2k+1` bands and is constructed in the \n"
    "LAPACK symmetrix banded storage: ``A[i, j] == ab[i-j, j]`` with `i >= j`. \n"
    "See the doctsring for `scipy.linalg.cholesky_banded` for more info. \n"
    " \n"
    "Parameters \n"
    "---------- \n"
    "x : ndarray, shape (n,) \n"
    "    sorted 1D array of x values \n"
    "t : ndarray, shape (nt + k + 1,) \n"
    "    sorted 1D array of knots \n"
    "k : int \n"
    "    spline order \n"
    "y : ndarray, shape (n, s) \n"
    "    a 2D array of y values. The second dimension contains all trailing \n"
    "    dimensions of the original array of ordinates. \n"
    "w : ndarray, shape(n,) \n"
    "    Weights. \n"
    "ab : ndarray, shape (k+1, n), in Fortran order. \n"
    "    This parameter is modified in-place. \n"
    "    On entry: should be zeroed out. \n"
    "    On exit: LHS of the normal equations. \n"
    "rhs : ndarray, shape (n, s), in C order. \n"
    "    This parameter is modified in-place. \n"
    "    On entry: should be zeroed out. \n"
    "    On exit: RHS of the normal equations. \n";

/*
 * def _norm_eq_lsq(const double[::1] x,
 *                  const double[::1] t,
 *                  int k,
 *                  const double[:, ::1] y,
 *                  const double[::1] w,
 *                  double[::1, :] ab,
 *                  double[:, ::1] rhs):
 */
static PyObject*
py_norm_eq_lsq(PyObject *self, PyObject *args)
{
    PyObject *py_x = NULL, *py_t = NULL, *py_y = NULL, *py_w = NULL, *py_abT = NULL, *py_rhs = NULL;
    int k;

    if(!PyArg_ParseTuple(args, "OOiOOOO", &py_x, &py_t, &k, &py_y, &py_w, &py_abT, &py_rhs)) {
        return NULL;
    }

    if (!(check_array(py_x, 1, NPY_DOUBLE) &&
          check_array(py_t, 1, NPY_DOUBLE) &&
          check_array(py_y, 2, NPY_DOUBLE) &&
          check_array(py_w, 1, NPY_DOUBLE) &&
          check_array(py_abT, 2, NPY_DOUBLE) &&
          check_array(py_rhs, 2, NPY_DOUBLE))) {
        return NULL;
    }
    PyArrayObject *a_x = (PyArrayObject *)py_x;
    PyArrayObject *a_t = (PyArrayObject *)py_t;
    PyArrayObject *a_y = (PyArrayObject *)py_y;
    PyArrayObject *a_w = (PyArrayObject *)py_w;
    PyArrayObject *a_abT = (PyArrayObject *)py_abT;
    PyArrayObject *a_rhs = (PyArrayObject *)py_rhs;

    // allocate temp storage
    std::vector<double> wrk(2*k + 2);

    // heavy lifting happens here
    try {
        fitpack::norm_eq_lsq(
            static_cast<const double*>(PyArray_DATA(a_x)), PyArray_DIM(a_x, 0),
            static_cast<const double*>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
            k,
            static_cast<const double*>(PyArray_DATA(a_y)), PyArray_DIM(a_y, 1),
            static_cast<const double*>(PyArray_DATA(a_w)),
            static_cast<double*>(PyArray_DATA(a_abT)),
            static_cast<double*>(PyArray_DATA(a_rhs)),
            wrk.data()
        );

        Py_RETURN_NONE;

    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }

}


static char doc_evaluate_spline[] =
   "Evaluate a spline in the B-spline basis. \n"
   "\n"
   "Parameters \n"
   "---------- \n"
   "t : ndarray, shape (n+k+1) \n"
   "     knots \n"
   " c : ndarray, shape (n, m) \n"
   "     B-spline coefficients \n"
   " xp : ndarray, shape (s,) \n"
   "     Points to evaluate the spline at. \n"
   " nu : int \n"
   "     Order of derivative to evaluate. \n"
   " extrapolate : int, optional \n"
   "     Whether to extrapolate to ouf-of-bounds points, or to return NaNs. \n"
   "\n"
   "Returns\n"
   "-------"
   " out : ndarray, shape (s, m) \n"
   "     Computed values of the spline at each of the input points. \n";
/*
 * def evaluate_spline(const double[::1] t,
 *            const double[:, ::1] c,
 *            int k,
 *            const double[::1] xp,
 *            int nu,
 *            bint extrapolate):
 */
static PyObject*
py_evaluate_spline(PyObject *self, PyObject *args)
{
    PyObject *py_t = NULL, *py_c = NULL, *py_xp = NULL;
    int k, nu, i_extrap;

    if(!PyArg_ParseTuple(args, "OOiOip", &py_t, &py_c, &k, &py_xp, &nu, &i_extrap)) {
        return NULL;
    }

    if (!(check_array(py_t, 1, NPY_DOUBLE) &&
          check_array(py_c, 2, NPY_DOUBLE) &&
          check_array(py_xp, 1, NPY_DOUBLE))) {
        return NULL;
    }
    PyArrayObject *a_t = (PyArrayObject *)py_t;
    PyArrayObject *a_c = (PyArrayObject *)py_c;
    PyArrayObject *a_xp = (PyArrayObject *)py_xp;

    // check derivative order
    if (nu < 0) {
        std::string msg = "Cannot do derivative order nu= " + std::to_string(nu);
        PyErr_SetString(PyExc_NotImplementedError, msg.c_str());
        return NULL;
    }

    // allocate temp storage
    std::vector<double> wrk(2*k + 2);

    // allocate the output array, shape (x.shape[0], c.shape[-1])
    npy_intp dims[2] = {PyArray_DIM(a_xp, 0), PyArray_DIM(a_c, 1)};
    PyArrayObject *a_out = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (a_out == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    // heavy lifting happens here
    try {
        fitpack::_evaluate_spline(
            static_cast<const double*>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
            static_cast<const double*>(PyArray_DATA(a_c)), PyArray_DIM(a_c, 0), PyArray_DIM(a_c, 1),
            k,
            static_cast<const double*>(PyArray_DATA(a_xp)), PyArray_DIM(a_xp, 0),
            nu,
            i_extrap,
            static_cast<double*>(PyArray_DATA(a_out)),
            wrk.data()
        );

        return (PyObject *)(a_out);
    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
};


static char doc_evaluate_all_bspl[] =
    "Evaluate the ``k+1`` B-splines which are non-zero on interval ``m``. \n"
    " \n"
    "Parameters \n"
    "---------- \n"
    "t : ndarray, shape (nt + k + 1,) \n"
    "    sorted 1D array of knots \n"
    "k : int \n"
    "    spline order \n"
    "xval: float \n"
    "    argument at which to evaluate the B-splines \n"
    "m : int \n"
    "    index of the left edge of the evaluation interval, ``t[m] <= x < t[m+1]`` \n"
    "nu : int, optional \n"
    "    Evaluate derivatives order `nu`. Default is zero. \n"
    " \n"
    "Returns \n"
    "------- \n"
    "ndarray, shape (k+1,) \n"
    "    The values of B-splines :math:`[B_{m-k}(xval), ..., B_{m}(xval)]` if \n"
    "    `nu` is zero, otherwise the derivatives of order `nu`. \n"
    " \n"
    "Examples \n"
    "-------- \n"
    " \n"
    "A textbook use of this sort of routine is plotting the ``k+1`` polynomial \n"
    "pieces which make up a B-spline of order `k`. \n"
    " \n"
    "Consider a cubic spline \n"
    " \n"
    ">>> k = 3 \n"
    ">>> t = [0., 1., 2., 3., 4.]   # internal knots \n"
    ">>> a, b = t[0], t[-1]    # base interval is [a, b) \n"
    ">>> t = np.array([a]*k + t + [b]*k)  # add boundary knots \n"
    " \n"
    ">>> import matplotlib.pyplot as plt \n"
    ">>> xx = np.linspace(a, b, 100) \n"
    ">>> plt.plot(xx, BSpline.basis_element(t[k:-k])(xx), \n"
    "...          lw=3, alpha=0.5, label='basis_element') \n"
    " \n"
    "Now we use slide an interval ``t[m]..t[m+1]`` along the base interval \n"
    "``a..b`` and use `evaluate_all_bspl` to compute the restriction of \n"
    "the B-spline of interest to this interval: \n"
    " \n"
    ">>> for i in range(k+1): \n"
    "...    x1, x2 = t[2*k - i], t[2*k - i + 1] \n"
    "...    xx = np.linspace(x1 - 0.5, x2 + 0.5) \n"
    "...    yy = [evaluate_all_bspl(t, k, x, 2*k - i)[i] for x in xx] \n"
    "...    plt.plot(xx, yy, '--', label=str(i)) \n"
    ">>> plt.grid(True) \n"
    ">>> plt.legend() \n"
    ">>> plt.show() \n";
/*
 * def evaluate_all_bspl(const double[::1] t, int k, double xval, int m, int nu=0):
 */
static PyObject*
py_evaluate_all_bspl(PyObject* self, PyObject* args)
{
    PyObject *py_t = NULL;
    int k, m;
    int nu = 0;
    double xval;

    if(!PyArg_ParseTuple(args, "Oidi|i", &py_t, &k, &xval, &m, &nu)) {
        return NULL;
    }
    if (!check_array(py_t, 1, NPY_DOUBLE)) {
        return NULL;
    }
    PyArrayObject *a_t = (PyArrayObject *)py_t;

    // allocate temp storage
    std::vector<double> wrk(2*k + 2);

    // compute non-zero bsplines
    fitpack::_deBoor_D(
        static_cast<const double*>(PyArray_DATA(a_t)),
        xval,
        k,
        m,
        nu,
        wrk.data()
    );

    // allocate and fill the output
    npy_intp dims[1] = {k+1};
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (arr == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    memcpy(PyArray_DATA(arr), wrk.data(), (k+1)*sizeof(double));
    return (PyObject *)arr;
}


static char doc_find_interval[] = 
    "Find an interval such that t[interval] <= xval < t[interval+1]. \n"
    "\n"
    "Uses a linear search with locality, see fitpack's splev. \n"
    "\n"
    "Parameters \n"
    "---------- \n"
    "t : ndarray, shape (nt,) \n"
    "    Knots \n"
    "k : int \n"
    "    B-spline degree \n"
    "xval : double \n"
    "    value to find the interval for \n"
    "prev_l : int \n"
    "    interval where the previous value was located. \n"
    "    if unknown, use any value < k to start the search. \n"
    "extrapolate : int \n"
    "    whether to return the last or the first interval if xval \n"
    "    is out of bounds. \n"
    "\n"
    "Returns \n"
    "------- \n"
    "interval : int \n"
    "    Suitable interval or -1 if xval was nan. \n";
/*
 * def _py_find_interval(const double[::1] t,
 *                       int k,
 *                       double xval,
 *                       int prev_l,
 *                       bint extrapolate):
 */
static PyObject*
py_find_interval(PyObject *self, PyObject *args)
{
    PyObject *py_t = NULL;
    double xval;
    int k, prev_l, i_extrap;

    if(!PyArg_ParseTuple(args, "Oidip", &py_t, &k, &xval, &prev_l, &i_extrap)) {
        return NULL;
    }

    if (!check_array(py_t, 1, NPY_DOUBLE)) {
        return NULL;
    }
    PyArrayObject *a_t = (PyArrayObject *)py_t;

    Py_ssize_t interval = fitpack::_find_interval(
        static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
        k,
        xval,
        prev_l,
        i_extrap
    );

    PyObject *py_interval = PyLong_FromSsize_t(interval);
    return py_interval;
}


/*** NDBspline ***/


static char doc_evaluate_ndbspline[] =
        "Evaluate an N-dim tensor product spline or its derivative.\n"
        "\n"
        "Parameters\n"
        "----------\n"
        "xi : ndarray, shape(npoints, ndim)\n"
        "    ``npoints`` values to evaluate the spline at, each value is\n"
        "    a point in an ``ndim``-dimensional space.\n"
        "t : ndarray, shape(ndim, max_len_t)\n"
        "    Array of knots for each dimension.\n"
        "    This array packs the tuple of knot arrays per dimension into a single\n"
        "    2D array. The array is ragged (knot lengths may differ), hence\n"
        "    the real knots in dimension ``d`` are ``t[d, :len_t[d]]``.\n"
        "len_t : ndarray, 1D, shape (ndim,)\n"
        "    Lengths of the knot arrays, per dimension.\n"
        "k : tuple of ints, len(ndim)\n"
        "    Spline degrees in each dimension.\n"
        "nu : ndarray of ints, shape(ndim,)\n"
        "    Orders of derivatives to compute, per dimension.\n"
        "extrapolate : int\n"
        "    Whether to extrapolate out of bounds or return nans.\n"
        "c1r: ndarray, one-dimensional\n"
        "    Flattened array of coefficients.\n"
        "    The original N-dimensional coefficient array ``c`` has shape\n"
        "    ``(n1, ..., nd, ...)`` where each ``ni == len(t[d]) - k[d] - 1``,\n"
        "    and the second '...' represents trailing dimensions of ``c``.\n"
        "    In code, given the C-ordered array ``c``, ``c1r`` is\n"
        "    ``c1 = c.reshape(c.shape[:ndim] + (-1,)); c1r = c1.ravel()``\n"
        "num_c_tr : int\n"
        "    The number of elements of ``c1r``, which correspond to the trailing\n"
        "    dimensions of ``c``. In code, this is\n"
        "    ``c1 = c.reshape(c.shape[:ndim] + (-1,)); num_c_tr = c1.shape[-1]``.\n"
        "strides_c1 : ndarray, one-dimensional\n"
        "    Pre-computed strides of the ``c1`` array.\n"
        "    Note: These are *data* strides, not numpy-style byte strides.\n"
        "    This array is equivalent to\n"
        "    ``[stride // s1.dtype.itemsize for stride in s1.strides]``.\n"
        "indices_k1d : ndarray, shape((k+1)**ndim, ndim)\n"
        "    Pre-computed mapping between indices for iterating over a flattened\n"
        "    array of shape ``[k[d] + 1) for d in range(ndim)`` and\n"
        "    ndim-dimensional indices of the ``(k+1,)*ndim`` dimensional array.\n"
        "    This is essentially a transposed version of\n"
        "    ``np.unravel_index(np.arange((k+1)**ndim), (k+1,)*ndim)``.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "out : ndarray, shape (npoints, num_c_tr)\n"
        "    Output values of the b-spline at given ``xi`` points.\n"
        "\n"
        "Notes\n"
        "-----\n"
        "\n"
        "This function is essentially equivalent to the following: given an\n"
        "N-dimensional vector ``x = (x1, x2, ..., xN)``, iterate over the\n"
        "dimensions, form linear combinations of products,\n"
        "B(x1) * B(x2) * ... B(xN) of (k+1)**N b-splines which are non-zero\n"
        "at ``x``.\n"
        "\n"
        "Since b-splines are localized, the sum has (k+1)**N non-zero elements.\n"
        "\n"
        "If ``i = (i1, i2, ..., iN)`` is a vector if intervals of the knot\n"
        "vectors, ``t[d, id] <= xd < t[d, id+1]``, for ``d=1, 2, ..., N``, then\n"
        "the core loop of this function is nothing but\n"
        "\n"
        "```\n"
        "result = 0\n"
        "iters = [range(i[d] - self.k[d], i[d] + 1) for d in range(ndim)]\n"
        "for idx in itertools.product(*iters):\n"
        "    term = self.c[idx] * np.prod([B(x[d], self.k[d], idx[d], self.t[d])\n"
        "                                  for d in range(ndim)])\n"
        "    result += term\n"
        "```\n"
        "\n"
        "For efficiency reasons, we iterate over the flattened versions of the arrays.\n";
/*
def evaluate_ndbspline(const double[:, ::1] xi,
                       const double[:, ::1] t,
                       const npy_int64[::1] len_t,
                       const npy_int64[::1] k,
                       npy_int64[::1] nu,
                       bint extrapolate,
                       const double[::1] c1r,
                       int num_c_tr,
                       const npy_int64[::1] strides_c1,
                       const npy_int64[:, ::] indices_k1d,
*/
static PyObject*
py_evaluate_ndbspline(PyObject *self, PyObject *args)
{
    PyObject *py_xi=NULL;
    PyObject *py_t=NULL, *py_c1r=NULL, *py_strides_c1=NULL, *py_indices_k1d=NULL;

    PyObject *py_len_t=NULL, *py_k=NULL, *py_nu=NULL;
    int num_c_tr;
    int i_extrap;

    if(!PyArg_ParseTuple(args, "OOOOOiOiOO",
                         &py_xi, &py_t, &py_len_t, &py_k, &py_nu, &i_extrap,
                         &py_c1r, &num_c_tr, &py_strides_c1, &py_indices_k1d)) {
        return NULL;
    }

    if (!(check_array(py_xi, 2, NPY_DOUBLE) &&
          check_array(py_t, 2, NPY_DOUBLE) &&
          check_array(py_len_t, 1, NPY_INT64) &&
          check_array(py_k, 1, NPY_INT64) &&
          check_array(py_nu, 1, NPY_INT64) &&
          check_array(py_c1r, 1, NPY_DOUBLE) &&
          check_array(py_strides_c1, 1, NPY_INT64) &&
          check_array(py_indices_k1d, 2, NPY_INT64))) {
        return NULL;
    }
    PyArrayObject *a_xi = (PyArrayObject *)py_xi;
    PyArrayObject *a_t = (PyArrayObject *)py_t;

    PyArrayObject *a_len_t = (PyArrayObject *)py_len_t;
    PyArrayObject *a_k = (PyArrayObject *)py_k;
    PyArrayObject *a_nu = (PyArrayObject *)py_nu;

    PyArrayObject *a_c1r = (PyArrayObject *)py_c1r;
    PyArrayObject *a_strides_c1 = (PyArrayObject *)py_strides_c1;
    PyArrayObject *a_indices_k1d = (PyArrayObject *)py_indices_k1d;

    // sanity checks
    int64_t ndim = PyArray_DIM(a_t, 0);
    if (PyArray_DIM(a_xi, 1) != ndim) {
        std::string msg = ("Expected data points in " + std::to_string(ndim) + "-D"
                           " space, got " + std::to_string(PyArray_DIM(a_xi, 1)) +
                           "-D points.");
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        return NULL;
    }

    // allocate the output
    npy_intp dims[2] = {PyArray_DIM(a_xi, 0), num_c_tr};
    PyArrayObject *a_out = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (a_out == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    // heavy lifting happens here
    try {
        fitpack::_evaluate_ndbspline(
            /* inputs */
            static_cast<const double *>(PyArray_DATA(a_xi)), PyArray_DIM(a_xi, 0), PyArray_DIM(a_xi, 1),
            static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 1),
            static_cast<const int64_t *>(PyArray_DATA(a_len_t)),
            static_cast<const int64_t *>(PyArray_DATA(a_k)),
            static_cast<const int64_t *>(PyArray_DATA(a_nu)),
            i_extrap,
            /* flattened coefficients */
            static_cast<const double *>(PyArray_DATA(a_c1r)), PyArray_DIM(a_c1r, 0),
            /* tabulated helpers */
            static_cast<const int64_t *>(PyArray_DATA(a_strides_c1)),
            static_cast<const int64_t *>(PyArray_DATA(a_indices_k1d)), PyArray_DIM(a_indices_k1d, 0),

            /* output */
            static_cast<double*>(PyArray_DATA(a_out)), num_c_tr
        );

        return (PyObject *)(a_out);
    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}


static char doc_coloc_nd[] =
    "Construct the N-D tensor product collocation matrix as a CSR array.\n"
    "\n"
    "In the dense representation, each row of the collocation matrix corresponds\n"
    "to a data point and contains non-zero b-spline basis functions which are\n"
    "non-zero at this data point.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "xvals : ndarray, shape(size, ndim)\n"
    "    Data points. ``xvals[j, :]`` gives the ``j``-th data point as an\n"
    "    ``ndim``-dimensional array.\n"
    "t : tuple of 1D arrays, length-ndim\n"
    "    Tuple of knot vectors\n"
    "k : ndarray, shape (ndim,)\n"
    "    Spline degrees\n"
    "\n"
    "Returns\n"
    "-------\n"
    "csr_data, csr_indices, csr_indptr\n"
    "    The collocation matrix in the CSR array format.\n"
    "\n"
    "Notes\n"
    "-----\n"
    "Algorithm: given `xvals` and the tuple of knots `t`, we construct a tensor\n"
    "product spline, i.e. a linear combination of\n"
    "\n"
    "   B(x1; i1, t1) * B(x2; i2, t2) * ... * B(xN; iN, tN)\n"
    "\n"
    "Here ``B(x; i, t)`` is the ``i``-th b-spline defined by the knot vector\n"
    "``t`` evaluated at ``x``.\n"
    "\n"
    "Since ``B`` functions are localized, for each point `(x1, ..., xN)` we\n"
    "loop over the dimensions, and\n"
    "- find the location in the knot array, `t[i] <= x < t[i+1]`,\n"
    "- compute all non-zero `B` values\n"
    "- place these values into the relevant row\n"
    "\n"
    "In the dense representation, the collocation matrix would have had a row per\n"
    "data point, and each row has the values of the basis elements (i.e., tensor\n"
    "products of B-splines) evaluated at this data point. Since the matrix is very\n"
    "sparse (has size = len(x)**ndim, with only (k+1)**ndim non-zero elements per\n"
    "row), we construct it in the CSR format.\n";
/*
def _colloc_nd(const double[:, ::1] xvals,
               const double[:, ::1] _t,
               const npy_int64[::1] len_t,
               const npy_int64[::1] k,
               const npy_int64[:, ::1] _indices_k1d,
               const npy_int64[::1] _cstrides):
*/
static PyObject*
py_coloc_nd(PyObject *self, PyObject *args)
{
    PyObject *py_xi, *py_t, *py_len_t, *py_k, *py_indices_k1d, *py_strides;

    if(!PyArg_ParseTuple(args, "OOOOOO",
                         &py_xi, &py_t, &py_len_t, &py_k,
                         &py_indices_k1d, &py_strides)) {
        return NULL;
    }

    if (!(check_array(py_xi, 2, NPY_DOUBLE) &&
          check_array(py_t, 2, NPY_DOUBLE) &&
          check_array(py_len_t, 1, NPY_INT64) &&
          check_array(py_k, 1, NPY_INT64) &&
          check_array(py_indices_k1d, 2, NPY_INT64) &&
          check_array(py_strides, 1, NPY_INT64))) {
        return NULL;
    }
    PyArrayObject *a_xi = (PyArrayObject *)py_xi;
    PyArrayObject *a_t = (PyArrayObject *)py_t;
    PyArrayObject *a_len_t = (PyArrayObject *)py_len_t;
    PyArrayObject *a_k = (PyArrayObject *)py_k;
    PyArrayObject *a_indices_k1d = (PyArrayObject *)py_indices_k1d;
    PyArrayObject *a_strides = (PyArrayObject *)py_strides;

    /* allocate the outputs */
    npy_intp npts = PyArray_DIM(a_xi, 0);
    npy_intp ndim = PyArray_DIM(a_xi, 1);

    // the number of non-zero b-splines at each data point
    npy_intp volume = 1;
    int64_t *k_data = static_cast<int64_t *>(PyArray_DATA(a_k));
    for (int d=0; d < ndim; d++) {
        volume *= k_data[d] + 1;
    }

    // Allocate the colocation matrix in the CSR format.
    npy_intp dims[1] = {npts*volume};
    PyObject *py_csr_data = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *py_csr_indices = PyArray_SimpleNew(1, dims, NPY_INT64);
    PyObject *py_csr_indptr = PyArray_Arange(0, volume*npts + 1, volume, NPY_INT64);

    if ((py_csr_data == NULL) || (py_csr_indices == NULL) || (py_csr_indptr == NULL)) {
        PyErr_NoMemory();
        return NULL;
    }

    PyArrayObject *a_csr_data = (PyArrayObject *)py_csr_data;
    PyArrayObject *a_csr_indices = (PyArrayObject *)py_csr_indices;

    // heavy lifting happens here
    try {
        int status = fitpack::_coloc_nd(
            /* inputs */
            static_cast<const double *>(PyArray_DATA(a_xi)), npts, ndim,
            static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 1),
            static_cast<const int64_t *>(PyArray_DATA(a_len_t)),
            static_cast<const int64_t *>(PyArray_DATA(a_k)),
            /* tabulated helpers */
            static_cast<const int64_t *>(PyArray_DATA(a_indices_k1d)), PyArray_DIM(a_indices_k1d, 0),
            static_cast<const int64_t *>(PyArray_DATA(a_strides)),
            /* outputs */
            static_cast<int64_t *>(PyArray_DATA(a_csr_indices)), volume,
            static_cast<double *>(PyArray_DATA(a_csr_data))
        );
        if (status < 0) {
            std::string mesg = ("Data point " + std::to_string(-status) + " is out of bounds");
            PyErr_SetString(PyExc_ValueError, mesg.c_str());
        }

        return Py_BuildValue("(NNN)", PyArray_Return(a_csr_data),
                                      PyArray_Return(a_csr_indices),
                                      py_csr_indptr
        );
    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}


/////////////////////////////////////

static PyMethodDef DierckxMethods[] = {
    /* FITPACK replacement helpers*/
    {"fpknot", py_fpknot, METH_VARARGS, 
     "fpknot replacement"},
    {"fpback", py_fpback, METH_VARARGS,
     "backsubstitution, triangular matrix"},
    {"qr_reduce", (PyCFunction)py_qr_reduce, METH_VARARGS | METH_KEYWORDS,
     "row-by-row QR triangularization"},
    {"data_matrix", py_data_matrix, METH_VARARGS,
     "(m, k+1) array of non-zero b-splines"},
    /* BSpline helpers */
    {"evaluate_spline", py_evaluate_spline, METH_VARARGS,
     doc_evaluate_spline},
    {"evaluate_all_bspl", py_evaluate_all_bspl, METH_VARARGS,
     doc_evaluate_all_bspl},
    {"find_interval", py_find_interval, METH_VARARGS,
     doc_find_interval},
    /* make_{interp,lsq}_spline helpers*/
    {"_coloc", py_coloc, METH_VARARGS,
      doc_coloc},
    {"_norm_eq_lsq", py_norm_eq_lsq, METH_VARARGS,
     doc_norm_eq_lsq},
    /* NdBSpline helpers */
    {"evaluate_ndbspline", py_evaluate_ndbspline, METH_VARARGS,
     doc_evaluate_ndbspline},
    {"_coloc_nd", py_coloc_nd, METH_VARARGS,
     doc_coloc_nd},
    //...
    {NULL, NULL, 0, NULL}        /* Sentinel */
};



static struct PyModuleDef dierckxmodule = {
    PyModuleDef_HEAD_INIT,
    "_dierckx",   /* name of module */
    NULL, //spam_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    DierckxMethods
};


PyMODINIT_FUNC
PyInit__dierckx(void)
{
    PyObject *module;

    import_array();

    module = PyModule_Create(&dierckxmodule);
    if (module == NULL) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}
