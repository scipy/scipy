#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <string>
#include <memory>
#include "numpy/arrayobject.h"
#include "__fitpack.h"

/*
 XXX:
   1. how to check for refleaks / refcounting errors?
   2. goto fail & crosses initialization errors: a better pattern or copy-paste?
   3. npy_intp vs Py_ssize_t : interchangeable or not?
   4. python int -> C int, robust pattern?
 */


/*
 * Check obj to be an ndim C contiguous array of specified typenum
 *
 *  sort of mimic checks cython does for `double[::1] obj`
 */
int
check_array(PyObject *obj, npy_intp ndim, int typenum) {

    int cond = (PyArray_CheckExact(obj) &&
                (PyArray_TYPE((PyArrayObject*)obj) == typenum) &&
                (PyArray_NDIM((PyArrayObject*)obj) == ndim) &&
                PyArray_CHKFLAGS((PyArrayObject*)obj, NPY_ARRAY_ALIGNED | NPY_ARRAY_C_CONTIGUOUS)
    );

    if(!cond) {
        std::string msg = "Expected a " + std::to_string(ndim) + "-dim C contiguous array";
        // XXX: name the dtype from typenum? Also type of arg if not array
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        return 0;
    }
    return 1;
}


/*
 * RAII malloc-d work arrays
 */
struct malloc_deleter {
    void operator()(void *p) const { std::free(p); }
};

template <typename T>
using unique_fptr = std::unique_ptr<T, malloc_deleter>;




/*
def _fpknot(const double[::1] x,
            const double[::1] t,
            int k,
            const double[::1] residuals):
*/
static PyObject*
py_fpknot(PyObject* self, PyObject *args)
{
    PyObject *py_x=NULL, *py_t=NULL, *py_residuals=NULL;
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

    // XXX: npy_intp vs Py_ssize_t vs ssize_t (== ptrdiff_t)?
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
    PyObject *py_R=NULL, *py_y=NULL;
    Py_ssize_t nc;
   // XXX: ssize_t in C++, is "n" a correct format? interchangeable with Py_ssize_t?

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
    PyObject *py_a=NULL, *py_offs=NULL, *py_y=NULL;
    Py_ssize_t nc;
    Py_ssize_t startrow=1; // The default is one (start from the 2nd row)

    // XXX: if the overhead is large, flip back to positional only arguments
    const char *kwlist[] = {"a", "offset", "nc", "y", "startrow", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OOnO|n", const_cast<char **>(kwlist),
                                    &py_a, &py_offs, &nc, &py_y, &startrow)) {
        return NULL;
    }

    if (!(check_array(py_a, 2, NPY_DOUBLE) &&
          check_array(py_offs, 1, NPY_INTP) && // XXX: typecode for ssize_t (==ptrdiff_t)?
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
            static_cast<ssize_t *>(PyArray_DATA(a_offs)),
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
 *                  const double[::1] w):
 */
static PyObject*
py_data_matrix(PyObject *self, PyObject *args)
{
    PyObject *py_x=NULL, *py_t=NULL, *py_w=NULL;
    npy_intp nc;
    int k;   // NB: declare as npy_intp, and it's garbage

    if(!PyArg_ParseTuple(args, "OOiO", &py_x, &py_t, &k, &py_w)) {
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
    PyArrayObject *a_offs = (PyArrayObject*)PyArray_ZEROS(1, dims, NPY_INTP, 0);

    unique_fptr<double> wrk( (double*)malloc((2*k+2)*sizeof(double)) );


    if ((a_A == NULL) || (a_offs == NULL) || (wrk.get() == NULL)) {
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
            static_cast<double *>(PyArray_DATA(a_A)),     // output: (A, offset, nc)
            static_cast<npy_intp*>(PyArray_DATA(a_offs)),   // XXX: callee expects ssize_t*
            &nc,
            wrk.get()
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
    PyObject *py_x=NULL, *py_t=NULL, *py_abT=NULL;
    int k, offset=0;

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
    unique_fptr<double> wrk( (double*)malloc((2*k+2)*sizeof(double)) );

    // heavy lifting happens here
    try {
        fitpack::_coloc_matrix(
            static_cast<const double *>(PyArray_DATA(a_x)), PyArray_DIM(a_x, 0),
            static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
            k,
            // abT.shape[1] is nbands because ab.shape == (nbands, nt) and abT is ab.T
            static_cast<double *>(PyArray_DATA(a_abT)), PyArray_DIM(a_abT, 1),
            offset,
            wrk.get()
        );

        Py_RETURN_NONE;
    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
}


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
    PyObject *py_x=NULL, *py_t=NULL, *py_y=NULL, *py_w=NULL, *py_abT=NULL, *py_rhs=NULL;
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
    unique_fptr<double> wrk( (double*)malloc((2*k+2)*sizeof(double)) );

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
            wrk.get()
        );

        Py_RETURN_NONE;

    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }

}


/*
 * def evaluate_spline(const double[::1] t,
 *            const double[:, ::1] c,
 *            int k,
 *            const double[::1] xp,
 *            int nu,
 *            bint extrapolate,
 *            double[:, ::1] out):
 */
static PyObject*
py_evaluate_spline(PyObject *self, PyObject *args)
{
    PyObject *py_t=NULL, *py_c=NULL, *py_xp=NULL, *py_out=NULL;
    int k, nu, i_extrap;

    if(!PyArg_ParseTuple(args, "OOiOipO", &py_t, &py_c, &k, &py_xp, &nu, &i_extrap, &py_out)) {
        return NULL;
    }

    if (!(check_array(py_t, 1, NPY_DOUBLE) &&
          check_array(py_c, 2, NPY_DOUBLE) &&
          check_array(py_xp, 1, NPY_DOUBLE) &&
          check_array(py_out, 2, NPY_DOUBLE))) {
        return NULL;
    }
    PyArrayObject *a_t = (PyArrayObject *)py_t;
    PyArrayObject *a_c = (PyArrayObject *)py_c;
    PyArrayObject *a_xp = (PyArrayObject *)py_xp;
    PyArrayObject *a_out = (PyArrayObject *)py_out;

    // check derivative order
    if (nu < 0) {
        std::string msg = "Cannot do derivative order nu= " + std::to_string(nu);
        PyErr_SetString(PyExc_NotImplementedError, msg.c_str());
        return NULL;
    }

    // sanity check sizes
    if (PyArray_DIM(a_out, 0) != PyArray_DIM(a_xp, 0)) {
        PyErr_SetString(PyExc_ValueError, "out and xp have incompatible shapes");
        return NULL;
    }
    if (PyArray_DIM(a_out, 1) != PyArray_DIM(a_c, 1)) {
        PyErr_SetString(PyExc_ValueError, "out and c have incompatible shapes");
        return NULL;
    }

    // allocate temp storage
    unique_fptr<double> wrk( (double*)malloc((2*k+2)*sizeof(double)) );

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
            wrk.get()
        );

        Py_RETURN_NONE;
    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
    

};


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


/////////////////////////////////////

static PyMethodDef DierckxMethods[] = {
    //...
    {"fpknot", py_fpknot, METH_VARARGS, 
     "fpknot replacement"},
    {"fpback", py_fpback, METH_VARARGS,
     "backsubstitution, triangular matrix"},
    {"qr_reduce", (PyCFunction)py_qr_reduce, METH_VARARGS | METH_KEYWORDS,
     "row-by-row QR triangularization"},
    {"data_matrix", py_data_matrix, METH_VARARGS,
     "(m, k+1) array of non-zero b-splines"},
    {"_coloc", py_coloc, METH_VARARGS,
     "colocation matrix in the F banded storage"},
    {"_norm_eq_lsq", py_norm_eq_lsq, METH_VARARGS,
     "lhs and rhs of the normal equations for a spline fit"},
    {"evaluate_spline", py_evaluate_spline, METH_VARARGS,
     "evaluate the spline function"},
    {"find_interval", py_find_interval, METH_VARARGS,
     doc_find_interval},
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
    import_array();

    return PyModule_Create(&dierckxmodule);
}
