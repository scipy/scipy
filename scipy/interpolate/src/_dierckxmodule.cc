#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <string>
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
def _fpknot(const double[::1] x,
            const double[::1] t,
            int k,
            const double[::1] residuals):
*/
static PyObject*
py_fpknot(PyObject* self, PyObject *args)
{
    PyObject *py_x=NULL, *py_t=NULL, *py_residuals=NULL;
    int k;   // XXX: python int => "i" in ParseTuple? 

    if(!PyArg_ParseTuple(args, "OOiO", &py_x, &py_t, &k, &py_residuals)) {
        return NULL;
    }

    // XXX: Do not need to DECREF py_x etc from PyArg_ParseTuple?


/*
    // XXX: overkill? Mimic cython's `double[::1] x` etc --- replaced by check_array stanza below
    PyArrayObject *a_x=NULL, *a_t=NULL, *a_residuals=NULL;
    a_x = (PyArrayObject *)PyArray_ContiguousFromObject(py_x, NPY_DOUBLE, 1, 1);
    a_t = (PyArrayObject *)PyArray_ContiguousFromObject(py_t, NPY_DOUBLE, 1, 1);
    a_residuals = (PyArrayObject *)PyArray_ContiguousFromObject(py_residuals, NPY_DOUBLE, 1, 1);
    if (a_x == NULL || a_t == NULL || a_residuals == NULL) {
        Py_XDECREF(a_x);
        Py_XDECREF(a_t);
        Py_XDECREF(a_residuals);
        return NULL;
    }
*/

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
    //    Py_XDECREF(a_x);
    //    Py_XDECREF(a_t);
    //    Py_XDECREF(a_residuals);
        return NULL;
    }

    // heavy lifting happens here
    double new_knot = fitpack::fpknot(
        static_cast<const double *>(PyArray_DATA(a_x)), PyArray_DIM(a_x, 0),
        static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
        k,
        static_cast<const double *>(PyArray_DATA(a_residuals))
    );

    // XXX: need to DECREF a_* variables? Apparently not unless PyArray_ContiguousFromObject
 //   Py_XDECREF(a_x);
 //   Py_XDECREF(a_t);
 //   Py_XDECREF(a_residuals);

    return PyFloat_FromDouble(new_knot);
}


// XXX: rework a la py_fpknot: no extra allocations, no DECREFs

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
        // XXX: 1) is the C contiguity guaranteed? 2) DECREF on failure?
        PyErr_NoMemory();
        return NULL;     
    }

    // heavy lifting happens here
    fitpack::fpback(static_cast<const double *>(PyArray_DATA(a_R)), m, nz,
                    nc,
                    static_cast<const double *>(PyArray_DATA(a_y)), PyArray_DIM(a_y, 1),
                    static_cast<double *>(PyArray_DATA(a_c))
    );

    return (PyObject *)a_c;  // XXX like this or incref?
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
    Py_ssize_t startrow=1;  // XXX: optional, keeps the value intact if not given?

    const char *kwlist[] = {"a", "offset", "nc", "y", "startrow", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OOnO|n", const_cast<char **>(kwlist),
                                    &py_a, &py_offs, &nc, &py_y, &startrow)) {
        return NULL;
    }

    if (!(check_array(py_a, 2, NPY_DOUBLE) &&
          check_array(py_offs, 1, NPY_LONG) && // XXX: typecode for ssize_t (==ptrdiff_t)?
          check_array(py_y, 2, NPY_DOUBLE))) {
        return NULL;
    }

    PyArrayObject *a_a = (PyArrayObject *)py_a;
    PyArrayObject *a_y = (PyArrayObject *)py_y;
    PyArrayObject *a_offs = (PyArrayObject *)py_offs;

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

    //Py_DECREF(a_offs);
    // XXX: a & y modified in-place: need to incref?
    Py_RETURN_NONE;
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
