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


static PyObject*
py_fpknot(PyObject* self, PyObject *args)
{
    PyObject *py_x=NULL, *py_t=NULL, *py_residuals=NULL;
    PyArrayObject *a_x=NULL, *a_t=NULL, *a_residuals=NULL;
    int k;   // XXX: python int => ? 

    if(!PyArg_ParseTuple(args, "OOiO", &py_x, &py_t, &k, &py_residuals)) {
        return NULL;
    }

    // XXX: overkill? Mimic cython's `double[::1] x` etc
    a_x = (PyArrayObject *)PyArray_ContiguousFromObject(py_x, NPY_DOUBLE, 1, 1);
    a_t = (PyArrayObject *)PyArray_ContiguousFromObject(py_t, NPY_DOUBLE, 1, 1);
    a_residuals = (PyArrayObject *)PyArray_ContiguousFromObject(py_residuals, NPY_DOUBLE, 1, 1);

    if (a_x == NULL || a_t == NULL || a_residuals == NULL) {
        Py_XDECREF(a_x);
        Py_XDECREF(a_t);
        Py_XDECREF(a_residuals);
        return NULL;
    }

    // XXX: npy_intp vs Py_ssize_t ?
    npy_intp len_x = PyArray_DIM(a_x, 0);
    npy_intp len_r = PyArray_DIM(a_residuals, 0);

    if (len_x != len_r) {
        std::string msg = ("len(x) = " + std::to_string(len_x) + " != " +
                          std::to_string(len_r) + " = len(residuals)");
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        Py_XDECREF(a_x);
        Py_XDECREF(a_t);
        Py_XDECREF(a_residuals);
        return NULL;
    }

    double new_knot = fitpack::fpknot(
        static_cast<const double *>(PyArray_DATA(a_x)), PyArray_DIM(a_x, 0),
        static_cast<const double *>(PyArray_DATA(a_t)), PyArray_DIM(a_t, 0),
        k,
        static_cast<const double *>(PyArray_DATA(a_residuals))
    );

    // XXX: need to DECREF a_* variables?
    Py_XDECREF(a_x);
    Py_XDECREF(a_t);
    Py_XDECREF(a_residuals);

    return PyFloat_FromDouble(new_knot);
}




static PyObject*
py_fpback(PyObject* self, PyObject *args)
{
    PyObject *py_R=NULL, *py_y=NULL;
    PyArrayObject *a_R=NULL, *a_y=NULL;
    Py_ssize_t nc;
   // XXX: ssize_t in C++, is "n" a correct format? interchangeable with Py_ssize_t?

    if(!PyArg_ParseTuple(args, "OnO", &py_R, &nc, &py_y)) {
        return NULL;
    }

    // XXX: double[:, ::1] R, double[:, ::1] y
    a_R = (PyArrayObject *)PyArray_ContiguousFromObject(py_R, NPY_DOUBLE, 2, 2);
    a_y = (PyArrayObject *)PyArray_ContiguousFromObject(py_y, NPY_DOUBLE, 2, 2);

    if (a_R == NULL || a_y == NULL) {
        Py_XDECREF(a_R);
        Py_XDECREF(a_y);
        return NULL;
    }

    // check consistency of array sizes
    Py_ssize_t m = PyArray_DIM(a_R, 0);
    Py_ssize_t nz = PyArray_DIM(a_R, 1);

    if (PyArray_DIM(a_y, 0) != m) {
        std::string msg = ("len(y) = " + std::to_string(PyArray_DIM(a_y, 0)) + " != " +
                  std::to_string(m) + " = m");
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        Py_XDECREF(a_R);
        Py_XDECREF(a_y);
        return NULL;
    }
    if (nc > m) {
        std::string msg = "nc = " + std::to_string(nc) + " > m = " + std::to_string(m);
        PyErr_SetString(PyExc_ValueError, msg.c_str());
        Py_XDECREF(a_R);
        Py_XDECREF(a_y);
        return NULL;        
    }

    // allocate the output buffer
    npy_intp dims[2] = {nc, PyArray_DIM(a_y, 1)};
    PyArrayObject *a_c = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (a_c == NULL) {
        // XXX: 1) need to guard PyArray_SimpleNew?   2) is the C contiguity guaranteed?
        PyErr_NoMemory();
        Py_XDECREF(a_R);
        Py_XDECREF(a_y);
        return NULL;     
    }

    // heavy lifting happens here
    fitpack::fpback(static_cast<const double *>(PyArray_DATA(a_R)), m, nz,
                    nc,
                    static_cast<const double *>(PyArray_DATA(a_y)), PyArray_DIM(a_y, 1),
                    static_cast<double *>(PyArray_DATA(a_c))
    );

    Py_DECREF(a_R);
    Py_DECREF(a_y);
    return (PyObject *)a_c;  // XXX like this? incref?
}



static PyObject*
py_qr_reduce(PyObject* self, PyObject *args, PyObject *kwargs)
{
    PyObject *py_a=NULL, *py_offs=NULL, *py_y=NULL;
    PyArrayObject *a_a=NULL, *a_offs=NULL, *a_y=NULL;
    Py_ssize_t nc;
    Py_ssize_t startrow=1;  // XXX: optional, keeps the value intact if not given?

    const char *kwlist[] = {"a", "offset", "nc", "y", "startrow", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OOnO|n", const_cast<char **>(kwlist),
                                    &py_a, &py_offs, &nc, &py_y, &startrow)) {
        return NULL;
    }

    // XXX: double[:, ::1] a, double[:, ::1] y; a&y are modified in in-place
    a_a = (PyArrayObject *)PyArray_ContiguousFromObject(py_a, NPY_DOUBLE, 2, 2);
    a_offs = (PyArrayObject *)PyArray_ContiguousFromObject(py_offs, NPY_LONG, 1, 1); // XXX: ssize_t typecode?
    a_y = (PyArrayObject *)PyArray_ContiguousFromObject(py_y, NPY_DOUBLE, 2, 2);

    if (a_a == NULL || a_offs == NULL || a_y == NULL) {
        Py_XDECREF(a_a);
        Py_XDECREF(a_offs);
        Py_XDECREF(a_y);
        return NULL;
    }

    // heavy lifting happens here, *in-place*
    fitpack::qr_reduce(
        static_cast<double *>(PyArray_DATA(a_a)), PyArray_DIM(a_a, 0), PyArray_DIM(a_a, 1), // a(m, nz), packed
        static_cast<ssize_t *>(PyArray_DATA(a_offs)),                                       // offset(m)
        nc,                                                                                 // dense would be a(m, nc)
        static_cast<double *>(PyArray_DATA(a_y)), PyArray_DIM(a_y, 1),                      // y(m, ydim2)
        startrow
    );

    Py_DECREF(a_offs);
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
