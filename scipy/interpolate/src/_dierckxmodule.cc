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



/////////////////////////////////////

static PyMethodDef DierckxMethods[] = {
    //...
    {"fpknot", py_fpknot, METH_VARARGS,
     "fpknot replacement"},
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
