#include "Python.h"
#include <stdlib.h>

#include "interpolate.h"
#include "numpy/arrayobject.h"

using namespace std;

extern "C" {

static PyObject* linear_method(PyObject*self, PyObject* args, PyObject* kywds)
{
    static char *kwlist[] = {"x","y","new_x","new_y", NULL};
    PyObject *py_x, *py_y, *py_new_x, *py_new_y;
    py_x = py_y = py_new_x = py_new_y = NULL;
    PyObject *arr_x, *arr_y, *arr_new_x, *arr_new_y;
    arr_x = arr_y = arr_new_x = arr_new_y = NULL;

    if(!PyArg_ParseTupleAndKeywords(args,kywds,"OOOO:linear_dddd",kwlist,&py_x, &py_y, &py_new_x, &py_new_y))
       return NULL;
    arr_x = PyArray_FROMANY(py_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_x) {
        PyErr_SetString(PyExc_ValueError, "x must be a 1-D array of floats");
        goto fail;
    }
    arr_y = PyArray_FROMANY(py_y, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_y) {
        PyErr_SetString(PyExc_ValueError, "y must be a 1-D array of floats");
        goto fail;
    }
    arr_new_x = PyArray_FROMANY(py_new_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_new_x) {
        PyErr_SetString(PyExc_ValueError, "new_x must be a 1-D array of floats");
        goto fail;
    }
    arr_new_y = PyArray_FROMANY(py_new_y, PyArray_DOUBLE, 1, 1, NPY_INOUT_ARRAY);
    if (!arr_new_y) {
        PyErr_SetString(PyExc_ValueError, "new_y must be a 1-D array of floats");
        goto fail;
    }

    linear((double*)PyArray_DATA(arr_x), (double*)PyArray_DATA(arr_y),
            PyArray_DIM(arr_x,0), (double*)PyArray_DATA(arr_new_x),
            (double*)PyArray_DATA(arr_new_y), PyArray_DIM(arr_new_x,0));

    Py_DECREF(arr_x);
    Py_DECREF(arr_y);
    Py_DECREF(arr_new_x);
    Py_DECREF(arr_new_y);

    Py_RETURN_NONE;

fail:
    Py_XDECREF(arr_x);
    Py_XDECREF(arr_y);
    Py_XDECREF(arr_new_x);
    Py_XDECREF(arr_new_y);
    return NULL;
}

static PyObject* loginterp_method(PyObject*self, PyObject* args, PyObject* kywds)
{
    static char *kwlist[] = {"x","y","new_x","new_y", NULL};
    PyObject *py_x, *py_y, *py_new_x, *py_new_y;
    py_x = py_y = py_new_x = py_new_y = NULL;
    PyObject *arr_x, *arr_y, *arr_new_x, *arr_new_y;
    arr_x = arr_y = arr_new_x = arr_new_y = NULL;

    if(!PyArg_ParseTupleAndKeywords(args,kywds,"OOOO:loginterp_dddd",kwlist,&py_x, &py_y, &py_new_x, &py_new_y))
       return NULL;
    arr_x = PyArray_FROMANY(py_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_x) {
        PyErr_SetString(PyExc_ValueError, "x must be a 1-D array of floats");
        goto fail;
    }
    arr_y = PyArray_FROMANY(py_y, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_y) {
        PyErr_SetString(PyExc_ValueError, "y must be a 1-D array of floats");
        goto fail;
    }
    arr_new_x = PyArray_FROMANY(py_new_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_new_x) {
        PyErr_SetString(PyExc_ValueError, "new_x must be a 1-D array of floats");
        goto fail;
    }
    arr_new_y = PyArray_FROMANY(py_new_y, PyArray_DOUBLE, 1, 1, NPY_INOUT_ARRAY);
    if (!arr_new_y) {
        PyErr_SetString(PyExc_ValueError, "new_y must be a 1-D array of floats");
        goto fail;
    }

    loginterp((double*)PyArray_DATA(arr_x), (double*)PyArray_DATA(arr_y),
            PyArray_DIM(arr_x,0), (double*)PyArray_DATA(arr_new_x),
            (double*)PyArray_DATA(arr_new_y), PyArray_DIM(arr_new_x,0));

    Py_DECREF(arr_x);
    Py_DECREF(arr_y);
    Py_DECREF(arr_new_x);
    Py_DECREF(arr_new_y);

    Py_RETURN_NONE;

fail:
    Py_XDECREF(arr_x);
    Py_XDECREF(arr_y);
    Py_XDECREF(arr_new_x);
    Py_XDECREF(arr_new_y);
    return NULL;
}

static PyObject* window_average_method(PyObject*self, PyObject* args, PyObject* kywds)
{
    static char *kwlist[] = {"x","y","new_x","new_y", NULL};
    PyObject *py_x, *py_y, *py_new_x, *py_new_y;
    py_x = py_y = py_new_x = py_new_y = NULL;
    PyObject *arr_x, *arr_y, *arr_new_x, *arr_new_y;
    arr_x = arr_y = arr_new_x = arr_new_y = NULL;
    double width;

    if(!PyArg_ParseTupleAndKeywords(args,kywds,"OOOOd:loginterp_dddd",kwlist,&py_x, &py_y, &py_new_x, &py_new_y, &width))
       return NULL;
    arr_x = PyArray_FROMANY(py_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_x) {
        PyErr_SetString(PyExc_ValueError, "x must be a 1-D array of floats");
        goto fail;
    }
    arr_y = PyArray_FROMANY(py_y, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_y) {
        PyErr_SetString(PyExc_ValueError, "y must be a 1-D array of floats");
        goto fail;
    }
    arr_new_x = PyArray_FROMANY(py_new_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_new_x) {
        PyErr_SetString(PyExc_ValueError, "new_x must be a 1-D array of floats");
        goto fail;
    }
    arr_new_y = PyArray_FROMANY(py_new_y, PyArray_DOUBLE, 1, 1, NPY_INOUT_ARRAY);
    if (!arr_new_y) {
        PyErr_SetString(PyExc_ValueError, "new_y must be a 1-D array of floats");
        goto fail;
    }

    window_average((double*)PyArray_DATA(arr_x), (double*)PyArray_DATA(arr_y),
            PyArray_DIM(arr_x,0), (double*)PyArray_DATA(arr_new_x),
            (double*)PyArray_DATA(arr_new_y), PyArray_DIM(arr_new_x,0), width);

    Py_DECREF(arr_x);
    Py_DECREF(arr_y);
    Py_DECREF(arr_new_x);
    Py_DECREF(arr_new_y);

    Py_RETURN_NONE;

fail:
    Py_XDECREF(arr_x);
    Py_XDECREF(arr_y);
    Py_XDECREF(arr_new_x);
    Py_XDECREF(arr_new_y);
    return NULL;
}

static PyObject* block_average_above_method(PyObject*self, PyObject* args, PyObject* kywds)
{
    static char *kwlist[] = {"x","y","new_x","new_y", NULL};
    PyObject *py_x, *py_y, *py_new_x, *py_new_y;
    py_x = py_y = py_new_x = py_new_y = NULL;
    PyObject *arr_x, *arr_y, *arr_new_x, *arr_new_y;
    arr_x = arr_y = arr_new_x = arr_new_y = NULL;

    if(!PyArg_ParseTupleAndKeywords(args,kywds,"OOOO:loginterp_dddd",kwlist,&py_x, &py_y, &py_new_x, &py_new_y))
       return NULL;
    arr_x = PyArray_FROMANY(py_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_x) {
        PyErr_SetString(PyExc_ValueError, "x must be a 1-D array of floats");
        goto fail;
    }
    arr_y = PyArray_FROMANY(py_y, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_y) {
        PyErr_SetString(PyExc_ValueError, "y must be a 1-D array of floats");
        goto fail;
    }
    arr_new_x = PyArray_FROMANY(py_new_x, PyArray_DOUBLE, 1, 1, NPY_IN_ARRAY);
    if (!arr_new_x) {
        PyErr_SetString(PyExc_ValueError, "new_x must be a 1-D array of floats");
        goto fail;
    }
    arr_new_y = PyArray_FROMANY(py_new_y, PyArray_DOUBLE, 1, 1, NPY_INOUT_ARRAY);
    if (!arr_new_y) {
        PyErr_SetString(PyExc_ValueError, "new_y must be a 1-D array of floats");
        goto fail;
    }

    block_average_above((double*)PyArray_DATA(arr_x), (double*)PyArray_DATA(arr_y),
            PyArray_DIM(arr_x,0), (double*)PyArray_DATA(arr_new_x),
            (double*)PyArray_DATA(arr_new_y), PyArray_DIM(arr_new_x,0));

    Py_DECREF(arr_x);
    Py_DECREF(arr_y);
    Py_DECREF(arr_new_x);
    Py_DECREF(arr_new_y);

    Py_RETURN_NONE;

fail:
    Py_XDECREF(arr_x);
    Py_XDECREF(arr_y);
    Py_XDECREF(arr_new_x);
    Py_XDECREF(arr_new_y);
    return NULL;
}

static PyMethodDef interpolate_methods[] = {
    {"linear_dddd", (PyCFunction)linear_method, METH_VARARGS|METH_KEYWORDS,
        ""},
    {"loginterp_dddd", (PyCFunction)loginterp_method, METH_VARARGS|METH_KEYWORDS,
        ""},
    {"window_average_ddddd", (PyCFunction)window_average_method, METH_VARARGS|METH_KEYWORDS,
        ""},
    {"block_average_above_dddd", (PyCFunction)block_average_above_method, METH_VARARGS|METH_KEYWORDS,
        ""},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC init_interpolate(void)
{
    PyObject* m;
    m = Py_InitModule3("_interpolate", interpolate_methods, 
        "A few interpolation routines.\n"
        );
    if (m == NULL)
        return;
    import_array();
}

} // extern "C"
