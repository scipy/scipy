/*
 * Last Change: Wed Jun 20 04:00 PM 2007 J
 *
 */
#include <Python.h>

#include <numpy/arrayobject.h>

#include "vq.h"

PyObject* compute_vq(PyObject*, PyObject*);

static PyMethodDef vqmethods [] = {
    {"vq", compute_vq, METH_VARARGS, "TODO docstring"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_vq(void) 
{
    Py_InitModule("_vq", vqmethods);
    import_array();
}

PyObject* compute_vq(PyObject* self, PyObject* args)
{
    PyObject *obs, *code, *out;
    PyArrayObject *obs_a, *code_a;
    PyArrayObject *index_a, *dist_a;
    int typenum1, typenum2;
    npy_intp nc, nd;
    npy_intp n, d;

    if ( !PyArg_ParseTuple(args, "OO", &obs, &code) ) {
        return NULL;
    }

    /* Check that obs and code both are arrays of same type, conformant
     * dimensions, etc...*/
    if (!(PyArray_Check(obs) && PyArray_Check(code))) {
		PyErr_Format(PyExc_ValueError,
			     "observation and code should be numpy arrays");
        return NULL;
    }

    typenum1 = PyArray_TYPE(obs);
    typenum2 = PyArray_TYPE(code);
    if (typenum1 != typenum1) {
		PyErr_Format(PyExc_ValueError,
			     "observation and code should have same type");
        return NULL;
    }
    obs_a = (PyArrayObject*)PyArray_FROM_OF(obs, 
                NPY_CONTIGUOUS | NPY_NOTSWAPPED | NPY_ALIGNED);
    if (obs_a == NULL) {
        return NULL;
    }

    code_a = (PyArrayObject*)PyArray_FROM_OF(code, 
                NPY_CONTIGUOUS | NPY_NOTSWAPPED | NPY_ALIGNED);
    if (code_a == NULL) {
        goto clean_obs_a;
    }

    if( !(obs_a->nd == code_a->nd)) {
		PyErr_Format(PyExc_ValueError,
			     "observation and code should have same shape");
        goto clean_code_a;
    }

    switch (obs_a->nd) {
        case 1:
            nd = 1;
            d = 1;
            n = PyArray_DIM(obs, 0);
            nc = PyArray_DIM(code, 0);
            break;
        case 2:
            nd = 2;
            n = PyArray_DIM(obs, 0);
            d = PyArray_DIM(obs, 1);
            nc = PyArray_DIM(code, 0);
            if (! (d == PyArray_DIM(code, 1)) ) {
                PyErr_Format(PyExc_ValueError,
                         "obs and code should have same number of "
                         " features (columns)");
                goto clean_code_a;
            }
            break;
        default:
            PyErr_Format(PyExc_ValueError,
                     "rank different than 1 or 2 are not supported");
            goto clean_code_a;
    }

    switch (PyArray_TYPE(obs)) {
        case NPY_FLOAT:
            dist_a = (PyArrayObject*)PyArray_EMPTY(1, &n, typenum1, 0);
            if (dist_a == NULL) {
                goto clean_code_a;
            }
            index_a = (PyArrayObject*)PyArray_EMPTY(1, &n, PyArray_INTP, 0);
            if (index_a == NULL) {
                goto clean_dist_a;
            }
            float_tvq((float*)obs_a->data, (float*)code_a->data, n, nc, d,
                    (npy_intp*)index_a->data, (float*)dist_a->data);
            break;
        case NPY_DOUBLE:
            dist_a = (PyArrayObject*)PyArray_EMPTY(1, &n, typenum1, 0);
            if (dist_a == NULL) {
                goto clean_code_a;
            }
            index_a = (PyArrayObject*)PyArray_EMPTY(1, &n, PyArray_INTP, 0);
            if (index_a == NULL) {
                goto clean_dist_a;
            }
            double_tvq((double*)obs_a->data, (double*)code_a->data, n, nc, d,
                    (npy_intp*)index_a->data, (double*)dist_a->data);
            break;
        default:
            PyErr_Format(PyExc_ValueError,
                     "type other than float or double not supported");
            goto clean_code_a;
    }

    /* Create output */
    out = PyTuple_New(2);
    if (out == NULL) {
        goto clean_index_a;
    }
    if (PyTuple_SetItem(out, 0, (PyObject*)index_a)) {
        goto clean_out;
    }
    if (PyTuple_SetItem(out, 1, (PyObject*)dist_a)) {
        goto clean_out;
    }

    /* Clean everything */
    Py_DECREF(code_a);
    Py_DECREF(obs_a);
    return out;

clean_out:
    Py_DECREF(out);
clean_dist_a:
    Py_DECREF(dist_a);
clean_index_a:
    Py_DECREF(index_a);
clean_code_a:
    Py_DECREF(code_a);
clean_obs_a:
    Py_DECREF(obs_a);
    return NULL;
}
