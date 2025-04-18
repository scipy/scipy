#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <iostream>
#include <climits>
#include <cassert>
#include "_batched_linalg.h"


/*
 * Basic checks: make sure `obj` is a LAPACK-compatible array
 */
static PyArrayObject*
ensure_array(PyObject *obj)
{
    if (!PyArray_CheckExact(obj)){
        return NULL;
    }

    PyArrayObject *arr = (PyArrayObject *)obj;
    int typenum = PyArray_TYPE(arr);

    bool dtype_ok = (typenum == NPY_FLOAT)
                     || (typenum == NPY_DOUBLE)
                     || (typenum == NPY_CFLOAT)
                     || (typenum == NPY_CDOUBLE);
    if(!dtype_ok || !PyArray_ISBEHAVED(arr)) {
        PyErr_SetString(PyExc_TypeError, "Expected a real or complex array.");
        return NULL;
    }

    return arr;
}


static PyObject*
py_inv(PyObject *self, PyObject *args)
{
    PyObject *py_a;
    int overwrite_a = 0;

    if(!PyArg_ParseTuple(args, "O|p", &py_a, &overwrite_a)) {
        return NULL;
    }

    PyArrayObject *a = ensure_array(py_a);
    if(a == NULL) {
        return NULL;
    }

    overwrite_a = overwrite_a
                  && PyArray_ISWRITEABLE(a)
                  && (PyArray_NDIM(a) == 2) && (PyArray_IS_F_CONTIGUOUS(a));

    PyArrayObject *a_inv;

    /*
        8. check LAPACK info, bail out : raise or warn or quiet (currently, the latter)
       11. Error handling. Make it talk to np.errstate?
     */

    npy_intp ndim = PyArray_NDIM(a);
    npy_intp *shape = PyArray_SHAPE(a);
    int typenum = PyArray_TYPE(a);

    if(!overwrite_a) {
        /* Allocate the result */
        a_inv = (PyArrayObject *)PyArray_SimpleNew(ndim, shape, typenum);
        if (!a_inv) {
            PyErr_NoMemory();
            return NULL;
        };
    } else {
        /* Reuse py_a's memory buffer. */
        a_inv = a;
        Py_INCREF(py_a);
    }

    long long status = -1;
    switch(typenum) {
    case(NPY_FLOAT) : status=inv_loop<float>(a, a_inv); break;
    case(NPY_DOUBLE) : status=inv_loop<double>(a, a_inv); break;
    case(NPY_CFLOAT) : status = inv_loop<npy_cfloat>(a, a_inv); break;
    case(NPY_CDOUBLE) : status = inv_loop<npy_cdouble>(a, a_inv); break;
    default:
        PyErr_SetString(PyExc_TypeError, "unreachable.");
        return NULL;
    }

    if (status == 0) {
        return (PyObject *)a_inv;
    }

    // some error occured
    if ((status == LLONG_MIN) && overwrite_a){
        // memory error; inverse was not computed
        Py_DECREF(py_a);
    }
    else {
        // XXX: LinalgError
        PyErr_Format(PyExc_ValueError,
                     "A singular matrix. LAPACK error code is %ll", status);
     }

     return NULL;
}




/////////////////////////////////////

static PyMethodDef BatchedLinalgMethods[] = {
    /* Batched linalg functions */
    {"inv", py_inv, METH_VARARGS, 
     "Invert a possibly batched matrix"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};



static struct PyModuleDef batched_linalg_module = {
    PyModuleDef_HEAD_INIT,
    "_batched_linalg",   /* name of module */
    NULL, //spam_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    BatchedLinalgMethods
};


PyMODINIT_FUNC
PyInit__batched_linalg(void)
{
    PyObject *module;

    import_array();

    module = PyModule_Create(&batched_linalg_module);
    if (module == NULL) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}
