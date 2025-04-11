#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <iostream>
#include <limits>
#include <cassert>
#include "_batched_linalg.h"


static PyObject*
py_inv(PyObject *self, PyObject *args)
{

    PyObject *py_a;
    int overwrite_a = 0;

    if(!PyArg_ParseTuple(args, "O|p", &py_a, &overwrite_a)) {
        return NULL;
    }

    PyArrayObject *a = (PyArrayObject *)py_a;

    /*
        8. check LAPACK info, bail out : raise or warn or quiet (currently, the latter)
       10. overwrite_a
       11. checks/asserts: m, n > 0, lda >= n etc

     */

    npy_intp ndim = PyArray_NDIM(a);
    npy_intp *shape = PyArray_SHAPE(a);
    int typenum = PyArray_TYPE(a);

    assert(ndim >= 3);  // XXX needs to make sure of this on the python side

    /* Allocate the result */
    PyArrayObject *a_inv = (PyArrayObject *)PyArray_SimpleNew(ndim, shape, typenum);
    if (!a_inv) {
        PyErr_NoMemory();
        return NULL;
    }

    int status = -1;
    switch(typenum) {
    case(NPY_FLOAT) : status=inv_loop<float>(a, a_inv); break;
    case(NPY_DOUBLE) : status=inv_loop<double>(a, a_inv); break;
    case(NPY_CFLOAT) : status = inv_loop<npy_cfloat>(a, a_inv); break;
    case(NPY_CDOUBLE) : status = inv_loop<npy_cdouble>(a, a_inv); break;
    default:
        PyErr_SetString(PyExc_TypeError, "unknown typenum");
        return NULL;
    }

    if (status != 0){
        PyErr_SetString(PyExc_RuntimeError, "Something went wrong.");
        return NULL;
    }

    return (PyObject *)a_inv;

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
