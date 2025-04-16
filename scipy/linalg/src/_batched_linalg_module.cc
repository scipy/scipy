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
    if (!PyArray_CheckExact(py_a)){
        return NULL;
    }
    // XXX: check NDIM, dtypes


    PyArrayObject *a = (PyArrayObject *)py_a;
    PyArrayObject *a_inv;

    /*
        8. check LAPACK info, bail out : raise or warn or quiet (currently, the latter)
        9. checks/asserts: m, n > 0, lda >= n etc
       10. overwrite_a related checks: flags etc, move from Python to here.
       11. Error handling. Make it talk to np.errstate?

     */

    npy_intp ndim = PyArray_NDIM(a);
    npy_intp *shape = PyArray_SHAPE(a);
    int typenum = PyArray_TYPE(a);

    bool dtype_ok = ((typenum == NPY_FLOAT)
                      || (typenum == NPY_DOUBLE)
                      || (typenum == NPY_CFLOAT)
                      || (typenum == NPY_CDOUBLE));
    if(!dtype_ok){
        PyErr_SetString(PyExc_TypeError, "Incompatible dtype.");
        return NULL;
    }

    if(!overwrite_a) {
        /* Allocate the result */
        a_inv = (PyArrayObject *)PyArray_SimpleNew(ndim, shape, typenum);
        if (!a_inv) {
            PyErr_NoMemory();
            return NULL;
        };
    } else {
        // XXX: check flags/strides
        a_inv = a;
    }

    int status = -1;
    switch(typenum) {
    case(NPY_FLOAT) : status=inv_loop<float>(a, a_inv); break;
    case(NPY_DOUBLE) : status=inv_loop<double>(a, a_inv); break;
    case(NPY_CFLOAT) : status = inv_loop<npy_cfloat>(a, a_inv); break;
    case(NPY_CDOUBLE) : status = inv_loop<npy_cdouble>(a, a_inv); break;
    default:
        PyErr_SetString(PyExc_TypeError, "unreachable.");
        return NULL;
    }

    if (status == -101){
        // memory error
        if(!overwrite_a) {
            // hey garbage collector
            Py_DECREF(a_inv);
        }
        else {
            // XXX: LinalgError
            PyErr_Format(PyExc_ValueError, "A singular matrix. LAPACK error code is %d", status);
        }
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
