#ifndef __SLSQPLIB_H
#define __SLSQPLIB_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>
#include "__nnls.h"


#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}
static PyObject* slsqp_error;

static char doc_nnls[] = ("Compute the nonnegative least squares solution.\n\n"
                           "    x, info = nnls(A)\n\n");

static PyObject* nnls(PyObject *dummy, PyObject *args);

static struct PyMethodDef slsqplib_module_methods[] = {
  {"nnls", nnls, METH_VARARGS, doc_nnls},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_slsqplib",
    NULL,
    -1,
    slsqplib_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__slsqplib(void)
{
    PyObject *module, *mdict;

    import_array();

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }

    mdict = PyModule_GetDict(module);
    if (mdict == NULL) {
        return NULL;
    }
    slsqp_error = PyErr_NewException("slsqplib.error", NULL, NULL);
    if (slsqp_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", slsqp_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}


#endif // __SLSQPLIB_H