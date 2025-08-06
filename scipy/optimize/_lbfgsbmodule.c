#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>
#include "src/lbfgsb.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}

static PyObject* lbfgsb_error;

static char doc_setulb[] = "setulb(m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,lsave,isave,dsave,maxls,ln_task)";

static PyObject*
lbfgsb_setulb(PyObject *self, PyObject *args)
{
    int m, n, maxls;
    double f, factr, pgtol;
    int *nbd, *iwa, *isave, *lsave, *taskptr, *ln_taskptr;
    double *x, *l, *u, *g, *wa, *dsave;

    PyArrayObject *ap_x=NULL, *ap_l=NULL, *ap_u=NULL, *ap_g=NULL, *ap_nbd=NULL;
    PyArrayObject *ap_wa=NULL, *ap_iwa=NULL, *ap_task=NULL, *ap_ln_task=NULL;
    PyArrayObject *ap_lsave=NULL, *ap_isave=NULL, *ap_dsave=NULL;

    // Check argument types
    // m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,lsave,isave,dsave,maxls,ln_task
    if (!(PyArg_ParseTuple(args,
                           "iO!O!O!O!dO!ddO!O!O!O!O!O!iO!",
                           &m,                                      // i
                           &PyArray_Type, (PyObject **)&ap_x,       // O!
                           &PyArray_Type, (PyObject **)&ap_l,       // O!
                           &PyArray_Type, (PyObject **)&ap_u,       // O!
                           &PyArray_Type, (PyObject **)&ap_nbd,     // O!
                           &f,                                      // d
                           &PyArray_Type, (PyObject **)&ap_g,       // O!
                           &factr,                                  // d
                           &pgtol,                                  // d
                           &PyArray_Type, (PyObject **)&ap_wa,      // O!
                           &PyArray_Type, (PyObject **)&ap_iwa,     // O!
                           &PyArray_Type, (PyObject **)&ap_task,    // O!
                           &PyArray_Type, (PyObject **)&ap_lsave,   // O!
                           &PyArray_Type, (PyObject **)&ap_isave,   // O!
                           &PyArray_Type, (PyObject **)&ap_dsave,   // O!
                           &maxls,                                  // i
                           &PyArray_Type, (PyObject **)&ap_ln_task  // O!
                           ))) { return NULL; }

    // Check if arrays are contiguous and with the right dtype.
    // All arrays are 1D hence both F and C contiguous flags are set to True.
    if (!(PyArray_IS_C_CONTIGUOUS(ap_x)) || (PyArray_TYPE(ap_x) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (x) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_l)) || (PyArray_TYPE(ap_l) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (l) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_u)) || (PyArray_TYPE(ap_u) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (u) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_wa)) || (PyArray_TYPE(ap_wa) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (wa) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_iwa)) || (PyArray_TYPE(ap_iwa) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (iwa) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_dsave)) || (PyArray_TYPE(ap_dsave) != NPY_FLOAT64))
    {
        PYERR(lbfgsb_error, " Argument (dsave) must be a contiguous array of type float64.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_nbd)) || (PyArray_TYPE(ap_nbd) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (nbd) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_task)) || (PyArray_TYPE(ap_task) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (task) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_lsave)) || (PyArray_TYPE(ap_lsave) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (lsave) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_isave)) || (PyArray_TYPE(ap_isave) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (isave) must be a contiguous array of type int32.");
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_ln_task)) || (PyArray_TYPE(ap_ln_task) != NPY_INT32))
    {
        PYERR(lbfgsb_error, " Argument (ln_task) must be a contiguous array of type int32.");
    }
    n = PyArray_DIMS(ap_x)[0];

    x = (double *)PyArray_DATA(ap_x);
    l = (double *)PyArray_DATA(ap_l);
    u = (double *)PyArray_DATA(ap_u);
    nbd = (int *)PyArray_DATA(ap_nbd);
    g = (double *)PyArray_DATA(ap_g);
    wa = (double *)PyArray_DATA(ap_wa);
    iwa = (int *)PyArray_DATA(ap_iwa);
    lsave = (int *)PyArray_DATA(ap_lsave);
    isave = (int *)PyArray_DATA(ap_isave);
    dsave = (double *)PyArray_DATA(ap_dsave);
    taskptr = (int *)PyArray_DATA(ap_task);
    ln_taskptr = (int *)PyArray_DATA(ap_ln_task);

    // Make the function call
    setulb(n, m, x, l, u, nbd, &f, g, factr, pgtol, wa, iwa, taskptr, lsave,
           isave, dsave, maxls, ln_taskptr);

    Py_RETURN_NONE;
}


static struct PyMethodDef lbfgsb_module_methods[] = {
  {"setulb", lbfgsb_setulb, METH_VARARGS, doc_setulb},
  {NULL,     NULL,          0,            NULL}
};


static int module_exec(PyObject *module) {

    if (_import_array() < 0) { return -1; }

    lbfgsb_error = PyErr_NewException("_lbfgsb.error", NULL, NULL);
    if (lbfgsb_error == NULL) { return -1; }

    if (PyModule_AddObject(module, "error", lbfgsb_error) < 0) {
        Py_DECREF(lbfgsb_error);
        return -1;
    }
#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif
    return 0;
}


static struct PyModuleDef_Slot _lbfgsb_slots[] = {
    {Py_mod_exec, module_exec},
#if PY_VERSION_HEX >= 0x030c00f0  // Python 3.12+
    // signal that this module can be imported in isolated subinterpreters
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    // signal that this module supports running without an active GIL
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_lbfgsb",
    .m_size = 0,
    .m_methods = lbfgsb_module_methods,
    .m_slots = _lbfgsb_slots
};


PyMODINIT_FUNC
PyInit__lbfgsb(void)
{
    return PyModuleDef_Init(&moduledef);
}
