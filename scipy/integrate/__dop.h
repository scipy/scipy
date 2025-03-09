#ifndef __DOP_H
#define __DOP_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "ccallback.h"
#include <math.h>


#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

// dopri_fcn(n, x, y, f)
// dopri_solout(nr, xold, x, y, n, con, icomp, nd)
typedef void dopri_fcn(int, double, double*, double*);
typedef int dopri_solout(int, double, double, double*, int, double*, int*, int);


static PyObject* dop_error;

static char doc_dopri853[] = "[result,abserr,infodict,ier] = _qagse(fun, a, b, | args, full_output, epsabs, epsrel, limit)";
static char doc_dopri5[] = "x,y,iwork,idid = dopri5(fcn,x,y,xend,rtol,atol,solout,iout,work,iwork,[fcn_extra_args,overwrite_y,solout_extra_args])";

void dopri853(const int n, dopri_fcn* fcn, double x, double* y, double* xend, double* rtol,
              double* atol, const int itol, dopri_solout* solout, const int iout, double* work,
              int* iwork, double* rpar, int* ipar, int* ierr);
void dopri5(int n, dopri_fcn* fcn, double x, double* y, double* xend, double* rtol,
            double* atol, const int itol, dopri_solout* solout, const int iout,
            double* work, int* iwork, int* ierr);


/*
 *
 * Handle the callback handler here for Python - C bridge.
 *
 */

static struct PyMethodDef dop_module_methods[] = {
  {"_dop853", dop_dopri853, METH_VARARGS, doc_dopri853},
  {"_dopri5", dop_dopri5, METH_VARARGS, doc_dopri5},
  {NULL    , NULL          , 0           , NULL     }
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_dop",
    NULL,
    -1,
    dop_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__dop(void)
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

    dop_error = PyErr_NewException ("_dop.error", NULL, NULL);
    if (dop_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", dop_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif
    return module;
}

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif