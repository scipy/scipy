#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif /* PY_SSIZE_T_CLEAN */
#include "Python.h"
#include "numpy/arrayobject.h"
#include "_direct/direct.h"

static PyObject *
direct(PyObject *self, PyObject *args)
{
    PyObject *f, *f_args, *lb, *ub;
    int dimension, max_feval, max_iter, force_stop, disp;
    const double *lower_bounds, *upper_bounds;
    double minf, magic_eps, magic_eps_abs;
    double volume_reltol, sigma_reltol;
    double fglobal, fglobal_reltol;
    FILE *logfile = NULL;
    direct_algorithm algorithm;
    direct_return_code ret_code;

    if (!PyArg_ParseTuple(args, "OOOOidiiidddd", 
                          &f, &lb, &ub, &f_args, &disp, &magic_eps,
                          &max_feval, &max_iter, (int*) &algorithm,
                          &fglobal, &fglobal_reltol,
                          &volume_reltol, &sigma_reltol))
    {
        return NULL;
    }

    dimension = PyArray_DIMS((PyArrayObject*)lb)[0];
    double x[dimension + 1];
    PyObject *x_seq = PyList_New(dimension);
    lower_bounds = (double*)PyArray_DATA((PyArrayObject*)lb);
    upper_bounds = (double*)PyArray_DATA((PyArrayObject*)ub);
    magic_eps_abs = abs(magic_eps);
    volume_reltol /= 100.0;
    sigma_reltol /= 100.0;
    fglobal_reltol /= 100.0;
    force_stop = 0;
    direct_return_info info;

    ret_code =  direct_optimize(f, x, x_seq, f_args, dimension, lower_bounds,
                                upper_bounds, &minf, max_feval, max_iter,
                                magic_eps, magic_eps_abs, volume_reltol,
                                sigma_reltol, &force_stop, fglobal, fglobal_reltol,
                                logfile, algorithm, &info);
    PyObject* ret_py = Py_BuildValue("Odiii", x_seq, minf, (int) ret_code,
                                     info.numfunc, info.numiter);
    return ret_py;
}

/*
 * Standard Python module interface
 */

static PyMethodDef
DIRECTMethods[] = {
    {"direct", direct, METH_VARARGS, "DIRECT Optimization Algorithm"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_direct_lib",
    NULL,
    -1,
    DIRECTMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit__directmodule(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);

    return m;
}
