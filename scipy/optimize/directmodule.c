#include "Python.h"
#include "numpy/arrayobject.h"
#include "directmodule.h"
#define MY_ALLOC(p, t, n, ret_code) p = (t *) malloc(sizeof(t) * (n)); \
                          if (!(p)) { *ret_code = DIRECT_OUT_OF_MEMORY; }
#define MY_FREE(p) if (p) free(p)

static PyObject *
direct(PyObject *self, PyObject *args)
{
    PyObject *f, *f_args, *lb, *ub, *callback;
    int dimension, max_feval, max_iter, force_stop, disp;
    const double *lower_bounds, *upper_bounds;
    double minf, magic_eps, magic_eps_abs, *x;
    double volume_reltol, sigma_reltol;
    double fglobal, fglobal_reltol;
    FILE *logfile = NULL;
    direct_algorithm algorithm;
    direct_return_code ret_code;

    if (!PyArg_ParseTuple(args, "OOOOidiiiddddO",
                          &f, &lb, &ub, &f_args, &disp, &magic_eps,
                          &max_feval, &max_iter, (int*) &algorithm,
                          &fglobal, &fglobal_reltol,
                          &volume_reltol, &sigma_reltol, &callback))
    {
        return NULL;
    }

    if (disp) {
        logfile = stdout;
    }

    dimension = PyArray_DIMS((PyArrayObject*)lb)[0];
    MY_ALLOC(x, double, dimension + 1, &ret_code);
    PyObject *x_seq = PyList_New(dimension);
    lower_bounds = (double*)PyArray_DATA((PyArrayObject*)lb);
    upper_bounds = (double*)PyArray_DATA((PyArrayObject*)ub);
    magic_eps_abs = 0.0;
    force_stop = 0;
    direct_return_info info;

    if (!direct_optimize(f, x, x_seq, f_args, dimension, lower_bounds,
                         upper_bounds, &minf, max_feval, max_iter,
                         magic_eps, magic_eps_abs, volume_reltol,
                         sigma_reltol, &force_stop, fglobal, fglobal_reltol,
                         logfile, algorithm, &info, &ret_code, callback)) {
        MY_FREE(x);
        return NULL;
    }
    PyObject* ret_py = Py_BuildValue("Odiii", x_seq, minf, (int) ret_code,
                                     info.numfunc, info.numiter);
    MY_FREE(x);
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
    import_array();

    return m;
}
