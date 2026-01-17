#include <Python.h>
#include <numpy/arrayobject.h>
#include "_directmodule.h"

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

    /*
     * Logging functionality is disable by default in the Python API. However
     * we still keep it here in C code, so that we don't have to re-implement
     * the whole functionality from scratch once we decide to expose it again
     * by default in Python API.
     */
    if (disp) {
        logfile = stdout;
    }

    dimension = PyArray_DIMS((PyArrayObject*)lb)[0];
    x = (double *) malloc(sizeof(double) * (dimension + 1));
    if (!x) {
        ret_code = DIRECT_OUT_OF_MEMORY;
    }
    PyObject *x_seq = PyList_New(dimension);
    lower_bounds = (double*)PyArray_DATA((PyArrayObject*)lb);
    upper_bounds = (double*)PyArray_DATA((PyArrayObject*)ub);
    magic_eps_abs = 0.0;
    force_stop = 0;
    direct_return_info info;

    PyObject *direct_ret = direct_optimize(f, x, x_seq, f_args, dimension, lower_bounds,
                         upper_bounds, &minf, max_feval, max_iter,
                         magic_eps, magic_eps_abs, volume_reltol,
                         sigma_reltol, &force_stop, fglobal, fglobal_reltol,
                         logfile, algorithm, &info, &ret_code, callback);
    if (!direct_ret) {
        Py_DECREF(x_seq);
        if (x)
            free(x);
        return NULL;
    }
    /* DECREF the return value from direct_optimize - we only needed it for error checking */
    Py_DECREF(direct_ret);
    PyObject* ret_py = Py_BuildValue("Odiii", x_seq, minf, (int) ret_code,
                                     info.numfunc, info.numiter);
    /* Py_BuildValue with "O" increments refcount. We need to DECREF our
       original reference since the tuple now owns it. */
    Py_DECREF(x_seq);
    if (x)
        free(x);
    return ret_py;
}

/*
 * Standard Python module interface
 */

static struct PyMethodDef direct_module_methods[] = {
    {"direct", direct, METH_VARARGS, "DIRECT Optimization Algorithm"},
    {NULL, NULL, 0, NULL}
};


static int module_exec(PyObject *module) {
    (void)module;  /* unused */

    if (_import_array() < 0) { return -1; }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif
    return 0;
}


static struct PyModuleDef_Slot direct_slots[] = {
    {Py_mod_exec, module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    /* signal that this module supports running without an active GIL */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_direct",
    .m_size = 0,
    .m_methods = direct_module_methods,
    .m_slots = direct_slots
};


PyMODINIT_FUNC
PyInit__direct(void)
{
    return PyModuleDef_Init(&moduledef);
}
