
/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modifications by Travis Oliphant to separate Python code from C
   routines */

#include "Python.h"
#include <setjmp.h>
#include "Zeros/zeros.h"

/*
 * Caller entry point functions
 */

#ifdef PYPY_VERSION
    /*
     * As described in http://doc.pypy.org/en/latest/cpython_differences.html#c-api-differences,
     * "assignment to a PyTupleObject is not supported after the tuple is used internally,
     * even by another C-API function call."
     */
    #define PyArgs(Operation) PyList_##Operation
#else
    /* Using a list in CPython raises "TypeError: argument list must be a tuple" */
    #define PyArgs(Operation) PyTuple_##Operation
#endif

typedef struct {
    PyObject *function;
    PyObject *xargs;
    jmp_buf env;
} scipy_zeros_parameters;



static double
scipy_zeros_functions_func(double x, void *params)
{
    scipy_zeros_parameters *myparams = params;
    PyObject *args, *xargs, *item, *f, *retval=NULL;
    Py_ssize_t i, len;
    double val;

    xargs = myparams->xargs;
    /* Need to create a new 'args' tuple on each call in case 'f' is
       stateful and keeps references to it (e.g. functools.lru_cache) */
    len = PyTuple_Size(xargs);
    /* Make room for the double as first argument */
    args = PyArgs(New)(len + 1);
    if (args == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to allocate arguments");
        longjmp(myparams->env, 1);
    }
    PyArgs(SET_ITEM)(args, 0, Py_BuildValue("d", x));
    for (i = 0; i < len; i++) {
        item = PyTuple_GetItem(xargs, i);
        if (item == NULL) {
            Py_DECREF(args);
            longjmp(myparams->env, 1);
        }
        Py_INCREF(item);
        PyArgs(SET_ITEM)(args, i+1, item);
    }

    f = myparams->function;
    retval = PyObject_CallObject(f,args);
    Py_DECREF(args);
    if (retval == NULL) {
        longjmp(myparams->env, 1);
    }
    val = PyFloat_AsDouble(retval);
    Py_XDECREF(retval);
    return val;
}


/*
 * Helper function that calls a Python function with extended arguments
 */

static PyObject *
call_solver(solver_type solver, PyObject *self, PyObject *args)
{
    double a, b, xtol, rtol, zero;
    int iter, fulloutput, disp=1, flag=0;
    scipy_zeros_parameters params;
    scipy_zeros_info solver_stats;
    PyObject *f, *xargs;

    if (!PyArg_ParseTuple(args, "OddddiOi|i",
                &f, &a, &b, &xtol, &rtol, &iter, &xargs, &fulloutput, &disp)) {
        PyErr_SetString(PyExc_RuntimeError, "Unable to parse arguments");
        return NULL;
    }
    if (xtol < 0) {
        PyErr_SetString(PyExc_ValueError, "xtol must be >= 0");
        return NULL;
    }
    if (iter < 0) {
        PyErr_SetString(PyExc_ValueError, "maxiter must be >= 0");
        return NULL;
    }

    params.function = f;
    params.xargs = xargs;

    if (!setjmp(params.env)) {
        /* direct return */
        solver_stats.error_num = 0;
        zero = solver(scipy_zeros_functions_func, a, b, xtol, rtol,
                      iter, (void*)&params, &solver_stats);
    } else {
        /* error return from Python function */
        return NULL;
    }

    if (solver_stats.error_num != CONVERGED) {
        if (solver_stats.error_num == SIGNERR) {
            PyErr_SetString(PyExc_ValueError,
                    "f(a) and f(b) must have different signs");
            return NULL;
        }
        if (solver_stats.error_num == CONVERR) {
            if (disp) {
                char msg[100];
                PyOS_snprintf(msg, sizeof(msg),
                        "Failed to converge after %d iterations.",
                        solver_stats.iterations);
                PyErr_SetString(PyExc_RuntimeError, msg);
                return NULL;
            }
            flag = CONVERR;
        }
    }
    else {
        flag = CONVERGED;
    }
    if (fulloutput) {
        return Py_BuildValue("diii",
                zero, solver_stats.funcalls, solver_stats.iterations, flag);
    }
    else {
        return Py_BuildValue("d", zero);
    }
}

/*
 * These routines interface with the solvers through call_solver
 */

static PyObject *
_bisect(PyObject *self, PyObject *args)
{
        return call_solver(bisect,self,args);
}

static PyObject *
_ridder(PyObject *self, PyObject *args)
{
        return call_solver(ridder,self,args);
}

static PyObject *
_brenth(PyObject *self, PyObject *args)
{
        return call_solver(brenth,self,args);
}

static PyObject *
_brentq(PyObject *self, PyObject *args)
{
        return call_solver(brentq,self,args);
}

/*
 * Standard Python module interface
 */

static PyMethodDef
Zerosmethods[] = {
	{"_bisect", _bisect, METH_VARARGS, "a"},
	{"_ridder", _ridder, METH_VARARGS, "a"},
	{"_brenth", _brenth, METH_VARARGS, "a"},
	{"_brentq", _brentq, METH_VARARGS, "a"},
	{NULL, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_zeros",
    NULL,
    -1,
    Zerosmethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__zeros(void)
{
    PyObject *module;

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return module;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}
