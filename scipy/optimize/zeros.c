
/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modifications by Travis Oliphant to separate Python code from C
   routines */

#include "Python.h"
#include <setjmp.h>
#include <numpy/npy_math.h>

#include <stdio.h>

#if PY_VERSION_HEX >= 0x03000000
    #define PyString_FromString PyBytes_FromString
    #define PyString_Concat PyBytes_Concat
    #define PyString_AsString PyBytes_AsString
    #define PyInt_FromLong PyLong_FromLong
    #define PyInt_AsLong PyLong_AsLong
#endif

#ifndef FALSE
#define FALSE 0
#endif /* FALSE */
#ifndef TRUE
#define TRUE 1
#endif /* TRUE */

#include "ccallback.h"
#include "Zeros/zeros.h"

/*
 * Caller entry point functions
 */

typedef enum {
    CB_D = 1,         // f(x)
    CB_D_D = 2,       // f(x, double)
    CB_D_VOIDSTAR = 100,     // f(x, void *) or f(x, const void *)
    CB_D_USERDATA = 1001     // f(x, void *) or f(x, const void *) using UserData
} zeros_signature_t;

static const ccallback_signature_t signatures_with_args[] = {
    {"double (double, double)", CB_D_D},
    {"double (double, void const *)", CB_D_VOIDSTAR},
    {"double (double, void *)", CB_D_VOIDSTAR},
    {NULL}
};

static const ccallback_signature_t signatures_no_args[] = {
    {"double (double)", CB_D},
    {"double (double, void const *)", CB_D_USERDATA},
    {"double (double, void *)", CB_D_USERDATA},
    {NULL}
};


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
    PyObject *args;
    jmp_buf env;
} scipy_zeros_parameters;

typedef struct {
    double a;
    int n;
} _x_to_the_n_minus_a_data_t;


static double
scipy_zeros_functions_lowlevelfunc(double x, void *func_data)
{
    double result;
    zeros_signature_t sigval;
    ccallback_t *callback = (ccallback_t*)func_data;
    assert(callback->c_function != NULL);
    assert(callback->signature != NULL);
    assert(callback->signature->value == CB_D || callback->info_p != NULL);
    sigval = callback->signature->value;
    switch (sigval) {
      case CB_D:
        result = ((double(*)(double))callback->c_function)(x);
        break;

      case CB_D_D:
        result = ((double(*)(double, double))callback->c_function)(
            x, *(double *)(callback->info_p));
        break;

      case CB_D_VOIDSTAR:
        result = ((double(*)(double, const void *))callback->c_function)(
            x, *(const void **)(callback->info_p));
        break;

      case CB_D_USERDATA:
        result = ((double(*)(double, const void *))callback->c_function)(
            x, (const void *)(callback->user_data));
        break;

      default:
         result = NPY_NAN;
         break;
    }
    return result;
}

static double
scipy_zeros_functions_func(double x, void *params)
{
    scipy_zeros_parameters *myparams = params;
    PyObject *args, *f, *retval=NULL;
    double val;

    args = myparams->args;
    f = myparams->function;
    PyArgs(SetItem)(args, 0, Py_BuildValue("d",x));
    retval = PyObject_CallObject(f,args);
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

static int
init_func_extra_args(ccallback_t *callback, PyObject *extra_arguments)
{
    double *p;
    Py_ssize_t i, num_xargs;

    callback->info_p = NULL;

    if (!PyTuple_Check(extra_arguments)) {
        return -1;
    }
    num_xargs = PyTuple_GET_SIZE(extra_arguments);

    p = (double *)calloc(num_xargs, sizeof(double));
    if (p == NULL) {
        PyErr_SetString(PyExc_MemoryError, "failed to allocate memory");
        return -1;
    }

    for (i = 0; i < num_xargs; ++i) {
        PyObject *item = PyTuple_GET_ITEM(extra_arguments, i);
        p[i] = PyFloat_AsDouble(item);
        if (PyErr_Occurred()) {
            free(p);
            return -1;
        }
    }

    callback->info_p = (void *)p;
    callback->info = (int)num_xargs;
    return 0;
}


static int fill_in_ccallback(PyObject *f, PyObject *xargs, ccallback_t *pcallback)
{
    int use_ccallback = FALSE;
    int bHasArgs = (xargs && PyTuple_Check(xargs) && PyTuple_GET_SIZE(xargs) > 0);

    /* f was filled in previously by PyArg_ParseTuple */
    const ccallback_signature_t *sigs =  (
        bHasArgs ? signatures_with_args : signatures_no_args);
    if (ccallback_prepare(pcallback, (ccallback_signature_t *)sigs,
                          f, CCALLBACK_DEFAULTS)) {
        return FALSE;
    }

    pcallback->info_p = NULL;
    if (pcallback->signature == NULL) {
        /* pure-Python */
        return FALSE;
    }
    switch(pcallback->signature->value) {
        case CB_D:
            /* extra_arguments are just ignored */
            pcallback->info_p = NULL;
            use_ccallback = TRUE;
            break;

        case CB_D_D:
            if (!init_func_extra_args(pcallback, xargs)) {
                /* Check for a single argument */
                if (pcallback->info != 1) {
                    PyErr_SetString(PyExc_ValueError,
                        "Wrong number of extra arguments.");
                } else {
                    use_ccallback = TRUE;
                }
            } else {
                PyErr_SetString(PyExc_ValueError,
                    "Unable to parse xargs as a double.");
            }
            break;

        case CB_D_VOIDSTAR:
            if (xargs && PyTuple_Check(xargs) && PyTuple_GET_SIZE(xargs) == 1) {
                PyObject *item = PyTuple_GET_ITEM(xargs, 0);
                // Expecting an object created by ctypes.pointer (E.g. c_void_p gives a CDataObject)
                // [ ctypes.byref returns CArgObject, more difficult to deal with.  It has a "_obj" member.]
                // Assume it is a CDataObject
                if (PyObject_HasAttrString(item, "value") ) {
                    void **p = NULL;
                    PyObject * p2 = PyObject_GetAttrString(item, "value");
                    assert(p2);
                    p = malloc(sizeof(void *));
                    *p = PyLong_AsVoidPtr(p2);
                    Py_DECREF(p2);
                    pcallback->info = 1;
                    pcallback->info_p = p;
                    use_ccallback = TRUE;
                } else {
                    PyErr_SetString(PyExc_ValueError,
                        "Unable to Extract value from xargs.");
                }
            } else {
                PyErr_SetString(PyExc_ValueError,
                    "Wrong number of extra arguments for void *.");
            }
            break;

        case CB_D_USERDATA:
            // Is the data passed in, or stored in user_data?
            if (xargs && PyTuple_Check(xargs) && PyTuple_GET_SIZE(xargs) > 0) {
                PyErr_SetString(PyExc_ValueError,
                        "Wrong number of extra arguments for user_data.");
            } else {
                use_ccallback = TRUE;
            }
            break;

        default:
        // FUTURE: Set an error string?
            break;
    }
    return use_ccallback;
}

static PyObject *
call_solver(solver_type solver, PyObject *self, PyObject *args)
{
    double a, b, xtol, rtol, zero;
    Py_ssize_t len;
    int iter, i, fulloutput, disp=1, flag=0;
    scipy_zeros_parameters params;
    scipy_zeros_info solver_stats;
    PyObject *f, *xargs, *item;
    volatile PyObject *fargs = NULL;
    int use_ccallback = FALSE, isllc = FALSE;
    ccallback_t callback;

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
        PyErr_SetString(PyExc_ValueError, "maxiter should be > 0");
        return NULL;
    }

    memset(&callback, sizeof(callback), '\0');
    callback.info_p = NULL;
    isllc = ccallback_is_lowlevelcallable(f);
    if (isllc) {
        use_ccallback = fill_in_ccallback(f, xargs, &callback);
        if (!use_ccallback) {
            /* error message already set inside ccallback_prepare/fill_in_ccallback */
//            PyErr_SetString(PyExc_RuntimeError, "Not a matching callback");
            return NULL;
        }
    } else if (!PyCallable_Check(f)) {
        PyErr_SetString(PyExc_RuntimeError, "Is not callable");
        return NULL;
    }


    if (use_ccallback) {
        solver_stats.error_num = 0;
        zero = solver(scipy_zeros_functions_lowlevelfunc, a, b,
                      xtol, rtol, iter, (void*)&callback, &solver_stats);
        free(callback.info_p);
        callback.info_p = NULL;
        ccallback_release(&callback);
    } else {
        len = PyTuple_Size(xargs);
        /* Make room for the double as first argument */
        fargs = PyArgs(New)(len + 1);
        if (fargs == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to allocate arguments");
            return NULL;
        }

        for (i = 0; i < len; i++) {
            item = PyTuple_GetItem(xargs, i);
            if (item == NULL) {
                Py_DECREF(fargs);
                return NULL;
            }
            Py_INCREF(item);
            PyArgs(SET_ITEM)(fargs, i+1, item);
        }

        params.function = f;
        params.args = (PyObject *)fargs;  /* Discard the volatile attribute */

        if (!setjmp(params.env)) {
            /* direct return */
            solver_stats.error_num = 0;
            zero = solver(scipy_zeros_functions_func, a, b, xtol, rtol,
                          iter, (void*)&params, &solver_stats);
            Py_DECREF(fargs);
            fargs = NULL;
        } else {
            /* error return from Python function */
            Py_DECREF(fargs);
            return NULL;
        }
    }

    if (solver_stats.error_num != 0) {
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
                flag = 1;
                return NULL;
            }
        }
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

#if PY_VERSION_HEX >= 0x03000000
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

PyObject *PyInit__zeros(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);

    return m;
}
#else
PyMODINIT_FUNC init_zeros(void)
{
        Py_InitModule("_zeros", Zerosmethods);
}
#endif
