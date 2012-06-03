/* MULTIPACK module by Travis Oliphant

Copyright (c) 1999 Travis Oliphant all rights reserved
oliphant.travis@ieee.org
Permission to use, modify, and distribute this software is given under the 
terms of the Scipy License

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
*/


/* This extension module is a collection of wrapper functions around
common FORTRAN code in the packages MINPACK, ODEPACK, and QUADPACK plus
some differential algebraic equation solvers.

The wrappers are meant to be nearly direct translations between the
FORTAN code and Python.  Some parameters like sizes do not need to be 
passed since they are available from the objects.  

It is anticipated that a pure Python module be written to call these lower
level routines and make a simpler user interface.  All of the routines define
default values for little-used parameters so that even the raw routines are
quite useful without a separate wrapper. 

FORTRAN Outputs that are not either an error indicator or the sought-after
results are placed in a dictionary and returned as an optional member of
the result tuple when the full_output argument is non-zero.
*/


#include "Python.h"
#include "numpy/npy_3kcompat.h"
#include "numpy/arrayobject.h"
#include <setjmp.h>


#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)


static PyObject *quadpack_python_function=NULL;
static PyObject *quadpack_extra_arguments=NULL;    /* a tuple */
static jmp_buf quadpack_jmpbuf;

static double (*quadpack_ctypes_function)(double) = NULL;

static PyObject *quadpack_error;

/* Stack Storage for re-entrant capability */
typedef struct {
    void *global0;
    void *global1;
    jmp_buf jmp;    
    PyObject *arg;
} QStorage;

typedef double (*_sp_double_func)(double);

typedef struct {
    PyObject_HEAD
    char *b_ptr;
} _sp_cfuncptr_object;

static _sp_double_func
get_ctypes_function_pointer(PyObject *obj) {
    return (*((void **)(((_sp_cfuncptr_object *)(obj))->b_ptr)));
}

static int 
quad_init_func(QStorage *store, PyObject *fun, PyObject *arg) {
    store->global0 = (void *)quadpack_python_function;
    store->global1 = (void *)quadpack_extra_arguments;
    memcpy(&(store->jmp), &quadpack_jmpbuf, sizeof(jmp_buf));
    store->arg = arg;
    if (store->arg == NULL) {
        if ((store->arg = PyTuple_New(0)) == NULL) 
            return NPY_FAIL;
    }
    else {
        Py_INCREF(store->arg);  /* We decrement on restore */
    }
    if (!PyTuple_Check(store->arg)) {
        PyErr_SetString(quadpack_error, "Extra Arguments must be in a tuple");
        Py_XDECREF(store->arg);
        return NPY_FAIL;
    }
    quadpack_python_function = fun;
    quadpack_extra_arguments = store->arg;
    return NPY_SUCCEED;
}

static void
quad_restore_func(QStorage *store, int *ierr) {
    quadpack_python_function = (PyObject *)store->global0;
    quadpack_extra_arguments = (PyObject *)store->global1;
    memcpy(&quadpack_jmpbuf, &(store->jmp), sizeof(jmp_buf));
    Py_XDECREF(store->arg);
    if (ierr != NULL) {
        if (PyErr_Occurred()) {
            *ierr = 80;             /* Python error */
            PyErr_Clear();
        }
    }
}

static int
init_ctypes_func(QStorage *store, PyObject *fun) {
    store->global0 = quadpack_ctypes_function;
    store->global1 = get_ctypes_function_pointer(fun);
    if (store->global1 == NULL) return NPY_FAIL;
    quadpack_ctypes_function = store->global1;
    return NPY_SUCCEED;
}

static void
restore_ctypes_func(QStorage *store) {
    quadpack_ctypes_function = store->global0;
}












