/* MULTIPACK module by Travis Oliphant

Copyright (c) 1999 Travis Oliphant all rights reserved
Oliphant.Travis@altavista.net
Permission to use, modify, and distribute this software is given under the 
terms of the LGPL.  See http://www.fsf.org

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
#include "Numeric/arrayobject.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

#define MAX(n1,n2) ((n1) > (n2))?(n1):(n2);
#define MIN(n1,n2) ((n1) > (n2))?(n2):(n1);

#define STORE_VARS() PyObject *store_multipack_globals[4]

#define INIT_FUNC(fun,arg,errobj) { /* Get extra arguments or set to zero length tuple */ \
  store_multipack_globals[0] = multipack_python_function; \
  store_multipack_globals[1] = multipack_extra_arguments; \
  if (arg == NULL) { \
    if ((arg = PyTuple_New(0)) == NULL) goto fail; \
  } \
  else \
    Py_INCREF(arg);   /* We decrement on exit. */ \
  if (!PyTuple_Check(arg))  \
    PYERR(errobj,"Extra Arguments must be in a tuple"); \
  /* Set up callback functions */ \
  if (!PyCallable_Check(fun)) \
    PYERR(errobj,"First argument must be a callable function."); \
  multipack_python_function = fun; \
  multipack_extra_arguments = arg; }

#define RESTORE_FUNC() multipack_python_function = store_multipack_globals[0]; \
  multipack_extra_arguments = store_multipack_globals[1];

static PyObject *multipack_python_function=NULL;
static PyObject *multipack_extra_arguments=NULL;    /* a tuple */













