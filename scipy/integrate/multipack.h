/* MULTIPACK module by Travis Oliphant

Copyright (c) 2002 Travis Oliphant all rights reserved
oliphant.travis@ieee.org
Permission to use, modify, and distribute this software is given under the 
terms of the SciPy (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

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

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

#define STORE_VARS() PyObject *store_multipack_globals[4]; int store_multipack_globals3;

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

#define INIT_JAC_FUNC(fun,Dfun,arg,col_deriv,errobj,jac_type) { \
  store_multipack_globals[0] = multipack_python_function; \
  store_multipack_globals[1] = multipack_extra_arguments; \
  store_multipack_globals[2] = multipack_python_jacobian; \
  store_multipack_globals3 = multipack_jac_transpose; \
  if (arg == NULL) { \
    if ((arg = PyTuple_New(0)) == NULL) goto fail; \
  } \
  else \
    Py_INCREF(arg);   /* We decrement on exit. */ \
  if (!PyTuple_Check(arg))  \
    PYERR(errobj,"Extra Arguments must be in a tuple"); \
  /* Set up callback functions */ \
  if (!PyCallable_Check(fun) || (Dfun != Py_None && !PyCallable_Check(Dfun))) \
    PYERR(errobj,"The function and its Jacobian must be callable functions."); \
  multipack_python_function = fun; \
  multipack_extra_arguments = arg; \
  multipack_python_jacobian = Dfun; \
  multipack_jac_transpose = !(col_deriv); \
  multipack_jac_type = jac_type;}

#define RESTORE_JAC_FUNC() multipack_python_function = store_multipack_globals[0]; \
  multipack_extra_arguments = store_multipack_globals[1]; \
  multipack_python_jacobian = store_multipack_globals[2]; \
  multipack_jac_transpose = store_multipack_globals3;

#define RESTORE_FUNC() multipack_python_function = store_multipack_globals[0]; \
  multipack_extra_arguments = store_multipack_globals[1];

#define SET_DIAG(ap_diag,o_diag,mode) { /* Set the diag vector from input */ \
  if (o_diag == NULL || o_diag == Py_None) { \
    ap_diag = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_DOUBLE); \
    if (ap_diag == NULL) goto fail; \
    diag = (double *)ap_diag -> data; \
    mode = 1; \
  } \
  else { \
    ap_diag = (PyArrayObject *)PyArray_ContiguousFromObject(o_diag, NPY_DOUBLE, 1, 1); \
    if (ap_diag == NULL) goto fail; \
    diag = (double *)ap_diag -> data; \
    mode = 2; } }

static PyObject *multipack_python_function=NULL;
static PyObject *multipack_python_jacobian=NULL;
static PyObject *multipack_extra_arguments=NULL;    /* a tuple */
static int multipack_jac_transpose=1;
static int multipack_jac_type;


static PyObject *call_python_function(PyObject *func, npy_intp n, double *x,
                                      PyObject *args, PyObject *error_obj)
{
  /*
    This is a generic function to call a python function that takes a 1-D
    sequence as a first argument and optional extra_arguments (should be a
    zero-length tuple if none desired).  The result of the function is 
    returned in a multiarray object.
        -- build sequence object from values in x.
	-- add extra arguments (if any) to an argument list.
	-- call Python callable object
 	-- check if error occurred:
	         if so return NULL
	-- if no error, place result of Python code into multiarray object.
  */

  PyArrayObject *sequence = NULL;
  PyObject *arglist = NULL;
  PyObject *arg1 = NULL, *str1 = NULL;
  PyObject *result = NULL;
  PyArrayObject *result_array = NULL;

  /* Build sequence argument from inputs */
  sequence = (PyArrayObject *)PyArray_SimpleNewFromData(1, &n, NPY_DOUBLE, (char *)x);
  if (sequence == NULL) {
    goto fail;
  }

  /* Build argument list */
  if ((arg1 = PyTuple_New(1)) == NULL) {
    Py_DECREF(sequence);
    return NULL;
  }
  PyTuple_SET_ITEM(arg1, 0, (PyObject *)sequence); 
                /* arg1 now owns sequence reference */
  if ((arglist = PySequence_Concat(arg1, args)) == NULL) {
    goto fail;
  }

  Py_DECREF(arg1);    /* arglist has a reference to sequence, now. */
  arg1 = NULL;

  /* Call function object --- variable passed to routine.  Extra
          arguments are in another passed variable.
   */
  if ((result = PyEval_CallObject(func, arglist))==NULL) {
    goto fail;
  }

  result_array = (PyArrayObject *)PyArray_ContiguousFromObject(
                                      result, NPY_DOUBLE, 0, 0);
  if (result_array == NULL) {
    goto fail;
  }

  Py_DECREF(result);
  Py_DECREF(arglist);
  return (PyObject *)result_array;

 fail:
  Py_XDECREF(arglist);
  Py_XDECREF(result);
  Py_XDECREF(arg1);
  return NULL;
}
