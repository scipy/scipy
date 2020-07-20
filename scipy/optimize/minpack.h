/* MULTIPACK module by Travis Oliphant

Copyright (c) 2002 Travis Oliphant all rights reserved
Oliphant.Travis@altavista.net
Permission to use, modify, and distribute this software is given under the 
terms of the SciPy (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
*/


/* This extension module is a collection of wrapper functions around
common FORTRAN code in the packages MINPACK, ODEPACK, and QUADPACK plus
some differential algebraic equation solvers.

The wrappers are meant to be nearly direct translations between the
FORTRAN code and Python.  Some parameters like sizes do not need to be 
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
#include "numpy/arrayobject.h"
#include "ccallback.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}

#define STORE_VARS() ccallback_t callback; int callback_inited = 0; jac_callback_info_t jac_callback_info;

#define INIT_FUNC(fun,arg,errobj) do { /* Get extra arguments or set to zero length tuple */ \
  if (arg == NULL) { \
    if ((arg = PyTuple_New(0)) == NULL) goto fail_free; \
  } \
  else \
    Py_INCREF(arg);   /* We decrement on exit. */ \
  if (!PyTuple_Check(arg))  \
    PYERR(errobj,"Extra Arguments must be in a tuple"); \
  /* Set up callback functions */ \
  if (!PyCallable_Check(fun)) \
    PYERR(errobj,"First argument must be a callable function."); \
  if (init_callback(&callback, fun, arg) != 0) \
    PYERR(errobj,"Could not init callback");\
  callback_inited = 1; \
  } while(0)

#define INIT_JAC_FUNC(fun,Dfun,arg,col_deriv,errobj) do { \
  if (arg == NULL) { \
    if ((arg = PyTuple_New(0)) == NULL) goto fail_free; \
  } \
  else \
    Py_INCREF(arg);   /* We decrement on exit. */ \
  if (!PyTuple_Check(arg))  \
    PYERR(errobj,"Extra Arguments must be in a tuple"); \
  /* Set up callback functions */ \
  if (!PyCallable_Check(fun) || (Dfun != Py_None && !PyCallable_Check(Dfun))) \
    PYERR(errobj,"The function and its Jacobian must be callable functions."); \
  if (init_jac_callback(&callback, &jac_callback_info, fun, Dfun, arg, col_deriv) != 0) \
    PYERR(errobj,"Could not init callback");\
  callback_inited = 1; \
} while(0)

#define RESTORE_JAC_FUNC() do { \
  if (callback_inited && release_callback(&callback) != 0) { \
    goto fail_free; \
  }} while(0)

#define RESTORE_FUNC() do { \
  if (callback_inited && release_callback(&callback) != 0) { \
    goto fail_free; \
  }} while(0)

#define SET_DIAG(ap_diag,o_diag,mode) { /* Set the diag vector from input */ \
  if (o_diag == NULL || o_diag == Py_None) { \
    ap_diag = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_DOUBLE); \
    if (ap_diag == NULL) goto fail; \
    diag = (double *)PyArray_DATA(ap_diag); \
    mode = 1; \
  } \
  else { \
    ap_diag = (PyArrayObject *)PyArray_ContiguousFromObject(o_diag, NPY_DOUBLE, 1, 1); \
    if (ap_diag == NULL) goto fail; \
    diag = (double *)PyArray_DATA(ap_diag); \
    mode = 2; } }

#define MATRIXC2F(jac,data,m,n) {double *p1=(double *)(jac), *p2, *p3=(double *)(data);\
int i,j;\
for (j=0;j<(m);p3++,j++) \
  for (p2=p3,i=0;i<(n);p2+=(m),i++,p1++) \
    *p1 = *p2; }

typedef struct {
  PyObject *Dfun;
  PyObject *extra_args;
  int jac_transpose;
} jac_callback_info_t;

static PyObject *call_python_function(PyObject *func, npy_intp n, double *x, PyObject *args, int dim, PyObject *error_obj, npy_intp out_size)
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
  PyObject *arg1 = NULL;
  PyObject *result = NULL;
  PyArrayObject *result_array = NULL;
  npy_intp fvec_sz = 0;

  /* Build sequence argument from inputs */
  sequence = (PyArrayObject *)PyArray_SimpleNewFromData(1, &n, NPY_DOUBLE, (char *)x);
  if (sequence == NULL) PYERR2(error_obj,"Internal failure to make an array of doubles out of first\n                 argument to function call.");

  /* Build argument list */
  if ((arg1 = PyTuple_New(1)) == NULL) {
    Py_DECREF(sequence);
    return NULL;
  }
  PyTuple_SET_ITEM(arg1, 0, (PyObject *)sequence); 
                /* arg1 now owns sequence reference */
  if ((arglist = PySequence_Concat( arg1, args)) == NULL)
    PYERR2(error_obj,"Internal error constructing argument list.");

  Py_DECREF(arg1);    /* arglist has a reference to sequence, now. */
  arg1 = NULL;

  /* Call function object --- variable passed to routine.  Extra
          arguments are in another passed variable.
   */
  if ((result = PyEval_CallObject(func, arglist))==NULL) {
      goto fail;
  }

  if ((result_array = (PyArrayObject *)PyArray_ContiguousFromObject(result, NPY_DOUBLE, dim-1, dim))==NULL) 
    PYERR2(error_obj,"Result from function call is not a proper array of floats.");

  fvec_sz = PyArray_SIZE(result_array);
  if(out_size != -1 && fvec_sz != out_size){
      PyErr_SetString(PyExc_ValueError, "The array returned by a function changed size between calls");
      Py_DECREF(result_array);
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
