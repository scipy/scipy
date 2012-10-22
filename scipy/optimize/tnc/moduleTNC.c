/* Python TNC module */

/*
 * Copyright (c) 2004-2005, Jean-Sebastien Roy (js@jeannot.org)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

static char const rcsid[] =
  "@(#) $Jeannot: moduleTNC.c,v 1.12 2005/01/28 18:27:31 js Exp $";

#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tnc.h"

typedef struct _pytnc_state
{
  PyObject *py_function;
  PyObject *py_callback;
  int n;
  int failed;
} pytnc_state;

static tnc_function function;
static PyObject *moduleTNC_minimize(PyObject *self, PyObject *args);
static int PyObject_AsDouble(PyObject *py_obj, double *x);
static PyObject *PyDoubleArray_AsList(int size, double *x);
static int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size);

int PyObject_AsDouble(PyObject *py_obj, double *x)
{
  PyObject *py_float;

  py_float = PyNumber_Float(py_obj);

  if (py_float == NULL) return -1;

  *x = PyFloat_AsDouble(py_float);

  Py_DECREF(py_float);
  return 0;
}

int PyList_IntoDoubleArray(PyObject *py_list, double *x, int size)
{
  int i;

  if (py_list == NULL) return 1;

  if (!PyList_Check(py_list)) return 1;

  if (size != PyList_Size(py_list)) return 1;

  for (i=0; i<size; i++)
  {
    PyObject *py_float = PyList_GetItem(py_list, i);
    if (py_float == NULL || PyObject_AsDouble(py_float, &(x[i])))
      return 1;
  }

  return 0;
}

PyObject *PyDoubleArray_AsList(int size, double *x)
{
  int i;
  PyObject *py_list;

  py_list = PyList_New(size);
  if (py_list == NULL) return NULL;

  for (i=0; i<size; i++)
  {
    PyObject *py_float;
    py_float = PyFloat_FromDouble(x[i]);
    if (py_float == NULL || PyList_SetItem(py_list, i, py_float))
    {
      Py_DECREF(py_list);
      return NULL;
    }
  }

  return py_list;
}

static int function(double x[], double *f, double g[], void *state)
{
  PyObject *py_list, *arglist, *py_grad, *result = NULL;
  pytnc_state *py_state = (pytnc_state *)state;

  py_list = PyDoubleArray_AsList(py_state->n, x);
  if (py_list == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "tnc: memory allocation failed.");
    goto failure;
  }

  arglist = Py_BuildValue("(N)", py_list);
  result = PyEval_CallObject(py_state->py_function, arglist);
  Py_DECREF(arglist);

  if (result == NULL)
    goto failure;

  if (result == Py_None)
  {
    Py_DECREF(result);
    return 1;
  }

  if (!PyArg_ParseTuple(result, "dO!", f, &PyList_Type, &py_grad))
  {
    PyErr_SetString(PyExc_ValueError,
      "tnc: invalid return value from minimized function.");
    goto failure;
  }

  if (PyList_IntoDoubleArray(py_grad, g, py_state->n))
    goto failure;

  Py_DECREF(result);

  return 0;

failure:
  py_state->failed = 1;
  Py_XDECREF(result);
  return 1;
}

static void callback(double x[], void *state)
{
  PyObject *py_list, *arglist, *result = NULL;
  pytnc_state *py_state = (pytnc_state *)state;

  py_list = PyDoubleArray_AsList(py_state->n, x);
  if (py_list == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "tnc: memory allocation failed.");
  }

  arglist = Py_BuildValue("(N)", py_list);
  result = PyEval_CallObject(py_state->py_callback, arglist);
  Py_DECREF(arglist);
  Py_XDECREF(result);
}

PyObject *moduleTNC_minimize(PyObject *self, PyObject *args)
{
  PyArrayObject *py_x0, *py_low, *py_up, *py_scale, *py_offset;
  PyObject *py_function = NULL;
  PyObject *py_callback = NULL;
  pytnc_state py_state;
  int n, n1, n2, n3, n4;
  tnc_callback *callback_function = NULL;

  int rc, msg, maxCGit, maxnfeval, nfeval = 0, niter = 0;
  double *x = NULL, *low = NULL, *up = NULL, *scale = NULL, *offset = NULL;
  double f, eta, stepmx, accuracy, fmin, ftol, xtol, pgtol, rescale;

  if (!PyArg_ParseTuple(args, "OO!O!O!O!O!iiiddddddddO",
    &py_function,
    &PyArray_Type, &py_x0,
    &PyArray_Type, &py_low,
    &PyArray_Type, &py_up,
    &PyArray_Type, &py_scale,
    &PyArray_Type, &py_offset,
    &msg, &maxCGit, &maxnfeval, &eta, &stepmx, &accuracy, &fmin, &ftol,
    &xtol, &pgtol,
    &rescale,
    &py_callback
    ))
    return NULL;

  if (!PyCallable_Check(py_function))
  {
    PyErr_SetString(PyExc_TypeError, "tnc: function must be callable");
    return NULL;
  }

  if ((n3 = PyArray_Size((PyObject *)py_scale)) != 0)
  {
    scale = (double *)PyArray_GETPTR1(py_scale, 0);
    if (scale == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "tnc: invalid scaling parameters.");
      return NULL;
    }
  }

  if ((n4 = PyArray_Size((PyObject *)py_offset)) != 0)
  {
    offset = (double *)PyArray_GETPTR1(py_offset, 0);
    if (offset == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "tnc: invalid offset parameters.");
      return NULL;
    }
  }

  if ((n = PyArray_Size((PyObject *)py_x0)) != 0)
  {
    x = (double *)PyArray_GETPTR1(py_x0, 0);
    if (x == NULL)
    {
      PyErr_SetString(PyExc_ValueError, "tnc: invalid initial vector.");
      return NULL;
    }
  }

  if ((n1 = PyArray_Size((PyObject *)py_low)) != 0)
  {
    low = (double *)PyArray_GETPTR1(py_low, 0);
  }
  if ((n2 = PyArray_Size((PyObject *)py_up)) != 0)
  {
    up = (double *)PyArray_GETPTR1(py_up, 0);
  }

  if ((n1 != 0 && low == NULL) || (n2 != 0 && up == NULL))
  {
    PyErr_SetString(PyExc_ValueError, "tnc: invalid bounds.");
    return NULL;
  }

  if (n1 != n2 || n != n1 || (scale != NULL && n != n3)
    || (offset != NULL && n != n4))
  {
    PyErr_SetString(PyExc_ValueError, "tnc: vector sizes must be equal.");
    return NULL;
  }

  py_state.py_function = py_function;
  py_state.n = n;
  py_state.failed = 0;

  Py_INCREF(py_function);

  if (py_callback != Py_None)
  {
      if (!PyCallable_Check(py_callback))
      {
        PyErr_SetString(PyExc_TypeError,
                        "tnc: callback must be callable or None.");
        return NULL;
      }
      py_state.py_callback = py_callback;
      Py_INCREF(py_callback);
      callback_function = callback;
  }

  rc = tnc(n, x, &f, NULL, function, &py_state, low, up, scale, offset, msg,
    maxCGit, maxnfeval, eta, stepmx, accuracy, fmin, ftol, xtol, pgtol, rescale,
    &nfeval, &niter, callback_function);

  Py_DECREF(py_function);
  if (py_callback != Py_None) Py_DECREF(py_callback);

  if (py_state.failed) return NULL;

  if (rc == TNC_ENOMEM)
  {
    PyErr_SetString(PyExc_MemoryError, "tnc: memory allocation failed.");
    return NULL;
  }

  return Py_BuildValue("(iiiO)", rc, nfeval, niter, PyArray_Return(py_x0));
}

static PyMethodDef moduleTNC_methods[] =
{
  {"minimize", moduleTNC_minimize, METH_VARARGS},
  {NULL, NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "moduleTNC",
    NULL,
    -1,
    moduleTNC_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit_moduleTNC(void)
{
    return PyModule_Create(&moduledef);
}
#else
PyMODINIT_FUNC initmoduleTNC(void)
{
  (void) Py_InitModule("moduleTNC", moduleTNC_methods);
  import_array();
}
#endif
