#include "vectorcall.h"

#ifdef PYPY_VERSION

/* PyPy doesn't have any support for Vectorcall/FastCall.
 * These helpers are for translating to PyObject_Call. */

static PyObject * build_arg_tuple(PyObject * const * args, Py_ssize_t nargs) {
  PyObject * tuple = PyTuple_New(nargs);
  if (!tuple) {
    return NULL;
  }

  for (Py_ssize_t i = 0; i < nargs; ++i) {
    Py_INCREF(args[i]); /* SET_ITEM steals a reference */
    PyTuple_SET_ITEM(tuple, i, args[i]);
  }
  return tuple;
}

static PyObject * build_kwarg_dict(
    PyObject * const * args, PyObject * names, Py_ssize_t nargs) {
  PyObject * dict = PyDict_New();
  if (!dict) {
    return NULL;
  }

  for (Py_ssize_t i = 0; i < nargs; ++i) {
    PyObject * key = PyTuple_GET_ITEM(names, i);
    int success = PyDict_SetItem(dict, key, args[i]);
    if (success == -1) {
      Py_DECREF(dict);
      return NULL;
    }
  }
  return dict;
}
#endif /* PYPY_VERSION */


Py_ssize_t Q_PyVectorcall_NARGS(size_t n) {
  return n & (~Q_PY_VECTORCALL_ARGUMENTS_OFFSET);
}

PyObject * Q_PyObject_Vectorcall(
    PyObject * callable, PyObject * const * args, size_t nargsf,
    PyObject * kwnames) {
#ifdef PYPY_VERSION
  PyObject * dict = NULL;
  Py_ssize_t nargs = Q_PyVectorcall_NARGS(nargsf);
  if (kwnames) {
    Py_ssize_t nkwargs = PyTuple_GET_SIZE(kwnames);
    dict = build_kwarg_dict(&args[nargs - nkwargs], kwnames, nkwargs);
    if (!dict) {
      return NULL;
    }
    nargs -= nkwargs;
  }
  PyObject * ret = Q_PyObject_VectorcallDict(callable, args, nargs, dict);
  Py_XDECREF(dict);
  return ret;
#else
  return PyObject_Vectorcall(callable, args, nargsf, kwnames);
#endif
}

PyObject * Q_PyObject_VectorcallDict(
    PyObject * callable, PyObject * const * args, size_t nargsf,
    PyObject * kwdict) {
#ifdef PYPY_VERSION
  Py_ssize_t nargs = Q_PyVectorcall_NARGS(nargsf);
  PyObject * tuple = build_arg_tuple(args, nargs);
  if (!tuple) {
    return NULL;
  }
  PyObject * ret = PyObject_Call(callable, tuple, kwdict);
  Py_DECREF(tuple);
  return ret;
#else
  return PyObject_VectorcallDict(callable, args, nargsf, kwdict);
#endif
}

PyObject * Q_PyObject_VectorcallMethod(
    PyObject * name, PyObject * const * args, size_t nargsf,
    PyObject * kwnames) {
#ifdef PYPY_VERSION
  PyObject * callable = PyObject_GetAttr(args[0], name);
  if (!callable) {
    return NULL;
  }
  PyObject * result =
      Q_PyObject_Vectorcall(callable, &args[1], nargsf - 1, kwnames);
  Py_DECREF(callable);
  return result;
#else
  return PyObject_VectorcallMethod(name, args, nargsf, kwnames);
#endif
}
