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
#elif (PY_VERSION_HEX < 0x03090000)
// clang-format off

static int is_method_descr(PyTypeObject* descr_tp) {
    return (
        (descr_tp->tp_flags & Q_Py_TPFLAGS_METHOD_DESCRIPTOR) != 0 ||
        (descr_tp == &PyFunction_Type) ||
        (descr_tp == &PyMethodDescr_Type));
}

/* The below code is derivative of CPython source, and is taken because
_PyObject_GetMethod is not exported from shared objects.

Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010,
2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 Python Software Foundation;
All Rights Reserved
*/

/* Specialized version of _PyObject_GenericGetAttrWithDict
   specifically for the LOAD_METHOD opcode.

   Return 1 if a method is found, 0 if it's a regular attribute
   from __dict__ or something returned by using a descriptor
   protocol.

   `method` will point to the resolved attribute or NULL.  In the
   latter case, an error will be set.
*/
static int _PyObject_GetMethod(PyObject *obj, PyObject *name, PyObject **method)
{
    PyTypeObject *tp = Py_TYPE(obj);
    PyObject *descr;
    descrgetfunc f = NULL;
    PyObject **dictptr, *dict;
    PyObject *attr;
    int meth_found = 0;

    assert(*method == NULL);

    if (Py_TYPE(obj)->tp_getattro != PyObject_GenericGetAttr
            || !PyUnicode_Check(name)) {
        *method = PyObject_GetAttr(obj, name);
        return 0;
    }

    if (tp->tp_dict == NULL && PyType_Ready(tp) < 0)
        return 0;

    descr = _PyType_Lookup(tp, name);
    if (descr != NULL) {
        Py_INCREF(descr);
        if (is_method_descr(Py_TYPE(descr))) {
            meth_found = 1;
        } else {
            f = Py_TYPE(descr)->tp_descr_get;
            if (f != NULL && PyDescr_IsData(descr)) {
                *method = f(descr, obj, (PyObject *)Py_TYPE(obj));
                Py_DECREF(descr);
                return 0;
            }
        }
    }

    dictptr = _PyObject_GetDictPtr(obj);
    if (dictptr != NULL && (dict = *dictptr) != NULL) {
        Py_INCREF(dict);
        attr = PyDict_GetItemWithError(dict, name);
        if (attr != NULL) {
            Py_INCREF(attr);
            *method = attr;
            Py_DECREF(dict);
            Py_XDECREF(descr);
            return 0;
        }
        else {
            Py_DECREF(dict);
            if (PyErr_Occurred()) {
                Py_XDECREF(descr);
                return 0;
            }
        }
    }

    if (meth_found) {
        *method = descr;
        return 1;
    }

    if (f != NULL) {
        *method = f(descr, obj, (PyObject *)Py_TYPE(obj));
        Py_DECREF(descr);
        return 0;
    }

    if (descr != NULL) {
        *method = descr;
        return 0;
    }

    PyErr_Format(PyExc_AttributeError,
                 "'%.50s' object has no attribute '%U'",
                 tp->tp_name, name);
    return 0;
}

// clang-format on
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
#elif (PY_VERSION_HEX >= 0x03090000)
  return PyObject_Vectorcall(callable, args, nargsf, kwnames);
#elif (PY_VERSION_HEX >= 0x03080000)
  return _PyObject_Vectorcall(callable, args, nargsf, kwnames);
#else
  Py_ssize_t nargs = Q_PyVectorcall_NARGS(nargsf);
  return _PyObject_FastCallKeywords(
      callable, (PyObject **)args, nargs, kwnames);
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
#elif (PY_VERSION_HEX >= 0x03090000)
  return PyObject_VectorcallDict(callable, args, nargsf, kwdict);
#else
  Py_ssize_t nargs = Q_PyVectorcall_NARGS(nargsf);
  return _PyObject_FastCallDict(callable, (PyObject **)args, nargs, kwdict);
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
#elif (PY_VERSION_HEX >= 0x03090000)
  return PyObject_VectorcallMethod(name, args, nargsf, kwnames);
#else
  /* Private CPython code for CALL_METHOD opcode */
  PyObject * callable = NULL;
  int unbound = _PyObject_GetMethod(args[0], name, &callable);
  if (callable == NULL) {
    return NULL;
  }

  if (unbound) {
    /* We must remove PY_VECTORCALL_ARGUMENTS_OFFSET since
     * that would be interpreted as allowing to change args[-1] */
    nargsf &= ~Q_PY_VECTORCALL_ARGUMENTS_OFFSET;
  } else {
    /* Skip "self". We can keep PY_VECTORCALL_ARGUMENTS_OFFSET since
     * args[-1] in the onward call is args[0] here. */
    args++;
    nargsf--;
  }
  PyObject * result = Q_PyObject_Vectorcall(callable, args, nargsf, kwnames);
  Py_DECREF(callable);
  return result;
#endif
}
