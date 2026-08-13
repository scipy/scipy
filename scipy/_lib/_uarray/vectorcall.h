#pragma once
#include <Python.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PYPY_VERSION
#  define Q_Py_TPFLAGS_HAVE_VECTORCALL 0
#  define Q_Py_TPFLAGS_METHOD_DESCRIPTOR 0
#  define Q_PY_VECTORCALL_ARGUMENTS_OFFSET                                     \
    ((size_t)1 << (8 * sizeof(size_t) - 1))
#else
#  define Q_Py_TPFLAGS_HAVE_VECTORCALL Py_TPFLAGS_HAVE_VECTORCALL
#  define Q_Py_TPFLAGS_METHOD_DESCRIPTOR Py_TPFLAGS_METHOD_DESCRIPTOR
#  define Q_PY_VECTORCALL_ARGUMENTS_OFFSET PY_VECTORCALL_ARGUMENTS_OFFSET
#endif

Py_ssize_t Q_PyVectorcall_NARGS(size_t n);

PyObject * Q_PyObject_Vectorcall(
    PyObject * callable, PyObject * const * args, size_t nargsf,
    PyObject * kwnames);
PyObject * Q_PyObject_VectorcallDict(
    PyObject * callable, PyObject * const * args, size_t nargsf,
    PyObject * kwdict);
PyObject * Q_PyObject_VectorcallMethod(
    PyObject * name, PyObject * const * args, size_t nargsf, PyObject * kwdict);


#ifdef __cplusplus
} // extern "C"
#endif
