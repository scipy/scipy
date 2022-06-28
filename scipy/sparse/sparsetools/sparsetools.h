#ifndef SPARSETOOLS_H
#define SPARSETOOLS_H

#include <Python.h>
#include "numpy/ndarrayobject.h"

#include <stdexcept>

#include "bool_ops.h"
#include "complex_ops.h"

typedef PY_LONG_LONG thunk_t(int I_typenum, int T_typenum, void **args);

NPY_VISIBILITY_HIDDEN PyObject *
call_thunk(char ret_spec, const char *spec, thunk_t *thunk, PyObject *args);

#endif
