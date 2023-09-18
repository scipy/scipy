#ifndef PYTHONIC_TYPES_BOOL_HPP
#define PYTHONIC_TYPES_BOOL_HPP

#include "pythonic/include/types/bool.hpp"

#ifdef ENABLE_PYTHON_MODULE

PYTHONIC_NS_BEGIN
inline PyObject *to_python<bool>::convert(bool b)
{
  if (b)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

inline bool from_python<bool>::is_convertible(PyObject *obj)
{
  return obj == Py_True || obj == Py_False;
}
inline bool from_python<bool>::convert(PyObject *obj)
{
  return obj == Py_True;
}

PYTHONIC_NS_END

#endif

#endif
