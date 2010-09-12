#ifndef __NPY_PY3H__
#define __NPY_PY3H__

#include <Python.h>

#if PY_VERSION_HEX >= 0x03000000
    #define PyUstring_FromString PyUnicode_FromString
    #define PyInt_FromLong PyLong_FromLong
    #define PyInt_Check PyLong_Check
#else
    #define PyUstring_FromString PyString_FromString
#endif

#endif
