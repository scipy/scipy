# -*- python -*- or rather like

from python_string cimport PyString_FromStringAndSize, \
    PyString_AS_STRING, PyString_Size


# Function to allocate, wrap memory via Python string creation
cdef inline object pyalloc_v(Py_ssize_t n, void **pp):
    cdef object ob = PyString_FromStringAndSize(NULL, n)
    pp[0] = <void*> PyString_AS_STRING(ob)
    return ob



