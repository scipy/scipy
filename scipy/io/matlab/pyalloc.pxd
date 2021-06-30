# -*- python -*- or rather like

from cpython cimport PyBytes_FromStringAndSize, \
    PyBytes_AS_STRING, PyBytes_Size


# Function to allocate, wrap memory via Python string creation
cdef inline object pyalloc_v(Py_ssize_t n, void **pp):
    cdef object ob = PyBytes_FromStringAndSize(NULL, n)
    pp[0] = <void*> PyBytes_AS_STRING(ob)
    return ob
