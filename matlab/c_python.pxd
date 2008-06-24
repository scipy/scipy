# -*- python -*-
# :Author:    Robert Kern 
# :Copyright: 2004, Enthought, Inc.
# :License:   BSD Style


cdef extern from "Python.h":
    # Not part of the Python API, but we might as well define it here.
    # Note that the exact type doesn't actually matter for Pyrex.
    ctypedef int size_t

    # Some type declarations we need
    ctypedef int Py_intptr_t


    # String API
    char* PyString_AsString(object string)
    char* PyString_AS_STRING(object string)
    object PyString_FromString(char* c_string)
    object PyString_FromStringAndSize(char* c_string, int length)
    object PyString_InternFromString(char *v)

    # Float API
    object PyFloat_FromDouble(double v)
    double PyFloat_AsDouble(object ob)
    long PyInt_AsLong(object ob)


    # Memory API
    void* PyMem_Malloc(size_t n)
    void* PyMem_Realloc(void* buf, size_t n)
    void PyMem_Free(void* buf)

    void Py_DECREF(object obj)
    void Py_XDECREF(object obj)
    void Py_INCREF(object obj)
    void Py_XINCREF(object obj)

    # CObject API
    ctypedef void (*destructor1)(void* cobj)
    ctypedef void (*destructor2)(void* cobj, void* desc)
    int PyCObject_Check(object p)
    object PyCObject_FromVoidPtr(void* cobj, destructor1 destr)
    object PyCObject_FromVoidPtrAndDesc(void* cobj, void* desc, 
        destructor2 destr)
    void* PyCObject_AsVoidPtr(object self)
    void* PyCObject_GetDesc(object self)
    int PyCObject_SetVoidPtr(object self, void* cobj)  

    # TypeCheck API
    int PyFloat_Check(object obj)
    int PyInt_Check(object obj)

    # Error API
    int PyErr_Occurred()
    void PyErr_Clear()
    int  PyErr_CheckSignals()

    # File API
    ctypedef struct FILE
    FILE* PyFile_AsFile(object)

cdef extern from "stdio.h":
    size_t fread(void *ptr, size_t size, size_t n, FILE *file)

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)

cdef extern from "math.h":
    double fabs(double x)

cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass

