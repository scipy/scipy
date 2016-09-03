from cpython.pycapsule cimport (
    PyCapsule_CheckExact, PyCapsule_New, PyCapsule_SetContext, PyCapsule_GetName, PyCapsule_GetPointer,
    PyCapsule_GetContext
)
from cpython.long cimport PyLong_AsVoidPtr
from libc.stdlib cimport malloc, free
from libc.string cimport strdup


#
# PyCapsule helpers
#

cdef void raw_capsule_destructor(object capsule):
    cdef char *name
    name = PyCapsule_GetName(capsule)
    free(name)

def get_raw_capsule(func_obj, char *name, context_obj):
    """
    get_raw_capsule(ptr, name, context)

    Create a new PyCapsule with given pointer, name, and context.

    Parameters
    ----------
    ptr : {PyCapsule, int}
        Memory address of the pointer.
    name : str
        Python string containing the signature.
    context : {PyCapsule, int}
        Memory address of the context.
        If NULL and ptr is a PyCapsule, use the one from the context of ptr.

    """
    cdef:
        void *func
        void *context
        char *capsule_name
        char *name_copy

    if PyCapsule_CheckExact(context_obj):
        capsule_name = PyCapsule_GetName(context_obj)
        context = PyCapsule_GetPointer(context_obj, capsule_name)
    elif context_obj is None:
        context = NULL
    else:
        context = PyLong_AsVoidPtr(long(context_obj))

    if PyCapsule_CheckExact(func_obj):
        capsule_name = PyCapsule_GetName(func_obj)
        func = PyCapsule_GetPointer(func_obj, capsule_name)

        if context == NULL:
            context = PyCapsule_GetContext(func_obj)

        if name[0] == '\0':
            name = capsule_name
    else:
        func = PyLong_AsVoidPtr(long(func_obj))

    if name == NULL:
        name_copy = name
    else:
        name_copy = strdup(name)

    capsule = PyCapsule_New(func, name_copy, &raw_capsule_destructor)
    if context != NULL:
        PyCapsule_SetContext(capsule, context)
    return capsule


def check_capsule(item):
    """
    check_capsule(item)

    Return True if the given object is a PyCapsule.
    """
    if PyCapsule_CheckExact(item):
        return True
    return False


#
# Test code for src/ccallback.h
#

DEF ERROR_VALUE = 2

from libc.math cimport sin

cdef double plus1_cython(double a, int *error_flag, void *user_data) nogil except *:
    if a == ERROR_VALUE:
        error_flag[0] = 1
        with gil:
            raise ValueError("failure...")

    if user_data == NULL:
        return a + 1
    else:
        return a + (<double *>user_data)[0]

cdef double plus1b_cython(double a, double b, int *error_flag, void *user_data) nogil except *:
    return plus1_cython(a, error_flag, user_data) + b

cdef double plus1bc_cython(double a, double b, double c, int *error_flag, void *user_data) nogil except *:
    return plus1_cython(a, error_flag, user_data) + b + c

cdef double sine(double x, void *user_data) nogil except *:
    return sin(x)


import ctypes

plus1_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
plus1_ctypes = ctypes.cast(<size_t>&plus1_cython, plus1_t)

plus1b_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
plus1b_ctypes = ctypes.cast(<size_t>&plus1b_cython, plus1b_t)

plus1bc_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
plus1bc_ctypes = ctypes.cast(<size_t>&plus1bc_cython, plus1bc_t)

sine_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_void_p)
sine_ctypes = ctypes.cast(<size_t>&sine, sine_t)
