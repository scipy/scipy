# cython: show_performance_hints=False

from cpython.pycapsule cimport (
    PyCapsule_CheckExact, PyCapsule_New, PyCapsule_SetContext, PyCapsule_GetName, PyCapsule_GetPointer,
    PyCapsule_GetContext
)
from cpython.long cimport PyLong_AsVoidPtr
from libc.stdlib cimport free
from libc.string cimport strdup
from libc.math cimport sin

from .ccallback cimport (ccallback_t, ccallback_prepare, ccallback_release, CCALLBACK_DEFAULTS,
                         ccallback_signature_t)


#
# PyCapsule helpers
#

cdef void raw_capsule_destructor(object capsule) noexcept:
    cdef const char *name
    name = PyCapsule_GetName(capsule)
    free(<char*>name)


def get_raw_capsule(func_obj, name_obj, context_obj):
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
        const char *capsule_name
        const char *name
        const char *name_copy

    if name_obj is None:
        name = NULL
    elif not isinstance(name_obj, bytes):
        name_obj = name_obj.encode('ascii')
        name = <char*>name_obj
    else:
        name = <char*>name_obj

    if PyCapsule_CheckExact(context_obj):
        capsule_name = PyCapsule_GetName(context_obj)
        context = PyCapsule_GetPointer(context_obj, capsule_name)
    elif context_obj is None:
        context = NULL
    else:
        context = PyLong_AsVoidPtr(int(context_obj))

    if PyCapsule_CheckExact(func_obj):
        capsule_name = PyCapsule_GetName(func_obj)
        func = PyCapsule_GetPointer(func_obj, capsule_name)

        if context == NULL:
            context = PyCapsule_GetContext(func_obj)

        if name == NULL:
            name = capsule_name
    else:
        func = PyLong_AsVoidPtr(int(func_obj))

    if name == NULL:
        name_copy = name
    else:
        name_copy = strdup(name)

    capsule = PyCapsule_New(func, name_copy, &raw_capsule_destructor)
    if context != NULL:
        PyCapsule_SetContext(capsule, context)
    return capsule


def get_capsule_signature(capsule_obj):
    cdef const char *name
    name = PyCapsule_GetName(capsule_obj)
    if name == NULL:
        raise ValueError("Capsule has no signature")
    return bytes(name).decode('ascii')


def check_capsule(item):
    """
    check_capsule(item)

    Return True if the given object is a PyCapsule.
    """
    if PyCapsule_CheckExact(item):
        return True
    return False


sigs = [
    (b"double (double, int *, void *)", 0),
    (b"double (double, double, int *, void *)", 1)
]

if sizeof(int) == sizeof(long):
    sigs.append((b"double (double, long *, void *)", 0))
    sigs.append((b"double (double, double, long *, void *)", 1))

if sizeof(int) == sizeof(short):
    sigs.append((b"double (double, short *, void *)", 0))
    sigs.append((b"double (double, double, short *, void *)", 1))

cdef ccallback_signature_t signatures[7]

for idx, sig in enumerate(sigs):
    signatures[idx].signature = sig[0]
    signatures[idx].value = sig[1]

signatures[idx + 1].signature = NULL


cdef double test_thunk_cython(double a, int *error_flag, void *data) except? -1.0 nogil:
    """
    Implementation of a thunk routine in Cython
    """
    cdef:
        ccallback_t *callback = <ccallback_t *>data
        double result = 0

    if callback.c_function != NULL:
        if callback.signature.value == 0:
            result = (<double(*)(double, int *, void *) nogil>callback.c_function)(
                a, error_flag, callback.user_data)
        else:
            result = (<double(*)(double, double, int *, void *) nogil>callback.c_function)(
                a, 0.0, error_flag, callback.user_data)

        if error_flag[0]:
            # Python error has been set by the callback
            return -1.0
    else:
        with gil:
            try:
                return float((<object>callback.py_function)(a))
            except:  # noqa: E722
                error_flag[0] = 1
                raise

    return result


def test_call_cython(callback_obj, double value):
    """
    Implementation of a caller routine in Cython
    """
    cdef:
        ccallback_t callback
        int error_flag = 0
        double result

    ccallback_prepare(&callback, signatures, callback_obj, CCALLBACK_DEFAULTS)

    with nogil:
        result = test_thunk_cython(value, &error_flag, <void *>&callback)

    ccallback_release(&callback)

    return result


cdef double plus1_cython(double a, int *error_flag, void *user_data) except * nogil:
    """
    Implementation of a callable in Cython
    """
    # 2.0 is ERROR_VALUE in the test code (src/_test_callback.c)
    if a == 2.0:
        error_flag[0] = 1
        with gil:
            raise ValueError("failure...")

    if user_data == NULL:
        return a + 1
    else:
        return a + (<double *>user_data)[0]

cdef double plus1b_cython(double a, double b, int *error_flag, void *user_data) except * nogil:
    return plus1_cython(a, error_flag, user_data) + b

cdef double plus1bc_cython(double a, double b, double c, int *error_flag, void *user_data) except * nogil:
    return plus1_cython(a, error_flag, user_data) + b + c

cdef double sine(double x, void *user_data) except * nogil:
    return sin(x)


# Ctypes declarations of the callables above
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
