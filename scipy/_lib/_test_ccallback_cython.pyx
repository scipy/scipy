#
# Test code for src/ccallback.h
#

DEF ERROR_VALUE = 2

cdef double plus1_cython(double a, int *error_flag, void *user_data) except *:
    if a == ERROR_VALUE:
        error_flag[0] = 1
        raise ValueError("failure...")

    if user_data == NULL:
        return a + 1
    else:
        return a + <object>user_data

import ctypes
plus1_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
plus1_ctypes = ctypes.cast(<size_t>&plus1_cython, plus1_t)

# Note that the PyCapsule approach in src/_test_ccallback.c:get_plus1_capsule can
# also be used to prepare Cython callback functions.
