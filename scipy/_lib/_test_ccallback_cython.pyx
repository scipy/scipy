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

cdef double plus1b_cython(double a, double b, int *error_flag, void *user_data) except *:
    return plus1_cython(a, error_flag, user_data) + b

cdef double plus1bc_cython(double a, double b, double c, int *error_flag, void *user_data) except *:
    return plus1_cython(a, error_flag, user_data) + b + c


import ctypes

plus1_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
plus1_ctypes = ctypes.cast(<size_t>&plus1_cython, plus1_t)

plus1b_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
plus1b_ctypes = ctypes.cast(<size_t>&plus1b_cython, plus1b_t)

plus1bc_t = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                            ctypes.POINTER(ctypes.c_int), ctypes.c_void_p)
plus1bc_ctypes = ctypes.cast(<size_t>&plus1bc_cython, plus1bc_t)

# Note that the PyCapsule approach in src/_test_ccallback.c:get_plus1_capsule can
# also be used to prepare Cython callback functions.
