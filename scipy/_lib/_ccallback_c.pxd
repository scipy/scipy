#
# Test function exports
#

cdef double plus1_cython(double a, int *error_flag, void *user_data) nogil except *
cdef double plus1b_cython(double a, double b, int *error_flag, void *user_data) nogil except *
cdef double plus1bc_cython(double a, double b, double c, int *error_flag, void *user_data) nogil except *

cdef double sine(double a, void *user_data) nogil except *
