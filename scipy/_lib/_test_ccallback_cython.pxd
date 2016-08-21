#
# Test code for src/ccallback.h
#

cdef double plus1_cython(double a, int *error_flag, void *user_data) except *
