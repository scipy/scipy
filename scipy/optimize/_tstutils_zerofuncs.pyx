"""
Test functions for optimize.zeros LowLevelCallable code
"""

from __future__ import absolute_import
cimport cython
from libc.math cimport pow

# 3 functions of increasing generality used to compute x^3 - 2,
# x**3 - 2
# x**n - 2
# x**3 - a
# x**n - a

# Function 1: Most specialized
cdef double xcubed_minus_2(double x) nogil:
    return x*x*x - 2


# Function 2: The exponent is specified at run-time
cdef double x_to_the_n_minus_2(double x, const void *user_data) nogil:
    cdef int ni = (<int *>user_data)[0]
    return pow(x, ni) - 2.0


# Now the function to compute x**n-a
cdef double x_to_the_3_minus_a(double x, const void *user_data) nogil:
    cdef double a = (<double *>user_data)[0]
    cdef double result
    result = pow(x, 3) - a
    return result


# Function 4:  Now the function to compute x**n-a
# The exponent and subtractand are specified at run-time via user_data
# First a struct to hold n, a
ctypedef struct _x_to_the_n_minus_a_data_t:
    double a
    int n

cdef double x_to_the_n_minus_a(double x, const void *user_data) nogil:
    cdef _x_to_the_n_minus_a_data_t *data = <_x_to_the_n_minus_a_data_t *>user_data
    cdef int n
    cdef double a, result
    n = data[0].n
    a = data[0].a
    result = pow(x, n) - a
    return result
