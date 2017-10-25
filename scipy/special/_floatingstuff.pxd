# Utilities for writing generic functions with floating-point
# arugments.

ctypedef long double long_double

ctypedef fused floating:
    float
    double
    long_double

cdef extern from "numpy/npy_math.h":
    float NPY_INFINITYF
    double NPY_INFINITY
    long double NPY_INFINITYL
    float NPY_NANF
    double NPY_NAN
    long double NPY_NANL
    float NPY_EULERf
    double NPY_EULER
    long double NPY_EULERl

cimport libc.math
cdef extern from "math.h" nogil:
    # These aren't in Cython's `math.pxd`
    float logf(float)
    long double logl(long double)
    float fabsf(float)
    long double fabsl(long double)
    float expf(float)
    long double expl(long double)

cdef extern from "float.h":
    float FLT_EPSILON
    double DBL_EPSILON
    long double LDBL_EPSILON


cdef inline floating eps(floating x) nogil:
    if floating is float:
        return FLT_EPSILON
    elif floating is double:
        return DBL_EPSILON
    else:
        return LDBL_EPSILON


cdef inline floating infinity(floating x) nogil:
    if floating is float:
        return NPY_INFINITYF
    elif floating is double:
        return NPY_INFINITY
    else:
        return NPY_INFINITYL


cdef inline floating nan(floating x) nogil:
    if floating is float:
        return NPY_NANF
    elif floating is double:
        return NPY_NAN
    else:
        return NPY_NANL


cdef inline floating euler(floating x) nogil:
    if floating is float:
        return NPY_EULERf
    elif floating is double:
        return NPY_EULER
    else:
        return NPY_EULERl


cdef inline floating fabs(floating x) nogil:
    if floating is float:
        return fabsf(x)
    elif floating is double:
        return libc.math.fabs(x)
    else:
        return fabsl(x)


cdef inline floating log(floating x) nogil:
    if floating is float:
        return logf(x)
    elif floating is double:
        return libc.math.log(x)
    else:
        return logl(x)


cdef inline floating exp(floating x) nogil:
    if floating is float:
        return expf(x)
    elif floating is double:
        return libc.math.exp(x)
    else:
        return expl(x)
