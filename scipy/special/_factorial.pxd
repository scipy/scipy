from ._cephes cimport Gamma


cdef extern from "numpy/npy_math.h" nogil:
    double npy_isnan(double)


cdef inline double _factorial(double n) nogil:
    if n < 0:
        return 0
    elif npy_isnan(n):
        return n
    elif <int>n != n:
        with gil:
            raise ValueError("factorial() only accepts integral values")
    else:
        return Gamma(n + 1)
