from ._cephes cimport Gamma


cdef inline double _factorial(double n) noexcept nogil:
    if n < 0:
        return 0
    else:
        return Gamma(n + 1)
