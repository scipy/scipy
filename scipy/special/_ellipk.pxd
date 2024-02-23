from ._cephes cimport ellpk


cdef inline double ellipk(double m) noexcept nogil:
    return ellpk(1.0 - m)
