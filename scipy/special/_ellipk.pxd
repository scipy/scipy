cdef extern from "special_c_wrappers.h" nogil:
    double cephes_ellpk_wrap(double x)


cdef inline double ellipk(double m) noexcept nogil:
    return cephes_ellpk_wrap(1.0 - m)
