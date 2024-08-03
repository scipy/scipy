cdef extern from "xsf_wrappers.h" nogil:
    double cephes_gamma_wrap(double x)


cdef inline double _factorial(double n) noexcept nogil:
    if n < 0:
        return 0
    else:
        return cephes_gamma_wrap(n + 1)
