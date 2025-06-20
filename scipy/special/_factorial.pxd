cdef extern from "xsf_wrappers.h" nogil:
    double xsf_gamma(double x)


cdef inline double _factorial(double n) noexcept nogil:
    if n < 0:
        return 0
    else:
        return xsf_gamma(n + 1)
