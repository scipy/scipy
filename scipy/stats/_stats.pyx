from libc cimport math
cimport cython
cimport numpy as np


ctypedef fused dtype:
    np.uint8_t
    np.uint16_t
    np.uint32_t
    np.uint64_t
    np.int8_t
    np.int16_t
    np.int32_t
    np.int64_t
    np.float32_t
    np.float64_t
    np.longdouble_t


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double ks_2samp(dtype[:] data1, dtype[:] data2):
    cdef:
        size_t i = 0, j = 0, n1 = data1.shape[0], n2 = data2.shape[0]
        dtype d1i, d2j
        double d = 0, maxd = 0, inv_n1 = 1. / n1, inv_n2 = 1. / n2
    while i < n1 and j < n2:
        d1i = data1[i]
        d2j = data2[j]
        if d1i <= d2j:
            d += inv_n1
            i += 1
        if d1i >= d2j:
            d -= inv_n2
            j += 1
        maxd = max(maxd, math.fabs(d))
    return maxd
