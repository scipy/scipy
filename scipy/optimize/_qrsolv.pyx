cdef extern from "alloca.h":
    void *alloca(int) nogil

cdef extern from "string.h":
    void *memset(void *, int, int) nogil

cdef extern from "math.h":
    double sqrt(double) nogil
    double fabs(double) nogil

cimport numpy
cimport cython
import numpy

@cython.boundscheck(False)
@cython.cdivision(True)
cdef void qrsolv_float(numpy.ndarray[float, ndim=2] s,
           numpy.ndarray[float] diag) nogil:
    cdef unsigned int N
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l
    cdef float sin
    cdef float cos
    cdef float tan
    cdef float cotan
    cdef float tmp
    cdef float *ta

    N = diag.shape[0]
    ta = <float *> alloca((N + 1) * sizeof(float))

    for 0 <= j < N:
        if diag[j] != 0:
            memset(ta, 0, (N + 1) * sizeof(float))
            ta[j] = diag[j]
            for j <= k < N:
                if ta[k] != 0:
                    if fabs(s[k, k]) > fabs(ta[k]):
                        tan = ta[k] / s[k, k]
                        cos = 1 / sqrt(1 + tan * tan)
                        sin = cos * tan
                    else:
                        cotan = s[k, k] / ta[k]
                        sin = 1 / sqrt(1 + cotan * cotan)
                        cos = sin * cotan
                    for k <= l <= N:
                        tmp = s[k, l]
                        s[k, l] = cos * tmp + sin * ta[l]
                        ta[l] = -sin * tmp + cos * ta[l]

@cython.boundscheck(False)
@cython.cdivision(True)
cdef void qrsolv_double(numpy.ndarray[double, ndim=2] s,
           numpy.ndarray[double] diag) nogil:
    cdef unsigned int N
    cdef unsigned int j
    cdef unsigned int k
    cdef unsigned int l
    cdef double sin
    cdef double cos
    cdef double tan
    cdef double cotan
    cdef double tmp
    cdef double *ta

    N = diag.shape[0]
    ta = <double *> alloca((N + 1) * sizeof(double))

    for 0 <= j < N:
        if diag[j] != 0:
            memset(ta, 0, (N + 1) * sizeof(double))
            ta[j] = diag[j]
            for j <= k < N:
                if ta[k] != 0:
                    if fabs(s[k, k]) > fabs(ta[k]):
                        tan = ta[k] / s[k, k]
                        cos = 1 / sqrt(1 + tan * tan)
                        sin = cos * tan
                    else:
                        cotan = s[k, k] / ta[k]
                        sin = 1 / sqrt(1 + cotan * cotan)
                        cos = sin * cotan
                    for k <= l <= N:
                        tmp = s[k, l]
                        s[k, l] = cos * tmp + sin * ta[l]
                        ta[l] = -sin * tmp + cos * ta[l]

def qrsolv(numpy.ndarray s, numpy.ndarray diag):
    if s.dtype == numpy.float64:
        qrsolv_double(s, diag)
    elif s.dtype == numpy.float32:
        qrsolv_float(s, diag)
    else:
        raise NotImplementedError(
            "qrsolv is only implemented for float23 and float64")
    return s
