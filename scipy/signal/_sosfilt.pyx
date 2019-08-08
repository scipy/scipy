import numpy as np
cimport numpy as np
cimport cython

ctypedef fused DTYPE_t:
    float
    float complex
    double
    double complex
    long double
    long double complex



@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def _sosfilt(DTYPE_t [:, :] sos,
             DTYPE_t [:, :] x,
             DTYPE_t [:, :, :] zi):
    # Modifies x and zi in place
    cdef Py_ssize_t n_signals = x.shape[0]
    cdef Py_ssize_t n_samples = x.shape[1]
    cdef Py_ssize_t n_sections = sos.shape[0]
    cdef DTYPE_t x_n
    cdef DTYPE_t [:, :] b = sos[:, :3]
    cdef DTYPE_t [:, :] a = sos[:, 4:]
    with nogil:
        for i in range(n_signals):
            for n in range(n_samples):
                for s in range(n_sections):
                    x_n = x[i, n]  # make a temporary copy
                    # Use direct II transposed structure:
                    x[i, n] = b[s, 0] * x_n + zi[i, s, 0]
                    zi[i, s, 0] = b[s, 1] * x_n - a[s, 0] * x[i, n] + zi[i, s, 1]
                    zi[i, s, 1] = b[s, 2] * x_n - a[s, 1] * x[i, n]
