cimport numpy as np
cimport cython

ctypedef fused DTYPE_t:
    float
    float complex
    double
    double complex
    long double
    long double complex
    object


# Once Cython 3.0 is out, we can just do the following below:
#
#     with nogil(DTYPE_t is not object):
#
# But until then, we'll neeed two copies of the loops, one with
# nogil and another with gil.

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def _sosfilt(DTYPE_t [:, ::1] sos,
             DTYPE_t [:, ::1] x,
             DTYPE_t [:, :, ::1] zi):
    # Modifies x and zi in place
    cdef Py_ssize_t n_signals = x.shape[0]
    cdef Py_ssize_t n_samples = x.shape[1]
    cdef Py_ssize_t n_sections = sos.shape[0]
    cdef Py_ssize_t i, n, s
    cdef DTYPE_t x_n
    cdef DTYPE_t [:, ::1] b = sos[:, :3]
    # We ignore sos[:, 3] here because _validate_sos guarantees it's 1.
    cdef DTYPE_t [:, ::1] a = sos[:, 4:]
    if DTYPE_t is object:
        for i in range(n_signals):
            for n in range(n_samples):
                for s in range(n_sections):
                    x_n = x[i, n]  # make a temporary copy
                    # Use direct II transposed structure:
                    x[i, n] = b[s, 0] * x_n + zi[i, s, 0]
                    zi[i, s, 0] = (
                        b[s, 1] * x_n - a[s, 0] * x[i, n] + zi[i, s, 1])
                    zi[i, s, 1] = (
                        b[s, 2] * x_n - a[s, 1] * x[i, n])
    else:
        with nogil:
            for i in range(n_signals):
                for n in range(n_samples):
                    for s in range(n_sections):
                        x_n = x[i, n]  # make a temporary copy
                        # Use direct II transposed structure:
                        x[i, n] = b[s, 0] * x_n + zi[i, s, 0]
                        zi[i, s, 0] = (
                            b[s, 1] * x_n - a[s, 0] * x[i, n] + zi[i, s, 1])
                        zi[i, s, 1] = (
                            b[s, 2] * x_n - a[s, 1] * x[i, n])
