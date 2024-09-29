cimport numpy as np
cimport cython

np.import_array()

ctypedef fused DTYPE_floating_t:
    float
    float complex
    double
    double complex
    long double
    long double complex

ctypedef fused DTYPE_t:
    DTYPE_floating_t
    object


# Once Cython 3.0 is out, we can just do the following below:
#
#     with nogil(DTYPE_t is not object):
#
# But until then, we'll need two copies of the loops, one with
# nogil and another with gil.

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _sosfilt_float(DTYPE_floating_t [:, ::1] sos,
                         DTYPE_floating_t [:, ::1] x,
                         DTYPE_floating_t [:, :, ::1] zi) noexcept nogil:
    # Modifies x and zi in place
    cdef Py_ssize_t n_signals = x.shape[0]
    cdef Py_ssize_t n_samples = x.shape[1]
    cdef Py_ssize_t n_sections = sos.shape[0]
    cdef Py_ssize_t i, n, s
    cdef DTYPE_floating_t x_new, x_cur
    cdef DTYPE_floating_t[:, ::1] zi_slice
    cdef DTYPE_floating_t const_1 = 1.0

    # jumping through a few memoryview hoops to reduce array lookups,
    # the original version is still in the gil version below.
    for i in xrange(n_signals):
        zi_slice = zi[i, :, :]
        for n in xrange(n_samples):

            x_cur = const_1 * x[i, n]  # make sure x_cur is a copy

            for s in xrange(n_sections):
                x_new = sos[s, 0] * x_cur + zi_slice[s, 0]
                zi_slice[s, 0] = (sos[s, 1] * x_cur - sos[s, 4] * x_new
                                  + zi_slice[s, 1])
                zi_slice[s, 1] = sos[s, 2] * x_cur - sos[s, 5] * x_new
                x_cur = x_new

            x[i, n] = x_cur


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def _sosfilt_object(object [:, ::1] sos,
                    object [:, ::1] x,
                    object [:, :, ::1] zi):
    # Modifies x and zi in place
    cdef Py_ssize_t n_signals = x.shape[0]
    cdef Py_ssize_t n_samples = x.shape[1]
    cdef Py_ssize_t n_sections = sos.shape[0]
    cdef Py_ssize_t i, n, s
    cdef object x_n
    for i in xrange(n_signals):
        for n in xrange(n_samples):
            for s in xrange(n_sections):
                x_n = x[i, n]  # make a temporary copy
                # Use direct II transposed structure:
                x[i, n] = sos[s, 0] * x_n + zi[i, s, 0]
                zi[i, s, 0] = (sos[s, 1] * x_n - sos[s, 4] * x[i, n] +
                               zi[i, s, 1])
                zi[i, s, 1] = (sos[s, 2] * x_n - sos[s, 5] * x[i, n])


def _sosfilt(DTYPE_t [:, ::1] sos,
             DTYPE_t [:, ::1] x,
             DTYPE_t [:, :, ::1] zi):
    if DTYPE_t is object:
        _sosfilt_object(sos, x, zi)
    else:
        with nogil:
            _sosfilt_float(sos, x, zi)
