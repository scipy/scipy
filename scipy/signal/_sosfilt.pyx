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

    # jumping through a few memoryview hoops to reduce array lookups
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


def _sosfilt(DTYPE_floating_t [:, ::1] sos,
             DTYPE_floating_t [:, ::1] x,
             DTYPE_floating_t [:, :, ::1] zi):
        with nogil:
            _sosfilt_float(sos, x, zi)
