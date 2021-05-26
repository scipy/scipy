# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

# distutils: language = c++

import numpy as np
cimport numpy as np
from libc.math cimport fabs, sqrt, pow

np.import_array()

cdef extern from "<thread>" namespace "std" nogil:
    cdef cppclass thread:
        thread()
        void thread[A, B, C, D, E, F](A, B, C, D, E, F)
        void join()

cdef extern from "<mutex>" namespace "std" nogil:
    cdef cppclass mutex:
        void lock()
        void unlock()

cdef extern from "<functional>" namespace "std" nogil:
    cdef cppclass reference_wrapper[T]:
        pass
    cdef reference_wrapper[T] ref[T](T&)

from libcpp.vector cimport vector


cdef mutex threaded_sum_mutex

def _cy_wrapper_centered_discrepancy(double[:, ::1] sample, bint iterative,
                                     workers):
    return centered_discrepancy(sample, iterative, workers)


def _cy_wrapper_wrap_around_discrepancy(double[:, ::1] sample,
                                        bint iterative, workers):
    return wrap_around_discrepancy(sample, iterative, workers)


def _cy_wrapper_mixture_discrepancy(double[:, ::1] sample,
                                    bint iterative, workers):
    return mixture_discrepancy(sample, iterative, workers)


def _cy_wrapper_l2_star_discrepancy(double[:, ::1] sample,
                                    bint iterative, workers):
    return l2_star_discrepancy(sample, iterative, workers)


cdef double centered_discrepancy(double[:, ::1] sample_view,
                                 bint iterative, unsigned int workers) nogil:
    cdef:
        Py_ssize_t n = sample_view.shape[0]
        Py_ssize_t d = sample_view.shape[1]
        Py_ssize_t i = 0, j = 0
        double prod, disc1 = 0

    for i in range(n):
        prod = 1
        for j in range(d):
            prod *= (
                1 + 0.5 * fabs(sample_view[i, j] - 0.5) - 0.5
                * fabs(sample_view[i, j] - 0.5) ** 2
            )
        disc1 += prod

    cdef double disc2 = threaded_loops(centered_discrepancy_loop, sample_view,
                                       workers)

    if iterative:
        n += 1

    return ((13.0 / 12.0) ** d - 2.0 / n * disc1
            + 1.0 / (n ** 2) * disc2)


cdef double centered_discrepancy_loop(double[:, ::1] sample_view,
                                      Py_ssize_t istart, Py_ssize_t istop) nogil:

    cdef:
        Py_ssize_t i, j, k
        double prod, disc2 = 0

    for i in range(istart, istop):
        for j in range(sample_view.shape[0]):
            prod = 1
            for k in range(sample_view.shape[1]):
                prod *= (
                    1 + 0.5 * fabs(sample_view[i, k] - 0.5)
                    + 0.5 * fabs(sample_view[j, k] - 0.5)
                    - 0.5 * fabs(sample_view[i, k] - sample_view[j, k])
                )
            disc2 += prod

    return disc2


cdef double wrap_around_discrepancy(double[:, ::1] sample_view,
                                    bint iterative, unsigned int workers) nogil:
    cdef:
        Py_ssize_t n = sample_view.shape[0]
        Py_ssize_t d = sample_view.shape[1]
        Py_ssize_t i = 0, j = 0, k = 0
        double x_kikj, prod = 1, disc

    disc = threaded_loops(wrap_around_loop, sample_view,
                          workers)

    if iterative:
        n += 1

    return - (4.0 / 3.0) ** d + 1.0 / (n ** 2) * disc


cdef double wrap_around_loop(double[:, ::1] sample_view,
                             Py_ssize_t istart, Py_ssize_t istop) nogil:

    cdef:
        Py_ssize_t i, j, k
        double prod, disc = 0

    for i in range(istart, istop):
        for j in range(sample_view.shape[0]):
            prod = 1
            for k in range(sample_view.shape[1]):
                x_kikj = fabs(sample_view[i, k] - sample_view[j, k])
                prod *= 3.0 / 2.0 - x_kikj + x_kikj ** 2
            disc += prod

    return disc


cdef double mixture_discrepancy(double[:, ::1] sample_view,
                                bint iterative, unsigned int workers) nogil:
    cdef:
        Py_ssize_t n = sample_view.shape[0]
        Py_ssize_t d = sample_view.shape[1]
        Py_ssize_t i = 0, j = 0, k = 0
        double prod = 1, disc = 0, disc1 = 0

    for i in range(n):
        for j in range(d):
            prod *= (
                5.0 / 3.0 - 0.25 * fabs(sample_view[i, j] - 0.5)
                - 0.25 * fabs(sample_view[i, j] - 0.5) ** 2
            )
        disc1 += prod
        prod = 1

    cdef double disc2 = threaded_loops(mixture_loop, sample_view, workers)

    if iterative:
        n += 1

    disc = (19.0 / 12.0) ** d
    disc1 = 2.0 / n * disc1
    disc2 = 1.0 / (n ** 2) * disc2

    return disc - disc1 + disc2


cdef double mixture_loop(double[:, ::1] sample_view, Py_ssize_t istart,
                         Py_ssize_t istop) nogil:

    cdef:
        Py_ssize_t i, j, k
        double prod, disc2 = 0

    for i in range(istart, istop):
        for j in range(sample_view.shape[0]):
            prod = 1
            for k in range(sample_view.shape[1]):
                prod *= (15.0 / 8.0
                         - 0.25 * fabs(sample_view[i, k] - 0.5)
                         - 0.25 * fabs(sample_view[j, k] - 0.5)
                         - 3.0 / 4.0 * fabs(sample_view[i, k]
                                            - sample_view[j, k])
                         + 0.5
                         * fabs(sample_view[i, k] - sample_view[j, k]) ** 2)
            disc2 += prod

    return disc2


cdef double l2_star_discrepancy(double[:, ::1] sample_view,
                                bint iterative, unsigned int workers) nogil:
    cdef:
        Py_ssize_t n = sample_view.shape[0]
        Py_ssize_t d = sample_view.shape[1]
        Py_ssize_t i = 0, j = 0, k = 0
        double prod = 1, disc1 = 0

    for i in range(n):
        for j in range(d):
            prod *= 1 - sample_view[i, j] ** 2

        disc1 += prod
        prod = 1

    cdef double disc2 = threaded_loops(l2_star_loop, sample_view, workers)

    if iterative:
        n += 1

    cdef double one_div_n = <double> 1 / n
    return sqrt(
        pow(3, -d) - one_div_n * pow(2, 1 - d) * disc1 + 1 / pow(n, 2) * disc2
    )


cdef double l2_star_loop(double[:, ::1] sample_view, Py_ssize_t istart,
                         Py_ssize_t istop) nogil:

    cdef:
        Py_ssize_t i, j, k
        double prod = 1, disc2 = 0, tmp_sum = 0

    for i in range(istart, istop):
        for j in range(sample_view.shape[0]):
            prod = 1
            for k in range(sample_view.shape[1]):
                prod *= (
                    1 - max(sample_view[i, k], sample_view[j, k])
                )
            tmp_sum += prod

        disc2 += tmp_sum
        tmp_sum = 0

    return disc2


def _cy_wrapper_update_discrepancy(double[::1] x_new_view,
                                   double[:, ::1] sample_view,
                                   double initial_disc):
    return c_update_discrepancy(x_new_view, sample_view, initial_disc)


cdef double c_update_discrepancy(double[::1] x_new_view,
                                 double[:, ::1] sample_view,
                                 double initial_disc):
    cdef:
        Py_ssize_t n = sample_view.shape[0] + 1
        Py_ssize_t xnew_nlines = x_new_view.shape[0]
        Py_ssize_t i = 0, j = 0, k = 0
        double prod = 1, tmp_sum= 0
        double  disc1 = 0, disc2 = 0, disc3 = 0
        double[::1] abs_ = np.zeros(n, dtype=np.float64)


    # derivation from P.T. Roy (@tupui)
    for i in range(xnew_nlines):
        abs_[i] = fabs(x_new_view[i] - 0.5)
        prod *= (
            1 + 0.5 * abs_[i]
            - 0.5 * pow(abs_[i], 2)
        )

    disc1 = (- 2 / <double> n) * prod

    prod = 1
    for i in range(n - 1):
        for j in range(xnew_nlines):
            prod *= (
                1 + 0.5 * abs_[j]
                + 0.5 * fabs(sample_view[i, j] - 0.5)
                - 0.5 * fabs(x_new_view[j] - sample_view[i, j])
            )
        disc2 += prod
        prod = 1

    disc2 *= 2 / pow(n, 2)

    for i in range(xnew_nlines):
        prod *= 1 + abs_[i]

    disc3 = 1 / pow(n, 2) * prod

    return initial_disc + disc1 + disc2 + disc3


ctypedef double (*func_type)(double[:, ::1], Py_ssize_t,
                             Py_ssize_t) nogil


cdef double threaded_loops(func_type loop_func,
                           double[:, ::1] sample_view,
                           unsigned int workers) nogil:
    cdef:
        Py_ssize_t n = sample_view.shape[0]
        double disc2 = 0

    if workers <= 1:
        return loop_func(sample_view, 0, n)

    cdef:
        vector[thread] threads
        unsigned int tid
        Py_ssize_t istart, istop

    for tid in range(workers):
        istart = <Py_ssize_t> (n / workers * tid)
        istop = <Py_ssize_t> (
            n / workers * (tid + 1)) if tid < workers - 1 else n
        threads.push_back(
            thread(one_thread_loop, loop_func, ref(disc2),
                   sample_view, istart, istop)
        )

    for tid in range(workers):
        threads[tid].join()

    return disc2


cdef void one_thread_loop(func_type loop_func, double& disc, double[:,
                         ::1] sample_view, Py_ssize_t istart, Py_ssize_t istop) nogil:

    cdef double tmp = loop_func(sample_view, istart, istop)

    threaded_sum_mutex.lock()
    (&disc)[0] += tmp # workaround to "disc += tmp", see cython issue #1863
    threaded_sum_mutex.unlock()
