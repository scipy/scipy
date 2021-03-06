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


def discrepancy(sample, bint iterative=False, method='CD', workers=1):
    """Discrepancy of a given sample.

    Parameters
    ----------
    sample : array_like (n, d)
        The sample to compute the discrepancy from.
    iterative : bool, optional
        Must be False if not using it for updating the discrepancy.
        Default is False. Refer to the notes for more details.
    method : str, optional
        Type of discrepancy, can be ``CD``, ``WD``, ``MD`` or ``L2-star``.
        Refer to the notes for more details. Default is ``CD``.
    workers : int, optional
        If `workers` is greater than 1, the population is subdivided into
        `workers` sections and evaluated in parallel (uses
        c++ threads).

    Returns
    -------
    discrepancy : float
        Discrepancy.

    Notes
    -----
    The discrepancy is a uniformity criterion used to assess the space
    filling
    of a number of samples in a hypercube. A discrepancy quantifies the
    distance between the continuous uniform distribution on a hypercube
    and the
    discrete uniform distribution on :math:`n` distinct sample points.

    The lower the value is, the better the coverage of the parameter
    space is.

    For a collection of subsets of the hypercube, the discrepancy is the
    difference between the fraction of sample points in one of those
    subsets and the volume of that subset. There are different
    definitions of
    discrepancy corresponding to different collections of subsets. Some
    versions take a root mean square difference over subsets instead of
    a maximum.

    A measure of uniformity is reasonable if it satisfies the following
    criteria [1]_:

    1. It is invariant under permuting factors and/or runs.
    2. It is invariant under rotation of the coordinates.
    3. It can measure not only uniformity of the sample over the hypercube,
       but also the projection uniformity of the sample over non-empty
       subset of lower dimension hypercubes.
    4. There is some reasonable geometric meaning.
    5. It is easy to compute.
    6. It satisfies the Koksma-Hlawka-like inequality.
    7. It is consistent with other criteria in experimental design.

    Four methods are available:

    * ``CD``: Centered Discrepancy - subspace involves a corner of the
      hypercube
    * ``WD``: Wrap-around Discrepancy - subspace can wrap around bounds
    * ``MD``: Mixture Discrepancy - mix between CD/WD covering more
    criteria
    * ``L2-star``: L2-star discrepancy - like CD BUT variant to rotation

    See [2]_ for precise definitions of each method.

    Lastly, using ``iterative=True``, it is possible to compute the
    discrepancy as if we had :math:`n+1` samples. This is useful if we want
    to add a point to a sampling and check the candidate which would
    give the
    lowest discrepancy. Then you could just update the discrepancy with
    each candidate using `update_discrepancy`. This method is faster than
    computing the discrepancy for a large number of candidates.

    References
    ----------
    .. [1] Fang et al. "Design and modeling for computer experiments".
       Computer Science and Data Analysis Series, 2006.
    .. [2] Zhou Y.-D. et al. Mixture discrepancy for quasi-random point
    sets.
       Journal of Complexity, 29 (3-4) , pp. 283-301, 2013.
    .. [3] T. T. Warnock. "Computational investigations of low discrepancy
       point sets". Applications of Number Theory to Numerical
       Analysis, Academic Press, pp. 319-343, 1972.

    Examples
    --------
    Calculate the quality of the sample using the discrepancy:

    >>> from scipy.stats import qmc
    >>> space = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
    >>> l_bounds = [0.5, 0.5]
    >>> u_bounds = [6.5, 6.5]
    >>> space = qmc.scale(space, l_bounds, u_bounds, reverse=True)
    >>> space
    array([[0.08333333, 0.41666667],
           [0.25      , 0.91666667],
           [0.41666667, 0.25      ],
           [0.58333333, 0.75      ],
           [0.75      , 0.08333333],
           [0.91666667, 0.58333333]])
    >>> qmc.discrepancy_vect(space)
    0.008142039609053464

    We can also compute iteratively the ``CD`` discrepancy by using
    ``iterative=True``.

    >>> disc_init = qmc.discrepancy_vect(space[:-1], iterative=True)
    >>> disc_init
    0.04769081147119336
    >>> qmc.update_discrepancy(space[-1], space[:-1], disc_init)
    0.008142039609053513

    """
    sample = np.asarray(sample, dtype=np.float64, order='C')

    # Checking that sample is within the hypercube and 2D
    if not sample.ndim == 2:
        raise ValueError('Sample is not a 2D array')

    if not (np.all(sample >= 0) and np.all(sample <= 1)):
        raise ValueError('Sample is not in unit hypercube')

    if not isinstance(workers, int):
        raise ValueError(f'workers param need to be an int. Found '
                         f'{type(workers)!r}')

    if method == 'CD':
        return centered_discrepancy(sample, iterative, workers=workers)

    elif method == 'WD':
        return wrap_around_discrepancy(sample, iterative, workers=workers)

    elif method == 'MD':
        return mixture_discrepancy(sample, iterative, workers=workers)

    elif method == 'L2-star':
        return l2_star_discrepancy(sample, iterative, workers=workers)
    else:
        raise ValueError('{!r} is not a valid method. Options are '
                     'CD, WD, MD, L2-star.'.format(method))


cdef double centered_discrepancy(double[:, ::1] sample_view,
                             bint iterative, unsigned int workers):
    cdef:
        Py_ssize_t n = sample_view.shape[0]
        Py_ssize_t d = sample_view.shape[1]
        Py_ssize_t i = 0, j = 0
        double prod, disc1 = 0

    for i in range(n):
        prod = 1
        for j in range(d):
                prod *= (
                        1 + 0.5 * fabs(sample_view[i, j] - 0.5) - 0.5 * fabs(
                        sample_view[i, j] - 0.5) ** 2
                )
        disc1 += prod

    cdef double disc2 = threaded_loops(centered_discrepancy_loop, sample_view,
                                       workers)

    if iterative:
        n += 1

    return ((13.0 / 12.0) ** d - 2.0 / n * disc1 +
                1.0 / (n ** 2) * disc2)


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
                             bint iterative, unsigned int workers):
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
                             bint iterative, unsigned int workers):
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
                             bint iterative, unsigned int workers):
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


def update_discrepancy(x_new, sample, double initial_disc):
    """Update the centered discrepancy with a new sample.

    Parameters
    ----------
    x_new : array_like (1, d)
        The new sample to add in `sample`.
    sample : array_like (n, d)
        The initial sample.
    initial_disc : float
        Centered discrepancy of the `sample`.

    Returns
    -------
    discrepancy : float
        Centered discrepancy of the sample composed of `x_new` and `sample`.

    Examples
    --------
    We can also compute iteratively the discrepancy by using
    ``iterative=True``.

    >>> from scipy.stats import qmc
    >>> space = np.array([[1, 3], [2, 6], [3, 2], [4, 5], [5, 1], [6, 4]])
    >>> l_bounds = [0.5, 0.5]
    >>> u_bounds = [6.5, 6.5]
    >>> space = qmc.scale(space, l_bounds, u_bounds, reverse=True)
    >>> disc_init = qmc.discrepancy(space[:-1], iterative=True)
    >>> disc_init
    0.04769081147119336
    >>> qmc.update_discrepancy(space[-1], space[:-1], disc_init)
    0.008142039609053513

    """
    sample = np.asarray(sample)
    x_new = np.asarray(x_new)

    # Checking that sample is within the hypercube and 2D
    if not sample.ndim == 2:
        raise ValueError('Sample is not a 2D array')

    if not (np.all(sample >= 0) and np.all(sample <= 1)):
        raise ValueError('Sample is not in unit hypercube')

    # Checking that x_new is within the hypercube and 1D
    if not x_new.ndim == 1:
        raise ValueError('x_new is not a 1D array')

    if not (np.all(x_new >= 0) and np.all(x_new <= 1)):
        raise ValueError('x_new is not in unit hypercube')

    if x_new.shape[0] != sample.shape[1]:
        raise ValueError('x_new and Sample must be broadcastable')

    return c_update_discrepancy(x_new, sample, initial_disc)


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
                1 + 0.5 * abs_[j] +
                0.5 * fabs(sample_view[i, j] - 0.5) -
                0.5 * fabs(x_new_view[j] - sample_view[i, j])
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


cdef double threaded_loops(func_type loop_func, double[:,
                         ::1] sample_view, unsigned int workers):
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
                   sample_view, istart, istop))

    for tid in range(workers):
        threads[tid].join()

    return disc2


cdef void one_thread_loop(func_type loop_func, double& disc, double[:,
                         ::1] sample_view, Py_ssize_t istart, Py_ssize_t istop) nogil:

    cdef double tmp = loop_func(sample_view, istart, istop)

    threaded_sum_mutex.lock()
    (&disc)[0] += tmp # workaround to "disc += tmp", see cython issue #1863
    threaded_sum_mutex.unlock()