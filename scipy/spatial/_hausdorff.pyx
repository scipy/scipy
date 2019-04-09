"""
Directed Hausdorff Code

.. versionadded:: 0.19.0

"""
#
# Copyright (C)  Tyler Reddy, Richard Gowers, and Max Linke, 2016
#
# For the supporting index shuffling algorithm:
# Copyright (C) Ben Pfaff, 2004.
#
# Distributed under the same BSD license as Scipy.
#

from __future__ import absolute_import

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt
from libc.stdlib cimport srand, rand, RAND_MAX
from libc.limits cimport UINT_MAX

__all__ = ['directed_hausdorff']

cdef extern from "numpy/npy_math.h":
    double inf "NPY_INFINITY"

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef void index_shuffler(long [::1] arr, size_t size):
    cdef:
        size_t i = 0
        size_t j
        int t

    while i < (size - 1):
        j = i + rand() / (RAND_MAX / (size - i) + 1)
        t = arr[j]
        arr[j] = arr[i]
        arr[i] = t
        i += 1

@cython.wraparound(False)
@cython.boundscheck(False)
def directed_hausdorff(const double[:,::1] ar1,
                       const double[:,::1] ar2,
                       seed=0):

    cdef double cmax, cmin, d
    cdef int N1 = ar1.shape[0]
    cdef int N2 = ar2.shape[0]
    cdef int data_dims = ar1.shape[1]
    cdef unsigned int i, j, k
    cdef unsigned int i_store = 0, j_store = 0, i_ret = 0, j_ret = 0
    cdef long[::1] resort1 = np.arange(N1)
    cdef long[::1] resort2 = np.arange(N2)

    # shuffling the points in each array generally increases the likelihood of
    # an advantageous break in the inner search loop and never decreases the
    # performance of the algorithm
    rng = np.random.RandomState(seed)
    srand(rng.randint(UINT_MAX, dtype=np.uint64))
    index_shuffler(resort1, N1)
    index_shuffler(resort2, N2)

    cmax = 0
    for i in range(N1):
        cmin = inf
        for j in range(N2):
            d = (ar1[resort1[i], 0] - ar2[resort2[j], 0])**2
	    # faster performance with square of distance
	    # avoid sqrt until very end
            for k in range(1, data_dims):
                d += (ar1[resort1[i], k] - ar2[resort2[j], k])**2
            if d < cmax: # break out of `for j` loop
                cmin = -inf
                break

            if d < cmin: # always true on first iteration of for-j loop
                cmin = d
                i_store = resort1[i]
                j_store = resort2[j]

        # always true on first iteration of for-j loop, after that only
        # if d >= cmax
        if cmin > cmax:
            cmax = cmin
            i_ret = i_store
            j_ret = j_store

    return (sqrt(cmax), i_ret, j_ret)
