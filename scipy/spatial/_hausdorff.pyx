"""
Directed Hausdorff Code

.. versionadded:: 0.19.0

"""
#
# Copyright (C)  Tyler Reddy, Richard Gowers, and Max Linke, 2016
#
# Distributed under the same BSD license as Scipy.
#

from __future__ import absolute_import

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt

__all__ = ['directed_hausdorff']

@cython.boundscheck(False)
def directed_hausdorff(double[:,::1] ar1, double[:,::1] ar2, seed=0):

    cdef double cmax, cmin, d
    cdef bint no_break_occurred
    cdef int N1 = ar1.shape[0]
    cdef int N2 = ar2.shape[0]
    cdef int data_dims = ar1.shape[1]
    cdef unsigned int i, j, k
    cdef unsigned int i_store = 0, j_store = 0, i_ret = 0, j_ret = 0
    cdef long[:] resort1, resort2

    # shuffling the points in each array generally increases the likelihood of
    # an advantageous break in the inner search loop and never decreases the
    # performance of the algorithm
    rng = np.random.RandomState(seed)
    resort1 = np.arange(N1)
    resort2 = np.arange(N2)
    rng.shuffle(resort1)
    rng.shuffle(resort2)
    ar1 = np.asarray(ar1)[resort1]
    ar2 = np.asarray(ar2)[resort2]

    cmax = 0
    for i in range(N1):
        no_break_occurred = True
        cmin = np.inf
        for j in range(N2):
            d = 0
	    # faster performance with square of distance
	    # avoid sqrt until very end
            for k in range(data_dims):
                d += (ar1[i, k] - ar2[j, k])**2
            if d < cmax: # break out of `for j` loop
                no_break_occurred = False
                break

            if d < cmin: # always true on first iteration of for-j loop
                cmin = d
                i_store = i
                j_store = j

        # always true on first iteration of for-j loop, after that only
        # if d >= cmax
        if cmin != np.inf and cmin > cmax and no_break_occurred == True:
            cmax = cmin
            i_ret = i_store
            j_ret = j_store

    return (sqrt(cmax), resort1[i_ret], resort2[j_ret])
