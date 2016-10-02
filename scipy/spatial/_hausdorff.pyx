"""
Directed Hausdorff Code

.. versionadded:: 0.19.0

"""
#
# Copyright (C)  Tyler Reddy and Richard Gowers, 2016
#
# Distributed under the same BSD license as Scipy.
#

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt

__all__ = ['directed_hausdorff']

@cython.boundscheck(False)
def directed_hausdorff(np.ndarray[np.float64_t, ndim =2] ar1,
                       np.ndarray[np.float64_t, ndim =2] ar2):

    cdef double cmax, cmin
    cdef int break_occurred
    cdef int N1 = ar1.shape[0]
    cdef int N2 = ar2.shape[0]
    cdef int data_dims = ar1.shape[1]
    cdef np.float64_t d
    cdef unsigned int i, j, k

    # shuffling the points in each array generally increases the likelihood of
    # an advantageous break in the inner search loop and never decreases the
    # performance of the algorithm
    np.random.shuffle(ar1)
    np.random.shuffle(ar2)
                                                                                                                                                                                                     
    cmax = 0 
    for i in range(N1):
        break_occurred = 0
        cmin = np.inf
        for j in range(N2):
            d = 0
	    # faster performance with square of distance
	    # avoid sqrt until very end
            for k in range(data_dims):
                d += (ar1[i, k] - ar2[j, k]) * (ar1[i, k] - ar2[j, k])
            if d < cmax: # early break
                break_occurred += 1
                break
            if d < cmin:
                cmin = d
        if cmin > cmax and cmin != np.inf and break_occurred == 0:
            cmax = cmin
    return sqrt(cmax)
