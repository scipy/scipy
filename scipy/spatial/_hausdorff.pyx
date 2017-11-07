"""
Directed Hausdorff Code

.. versionadded:: 0.19.0

"""
#
# Copyright (C)  Tyler Reddy, Richard Gowers, and Max Linke, 2016
#
# Distributed under the same BSD license as Scipy.
#

from __future__ import absolute_import, print_function

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt

__all__ = ['directed_hausdorff']

cdef extern from "hausdorff_util.h":
    struct return_values:
        double cmax
        int index_1
        int index_2

    void hausdorff_loop(int data_dims,
                    double ar1[],
                    double ar2[],
                    int N1,
                    int N2,
                    return_values *ret_vals)

@cython.boundscheck(False)
def directed_hausdorff(double[:,::1] ar1, double[:,::1] ar2, seed=0):

    cdef int N1 = ar1.shape[0]
    cdef int N2 = ar2.shape[0]
    cdef int data_dims = ar1.shape[1]
    cdef long[:] resort1, resort2
    cdef return_values ret_vals

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

    hausdorff_loop(data_dims, &ar1[0,0], &ar2[0,0],
                          N1, N2, &ret_vals)

    return (sqrt(ret_vals.cmax),
            resort1[ret_vals.index_1],
            resort2[ret_vals.index_2])

