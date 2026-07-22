"""
Spherical Voronoi Cython Code

.. versionadded:: 0.19.0

"""
#
# Copyright (C)  Tyler Reddy, 2016
#
# Distributed under the same BSD license as Scipy.
#

import numpy as np
cimport numpy as np
cimport cython

np.import_array()

__all__ = ['sort_vertices_of_regions']


@cython.boundscheck(False)
def sort_vertices_of_regions(const int[:,::1] simplices, list regions):
    cdef np.npy_intp n, k, s, i, max_len
    cdef np.npy_intp num_regions = len(regions)
    cdef np.npy_intp current_simplex, current_vertex
    cdef np.npy_intp remaining_size
    cdef np.npy_intp[:] remaining
    cdef np.ndarray[np.intp_t, ndim=1] sorted_vertices

    max_len = 0
    for region in regions:
        max_len = max(max_len, len(region))
    sorted_vertices = np.empty(max_len, dtype=np.intp)

    for n in range(num_regions):
        remaining = np.asarray(regions[n][:])
        remaining_size = remaining.shape[0]
        current_simplex = remaining[0]
        for i in range(3):
            k = simplices[current_simplex, i]
            if k != n:
                current_vertex = k
                break
        for k in range(remaining_size):
            sorted_vertices[k] = current_simplex
            for i in range(remaining_size):
                if remaining[i] == sorted_vertices[k]:
                    continue
                s = remaining[i]
                for j in range(3):
                    if current_vertex == simplices[s, j]:
                        current_simplex = s
            for i in range(3):
                s = simplices[current_simplex, i]
                if s != n and s != current_vertex:
                    current_vertex = s
                    break
        regions[n] = list(sorted_vertices[:remaining_size])
