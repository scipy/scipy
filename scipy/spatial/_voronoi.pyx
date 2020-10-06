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

__all__ = ['sort_vertices_of_regions']

# array-filling placeholder that can never occur
DEF ARRAY_FILLER = -2

@cython.boundscheck(False)
cdef void remaining_filter(np.npy_intp[:] remaining,
                           np.npy_intp current_simplex):
    cdef np.npy_intp i
    for i in range(remaining.shape[0]):
        if remaining[i] == current_simplex:
            remaining[i] = ARRAY_FILLER


@cython.boundscheck(False)
def sort_vertices_of_regions(int[:,::1] simplices, regions):
    cdef np.npy_intp n, k, s, i
    cdef np.npy_intp num_regions = len(regions)
    cdef np.npy_intp current_simplex, current_vertex, remaining_count
    cdef np.npy_intp cs_identified, remaining_size
    cdef np.npy_intp[:] remaining
    cdef np.npy_intp[:] sorted_vertices = np.empty(max([len(region) for region
                                                   in regions]),
                                                   dtype=np.intp)

    for n in range(num_regions):
        remaining_count = 0
        remaining = np.asarray(regions[n][:])
        remaining_size = remaining.shape[0]
        sorted_vertices[:] = ARRAY_FILLER
        current_simplex = remaining[0]
        for i in range(3):
            k = simplices[current_simplex, i]
            if k != n:
                current_vertex = k
                break
        sorted_vertices[remaining_count] = current_simplex
        remaining_count += 1
        remaining_filter(remaining, current_simplex)
        while remaining_size > remaining_count:
            cs_identified = 0
            for i in range(remaining_size):
                if remaining[i] == ARRAY_FILLER:
                    continue
                s = remaining[i]
                for j in range(3):
                    if current_vertex == simplices[s, j]:
                        current_simplex = remaining[i]
                        cs_identified += 1
                        break
                if cs_identified > 0:
                    break
            for i in range(3):
                s = simplices[current_simplex, i]
                if s != n and s != current_vertex:
                    current_vertex = s
                    break
            sorted_vertices[remaining_count] = current_simplex
            remaining_count += 1
            remaining_filter(remaining, current_simplex)
        regions_arr = np.asarray(sorted_vertices)
        regions[n] = regions_arr[regions_arr > ARRAY_FILLER].tolist()
