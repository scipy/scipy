# cython: language_level=3
"""
Cythonized routines for the RegularGridInterpolator.
"""

from libc.math cimport NAN
import numpy as np
cimport numpy as np
cimport cython

include "_poly_common.pxi"

np.import_array()


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.initializedcheck(False)
def evaluate_linear_2d(double_or_complex[:, :] values, # cannot declare as ::1
                       const long[:, :] indices,       # unless prior
                       double[:, :] norm_distances,    # np.ascontiguousarray
                       tuple grid not None,
                       double_or_complex[:] out):
    cdef:
        long num_points = indices.shape[1]      # XXX: npy_intp?
        long i0, i1, point
        double_or_complex y0, y1, result
    assert out.shape[0] == num_points

    if grid[1].shape[0] == 1:
        # linear interpolation along axis=0
        for point in range(num_points):
            i0 = indices[0, point]
            if i0 >= 0:
                y0 = norm_distances[0, point]
                result = values[i0, 0]*(1 - y0) + values[i0+1, 0]*y0
                out[point] = result
            else:
                # xi was nan: find_interval returns -1
                out[point] = NAN
    elif grid[0].shape[0] == 1:
        # linear interpolation along axis=1
        for point in range(num_points):
            i1 = indices[1, point]
            if i1 >= 0:
                y1 = norm_distances[1, point]
                result = values[0, i1]*(1 - y1) + values[0, i1+1]*y1
                out[point] = result
            else:
                # xi was nan: find_interval returns -1
                out[point] = NAN
    else:
        for point in range(num_points):
            i0, i1 = indices[0, point], indices[1, point]
            if i0 >=0 and i1 >=0:
                y0, y1 = norm_distances[0, point], norm_distances[1, point]

                result = 0.0
                result = result + values[i0, i1] * (1 - y0) * (1 - y1)
                result = result + values[i0, i1+1] * (1 - y0) * y1
                result = result + values[i0+1, i1] * y0 * (1 - y1)
                result = result + values[i0+1, i1+1] * y0 * y1
                out[point] = result
            else:
                # xi was nan
                out[point] = NAN

    return np.asarray(out)


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.initializedcheck(False)
def find_indices(tuple grid not None, const double[:, :] xi):
    # const is required for xi above in case xi is read-only
    cdef:
        long i, j, grid_i_size
        double denom, value
        # const is required in case grid is read-only
        const double[::1] grid_i

        # Axes to iterate over
        long I = xi.shape[0]
        long J = xi.shape[1]

        int index = 0

        # Indices of relevant edges between which xi are situated
        long[:,::1] indices = np.empty_like(xi, dtype=int)

        # Distances to lower edge in unity units
        double[:,::1] norm_distances = np.zeros_like(xi, dtype=float)

    # iterate through dimensions
    for i in range(I):
        grid_i = grid[i]
        grid_i_size = grid_i.shape[0]

        if grid_i_size == 1:
            # special case length-one dimensions
            for j in range(J):
                # Should equal 0. Setting it to -1 is a hack: evaluate_linear 
                # looks at indices [i, i+1] which both end up =0 with wraparound. 
                # Conclusion: change -1 to 0 here together with refactoring
                # evaluate_linear, which will also need to special-case
                # length-one axes
                indices[i, j] = -1
                # norm_distances[i, j] is already zero
        else:
            for j in range(J):
                value = xi[i, j]
                index = find_interval_ascending(&grid_i[0],
                                                grid_i_size,
                                                value,
                                                prev_interval=index,
                                                extrapolate=1)
                indices[i, j] = index

                if value == value:
                    denom = grid_i[index + 1] - grid_i[index]
                    norm_distances[i, j] = (value - grid_i[index]) / denom
                else:
                    # xi[i, j] is nan
                    norm_distances[i, j] = NAN


    return np.asarray(indices), np.asarray(norm_distances)
