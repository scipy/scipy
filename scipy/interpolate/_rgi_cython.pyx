import numpy as np
cimport numpy as np
cimport cython

np.import_array()

import itertools


def evaluate_linear(values, indices, norm_distances, out_of_bounds):
    # slice for broadcasting over trailing dimensions in self.values
    vslice = (slice(None),) + (None,) * (values.ndim - len(indices))

    # Compute shifting up front before zipping everything together
    shift_norm_distances = [1 - yi for yi in norm_distances]
    shift_indices = [i + 1 for i in indices]

    # The formula for linear interpolation in 2d takes the form:
    # values = self.values[(i0, i1)] * (1 - y0) * (1 - y1) + \
    #          self.values[(i0, i1 + 1)] * (1 - y0) * y1 + \
    #          self.values[(i0 + 1, i1)] * y0 * (1 - y1) + \
    #          self.values[(i0 + 1, i1 + 1)] * y0 * y1
    # We pair i with 1 - yi (zipped1) and i + 1 with yi (zipped2)
    zipped1 = zip(indices, shift_norm_distances)
    zipped2 = zip(shift_indices, norm_distances)

    # Take all products of zipped1 and zipped2 and iterate over them
    # to get the terms in the above formula. This corresponds to iterating
    # over the vertices of a hypercube.
    hypercube = itertools.product(*zip(zipped1, zipped2))
    value = np.array([0.])
    for h in hypercube:
        edge_indices, weights = zip(*h)
        weight = np.array([1.])
        for w in weights:
            weight = weight * w
        value = value + np.asarray(values[edge_indices]) * weight[vslice]
    return value


@cython.wraparound(False)
@cython.boundscheck(False)
def find_indices(grid, bounds_error, xi):
    cdef int i
    cdef int j
    cdef np.ndarray[long, ndim=2] indices
    cdef np.ndarray[double, ndim=2] norm_distances
    cdef np.ndarray denom

    # find relevant edges between which xi are situated
    indices = np.empty_like(xi, dtype=int)
    # compute distance to lower edge in unity units
    norm_distances = np.empty_like(xi, dtype=np.float64)
    # check for out of bounds xi
    out_of_bounds = np.zeros(xi.shape[1], dtype=bool)

    # iterate through dimensions
    for i in range(xi.shape[0]):
        index = np.searchsorted(grid[i], xi[i]) - 1
        index[index < 0] = 0
        index[index > grid[i].size - 2] = grid[i].size - 2
        indices[i] = index

        # compute norm_distances, incl length-1 grids,
        # where `grid[index+1] == grid[index]`
        denom = grid[i][index + 1] - grid[i][index]
        for j in range(denom.shape[0]):
            if denom[j] == 0:
                norm_distances[i][j] = 0
            else:
                norm_distances[i][j] = (xi[i][j] - grid[i][index][j]) / denom[j]

        if not bounds_error:
            out_of_bounds += xi[i] < grid[i][0]
            out_of_bounds += xi[i] > grid[i][grid[i].shape[0] - 1]
    return indices, norm_distances, out_of_bounds
