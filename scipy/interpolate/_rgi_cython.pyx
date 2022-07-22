import numpy as np
cimport numpy as np
cimport cython

np.import_array()

import itertools

ctypedef fused double_or_long:
    double
    long


cdef extern from "numpy/npy_math.h":
    double nan "NPY_NAN"


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef int find_interval_ascending(const double *x,
                                 size_t nx,
                                 double xval,
                                 int prev_interval=0,
                                 bint extrapolate=1) nogil:
    """
    Find an interval such that x[interval] <= xval < x[interval+1]. Assuming
    that x is sorted in the ascending order.
    If xval < x[0], then interval = 0, if xval > x[-1] then interval = n - 2.

    Parameters
    ----------
    x : array of double, shape (m,)
        Piecewise polynomial breakpoints sorted in ascending order.
    xval : double
        Point to find.
    prev_interval : int, optional
        Interval where a previous point was found.
    extrapolate : bint, optional
        Whether to return the last of the first interval if the
        point is out-of-bounds.

    Returns
    -------
    interval : int
        Suitable interval or -1 if nan.

    """
    cdef int interval, high, low, mid
    cdef double a, b

    a = x[0]
    b = x[nx-1]

    interval = prev_interval
    if interval < 0 or interval >= nx:
        interval = 0

    if not (a <= xval <= b):
        # Out-of-bounds (or nan)
        if xval < a and extrapolate:
            # below
            interval = 0
        elif xval > b and extrapolate:
            # above
            interval = nx - 2
        else:
            # nan or no extrapolation
            interval = -1
    elif xval == b:
        # Make the interval closed from the right
        interval = nx - 2
    else:
        # Find the interval the coordinate is in
        # (binary search with locality)
        if xval >= x[interval]:
            low = interval
            high = nx - 2
        else:
            low = 0
            high = interval

        if xval < x[low+1]:
            high = low

        while low < high:
            mid = (high + low)//2
            if xval < x[mid]:
                # mid < high
                high = mid
            elif xval >= x[mid + 1]:
                low = mid + 1
            else:
                # x[mid] <= xval < x[mid+1]
                low = mid
                break

        interval = low

    return interval


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
@cython.boundscheck(True)
@cython.cdivision(True)
def find_indices(tuple grid not None, double[:, :] xi):
    cdef long i, j, I, J, grid_i_size
    cdef double denom, value
    cdef long[:,::1] indices
    cdef double[:,::1] norm_distances
    cdef double[::1] grid_i

    cdef int index = 0
    # find relevant edges between which xi are situated
    indices = np.empty_like(xi, dtype=int)
    # compute distance to lower edge in unity units
    norm_distances = np.zeros_like(xi, dtype=float)

    # iterate through dimensions
    for i in range(I):
        grid_i = grid[i]
        grid_i_size = grid_i.shape[0]

        if grid_i_size == 1:
            # special case length-one dimensions
            for j in range(J):
                indices[i, j] = -1    # Should equal 0. Setting it to -1 is a hack: evaluate_linear looks at indices [i, i+1]
                                      # which both end up =0 with wraparound. 
                                      # Conclusion: change -1 to 0 here together with refactoring evaluate_linear, which
                                      # will also need to special-case length-one axes
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
                    norm_distances[i, j] = nan


    return np.asarray(indices), np.asarray(norm_distances)
