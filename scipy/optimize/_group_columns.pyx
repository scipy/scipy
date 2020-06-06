"""
Cython implementation of columns grouping for finite difference Jacobian
estimation. Used by ._numdiff.group_columns.
"""

cimport cython

import numpy as np

cimport numpy as np
from cpython cimport bool


@cython.boundscheck(False)
@cython.wraparound(False)
def group_dense(int m, int n, int [:, :] A):
    cdef int [:, :] B = A.T  # Transposed view for convenience.

    cdef int [:] groups = np.full(n, -1, dtype=np.int32)
    cdef int current_group = 0

    cdef int i, j, k

    union = np.empty(m, dtype=np.int32)
    cdef int [:] union_v = union

    # Loop through all the columns.
    for i in range(n):
        if groups[i] >= 0:  # A group was already assigned.
            continue

        groups[i] = current_group
        all_grouped = True

        union_v[:] = B[i]  # Here we store the union of grouped columns.

        for j in range(groups.shape[0]):
            if groups[j] < 0:
                all_grouped = False
            else:
                continue

            # Determine if j-th column intersects with the union.
            intersect = False
            for k in range(m):
                if union_v[k] > 0 and B[j, k] > 0:
                    intersect = True
                    break

            # If not, add it to the union and assign the group to it.
            if not intersect:
                union += B[j]
                groups[j] = current_group

        if all_grouped:
            break

        current_group += 1

    return groups.base


@cython.wraparound(False)
def group_sparse(int m, int n, int [:] indices, int [:] indptr):
    cdef int [:] groups = np.full(n, -1, dtype=np.int32)
    cdef int current_group = 0

    cdef int i, j, k

    union = np.empty(m, dtype=np.int32)
    cdef int [:] union_v = union

    for i in range(n):
        if groups[i] >= 0:
            continue

        groups[i] = current_group
        all_grouped = True

        union.fill(0)
        for k in range(indptr[i], indptr[i + 1]):
            union_v[indices[k]] = 1

        for j in range(groups.shape[0]):
            if groups[j] < 0:
                all_grouped = False
            else:
                continue

            intersect = False
            for k in range(indptr[j], indptr[j + 1]):
                if union_v[indices[k]] == 1:
                    intersect = True
                    break
            if not intersect:
                for k in range(indptr[j], indptr[j + 1]):
                    union_v[indices[k]] = 1
                groups[j] = current_group

        if all_grouped:
            break

        current_group += 1

    return groups.base
