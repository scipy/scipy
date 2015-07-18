"""
Cython implementation of columns grouping for finite difference Jacobian
estimation. Used by ._numdiff.group_columns.
"""

import numpy as np

cimport numpy as np
from cpython cimport bool


def group_dense(int m, int n, int [:, :] A):
    cdef int [:, :] B = A.T

    groups = -np.ones(n, dtype=np.int32)
    cdef int [:] groups_v = groups
    cdef int current_group = 0

    cdef int i, j, k
    cdef bool intersect

    union = np.empty(m, dtype=np.int32)
    cdef int [:] union_v = union

    for i in range(n):
        if groups_v[i] >= 0:
            continue

        groups_v[i] = current_group
        non_grouped, = np.where(groups < 0)

        union_v[:] = B[i]
        if non_grouped.size == 0:
            break

        for j in non_grouped:
            intersect = False
            for k in range(m):
                if union_v[k] > 0 and B[j, k] > 0:
                    intersect = True
                    break

            if not intersect:
                union += B[j]
                groups_v[j] = current_group

        current_group += 1

    return groups


def group_sparse(int m, int n, int [:] indices, int [:] indptr):
    groups = -np.ones(n, dtype=np.int32)
    cdef int [:] groups_v = groups
    cdef int current_group = 0

    cdef int i, j, k
    cdef bool intersect

    union = np.empty(m, dtype=np.int32)
    cdef int [:] union_v = union

    for i in range(n):
        if groups_v[i] >= 0:
            continue

        groups_v[i] = current_group
        non_grouped, = np.where(groups < 0)

        if non_grouped.size == 0:
            break

        union.fill(0)
        for k in range(indptr[i], indptr[i + 1]):
            union_v[indices[k]] = 1

        for j in non_grouped:
            intersect = False
            for k in range(indptr[j], indptr[j + 1]):
                if union_v[indices[k]] == 1:
                    intersect = True
                    break
            if not intersect:
                for k in range(indptr[j], indptr[j + 1]):
                    union_v[indices[k]] = 1
                groups_v[j] = current_group

        current_group += 1

    return groups
