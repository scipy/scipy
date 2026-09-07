"""
Pythran implementation of columns grouping for finite difference Jacobian
estimation. Used by ._numdiff.group_columns and based on the Cython version.
"""

import numpy as np

#pythran export group_dense(int, int, intc[:,:])
#pythran export group_dense(int, int, int[:,:])
def group_dense(m, n, A):
    B = A.T  # Transposed view for convenience.

    groups = np.full(n, -1, dtype=np.intp)
    current_group = 0

    union = np.empty(m, dtype=A.dtype)

    # Loop through all the columns.
    for i in range(n):
        if groups[i] >= 0:  # A group was already assigned.
            continue

        groups[i] = current_group
        all_grouped = True

        union[:] = B[i]  # Here we store the union of grouped columns.

        for j in range(i + 1, groups.shape[0]):
            if groups[j] < 0:
                all_grouped = False
            else:
                continue

            # Determine if j-th column intersects with the union.
            for k in range(m):
                if union[k] > 0 and B[j, k] > 0:
                    break
            # If not, add it to the union and assign the group to it.
            else:
                union += B[j]
                groups[j] = current_group

        if all_grouped:
            break

        current_group += 1

    return groups


#pythran export group_sparse(int, int, int32[], int32[])
#pythran export group_sparse(int, int, int64[], int64[])
#pythran export group_sparse(int, int, int32[::], int32[::])
#pythran export group_sparse(int, int, int64[::], int64[::])
def group_sparse(m, n, indices, indptr):
    groups = np.full(n, -1, dtype=np.intp)
    current_group = 0

    union = np.empty(m, dtype=bool)

    for i in range(n):
        if groups[i] >= 0:
            continue

        groups[i] = current_group
        all_grouped = True

        union.fill(False)
        for k in range(indptr[i], indptr[i + 1]):
            union[indices[k]] = True

        for j in range(i + 1, groups.shape[0]):
            if groups[j] < 0:
                all_grouped = False
            else:
                continue

            for k in range(indptr[j], indptr[j + 1]):
                if union[indices[k]]:
                    break
            else:
                for k in range(indptr[j], indptr[j + 1]):
                    union[indices[k]] = True
                groups[j] = current_group

        if all_grouped:
            break

        current_group += 1

    return groups
