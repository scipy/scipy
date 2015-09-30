from scipy.linalg.cython_lapack cimport dlartg
from scipy.linalg.cython_blas cimport drot

import numpy as np


def givens_elimination(double [:, ::1] S, double [:] v, double [:] diag):
    """Zero out a diagonal block of a matrix by series of Givens rotations.

    The matrix has the structure::

        [ S ]
        [ D ]

    Where S is an upper triangular matrix with shape (n, n) and D is a
    diagonal matrix with shape (n, n) with elements from `diag`. This function
    applies Givens rotations to it such that the resulting matrix has zeros
    in place of D.

    Array `S` will be modified in-place.

    Array `v` of shape (n,) is the part of the full vector with shape (2*n,)::

        [ v ]
        [ 0 ]

    to which Givens rotations are applied. This array is modified in place,
    such that on exit it contains the first n components of the above
    mentioned vector after rotations were applied.
    """
    cdef int n = diag.size
    cdef int k

    cdef int i, j

    cdef double f, g, r
    cdef double cs, sn
    cdef int one = 1

    cdef double [:] diag_row = np.empty(n)
    cdef double u  # For `v` rotations.

    for i in range(n):
        if diag[i] == 0:
            continue

        diag_row[:] = 0
        diag_row[i] = diag[i]
        u = 0

        for j in range(i, n):
            if diag_row[j] != 0:
                f = S[j, j]
                g = diag_row[j]

                # Compute cosine and sine of rotation angle.
                dlartg(&f, &g, &cs, &sn, &r)
                S[j, j] = r
                # diag_row[j] is implicitly 0 now.

                # Now rotate the remaining elements in rows.
                k = n - j - 1
                if k > 0:
                    drot(&k, &S[j, j+1], &one, &diag_row[j+1], &one, &cs, &sn)

                # Some custom code for rotating `v`.
                f = v[j]
                v[j] = cs * f + sn * u
                u = -sn * f + cs * u
