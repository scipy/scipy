#!/usr/bin/env python3
# =============================================================================
#     File: check_normcond_behavior.py
#  Created: 2025-07-02 12:20
#   Author: Bernie Roesler
#
"""Development script to test the behavior of `np.linalg.norm`,
`np.linalg.cond`, `scipy.sparse.linalg.norm`, and
`scipy.sparse.linalg.onenormest` on various sparse matrices, so that we can
make `sicpy.sparse.linalg.cond1est` behave consistently with these existing
functions.

This script is not part of the test suite. Run manually.
"""
# =============================================================================

import numpy as np
import numpy.linalg as nla
import scipy.linalg as sla

from scipy import sparse
from scipy.sparse import linalg as spla

from numpy.testing import assert_warns, assert_allclose, assert_raises_regex


# -----------------------------------------------------------------------------
#         Empty Matrix
# -----------------------------------------------------------------------------
A = sparse.csc_array((0, 0))

# ----- norm(A)
assert(nla.norm(A.toarray(), ord=1) == 0)

with assert_raises_regex(ValueError, "Improper number of dimensions to norm"):
    sla.norm(A, ord=1)

# raises ValueError: abs(A).sum(axis=0).max() fails because abs(A).sum(axis=0)
# produces an empty matrix
with assert_raises_regex(
    ValueError,
    "zero-size array to reduction operation maximum which has no identity"
):
    spla.norm(A, ord=1)

with assert_raises_regex(ValueError, "attempt to get argmax of an empty sequence"):
    spla.onenormest(A)

# ----- norm(inv(A))
assert(nla.inv(A.toarray()).shape == (0, 0))
assert(sla.inv(A.toarray()).shape == (0, 0))

assert(nla.norm(nla.inv(A.toarray()), ord=1) == 0)
assert(sla.norm(sla.inv(A.toarray()), ord=1) == 0)

with assert_raises_regex(ValueError, "need at least one array to concatenate"):
    # spla.inv() uses SuperLU, which raises this error
    spla.inv(A)

# ----- cond(A)
with assert_raises_regex(nla.LinAlgError, "cond is not defined on empty arrays"):
    nla.cond(A.toarray(), p=1)


# -----------------------------------------------------------------------------
#         All Zero Matrix
# -----------------------------------------------------------------------------
A = sparse.csc_array((3, 3))

# ----- norm(A)
assert(nla.norm(A.toarray(), ord=1) == 0)
assert(sla.norm(A.toarray(), ord=1) == 0)
assert(spla.norm(A, ord=1) == 0)
assert(spla.onenormest(A) == 0)

# ----- norm(inv(A))
with assert_raises_regex(nla.LinAlgError, "Singular matrix"):
    nla.inv(A.toarray())

with assert_raises_regex(sla.LinAlgError, "A singular matrix detected"):
    sla.inv(A.toarray())

with assert_raises_regex(RuntimeError, "Factor is exactly singular"):
    # spla.inv() uses SuperLU, which raises this error
    spla.inv(A)

# ----- cond(A)
assert(nla.cond(A.toarray(), p=1) == np.inf)


# -----------------------------------------------------------------------------
#         Singleton Matrix
# -----------------------------------------------------------------------------
A = sparse.csc_array([[2.0]])

# ----- norm(A)
assert_allclose(nla.norm(A.toarray(), ord=1), 2.0)
assert_allclose(sla.norm(A.toarray(), ord=1), 2.0)
assert_allclose(spla.norm(A, ord=1), 2.0)
assert_allclose(spla.onenormest(A), 2.0)

# ----- norm(inv(A))
assert_allclose(nla.norm(nla.inv(A.toarray()), ord=1), 0.5)
assert_allclose(sla.norm(sla.inv(A.toarray()), ord=1), 0.5)

# The inverse is a *dense* array!
sp_A_inv = spla.inv(A)
assert_allclose(sp_A_inv, np.array([0.5]), strict=True)

with assert_raises_regex(TypeError, "input is not sparse"):
    spla.norm(spla.inv(A), ord=1)

assert_allclose(spla.onenormest(spla.inv(A)), 0.5)


# -----------------------------------------------------------------------------
#         Exactly Singular Matrix
# -----------------------------------------------------------------------------
A = sparse.random_array((5, 5), density=0.5, format="lil", random_state=5656)
A[0] = 0       # make it exactly singular
A = A.tocsc()  # avoid SparseEfficiencyWarning

# ----- norm(A)
the_norm = nla.norm(A.toarray(), ord=1)
assert_allclose(sla.norm(A.toarray(), ord=1), the_norm)
assert_allclose(spla.norm(A, ord=1), the_norm)
assert_allclose(spla.onenormest(A), the_norm)

# ----- norm(inv(A))
with assert_raises_regex(nla.LinAlgError, "Singular matrix"):
    nla.inv(A.toarray())

with assert_raises_regex(sla.LinAlgError, "A singular matrix detected"):
    sla.inv(A.toarray())

with assert_raises_regex(RuntimeError, "Factor is exactly singular"):
    # spla.inv() uses SuperLU, which raises this error
    spla.inv(A)

# ----- cond(A)
assert(nla.cond(A.toarray(), p=1) == np.inf)


# -----------------------------------------------------------------------------
#         Nearly Singular Matrix
# -----------------------------------------------------------------------------
A = sparse.random_array((5, 5), density=0.5, format="lil", random_state=5656)
eps = np.finfo(float).eps
A[0] = eps * eps  # make it nearly singular
A = A.tocsc()     # avoid SparseEfficiencyWarning

# ----- norm(A)
the_norm = nla.norm(A.toarray(), ord=1)
assert_allclose(sla.norm(A.toarray(), ord=1), the_norm)
assert_allclose(spla.norm(A, ord=1), the_norm)
assert_allclose(spla.onenormest(A), the_norm)

# ----- norm(inv(A))
the_norm = nla.norm(nla.inv(A.toarray()), ord=1)

with assert_warns(sla.LinAlgWarning):
    # "An ill-conditioned matrix is detected"
    assert_allclose(sla.norm(sla.inv(A.toarray()), ord=1), the_norm)

assert_allclose(spla.norm(spla.inv(A), ord=1), the_norm)
assert_allclose(spla.onenormest(spla.inv(A)), the_norm)


# -----------------------------------------------------------------------------
#         1-D array
# -----------------------------------------------------------------------------
A = sparse.coo_array(np.arange(5))  # CSC only for 2D

# ----- norm(A)
the_norm = nla.norm(A.toarray(), ord=1)
assert_allclose(sla.norm(A.toarray(), ord=1), the_norm)
assert_allclose(spla.norm(A, ord=1), the_norm)

with assert_raises_regex(ValueError, r"invalid shape.*(must be 2-d)"):
    spla.onenormest(A)

# ----- norm(inv(A))
with assert_raises_regex(nla.LinAlgError, "Array must be at least two-dimensional"):
    nla.inv(A.toarray())

with assert_raises_regex(ValueError, "Expected at least ndim=2"):
    sla.inv(A.toarray())

with assert_raises_regex(IndexError, "tuple index out of range"):
    spla.inv(A)

# ----- cond(A)
with assert_raises_regex(nla.LinAlgError, "Array must be at least two-dimensional"):
    nla.cond(A.toarray(), p=1)


# -----------------------------------------------------------------------------
#         Non-Square Array (Non-Singular)
# -----------------------------------------------------------------------------
A = sparse.random_array((3, 5), density=0.5, random_state=5656)
A.setdiag(1.0)  # make it non-singular
A = A.tocsc()   # avoid SparseEfficiencyWarning

# ----- norm(A)
the_norm = nla.norm(A.toarray(), ord=1)
assert_allclose(sla.norm(A.toarray(), ord=1), the_norm)
assert_allclose(spla.norm(A, ord=1), the_norm)

with assert_raises_regex(
    ValueError,
    "expected the operator to act like a square matrix"
):
    spla.onenormest(A)

# ----- norm(inv(A))
with assert_raises_regex(
    nla.LinAlgError,
    "Last 2 dimensions of the array must be square"
):
    nla.inv(A.toarray())

with assert_raises_regex(ValueError, "Expected square matrix"):
    sla.inv(A.toarray())

with assert_raises_regex(ValueError, "matrix must be square"):
    spla.inv(A)

# ----- cond(A)
with assert_raises_regex(
    nla.LinAlgError,
    "Last 2 dimensions of the array must be square"
):
    nla.cond(A.toarray(), p=1)


# -----------------------------------------------------------------------------
#         Multi-Dimensional Array
# -----------------------------------------------------------------------------
# Define a dense array, since COO is not subscriptable
M, N = 3, 5
A = sparse.random_array((M, N, N), density=0.5, random_state=5656).toarray()
# make each slice non-singular
for i in range(M):
    for j in range(N):
        A[i, j, j] = 1.0

A = sparse.coo_array(A)  # convert to sparse format for testing

# ----- norm(A)
with assert_raises_regex(ValueError, "Improper number of dimensions to norm"):
    nla.norm(A.toarray(), ord=1)

# Pass additional axis argument to specify which 2 axes are square
axis = (1, 2)
the_norm = nla.norm(A.toarray(), ord=1, axis=axis)
assert(the_norm.shape == (M,))

assert_allclose(sla.norm(A.toarray(), ord=1, axis=axis), the_norm)

# NOTE likely a bug in spla.norm, because the docs say it should accept the
# axis argument:
#     axis : {int, 2-tuple of ints, None}, optional
#         If `axis` is an integer, it specifies the axis of `x` along which to
#         compute the vector norms.  If `axis` is a 2-tuple, it specifies the
#         axes that hold 2-D matrices, and the matrix norms of these matrices
#         are computed.  If `axis` is None then either a vector norm (when `x`
#         is 1-D) or a matrix norm (when `x` is 2-D) is returned.
# The behavior should match that of `np.linalg.norm`, but it does not. The
# error is raised by the line ``x = x.tocsr()``, which occurs before the
# ``axis`` argument is checked. This line needs to be changed for the
# multi-dimensional array interface to work.
with assert_raises_regex(ValueError, "Cannot convert. CSR must be 1D or 2D. Got 3D"):
    spla.norm(A, ord=1, axis=axis)

with assert_raises_regex(ValueError, "invalid shape.*(must be 2-d)"):
    spla.onenormest(A)

# ----- norm(inv(A))
# Automatically uses last 2 dimensions
the_inv = nla.inv(A.toarray())
assert_allclose(sla.inv(A.toarray()), the_inv)

with assert_raises_regex(ValueError, "Cannot convert. CSC format must be 2D"):
    # spla.inv() uses SuperLU, which requires 2D input
    spla.inv(A)

# ----- cond(A)
# Automatically uses last 2 dimensions
the_cond = nla.cond(A.toarray(), p=1)
assert(the_cond.shape == (M,))
assert(np.all((0 < the_cond) & (the_cond < np.inf)))


# =============================================================================
# =============================================================================
