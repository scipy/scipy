"""Test the sparse.linalg._cond1est module."""

import numpy as np
import pytest

from numpy.testing import assert_allclose, assert_raises
from scipy import sparse
from scipy.sparse import linalg as spla

# TODO try different matrix sizes and densities (parametrize entire class by N)
# TODO run many trials for a given seed
# TODO try both invertible and non-invertible (ensure cond1est == np.inf)
# TODO test with empty matrices
# TODO test with single element matrices
# TODO test with identity matrix (should return 1.0)
# TODO test for failure on dense matrices
# TODO test with complex matrices

# Existing Behavior
# -----------------
# 
# - empty matrix (shape (0, 0) with 0 elements)
#
#     A = sparse.csc_array((0, 0))
# 
#     NORM(A)
# 
#     np.linalg.norm(A.toarray(), 1) = 0.0
#     spla.norm(A, 1) -> ValueError (abs(A).sum(axis=0).max() fails because
#       abs(A).sum(axis=0) produces an empty matrix)
#     spla.onenormest(A) -> ValueError (attempt to get argmax of empty
#       sequence)
# 
#     NORM(INV(A))
# 
#     np.linalg.norm(np.linalg.inv(A.toarray()), 1) = 0.0 (inverse is also an
#       empty array).
#     spla.norm(spla.inv(A), 1) -> ValueError (need at least one array to
#       concatenate in call to spsolve(A, I))
# 
#     COND(A)
# 
#     np.linalg.cond(A, p=1) -> LinAlgError 'cond is not defined on empty arrays'
# 
# - all zero matrix (shape (N, N) with 0 non-zero elements)
# 
#     A = sparse.csc_array((N, N))
#
#     NORM(A)
# 
#     np.linalg.norm(A.toarray(), 1) = 0.0
#     spla.norm(A, 1) = 0.0
#     spla.onenormest(A) = 0.0
# 
#     NORM(INV(A))
# 
#     np.linalg.norm(np.linalg.inv(A.toarray()), 1) -> LinAlgError "Singular matrix"
#     spla.norm(spla.inv(A), 1) -> RuntimeError "Factor is exactly singular"
#     spla.onenormest(spla.inv(A)) -> RuntimeError "Factor is exactly singular"
#     (spla.inv() uses SuperLU.)
# 
#     COND(A)
# 
#     np.linalg.cond(A.toarray(), p=1) = np.float64(inf)
# 
# - singleton matrix (A = [2.0])
#
#     A = sparse.csc_array([[2.0]])
# 
#     np.linalg.norm(A.toarray(), ord=1) = 2.0
#     np.linalg.norm(np.linalg.inv(A.toarray()), 1) = 0.5
#     spla.onenormest(A) = 2.0
#     spla.onenormest(spla.inv(A)) = 0.5
#     spla.norm(spla.inv(A), 1) -> TypeError input is not sparse!
#       spla.inv(A) = array([0.5]).
# 
# - exactly singular matrix (row of 0's)
# 
#     A = sparse.random_array((N, N), density=0.5, format="lil")
#     A[0] = 0  # make it exactly singular
#     A = A.tocsc()
# 
#     norms: same as above
#     np.linalg.cond(A.toarray(), p=1) = np.float64(inf)
# 
# - nearly singular matrix (row of np.finfo(float).eps / 10)
# 
#     A = sparse.random_array((N, N), density=0.5, format="lil")
#     eps = np.finfo(float).eps
#     A[0] = eps * eps  # make it nearly singular
#     A = A.tocsc()
#
#     all give same answer: np.float64(1.1716665819012448e+17)
#     np.linalg.cond(A.toarray(), p=1) = np.float64(1.6716633826598544e+17)
#


def generate_matrix(N, dtype=None, singular=False):
    """Generate a random sparse matrix of size N x N.

    Parameters
    ----------
    N : int
        Size of the matrix (N x N).
    dtype : data-type, optional
        Data type of the matrix elements (default is None, which uses
        float64).
    singular : bool, optional
        If True, the matrix will be singular (default is False, which
        creates an invertible matrix).

    Returns
    -------
    result : (N, N) sparse.csc_matrix
        A random sparse matrix in Compressed Sparse Column (CSC) format.
    """
    rng = np.random.default_rng(789002319)

    A = sparse.random_array(
        (N, N),
        density=0.5,
        format="lil",
        dtype=dtype,
        random_state=rng,
    )

    if singular:
        # NOTE we can vary this value to control the condition number.
        # * 0 would make it *exactly* singular, so cond = np.inf.
        # * ~1e±10 would make it very ill-conditioned.
        A[0] = 0  # make it exactly singular
    else:
        A.setdiag(N * (1 + np.arange(N)))  # make it invertible

    return A.tocsc()


class TestNormEstInv:
    """Test the normest_inv method of SuperLU."""

    def test_error_unsupported_norm(self):
        """Test that an error is raised for unsupported inputs."""
        A = generate_matrix(5, dtype=np.float64)
        with assert_raises(ValueError):
            spla.splu(A).normest_inv(ord=2)

    # @pytest.mark.parametrize("dtype",
    #     [np.float32, np.float64, np.complex64, np.complex128]
    # )  # FIXME complex fails
    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    @pytest.mark.parametrize("ord", [1, np.inf])
    def test_normest_inv(self, ord, dtype):
        A = generate_matrix(5, dtype)
        true_invnorm = np.linalg.norm(np.linalg.inv(A.toarray()), ord=ord)
        est_invnorm = spla.splu(A).normest_inv(ord=ord)
        assert_allclose(est_invnorm, true_invnorm)



class TestCond1Est:
    """Test the cond1est function."""

    def test_empty_matrix(self):
        """Test that an error is raised for an empty matrix."""
        pass

    # @pytest.mark.parametrize("dtype",
    #     [np.float32, np.float64, np.complex64, np.complex128]
    # )  # FIXME complex fails
    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    @pytest.mark.parametrize("singular", [True, False])
    def test_cond1est(self, dtype, singular):
        A = generate_matrix(5, dtype=dtype, singular=singular)
        true_cond1norm = np.linalg.cond(A.toarray(), p=1)
        est_cond1norm = spla.cond1est(A)
        assert_allclose(est_cond1norm, true_cond1norm)
