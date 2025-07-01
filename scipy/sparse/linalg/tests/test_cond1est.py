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

class TestCond1Est:
    @staticmethod
    def generate_matrix(n, dtype=None, singular=False):
        """Generate a random sparse matrix of size n x n.

        Parameters
        ----------
        n : int
            Size of the matrix (n x n).
        dtype : data-type, optional
            Data type of the matrix elements (default is None, which uses
            float64).
        singular : bool, optional
            If True, the matrix will be singular (default is False, which
            creates an invertible matrix).

        Returns
        -------
        result : (n, n) sparse.csc_matrix
            A random sparse matrix in Compressed Sparse Column (CSC) format.
        """
        rng = np.random.default_rng(789002319)

        A = sparse.random_array(
            (n, n),
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
            A.setdiag(n * (1 + np.arange(n)))  # make it invertible

        return A.tocsc()

    # @pytest.mark.parametrize("dtype",
    #     [np.float32, np.float64, np.complex64, np.complex128]
    # )  # FIXME complex fails
    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    @pytest.mark.parametrize("ord", [1, np.inf])
    def test_invnormest(self, ord, dtype):
        A = self.generate_matrix(5, dtype)
        true_invnorm = np.linalg.norm(np.linalg.inv(A.toarray()), ord=ord)
        est_invnorm = spla.splu(A).normest_inv(ord=ord)
        assert_allclose(est_invnorm, true_invnorm)

    # @pytest.mark.parametrize("dtype",
    #     [np.float32, np.float64, np.complex64, np.complex128]
    # )  # FIXME complex fails
    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    @pytest.mark.parametrize("singular", [True, False])
    def test_cond1est(self, dtype, singular):
        A = self.generate_matrix(5, dtype=dtype, singular=singular)
        true_cond1norm = np.linalg.cond(A.toarray(), p=1)
        est_cond1norm = spla.cond1est(A)
        assert_allclose(est_cond1norm, true_cond1norm)

    def test_error_unsupported_norm(self):
        A = self.generate_matrix(5, np.float64)
        with assert_raises(ValueError):
            spla.splu(A).normest_inv(ord=2)
        
