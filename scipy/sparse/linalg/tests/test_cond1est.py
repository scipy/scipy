"""Test functions for the sparse.linalg._cond1est module
"""

import numpy as np
from numpy.testing import assert_allclose, assert_raises
import pytest
import scipy.linalg
import scipy.sparse.linalg

class TestCond1Est:
    def generate_matrix(self, n, dtype):
        rng = np.random.default_rng(789002319)
        rvs = rng.random
        A = scipy.sparse.random(n, n, density=0.3, format="lil", dtype=dtype,
            random_state=rng, data_rvs=rvs)
        if dtype==np.complex64 or dtype==np.complex128:
          A += 1j * scipy.sparse.random(n, n, density=0.3, format="lil", dtype=dtype,
              random_state=rng, data_rvs=rvs)
        A = A + scipy.sparse.eye(n, format="lil", dtype=dtype) # make it likely invertible
        A = A.tocsc()
        return A

    @pytest.mark.parametrize("norm", [1, np.inf])
    @pytest.mark.parametrize("dtype", [np.float32, np.float64, np.complex64, np.complex128])
    def test_invnormest(self, norm, dtype):
        A = self.generate_matrix(5, dtype)
        true_invnorm = np.linalg.norm(np.linalg.inv(A.toarray()), ord=norm)
        est_invnorm = scipy.sparse.linalg.invnormest(A, norm="1" if norm==1 else "I")
        assert_allclose(est_invnorm, true_invnorm)

    @pytest.mark.parametrize("dtype", [np.float32, np.float64, np.complex64, np.complex128])
    def test_cond1normest(self, dtype):
        A = self.generate_matrix(5, dtype)
        true_cond1norm = np.linalg.cond(A.toarray(), p=1)
        est_cond1norm = scipy.sparse.linalg.cond1est(A)
        assert_allclose(est_cond1norm, true_cond1norm, rtol=1e-4)

    def test_error_unsupported_norm(self):
        A = self.generate_matrix(5, np.float64)
        with assert_raises(ValueError):
            scipy.sparse.linalg.splu(A).invnormest(norm="2")
        


