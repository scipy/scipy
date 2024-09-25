"""Test functions for the sparse.linalg._rinvnormest module
"""

import numpy as np
from numpy.testing import assert_allclose, assert_raises
import pytest
import scipy.linalg
import scipy.sparse.linalg

class TestRinvnormest:
    def generate_matrix(self, n, dtype):
        rng = np.random.default_rng(789002319)
        rvs = rng.random
        A = scipy.sparse.random(n, n, density=0.3, format="lil", dtype=dtype,
            random_state=rng, data_rvs=rvs)
        A = A + scipy.sparse.eye(n, format="lil", dtype=dtype) # make it likely invertible
        A = A.tocsc()
        return A

    @pytest.mark.parametrize("norm", [1, np.inf])
    @pytest.mark.parametrize("dtype", [np.float32, np.float64, np.complex64, np.complex128])
    def test_rinvnormest(self, norm, dtype):
        A = self.generate_matrix(5, dtype)
        true_rinvnorm = 1/np.linalg.norm(np.linalg.inv(A.toarray()), ord=norm)
        est_rinvnorm = scipy.sparse.linalg.rinvnormest(A, norm="1" if norm==1 else "I")
        assert_allclose(est_rinvnorm, true_rinvnorm)

    @pytest.mark.parametrize("dtype", [np.float32, np.float64, np.complex64, np.complex128])
    def test_cond1normest(self, dtype):
        A = self.generate_matrix(5, dtype)
        true_cond1norm = np.linalg.cond(A.toarray(), p=1)
        est_cond1norm = scipy.sparse.linalg.cond1est(A)
        assert_allclose(est_cond1norm, true_cond1norm)

    def test_error_unsupported_norm(self):
        A = self.generate_matrix(5, np.float64)
        with assert_raises(ValueError):
            scipy.sparse.linalg.splu(A).rinvnormest(norm="2")
        


