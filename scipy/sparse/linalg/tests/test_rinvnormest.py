"""Test functions for the sparse.linalg._rinvnormest module
"""

import numpy as np
from numpy.testing import assert_allclose, assert_equal, assert_
import pytest
import scipy.linalg
import scipy.sparse.linalg

class TestRinvnormest:
    @pytest.mark.parametrize("norm", [1, np.inf])
    @pytest.mark.parametrize("dtype", [np.float32, np.float64, np.complex64, np.complex128])
    def test_rinvnormest(self, norm, dtype):
        n = 5
        rng = np.random.default_rng(789002319)
        rvs = rng.random
        A = scipy.sparse.random(n, n, density=0.3, format="lil", dtype=dtype,
            random_state=rng, data_rvs=rvs)
        A = A + scipy.sparse.eye(n, format="lil") # make it likely invertible
        A = A.tocsc()
        true_rinvnorm = 1/np.linalg.norm(np.linalg.inv(A.toarray()), ord=norm)
        est_rinvnorm = scipy.sparse.linalg.rinvnormest(A, norm="1" if norm==1 else "I")
        assert_allclose(est_rinvnorm, true_rinvnorm)



