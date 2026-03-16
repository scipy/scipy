import numpy as np
from scipy.sparse.linalg import tfqmr
import pytest

def test_tfqmr_zero_iter_show_true_basic():
    # Dense matrix case: expect ValueError when maxiter=0 and show=True (buggy behavior)
    A = np.eye(3)
    b = np.ones(3)
    with pytest.raises(ValueError, match='maxiter must be an integer not less than'):
        tfqmr(A, b, maxiter=0, show=True)