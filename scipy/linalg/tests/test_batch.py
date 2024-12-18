import numpy as np
from numpy.testing import assert_allclose
from scipy import linalg


class TestMatrixInScalarOut:

    def batch_test(self, fun, args=(), kwargs=None, dtype=np.float64,
                   batch_shape=(3, 2), core_shape=(4, 4), seed=8342310302941288912051):
        kwargs = {} if kwargs is None else kwargs
        rng = np.random.default_rng(seed)
        A = rng.random(batch_shape + core_shape).astype(dtype)
        res = fun(A, *args, **kwargs)

        for i in range(batch_shape[0]):
            for j in range(batch_shape[1]):
                ref = fun(A[i, j], *args, **kwargs)
                assert_allclose(res[i, j], ref)

    def test_expm_cond(self):
        self.batch_test(linalg.expm_cond)
