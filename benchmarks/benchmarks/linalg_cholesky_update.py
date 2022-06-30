from .common import Benchmark, safe_import
import numpy as np
import scipy.linalg

with safe_import():
    from scipy.linalg._cholesky_update import cholesky_update


class CholUpdate(Benchmark):
    sizes = [100, 500, 1000, 2000]
    shapes = [(i, i) for i in sizes]
    orders = ['C', 'F']
    param_names = ['shape', 'order']
    params = [shapes, orders]

    def setup(self, shape, order):
        np.random.seed(29382)
        a = np.random.random(shape)
        a = a + a.conj().T + a.shape[0] * np.eye(a.shape[0])
        u1 = np.random.random((a.shape[0], 1))
        u2 = np.random.random((a.shape[0], 2))
        u20 = np.random.random((a.shape[0], 20))

        self.a = a.copy(order=order)
        self.r = scipy.linalg.cholesky(a, lower=False).copy(order=order)
        self.u1 = u1.copy(order=order)
        self.u2 = u2.copy(order=order)
        self.u20 = u20.copy(order=order)

    def time_update1(self, *args):
        cholesky_update(self.r, self.u1, downdate=False, lower=False,
                        overwrite_rz=False, check_finite=False)

    def time_update20(self, *args):
        cholesky_update(self.r, self.u2, downdate=False, lower=False,
                        overwrite_rz=False, check_finite=False)

    def time_cholesky(self, *args):
        # overwrite = True will fail with LinAlgErrors
        # unclear why since setup should always be called.
        scipy.linalg.cholesky(self.a, lower=False, check_finite=False,
                              overwrite_a=False)
