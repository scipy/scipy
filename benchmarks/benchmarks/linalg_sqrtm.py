""" Benchmark linalg.sqrtm for various blocksizes.

"""
import numpy as np

try:
    import scipy.linalg
except ImportError:
    pass

from .common import Benchmark


class Sqrtm(Benchmark):
    params = [
        ['float64', 'complex128'],
        [64, 256],
        [32, 64, 256]
    ]
    param_names = ['dtype', 'n', 'blocksize']

    def setup(self, dtype, n, blocksize):
        n = int(n)
        dtype = np.dtype(dtype)
        blocksize = int(blocksize)
        A = np.random.rand(n, n)
        if dtype == np.complex128:
            A = A + 1j*np.random.rand(n, n)
        self.A = A

        if blocksize > n:
            raise NotImplementedError()

    def time_sqrtm(self, dtype, n, blocksize):
        scipy.linalg.sqrtm(self.A, disp=False, blocksize=blocksize)
