import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    import scipy.linalg.lapack as la
    import scipy.linalg.blas as bla


class GetBlasLapackFuncs(Benchmark):
    """
    Test the speed of grabbing the correct BLAS/LAPACK routine flavor.

    In particular, upon receiving strange dtype arrays the results shouldn't
    diverge too much. Hence the results here should be comparable
    """

    param_names = ['dtype1', 'dtype2',
                   'dtype1_ord', 'dtype2_ord',
                   'size']
    params = [
        ['b', 'G', 'd'],
        ['d', 'F', '?'],
        ['C', 'F'],
        ['C', 'F'],
        [10, 100, 1000]
    ]

    def setup(self, dtype1, dtype2, dtype1_ord, dtype2_ord, size):
        self.arr1 = np.empty(size, dtype=dtype1, order=dtype1_ord)
        self.arr2 = np.empty(size, dtype=dtype2, order=dtype2_ord)

    def time_find_best_blas_type(self, dtype1, dtype2, dtype1_ord, dtype2_ord, size):
        prefix, dtype, prefer_fortran = bla.find_best_blas_type((self.arr1, self.arr2))
