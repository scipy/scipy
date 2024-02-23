import numpy as np

import scipy.sparse
from .common import Benchmark, safe_import

with safe_import():
    from scipy.sparse.csgraph import maximum_flow


class MaximumFlow(Benchmark):
    params = [[200, 500, 1500], [0.1, 0.3, 0.5]]
    param_names = ['n', 'density']

    def setup(self, n, density):
        # Create random matrices whose values are integers between 0 and 100.
        data = (scipy.sparse.rand(n, n, density=density, format='lil',
                                  random_state=42)*100).astype(np.int32)
        data.setdiag(np.zeros(n, dtype=np.int32))
        self.data = scipy.sparse.csr_matrix(data)

    def time_maximum_flow(self, n, density):
        maximum_flow(self.data, 0, n - 1)
