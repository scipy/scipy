"""benchmarks for the scipy.sparse.csgraph module"""
import numpy as np
import scipy.sparse

from .common import Benchmark, safe_import

with safe_import():
    from scipy.sparse.csgraph import dijkstra


class Dijkstra(Benchmark):
    params = [
        [30, 300, 900],
        [True, False],
        ['random', 'star']
    ]
    param_names = ['n', 'min_only', 'format']

    def setup(self, n, min_only, format):
        rng = np.random.default_rng(1234)
        if format == 'random':
            # make a random connectivity matrix
            data = scipy.sparse.rand(n, n, density=0.2, format='csc',
                                     random_state=42, dtype=np.bool_)
            data.setdiag(np.zeros(n, dtype=np.bool_))
            self.data = data
        elif format == 'star':
            rows = [0 for i in range(n - 1)] + [i + 1 for i in range(n - 1)]
            cols = [i + 1 for i in range(n - 1)] + [0 for i in range(n - 1)]
            weights = [i + 1 for i in range(n - 1)] * 2
            self.data = scipy.sparse.csr_matrix((weights, (rows, cols)),
                                                shape=(n, n))
        # choose some random vertices
        v = np.arange(n)
        rng.shuffle(v)
        self.indices = v[:int(n*.1)]

    def time_dijkstra_multi(self, n, min_only, format):
        dijkstra(self.data,
                 directed=False,
                 indices=self.indices,
                 min_only=min_only)
