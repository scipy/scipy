"""benchmarks for the scipy.sparse.csgraph module"""
import numpy as np
import scipy.sparse

from .common import Benchmark, safe_import

with safe_import():
    from scipy.sparse.csgraph import dijkstra


class Dijkstra(Benchmark):
    params = [
        [30, 300, 900],
        [True, False]
    ]
    param_names = ['n', 'min_only']

    def setup(self, n, min_only):
        np.random.seed(1234)
        # make a random connectivity matrix
        data = scipy.sparse.rand(n, n, density=0.2, format='csc',
                                 random_state=42, dtype=np.bool_)
        data.setdiag(np.zeros(n, dtype=np.bool_))
        self.data = data
        # choose some random vertices
        v = np.arange(n)
        np.random.shuffle(v)
        self.indices = v[:int(n*.1)]

    def time_dijkstra_multi(self, n, min_only):
        dijkstra(self.data,
                 directed=False,
                 indices=self.indices,
                 min_only=min_only)
