"""benchmarks for the scipy.sparse.csgraph module"""
from __future__ import print_function, absolute_import
import numpy as np
import scipy.sparse

try:
    from scipy.sparse.csgraph import dijkstra
except ImportError:
    pass

from .common import Benchmark


class Dijkstra(Benchmark):
    params = [
        [30, 300, 900]
    ]
    param_names = ['n', 'format', 'normed']

    def setup(self, n):
        # make a random connectivity matrix
        data = scipy.sparse.rand(n, n, density=0.2, format='csc', random_state=42, dtype=np.bool)
        data.setdiag(np.zeros(n, dtype=np.bool))
        self.data = data
        # choose some random vertices
        v = np.arange(n)
        np.random.shuffle(v)
        self.indices = v[:int(n*.1)]

    def time_dijkstra_single(self, n):
        dm = dijkstra(self.data,
                      directed=False,
                      indices=self.indices,
                      min_only=False)
        ds = np.min(dm, axis=1) 

    def time_dijkstra_multi(self, n):
        ds = dijkstra(self.data,
                      directed=False,
                      indices=self.indices,
                      min_only=True)
