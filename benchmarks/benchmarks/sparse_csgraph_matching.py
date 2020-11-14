import numpy as np
import scipy.sparse

from .common import Benchmark, safe_import

with safe_import():
    from scipy.sparse.csgraph import maximum_bipartite_matching


class MaximumBipartiteMatching(Benchmark):
    params = [[5000, 7500, 10000], [0.0001, 0.0005, 0.001]]
    param_names = ['n', 'density']

    def setup(self, n, density):
        # Create random sparse matrices. Note that we could use
        # scipy.sparse.rand for this purpose, but simply using np.random and
        # disregarding duplicates is quite a bit faster.
        np.random.seed(42)
        d = np.random.randint(0, n, size=(int(n*n*density), 2))
        graph = scipy.sparse.csr_matrix((np.ones(len(d)), (d[:, 0], d[:, 1])),
                                        shape=(n, n))
        self.graph = graph

    def time_maximum_bipartite_matching(self, n, density):
        maximum_bipartite_matching(self.graph)
