import scipy.sparse

try:
    from scipy.sparse.csgraph import maximum_bipartite_matching
except ImportError:
    pass

from .common import Benchmark


class MaximumBipartiteMatching(Benchmark):
    params = [[1000, 2500], [0.01, 0.1]]
    param_names = ['n', 'density']

    def setup(self, n, density):
        graph = scipy.sparse.rand(n, n, density=density,
                                  format='csr', random_state=42)
        self.graph = graph

    def time_maximum_bipartite_matching(self, n, density):
        maximum_bipartite_matching(self.graph)
