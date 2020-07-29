import numpy as np

try:
    from scipy.sparse.csgraph import DisjointSet
except ImportError:
    pass

from .common import Benchmark


class Bench(Benchmark):
    params = [[100, 1000, 10000]]
    param_names = ['n']

    def setup(self, n):
        # Create random edges
        rng = np.random.RandomState(seed=0)
        self.edges = rng.randint(0, 10 * n, (n, 2))
        self.nodes = np.unique(self.edges)
        self.disjoint_set = DisjointSet()

    def time_merge_find(self, n):
        dis = self.disjoint_set
        for a, b in self.edges:
            dis.merge(a, b)

        return [dis[i] for i in self.nodes]
