import numpy as np

try:
    from scipy.cluster.hierarchy import DisjointSet
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
        self.disjoint_set = DisjointSet(self.nodes)

        self.pre_merged = DisjointSet(self.nodes)
        for a, b in self.edges:
            self.pre_merged.merge(a, b)

        self.pre_merged_found = DisjointSet(self.nodes)
        for a, b in self.edges:
            self.pre_merged_found.merge(a, b)
        for x in self.nodes:
            self.pre_merged_found[x]

    def time_merge(self, n):
        dis = self.disjoint_set
        for a, b in self.edges:
            dis.merge(a, b)

    def time_merge_already_merged(self, n):
        dis = self.pre_merged
        for a, b in self.edges:
            dis.merge(a, b)

    def time_find(self, n):
        dis = self.pre_merged
        return [dis[i] for i in self.nodes]

    def time_find_already_found(self, n):
        dis = self.pre_merged_found
        return [dis[i] for i in self.nodes]

    def time_contains(self, n):
        assert self.nodes[0] in self.pre_merged
        assert self.nodes[n // 2] in self.pre_merged
        assert self.nodes[-1] in self.pre_merged

    def time_absence(self, n):
        # Test for absence
        assert None not in self.pre_merged
        assert "dummy" not in self.pre_merged
        assert (1, 2, 3) not in self.pre_merged
