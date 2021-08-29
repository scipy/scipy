import numpy as np
from numpy.testing import suppress_warnings

from .common import Benchmark, safe_import

with safe_import():
    from scipy.cluster.hierarchy import linkage
    from scipy.cluster.vq import kmeans, kmeans2, vq


class HierarchyLinkage(Benchmark):
    params = ['single', 'complete', 'average', 'weighted', 'centroid',
              'median', 'ward']
    param_names = ['method']

    def __init__(self):
        rnd = np.random.RandomState(0)
        self.X = rnd.randn(2000, 2)

    def time_linkage(self, method):
        linkage(self.X, method=method)


class KMeans(Benchmark):
    params = [2, 10, 50]
    param_names = ['k']

    def __init__(self):
        rnd = np.random.RandomState(0)
        self.obs = rnd.rand(1000, 5)

    def time_kmeans(self, k):
        kmeans(self.obs, k, iter=10)


class KMeans2(Benchmark):
    params = [[2, 10, 50], ['random', 'points', '++']]
    param_names = ['k', 'init']

    def __init__(self):
        rnd = np.random.RandomState(0)
        self.obs = rnd.rand(1000, 5)

    def time_kmeans2(self, k, init):
        with suppress_warnings() as sup:
            sup.filter(UserWarning,
                       "One of the clusters is empty. Re-run kmeans with a "
                       "different initialization")
            kmeans2(self.obs, k, minit=init, iter=10)


class VQ(Benchmark):
    params = [[2, 10, 50], ['float32', 'float64', 'float128']]
    param_names = ['k', 'dtype']

    def __init__(self):
        rnd = np.random.RandomState(0)
        self.data = rnd.rand(5000, 5)
        self.cbook_source = rnd.rand(50, 5)

    def setup(self, k, dtype):
        self.obs = self.data.astype(dtype)
        self.cbook = self.cbook_source[:k].astype(dtype)

    def time_vq(self, k, dtype):
        vq(self.obs, self.cbook)
