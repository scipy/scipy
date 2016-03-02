import numpy as np

from .common import Benchmark

try:
    from scipy.cluster.hierarchy import linkage
except ImportError:
    pass


class HierarchyLinkage(Benchmark):
    params = ['single', 'complete', 'average', 'weighted', 'centroid',
              'median', 'ward']
    param_names = ['method']

    def __init__(self):
        rnd = np.random.RandomState(0)
        self.X = rnd.randn(2000, 2)

    def time_linkage(self, method):
        linkage(self.X, method=method)
