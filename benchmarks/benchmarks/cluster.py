import numpy as np
from numpy.testing import suppress_warnings

from .common import Benchmark, XPBenchmark, is_xslow, safe_import

with safe_import():
    from scipy.cluster.hierarchy import linkage, is_isomorphic
    from scipy.cluster.vq import kmeans, kmeans2, vq, whiten


class HierarchyLinkage(Benchmark):
    params = ['single', 'complete', 'average', 'weighted', 'centroid',
              'median', 'ward']
    param_names = ['method']

    def __init__(self):
        rnd = np.random.RandomState(0)
        self.X = rnd.randn(2000, 2)

    def time_linkage(self, method):
        linkage(self.X, method=method)


class IsIsomorphic(XPBenchmark):
    NCLUSTERS = 5
    # This is very slow and memory intensive, but necessary to
    # let _most_ backends approach O(n*logn) behaviour.
    # Note: memory usage = 16 * nobs
    if is_xslow():
        nobs = [100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]
    else:
        nobs = [100, 100_000]
    # Skip cpu backends for nobs greater than this. 
    # They should all have reached O(n*logn) behaviour by then.
    CPU_MAX_OBS = 1_000_000

    param_names = (*XPBenchmark.param_names, "nobs")
    params = (*XPBenchmark.params, nobs)

    def setup(self, backend, nobs):
        use_cuda = backend == "cupy" or backend.endswith(":cuda")
        if not use_cuda and nobs > self.CPU_MAX_OBS:
            raise NotImplementedError("Skipping huge size on CPU")

        super().setup(backend, is_isomorphic)

        rng = np.random.default_rng(0)
        a = self.xp.asarray(rng.integers(0, self.NCLUSTERS, size=nobs))
        p = self.xp.asarray(rng.permutation(self.NCLUSTERS))
        b = self.xp.take(p, a)
        self.a, self.b = self.synchronize(a, b)

        if self.warmup:
            self.func(self.a, self.b)

    def time_is_isomorphic(self, backend, nobs):
        self.func(self.a, self.b)


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
    params = [[2, 10, 50], ['float32', 'float64']]
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


class Whiten(XPBenchmark):
    # This is very slow and memory intensive, but necessary to
    # let _most_ backends approach O(n*logn) behaviour.
    # Note: memory usage = 8 * shape[0] * shape[1]
    if is_xslow():
        shape = [(10, 10), (100, 100), (1000, 1000), (10_000, 10_000)]
    else:
        shape = [(10, 10), (100, 100)]
    # Skip cpu backends for shape greater than this. 
    # They should all have reached O(n*logn) behaviour by then.
    CPU_MAX_OBS = (1000, 1000)

    param_names = (*XPBenchmark.param_names, "shape")
    params = (*XPBenchmark.params, shape)

    def setup(self, backend, shape):
        super().setup(backend, whiten, static_argnames="check_finite")

        rng = np.random.default_rng(0)
        obs = self.xp.asarray(rng.uniform(0, 100.0, size=shape))
        self.obs = self.synchronize(obs)

        if self.warmup:
            self.func(self.obs, check_finite=False)

    def time_whiten(self, backend, shape):
        self.func(self.obs, check_finite=False)
