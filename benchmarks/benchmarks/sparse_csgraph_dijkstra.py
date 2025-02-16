"""benchmarks for the scipy.sparse.csgraph module"""
import numpy as np
import scipy.sparse

from .common import Benchmark, safe_import, is_xslow

with safe_import():
    from scipy.sparse.csgraph import dijkstra, shortest_path


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
            data = scipy.sparse.rand(n, n, density=0.2, format='lil',
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


class DijkstraDensity(Benchmark):
    """
    Benchmark performance of Dijkstra, adapted from [^1]

    [^1]: https://github.com/scipy/scipy/pull/20717#issuecomment-2562795171
    """
    params = [
        [10, 100, 1000],
        [0.1, 0.3, 0.5, 0.9],
    ]
    param_names = ["n", "density"]

    def setup(self, n, density):
        if n >= 1000 and not is_xslow():
            raise NotImplementedError("skipped")

        rng = np.random.default_rng(42)
        self.graph = scipy.sparse.random_array(
            shape=(n, n),
            density=density,
            format='csr',
            rng=rng,
            data_sampler=lambda size: rng.integers(100, size=size, dtype=np.uint32),
        )


    def time_test_shortest_path(self, n, density):
        shortest_path(self.graph, method="D", directed=False)
