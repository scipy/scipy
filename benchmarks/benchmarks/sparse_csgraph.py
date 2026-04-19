"""benchmarks for the scipy.sparse.csgraph module"""
import numpy as np
import scipy.sparse

from .common import Benchmark, safe_import

with safe_import():
    from scipy.sparse.csgraph import laplacian, connected_components


class Laplacian(Benchmark):
    params = [
        [30, 300, 900],
        ['dense', 'coo', 'csc', 'csr', 'dia'],
        [True, False]
    ]
    param_names = ['n', 'format', 'normed']

    def setup(self, n, format, normed):
        data = scipy.sparse.rand(9, n, density=0.5, random_state=42).toarray()
        data = np.vstack((data, data))
        diags = list(range(-9, 0)) + list(range(1, 10))
        A = scipy.sparse.spdiags(data, diags, n, n)
        if format == 'dense':
            self.A = A.toarray()
        else:
            self.A = A.asformat(format)

    def time_laplacian(self, n, format, normed):
        laplacian(self.A, normed=normed)

class StronglyConnectedComponents(Benchmark):
    params = [["random", "single_scc", "chain"]]
    param_names = ["kind"]

    def setup(self, kind):
        n = 1_000_000
        rng = np.random.default_rng(42)
        if kind == "random":
            self.G = scipy.sparse.random_array(
                shape=(n, n),
                density=100 / n,
                format="csr",
                rng=rng,
            )
        elif kind == "single_scc":
            # Hamiltonian cycle (one giant SCC) plus random edges.
            perm = rng.permutation(n)
            row = np.concatenate([perm, rng.integers(0, n, size=99 * n)])
            col = np.concatenate([np.roll(perm, -1), rng.integers(0, n, size=99 * n)])
            self.G = scipy.sparse.csr_array(
                (np.ones(len(row)), (row, col)), shape=(n, n)
            )
        elif kind == "chain":
            row = np.arange(n - 1)
            self.G = scipy.sparse.csr_array(
                (np.ones(n - 1), (row, row + 1)), shape=(n, n)
            )

    def time_strongly_connected_components(self, kind):
        connected_components(self.G, directed=True, connection="strong")
