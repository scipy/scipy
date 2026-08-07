"""benchmarks for the scipy.sparse.csgraph module"""
import numpy as np
import scipy.sparse

from .common import Benchmark, get_mem_info, safe_import

with safe_import():
    from scipy.sparse.csgraph import laplacian, connected_components


def _require_memory(nbytes):
    """Skip rather than get OOM-killed on machines that can't fit the benchmark."""
    try:
        available = get_mem_info()["memavailable"]
    except Exception:
        # without psutil we can't tell; run it rather than skipping everywhere
        return
    if nbytes > available:
        raise NotImplementedError(
            f"needs ~{nbytes / 1e9:.1f} GB, only {available / 1e9:.1f} GB available"
        )


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
            # ~4.0 GB peak RSS measured; the headroom is needed because asv itself
            # shares the cgroup, so CircleCI's 4.3 GB "available" still gets OOM-killed
            _require_memory(6 * (100 * n) * 8)
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
