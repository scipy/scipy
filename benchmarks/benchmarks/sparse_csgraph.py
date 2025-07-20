"""benchmarks for the scipy.sparse.csgraph module"""
import numpy as np
import scipy.sparse

from .common import Benchmark, safe_import

with safe_import():
    from scipy.sparse.csgraph import laplacian


class Laplacian(Benchmark):
    params = [
        [30, 300, 900],
        ['dense', 'coo', 'csc', 'csr', 'dia'],
        [True, False],
        ['standard', 'opposing', 'unsigned'],
    ]
    param_names = ['n', 'format', 'normed', 'signed_graph_variant']

    def setup(self, n, format, normed, signed_graph_variant):
        data = scipy.sparse.rand(9, n, density=0.5, random_state=42).toarray()
        data = np.vstack((data, data))
        diags = list(range(-9, 0)) + list(range(1, 10))
        A = scipy.sparse.spdiags(data, diags, n, n)
        if format == 'dense':
            self.A = A.toarray()
        else:
            self.A = A.asformat(format)

    def time_laplacian(self, n, format, normed, signed_graph_variant):
        laplacian(self.A, normed=normed, signed_graph_variant=signed_graph_variant)
