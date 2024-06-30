"""
Simple benchmarks for the sparse module
"""
import warnings
import time
import timeit
import pickle

import numpy as np
from numpy import ones, array, asarray, empty

from .common import Benchmark, safe_import

with safe_import():
    from scipy import sparse
    from scipy.sparse import (coo_matrix, dia_matrix, lil_matrix,
                              dok_matrix, rand, SparseEfficiencyWarning)


def random_sparse(m, n, nnz_per_row):
    rows = np.arange(m).repeat(nnz_per_row)
    cols = np.random.randint(0, n, size=nnz_per_row*m)
    vals = np.random.random_sample(m*nnz_per_row)
    return coo_matrix((vals, (rows, cols)), (m, n)).tocsr()


# TODO move this to a matrix gallery and add unittests
def poisson2d(N, dtype='d', format=None):
    """
    Return a sparse matrix for the 2D Poisson problem
    with standard 5-point finite difference stencil on a
    square N-by-N grid.
    """
    if N == 1:
        diags = asarray([[4]], dtype=dtype)
        return dia_matrix((diags, [0]), shape=(1, 1)).asformat(format)

    offsets = array([0, -N, N, -1, 1])

    diags = empty((5, N**2), dtype=dtype)

    diags[0] = 4  # main diagonal
    diags[1:] = -1  # all offdiagonals

    diags[3, N-1::N] = 0  # first lower diagonal
    diags[4, N::N] = 0  # first upper diagonal

    return dia_matrix((diags, offsets), shape=(N**2, N**2)).asformat(format)


class Arithmetic(Benchmark):
    param_names = ['format', 'XY', 'op']
    params = [
        ['csr', 'csc', 'coo', 'dia'],
        ['AA', 'AB', 'BA', 'BB'],
        ['__add__', '__sub__', 'multiply', '__mul__']
    ]

    def setup(self, format, XY, op):
        matrices = dict(A=poisson2d(250, format=format),
                        B=poisson2d(250, format=format)**2)

        x = matrices[XY[0]]
        self.y = matrices[XY[1]]
        self.fn = getattr(x, op)
        self.fn(self.y)  # warmup

    def time_arithmetic(self, format, XY, op):
        self.fn(self.y)


class Sort(Benchmark):
    params = ['Rand10', 'Rand25', 'Rand50', 'Rand100', 'Rand200']
    param_names = ['matrix']

    def setup(self, matrix):
        n = 10000
        if matrix.startswith('Rand'):
            k = int(matrix[4:])
            self.A = random_sparse(n, n, k)
            self.A.has_sorted_indices = False
            self.A.indices[:2] = 2, 1
        else:
            raise NotImplementedError()

    def time_sort(self, matrix):
        """sort CSR column indices"""
        self.A.sort_indices()


class Matvec(Benchmark):
    params = [
        ['Identity', 'Poisson5pt', 'Block2x2', 'Block3x3'],
        ['dia', 'csr', 'csc', 'dok', 'lil', 'coo', 'bsr']
    ]
    param_names = ['matrix', 'format']

    def setup(self, matrix, format):
        if matrix == 'Identity':
            if format in ('lil', 'dok'):
                raise NotImplementedError()
            self.A = sparse.eye(10000, 10000, format=format)
        elif matrix == 'Poisson5pt':
            self.A = poisson2d(300, format=format)
        elif matrix == 'Block2x2':
            if format not in ('csr', 'bsr'):
                raise NotImplementedError()
            b = (2, 2)
            self.A = sparse.kron(poisson2d(150),
                                 ones(b)).tobsr(blocksize=b).asformat(format)
        elif matrix == 'Block3x3':
            if format not in ('csr', 'bsr'):
                raise NotImplementedError()
            b = (3, 3)
            self.A = sparse.kron(poisson2d(100),
                                 ones(b)).tobsr(blocksize=b).asformat(format)
        else:
            raise NotImplementedError()

        self.x = ones(self.A.shape[1], dtype=float)

    def time_matvec(self, matrix, format):
        self.A * self.x


class Matvecs(Benchmark):
    params = ['dia', 'coo', 'csr', 'csc', 'bsr']
    param_names = ["format"]

    def setup(self, format):
        self.A = poisson2d(300, format=format)
        self.x = ones((self.A.shape[1], 10), dtype=self.A.dtype)

    def time_matvecs(self, format):
        self.A * self.x


class Matmul(Benchmark):
    def setup(self):
        H1, W1 = 1, 100000
        H2, W2 = W1, 1000
        C1 = 10
        C2 = 1000000

        rng = np.random.default_rng(0)

        i = rng.integers(H1, size=C1)
        j = rng.integers(W1, size=C1)
        data = rng.random(C1)
        self.matrix1 = coo_matrix((data, (i, j)), shape=(H1, W1)).tocsr()

        i = rng.integers(H2, size=C2)
        j = rng.integers(W2, size=C2)
        data = rng.random(C2)
        self.matrix2 = coo_matrix((data, (i, j)), shape=(H2, W2)).tocsr()

    def time_large(self):
        for i in range(100):
            self.matrix1 * self.matrix2

    # Retain old benchmark results (remove this if changing the benchmark)
    time_large.version = (
        "33aee08539377a7cb0fabaf0d9ff9d6d80079a428873f451b378c39f6ead48cb"
    )


class Construction(Benchmark):
    params = [
        ['Empty', 'Identity', 'Poisson5pt'],
        ['lil', 'dok']
    ]
    param_names = ['matrix', 'format']

    def setup(self, name, format):
        if name == 'Empty':
            self.A = coo_matrix((10000, 10000))
        elif name == 'Identity':
            self.A = sparse.eye(10000, format='coo')
        else:
            self.A = poisson2d(100, format='coo')

        formats = {'lil': lil_matrix, 'dok': dok_matrix}
        self.cls = formats[format]

    def time_construction(self, name, format):
        T = self.cls(self.A.shape)
        for i, j, v in zip(self.A.row, self.A.col, self.A.data):
            T[i, j] = v


class BlockDiagDenseConstruction(Benchmark):
    param_names = ['num_matrices']
    params = [1000, 5000, 10000, 15000, 20000]

    def setup(self, num_matrices):
        self.matrices = []
        for i in range(num_matrices):
            rows = np.random.randint(1, 4)
            columns = np.random.randint(1, 4)
            mat = np.random.randint(0, 10, (rows, columns))
            self.matrices.append(mat)

    def time_block_diag(self, num_matrices):
        sparse.block_diag(self.matrices)


class BlockDiagSparseConstruction(Benchmark):
    param_names = ['num_matrices']
    params = [100, 500, 1000, 1500, 2000]

    def setup(self, num_matrices):
        self.matrices = []
        for i in range(num_matrices):
            rows = np.random.randint(1, 20)
            columns = np.random.randint(1, 20)
            mat = np.random.randint(0, 10, (rows, columns))
            self.matrices.append(mat)

    def time_block_diag(self, num_matrices):
        sparse.block_diag(self.matrices)


class CsrHstack(Benchmark):
    param_names = ['num_rows']
    params = [10000, 25000, 50000, 100000, 250000]

    def setup(self, num_rows):
        num_cols = int(1e5)
        density = 2e-3
        nnz_per_row = int(density*num_cols)
        self.mat = random_sparse(num_rows, num_cols, nnz_per_row)

    def time_csr_hstack(self, num_rows):
        sparse.hstack([self.mat, self.mat])


class Conversion(Benchmark):
    params = [
        ['csr', 'csc', 'coo', 'dia', 'lil', 'dok', 'bsr'],
        ['csr', 'csc', 'coo', 'dia', 'lil', 'dok', 'bsr'],
    ]
    param_names = ['from_format', 'to_format']

    def setup(self, fromfmt, tofmt):
        base = poisson2d(100, format=fromfmt)

        try:
            self.fn = getattr(base, 'to' + tofmt)
        except Exception:
            def fn():
                raise RuntimeError()
            self.fn = fn

    def time_conversion(self, fromfmt, tofmt):
        self.fn()


class Getset(Benchmark):
    params = [
        [1, 10, 100, 1000, 10000],
        ['different', 'same'],
        ['csr', 'csc', 'lil', 'dok']
    ]
    param_names = ['N', 'sparsity pattern', 'format']
    unit = "seconds"

    def setup(self, N, sparsity_pattern, format):
        if format == 'dok' and N > 500:
            raise NotImplementedError()

        self.A = rand(1000, 1000, density=1e-5)

        A = self.A
        N = int(N)

        # indices to assign to
        i, j = [], []
        while len(i) < N:
            n = N - len(i)
            ip = np.random.randint(0, A.shape[0], size=n)
            jp = np.random.randint(0, A.shape[1], size=n)
            i = np.r_[i, ip]
            j = np.r_[j, jp]
        v = np.random.rand(n)

        if N == 1:
            i = int(i)
            j = int(j)
            v = float(v)

        base = A.asformat(format)

        self.m = base.copy()
        self.i = i
        self.j = j
        self.v = v

    def _timeit(self, kernel, recopy):
        min_time = 1e99
        if not recopy:
            kernel(self.m, self.i, self.j, self.v)

        number = 1
        start = time.time()
        while time.time() - start < 0.1:
            if recopy:
                m = self.m.copy()
            else:
                m = self.m
            while True:
                duration = timeit.timeit(
                    lambda: kernel(m, self.i, self.j, self.v), number=number)
                if duration > 1e-5:
                    break
                else:
                    number *= 10
            min_time = min(min_time, duration/number)
        return min_time

    def track_fancy_setitem(self, N, sparsity_pattern, format):
        def kernel(A, i, j, v):
            A[i, j] = v

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SparseEfficiencyWarning)
            return self._timeit(kernel, sparsity_pattern == 'different')

    def time_fancy_getitem(self, N, sparsity_pattern, format):
        self.m[self.i, self.j]


class NullSlice(Benchmark):
    params = [[0.05, 0.01], ['csr', 'csc', 'lil']]
    param_names = ['density', 'format']

    def _setup(self, density, format):
        n = 100000
        k = 1000

        # faster version of rand(n, k, format=format, density=density),
        # with non-exact nnz
        nz = int(n*k * density)
        row = np.random.randint(0, n, size=nz)
        col = np.random.randint(0, k, size=nz)
        data = np.ones(nz, dtype=np.float64)
        X = coo_matrix((data, (row, col)), shape=(n, k))
        X.sum_duplicates()
        X = X.asformat(format)
        with open(f'{density}-{format}.pck', 'wb') as f:
            pickle.dump(X, f, protocol=pickle.HIGHEST_PROTOCOL)

    def setup_cache(self):
        for density in self.params[0]:
            for fmt in self.params[1]:
                self._setup(density, fmt)

    setup_cache.timeout = 120

    def setup(self, density, format):
        # Unpickling is faster than computing the random matrix...
        with open(f'{density}-{format}.pck', 'rb') as f:
            self.X = pickle.load(f)

    def time_getrow(self, density, format):
        self.X.getrow(100)

    def time_getcol(self, density, format):
        self.X.getcol(100)

    def time_3_rows(self, density, format):
        self.X[[0, 100, 105], :]

    def time_10000_rows(self, density, format):
        self.X[np.arange(10000), :]

    def time_3_cols(self, density, format):
        self.X[:, [0, 100, 105]]

    def time_100_cols(self, density, format):
        self.X[:, np.arange(100)]

    # Retain old benchmark results (remove this if changing the benchmark)
    time_10000_rows.version = (
        "dc19210b894d5fd41d4563f85b7459ef5836cddaf77154b539df3ea91c5d5c1c"
    )
    time_100_cols.version = (
        "8d43ed52084cdab150018eedb289a749a39f35d4dfa31f53280f1ef286a23046"
    )
    time_3_cols.version = (
        "93e5123910772d62b3f72abff56c2732f83d217221bce409b70e77b89c311d26"
    )
    time_3_rows.version = (
        "a9eac80863a0b2f4b510269955041930e5fdd15607238257eb78244f891ebfe6"
    )
    time_getcol.version = (
        "291388763b355f0f3935db9272a29965d14fa3f305d3306059381e15300e638b"
    )
    time_getrow.version = (
        "edb9e4291560d6ba8dd58ef371b3a343a333bc10744496adb3ff964762d33c68"
    )


class Diagonal(Benchmark):
    params = [[0.01, 0.1, 0.5], ['csr', 'csc', 'coo', 'lil', 'dok', 'dia']]
    param_names = ['density', 'format']

    def setup(self, density, format):
        n = 1000
        if format == 'dok' and n * density >= 500:
            raise NotImplementedError()

        warnings.simplefilter('ignore', SparseEfficiencyWarning)

        self.X = sparse.rand(n, n, format=format, density=density)

    def time_diagonal(self, density, format):
        self.X.diagonal()

    # Retain old benchmark results (remove this if changing the benchmark)
    time_diagonal.version = (
        "d84f53fdc6abc208136c8ce48ca156370f6803562f6908eb6bd1424f50310cf1"
    )


class Sum(Benchmark):
    params = [[0.01, 0.1, 0.5], ['csr', 'csc', 'coo', 'lil', 'dok', 'dia']]
    param_names = ['density', 'format']

    def setup(self, density, format):
        n = 1000
        if format == 'dok' and n * density >= 500:
            raise NotImplementedError()

        warnings.simplefilter('ignore', SparseEfficiencyWarning)
        self.X = sparse.rand(n, n, format=format, density=density)

    def time_sum(self, density, format):
        self.X.sum()

    def time_sum_axis0(self, density, format):
        self.X.sum(axis=0)

    def time_sum_axis1(self, density, format):
        self.X.sum(axis=1)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_sum.version = (
        "05c305857e771024535e546360203b17f5aca2b39b023a49ab296bd746d6cdd3"
    )
    time_sum_axis0.version = (
        "8aca682fd69aa140c69c028679826bdf43c717589b1961b4702d744ed72effc6"
    )
    time_sum_axis1.version = (
        "1a6e05244b77f857c61f8ee09ca3abd006a10ba07eff10b1c5f9e0ac20f331b2"
    )


class Iteration(Benchmark):
    params = [[0.05, 0.01], ['csr', 'csc', 'lil']]
    param_names = ['density', 'format']

    def setup(self, density, format):
        n = 500
        k = 1000
        self.X = sparse.rand(n, k, format=format, density=density)

    def time_iteration(self, density, format):
        for row in self.X:
            pass


class Densify(Benchmark):
    params = [
        ['dia', 'csr', 'csc', 'dok', 'lil', 'coo', 'bsr'],
        ['C', 'F'],
    ]
    param_names = ['format', 'order']

    def setup(self, format, order):
        warnings.simplefilter('ignore', SparseEfficiencyWarning)
        self.X = sparse.rand(1000, 1000, format=format, density=0.01)

    def time_toarray(self, format, order):
        self.X.toarray(order=order)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_toarray.version = (
        "2fbf492ec800b982946a62785beda803460b913cc80080043a5d407025893b2b"
    )


class Random(Benchmark):
    params = [
        np.arange(0, 1.1, 0.1).tolist()
    ]
    param_names = ['density']

    def setup(self, density):
        warnings.simplefilter('ignore', SparseEfficiencyWarning)
        self.nrows = 1000
        self.ncols = 1000
        self.format = 'csr'

    def time_rand(self, density):
        sparse.rand(self.nrows, self.ncols,
                    format=self.format, density=density)
