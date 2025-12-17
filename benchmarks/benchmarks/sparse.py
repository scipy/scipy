"""
Simple benchmarks for the sparse module
"""
import warnings
import time
import timeit
import pickle

import numpy as np

from .common import Benchmark, safe_import

with safe_import():
    from scipy import sparse


def random_sparse(m, n, nnz_per_row, sparse_type):
    rows = np.arange(m).repeat(nnz_per_row)
    cols = np.random.randint(0, n, size=rows.size)
    vals = np.random.random_sample(rows.size)
    coo = sparse.coo_array if sparse_type == "sparray" else sparse.coo_matrix
    return coo((vals, (rows, cols)), (m, n)).tocsr()


# TODO move this to a matrix gallery and add unittests
def poisson2d(N, dtype='d', format=None, sparse_type="sparray"):
    """
    Return a sparse matrix for the 2D Poisson problem
    with standard 5-point finite difference stencil on a
    square N-by-N grid.
    """
    dia = sparse.dia_array if sparse_type == "sparray" else sparse.dia_matrix
    if N == 1:
        diags = np.asarray([[4]], dtype=dtype)
        return dia((diags, [0]), shape=(1, 1)).asformat(format)

    offsets = np.array([0, -N, N, -1, 1])

    diags = np.empty((5, N**2), dtype=dtype)

    diags[0] = 4  # main diagonal
    diags[1:] = -1  # all offdiagonals

    diags[3, N-1::N] = 0  # first lower diagonal
    diags[4, N::N] = 0  # first upper diagonal

    return dia((diags, offsets), shape=(N**2, N**2)).asformat(format)


class Arithmetic(Benchmark):
    param_names = ['sparse_type', 'format', 'XY', 'op']
    params = [
        ['spmatrix', 'sparray'],
        ['csr', 'csc', 'coo', 'dia'],
        ['AA', 'AB', 'BA', 'BB'],
        ['__add__', '__sub__', 'multiply', '__mul__']
    ]

    def setup(self, sparse_type, format, XY, op):
        matrices = dict(A=poisson2d(250, format=format, sparse_type=sparse_type))
        matrices['B'] = (matrices['A']**2).asformat(format)

        x = matrices[XY[0]]
        self.y = matrices[XY[1]]
        self.fn = getattr(x, op)
        self.fn(self.y)  # warmup

    def time_arithmetic(self, sparse_type, format, XY, op):
        self.fn(self.y)


class Sort(Benchmark):
    param_names = ['sparse_type', 'name',]
    params = [
        ['spmatrix', 'sparray'],
        ['Rand10', 'Rand25', 'Rand50', 'Rand100', 'Rand200'],
    ]

    def setup(self, sparse_type, name):
        n = 10000
        if name.startswith('Rand'):
            k = int(name[4:])
            self.A = random_sparse(n, n, k, sparse_type)
            self.A.has_sorted_indices = False
            self.A.indices[:2] = 2, 1
        else:
            raise NotImplementedError()

    def time_sort(self, sparse_type, name):
        """sort CSR column indices"""
        self.A.sort_indices()


class Matvec(Benchmark):
    param_names = ['sparse_type', 'name', 'format']
    params = [
        ['spmatrix', 'sparray'],
        ['Identity', 'Poisson5pt', 'Block2x2', 'Block3x3'],
        ['dia', 'csr', 'csc', 'dok', 'lil', 'coo', 'bsr'],
    ]

    def setup(self, sparse_type, name, format):
        if name == 'Identity':
            if format in ('lil', 'dok'):
                raise NotImplementedError()
            if sparse_type == "sparray":
                self.A = sparse.eye_array(10000, format=format)
            else:
                self.A = sparse.eye(10000, format=format)
        elif name == 'Poisson5pt':
            self.A = poisson2d(300, format=format, sparse_type=sparse_type)
        elif name == 'Block2x2':
            if format not in ('csr', 'bsr'):
                raise NotImplementedError()
            b = (2, 2)
            self.A = sparse.kron(poisson2d(150, sparse_type=sparse_type),
                                 np.ones(b)).tobsr(blocksize=b).asformat(format)
        elif name == 'Block3x3':
            if format not in ('csr', 'bsr'):
                raise NotImplementedError()
            b = (3, 3)
            self.A = sparse.kron(poisson2d(100, sparse_type=sparse_type),
                                 np.ones(b)).tobsr(blocksize=b).asformat(format)
        else:
            raise NotImplementedError()

        self.x = np.ones(self.A.shape[1], dtype=float)

    def time_matvec(self, sparse_type, name, format):
        self.A @ self.x


class Matvecs(Benchmark):
    param_names = ['sparse_type', "format"]
    params = [
        ['spmatrix', 'sparray'],
        ['dia', 'coo', 'csr', 'csc', 'bsr'],
    ]

    def setup(self, sparse_type, format):
        self.A = poisson2d(300, format=format, sparse_type=sparse_type)
        self.x = np.ones((self.A.shape[1], 10), dtype=self.A.dtype)

    def time_matvecs(self, sparse_type, format):
        self.A @ self.x


class Matmul(Benchmark):
    param_names = ['sparse_type']
    params = [
        ['spmatrix', 'sparray']
    ]

    def setup(self, sparse_type):
        coo = sparse.coo_array if sparse_type == "sparray" else sparse.coo_matrix
        H1, W1 = 1, 100000
        H2, W2 = W1, 1000
        C1 = 10
        C2 = 1000000

        rng = np.random.default_rng(0)

        i = rng.integers(H1, size=C1)
        j = rng.integers(W1, size=C1)
        data = rng.random(C1)
        self.matrix1 = coo((data, (i, j)), shape=(H1, W1)).tocsr()

        i = rng.integers(H2, size=C2)
        j = rng.integers(W2, size=C2)
        data = rng.random(C2)
        self.matrix2 = coo((data, (i, j)), shape=(H2, W2)).tocsr()

    def time_large(self, sparse_type):
        for i in range(100):
            self.matrix1 @ self.matrix2

    # Retain old benchmark results (remove this if changing the benchmark)
    time_large.version = (
        "33aee08539377a7cb0fabaf0d9ff9d6d80079a428873f451b378c39f6ead48cb"
    )


class Construction(Benchmark):
    param_names = ['sparse_type', 'name', 'format']
    params = [
        ['spmatrix', 'sparray'],
        ['Empty', 'Identity', 'Poisson5pt'],
        ['lil', 'dok'],
    ]

    def setup(self, sparse_type, name, format):
        if name == 'Empty':
            self.A = sparse.coo_array((10000, 10000))
        elif name == 'Identity':
            self.A = sparse.eye_array(10000, format='coo')
        else:
            self.A = poisson2d(100, format='coo')
        if sparse_type == "spmatrix":
            self.A = sparse.coo_matrix(self.A)

        formats = {'lil': sparse.lil_matrix, 'dok': sparse.dok_matrix}
        self.cls = formats[format]

    def time_construction(self, sparse_type, name, format):
        T = self.cls(self.A.shape)
        for v, i, j in zip(self.A.data, *self.A.coords):
            T[i, j] = v


class BlockDiagDenseConstruction(Benchmark):
    param_names = ['sparse_type', 'num_matrices']
    params = [
        ['spmatrix', 'sparray'],
        [1000, 5000, 10000, 15000, 20000],
    ]

    def setup(self, sparse_type, num_matrices):
        coo = sparse.coo_array if sparse_type == "sparray" else sparse.coo_matrix
        self.matrices = []
        for i in range(num_matrices):
            rows = np.random.randint(1, 4)
            columns = np.random.randint(1, 4)
            mat = np.random.randint(0, 10, (rows, columns))
            if i == 0:
                self.matrices.append(coo(mat))  # make 1st requested sparse_type
            else:
                self.matrices.append(mat)

    def time_block_diag(self, sparse_type, num_matrices):
        sparse.block_diag(self.matrices)


class BlockDiagSparseConstruction(Benchmark):
    param_names = ['sparse_type', 'num_matrices']
    params = [
        ['spmatrix', 'sparray'],
        [1000, 5000, 10000, 15000, 20000],
    ]

    def setup(self, sparse_type, num_matrices):
        self.matrices = []
        for i in range(num_matrices):
            rows = np.random.randint(1, 20)
            columns = np.random.randint(1, 20)
            density = 2e-3
            nnz_per_row = int(density*columns)

            mat = random_sparse(rows, columns, nnz_per_row, sparse_type)
            self.matrices.append(mat)

    def time_block_diag(self, sparse_type, num_matrices):
        sparse.block_diag(self.matrices)


class CsrHstack(Benchmark):
    param_names = ['sparse_type', 'num_rows']
    params = [
        ['spmatrix', 'sparray'],
        [10000, 25000, 50000, 100000, 250000],
    ]

    def setup(self, sparse_type, num_rows):
        num_cols = int(1e5)
        density = 2e-3
        nnz_per_row = int(density*num_cols)
        self.mat = random_sparse(num_rows, num_cols, nnz_per_row, sparse_type)

    def time_csr_hstack(self, sparse_type, num_rows):
        sparse.hstack([self.mat, self.mat])


class Conversion(Benchmark):
    param_names = ['sparse_type', 'from_format', 'to_format']
    params = [
        ['spmatrix', 'sparray'],
        ['csr', 'csc', 'coo', 'dia', 'lil', 'dok', 'bsr'],
        ['csr', 'csc', 'coo', 'dia', 'lil', 'dok', 'bsr'],
    ]

    def setup(self, sparse_type, from_format, to_format):
        base = poisson2d(100, format=from_format, sparse_type=sparse_type)

        try:
            self.fn = getattr(base, 'to' + to_format)
        except Exception:
            def fn():
                raise RuntimeError()
            self.fn = fn

    def time_conversion(self, sparse_type, from_format, to_format):
        self.fn()


class Getset(Benchmark):
    param_names = ['sparse_type', 'N', 'sparsity pattern', 'format']
    params = [
        ['spmatrix', 'sparray'],
        [1, 10, 100, 1000, 10000],
        ['different', 'same'],
        ['csr', 'csc', 'lil', 'dok'],
    ]
    unit = "seconds"

    def setup(self, sparse_type, N, sparsity_pattern, format):
        if format == 'dok' and N > 500:
            raise NotImplementedError()

        if sparse_type == "sparray":
            A = self.A = sparse.random_array((1000, 1000), density=1e-5)
        else:
            A = self.A = sparse.random(1000, 1000, density=1e-5)

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

    def track_fancy_setitem(self, sparse_type, N, sparsity_pattern, format):
        def kernel(A, i, j, v):
            A[i, j] = v

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', sparse.SparseEfficiencyWarning)
            return self._timeit(kernel, sparsity_pattern == 'different')

    def time_fancy_getitem(self, sparse_type, N, sparsity_pattern, format):
        self.m[self.i, self.j]


class NullSlice(Benchmark):
    param_names = ['sparse_type', 'density', 'format']
    params = [
        ['spmatrix', 'sparray'],
        [0.05, 0.01],
        ['csr', 'csc', 'lil'],
    ]

    def _setup(self, sparse_type, density, format):
        n = 100000
        k = 1000

        # faster version of sparse.rand(n, k, format=format, density=density),
        # with non-exact nnz
        nz = int(n*k * density)
        row = np.random.randint(0, n, size=nz)
        col = np.random.randint(0, k, size=nz)
        data = np.ones(nz, dtype=np.float64)
        coo = sparse.coo_array if sparse_type == "sparray" else sparse.coo_matrix
        X = coo((data, (row, col)), shape=(n, k))
        X.sum_duplicates()
        X = X.asformat(format)
        with open(f'{sparse_type}-{density}-{format}.pck', 'wb') as f:
            pickle.dump(X, f, protocol=pickle.HIGHEST_PROTOCOL)

    def setup_cache(self):
        for sparse_type in self.params[0]:
            for density in self.params[1]:
                for fmt in self.params[2]:
                    self._setup(sparse_type, density, fmt)

    setup_cache.timeout = 120

    def setup(self, sparse_type, density, format):
        # Unpickling is faster than computing the random matrix...
        with open(f'{sparse_type}-{density}-{format}.pck', 'rb') as f:
            self.X = pickle.load(f)

    def time_getrow(self, sparse_type, density, format):
        if sparse_type == "sparray":
            self.X[100]
        else:
            self.X.getrow(100)

    def time_getcol(self, sparse_type, density, format):
        if sparse_type == "sparray":
            self.X[:, 100]
        else:
            self.X.getcol(100)

    def time_3_rows(self, sparse_type, density, format):
        self.X[[0, 100, 105], :]

    def time_10000_rows(self, sparse_type, density, format):
        self.X[np.arange(10000), :]

    def time_3_cols(self, sparse_type, density, format):
        self.X[:, [0, 100, 105]]

    def time_100_cols(self, sparse_type, density, format):
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
    param_names = ['sparse_type', 'density', 'format']
    params = [
        ['spmatrix', 'sparray'],
        [0.01, 0.1, 0.5],
        ['csr', 'csc', 'coo', 'lil', 'dok', 'dia'],
    ]

    def setup(self, sparse_type, density, format):
        n = 1000
        if format == 'dok' and n * density >= 500:
            raise NotImplementedError()

        warnings.simplefilter('ignore', sparse.SparseEfficiencyWarning)

        if sparse_type == "sparray":
            self.X = sparse.random_array((n, n), format=format, density=density)
        else:
            self.X = sparse.random(n, n, format=format, density=density)

    def time_diagonal(self, sparse_type, density, format):
        self.X.diagonal()

    # Retain old benchmark results (remove this if changing the benchmark)
    time_diagonal.version = (
        "d84f53fdc6abc208136c8ce48ca156370f6803562f6908eb6bd1424f50310cf1"
    )


class Sum(Benchmark):
    param_names = ['sparse_type', 'density', 'format']
    params = [
        ['spmatrix', 'sparray'],
        [0.01, 0.1, 0.5],
        ['csr', 'csc', 'coo', 'lil', 'dok', 'dia'],
    ]

    def setup(self, sparse_type, density, format):
        n = 1000
        if format == 'dok' and n * density >= 500:
            raise NotImplementedError()

        warnings.simplefilter('ignore', sparse.SparseEfficiencyWarning)
        if sparse_type == "sparray":
            self.X = sparse.random_array((n, n), format=format, density=density)
        else:
            self.X = sparse.random(n, n, format=format, density=density)

    def time_sum(self, sparse_type, density, format):
        self.X.sum()

    def time_sum_axis0(self, sparse_type, density, format):
        self.X.sum(axis=0)

    def time_sum_axis1(self, sparse_type, density, format):
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
    param_names = ['sparse_type', 'density', 'format']
    params = [
        ['spmatrix', 'sparray'],
        [0.05, 0.01],
        ['csr', 'csc', 'lil'],
    ]

    def setup(self, sparse_type, density, format):
        n = 500
        k = 1000
        if sparse_type == "sparray":
            self.X = sparse.random_array((n, k), format=format, density=density)
        else:
            self.X = sparse.random(n, k, format=format, density=density)

    def time_iteration(self, sparse_type, density, format):
        for row in self.X:
            pass


class Densify(Benchmark):
    param_names = ['sparse_type', 'format', 'order']
    params = [
        ['spmatrix', 'sparray'],
        ['dia', 'csr', 'csc', 'dok', 'lil', 'coo', 'bsr'],
        ['C', 'F'],
    ]

    def setup(self, sparse_type, format, order):
        warnings.simplefilter('ignore', sparse.SparseEfficiencyWarning)
        if sparse_type == "sparray":
            self.X = sparse.random_array((1000, 1000), format=format, density=0.01)
        else:
            self.X = sparse.random(1000, 1000, format=format, density=0.01)

    def time_toarray(self, sparse_type, format, order):
        self.X.toarray(order=order)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_toarray.version = (
        "2fbf492ec800b982946a62785beda803460b913cc80080043a5d407025893b2b"
    )


class Random(Benchmark):
    param_names = ['sparse_type', 'density']
    params = [
        ['spmatrix', 'sparray'],
        np.arange(0, 1.1, 0.1).tolist(),
    ]

    def setup(self, sparse_type, density):
        warnings.simplefilter('ignore', sparse.SparseEfficiencyWarning)
        self.nrows = 1000
        self.ncols = 1000
        self.format = 'csr'

    def time_rand(self, sparse_type, density):
        if sparse_type == "sparray":
            self.X = sparse.random_array(
                (self.nrows, self.ncols), format=self.format, density=density
            )
        else:
            self.X = sparse.random(
                self.nrows, self.ncols, format=self.format, density=density
            )


class Argmax(Benchmark):
    param_names = ['sparse_type', 'density', 'format', 'explicit']
    params = [
        ['spmatrix', 'sparray'],
        [0.01, 0.1, 0.5],
        ['csr', 'csc', 'coo'],
        [True, False],
    ]

    def setup(self, sparse_type, density, format, explicit):
        n = 1000

        warnings.simplefilter('ignore', sparse.SparseEfficiencyWarning)
        if sparse_type == "sparray":
            self.X = sparse.random_array((n, n), format=format, density=density)
        else:
            self.X = sparse.random(n, n, format=format, density=density)

    def time_argmax(self, sparse_type, density, format, explicit):
        self.X.argmax(explicit=explicit)
