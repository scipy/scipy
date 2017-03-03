"""
Simple benchmarks for the sparse module
"""
from __future__ import division, print_function, absolute_import

import warnings
import time
import timeit

import numpy
import numpy as np
from numpy import ones, array, asarray, empty, random, zeros
from numpy.testing import assert_

try:
    from scipy import sparse
    from scipy.sparse import (csr_matrix, coo_matrix, dia_matrix, lil_matrix,
                              dok_matrix, rand, SparseEfficiencyWarning)
except ImportError:
    pass

from .common import Benchmark


def random_sparse(m, n, nnz_per_row):
    rows = numpy.arange(m).repeat(nnz_per_row)
    cols = numpy.random.randint(0, n, size=nnz_per_row*m)
    vals = numpy.random.random_sample(m*nnz_per_row)
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

def has_inplace_matvec():
    A = sparse.eye(3, format='csr')
    x = ones(3)
    y = ones(3)
    try:
        sparse._sparsetools.csr_matvec(3, 3, A.indptr, A.indices, A.data, 2, x, 3, y)
        return True
    except:
        return False

if has_inplace_matvec():
    def sparse_gemm(a, b, c, alpha, beta):
        fmt = a.format
        M, N = a.shape
        n_vecs = 1
        if b.ndim == 2:
            n_vecs = b.shape[1]

        s = sparse._sparsetools
        if n_vecs == 1:
            if fmt == 'csr':
                s.csr_matvec(M, N, a.indptr, a.indices, a.data, alpha, b.ravel(), beta, c.ravel())
            elif fmt == 'csc':
                s.csc_matvec(M, N, a.indptr, a.indices, a.data, alpha, b.ravel(), beta, c.ravel())
            elif fmt == 'bsr':
                R, C = a.blocksize
                s.bsr_matvec(M//R, N//R, R, C, a.indptr, a.indices, a.data, alpha, b.ravel(),
                             beta, c.ravel())
            elif fmt == 'dia':
                L = a.data.shape[1]
                s.dia_matvec(M,N, len(a.offsets), L, a.offsets, a.data, alpha, b.ravel(), beta, c.ravel())
            elif fmt == 'coo':
                s.coo_matvec(M, a.nnz, a.row, a.col, a.data, alpha, b.ravel(), beta, c.ravel())
            else:
                q = alpha * a*b
                c *= beta
                c += q
        else:
            if fmt == 'csr':
                s.csr_matvecs(M, N, n_vecs, a.indptr, a.indices, a.data, alpha, b.ravel(), beta, c.ravel())
            elif fmt == 'csc':
                s.csc_matvecs(M, N, n_vecs, a.indptr, a.indices, a.data, alpha, b.ravel(), beta, c.ravel())
            elif fmt == 'bsr':
                R, C = a.blocksize
                s.bsr_matvecs(M//R, N//R, n_vecs, R, C, a.indptr, a.indices, a.data, alpha, b.ravel(),
                             beta, c.ravel())
            else:
                raise ValueError('unsupported type')
else:
    def sparse_gemm(a, b, c, alpha, beta):
        q = a * b
        q *= alpha
        c *= beta
        c += q


class Arithmetic(Benchmark):
    param_names = ['format', 'XY', 'op']
    params = [
        ['csr'],
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
        self.y = ones(self.A.shape[0], dtype=float)

    def time_matvec(self, matrix, format):
        self.A * self.x

    def time_matvec_inplace_1_0(self, matrix, format):
        sparse_gemm(self.A, self.x, self.y, 1, 0)

    def time_matvec_inplace_1_1(self, matrix, format):
        sparse_gemm(self.A, self.x, self.y, 1, 1)

    def time_matvec_inplace_2_3(self, matrix, format):
        sparse_gemm(self.A, self.x, self.y, 2, 3)


class Matvecs(Benchmark):
    params = ['dia', 'coo', 'csr', 'csc', 'bsr']
    param_names = ["format"]

    def setup(self, format):
        self.A = poisson2d(300, format=format)
        self.x = ones((self.A.shape[1], 10), dtype=self.A.dtype)
        self.y = ones((self.A.shape[0], 10), dtype=self.A.dtype)

    def time_matvecs(self, format):
        self.A * self.x

    def time_matvecs_inplace_1_0(self, fmt):
        A = self.matrices[fmt]
        sparse_gemm(self.A, self.x, self.y, 1, 0)

    def time_matvecs_inplace_1_1(self, fmt):
        A = self.matrices[fmt]
        sparse_gemm(self.A, self.x, self.y, 1, 1)

    def time_matvecs_inplace_2_3(self, fmt):
        A = self.matrices[fmt]
        sparse_gemm(self.A, self.x, self.y, 2, 3)


class Matmul(Benchmark):
    def setup(self):
        H1, W1 = 1, 100000
        H2, W2 = W1, 1000
        C1 = 10
        C2 = 1000000

        random.seed(0)

        matrix1 = lil_matrix(zeros((H1, W1)))
        matrix2 = lil_matrix(zeros((H2, W2)))
        for i in range(C1):
            matrix1[random.randint(H1), random.randint(W1)] = random.rand()
        for i in range(C2):
            matrix2[random.randint(H2), random.randint(W2)] = random.rand()
        self.matrix1 = matrix1.tocsr()
        self.matrix2 = matrix2.tocsr()

    def time_large(self):
        for i in range(100):
            self.matrix1 * self.matrix2


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


class Conversion(Benchmark):
    params = [
        ['csr', 'csc', 'coo', 'dia', 'lil', 'dok'],
        ['csr', 'csc', 'coo', 'dia', 'lil', 'dok'],
    ]
    param_names = ['from_format', 'to_format']

    def setup(self, fromfmt, tofmt):
        base = poisson2d(100, format=fromfmt)

        try:
            self.fn = getattr(base, 'to' + tofmt)
        except:
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
            ip = numpy.random.randint(0, A.shape[0], size=n)
            jp = numpy.random.randint(0, A.shape[1], size=n)
            i = numpy.r_[i, ip]
            j = numpy.r_[j, jp]
        v = numpy.random.rand(n)

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

    def setup(self, density, format):
        n = 100000
        k = 1000
        self.X = sparse.rand(n, k, format=format, density=density)

    def time_3_rows(self, density, format):
        self.X[[0, 100, 105], :]

    def time_10000_rows(self, density, format):
        self.X[np.arange(10000), :]

    def time_3_cols(self, density, format):
        self.X[:, [0, 100, 105]]

    def time_100_cols(self, density, format):
        self.X[:, np.arange(100)]


class Diagonal(Benchmark):
    params = [[0.01, 0.1, 0.5], ['csr', 'csc', 'coo', 'lil', 'dok', 'dia']]
    param_names = ['density', 'format']

    def setup(self, density, format):
        n = 1000
        if format == 'dok' and n * density >= 500:
            raise NotImplementedError()

        self.X = sparse.rand(n, n, format=format, density=density)

    def time_diagonal(self, density, format):
        self.X.diagonal()


class Sum(Benchmark):
    params = [[0.01, 0.1, 0.5], ['csr', 'csc', 'coo', 'lil', 'dok', 'dia']]
    param_names = ['density', 'format']

    def setup(self, density, format):
        n = 1000
        if format == 'dok' and n * density >= 500:
            raise NotImplementedError()

        self.X = sparse.rand(n, n, format=format, density=density)

    def time_sum(self, density, format):
        self.X.sum()

    def time_sum_axis0(self, density, format):
        self.X.sum(axis=0)

    def time_sum_axis1(self, density, format):
        self.X.sum(axis=1)
