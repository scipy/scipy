"""
Simple benchmarks for the sparse module
"""
from __future__ import division, print_function, absolute_import

import warnings
import time
import collections
import sys
import timeit

import numpy
import numpy as np
from numpy import ones, array, asarray, empty, random, zeros

try:
    import scipy
    from scipy import sparse
    from scipy.sparse import (csr_matrix, coo_matrix, dia_matrix, lil_matrix,
                              dok_matrix, rand, SparseEfficiencyWarning)
except ImportError:
    pass

from .common import Benchmark


def random_sparse(m,n,nnz_per_row):
    rows = numpy.arange(m).repeat(nnz_per_row)
    cols = numpy.random.random_integers(low=0,high=n-1,size=nnz_per_row*m)
    vals = numpy.random.random_sample(m*nnz_per_row)
    return coo_matrix((vals,(rows,cols)),(m,n)).tocsr()


# TODO move this to a matrix gallery and add unittests
def poisson2d(N,dtype='d',format=None):
    """
    Return a sparse matrix for the 2D Poisson problem
    with standard 5-point finite difference stencil on a
    square N-by-N grid.
    """
    if N == 1:
        diags = asarray([[4]],dtype=dtype)
        return dia_matrix((diags,[0]), shape=(1,1)).asformat(format)

    offsets = array([0,-N,N,-1,1])

    diags = empty((5,N**2),dtype=dtype)

    diags[0] = 4  # main diagonal
    diags[1:] = -1  # all offdiagonals

    diags[3,N-1::N] = 0  # first lower diagonal
    diags[4,N::N] = 0  # first upper diagonal

    return dia_matrix((diags,offsets),shape=(N**2,N**2)).asformat(format)


class Arithmetic(Benchmark):
    param_names = ['format', 'XY', 'op']
    params = [
        ['csr'],
        ['AA', 'AB', 'BA', 'BB'],
        ['__add__', '__sub__', 'multiply', '__mul__']
    ]

    def setup(self, format, XY, op):
        self.matrices = {}
        # matrices.append( ('A','Identity', sparse.eye(500**2,format='csr')) )
        self.matrices['A'] = poisson2d(250,format='csr')
        self.matrices['B'] = poisson2d(250,format='csr')**2

        X, Y = XY
        vars = dict([(var, mat.asformat(format)) 
                     for (var, mat) in self.matrices.items()])
        self.x, self.y = vars[X], vars[Y]
        self.fn = getattr(self.x, op)
        self.fn(self.y)  # warmup

    def time_arithmetic(self, format, XY, op):
        self.fn(self.y)


class Sort(Benchmark):
    params = ['Rand10', 'Rand25', 'Rand50', 'Rand100', 'Rand200']
    param_names = ['matrix']

    def setup(self, matrix):
        matrices = []
        matrices.append(('Rand10', (1e4, 10)))
        matrices.append(('Rand25', (1e4, 25)))
        matrices.append(('Rand50', (1e4, 50)))
        matrices.append(('Rand100', (1e4, 100)))
        matrices.append(('Rand200', (1e4, 200)))
        self.matrices = dict(matrices)

        N, K = self.matrices[matrix]
        N = int(float(N))
        K = int(float(K))
        self.A = random_sparse(N,N,K)

    def time_sort(self, matrix):
        """sort CSR column indices"""
        self.A.has_sorted_indices = False
        self.A.indices[:2] = 2,1
        self.A.sort_indices()


class Matvec(Benchmark):
    param_names = ['matrix']

    @property
    def params(self):
        return list(sorted(self._get_matrices().keys()))

    def _get_matrices(self):
        matrices = collections.OrderedDict()

        matrices['Identity_dia'] = sparse.eye(10**4,10**4,format='dia')
        matrices['Identity_csr'] = sparse.eye(10**4,10**4,format='csr')
        matrices['Poisson5pt_lil'] = poisson2d(300,format='lil')
        matrices['Poisson5pt_dok'] = poisson2d(300,format='dok')
        matrices['Poisson5pt_dia'] = poisson2d(300,format='dia')
        matrices['Poisson5pt_coo'] = poisson2d(300,format='coo')
        matrices['Poisson5pt_csr'] = poisson2d(300,format='csr')
        matrices['Poisson5pt_csc'] = poisson2d(300,format='csc')
        matrices['Poisson5pt_bsr'] = poisson2d(300,format='bsr')

        A = sparse.kron(poisson2d(150),ones((2,2))).tobsr(blocksize=(2,2))
        matrices['Block2x2_csr'] = A.tocsr()
        matrices['Block2x2_bsr'] = A

        A = sparse.kron(poisson2d(100),ones((3,3))).tobsr(blocksize=(3,3))
        matrices['Block3x3_csr'] = A.tocsr()
        matrices['Block3x3_bsr'] = A
        return matrices

    def setup(self, matrix):
        self.matrices = self._get_matrices()
        self.x = ones(max(A.shape[1] for A in self.matrices.values()), 
                      dtype=float)

        self.A = self.matrices[matrix]
        self.x = ones(self.A.shape[1], dtype=float)

    def time_matvec(self, matrix):
        self.A * self.x


class Matvecs(Benchmark):
    params = ['dia', 'coo', 'csr', 'csc', 'bsr']
    param_names = ["format"]

    def setup(self, *args):
        self.matrices = {}
        self.matrices['dia'] = poisson2d(300,format='dia')
        self.matrices['coo'] = poisson2d(300,format='coo')
        self.matrices['csr'] = poisson2d(300,format='csr')
        self.matrices['csc'] = poisson2d(300,format='csc')
        self.matrices['bsr'] = poisson2d(300,format='bsr')
        A = self.matrices['dia']
        self.x = ones((A.shape[1], 10), dtype=A.dtype)

    def time_matvecs(self, fmt):
        A = self.matrices[fmt]
        y = A*self.x


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
            matrix3 = self.matrix1 * self.matrix2


class Construction(Benchmark):
    params = [
        ['Empty', 'Identity', 'Poisson5pt'],
        ['lil', 'dok']
    ]
    param_names = ['matrix', 'format']

    def setup(self, name, format):
        self.matrices = {}
        self.matrices['Empty'] = csr_matrix((10000,10000))
        self.matrices['Identity'] = sparse.eye(10000)
        self.matrices['Poisson5pt'] = poisson2d(100)
        self.formats = {'lil': lil_matrix, 'dok': dok_matrix}

        A = self.matrices[name]
        self.cls = self.formats[format]
        self.A = A.tocoo()

    def time_construction(self, name, format):
        T = self.cls(self.A.shape)
        for i,j,v in zip(self.A.row,self.A.col,self.A.data):
            T[i,j] = v


class Conversion(Benchmark):
    params = [
        ['csr','csc','coo','dia','lil','dok'],
        ['csr','csc','coo','dia','lil','dok'],
    ]
    param_names = ['from_format', 'to_format']

    def setup(self, fromfmt, tofmt):
        self.A = poisson2d(100)

        A = self.A
        base = getattr(A,'to' + fromfmt)()

        result = np.nan
        try:
            self.fn = getattr(base, 'to' + tofmt)
        except:
            def fn():
                raise RuntimeError()
            self.fn = fn

    def time_conversion(self, fromfmt, tofmt):
        x = self.fn()


class Getset(Benchmark):
    params = [
        [1, 10, 100, 1000, 10000],
        ['different', 'same'],
        ['csr', 'csc', 'lil', 'dok']
    ]
    param_names = ['N', 'sparsity pattern', 'format']

    def setup(self, N, sparsity_pattern, format):
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
        while time.time() - start < 0.5:
            if recopy:
                m = self.m.copy()
            else:
                m = self.m
            while True:
                duration = timeit.timeit(lambda: kernel(m, self.i, self.j, self.v),
                                         number=number)
                if duration > 1e-5:
                    break
                else:
                    number *= 10
            min_time = min(min_time, duration/number)
        return min_time

    def track_fancy_setitem(self, N, sparsity_pattern, format):
        if format == 'dok' and N > 500:
            return np.nan

        def kernel(A, i, j, v):
            A[i, j] = v

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SparseEfficiencyWarning)
            return self._timeit(kernel, sparsity_pattern=='different')

    def track_fancy_getitem(self, N, sparsity_pattern, format):
        if format == 'dok' and N > 500:
            return np.nan

        def kernel(A, i, j, v):
            A[i, j]

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SparseEfficiencyWarning)
            return self._timeit(kernel, sparsity_pattern=='different')
