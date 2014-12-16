"""general tests and simple benchmarks for the sparse module"""
from __future__ import division, print_function, absolute_import

import time
import warnings

import numpy
from numpy import ones, array, asarray, empty, random, zeros

from numpy.testing import Tester, TestCase

import scipy
from scipy import sparse
from scipy.sparse import (csr_matrix, coo_matrix, dia_matrix, lil_matrix,
                          dok_matrix, rand, SparseEfficiencyWarning)


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


class BenchmarkSparse(TestCase):
    """Simple benchmarks for sparse matrix module"""

    def bench_arithmetic(self):
        matrices = []
        # matrices.append( ('A','Identity', sparse.eye(500**2,format='csr')) )
        matrices.append(('A','Poisson5pt', poisson2d(250,format='csr')))
        matrices.append(('B','Poisson5pt^2', poisson2d(250,format='csr')**2))

        print()
        print('                 Sparse Matrix Arithmetic')
        print('====================================================================')
        print(' var |     name       |         shape        |   dtype   |    nnz   ')
        print('--------------------------------------------------------------------')
        fmt = '  %1s  | %14s | %20s | %9s | %8d '

        for var,name,mat in matrices:
            name = name.center(14)
            shape = ("%s" % (mat.shape,)).center(20)
            dtype = mat.dtype.name.center(9)
            print(fmt % (var,name,shape,dtype,mat.nnz))

        space = ' ' * 10
        print()
        print(space+'              Timings')
        print(space+'==========================================')
        print(space+' format |     operation     | time (msec) ')
        print(space+'------------------------------------------')
        fmt = space+'   %3s  | %17s |  %7.1f  '

        for format in ['csr']:
            vars = dict([(var, mat.asformat(format)) for (var, _, mat) in matrices])
            for X,Y in [('A','A'),('A','B'),('B','A'),('B','B')]:
                x,y = vars[X],vars[Y]
                for op in ['__add__','__sub__','multiply','__mul__']:
                    fn = getattr(x,op)
                    fn(y)  # warmup

                    start = time.clock()
                    iter = 0
                    while iter < 3 or time.clock() < start + 0.5:
                        fn(y)
                        iter += 1
                    end = time.clock()

                    msec_per_it = 1000*(end - start)/float(iter)
                    operation = (X + '.' + op + '(' + Y + ')').center(17)
                    print(fmt % (format,operation,msec_per_it))

    def bench_sort(self):
        """sort CSR column indices"""
        matrices = []
        matrices.append(('Rand10', 1e4, 10))
        matrices.append(('Rand25', 1e4, 25))
        matrices.append(('Rand50', 1e4, 50))
        matrices.append(('Rand100', 1e4, 100))
        matrices.append(('Rand200', 1e4, 200))

        print()
        print('                    Sparse Matrix Index Sorting')
        print('=====================================================================')
        print(' type |    name      |         shape        |    nnz   | time (msec) ')
        print('---------------------------------------------------------------------')
        fmt = '  %3s | %12s | %20s | %8d |   %6.2f  '

        for name,N,K in matrices:
            N = int(N)
            A = random_sparse(N,N,K)

            start = time.clock()
            iter = 0
            while iter < 5 and time.clock() - start < 1:
                A.has_sorted_indices = False
                A.indices[:2] = 2,1
                A.sort_indices()
                iter += 1
            end = time.clock()

            name = name.center(12)
            shape = ("%s" % (A.shape,)).center(20)

            print(fmt % (A.format,name,shape,A.nnz,1e3*(end-start)/float(iter)))

    def bench_matvec(self):
        matrices = []
        matrices.append(('Identity', sparse.eye(10**4,format='dia')))
        matrices.append(('Identity', sparse.eye(10**4,format='csr')))
        matrices.append(('Poisson5pt', poisson2d(300,format='lil')))
        matrices.append(('Poisson5pt', poisson2d(300,format='dok')))
        matrices.append(('Poisson5pt', poisson2d(300,format='dia')))
        matrices.append(('Poisson5pt', poisson2d(300,format='coo')))
        matrices.append(('Poisson5pt', poisson2d(300,format='csr')))
        matrices.append(('Poisson5pt', poisson2d(300,format='csc')))
        matrices.append(('Poisson5pt', poisson2d(300,format='bsr')))

        A = sparse.kron(poisson2d(150),ones((2,2))).tobsr(blocksize=(2,2))
        matrices.append(('Block2x2', A.tocsr()))
        matrices.append(('Block2x2', A))

        A = sparse.kron(poisson2d(100),ones((3,3))).tobsr(blocksize=(3,3))
        matrices.append(('Block3x3', A.tocsr()))
        matrices.append(('Block3x3', A))

        print()
        print('                 Sparse Matrix Vector Product')
        print('==================================================================')
        print(' type |    name      |         shape        |    nnz   |  MFLOPs  ')
        print('------------------------------------------------------------------')
        fmt = '  %3s | %12s | %20s | %8d |  %6.1f '

        for name,A in matrices:
            x = ones(A.shape[1],dtype=A.dtype)

            y = A*x  # warmup

            start = time.clock()
            iter = 0
            while iter < 5 or time.clock() < start + 1:
                y = A*x
                iter += 1
            end = time.clock()

            del y

            name = name.center(12)
            shape = ("%s" % (A.shape,)).center(20)
            MFLOPs = (2*A.nnz*iter/(end-start))/float(1e6)

            print(fmt % (A.format,name,shape,A.nnz,MFLOPs))

    def bench_matvecs(self):
        matrices = []
        matrices.append(('Poisson5pt', poisson2d(300,format='dia')))
        matrices.append(('Poisson5pt', poisson2d(300,format='coo')))
        matrices.append(('Poisson5pt', poisson2d(300,format='csr')))
        matrices.append(('Poisson5pt', poisson2d(300,format='csc')))
        matrices.append(('Poisson5pt', poisson2d(300,format='bsr')))

        n_vecs = 10

        print()
        print('             Sparse Matrix (Block) Vector Product')
        print('                       Blocksize = %d' % (n_vecs,))
        print('==================================================================')
        print(' type |    name      |         shape        |    nnz   |  MFLOPs  ')
        print('------------------------------------------------------------------')
        fmt = '  %3s | %12s | %20s | %8d |  %6.1f '

        for name,A in matrices:
            x = ones((A.shape[1],10),dtype=A.dtype)

            y = A*x  # warmup

            start = time.clock()
            iter = 0
            while iter < 5 or time.clock() < start + 1:
                y = A*x
                iter += 1
            end = time.clock()

            del y

            name = name.center(12)
            shape = ("%s" % (A.shape,)).center(20)
            MFLOPs = (2*n_vecs*A.nnz*iter/(end-start))/float(1e6)

            print(fmt % (A.format,name,shape,A.nnz,MFLOPs))

    def bench_construction(self):
        """build matrices by inserting single values"""
        matrices = []
        matrices.append(('Empty',csr_matrix((10000,10000))))
        matrices.append(('Identity',sparse.eye(10000)))
        matrices.append(('Poisson5pt', poisson2d(100)))

        print()
        print('                    Sparse Matrix Construction')
        print('====================================================================')
        print(' type |    name      |         shape        |    nnz   | time (sec) ')
        print('--------------------------------------------------------------------')
        fmt = '  %3s | %12s | %20s | %8d |   %6.4f '

        for name,A in matrices:
            A = A.tocoo()

            for format, cls in [('lil', lil_matrix), ('dok', dok_matrix)]:

                start = time.clock()

                iter = 0
                while time.clock() < start + 0.5:
                    T = cls(A.shape)
                    for i,j,v in zip(A.row,A.col,A.data):
                        T[i,j] = v
                    iter += 1
                end = time.clock()

                del T
                name = name.center(12)
                shape = ("%s" % (A.shape,)).center(20)

                print(fmt % (format,name,shape,A.nnz,(end-start)/float(iter)))

    def bench_conversion(self):
        A = poisson2d(100)

        formats = ['csr','csc','coo','dia','lil','dok']

        print()
        print('                     Sparse Matrix Conversion')
        print('====================================================================')
        print(' format | tocsr() | tocsc() | tocoo() | todia() | tolil() | todok() ')
        print('--------------------------------------------------------------------')

        for fromfmt in formats:
            base = getattr(A,'to' + fromfmt)()

            times = []

            for tofmt in formats:
                try:
                    fn = getattr(base,'to' + tofmt)
                except:
                    times.append(None)
                else:
                    x = fn()  # warmup
                    start = time.clock()
                    iter = 0
                    while time.clock() < start + 0.2:
                        x = fn()
                        iter += 1
                    end = time.clock()
                    del x
                    times.append((end - start)/float(iter))

            output = "  %3s   " % fromfmt
            for t in times:
                if t is None:
                    output += '|    n/a    '
                else:
                    output += '| %5.1fms ' % (1000*t)
            print(output)

    def _getset_bench(self, kernel, formats):
        print('==========================================================')
        print('      N | s.patt. |' + ''.join(' %7s |' % fmt for fmt in formats))
        print('----------------------------------------------------------')

        A = rand(1000, 1000, density=1e-5)

        for N in [1, 10, 100, 1000, 10000]:
            for spat in [False, True]:
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

                times = []

                for fmt in formats:
                    if fmt == 'dok' and N > 500:
                        times.append(None)
                        continue

                    base = A.asformat(fmt)

                    m = base.copy()
                    if spat:
                        kernel(m, i, j, v)

                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', SparseEfficiencyWarning)

                        iter = 0
                        total_time = 0
                        while total_time < 0.2 and iter < 5000:
                            if not spat:
                                m = base.copy()
                            a = time.clock()
                            kernel(m, i, j, v)
                            total_time += time.clock() - a
                            iter += 1

                    times.append(total_time/float(iter))

                output = " %6d | %7s " % (N, "same" if spat else "change")
                for t in times:
                    if t is None:
                        output += '|    n/a    '
                    else:
                        output += '| %5.2fms ' % (1e3*t)
                print(output)

    def bench_setitem(self):
        def kernel(A, i, j, v):
            A[i, j] = v
        print()
        print('           Sparse Matrix fancy __setitem__')
        self._getset_bench(kernel, ['csr', 'csc', 'lil', 'dok'])

    def bench_getitem(self):
        def kernel(A, i, j, v=None):
            A[i, j]
        print()
        print('           Sparse Matrix fancy __getitem__')
        self._getset_bench(kernel, ['csr', 'csc', 'lil'])

    def bench_large(self):
        H1, W1 = 1, 100000
        H2, W2 = W1, 1000
        C1 = 10
        C2 = 1000000

        print()

        random.seed(0)

        print('                  Sparse Matrix Large Matrix Multiplication')
        start = time.time()
        matrix1 = lil_matrix(zeros((H1, W1)))
        matrix2 = lil_matrix(zeros((H2, W2)))
        for i in xrange(C1):
            matrix1[random.randint(H1), random.randint(W1)] = random.rand()
        for i in xrange(C2):
            matrix2[random.randint(H2), random.randint(W2)] = random.rand()
        matrix1 = matrix1.tocsr()
        matrix2 = matrix2.tocsr()
        end = time.time()

        start = time.time()
        for i in range(100):
            matrix3 = matrix1 * matrix2
        end = time.time()
        print('==============================================================================')
        print('Matrix 1 shape | Matrix 2 shape | Matrix 1 count | Matrix 2 count | time (sec)')
        print('------------------------------------------------------------------------------')
        print('%14s | %14s | %14d | %14d | %10.3f' % (matrix1.shape, matrix2.shape, C1, C2, end-start))

if __name__ == "__main__":
    test = Tester(__file__).bench()
