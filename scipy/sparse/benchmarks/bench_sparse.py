"""general tests and simple benchmarks for the sparse module"""
from __future__ import division, print_function, absolute_import

import time

import numpy
from numpy import ones, array, asarray, empty

from numpy.testing import *

from scipy import sparse
from scipy.lib.six.moves import xrange
from scipy.sparse import csr_matrix, coo_matrix, dia_matrix, lil_matrix, \
        dok_matrix


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
            vars = dict([(var,mat.asformat(format)) for (var,name,mat) in matrices])
            for X,Y in [('A','A'),('A','B'),('B','A'),('B','B')]:
                x,y = vars[X],vars[Y]
                for op in ['__add__','__sub__','multiply','__div__','__mul__']:
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

            for format in ['lil','dok']:

                start = time.clock()

                iter = 0
                while time.clock() < start + 0.5:
                    T = eval(format + '_matrix')(A.shape)
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


# class TestLarge(TestCase):
#    def bench_large(self):
#        # Create a 100x100 matrix with 100 non-zero elements
#        # and play around with it
#        #TODO move this out of Common since it doesn't use spmatrix
#        random.seed(0)
#        A = dok_matrix((100,100))
#        for k in xrange(100):
#            i = random.randrange(100)
#            j = random.randrange(100)
#            A[i,j] = 1.
#        csr = A.tocsr()
#        csc = A.tocsc()
#        csc2 = csr.tocsc()
#        coo = A.tocoo()
#        csr2 = coo.tocsr()
#        assert_array_equal(A.transpose().todense(), csr.transpose().todense())
#        assert_array_equal(csc.todense(), csr.todense())
#        assert_array_equal(csr.todense(), csr2.todense())
#        assert_array_equal(csr2.todense().transpose(), coo.transpose().todense())
#        assert_array_equal(csr2.todense(), csc2.todense())
#        csr_plus_csc = csr + csc
#        csc_plus_csr = csc + csr
#        assert_array_equal(csr_plus_csc.todense(), (2*A).todense())
#        assert_array_equal(csr_plus_csc.todense(), csc_plus_csr.todense())


if __name__ == "__main__":
    run_module_suite()
