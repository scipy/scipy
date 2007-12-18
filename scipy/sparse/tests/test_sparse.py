"""general tests and simple benchmarks for the sparse module"""

import numpy
from numpy import ones, array, asarray, empty

import random
from numpy.testing import *
set_package_path()
from scipy.sparse import csc_matrix, csr_matrix, dok_matrix, \
        coo_matrix, lil_matrix, dia_matrix, spidentity, spdiags
from scipy.linsolve import splu
restore_path()

#TODO move this to a matrix gallery and add unittests
def poisson2d(N,dtype='d',format=None):
    """
    Return a sparse matrix for the 2d poisson problem
    with standard 5-point finite difference stencil on a
    square N-by-N grid.
    """
    if N == 1:
        diags   = asarray( [[4]],dtype=dtype)
        return dia_matrix((diags,[0]), shape=(1,1)).asformat(format)

    offsets = array([0,-N,N,-1,1])

    diags = empty((5,N**2),dtype=dtype)

    diags[0]  =  4 #main diagonal
    diags[1:] = -1 #all offdiagonals

    diags[3,N-1::N] = 0  #first lower diagonal
    diags[4,N::N] = 0  #first upper diagonal

    return dia_matrix((diags,offsets),shape=(N**2,N**2)).asformat(format)

import time
class TestSparseTools(NumpyTestCase):
    """Simple benchmarks for sparse matrix module"""

    def test_matvec(self,level=5):
        matrices = []
        matrices.append(('Identity',   spidentity(10**5,format='csr')))
        matrices.append(('Poisson5pt', poisson2d(1000,format='csr')))
        matrices.append(('Poisson5pt', poisson2d(1000,format='dia')))

        print
        print '                 Sparse Matrix Vector Product'
        print '=================================================================='
        print ' type |    name      |         shape        |    nnz   |  MFLOPs  '
        print '------------------------------------------------------------------'
        fmt = '  %3s | %12s | %20s | %8d |  %6.1f '

        for name,A in matrices:
            x = ones(A.shape[1],dtype=A.dtype)

            y = A*x  #warmup

            start = time.clock()
            iter = 0
            while iter < 5 or time.clock() < start + 1:
                y = A*x
                iter += 1
            end = time.clock()

            name = name.center(12)
            shape = ("%s" % (A.shape,)).center(20)
            MFLOPs = (2*A.nnz*iter/(end-start))/float(1e6)

            print fmt % (A.format,name,shape,A.nnz,MFLOPs)
            
    def bench_construction(self,level=5):
        """build matrices by inserting single values"""
        matrices = []
        matrices.append( ('Empty',csr_matrix((10000,10000))) )
        matrices.append( ('Identity',spidentity(10000)) )
        matrices.append( ('Poisson5pt', poisson2d(100)) )
        
        print
        print '                    Sparse Matrix Construction'
        print '===================================================================='
        print ' type |    name      |         shape        |    nnz   | time (sec) '
        print '--------------------------------------------------------------------'
        fmt = '  %3s | %12s | %20s | %8d |   %6.4f '

        for name,A in matrices:
            A = A.tocoo()
             
            for format in ['lil','dok']: 

                start = time.clock()
                
                iter = 0
                while time.clock() < start + 0.1:
                    T = eval(format + '_matrix')(A.shape)
                    for i,j,v in zip(A.row,A.col,A.data):
                        T[i,j] = v
                    iter += 1
                end = time.clock()

                del T
                name = name.center(12)
                shape = ("%s" % (A.shape,)).center(20)

                print fmt % (format,name,shape,A.nnz,(end-start)/float(iter))


    def bench_conversion(self,level=5):
        A = poisson2d(100)

        formats = ['csr','csc','coo','lil','dok']
       
        print
        print '                Sparse Matrix Conversion'
        print '=========================================================='
        print ' format | tocsr() | tocsc() | tocoo() | tolil() | todok() '
        print '----------------------------------------------------------'
        
        for fromfmt in formats:
            base = getattr(A,'to' + fromfmt)()
 
            times = []

            for tofmt in formats:
                try:
                    fn = getattr(base,'to' + tofmt)
                except:
                    times.append(None)
                else:
                    x = fn() #warmup
                    start = time.clock()
                    iter = 0
                    while time.clock() < start + 0.2:
                        x = fn()
                        iter += 1
                    end = time.clock()
                    del x 
                    times.append( (end - start)/float(iter))

            output = "  %3s   " % fromfmt
            for t in times:
                if t is None:
                    output += '|    n/a    '
                else:
                    output += '| %5.1fms ' % (1000*t) 
            print output


                
if __name__ == "__main__":
    NumpyTest().run()

