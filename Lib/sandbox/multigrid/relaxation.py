import multigridtools
import numpy

def gauss_seidel(A,x,b,iterations=1,sweep="forward"):
    """
    Perform Gauss-Seidel iteration on the linear system Ax=b
 
     Input:
         A - NxN csr_matrix
         x - rank 1 ndarray of length N
         b - rank 1 ndarray of length N
     Optional:
         iterations - number of iterations to perform (default: 1)
         sweep      - slice of unknowns to relax (default: all in forward direction)
    """ 
    if sweep == 'forward':
        row_start,row_stop,row_step = 0,len(x),1
    elif sweep == 'backward':
        row_start,row_stop,row_step = len(x)-1,-1,-1
    else:
       raise ValueError,'valid sweep directions are \'forward\' and \'backward\''

    for iter in xrange(iterations):
        multigridtools.gauss_seidel(A.shape[0],
                                    A.indptr, A.indices, A.data,
                                    x, b,
                                    row_start, row_stop, row_step)

def jacobi(A,x,b,iterations=1,omega=1.0):
    """
    Perform Gauss-Seidel iteration on the linear system Ax=b
 
     Input:
         A - NxN csr_matrix
         x - rank 1 ndarray of length N
         b - rank 1 ndarray of length N
     Optional:
         iterations - number of iterations to perform (default: 1)
         sweep      - slice of unknowns to relax (default: all in forward direction)
    """ 
    sweep = slice(None)
    (row_start,row_stop,row_step) = sweep.indices(A.shape[0])
    
    if (row_stop - row_start) * row_step <= 0:  #no work to do
        return

    temp = numpy.empty_like(x)
    
    for iter in xrange(iterations):
        multigridtools.jacobi(A.shape[0],
                              A.indptr, A.indices, A.data,
                              x, b, temp,
                              row_start, row_stop, row_step,
                              omega)


