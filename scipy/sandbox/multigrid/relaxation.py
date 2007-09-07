import multigridtools
import numpy

def gauss_seidel(A,x,b,iterations=1,sweep='forward'):
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
    if A.shape[0] != A.shape[1]:
        raise ValueError,'expected symmetric matrix'

    if A.shape[1] != len(x) or len(x) != len(b):
        raise ValueError,'unexpected number of unknowns'

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
    Perform Jacobi iteration on the linear system Ax=b
 
       x <- (1 - omega) x  +  omega * D^-1 (b - (A - D) x)
    
    where D is the diagonal of A.
    
    Input:
         A - NxN csr_matrix
         x - rank 1 ndarray of length N
         b - rank 1 ndarray of length N
     Optional:
         iterations - number of iterations to perform (default: 1)
         omega      - damping parameter (default: 1.0)
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


def polynomial_smoother(A,x,b,coeffs):
    """
    Apply a polynomial smoother to the system Ax=b

    The smoother has the form:
      x_new = x + p(A) (b - A*x)
    where p(A) is a polynomial in A whose scalar coeffients
    are specified (in decending order) by argument coeffs.

    Eg.

      Richardson iteration p(A) = c_0:
         polynomial_smoother(A,x,b,[c_0])

      Linear smoother p(A) = c_1*A + c_0:
         polynomial_smoother(A,x,b,[c_1,c_0])

      Quadratic smoother p(A) = c_2*A^2 + c_1*A + c_0:
         polynomial_smoother(A,x,b,[c_2,c_1,c_0])


    Note: Horner's Rule is applied to avoid computing A^k directly.
    """

    residual = (b - A*x)
    h = coeffs[0]*residual
    
    for c in coeffs[1:]:
        h = c*residual + A*h
        
    x += h
    

