import multigridtools
from numpy import empty_like


def sor(A,x,b,omega,iterations=1,sweep='forward'):
    """
    Perform SOR iteration on the linear system Ax=b
    """
    x_old = empty_like(x)
        
    for i in range(iterations):
        x_old[:] = x
        gauss_seidel(A,x,b,iterations=1,sweep=sweep)
        
        x     *= omega
        x_old *= (1-omega)
        x     += x_old
        
        
def gauss_seidel(A,x,b,iterations=1,sweep='forward'):
    """
    Perform Gauss-Seidel iteration on the linear system Ax=b
 
     Input:
         A - NxN csr_matrix
         x - rank 1 ndarray of length N
         b - rank 1 ndarray of length N
     Optional:
         iterations - number of iterations to perform (default: 1)
         sweep      - direction of sweep:
                        'forward' (default), 'backward', or 'symmetric'
    """ 
    if A.shape[0] != A.shape[1]:
        raise ValueError,'expected symmetric matrix'

    if A.shape[1] != len(x) or len(x) != len(b):
        raise ValueError,'unexpected number of unknowns'

    if sweep == 'forward':
        row_start,row_stop,row_step = 0,len(x),1
        for iter in xrange(iterations):
            multigridtools.gauss_seidel(A.shape[0],
                                        A.indptr, A.indices, A.data,
                                        x, b,
                                        row_start, row_stop, row_step)
    elif sweep == 'backward':
        row_start,row_stop,row_step = len(x)-1,-1,-1
        for iter in xrange(iterations):
            multigridtools.gauss_seidel(A.shape[0],
                                        A.indptr, A.indices, A.data,
                                        x, b,
                                        row_start, row_stop, row_step)
    elif sweep == 'symmetric':
        for iter in xrange(iterations):
            gauss_seidel(A,x,b,iterations=1,sweep='forward')
            gauss_seidel(A,x,b,iterations=1,sweep='backward')
    else:
       raise ValueError,'valid sweep directions are \'forward\', \'backward\', and \'symmetric\''


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

    temp = empty_like(x)
    
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

    #TODO skip first matvec if x is all zero
    
    residual = (b - A*x)
    h = coeffs[0]*residual
    
    for c in coeffs[1:]:
        h = c*residual + A*h
        
    x += h
    

