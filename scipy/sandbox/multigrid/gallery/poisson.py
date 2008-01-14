__all__ = ['poisson']

from scipy import array, empty 
from scipy.sparse import dia_matrix

def poisson(N, stencil='5pt', dtype=float, format=None):
    """Finite Difference approximations to the Poisson problem

    TheDirichlet boundary conditions are
   
    
    Parameters
    ==========
        - N : integer 
            - grid size
        - stencil : one of the following strings
            - '3pt' : 3-point Finite Difference stencil in 1 dimension
            - '5pt' : 5-point Finite Difference stencil in 2 dimensions
            - '8pt' : NotImplemented
            - '27pt' : NotImplemented

    """
    if N < 1:
        raise ValueError,'invalid grid size %s' % N

    if stencil == '3pt':
        if N == 1:
            diags   = array( [[2]], dtype=dtype)
            return dia_matrix((diags,[0]), shape=(1,1)).asformat(format)
        else:
            data = empty((3,N),dtype=dtype)
            data[0,:] = 2 #main diagonal
            data[1,:] = -1 
            data[2,:] = -1
    
            return dia_matrix((data,[0,-1,1]),shape=(N,N)).asformat(format)
    elif stencil == '5pt':
        if N == 1:
            data = array( [[4]], dtype=dtype)
            return dia_matrix((diags,[0]), shape=(1,1)).asformat(format)
        else:
            diags = array([0,-N,N,-1,1])

            data = empty((5,N**2),dtype=dtype)

            data[0]  =  4 #main diagonal
            data[1::,:] = -1
            data[3,N-1::N] = 0 
            data[4,N::N]   = 0    
            
            return dia_matrix((data,diags),shape=(N**2,N**2)).asformat(format)
    else:
        raise NotImplementedError,'unsupported stencil=%s' % stencil


