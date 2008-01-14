__all__ = ['poisson']

from scipy import arange, empty, intc, ravel, prod
from scipy.sparse import coo_matrix


def poisson( grid, spacing=None, dtype=float, format=None):
    """Finite Difference approximation to the Poisson problem on a 
    regular n-dimensional grid with Dirichlet boundary conditions.
   
    
    Parameters
    ==========
        - grid : tuple
            - grid dimensions e.g. (100,100)


    Examples
    ========

    >>> # 4 nodes in one dimension
    >>> poisson( (4,) ).todense()
    matrix([[ 2., -1.,  0.,  0.],
            [-1.,  2., -1.,  0.],
            [ 0., -1.,  2., -1.],
            [ 0.,  0., -1.,  2.]])

    >>> # rectangular two dimensional grid 
    >>> poisson( (2,3) ).todense()
    matrix([[ 4., -1.,  0., -1.,  0.,  0.],
            [-1.,  4., -1.,  0., -1.,  0.],
            [ 0., -1.,  4.,  0.,  0., -1.],
            [-1.,  0.,  0.,  4., -1.,  0.],
            [ 0., -1.,  0., -1.,  4., -1.],
            [ 0.,  0., -1.,  0., -1.,  4.]])

    """
    grid = tuple(grid)

    D = len(grid) # grid dimension

    if D < 1 or min(grid) < 1:
        raise ValueError,'invalid grid shape: %s' % str(grid)

    nodes = arange(prod(grid)).reshape(*grid)

    nnz = nodes.size 
    for i in range(D):
        nnz += 2 * prod( grid[:i] + grid[i+1:] ) * (grid[i] - 1)
    
    row  = empty(nnz, dtype=intc)
    col  = empty(nnz, dtype=intc)
    data = empty(nnz, dtype=dtype)
    
    row[:nodes.size]  = ravel(nodes)
    col[:nodes.size]  = ravel(nodes)
    data[:nodes.size] = 2*D
    data[nodes.size:] = -1
    
    ptr = nodes.size
    
    for i in range(D):
        s0 = [slice(None)] * i + [slice(0,-1)  ] + [slice(None)] * (D - i - 1)
        s1 = [slice(None)] * i + [slice(1,None)] + [slice(None)] * (D - i - 1)
    
        n0 = nodes[s0]
        n1 = nodes[s1]
    
        row0 = row[ ptr:ptr + n0.size].reshape(n0.shape)
        col0 = col[ ptr:ptr + n0.size].reshape(n0.shape)
        ptr += n0.size
    
        row1 = row[ ptr:ptr + n0.size].reshape(n0.shape)
        col1 = col[ ptr:ptr + n0.size].reshape(n0.shape)
        ptr += n0.size
    
        row0[:] = n0
        col0[:] = n1
    
        row1[:] = n1
        col1[:] = n0
    
    return coo_matrix((data,(row,col)),shape=(nodes.size,nodes.size)).asformat(format)

