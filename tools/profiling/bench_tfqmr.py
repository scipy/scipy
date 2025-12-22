import numpy as np
import scipy

def _create_sparse_poisson1d(n):
    # Make Gilbert Strang's favorite matrix
    # http://www-math.mit.edu/~gs/PIX/cupcakematrix.jpg
    P1d = scipy.sparse.diags(
        [[-1]*(n-1), [2]*n, [-1]*(n-1)], [-1, 0, 1],
        dtype=np.float64
    )
    return P1d

def _create_sparse_poisson2d(n):
    P1d = _create_sparse_poisson1d(n)
    P2d = scipy.sparse.kronsum(P1d, P1d)
    return P2d.tocsr()

n = 500
b = np.ones(n*n)
P_sparse = _create_sparse_poisson2d(n)

scipy.sparse.linalg.tfqmr(P_sparse, b)
