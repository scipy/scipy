import numpy as np
from scipy.sparse.sputils import isshape
from scipy.sparse import isspmatrix

__all__ = ['LinearOperator', 'aslinearoperator']

class LinearOperator:
    def __init__( self, shape, matvec, rmatvec=None, dtype=None ):
        """Common interface for performing matrix vector products
        """
        shape = tuple(shape)

        if not isshape(shape):
            raise ValueError('invalid shape')

        self.shape  = shape
        self.matvec = matvec

        if rmatvec is None:
            def rmatvec(x):
                raise NotImplementedError('rmatvec is not defined')
            self.rmatvec = rmatvec
        else:
            self.rmatvec = rmatvec

        if dtype is not None:
            self.dtype = np.dtype(dtype)

    def __repr__(self):
        M,N = self.shape
        if hasattr(self,'dtype'):
            dt = 'dtype=' + str(self.dtype)
        else:
            dt = 'unspecified dtype'

        return '<%dx%d LinearOperator with %s>' % (M,N,dt)

def aslinearoperator(A):
    """Return A as a LinearOperator

    'A' may be any of the following types
        - ndarray
        - matrix
        - sparse matrix (e.g. csr_matrix, lil_matrix, etc.)
        - LinearOperator
        - An object with .shape and .matvec attributes

    See the LinearOperator documentation for additonal information.

    Examples
    ========
    
    >>> from scipy import matrix
    >>> M = matrix( [[1,2,3],[4,5,6]], dtype='int32' )
    >>> aslinearoperator( M )
    <2x3 LinearOperator with dtype=int32>

    """

    if isinstance(A, LinearOperator):
        return A

    elif isinstance(A, np.ndarray) or isinstance(A,np.matrix):
        def matvec(x):
            return np.dot(np.asarray(A),x)
        def rmatvec(x):
            return np.dot(x,np.asarray(A))
        return LinearOperator( A.shape, matvec, rmatvec=rmatvec, dtype=A.dtype )

    elif isspmatrix(A):
        return LinearOperator( A.shape, A.matvec, rmatvec=A.rmatvec, dtype=A.dtype ) 

    else:
        if hasattr(A,'shape') and hasattr(A,'matvec'):
            rmatvec = None
            dtype = None

            if hasattr(A,'rmatvec'):
                rmatvec = A.rmatvec
            if hasattr(A,'dtype'):
                dtype = A.dtype
            return LinearOperator(A.shape, A.matvec, rmatvec=rmatvec, dtype=dtype)

        else:
            raise TypeError('type not understood')


