import numpy as np
from scipy.sparse.sputils import isshape
from scipy.sparse import isspmatrix

__all__ = ['LinearOperator', 'aslinearoperator']

class LinearOperator:
    """Common interface for performing matrix vector products

    Many iterative methods (e.g. cg, gmres) do not need to know the
    individual entries of a matrix to solve a linear system A*x=b.
    Such solvers only require the computation of matrix vector
    products, A*v where v is a dense vector.  This class serves as
    an abstract interface between iterative solvers and matrix-like
    objects.

    Parameters
    ----------
    shape : tuple
        Matrix dimensions (M,N)
    matvec : callable f(v)
        Returns returns A * v.

    Optional Parameters
    -------------------
    rmatvec : callable f(v)
        Returns A^H * v, where A^H is the conjugate transpose of A.
    matmat : callable f(V)
        Returns A * V, where V is a dense matrix with dimensions (N,K).
    dtype : dtype
        Data type of the matrix.

    See Also
    --------
    aslinearoperator : Construct LinearOperators

    Notes
    -----
    The user-defined matvec() function must properly handle the case
    where v has shape (N,) as well as the (N,1) case.  The shape of
    the return type is handled internally by LinearOperator.

    Examples
    --------
    >>> from scipy.sparse.linalg import LinearOperator
    >>> from scipy import *
    >>> def mv(v):
    ...     return array([ 2*v[0], 3*v[1]])
    ...
    >>> A = LinearOperator( (2,2), matvec=mv )
    >>> A
    <2x2 LinearOperator with unspecified dtype>
    >>> A.matvec( ones(2) )
    array([ 2.,  3.])
    >>> A * ones(2)
    array([ 2.,  3.])

    """
    def __init__(self, shape, matvec, rmatvec=None, matmat=None, dtype=None):

        shape = tuple(shape)

        if not isshape(shape):
            raise ValueError('invalid shape')

        self.shape  = shape
        self._matvec = matvec

        if rmatvec is None:
            def rmatvec(v):
                raise NotImplementedError('rmatvec is not defined')
            self.rmatvec = rmatvec
        else:
            self.rmatvec = rmatvec

        if matmat is not None:
            # matvec each column of V
            self._matmat = matmat

        if dtype is not None:
            self.dtype = np.dtype(dtype)


    def _matmat(self, X):
        """Default matrix-matrix multiplication handler.  Falls back on
        the user-defined matvec() routine, which is always provided.
        """

        return np.hstack( [ self.matvec(col.reshape(-1,1)) for col in X.T ] )


    def matvec(self, x):
        """Matrix-vector multiplication

        Performs the operation y=A*x where A is an MxN linear
        operator and x is a column vector or rank-1 array.

        Parameters
        ----------
        x : {matrix, ndarray}
            An array with shape (N,) or (N,1).

        Returns
        -------
        y : {matrix, ndarray}
            A matrix or ndarray with shape (M,) or (M,1) depending
            on the type and shape of the x argument.

        Notes
        -----
        This matvec wraps the user-specified matvec routine to ensure that
        y has the correct shape and type.

        """

        x = np.asanyarray(x)

        M,N = self.shape

        if x.shape != (N,) and x.shape != (N,1):
            raise ValueError('dimension mismatch')

        y = self._matvec(x)

        if isinstance(x, np.matrix):
            y = np.asmatrix(y)
        else:
            y = np.asarray(y)

        if x.ndim == 1:
            y = y.reshape(M)
        elif x.ndim == 2:
            y = y.reshape(M,1)
        else:
            raise ValueError('invalid shape returned by user-defined matvec()')


        return y


    def matmat(self, X):
        """Matrix-matrix multiplication

        Performs the operation y=A*X where A is an MxN linear
        operator and X dense N*K matrix or ndarray.

        Parameters
        ----------
        X : {matrix, ndarray}
            An array with shape (N,K).

        Returns
        -------
        Y : {matrix, ndarray}
            A matrix or ndarray with shape (M,K) depending on
            the type of the X argument.

        Notes
        -----
        This matmat wraps any user-specified matmat routine to ensure that
        y has the correct type.

        """

        X = np.asanyarray(X)

        if X.ndim != 2:
            raise ValueError('expected rank-2 ndarray or matrix')

        M,N = self.shape

        if X.shape[0] != N:
            raise ValueError('dimension mismatch')

        Y = self._matmat(X)

        if isinstance(Y, np.matrix):
            Y = np.asmatrix(Y)

        return Y


    def __mul__(self,x):
        x = np.asarray(x)

        if x.ndim == 1 or x.ndim == 2 and x.shape[1] == 1:
            return self.matvec(x)
        elif x.ndim == 2:
            return self.matmat(x)
        else:
            raise ValueError('expected rank-1 or rank-2 array or matrix')


    def __repr__(self):
        M,N = self.shape
        if hasattr(self,'dtype'):
            dt = 'dtype=' + str(self.dtype)
        else:
            dt = 'unspecified dtype'

        return '<%dx%d LinearOperator with %s>' % (M,N,dt)

def aslinearoperator(A):
    """Return A as a LinearOperator.

    'A' may be any of the following types:
     - ndarray
     - matrix
     - sparse matrix (e.g. csr_matrix, lil_matrix, etc.)
     - LinearOperator
     - An object with .shape and .matvec attributes

    See the LinearOperator documentation for additonal information.

    Examples
    --------
    >>> from scipy import matrix
    >>> M = matrix( [[1,2,3],[4,5,6]], dtype='int32' )
    >>> aslinearoperator( M )
    <2x3 LinearOperator with dtype=int32>

    """
    if isinstance(A, LinearOperator):
        return A

    elif isinstance(A, np.ndarray) or isinstance(A, np.matrix):
        if A.ndim > 2:
            raise ValueError('array must have rank <= 2')

        A = np.atleast_2d(np.asarray(A))

        def matvec(v):
            return np.dot(A, v)
        def rmatvec(v):
            return np.dot(A.conj().transpose(), v)
        def matmat(V):
            return np.dot(A, V)
        return LinearOperator(A.shape, matvec, rmatvec=rmatvec,
                              matmat=matmat, dtype=A.dtype)

    elif isspmatrix(A):
        def matvec(v):
            return A * v
        def rmatvec(v):
            return A.conj().transpose() * v
        def matmat(V):
            return A * V
        return LinearOperator(A.shape, matvec, rmatvec=rmatvec,
                              matmat=matmat, dtype=A.dtype)

    else:
        if hasattr(A, 'shape') and hasattr(A, 'matvec'):
            rmatvec = None
            dtype = None

            if hasattr(A, 'rmatvec'):
                rmatvec = A.rmatvec
            if hasattr(A, 'dtype'):
                dtype = A.dtype
            return LinearOperator(A.shape, A.matvec,
                                  rmatvec=rmatvec, dtype=dtype)

        else:
            raise TypeError('type not understood')
