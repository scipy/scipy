from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.sparse.sputils import isshape, isintlike
from scipy.sparse import isspmatrix

__all__ = ['LinearOperator', 'aslinearoperator']


class LinearOperator(object):
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

    Other Parameters
    ----------------
    rmatvec : callable f(v)
        Returns A^H * v, where A^H is the conjugate transpose of A.
    matmat : callable f(V)
        Returns A * V, where V is a dense matrix with dimensions (N,K).
    dtype : dtype
        Data type of the matrix.

    Attributes
    ----------
    args : tuple
        For linear operators describing products etc. of other linear
        operators, the operands of the binary operation.

    See Also
    --------
    aslinearoperator : Construct LinearOperators

    Notes
    -----
    The user-defined matvec() function must properly handle the case
    where v has shape (N,) as well as the (N,1) case.  The shape of
    the return type is handled internally by LinearOperator.

    LinearOperator instances can also be multiplied, added with each
    other and exponentiated, to produce a new linear operator.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse.linalg import LinearOperator
    >>> def mv(v):
    ...     return np.array([2*v[0], 3*v[1]])
    ...
    >>> A = LinearOperator((2,2), matvec=mv)
    >>> A
    <2x2 LinearOperator with unspecified dtype>
    >>> A.matvec(np.ones(2))
    array([ 2.,  3.])
    >>> A * np.ones(2)
    array([ 2.,  3.])

    """
    def __init__(self, shape, matvec, rmatvec=None, matmat=None, dtype=None):

        shape = tuple(shape)

        if not isshape(shape):
            raise ValueError('invalid shape')

        self.shape = shape
        self._matvec = matvec
        self.args = ()

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

        return np.hstack([self.matvec(col.reshape(-1,1)) for col in X.T])

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

    def __call__(self, x):
        return self*x

    def __mul__(self, x):
        if isinstance(x, LinearOperator):
            return _ProductLinearOperator(self, x)
        elif np.isscalar(x):
            return _ScaledLinearOperator(self, x)
        else:
            x = np.asarray(x)

            if x.ndim == 1 or x.ndim == 2 and x.shape[1] == 1:
                return self.matvec(x)
            elif x.ndim == 2:
                return self.matmat(x)
            else:
                raise ValueError('expected rank-1 or rank-2 array or matrix')

    def dot(self, other):
        # modeled after scipy.sparse.base.dot
        return self * other

    def __rmul__(self, x):
        if np.isscalar(x):
            return _ScaledLinearOperator(self, x)
        else:
            return NotImplemented

    def __pow__(self, p):
        if np.isscalar(p):
            return _PowerLinearOperator(self, p)
        else:
            return NotImplemented

    def __add__(self, x):
        if isinstance(x, LinearOperator):
            return _SumLinearOperator(self, x)
        else:
            return NotImplemented

    def __neg__(self):
        return _ScaledLinearOperator(self, -1)

    def __sub__(self, x):
        return self.__add__(-x)

    def __repr__(self):
        M,N = self.shape
        if hasattr(self,'dtype'):
            dt = 'dtype=' + str(self.dtype)
        else:
            dt = 'unspecified dtype'

        return '<%dx%d %s with %s>' % (M, N, self.__class__.__name__, dt)


def _get_dtype(operators, dtypes=[]):
    for obj in operators:
        if obj is not None and hasattr(obj, 'dtype'):
            dtypes.append(obj.dtype)
    return np.find_common_type(dtypes, [])


class _SumLinearOperator(LinearOperator):
    def __init__(self, A, B):
        if not isinstance(A, LinearOperator) or \
                not isinstance(B, LinearOperator):
            raise ValueError('both operands have to be a LinearOperator')
        if A.shape != B.shape:
            raise ValueError('shape mismatch')
        super(_SumLinearOperator, self).__init__(A.shape,
                self.matvec, self.rmatvec, self.matmat, _get_dtype([A,B]))
        self.args = (A, B)

    def matvec(self, x):
        return self.args[0].matvec(x) + self.args[1].matvec(x)

    def rmatvec(self, x):
        return self.args[0].rmatvec(x) + self.args[1].rmatvec(x)

    def matmat(self, x):
        return self.args[0].matmat(x) + self.args[1].matmat(x)


class _ProductLinearOperator(LinearOperator):
    def __init__(self, A, B):
        if not isinstance(A, LinearOperator) or \
                not isinstance(B, LinearOperator):
            raise ValueError('both operands have to be a LinearOperator')
        if A.shape[1] != B.shape[0]:
            raise ValueError('shape mismatch')
        super(_ProductLinearOperator, self).__init__((A.shape[0], B.shape[1]),
                self.matvec, self.rmatvec, self.matmat, _get_dtype([A,B]))
        self.args = (A, B)

    def matvec(self, x):
        return self.args[0].matvec(self.args[1].matvec(x))

    def rmatvec(self, x):
        return self.args[1].rmatvec(self.args[0].rmatvec(x))

    def matmat(self, x):
        return self.args[0].matmat(self.args[1].matmat(x))


class _ScaledLinearOperator(LinearOperator):
    def __init__(self, A, alpha):
        if not isinstance(A, LinearOperator):
            raise ValueError('LinearOperator expected as A')
        if not np.isscalar(alpha):
            raise ValueError('scalar expected as alpha')
        super(_ScaledLinearOperator, self).__init__(A.shape,
                self.matvec, self.rmatvec, self.matmat,
                _get_dtype([A], [type(alpha)]))
        self.args = (A, alpha)

    def matvec(self, x):
        return self.args[1] * self.args[0].matvec(x)

    def rmatvec(self, x):
        return np.conj(self.args[1]) * self.args[0].rmatvec(x)

    def matmat(self, x):
        return self.args[1] * self.args[0].matmat(x)


class _PowerLinearOperator(LinearOperator):
    def __init__(self, A, p):
        if not isinstance(A, LinearOperator):
            raise ValueError('LinearOperator expected as A')
        if A.shape[0] != A.shape[1]:
            raise ValueError('square LinearOperator expected as A')
        if not isintlike(p):
            raise ValueError('integer expected as p')
        super(_PowerLinearOperator, self).__init__(A.shape,
                self.matvec, self.rmatvec, self.matmat,
                _get_dtype([A]))
        self.args = (A, p)

    def _power(self, fun, x):
        res = np.array(x, copy=True)
        for i in range(self.args[1]):
            res = fun(res)
        return res

    def matvec(self, x):
        return self._power(self.args[0].matvec, x)

    def rmatvec(self, x):
        return self._power(self.args[0].rmatvec, x)

    def matmat(self, x):
        return self._power(self.args[0].matmat, x)


class MatrixLinearOperator(LinearOperator):
    def __init__(self, A):
        super(MatrixLinearOperator, self).__init__(shape=A.shape,
                dtype=A.dtype, matvec=None, rmatvec=self.rmatvec)
        self.matvec = A.dot
        self.matmat = A.dot
        self.__mul__ = A.dot
        self.A = A
        self.A_conj = None
        self.args = (A,)

    def rmatvec(self, x):
        if self.A_conj is None:
            self.A_conj = self.A.T.conj()
        return self.A_conj.dot(x)


class IdentityOperator(LinearOperator):
    def __init__(self, shape, dtype):
        super(IdentityOperator, self).__init__(shape=shape, dtype=dtype,
                matvec=None, rmatvec=self.rmatvec)

    def matvec(self, x):
        return x

    def rmatvec(self, x):
        return x

    def matmat(self, x):
        return x

    def __mul__(self, x):
        return x


def aslinearoperator(A):
    """Return A as a LinearOperator.

    'A' may be any of the following types:
     - ndarray
     - matrix
     - sparse matrix (e.g. csr_matrix, lil_matrix, etc.)
     - LinearOperator
     - An object with .shape and .matvec attributes

    See the LinearOperator documentation for additional information.

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
        return MatrixLinearOperator(A)

    elif isspmatrix(A):
        return MatrixLinearOperator(A)

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
