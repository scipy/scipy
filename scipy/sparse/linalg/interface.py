from __future__ import division, print_function, absolute_import

from abc import ABCMeta, abstractmethod

import numpy as np

from scipy.lib.six import with_metaclass
from scipy.sparse import isspmatrix
from scipy.sparse.sputils import isshape, isintlike

__all__ = ['LinearOperator', 'aslinearoperator']


class LinearOperator(with_metaclass(ABCMeta, object)):
    """Common interface for performing matrix vector products

    Many iterative methods (e.g. cg, gmres) do not need to know the
    individual entries of a matrix to solve a linear system A*x=b.
    Such solvers only require the computation of matrix vector
    products, A*v where v is a dense vector.  This class serves as
    an abstract interface between iterative solvers and matrix-like
    objects.

    To construct a concrete LinearOperator, either pass appropriate
    callables to the constructor of this class, or subclass it and
    implement at least the method ``_matvec`` and the attributes or
    properties ``shape`` (pair of integers) and dtype (may be None),
    and optionally the methods ``_rmatvec``, ``_matmat`` and
    ``_adjoint``.

    Parameters
    ----------
    shape : tuple
        Matrix dimensions (M,N).
    matvec : callable f(v)
        Returns returns A * v.
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
    def __new__(cls, *args, **kwargs):
        if cls is LinearOperator:
            # Operate as _CustomLinearOperator factory.
            return _CustomLinearOperator(*args, **kwargs)
        else:
            obj = super(LinearOperator, cls).__new__(cls)
            obj.__init__(*args, **kwargs)
            return obj

    def _matmat(self, X):
        """Default matrix-matrix multiplication handler.  Falls back on
        the user-defined matvec() routine, which is always provided.
        """

        return np.hstack([self.matvec(col.reshape(-1,1)) for col in X.T])

    @abstractmethod
    def _matvec(self, x):
        """Matrix-vector multiplication.

        If self is a linear operator of shape (M, N), then this method will
        be called on a shape (N,) or (N, 1) ndarray, and should return a
        shape (M,) or (M, 1) ndarray.
        """

    def matvec(self, x):
        """Matrix-vector multiplication.

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
        This matvec wraps the user-specified matvec routine or overridden
        _matvec method to ensure that y has the correct shape and type.

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

    def rmatvec(self, x):
        """Adjoint matrix-vector multiplication.

        Performs the operation y = A^H * x where A is an MxN linear
        operator and x is a column vector or rank-1 array.

        Parameters
        ----------
        x : {matrix, ndarray}
            An array with shape (M,) or (M,1).

        Returns
        -------
        y : {matrix, ndarray}
            A matrix or ndarray with shape (N,) or (N,1) depending
            on the type and shape of the x argument.

        Notes
        -----
        This rmatvec wraps the user-specified rmatvec routine or overridden
        _rmatvec method to ensure that y has the correct shape and type.

        """

        x = np.asanyarray(x)

        M,N = self.shape

        if x.shape != (M,) and x.shape != (M,1):
            raise ValueError('dimension mismatch')

        y = self._rmatvec(x)

        if isinstance(x, np.matrix):
            y = np.asmatrix(y)
        else:
            y = np.asarray(y)

        if x.ndim == 1:
            y = y.reshape(N)
        elif x.ndim == 2:
            y = y.reshape(N,1)
        else:
            raise ValueError('invalid shape returned by user-defined rmatvec()')

        return y

    def matmat(self, X):
        """Matrix-matrix multiplication.

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
        This matmat wraps any user-specified matmat routine or overridden
        _matmat method to ensure that y has the correct type.

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
        if self.dtype is None:
            dt = 'unspecified dtype'
        else:
            dt = 'dtype=' + str(self.dtype)

        return '<%dx%d %s with %s>' % (M, N, self.__class__.__name__, dt)

    def adjoint(self):
        """Hermitian adjoint.

        Returns the Hermitian adjoint of self, aka the Hermitian
        conjugate or Hermitian transpose. For a complex matrix, the
        Hermitian adjoint is equal to the conjugate transpose.

        Can be abbreviated self.H instead of self.adjoint().

        Returns
        -------
        A_H : LinearOperator
            Hermitian adjoint of self.
        """
        return self._adjoint()

    H = property(adjoint)


class _CustomLinearOperator(LinearOperator):
    """Linear operator defined in terms of user-specified operations."""

    def __init__(self, shape, matvec, rmatvec=None, matmat=None, dtype=None):
        shape = tuple(shape)

        if not isshape(shape):
            raise ValueError('invalid shape')

        self.args = ()
        self.shape = shape

        self.__matvec_impl = matvec
        self.__rmatvec_impl = rmatvec
        self.__matmat_impl = matmat

        if dtype is None:
            self.dtype = None
        else:
            self.dtype = np.dtype(dtype)

    def _matmat(self, X):
        if self.__matmat_impl is not None:
            return self.__matmat_impl(X)
        else:
            return super(_CustomLinearOperator, self)._matmat(X)

    def _matvec(self, x):
        return self.__matvec_impl(x)

    def _rmatvec(self, x):
        func = self.__rmatvec_impl
        if func is None:
            raise NotImplemented("rmatvec is not defined")
        return self.__rmatvec_impl(x)

    def _adjoint(self):
        return _CustomLinearOperator(shape=(self.shape[1], self.shape[0]),
                                     matvec=self.__rmatvec_impl,
                                     rmatvec=self.__matvec_impl,
                                     dtype=self.dtype)


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
            raise ValueError('cannot add %r and %r: shape mismatch'
                             % (A, B))
        self.args = (A, B)
        self.dtype = _get_dtype([A, B])
        self.shape = A.shape

    def _matvec(self, x):
        return self.args[0].matvec(x) + self.args[1].matvec(x)

    def _rmatvec(self, x):
        return self.args[0].rmatvec(x) + self.args[1].rmatvec(x)

    def _matmat(self, x):
        return self.args[0].matmat(x) + self.args[1].matmat(x)

    def _adjoint(self):
        A, B = self.args
        return A.H + B.H


class _ProductLinearOperator(LinearOperator):
    def __init__(self, A, B):
        if not isinstance(A, LinearOperator) or \
                not isinstance(B, LinearOperator):
            raise ValueError('both operands have to be a LinearOperator')
        if A.shape[1] != B.shape[0]:
            raise ValueError('cannot multiply %r and %r: shape mismatch'
                             % (A, B))
        self.shape = (A.shape[0], B.shape[1])
        self.args = (A, B)
        self.dtype = _get_dtype([A, B])

    def _matvec(self, x):
        return self.args[0].matvec(self.args[1].matvec(x))

    def _rmatvec(self, x):
        return self.args[1].rmatvec(self.args[0].rmatvec(x))

    def _matmat(self, x):
        return self.args[0].matmat(self.args[1].matmat(x))

    def _adjoint(self):
        A, B = self.args
        return B.H * A.H


class _ScaledLinearOperator(LinearOperator):
    def __init__(self, A, alpha):
        if not isinstance(A, LinearOperator):
            raise ValueError('LinearOperator expected as A')
        if not np.isscalar(alpha):
            raise ValueError('scalar expected as alpha')
        self.shape = A.shape
        self.dtype = _get_dtype([A], [type(alpha)])
        self.args = (A, alpha)

    def _matvec(self, x):
        return self.args[1] * self.args[0].matvec(x)

    def _rmatvec(self, x):
        return np.conj(self.args[1]) * self.args[0].rmatvec(x)

    def _matmat(self, x):
        return self.args[1] * self.args[0].matmat(x)

    def _adjoint(self):
        A, alpha = self.args
        return A.H * alpha


class _PowerLinearOperator(LinearOperator):
    def __init__(self, A, p):
        if not isinstance(A, LinearOperator):
            raise ValueError('LinearOperator expected as A')
        if A.shape[0] != A.shape[1]:
            raise ValueError('square LinearOperator expected, got %r' % A)
        if not isintlike(p):
            raise ValueError('integer expected as p')

        self.args = (A, p)
        self.dtype = _get_dtype([A])
        self.shape = A.shape

    def _power(self, fun, x):
        res = np.array(x, copy=True)
        for i in range(self.args[1]):
            res = fun(res)
        return res

    def _matvec(self, x):
        return self._power(self.args[0].matvec, x)

    def _rmatvec(self, x):
        return self._power(self.args[0].rmatvec, x)

    def _matmat(self, x):
        return self._power(self.args[0].matmat, x)

    def _adjoint(self):
        A, p = self.args
        return A.H ** p


class MatrixLinearOperator(LinearOperator):
    def __init__(self, A):
        self.A = A
        self.__adj = None
        self.args = (A,)
        self.dtype = A.dtype
        self.shape = A.shape

    def _matmat(self, X):
        return self.A.dot(X)

    def _matvec(self, x):
        return self.A.dot(x)

    def _rmatvec(self, x):
        return self._adjoint().matvec(x)

    def _adjoint(self):
        if self.__adj is None:
            self.__adj = _AdjointMatrixOperator(self)
        return self.__adj


class _AdjointMatrixOperator(MatrixLinearOperator):
    def __init__(self, adjoint):
        self.A = adjoint.A.T.conj()
        self.__adjoint = adjoint
        self.args = (adjoint,)
        self.shape = adjoint.shape[1], adjoint.shape[0]

    @property
    def dtype(self):
        return self.__adjoint.dtype

    def _adjoint(self):
        return self.__adjoint


class IdentityOperator(LinearOperator):
    def __init__(self, shape, dtype=None):
        self.dtype = dtype
        self.shape = shape

    def _matvec(self, x):
        return x

    def _rmatvec(self, x):
        return x

    def _matmat(self, x):
        return x

    def _adjoint(self):
        return self


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
