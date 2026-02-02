"""Abstract linear algebra library.

This module defines a class hierarchy that implements a kind of "lazy"
matrix representation, called the ``LinearOperator``. It can be used to do
linear algebra with extremely large sparse or structured matrices, without
representing those explicitly in memory. Such matrices can be added,
multiplied, transposed, etc.

As a motivating example, suppose you want have a matrix where almost all of
the elements have the value one. The standard sparse matrix representation
skips the storage of zeros, but not ones. By contrast, a LinearOperator is
able to represent such matrices efficiently. First, we need a compact way to
represent an all-ones matrix::

    >>> import numpy as np
    >>> from scipy.sparse.linalg._interface import LinearOperator
    >>> class Ones(LinearOperator):
    ...     def __init__(self, shape):
    ...         super().__init__(dtype=None, shape=shape)
    ...     def _matvec(self, x):
    ...         return np.repeat(x.sum(), self.shape[0])

Instances of this class emulate ``np.ones(shape)``, but using a constant
amount of storage, independent of ``shape``. The ``_matvec`` method specifies
how this linear operator multiplies with (operates on) a vector. We can now
add this operator to a sparse matrix that stores only offsets from one::

    >>> from scipy.sparse.linalg._interface import aslinearoperator
    >>> from scipy.sparse import csr_array
    >>> offsets = csr_array([[1, 0, 2], [0, -1, 0], [0, 0, 3]])
    >>> A = aslinearoperator(offsets) + Ones(offsets.shape)
    >>> A.dot([1, 2, 3])
    array([13,  4, 15])

The result is the same as that given by its dense, explicitly-stored
counterpart::

    >>> (np.ones(A.shape, A.dtype) + offsets.toarray()).dot([1, 2, 3])
    array([13,  4, 15])

Several algorithms in the ``scipy.sparse`` library are able to operate on
``LinearOperator`` instances.
"""

import os
import types
import warnings

import numpy as np

from scipy import sparse
from scipy.sparse import issparse
from scipy.sparse._sputils import isshape, isintlike, asmatrix, is_pydata_spmatrix

__all__ = ['LinearOperator', 'aslinearoperator']


class LinearOperator:
    """Common interface for performing matrix vector products.

    Many iterative methods (e.g. `cg`, `gmres`) do not need to know the
    individual entries of a matrix to solve a linear system ``A@x = b``.
    Such solvers only require the computation of matrix vector
    products, ``A@v``, where ``v`` is a dense vector.  This class serves as
    an abstract interface between iterative solvers and matrix-like
    objects.

    To construct a concrete `LinearOperator`, either pass appropriate
    callables to the constructor of this class, or subclass it.

    A subclass must implement either one of the methods ``_matvec``
    and ``_matmat``, and the attributes/properties ``shape`` (pair of
    integers, optionally with additional batch dimensions at the front)
    and ``dtype`` (may be None). It may call the ``__init__``
    on this class to have these attributes validated. Implementing
    ``_matvec`` automatically implements ``_matmat`` (using a naive
    algorithm) and vice-versa.

    Optionally, a subclass may implement ``_rmatvec`` or ``_adjoint``
    to implement the Hermitian adjoint (conjugate transpose). As with
    ``_matvec`` and ``_matmat``, implementing either ``_rmatvec`` or
    ``_adjoint`` implements the other automatically. Implementing
    ``_adjoint`` is preferable; ``_rmatvec`` is mostly there for
    backwards compatibility.

    The defined operator may have additional "batch" dimensions
    prepended to the core shape, to represent a batch of 2-D operators;
    see :ref:`linalg_batch` for details.

    Parameters
    ----------
    shape : tuple
        Matrix dimensions ``(..., M, N)``,
        where ``...`` represents any additional batch dimensions.
    matvec : callable f(v)
        Returns ``A @ v``, where ``v`` is a dense vector
        with shape ``(..., N)``.
    rmatvec : callable f(v)
        Returns ``A^H @ v``, where ``A^H`` is the conjugate transpose of ``A``,
        and ``v`` is a dense vector of shape ``(..., M)``.
    matmat : callable f(V)
        Returns ``A @ V``, where ``V`` is a dense matrix
        with dimensions ``(..., N, K)``.
    rmatmat : callable f(V)
        Returns ``A^H @ V``, where ``V`` is a dense matrix
        with dimensions ``(..., M, K)``.
    dtype : dtype
        Data type of the matrix or matrices.

    Attributes
    ----------
    args : tuple
        For linear operators describing products etc. of other linear
        operators, the operands of the binary operation.
    ndim : int
        Number of dimensions (greater than 2 in the case of batch dimensions).
    T : LinearOperator
        Transpose.
    H : LinearOperator
        Hermitian adjoint.
    
    Methods
    -------
    matvec
    matmat
    adjoint
    transpose
    rmatvec
    rmatmat
    dot
    rdot
    __mul__
    __matmul__
    __call__
    __add__
    __truediv__
    __rmul__
    __rmatmul__

    See Also
    --------
    aslinearoperator : Construct a `LinearOperator`.

    Notes
    -----
    The user-defined `matvec` function must properly handle the case
    where ``v`` has shape ``(..., N)``.

    It is highly recommended to explicitly specify the `dtype`, otherwise
    it is determined automatically at the cost of a single matvec application
    on ``int8`` zero vector using the promoted `dtype` of the output.
    It is assumed that `matmat`, `rmatvec`, and `rmatmat` would result in
    the same dtype of the output given an ``int8`` input as `matvec`.

    `LinearOperator` instances can also be multiplied, added with each
    other, and raised to integral powers, all lazily: the result of these
    operations
    is always a new, composite `LinearOperator`, that defers linear
    operations to the original operators and combines the results.

    More details regarding how to subclass a `LinearOperator` and several
    examples of concrete `LinearOperator` instances can be found in the
    external project `PyLops <https://pylops.readthedocs.io>`_.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse.linalg import LinearOperator
    >>> def mv(v):
    ...     return np.array([2*v[0], 3*v[1]])
    ...
    >>> A = LinearOperator((2,2), matvec=mv)
    >>> A
    <2x2 _CustomLinearOperator with dtype=int8>
    >>> A.matvec(np.ones(2))
    array([ 2.,  3.])
    >>> A @ np.ones(2)
    array([ 2.,  3.])

    """

    # Necessary for right matmul with numpy arrays.
    __array_ufunc__ = None

    # generic type compatibility with scipy-stubs
    __class_getitem__ = classmethod(types.GenericAlias)

    ndim: int

    def __new__(cls, *args, **kwargs):
        if cls is LinearOperator:
            # Operate as _CustomLinearOperator factory.
            return super().__new__(_CustomLinearOperator)
        else:
            obj = super().__new__(cls)

            if (type(obj)._matvec == LinearOperator._matvec
                    and type(obj)._matmat == LinearOperator._matmat):
                warnings.warn("LinearOperator subclass should implement"
                              " at least one of _matvec and _matmat.",
                              category=RuntimeWarning, stacklevel=2)

            return obj

    def __init__(self, dtype, shape):
        """Initialize this LinearOperator.

        To be called by subclasses. ``dtype`` may be None; ``shape`` should
        be convertible to a length >=2 tuple.
        """
        if dtype is not None:
            dtype = np.dtype(dtype)

        shape = tuple(shape)
        if not isshape(shape, check_nd=False) or len(shape) < 2:
            raise ValueError(f"invalid shape {shape!r} (must be at least 2-d)")

        self.dtype = dtype
        self.shape = shape
        self.ndim = len(shape)

    def _init_dtype(self):
        """Determine the dtype by executing `matvec` on an `int8` test vector.

        In `np.promote_types` hierarchy, the type `int8` is the smallest,
        so we call `matvec` on `int8` and use the promoted dtype of the output
        to set the default `dtype` of the `LinearOperator`.
        We assume that `matmat`, `rmatvec`, and `rmatmat` would result in
        the same dtype of the output given an `int8` input as `matvec`.

        Called from subclasses at the end of the __init__ routine.
        """
        if self.dtype is None:
            v = np.zeros(self.shape[-1], dtype=np.int8)
            try:
                matvec_v = np.asarray(self.matvec(v))
            except OverflowError:
                # Python large `int` promoted to `np.int64`or `np.int32`
                self.dtype = np.dtype(int)
            else:
                self.dtype = matvec_v.dtype

    def _matmat(self, X):
        """Default matrix-matrix multiplication handler.
        
        If ``self`` is a linear operator of shape ``(..., M, N)``,
        then this method will be called on a shape ``(..., N, K)`` array,
        and should return a shape ``(..., M, K)`` array.

        Falls back to `matvec`, so defining that will
        define matrix multiplication too (though in a very suboptimal way).
        """
        # Maintain backwards-compatibility for 1-D input
        if X.ndim == 1:
            X = X[np.newaxis]

        # NOTE: we can't use `_matvec` directly (for the unbatched case)
        # as we can't assume that user-defined `matvec` functions support batching.
        return np.stack(
            [self._matvec(X[..., :, i]) for i in range(X.shape[-1])],
            axis=-1
        )

    def _matvec(self, x):
        """Default matrix-vector multiplication handler.

        If ``self`` is a linear operator of shape ``(..., M, N)``,
        then this method will be called on a shape
        ``(..., N)`` array,
        and should return a shape ``(..., M)`` array.

        Falls back to `matmat`, so defining that
        will define matrix-vector multiplication as well.
        """
        return np.squeeze(self._matmat(x[..., np.newaxis]), axis=-1)

    def matvec(self, x):
        """Matrix-vector multiplication.

        Performs the operation ``A @ x`` where ``A`` is an ``M`` x ``N``
        linear operator (or batch of linear operators)
        and `x` is a row or column vector (or batch of such vectors).

        Parameters
        ----------
        x : {matrix, ndarray}
            An array with shape ``(..., N)`` representing a row vector
            (or batch of row vectors),
            or an array with shape ``(N, 1)`` representing a column vector.

            .. versionadded:: 1.18.0
                A ``FutureWarning`` is now emitted for column vector input of shape
                ``(N, 1)``. `matmat` can be called instead for identical behaviour
                on such input.

        Returns
        -------
        y : {matrix, ndarray}
            An array with shape ``(..., M)`` or ``(M, 1)`` depending
            on the type and shape of `x`.

        Notes
        -----
        This method wraps the user-specified ``matvec`` routine or overridden
        ``_matvec`` method to ensure that `y` has the correct shape and type.
        """

        x = np.asanyarray(x)

        *self_broadcast_dims, M, N = self.shape

        # TODO: deprecate `np.matrix` support
        if isinstance(x, np.matrix):
            y = self._matvec(x)
            if x.shape == (N, 1):
                y = y.reshape(M, 1)
            return asmatrix(y)

        x_broadcast_dims: tuple[int, ...] = ()
        row_vector: bool = False
        if column_vector := x.shape == (N, 1):
            msg = (
                "In the future, calling `matvec` on 'column vectors' of shape "
                "`(N, 1)` will be deprecated. Please call `matmat` instead "
                "for identical behaviour."
            )
            warnings.warn(
                msg, FutureWarning,
                skip_file_prefixes=(os.path.dirname(__file__),)
            )
            x = np.reshape(x, (N,))
        elif x.ndim >= 1 and (row_vector := x.shape[-1] == N):
            x_broadcast_dims = x.shape[:-1]

        if not (row_vector or column_vector):
            msg = (
                f"Dimension mismatch: `x` must have a shape ending in "
                f"`({N},)`, or shape `({N}, 1)`. Given shape: {x.shape}"
            )
            raise ValueError(msg)

        y = self._matvec(x)

        broadcasted_dims = np.broadcast_shapes(self_broadcast_dims, x_broadcast_dims)
        if row_vector:
            y = y.reshape(*broadcasted_dims, M)
        elif column_vector:
            y = y.reshape(*broadcasted_dims, M, 1)

        return y

    def rmatvec(self, x):
        """Adjoint matrix-vector multiplication.

        Performs the operation ``A^H @ x`` where ``A`` is an
        ``M`` x ``N`` linear operator (or batch of linear operators)
        and `x` is a row or column vector (or batch of such vectors).

        Parameters
        ----------
        x : {matrix, ndarray}
            An array with shape ``(..., M)`` representing a row vector
            (or batch of row vectors),
            or an array with shape ``(M, 1)`` representing a column vector.

            .. versionadded:: 1.18.0
                A ``FutureWarning`` is now emitted for column vector input of
                shape ``(M, 1)``.
                `rmatmat` can be called instead for identical behaviour
                on such input.

        Returns
        -------
        y : {matrix, ndarray}
            An array with shape ``(..., N)`` or ``(N, 1)`` depending
            on the shape of `x`.

        Notes
        -----
        This method wraps the user-specified ``rmatvec`` routine or overridden
        ``_rmatvec`` method to ensure that `y` has the correct shape and type.
        """

        x = np.asanyarray(x)

        *self_broadcast_dims, M, N = self.shape

        # TODO: deprecate `np.matrix` support
        if isinstance(x, np.matrix):
            y = self._rmatvec(x)
            if x.shape == (M, 1):
                y = y.reshape(N, 1)
            return asmatrix(y)

        x_broadcast_dims: tuple[int, ...] = ()
        row_vector: bool = False
        if column_vector := x.shape == (M, 1):
            msg = (
                "In the future, calling `rmatvec` on 'column vectors' of shape "
                "`(M, 1)` will be deprecated. Please call `rmatmat` instead "
                "for identical behaviour."
            )
            warnings.warn(
                msg, FutureWarning,
                skip_file_prefixes=(os.path.dirname(__file__),)
            )
            x = np.reshape(x, (M,))
        elif x.ndim >= 1 and (row_vector := x.shape[-1] == M):
            x_broadcast_dims = x.shape[:-1]
        if not (row_vector or column_vector):
            msg = (
                f"Dimension mismatch: `x` must have a shape ending in "
                f"`({M},)`, or shape `({M}, 1)`. Given shape: {x.shape}"
            )
            raise ValueError(msg)

        y = self._rmatvec(x)

        broadcasted_dims = np.broadcast_shapes(self_broadcast_dims, x_broadcast_dims)
        if row_vector:
            y = y.reshape(*broadcasted_dims, N)
        elif column_vector:
            y = y.reshape(*broadcasted_dims, N, 1)

        return y

    def _rmatvec(self, x):
        """Default implementation of `_rmatvec`.
        Defers to `_rmatmat` or `adjoint`."""
        if type(self)._adjoint == LinearOperator._adjoint:
            # _adjoint not overridden, prevent infinite recursion
            if (hasattr(self, "_rmatmat")
                    and type(self)._rmatmat != LinearOperator._rmatmat):
                # Try to use _rmatmat as a fallback
                return np.squeeze(self._rmatmat(x[..., np.newaxis]), axis=-1)
            raise NotImplementedError
        else:
            return self.H.matvec(x)

    def matmat(self, X):
        """Matrix-matrix multiplication.

        Performs the operation ``A @ X`` where ``A`` is an ``M`` x ``N``
        linear operator (or batch of linear operators)
        and ``X`` is a dense ``N`` x ``K`` matrix
        (or batch of dense matrices).

        Parameters
        ----------
        X : {matrix, ndarray}
            An array with shape ``(..., N, K)`` representing the dense matrix
            (or batch of dense matrices).

        Returns
        -------
        Y : {matrix, ndarray}
            An array with shape ``(..., M, K)`` depending on
            the type of `X`.

        Notes
        -----
        This method wraps any user-specified ``matmat`` routine or overridden
        ``_matmat`` method to ensure that `Y` has the correct type.
        """
        if not (issparse(X) or is_pydata_spmatrix(X)):
            X = np.asanyarray(X)

        if X.ndim < 2:
            raise ValueError(f'expected at least 2-d ndarray or matrix, not {X.ndim}-d')

        if X.shape[-2] != self.shape[-1]:
            raise ValueError(f'dimension mismatch: {self.shape}, {X.shape}')

        try:
            Y = self._matmat(X)
        except Exception as e:
            if issparse(X) or is_pydata_spmatrix(X):
                raise TypeError(
                    "Unable to multiply a LinearOperator with a sparse matrix."
                    " Wrap the matrix with `aslinearoperator` first."
                ) from e
            raise

        if isinstance(Y, np.matrix):
            Y = asmatrix(Y)

        return Y

    def rmatmat(self, X):
        """Adjoint matrix-matrix multiplication.

        Performs the operation ``A^H @ X`` where ``A`` is an ``M`` x ``N``
        linear operator (or batch of linear operators)
        and `X` is a dense ``M`` x ``K`` matrix
        (or batch of dense matrices).
        The default implementation defers to the adjoint.

        Parameters
        ----------
        X : {matrix, ndarray}
            An array with shape ``(..., M, K)`` representing the dense matrix
            (or batch of dense matrices).

        Returns
        -------
        Y : {matrix, ndarray}
            An array with shape ``(..., N, K)`` depending on the type of `X`.

        Notes
        -----
        This method wraps any user-specified ``rmatmat`` routine or overridden
        ``_rmatmat`` method to ensure that `Y` has the correct type.

        """
        if not (issparse(X) or is_pydata_spmatrix(X)):
            X = np.asanyarray(X)

        if X.ndim < 2:
            raise ValueError(f'expected at least 2-d ndarray or matrix, not {X.ndim}-d')

        if X.shape[-2] != self.shape[-2]:
            raise ValueError(f'dimension mismatch: {self.shape}, {X.shape}')

        try:
            Y = self._rmatmat(X)
        except Exception as e:
            if issparse(X) or is_pydata_spmatrix(X):
                raise TypeError(
                    "Unable to multiply a LinearOperator with a sparse matrix."
                    " Wrap the matrix in aslinearoperator() first."
                ) from e
            raise

        if isinstance(Y, np.matrix):
            Y = asmatrix(Y)
        return Y

    def _rmatmat(self, X):
        """Default implementation of `_rmatmat`; defers to `rmatvec` or `adjoint`."""
        if type(self)._adjoint == LinearOperator._adjoint:
            # Maintain backwards-compatibility for 1-D input
            if X.ndim == 1:
                X = X[np.newaxis]

            # NOTE: we can't use `_rmatvec` directly as we can't assume that
            # user-defined `rmatvec` functions support batching.
            return np.stack(
                [self._rmatvec(X[..., :, i]) for i in range(X.shape[-1])],
                axis=-1
            )
        else:
            return self.H.matmat(X)

    def __call__(self, x):
        """Apply this linear operator.
        
        Equivalent to `__matmul__`.
        """
        return self@x

    def __mul__(self, x):
        """Multiplication.
        
        Used by the ``*`` operator. Equivalent to `dot`.
        """
        return self.dot(x)

    def __truediv__(self, other):
        """Scalar Division.
        
        Returns a lazily scaled linear operator.
        """
        if not np.isscalar(other):
            raise ValueError("Can only divide a linear operator by a scalar.")

        return _ScaledLinearOperator(self, 1.0/other)

    def dot(self, x):
        """Multi-purpose multiplication method.

        Parameters
        ----------
        x : array_like or `LinearOperator` or scalar
            Array-like input will be interpreted as a 1-D row vector or
            2-D matrix (or batch of matrices)
            depending on its shape. See the Returns section for details.

        Returns
        -------
        Ax : array or `LinearOperator`
            - For `LinearOperator` input, operator composition is performed.

            - For scalar input, a lazily scaled operator is returned.

            - Otherwise, the input is expected to take the form of a dense
              1-D vector or 2-D matrix (or batch of matrices),
              interpreted as follows
              (where ``self`` is an ``M`` by ``N`` linear operator):
                
              - If `x` has shape ``(N,)``
                it is interpreted as a row vector
                and `matvec` is called.
              - If `x` has shape ``(..., N, K)`` for some
                integer ``K``, it is interpreted as a matrix
                (or batch of matrices if there are batch dimensions)
                and `matmat` is called.

        Notes
        -----
        To perform matrix-vector multiplication on batches of vectors,
        use `matvec`.

        For clarity, it is recommended to use the `matvec` or
        `matmat` methods directly instead of this method
        when interacting with dense vectors and matrices.

        See Also
        --------
        __mul__ : Equivalent method used by the ``*`` operator.
        __matmul__ :
            Method used by the ``@`` operator which rejects scalar
            input before calling this method.

        """
        if isinstance(x, LinearOperator):
            return _ProductLinearOperator(self, x)
        elif np.isscalar(x):
            return _ScaledLinearOperator(self, x)
        else:
            if not issparse(x) and not is_pydata_spmatrix(x):
                # Sparse matrices shouldn't be converted to numpy arrays.
                x = np.asarray(x)
                
            N = self.shape[-1]
    
            # treat 1-D input as a vector and >1-D input as a matrix, if the shape fits
            vector = x.shape == (N,)
            matrix = x.ndim >=2 and x.shape[-2] == N

            if not (vector or matrix):
                msg = (
                    f"Dimension mismatch: `x` must have shape `({N},)` "
                    f"or a shape ending in `({N}, K)` for some integer `K`. "
                    f"Given shape: {x.shape}"
                )
                raise ValueError(msg)
            
            if vector:
                return self.matvec(x)
            elif matrix:
                return self.matmat(x)

    def __matmul__(self, other):
        """Matrix Multiplication.
        
        Used by the ``@`` operator.
        Rejects scalar input.
        Otherwise, equivalent to `dot`.
        """
        if np.isscalar(other):
            raise ValueError("Scalar operands are not allowed, "
                             "use '*' instead")
        return self.__mul__(other)

    def __rmatmul__(self, other):
        """Matrix Multiplication from the right.
        
        Used by the ``@`` operator from the right.
        Rejects scalar input.
        Otherwise, equivalent to `rdot`.
        """
        if np.isscalar(other):
            raise ValueError("Scalar operands are not allowed, "
                             "use '*' instead")
        return self.__rmul__(other)

    def __rmul__(self, x):
        """Multiplication from the right.
        
        Used by the ``*`` operator from the right. Equivalent to `rdot`.
        """
        return self.rdot(x)

    def rdot(self, x):
        """Multi-purpose multiplication method from the right.

        .. note ::

            For complex data, this does not perform conjugation,
            returning ``xA`` rather than ``x A^H``.
            To calculate adjoint multiplication instead, use one of
            `rmatvec` or `rmatmat`, or take the adjoint first,
            like ``self.H.rdot(x)`` or ``x * self.H``.

        .. note ::
            
            Batched (>2-D) input to this function is unsupported.
            It is recommended to transpose data separately
            and then use forward operations like `matvec` and `matmat` directly.

        Parameters
        ----------
        x : array_like or `LinearOperator` or scalar
            Array-like input will be interpreted as a dense row vector
            or matrix depending on its shape.
            See the Returns section for details.

        Returns
        -------
        xA : array or `LinearOperator`
            - For `LinearOperator` input, operator composition is performed.

            - For scalar input, a lazily scaled operator is returned.

            - Otherwise, the input is expected to take the form of a dense
              vector or matrix, interpreted as follows
              (where ``self`` is an ``M`` by ``N`` linear operator):

              - If `x` has shape ``(M,)``
                it is interpreted as a row vector.
              - Otherwise, if `x` has shape ``(K, M)`` for some
                integer ``K``, it is interpreted as a matrix.

        See Also
        --------
        dot : Multi-purpose multiplication method from the left.
        __rmul__ :
            Equivalent method, used by the ``*`` operator from the right.
        __rmatmul__ :
            Method used by the ``@`` operator from the right
            which rejects scalar input before calling this method.
        """
        if isinstance(x, LinearOperator):
            return _ProductLinearOperator(x, self)
        elif np.isscalar(x):
            return _ScaledLinearOperator(self, x)
        else:
            if x.ndim > 2:
                msg = (
                    "Batched (>2-D) input is unsupported by `rdot`.\n"
                    "It is recommended to transpose data separately and then"
                    "use forward operations like `matvec` and `matmat` directly."
                )
                raise ValueError(msg)
            if not issparse(x) and not is_pydata_spmatrix(x):
                # Sparse matrices shouldn't be converted to numpy arrays.
                x = np.asarray(x)

            M = self.shape[-2]

            # treat 1-D input as a vector and 2-D input as a matrix, if the shape fits
            vector = x.shape == (M,)
            matrix = x.ndim == 2 and x.shape[-1] == M

            if not (vector or matrix):
                msg = (
                    f"Dimension mismatch: `x` must have shape `({M},)` "
                    f"or `(K, {M})` for some integer `K`. "
                    f"Given shape: {x.shape}"
                )
                raise ValueError(msg)
            
            # We use transpose instead of rmatvec/rmatmat to avoid
            # unnecessary complex conjugation if possible.
            if vector:
                return self.T.matvec(x.T).T
            elif matrix:
                return self.T.matmat(x.T).T

    def __pow__(self, p):
        if np.isscalar(p):
            return _PowerLinearOperator(self, p)
        else:
            return NotImplemented

    def __add__(self, x):
        """Linear operator addition.
        
        The input must be a `LinearOperator`.
        A lazily summed linear operator is returned.
        """
        if isinstance(x, LinearOperator):
            return _SumLinearOperator(self, x)
        else:
            return NotImplemented

    def __neg__(self):
        return _ScaledLinearOperator(self, -1)

    def __sub__(self, x):
        return self.__add__(-x)

    def __repr__(self):
        if self.dtype is None:
            dt = 'unspecified dtype'
        else:
            dt = 'dtype=' + str(self.dtype)

        shape = 'x'.join(str(dim) for dim in self.shape)
        return f'<{shape} {self.__class__.__name__} with {dt}>'

    def adjoint(self):
        """Hermitian adjoint.

        Returns the Hermitian adjoint of this linear operator,
        also known as the Hermitian
        conjugate or Hermitian transpose. For a complex matrix, the
        Hermitian adjoint is equal to the conjugate transpose.

        Returns
        -------
        `LinearOperator`
            Hermitian adjoint of self.
        
        See Also
        --------
        :attr:`~scipy.sparse.linalg.LinearOperator.H` : Equivalent attribute.
        """
        return self._adjoint()

    @property
    def H(self):
        """Hermitian adjoint.

        See Also
        --------
        scipy.sparse.linalg.LinearOperator.adjoint : Equivalent method.
        """
        return self.adjoint()

    def transpose(self):
        """Transpose.

        Returns
        -------
        `LinearOperator`
            Transpose of the linear operator.   
        
        See Also
        --------
        :attr:`~scipy.sparse.linalg.LinearOperator.T` : Equivalent attribute.
        """
        return self._transpose()

    @property
    def T(self):
        """Transpose.

        See Also
        --------
        scipy.sparse.linalg.LinearOperator.transpose : Equivalent method.
        """
        return self.transpose()

    def _adjoint(self):
        """Default implementation of `_adjoint`.
        Defers to adjoint functions, e.g. `_rmatvec` for `_matvec`."""
        return _AdjointLinearOperator(self)

    def _transpose(self):
        """Default implementation of `_transpose`.
        For `_matvec`, defers to `_rmatvec` + `np.conj`."""
        return _TransposedLinearOperator(self)


class _CustomLinearOperator(LinearOperator):
    """Linear operator defined in terms of user-specified operations."""

    def __init__(self, shape, matvec, rmatvec=None, matmat=None,
                 dtype=None, rmatmat=None):
        super().__init__(dtype, shape)

        self.args = ()

        self.__matvec_impl = matvec
        self.__rmatvec_impl = rmatvec
        self.__rmatmat_impl = rmatmat
        self.__matmat_impl = matmat

        self._init_dtype()

    def _matmat(self, X):
        if self.__matmat_impl is not None:
            return self.__matmat_impl(X)
        else:
            return super()._matmat(X)

    def _matvec(self, x):
        return self.__matvec_impl(x)

    def _rmatvec(self, x):
        func = self.__rmatvec_impl
        if func is None:
            raise NotImplementedError("rmatvec is not defined")
        return self.__rmatvec_impl(x)

    def _rmatmat(self, X):
        if self.__rmatmat_impl is not None:
            return self.__rmatmat_impl(X)
        else:
            return super()._rmatmat(X)

    def _adjoint(self):
        return _CustomLinearOperator(
            shape=(*self.shape[:-2], self.shape[-1], self.shape[-2]),
            matvec=self.__rmatvec_impl,
            rmatvec=self.__matvec_impl,
            matmat=self.__rmatmat_impl,
            rmatmat=self.__matmat_impl,
            dtype=self.dtype
        )


class _AdjointLinearOperator(LinearOperator):
    """Adjoint of arbitrary Linear Operator"""

    def __init__(self, A):
        shape = (*A.shape[:-2], A.shape[-1], A.shape[-2])
        super().__init__(dtype=A.dtype, shape=shape)
        self.A = A
        self.args = (A,)

    def _matvec(self, x):
        return self.A._rmatvec(x)

    def _rmatvec(self, x):
        return self.A._matvec(x)

    def _matmat(self, x):
        return self.A._rmatmat(x)

    def _rmatmat(self, x):
        return self.A._matmat(x)

class _TransposedLinearOperator(LinearOperator):
    """Transposition of arbitrary Linear Operator"""

    def __init__(self, A):
        shape = (*A.shape[:-2], A.shape[-1], A.shape[-2])
        super().__init__(dtype=A.dtype, shape=shape)
        self.A = A
        self.args = (A,)

    def _matvec(self, x):
        # NB. np.conj works also on sparse matrices
        return np.conj(self.A._rmatvec(np.conj(x)))

    def _rmatvec(self, x):
        return np.conj(self.A._matvec(np.conj(x)))

    def _matmat(self, x):
        # NB. np.conj works also on sparse matrices
        return np.conj(self.A._rmatmat(np.conj(x)))

    def _rmatmat(self, x):
        return np.conj(self.A._matmat(np.conj(x)))

def _get_dtype(operators, dtypes=None):
    """Returns the promoted dtype from input dtypes and operators."""
    if dtypes is None:
        dtypes = []
    for obj in operators:
        if obj is not None and hasattr(obj, 'dtype'):
            dtypes.append(obj.dtype)
    return np.result_type(*dtypes)


class _SumLinearOperator(LinearOperator):
    """Representing ``A + B``"""
    def __init__(self, A, B):
        if not isinstance(A, LinearOperator) or not isinstance(B, LinearOperator):
            raise ValueError('both operands have to be a LinearOperator')
        *A_broadcast_dims, A_M, A_N = A.shape
        *B_broadcast_dims, B_M, B_N = B.shape
        if (A_M, A_N) != (B_M, B_N):
            raise ValueError(f'cannot add {A} and {B}: shape mismatch')
        broadcasted_dims = np.broadcast_shapes(A_broadcast_dims, B_broadcast_dims)
        self.args = (A, B)
        super().__init__(_get_dtype([A, B]), (*broadcasted_dims, A_M, A_N))

    def _matvec(self, x):
        return self.args[0].matvec(x) + self.args[1].matvec(x)

    def _rmatvec(self, x):
        return self.args[0].rmatvec(x) + self.args[1].rmatvec(x)

    def _rmatmat(self, x):
        return self.args[0].rmatmat(x) + self.args[1].rmatmat(x)

    def _matmat(self, x):
        return self.args[0].matmat(x) + self.args[1].matmat(x)

    def _adjoint(self):
        A, B = self.args
        return A.H + B.H


class _ProductLinearOperator(LinearOperator):
    """Representing ``A @ B``"""
    def __init__(self, A, B):
        if not isinstance(A, LinearOperator) or not isinstance(B, LinearOperator):
            raise ValueError('both operands have to be a LinearOperator')
        *A_broadcast_dims, A_M, A_N = A.shape
        *B_broadcast_dims, B_M, B_N = B.shape
        if A_N != B_M:
            raise ValueError(f'cannot multiply {A} and {B}: shape mismatch')
        broadcasted_dims = np.broadcast_shapes(A_broadcast_dims, B_broadcast_dims)
        super().__init__(_get_dtype([A, B]), (*broadcasted_dims, A_M, B_N))
        self.args = (A, B)

    def _matvec(self, x):
        return self.args[0].matvec(self.args[1].matvec(x))

    def _rmatvec(self, x):
        return self.args[1].rmatvec(self.args[0].rmatvec(x))

    def _rmatmat(self, x):
        return self.args[1].rmatmat(self.args[0].rmatmat(x))

    def _matmat(self, x):
        return self.args[0].matmat(self.args[1].matmat(x))

    def _adjoint(self):
        A, B = self.args
        return B.H @ A.H


class _ScaledLinearOperator(LinearOperator):
    """Representing ``alpha * A``"""
    def __init__(self, A, alpha):
        if not isinstance(A, LinearOperator):
            raise ValueError('LinearOperator expected as A')
        if not np.isscalar(alpha):
            raise ValueError('scalar expected as alpha')
        if isinstance(A, _ScaledLinearOperator):
            A, alpha_original = A.args
            # Avoid in-place multiplication so that we don't accidentally mutate
            # the original prefactor.
            alpha = alpha * alpha_original

        dtype = _get_dtype([A], [type(alpha)])
        super().__init__(dtype, A.shape)
        self.args = (A, alpha)
        # Note: args[1] is alpha (a scalar), so use `*` below, not `@`

    def _matvec(self, x):
        return self.args[1] * self.args[0].matvec(x)

    def _rmatvec(self, x):
        return np.conj(self.args[1]) * self.args[0].rmatvec(x)

    def _rmatmat(self, x):
        return np.conj(self.args[1]) * self.args[0].rmatmat(x)

    def _matmat(self, x):
        return self.args[1] * self.args[0].matmat(x)

    def _adjoint(self):
        A, alpha = self.args
        return A.H * np.conj(alpha)


class _PowerLinearOperator(LinearOperator):
    """Representing ``A ** p``"""
    def __init__(self, A, p):
        if not isinstance(A, LinearOperator):
            raise ValueError('LinearOperator expected as A')
        if A.shape[-2] != A.shape[-1]:
            msg = f'square core-dimensions of LinearOperator expected, got {A!r}'
            raise ValueError(msg)
        if not isintlike(p) or p < 0:
            raise ValueError('non-negative integer expected as p')

        super().__init__(_get_dtype([A]), A.shape)
        self.args = (A, p)

    def _power(self, fun, x):
        res = np.array(x, copy=True)
        for i in range(self.args[1]):
            res = fun(res)
        return res

    def _matvec(self, x):
        return self._power(self.args[0].matvec, x)

    def _rmatvec(self, x):
        return self._power(self.args[0].rmatvec, x)

    def _rmatmat(self, x):
        return self._power(self.args[0].rmatmat, x)

    def _matmat(self, x):
        return self._power(self.args[0].matmat, x)

    def _adjoint(self):
        A, p = self.args
        return A.H ** p


class MatrixLinearOperator(LinearOperator):
    """Operator defined by a matrix `A` which implements ``@``."""
    def __init__(self, A):
        super().__init__(A.dtype, A.shape)
        self.A = A
        self.__adj = None
        self.args = (A,)

    def _matmat(self, X):
        return self.A @ X

    def _adjoint(self):
        if self.__adj is None:
            self.__adj = _AdjointMatrixOperator(self.A)
        return self.__adj


class _AdjointMatrixOperator(MatrixLinearOperator):
    """Representing ``A.H``, for `MatrixLinearOperator` `A`."""
    def __init__(self, A):
        if A.ndim > 2:
            if issparse(A):
                A_T = sparse.swapaxes(A, -1, -2)
            else:
                A_T = A.mT
        else:
            A_T = A.T
        self.A = A_T.conj()
        self.args = (A,)
        self.shape = (
            *A.shape[:-2], A.shape[-1], A.shape[-2]
        )
        self.ndim = A.ndim

    @property
    def dtype(self):
        return self.args[0].dtype

    def _adjoint(self):
        return MatrixLinearOperator(self.args[0])


class IdentityOperator(LinearOperator):
    def __init__(self, shape, dtype=None):
        super().__init__(dtype, shape)

    def _matvec(self, x):
        return x

    def _rmatvec(self, x):
        return x

    def _rmatmat(self, x):
        return x

    def _matmat(self, x):
        return x

    def _adjoint(self):
        return self


def aslinearoperator(A):
    """Return `A` as a `LinearOperator`.

    See the `LinearOperator` documentation for additional information.

    Parameters
    ----------
    A : object
        Object to convert to a `LinearOperator`. May be any one of the following types:

        - `numpy.ndarray`
        - `numpy.matrix`
        - `scipy.sparse` array
          (e.g. `~scipy.sparse.csr_array`, `~scipy.sparse.lil_array`, etc.)
        - `LinearOperator`
        - An object with ``.shape`` and ``.matvec`` attributes

    Returns
    -------
    B : LinearOperator
        A `LinearOperator` corresponding with `A`

    Notes
    -----
    If `A` has no ``.dtype`` attribute, the data type is determined by calling
    :func:`LinearOperator.matvec()` - set the ``.dtype`` attribute to prevent this
    call upon the linear operator creation.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse.linalg import aslinearoperator
    >>> M = np.array([[1,2,3],[4,5,6]], dtype=np.int32)
    >>> aslinearoperator(M)
    <2x3 MatrixLinearOperator with dtype=int32>
    """
    if isinstance(A, LinearOperator):
        return A

    elif isinstance(A, np.ndarray) or isinstance(A, np.matrix):
        A = np.atleast_2d(np.asarray(A))
        return MatrixLinearOperator(A)

    elif issparse(A) or is_pydata_spmatrix(A):
        return MatrixLinearOperator(A)

    else:
        if hasattr(A, 'shape') and hasattr(A, 'matvec'):
            rmatvec = None
            rmatmat = None
            dtype = None

            if hasattr(A, 'rmatvec'):
                rmatvec = A.rmatvec
            if hasattr(A, 'rmatmat'):
                rmatmat = A.rmatmat
            if hasattr(A, 'dtype'):
                dtype = A.dtype
            return LinearOperator(A.shape, A.matvec, rmatvec=rmatvec,
                                  rmatmat=rmatmat, dtype=dtype)

        else:
            raise TypeError('type not understood')
