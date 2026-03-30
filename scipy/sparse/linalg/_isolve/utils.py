__docformat__ = "restructuredtext en"

__all__ = []

import numpy as np

from scipy.sparse.linalg._interface import (
    aslinearoperator, LinearOperator, IdentityOperator
)
from scipy._lib._array_api import (
    _asarray,
    array_namespace,
    xp_copy,
    xp_ravel,
    xp_result_type,
)

def id(x):
    return x

def make_system(A, M, x0, b, nd_support=False):
    """Make a linear system Ax=b

    Parameters
    ----------
    A : LinearOperator
        sparse or dense matrix (or any valid input to aslinearoperator)
    M : {LinearOperator, None}
        preconditioner
        sparse or dense matrix (or any valid input to aslinearoperator)
    x0 : {array_like, str, None}
        initial guess to iterative method.
        ``x0 = 'Mb'`` means using the nonzero initial guess ``M @ b``.
        Default is `None`, which means using the zero initial guess.
    b : array_like
        right hand side
    nd_support: bool, optional
        Whether or not the calling algorithm supports n-dimensional
        batchesfor input `A` and `b`. Default: ``False``.

    Returns
    -------
    (A, M, x, b, xp)
        A : LinearOperator
            matrix of the linear system
        M : LinearOperator
            preconditioner
        x : rank 1 ndarray
            initial guess
        b : rank 1 ndarray
            right hand side
        xp : compatible array namespace
        batched : bool
            True if input is batched, i.e. ``A.ndim > 2 or b.ndim > 1``.

    """
    A_ = A
    A = aslinearoperator(A)

    if (N := A.shape[-2]) != A.shape[-1]:
        raise ValueError(
            f"expected square matrix or batch of square matrices, "
            f"but got shape={(A.shape,)}"
        )

    xp = A._xp
    if (xp_b := array_namespace(b)) != xp:
        msg = (
            f"Mismatched array namespaces. "
            f"Namespace for A is {xp}, namespace for b is {xp_b}."
        )
        raise TypeError(msg)

    b = _asarray(b, subok=True, xp=xp)

    # maintain column vector backwards-compatibility in 2-D case
    column_vector = b.ndim == 2 and b.shape[-2:] == (N, 1) 
    # otherwise treat as a row-vector
    row_vector = b.shape[-1] == N

    if not (column_vector or row_vector):
        raise ValueError(f'shapes of A {A.shape} and b {b.shape} are '
                         'incompatible')
    if column_vector:
        b = xp_ravel(b, xp=xp)
    
    # NOTE: unbatched column vectors are ravelled and hence are 1-D.
    batched = A.ndim > 2 or b.ndim > 1
    if batched and not nd_support:
        raise ValueError(
            f"{A.ndim}-dimensional `A` and `{b.ndim}-dimensional `b` "
            f"is unsupported, expected 2-D `A` and 1-D `b`."
        )

    dtype = xp_result_type(A.dtype, b.dtype, force_floating=True, xp=xp)

    b = xp.astype(b, dtype)  # make b the same type as x

    # process preconditioner
    if M is None:
        if hasattr(A_,'psolve'):
            psolve = A_.psolve
        else:
            psolve = id
        if hasattr(A_,'rpsolve'):
            rpsolve = A_.rpsolve
        else:
            rpsolve = id
        if psolve is id and rpsolve is id:
            M = IdentityOperator(shape=A.shape, dtype=A.dtype, xp=xp)
        else:
            M = LinearOperator(A.shape, matvec=psolve, rmatvec=rpsolve,
                               dtype=A.dtype, xp=xp)
    else:
        M = aslinearoperator(M)
        if (xp_M := M._xp) != xp:
            msg = (
                f"Mismatched array namespaces. "
                f"Namespace for A is {xp}, namespace for M is {xp_M}."
            )
            raise TypeError(msg)
        if A.shape != M.shape:
            raise ValueError('matrix and preconditioner have different shapes')

    # set initial guess
    if x0 is None:
        x = xp.zeros((*M.shape[:-2], N), dtype=dtype)
    elif isinstance(x0, str):
        if x0 == 'Mb':  # use nonzero initial guess ``M @ b``
            bCopy = xp_copy(b, xp=xp)
            x = M.matvec(bCopy)
        else:
            raise ValueError(f"invalid input for x0: {x0}")
    else:
        x = xp.asarray(x0, dtype=dtype, copy=True)

        # maintain column vector backwards-compatibility in 2-D case
        column_vector = x.ndim == 2 and x.shape[-2:] == (N, 1)
        # otherwise treat as a row-vector
        row_vector = x.shape[-1] == N

        if not (row_vector or column_vector):
            raise ValueError(f'shapes of A {A.shape} and '
                             f'x0 {x.shape} are incompatible')
        if column_vector:
            x = xp_ravel(x, xp=xp)

    result_shape = np.broadcast_shapes(A.shape[:-1], b.shape)
    def _broadcast_if_needed(arr):
        if arr.shape != result_shape:
            arr = xp.broadcast_to(arr, shape=result_shape)
            arr = xp_copy(arr, xp=xp) # avoid read-only arrays
        return arr
    
    x, b = map(_broadcast_if_needed, [x, b])

    return A, M, x, b, xp, batched
