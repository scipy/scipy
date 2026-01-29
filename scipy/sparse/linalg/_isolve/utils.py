__docformat__ = "restructuredtext en"

__all__ = []


from scipy.sparse.linalg._interface import (
    _xp_aslinearoperator, LinearOperator, IdentityOperator
)
from scipy._lib._array_api import (
    array_namespace, np_compat, xp_copy, xp_ravel, xp_result_type, _asarray,
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

    """
    A_ = A
    A, xp = _xp_aslinearoperator(A)

    if not nd_support and A.ndim > 2:
        raise ValueError(f"{A.ndim}-dimensional `A` is unsupported, expected 2-D.")

    # lazy = is_lazy_array(xp.empty(0))
    lazy = False
    if not lazy and (N := A.shape[-2]) != A.shape[-1]:
        raise ValueError(
            f"expected square matrix or batch of square matrices, "
            f"but got shape={(A.shape,)}"
        )

    xp_b = array_namespace(b)
    if xp_b != xp and not (xp.__name__ == 'sparse' and xp_b == np_compat):
        msg = (
            f"Mismatched array namespaces. "
            f"Namespace for A is {xp}, namespace for b is {xp_b}."
        )
        raise TypeError(msg)
    
    b = _asarray(b, subok=True, xp=xp)

    # maintain column vector backwards-compatibility in 2-D case
    column_vector = not lazy and b.ndim == 2 and b.shape[-2:] == (N, 1) 
    # otherwise treat as a row-vector
    row_vector = b.shape[-1] == N

    if not lazy and not (column_vector or row_vector):
        raise ValueError(f'shapes of A {A.shape} and b {b.shape} are '
                         'incompatible')

    if not xp.isdtype(b.dtype, ("real floating", "complex floating")):
        b = xp.astype(b, xp.float64)  # upcast non-FP types to float64

    if hasattr(A, 'dtype'):
        x_dtype = A.dtype
    else:
        x_dtype = A.matvec(b).dtype
    # XXX: does this match the previous coercion?
    x_dtype = xp_result_type(x_dtype, b.dtype, force_floating=True, xp=xp)

    b = xp.astype(b, x_dtype)  # make b the same type as x
    if column_vector:
        b = xp_ravel(b)

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
        M, xp_M = _xp_aslinearoperator(M)
        if xp_M != xp:
            msg = (
                f"Mismatched array namespaces. "
                f"Namespace for A is {xp}, namespace for M is {xp_M}."
            )
            raise TypeError(msg)
        if A.shape != M.shape:
            raise ValueError('matrix and preconditioner have different shapes')

    # set initial guess
    if x0 is None:
        x = xp.zeros((*M.shape[:-2], N), dtype=x_dtype)
    # XXX: proper error handling for `x0` of type `str` but not equal to `'Mb'`?
    elif isinstance(x0, str):
        if x0 == 'Mb':  # use nonzero initial guess ``M @ b``
            bCopy = xp_copy(b)
            x = M.matvec(bCopy)
    else:
        x = xp.asarray(x0, dtype=x_dtype, copy=True)
        
        # maintain column vector backwards-compatibility in 2-D case
        column_vector = x.ndim == 2 and x.shape[-2:] == (N, 1)
        # otherwise treat as a row-vector
        row_vector = x.shape[-1] == N
        
        if not (row_vector or column_vector):
            raise ValueError(f'shapes of A {A.shape} and '
                             f'x0 {x.shape} are incompatible')
        if column_vector:
            x = xp_ravel(x)

    return A, M, x, b, xp
