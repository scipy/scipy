"""
Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG).

References
----------
.. [1] A. V. Knyazev (2001),
       Toward the Optimal Preconditioned Eigensolver: Locally Optimal
       Block Preconditioned Conjugate Gradient Method.
       SIAM Journal on Scientific Computing 23, no. 2,
       pp. 517-541. :doi:`10.1137/S1064827500366124`

.. [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov (2007),
       Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX)
       in hypre and PETSc.  :arxiv:`0705.2626`

.. [3] A. V. Knyazev's C and MATLAB implementations:
       https://github.com/lobpcg/blopex
"""

import warnings
import numpy as np
from scipy.linalg import (inv, eigh, cho_factor, cho_solve,
                          cholesky, LinAlgError)
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import issparse

__all__ = ["lobpcg"]


def _report_nonhermitian(M, name):
    """
    Report if `M` is not a Hermitian matrix given its type.
    """
    from scipy.linalg import norm

    md = M - M.T.conj()
    nmd = norm(md, 1)
    tol = 10 * np.finfo(M.dtype).eps
    tol = max(tol, tol * norm(M, 1))
    if nmd > tol:
        warnings.warn(
              f"Matrix {name} of the type {M.dtype} is not Hermitian: "
              f"condition: {nmd} < {tol} fails.",
              UserWarning, stacklevel=4
         )

def _as2d(ar):
    """
    If the input array is 2D return it, if it is 1D, append a dimension,
    making it a column vector.
    """
    if ar.ndim == 2:
        return ar
    else:  # Assume 1!
        aux = np.asarray(ar).reshape((ar.shape[0], 1))
        return aux


def _makeMatMat(m):
    if m is None:
        return None
    elif callable(m):
        return lambda v: m(v)
    else:
        return lambda v: m @ v


def _matmul_inplace(x, y, verbosityLevel=0):
    """Perform 'np.matmul' in-place if possible.

    If some sufficient conditions for inplace matmul are met, do so.
    Otherwise try inplace update and fall back to overwrite if that fails.
    """
    if x.flags["CARRAY"] and x.shape[1] == y.shape[1] and x.dtype == y.dtype:
        # conditions where we can guarantee that inplace updates will work;
        # i.e. x is not a view/slice, x & y have compatible dtypes, and the
        # shape of the result of x @ y matches the shape of x.
        np.matmul(x, y, out=x)
    else:
        # ideally, we'd have an exhaustive list of conditions above when
        # inplace updates are possible; since we don't, we opportunistically
        # try if it works, and fall back to overwriting if necessary
        try:
            np.matmul(x, y, out=x)
        except Exception:
            if verbosityLevel:
                warnings.warn(
                    "Inplace update of x = x @ y failed, "
                    "x needs to be overwritten.",
                    UserWarning, stacklevel=3
                )
            x = x @ y
    return x


def _applyConstraints(blockVectorV, factYBY, blockVectorBY, blockVectorY):
    """Changes blockVectorV in-place."""
    YBV = blockVectorBY.T.conj() @ blockVectorV
    tmp = cho_solve(factYBY, YBV)
    blockVectorV -= blockVectorY @ tmp


def _b_orthonormalize(B, blockVectorV, blockVectorBV=None,
                      verbosityLevel=0):
    """in-place B-orthonormalize the given block vector using Cholesky."""
    if blockVectorBV is None:
        if B is None:
            blockVectorBV = blockVectorV
        else:
            try:
                blockVectorBV = B(blockVectorV)
            except Exception as e:
                if verbosityLevel:
                    warnings.warn(
                        f"Secondary MatMul call failed with error\n"
                        f"{e}\n",
                        UserWarning, stacklevel=3
                    )
                    return None, None, None
            if blockVectorBV.shape != blockVectorV.shape:
                raise ValueError(
                    f"The shape {blockVectorV.shape} "
                    f"of the orthogonalized matrix not preserved\n"
                    f"and changed to {blockVectorBV.shape} "
                    f"after multiplying by the secondary matrix.\n"
                )

    VBV = blockVectorV.T.conj() @ blockVectorBV
    try:
        # VBV is a Cholesky factor from now on...
        VBV = cholesky(VBV, overwrite_a=True)
        VBV = inv(VBV, overwrite_a=True)
        blockVectorV = _matmul_inplace(
            blockVectorV, VBV,
            verbosityLevel=verbosityLevel
        )
        if B is not None:
            blockVectorBV = _matmul_inplace(
                blockVectorBV, VBV,
                verbosityLevel=verbosityLevel
            )
        return blockVectorV, blockVectorBV, VBV
    except LinAlgError:
        if verbosityLevel:
            warnings.warn(
                "Cholesky has failed.",
                UserWarning, stacklevel=3
            )
        return None, None, None


def _get_indx(_lambda, num, largest):
    """Get `num` indices into `_lambda` depending on `largest` option."""
    ii = np.argsort(_lambda)
    if largest:
        ii = ii[:-num - 1:-1]
    else:
        ii = ii[:num]

    return ii


def _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel):
    if verbosityLevel:
        _report_nonhermitian(gramA, "gramA")
        _report_nonhermitian(gramB, "gramB")


def _check_X_Y(blockVectorX, blockVectorY):
    """Validate the initial approximations X and the constraints Y."""
    if blockVectorY is not None and len(blockVectorY.shape) != 2:
        warnings.warn(
            f"Expected rank-2 array for argument Y, instead got "
            f"{len(blockVectorY.shape)}, "
            f"so ignore it and use no constraints.",
            UserWarning, stacklevel=3
        )
        blockVectorY = None

    # Block size.
    if blockVectorX is None:
        raise ValueError("The mandatory initial matrix X cannot be None")
    if len(blockVectorX.shape) != 2:
        raise ValueError("expected rank-2 array for argument X")

    # Data type of iterates, determined by X, must be inexact
    if not np.issubdtype(blockVectorX.dtype, np.inexact):
        warnings.warn(
            f"Data type for argument X is {blockVectorX.dtype}, "
            f"which is not inexact, so casted to np.float32.",
            UserWarning, stacklevel=3
        )
        blockVectorX = np.asarray(blockVectorX, dtype=np.float32)
    return blockVectorX, blockVectorY


def _apply_checked(op, blockVector, name, which):
    """Apply `op` to `blockVector` and check that its shape is preserved."""
    result = op(blockVector)
    if result.shape != blockVector.shape:
        raise ValueError(
            f"The shape {blockVector.shape} "
            f"of the {name} not preserved\n"
            f"and changed to {result.shape} "
            f"after multiplying by the {which} matrix.\n"
        )
    return result


def _residuals(B, blockVectorX, blockVectorAX, blockVectorBX, _lambda):
    """Return the block of residuals and their column norms."""
    if B is not None:
        aux = blockVectorBX * _lambda[np.newaxis, :]
    else:
        aux = blockVectorX * _lambda[np.newaxis, :]

    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
    residualNorms = np.sqrt(np.abs(aux))
    return blockVectorR, residualNorms


def _print_problem_info(B, M, Y, n, sizeX, sizeY):
    aux = "Solving "
    if B is None:
        aux += "standard"
    else:
        aux += "generalized"
    aux += " eigenvalue problem with"
    if M is None:
        aux += "out"
    aux += " preconditioning\n\n"
    aux += f"matrix size {n}\n"
    aux += f"block size {sizeX}\n\n"
    if Y is None:
        aux += "No constraints\n\n"
    elif sizeY > 1:
        aux += f"{sizeY} constraints\n\n"
    else:
        aux += f"{sizeY} constraint\n\n"
    print(aux)


def _as_dense(A, n, which):
    """Convert the matrix `A` of any supported type to a dense array."""
    try:
        if isinstance(A, LinearOperator):
            A = A(np.eye(n, dtype=int))
        elif callable(A):
            A = A(np.eye(n, dtype=int))
            if A.shape != (n, n):
                raise ValueError(
                    f"The shape {A.shape} of the {which} matrix\n"
                    f"defined by a callable object is wrong.\n"
                    f"Expected {(n, n)}."
                )
        elif issparse(A):
            A = A.toarray()
        else:
            A = np.asarray(A)
    except Exception as e:
        raise Exception(
            f"{which.capitalize()} MatMul call failed with error\n"
            f"{e}\n")
    return A


def _dense_eigh(A, B, n, sizeX, largest):
    """Solve the eigenproblem with a dense eigensolver instead of LOBPCG."""
    # Define the closed range of indices of eigenvalues to return.
    if largest:
        eigvals = (n - sizeX, n - 1)
    else:
        eigvals = (0, sizeX - 1)

    A = _as_dense(A, n, "primary")
    if B is not None:
        B = _as_dense(B, n, "secondary")

    try:
        vals, vecs = eigh(A,
                          B,
                          subset_by_index=eigvals,
                          check_finite=False)
        if largest:
            # Reverse order to be compatible with eigs() in 'LM' mode.
            vals = vals[::-1]
            vecs = vecs[:, ::-1]

        return vals, vecs
    except Exception as e:
        raise Exception(
            f"Dense eigensolver failed with error\n"
            f"{e}\n"
        )


def _split_eig_block(eigBlockVector, sizeX, currentBlockSize, restart):
    """Split Ritz coefficients into the parts acting on X, R and P."""
    eigBlockVectorX = eigBlockVector[:sizeX]
    if restart:
        return eigBlockVectorX, eigBlockVector[sizeX:], None
    eigBlockVectorR = eigBlockVector[sizeX:sizeX + currentBlockSize]
    eigBlockVectorP = eigBlockVector[sizeX + currentBlockSize:]
    return eigBlockVectorX, eigBlockVectorR, eigBlockVectorP


def _new_direction(blockVectorR, blockVectorP,
                   eigBlockVectorR, eigBlockVectorP, restart):
    """Combine the R (and, unless restarting, P) blocks into new directions."""
    pp = np.dot(blockVectorR, eigBlockVectorR)
    if not restart:
        pp += np.dot(blockVectorP, eigBlockVectorP)
    return pp


def _b_orthogonalize_to_X(B, blockVectorX, blockVectorBX, blockVectorR):
    """B-orthogonalize the block vector R to the B-orthonormal X."""
    if B is not None:
        return blockVectorR - (
            blockVectorX @
            (blockVectorBX.T.conj() @ blockVectorR)
        )
    return blockVectorR - (
        blockVectorX @
        (blockVectorX.T.conj() @ blockVectorR)
    )


def _explicit_gram_flag(explicitGramFlag, residualNorms,
                        blockVectorR, blockVectorAR):
    """Decide whether the Gram matrices must be computed explicitly."""
    if blockVectorAR.dtype == "float32":
        myeps = 1
    else:
        myeps = np.sqrt(np.finfo(blockVectorR.dtype).eps)

    if residualNorms.max() > myeps and not explicitGramFlag:
        return False
    # Once explicitGramFlag, forever explicitGramFlag.
    return True


def _rayleigh_ritz(blockVectorX, blockVectorAX, blockVectorBX,
                   blockVectorR, blockVectorAR, blockVectorBR,
                   blockVectorP, blockVectorAP, blockVectorBP,
                   _lambda, sizeX, currentBlockSize, explicitGramFlag,
                   restart, iterationNumber, verbosityLevel):
    """Rayleigh-Ritz procedure on the trial subspace spanned by X, R and P.

    The directions P are dropped when restarting, or when the eigensolver
    fails on the full trial subspace. Returns ``(_lambda, eigBlockVector,
    restart)``, or None if the eigensolver fails on the reduced subspace.
    """
    # Compute symmetric Gram matrices:
    # Common submatrices:
    gramXAR = np.dot(blockVectorX.T.conj(), blockVectorAR)
    gramRAR = np.dot(blockVectorR.T.conj(), blockVectorAR)

    gramDtype = blockVectorAR.dtype
    if explicitGramFlag:
        gramRAR = (gramRAR + gramRAR.T.conj()) / 2
        gramXAX = np.dot(blockVectorX.T.conj(), blockVectorAX)
        gramXAX = (gramXAX + gramXAX.T.conj()) / 2
        gramXBX = np.dot(blockVectorX.T.conj(), blockVectorBX)
        gramRBR = np.dot(blockVectorR.T.conj(), blockVectorBR)
        gramXBR = np.dot(blockVectorX.T.conj(), blockVectorBR)
    else:
        gramXAX = np.diag(_lambda).astype(gramDtype)
        gramXBX = np.eye(sizeX, dtype=gramDtype)
        gramRBR = np.eye(currentBlockSize, dtype=gramDtype)
        gramXBR = np.zeros((sizeX, currentBlockSize), dtype=gramDtype)

    if not restart:
        gramXAP = np.dot(blockVectorX.T.conj(), blockVectorAP)
        gramRAP = np.dot(blockVectorR.T.conj(), blockVectorAP)
        gramPAP = np.dot(blockVectorP.T.conj(), blockVectorAP)
        gramXBP = np.dot(blockVectorX.T.conj(), blockVectorBP)
        gramRBP = np.dot(blockVectorR.T.conj(), blockVectorBP)
        if explicitGramFlag:
            gramPAP = (gramPAP + gramPAP.T.conj()) / 2
            gramPBP = np.dot(blockVectorP.T.conj(), blockVectorBP)
        else:
            gramPBP = np.eye(currentBlockSize, dtype=gramDtype)

        gramA = np.block(
            [
                [gramXAX, gramXAR, gramXAP],
                [gramXAR.T.conj(), gramRAR, gramRAP],
                [gramXAP.T.conj(), gramRAP.T.conj(), gramPAP],
            ]
        )
        gramB = np.block(
            [
                [gramXBX, gramXBR, gramXBP],
                [gramXBR.T.conj(), gramRBR, gramRBP],
                [gramXBP.T.conj(), gramRBP.T.conj(), gramPBP],
            ]
        )

        _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel)

        try:
            _lambda, eigBlockVector = eigh(gramA,
                                           gramB,
                                           check_finite=False)
            return _lambda, eigBlockVector, restart
        except LinAlgError as e:
            # raise ValueError("eigh failed in lobpcg iterations") from e
            if verbosityLevel:
                warnings.warn(
                    f"eigh failed at iteration {iterationNumber} \n"
                    f"with error {e} causing a restart.\n",
                    UserWarning, stacklevel=3
                )
            # try again after dropping the direction vectors P from RR
            restart = True

    gramA = np.block([[gramXAX, gramXAR], [gramXAR.T.conj(), gramRAR]])
    gramB = np.block([[gramXBX, gramXBR], [gramXBR.T.conj(), gramRBR]])

    _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel)

    try:
        _lambda, eigBlockVector = eigh(gramA,
                                       gramB,
                                       check_finite=False)
    except LinAlgError as e:
        # raise ValueError("eigh failed in lobpcg iterations") from e
        warnings.warn(
            f"eigh failed at iteration {iterationNumber} with error\n"
            f"{e}\n",
            UserWarning, stacklevel=3
        )
        return None
    return _lambda, eigBlockVector, restart


def _postprocess(A, B, blockVectorX, sizeX, largest, verbosityLevel):
    """Final "exact" Rayleigh-Ritz on X to make the eigenvectors
    "exactly" orthonormalized.

    Returns the eigenvalues, the eigenvectors and their residual norms.
    """
    blockVectorAX = _apply_checked(A, blockVectorX,
                                   "postprocessing iterate", "primary")
    gramXAX = np.dot(blockVectorX.T.conj(), blockVectorAX)

    blockVectorBX = blockVectorX
    if B is not None:
        blockVectorBX = _apply_checked(B, blockVectorX,
                                       "postprocessing iterate", "secondary")

    gramXBX = np.dot(blockVectorX.T.conj(), blockVectorBX)
    _handle_gramA_gramB_verbosity(gramXAX, gramXBX, verbosityLevel)
    gramXAX = (gramXAX + gramXAX.T.conj()) / 2
    gramXBX = (gramXBX + gramXBX.T.conj()) / 2
    try:
        _lambda, eigBlockVector = eigh(gramXAX,
                                       gramXBX,
                                       check_finite=False)
    except LinAlgError as e:
        raise ValueError("eigh has failed in lobpcg postprocessing") from e

    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    eigBlockVector = np.asarray(eigBlockVector[:, ii])

    blockVectorX = np.dot(blockVectorX, eigBlockVector)
    blockVectorAX = np.dot(blockVectorAX, eigBlockVector)

    if B is not None:
        blockVectorBX = np.dot(blockVectorBX, eigBlockVector)

    _, residualNorms = _residuals(
        B, blockVectorX, blockVectorAX, blockVectorBX, _lambda)
    return _lambda, blockVectorX, residualNorms


def _finalize_history(history, bestIterationNumber, values):
    """Record the final `values`, truncate and split the history into a list."""
    history[bestIterationNumber + 1, :] = values
    history = history[: bestIterationNumber + 2, :]
    history = np.vsplit(history, np.shape(history)[0])
    return [np.squeeze(i) for i in history]


def lobpcg(
    A,
    X,
    B=None,
    M=None,
    Y=None,
    tol=None,
    maxiter=None,
    largest=True,
    verbosityLevel=0,
    retLambdaHistory=False,
    retResidualNormsHistory=False,
    restartControl=20,
):
    """Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG).

    LOBPCG is a preconditioned eigensolver for large real symmetric and complex
    Hermitian definite generalized eigenproblems.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator, callable object}
        The Hermitian linear operator of the problem, usually given by a
        sparse matrix.  Often called the "stiffness matrix".
    X : ndarray, float32 or float64
        Initial approximation to the ``k`` eigenvectors (non-sparse).
        If `A` has ``shape=(n,n)`` then `X` must have ``shape=(n,k)``.
    B : {sparse matrix, ndarray, LinearOperator, callable object}
        Optional. By default ``B = None``, which is equivalent to identity.
        The right hand side operator in a generalized eigenproblem if present.
        Often called the "mass matrix". Must be Hermitian positive definite.
    M : {sparse matrix, ndarray, LinearOperator, callable object}
        Optional. By default ``M = None``, which is equivalent to identity.
        Preconditioner aiming to accelerate convergence.
    Y : ndarray, float32 or float64, default: None
        An ``n-by-sizeY`` ndarray of constraints with ``sizeY < n``.
        The iterations will be performed in the ``B``-orthogonal complement
        of the column-space of `Y`. `Y` must be full rank if present.
    tol : scalar, optional
        The default is ``tol=n*sqrt(eps)``.
        Solver tolerance for the stopping criterion.
    maxiter : int, default: 20
        Maximum number of iterations.
    largest : bool, default: True
        When True, solve for the largest eigenvalues, otherwise the smallest.
    verbosityLevel : int, optional
        By default ``verbosityLevel=0`` no output.
        Controls the solver standard/screen output.
    retLambdaHistory : bool, default: False
        Whether to return iterative eigenvalue history.
    retResidualNormsHistory : bool, default: False
        Whether to return iterative history of residual norms.
    restartControl : int, optional
        Iterations restart if the residuals jump ``2**restartControl`` times
        compared to the smallest recorded in ``retResidualNormsHistory``.
        The default is ``restartControl=20``, making the restarts rare for
        backward compatibility.

    Returns
    -------
    lambda : ndarray of the shape ``(k, )``.
        Array of ``k`` approximate eigenvalues.
    v : ndarray of the same shape as ``X.shape``.
        An array of ``k`` approximate eigenvectors.
    lambdaHistory : ndarray, optional.
        The eigenvalue history, if `retLambdaHistory` is ``True``.
    ResidualNormsHistory : ndarray, optional.
        The history of residual norms, if `retResidualNormsHistory`
        is ``True``.

    Notes
    -----
    The iterative loop runs ``maxit=maxiter`` (20 if ``maxit=None``)
    iterations at most and finishes earlier if the tolerance is met.
    Breaking backward compatibility with the previous version, LOBPCG
    now returns the block of iterative vectors with the best accuracy rather
    than the last one iterated, as a cure for possible divergence.

    If ``X.dtype == np.float32`` and user-provided operations/multiplications
    by `A`, `B`, and `M` all preserve the ``np.float32`` data type,
    all the calculations and the output are in ``np.float32``.

    The size of the iteration history output equals to the number of the best
    (limited by `maxit`) iterations plus 3: initial, final, and postprocessing.

    If both `retLambdaHistory` and `retResidualNormsHistory` are ``True``,
    the return tuple has the following format
    ``(lambda, V, lambda history, residual norms history)``.

    In the following ``n`` denotes the matrix size and ``k`` the number
    of required eigenvalues (smallest or largest).

    The LOBPCG code internally solves eigenproblems of the size ``3k`` on every
    iteration by calling the dense eigensolver `eigh`, so if ``k`` is not
    small enough compared to ``n``, it makes no sense to call the LOBPCG code.
    Moreover, if one calls the LOBPCG algorithm for ``5k > n``, it would likely
    break internally, so the code calls the standard function `eigh` instead.
    It is not that ``n`` should be large for the LOBPCG to work, but rather the
    ratio ``n / k`` should be large. It you call LOBPCG with ``k=1``
    and ``n=10``, it works though ``n`` is small. The method is intended
    for extremely large ``n / k``.

    The convergence speed depends basically on three factors:

    1. Quality of the initial approximations `X` to the seeking eigenvectors.
       Randomly distributed around the origin vectors work well if no better
       choice is known.

    2. Relative separation of the desired eigenvalues from the rest
       of the eigenvalues. One can vary ``k`` to improve the separation.

    3. Proper preconditioning to shrink the spectral spread.
       For example, a rod vibration test problem (under tests
       directory) is ill-conditioned for large ``n``, so convergence will be
       slow, unless efficient preconditioning is used. For this specific
       problem, a good simple preconditioner function would be a linear solve
       for `A`, which is easy to code since `A` is tridiagonal.

    References
    ----------
    .. [1] A. V. Knyazev (2001),
           Toward the Optimal Preconditioned Eigensolver: Locally Optimal
           Block Preconditioned Conjugate Gradient Method.
           SIAM Journal on Scientific Computing 23, no. 2,
           pp. 517-541. :doi:`10.1137/S1064827500366124`

    .. [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov
           (2007), Block Locally Optimal Preconditioned Eigenvalue Xolvers
           (BLOPEX) in hypre and PETSc. :arxiv:`0705.2626`

    .. [3] A. V. Knyazev's C and MATLAB implementations:
           https://github.com/lobpcg/blopex

    Examples
    --------
    Our first example is minimalistic - find the largest eigenvalue of
    a diagonal matrix by solving the non-generalized eigenvalue problem
    ``A x = lambda x`` without constraints or preconditioning.

    >>> import numpy as np
    >>> from scipy.sparse import diags_array
    >>> from scipy.sparse.linalg import LinearOperator, aslinearoperator
    >>> from scipy.sparse.linalg import lobpcg

    The square matrix size is

    >>> n = 100

    and its diagonal entries are 1, ..., 100 defined by

    >>> vals = np.arange(1, n + 1).astype(np.int16)

    The first mandatory input parameter in this test is
    the sparse diagonal matrix `A`
    of the eigenvalue problem ``A x = lambda x`` to solve.

    >>> A = diags_array(vals, offsets=0, shape=(n, n), dtype=None)
    >>> A.toarray()
    array([[  1,   0,   0, ...,   0,   0,   0],
           [  0,   2,   0, ...,   0,   0,   0],
           [  0,   0,   3, ...,   0,   0,   0],
           ...,
           [  0,   0,   0, ...,  98,   0,   0],
           [  0,   0,   0, ...,   0,  99,   0],
           [  0,   0,   0, ...,   0,   0, 100]], shape=(100, 100), dtype=int16)

    The second mandatory input parameter `X` is a 2D array with the
    row dimension determining the number of requested eigenvalues.
    `X` is an initial guess for targeted eigenvectors.
    `X` must have linearly independent columns.
    If no initial approximations available, randomly oriented vectors
    commonly work best, e.g., with components normally distributed
    around zero or uniformly distributed on the interval [-1 1].
    Setting the initial approximations to dtype ``np.float32``
    forces all iterative values to dtype ``np.float32`` speeding up
    the run while still allowing accurate eigenvalue computations.

    >>> k = 1
    >>> rng = np.random.default_rng()
    >>> X = rng.normal(size=(n, k))
    >>> X = X.astype(np.float32)

    >>> eigenvalues, _ = lobpcg(A, X, maxiter=60)
    >>> eigenvalues
    array([100.], dtype=float32)

    `lobpcg` needs only access the matrix product with `A` rather
    then the matrix itself. Since the matrix `A` is diagonal in
    this example, one can write a function of the matrix product
    ``A @ X`` using the diagonal values ``vals`` only, e.g., by
    element-wise multiplication with broadcasting in the lambda-function

    >>> A_lambda = lambda X: vals[:, np.newaxis] * X

    or the regular function

    >>> def A_matmat(X):
    ...     return vals[:, np.newaxis] * X

    and use the handle to one of these callables as an input

    >>> eigenvalues, _ = lobpcg(A_lambda, X, maxiter=60)
    >>> eigenvalues
    array([100.], dtype=float32)
    >>> eigenvalues, _ = lobpcg(A_matmat, X, maxiter=60)
    >>> eigenvalues
    array([100.], dtype=float32)

    The traditional callable `LinearOperator` is no longer
    necessary but still supported as the input to `lobpcg`.
    Specifying ``matmat=A_matmat`` explicitly improves performance.

    >>> A_lo = LinearOperator((n, n), matvec=A_matmat, matmat=A_matmat, dtype=np.int16)
    >>> eigenvalues, _ = lobpcg(A_lo, X, maxiter=80)
    >>> eigenvalues
    array([100.], dtype=float32)

    The least efficient callable option is `aslinearoperator`:

    >>> eigenvalues, _ = lobpcg(aslinearoperator(A), X, maxiter=80)
    >>> eigenvalues
    array([100.], dtype=float32)

    We now switch to computing the three smallest eigenvalues specifying

    >>> k = 3
    >>> X = np.random.default_rng().normal(size=(n, k))

    and ``largest=False`` parameter

    >>> eigenvalues, _ = lobpcg(A, X, largest=False, maxiter=90)
    >>> print(eigenvalues)
    [1. 2. 3.]

    The next example illustrates computing 3 smallest eigenvalues of
    the same matrix `A` given by the function handle ``A_matmat`` but
    with constraints and preconditioning.

    Constraints - an optional input parameter is a 2D array comprising
    of column vectors that the eigenvectors must be orthogonal to

    >>> Y = np.eye(n, 3)

    The preconditioner acts as the inverse of `A` in this example, but
    in the reduced precision ``np.float32`` even though the initial `X`
    and thus all iterates and the output are in full ``np.float64``.

    >>> inv_vals = 1./vals
    >>> inv_vals = inv_vals.astype(np.float32)
    >>> M = lambda X: inv_vals[:, np.newaxis] * X

    Let us now solve the eigenvalue problem for the matrix `A` first
    without preconditioning requesting 80 iterations

    >>> eigenvalues, _ = lobpcg(A_matmat, X, Y=Y, largest=False, maxiter=80)
    >>> eigenvalues
    array([4., 5., 6.])
    >>> eigenvalues.dtype
    dtype('float64')

    With preconditioning we need only 20 iterations from the same `X`

    >>> eigenvalues, _ = lobpcg(A_matmat, X, Y=Y, M=M, largest=False, maxiter=20)
    >>> eigenvalues
    array([4., 5., 6.])

    Note that the vectors passed in `Y` are the eigenvectors of the 3
    smallest eigenvalues. The results returned above are orthogonal to those.

    The primary matrix `A` may be indefinite, e.g., after shifting
    ``vals`` by 50 from 1, ..., 100 to -49, ..., 50, we still can compute
    the 3 smallest or largest eigenvalues.

    >>> vals = vals - 50
    >>> X = rng.normal(size=(n, k))
    >>> eigenvalues, _ = lobpcg(A_matmat, X, largest=False, maxiter=99)
    >>> eigenvalues
    array([-49., -48., -47.])
    >>> eigenvalues, _ = lobpcg(A_matmat, X, largest=True, maxiter=99)
    >>> eigenvalues
    array([50., 49., 48.])

    """
    blockVectorX = X
    bestblockVectorX = blockVectorX
    blockVectorY = Y
    residualTolerance = tol
    if maxiter is None:
        maxiter = 20

    bestIterationNumber = maxiter

    blockVectorX, blockVectorY = _check_X_Y(blockVectorX, blockVectorY)
    n, sizeX = blockVectorX.shape
    sizeY = 0 if blockVectorY is None else blockVectorY.shape[1]

    # The histories are always recorded, but only returned if requested.
    lambdaHistory = np.zeros((max(maxiter, 0) + 3, sizeX),
                             dtype=blockVectorX.dtype)
    residualNormsHistory = np.zeros_like(lambdaHistory)

    if verbosityLevel:
        _print_problem_info(B, M, blockVectorY, n, sizeX, sizeY)

    if (n - sizeY) < (5 * sizeX):
        warnings.warn(
            f"The problem size {n} minus the constraints size {sizeY} "
            f"is too small relative to the block size {sizeX}. "
            f"Using a dense eigensolver instead of LOBPCG iterations."
            f"No output of the history of the iterations.",
            UserWarning, stacklevel=2
        )

        if blockVectorY is not None:
            raise NotImplementedError(
                "The dense eigensolver does not support constraints."
            )

        return _dense_eigh(A, B, n, min(sizeX, n), largest)

    if (residualTolerance is None) or (residualTolerance <= 0.0):
        residualTolerance = np.sqrt(np.finfo(blockVectorX.dtype).eps) * n

    A = _makeMatMat(A)
    B = _makeMatMat(B)
    M = _makeMatMat(M)

    # Apply constraints to X.
    if blockVectorY is not None:

        if B is not None:
            blockVectorBY = _apply_checked(B, blockVectorY,
                                           "constraint", "secondary")
        else:
            blockVectorBY = blockVectorY

        # gramYBY is a dense array.
        gramYBY = blockVectorY.T.conj() @ blockVectorBY
        try:
            # gramYBY is a Cholesky factor from now on...
            gramYBY = cho_factor(gramYBY, overwrite_a=True)
        except LinAlgError as e:
            raise ValueError("Linearly dependent constraints") from e

        _applyConstraints(blockVectorX, gramYBY, blockVectorBY, blockVectorY)

    ##
    # B-orthonormalize X.
    blockVectorX, blockVectorBX, _ = _b_orthonormalize(
        B, blockVectorX, verbosityLevel=verbosityLevel)
    if blockVectorX is None:
        raise ValueError("Linearly dependent initial approximations")

    ##
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = _apply_checked(A, blockVectorX,
                                   "initial approximations", "primary")

    gramXAX = blockVectorX.T.conj() @ blockVectorAX

    _lambda, eigBlockVector = eigh(gramXAX, check_finite=False)
    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    lambdaHistory[0, :] = _lambda

    eigBlockVector = np.asarray(eigBlockVector[:, ii])
    blockVectorX = _matmul_inplace(
        blockVectorX, eigBlockVector,
        verbosityLevel=verbosityLevel
    )
    blockVectorAX = _matmul_inplace(
        blockVectorAX, eigBlockVector,
        verbosityLevel=verbosityLevel
    )
    if B is not None:
        blockVectorBX = _matmul_inplace(
            blockVectorBX, eigBlockVector,
            verbosityLevel=verbosityLevel
        )

    ##
    # Active index set.
    activeMask = np.ones((sizeX,), dtype=bool)

    ##
    # Main iteration loop.

    blockVectorP = None  # set during iteration
    blockVectorAP = None
    blockVectorBP = None
    activeBlockVectorP = None
    activeBlockVectorAP = None
    activeBlockVectorBP = None

    smallestResidualNorm = np.abs(np.finfo(blockVectorX.dtype).max)

    iterationNumber = -1
    restart = True
    forcedRestart = False
    explicitGramFlag = False
    while iterationNumber < maxiter:
        iterationNumber += 1

        blockVectorR, residualNorms = _residuals(
            B, blockVectorX, blockVectorAX, blockVectorBX, _lambda)
        residualNormsHistory[iterationNumber, :] = residualNorms
        residualNorm = np.sum(np.abs(residualNorms)) / sizeX

        if residualNorm < smallestResidualNorm:
            smallestResidualNorm = residualNorm
            bestIterationNumber = iterationNumber
            bestblockVectorX = blockVectorX
        elif residualNorm > 2**restartControl * smallestResidualNorm:
            forcedRestart = True
            blockVectorAX = _apply_checked(A, blockVectorX,
                                           "restarted iterate", "primary")
            if B is not None:
                blockVectorBX = _apply_checked(B, blockVectorX,
                                               "restarted iterate",
                                               "secondary")

        ii = np.where(residualNorms > residualTolerance, True, False)
        activeMask = activeMask & ii
        currentBlockSize = activeMask.sum()

        if verbosityLevel:
            print(f"iteration {iterationNumber}")
            print(f"current block size: {currentBlockSize}")
            print(f"eigenvalue(s):\n{_lambda}")
            print(f"residual norm(s):\n{residualNorms}")

        if currentBlockSize == 0:
            break

        activeBlockVectorR = _as2d(blockVectorR[:, activeMask])

        if iterationNumber > 0:
            activeBlockVectorP = _as2d(blockVectorP[:, activeMask])
            activeBlockVectorAP = _as2d(blockVectorAP[:, activeMask])
            if B is not None:
                activeBlockVectorBP = _as2d(blockVectorBP[:, activeMask])

        if M is not None:
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = M(activeBlockVectorR)

        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            _applyConstraints(activeBlockVectorR,
                              gramYBY,
                              blockVectorBY,
                              blockVectorY)

        ##
        # B-orthogonalize the preconditioned residuals to X.
        activeBlockVectorR = _b_orthogonalize_to_X(
            B, blockVectorX, blockVectorBX, activeBlockVectorR)

        ##
        # B-orthonormalize the preconditioned residuals.
        aux = _b_orthonormalize(
            B, activeBlockVectorR, verbosityLevel=verbosityLevel)
        activeBlockVectorR, activeBlockVectorBR, _ = aux

        if activeBlockVectorR is None:
            warnings.warn(
                f"Failed at iteration {iterationNumber} with accuracies "
                f"{residualNorms}\n not reaching the requested "
                f"tolerance {residualTolerance}.",
                UserWarning, stacklevel=2
            )
            break
        activeBlockVectorAR = A(activeBlockVectorR)

        if iterationNumber > 0:
            if B is not None:
                aux = _b_orthonormalize(
                    B, activeBlockVectorP, activeBlockVectorBP,
                    verbosityLevel=verbosityLevel
                )
                activeBlockVectorP, activeBlockVectorBP, invR = aux
            else:
                aux = _b_orthonormalize(B, activeBlockVectorP,
                                        verbosityLevel=verbosityLevel)
                activeBlockVectorP, _, invR = aux
            # Function _b_orthonormalize returns None if Cholesky fails
            if activeBlockVectorP is not None:
                activeBlockVectorAP = _matmul_inplace(
                    activeBlockVectorAP, invR,
                    verbosityLevel=verbosityLevel
                )
                restart = forcedRestart
            else:
                restart = True

        ##
        # Perform the Rayleigh Ritz Procedure.
        explicitGramFlag = _explicit_gram_flag(
            explicitGramFlag, residualNorms,
            activeBlockVectorR, activeBlockVectorAR)

        # Shared memory assignments to simplify the code
        if B is None:
            blockVectorBX = blockVectorX
            activeBlockVectorBR = activeBlockVectorR
            if not restart:
                activeBlockVectorBP = activeBlockVectorP

        aux = _rayleigh_ritz(
            blockVectorX, blockVectorAX, blockVectorBX,
            activeBlockVectorR, activeBlockVectorAR, activeBlockVectorBR,
            activeBlockVectorP, activeBlockVectorAP, activeBlockVectorBP,
            _lambda, sizeX, currentBlockSize, explicitGramFlag, restart,
            iterationNumber, verbosityLevel)
        if aux is None:
            # Keep the previous _lambda for the postprocessing below.
            break
        _lambda, eigBlockVector, restart = aux

        ii = _get_indx(_lambda, sizeX, largest)
        _lambda = _lambda[ii]
        eigBlockVector = eigBlockVector[:, ii]
        lambdaHistory[iterationNumber + 1, :] = _lambda

        # Compute Ritz vectors.
        eigBlockVectorX, eigBlockVectorR, eigBlockVectorP = _split_eig_block(
            eigBlockVector, sizeX, currentBlockSize, restart)

        pp = _new_direction(activeBlockVectorR, activeBlockVectorP,
                            eigBlockVectorR, eigBlockVectorP, restart)
        app = _new_direction(activeBlockVectorAR, activeBlockVectorAP,
                             eigBlockVectorR, eigBlockVectorP, restart)

        blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
        blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app
        blockVectorP, blockVectorAP = pp, app

        if B is not None:
            bpp = _new_direction(activeBlockVectorBR, activeBlockVectorBP,
                                 eigBlockVectorR, eigBlockVectorP, restart)
            blockVectorBX = np.dot(blockVectorBX, eigBlockVectorX) + bpp
            blockVectorBP = bpp

    _, residualNorms = _residuals(
        B, blockVectorX, blockVectorAX, blockVectorBX, _lambda)
    # Use old lambda in case of early loop exit.
    lambdaHistory[iterationNumber + 1, :] = _lambda
    residualNormsHistory[iterationNumber + 1, :] = residualNorms
    residualNorm = np.sum(np.abs(residualNorms)) / sizeX
    if residualNorm < smallestResidualNorm:
        smallestResidualNorm = residualNorm
        bestIterationNumber = iterationNumber + 1
        bestblockVectorX = blockVectorX

    if np.max(np.abs(residualNorms)) > residualTolerance:
        warnings.warn(
            f"Exited at iteration {iterationNumber} with accuracies \n"
            f"{residualNorms}\n"
            f"not reaching the requested tolerance {residualTolerance}.\n"
            f"Use iteration {bestIterationNumber} instead with accuracy \n"
            f"{smallestResidualNorm}.\n",
            UserWarning, stacklevel=2
        )

    if verbosityLevel:
        print(f"Final iterative eigenvalue(s):\n{_lambda}")
        print(f"Final iterative residual norm(s):\n{residualNorms}")

    blockVectorX = bestblockVectorX
    # Making eigenvectors "exactly" satisfy the blockVectorY constrains
    if blockVectorY is not None:
        _applyConstraints(blockVectorX,
                          gramYBY,
                          blockVectorBY,
                          blockVectorY)

    _lambda, blockVectorX, residualNorms = _postprocess(
        A, B, blockVectorX, sizeX, largest, verbosityLevel)

    result = (_lambda, blockVectorX)
    if retLambdaHistory:
        result += (_finalize_history(lambdaHistory, bestIterationNumber,
                                     _lambda),)
    if retResidualNormsHistory:
        result += (_finalize_history(residualNormsHistory,
                                     bestIterationNumber, residualNorms),)

    if np.max(np.abs(residualNorms)) > residualTolerance:
        warnings.warn(
            f"Exited postprocessing with accuracies \n"
            f"{residualNorms}\n"
            f"not reaching the requested tolerance {residualTolerance}.",
            UserWarning, stacklevel=2
        )

    if verbosityLevel:
        print(f"Final postprocessing eigenvalue(s):\n{_lambda}")
        print(f"Final residual norm(s):\n{residualNorms}")

    return result
