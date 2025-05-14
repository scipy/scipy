"""Restared Krylov method for evaluating f(A)b"""

import numpy as np
from ._isolve.iterative import _get_atol_rtol

__all__ = ['krylov_funmv']

def _krylov_funmv_arnoldi(A, b, bnorm, V, H, m):
    """
    The Arnoldi iteration for constructing the basis V and the projection H = V * A V
    for the Krylov subspace Km(A, b) of order m.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose matrix function is of interest.
    b : ndarray
        The vector b to multiply the f(A) with.
    V : ndarray
        The n x (m + 1) matrix whose columns determines the basis for
        Krylov subspace Km(A, b).
    H : ndarray
        A (m + 1) x m upper Hessenberg matrix representing the projection of A
        onto Km(A, b).
    m : int
        The order of the Krylov subspace.

    """

    dotprod = np.vdot if np.iscomplexobj(b) else np.dot
    V[:, 0] = b / bnorm

    for k in range(0, m):
        V[:, k + 1] = A.dot(V[:, k])

        # Uses the modified Gram-Schmift process to orthogonalize V[:, k + 1]
        # against the previous basis vectors
        for i in range(0, k + 1):
            H[i, k] = dotprod(V[:, i], V[:, k + 1])
            V[:, k + 1] = V[:, k + 1] - H[i, k] * V[:, i]

        H[k + 1, k] = np.linalg.norm(V[:, k + 1])
        if H[k + 1, k] != 0:
            V[:, k + 1] = V[:, k + 1] / H[k + 1, k]
        else:
            raise RuntimeError("krylov_funmv: The Arnoldi iteration broke down "
                               f"at k = {k}, i.e., V[:, k + 1] = 0")

def _krylov_funmv_lanczos(A, b, bnorm, V, H, m):
    """
    The Lanczos iteration for constructing the basis V and the projection H = V * A V
    for the Krylov subspace Km(A, b) of order m. A must be Hermitian.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose matrix function is of interest.
    b : ndarray
        The vector b to multiply the f(A) with.
    V : ndarray
        The n x (m + 1) matrix whose columns determines the basis for
        Krylov subspace Km(A, b).
    H : ndarray
        A (m + 1) x m upper Hessenberg matrix representing the projection of A
        onto Km(A, b).
    m : int
        The order of the Krylov subspace.

    """
    dotprod = np.vdot if np.iscomplexobj(b) else np.dot
    V[:, 0] = b / bnorm

    for k in range(0, m):
        if k > 0:
            V[:, k + 1] = A.dot(V[:, k]) - H[k, k - 1] * V[:, k - 1]
        else:
            V[:, k + 1] = A.dot(V[:, k])

        H[k, k] = dotprod(V[:, k + 1], V[:, k])
        V[:, k + 1] = V[:, k + 1] - H[k, k] * V[:, k]

        H[k + 1, k] = np.linalg.norm(V[:, k + 1])
        if H[k + 1, k] != 0:
            V[:, k + 1] = V[:, k + 1] / H[k + 1, k]
            if k < m - 1:
                H[k, k + 1] = H[k + 1, k]
        else:
            raise RuntimeError("krylov_funmv: The Lanczos iteration broke down "
                               f"at k = {k}, i.e., V[:, k + 1] = 0")


def krylov_funmv(f, t, A, b, atol = 0.0, btol = 1e-6, restart_length = None,
                 max_restarts = 20, ortho_method = "arnoldi", verbose = False):
    """
    A restarted Krylov method for evaluating ``y = f(tA) b``.

    Parameters
    ----------
    f : function
        Callable object that computes the matrix function ``F = f(X)``.

    t : float
        The value to scale the matrix ``A`` with.

    A : {sparse array, ndarray, LinearOperator}
        A real or complex N-by-N matrix.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g., ``scipy.sparse.linalg.LinearOperator``.

    b : ndarray, sparse array
        A vector to multiply the ``f(tA)`` with.

    atol, btol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(||y_k - y_k-1||) <= max(btol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``btol=1e-6``.

    restart_length : integer
        Number of iterations between restarts. Larger values increase
        iteration cost, but may be necessary for convergence.
        If omitted, ``min(20, n)`` is used.

    max_restarts : int, optional
        Maximum number of restart cycles.  The algorithm will stop
        after max_restarts cycles even if the specified tolerance has not been
        achieved. The default is ``max_restarts=20``

    ortho_method : string, optional
        Orthogonalization method to use. ``lanczos`` is faster, but it requires
        that ``A`` is Hermitian. ``arnoldi`` is the default option.

    verbose : bool
        Prints iteration logs if ``verbose=True``. Default is ``False``.

    Returns
    -------
    y : ndarray
        The resulting vector from the operation

    Notes
    ______

    The convergence of the Krylov method heavily depends on the spectrum
    of ``A`` and the function ``f``. With restarting, there are only formal
    proofs for functions of order 1 (e.g., ``exp``, ``sin``, ``cos``) and
    Stieltjes functions [2, 3], while the general case remains an open problem.

    Examples
    ________

    >>> import numpy as np
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.linalg import krylov_funmv
    >>> from scipy.linalg import expm, solve
    >>> A = csr_array([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> t = 0.1

    Compute ``y = exp(tA) b``.

    >>> y = krylov_funmv(expm, t, A, b)
    [3.6164913 3.88421511 0.96073457]

    >>> ref = expm(t * A.todense()) @ b
    >>> y - ref
    [4.44089210e-16 0.00000000e+00 2.22044605e-16]

    Compute ``y = phi_1(tA) b``, where

    .. math::

        \phi_1(A) = A^{-1}(e^{A} - I).


    >>> def phim_1(X):
    >>>     return solve(X, expm(X) - np.eye(X.shape[0]))
    >>> y = krylov_funmv(phim_1, t, A, b)
    [ 2.76984306  3.92769192 -0.03111392]

    >>> ref = phim_1(t * A.todense()) @ b
    >>> y - ref
    [ 0.00000000e+00  8.88178420e-16 -4.60742555e-15]

    References
    ----------
    .. [1] M. Afanasjew, M. Eiermann, O. G. Ernst, and S. Güttel,
          “Implementation of a restarted Krylov subspace method for the
          evaluation of matrix functions,” Linear Algebra and its Applications,
           vol. 429, no. 10, pp. 2293–2314, Nov. 2008, doi: 10.1016/j.laa.2008.06.029.

    .. [2] M. Eiermann and O. G. Ernst, “A Restarted Krylov Subspace Method
           for the Evaluation of Matrix Functions,” SIAM J. Numer. Anal., vol. 44,
           no. 6, pp. 2481–2504, Jan. 2006, doi: 10.1137/050633846.

    .. [3] A. Frommer, S. Güttel, and M. Schweitzer, “Convergence of Restarted
           Krylov Subspace Methods for Stieltjes Functions of Matrices,” SIAM J.
           Matrix Anal. Appl., vol. 35, no. 4, pp. 1602–1624,
           Jan. 2014, doi: 10.1137/140973463.

    """

    if len(b.shape) != 1:
        raise RuntimeError("krylov_funmv: b must be a vector.")

    n = b.shape[0]

    if restart_length is None:
        restart_length = min(20, n)
    m = restart_length
    max_restarts = min(max_restarts, int(n / restart_length) + 1)
    mmax = restart_length * max_restarts

    bnorm = np.linalg.norm(b)
    atol, _ = _get_atol_rtol("krylov_funmv", bnorm, atol, btol)

    if bnorm == 0:
        y = b
        return y

    # Pre-allocate the maximum memory space.
    # Using the column major order here since we work with
    # each individual column separately.
    V = np.zeros((n, m + 1), dtype = b.dtype, order = 'F')
    H = np.zeros((mmax + 1, mmax), dtype = b.dtype, order = 'F')
    y = np.zeros_like(b)

    if verbose:
        print(f"N = {n}")
        print(f"f = {f.__name__}")
        print(f"atol = {atol:8.2e}, btol = {btol:8.2e}")
        print(f"restart_length = {restart_length}, max_restarts = {max_restarts}")
        print(f"ortho_method = {ortho_method}")
        print(f"||b|| = {bnorm}")

    restart = 1
    if ortho_method == "lanczos":
        _krylov_funmv_lanczos(A, b, bnorm, V, H[:m + 1, :m], m)
    elif ortho_method == "arnoldi":
        _krylov_funmv_arnoldi(A, b, bnorm, V, H[:m + 1, :m], m)
    else:
        raise RuntimeError("krylov_funmv: Invalid orthogonalization method. "
                           "Available options: 'arnoldi' or 'lanczos'")

    fH = f(t * H[:m, :m])
    y = bnorm * V[:, :m].dot(fH[:, 0])
    update_norm = np.linalg.norm(bnorm * fH[:, 0])

    if verbose:
        print("{:^10}{:^10}{:^10}".format('restart', '||y_k||', '||y_k - y_k-1||'))
        print(f"{restart:^10}{np.linalg.norm(y):^10.3g}{update_norm:^10.3g}")


    while restart < max_restarts and update_norm > atol:
        begin = restart * m
        end = (restart + 1) * m

        if ortho_method == "lanczos":
            _krylov_funmv_lanczos(A, V[:, m], 1, V, H[begin:end + 1, begin:end], m)
        elif ortho_method == "arnoldi":
            _krylov_funmv_arnoldi(A, V[:, m], 1, V, H[begin:end + 1, begin:end], m)

        fH = f(t * H[:end, :end])
        y = y + bnorm * V[:, :m].dot(fH[begin:end, 0])
        update_norm = np.linalg.norm(bnorm * fH[begin:end, 0])
        restart += 1

        if verbose:
            print(f"{restart:^10}{np.linalg.norm(y):^10.3g}{update_norm:^10.3g}")

    if verbose:
        print("\n")

    return y
