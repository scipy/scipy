"""Restart GMRES method for solving asymmetric linear systems"""

import numpy as np
from .utils import make_system

__all__ = ['gmresk']


def gmresk(A, b, x0=None, restart=None, tol=1e-5, maxiter=None, M=None,
           callback=None, atol=None, callback_type=None, show=False):
    """
    Use Restarted Generalized Minimal RESidual based on modified Gram-Schmidt
    iteration to solve ``Ax = b``. This version is re-entrant and faster than
    the old version.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        `scipy.sparse.linalg.LinearOperator`.
    b : {ndarray}
        Right hand side of the linear system. Has shape (N,) or (N,1).

    Returns
    -------
    x : {ndarray}
        The converged solution.
    info : int
        Provides convergence information:
          * 0  : successful exit
          * >0 : convergence to tolerance not achieved, number of iterations
          * <0 : illegal input or breakdown

    Other parameters
    ----------------
    x0 : {ndarray}
        Starting guess for the solution (a vector of zeros by default).
    restart : int, optional
        Number of iterations between restarts. Larger values increase
        iteration cost, but may be necessary for convergence.
        Default is 20.
    tol : float, optional
        Tolerances for convergence, ``norm(residual) <= tol*norm(b)``.
        Default is 1.0e-5.
    maxiter : int, optional
        Maximum number of iterations (restart cycles).  Iteration will stop
        after maxiter steps even if the specified tolerance has not been
        achieved.
    M : {sparse matrix, ndarray, LinearOperator}
        Inverse of the preconditioner of A.  M should approximate the
        inverse of A and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.
    callback : function
        User-supplied function to call after each iteration.
    callback_type : {None, 'x', 'prnorm'}, optional
        Callback function argument requested:
          - ``None``: preconditioned residual norm (float), called on every
            inner iteration.
          - ``x``: current iterate (ndarray), called on every restart.
          - ``prnorm``: (absolute) preconditioned residual norm (float),
            called on every inner iteration.
    show : bool, optional
        Specify ``show = True`` to show the convergence, ``show = False`` is
        to close the output of the convergence.
        Default is `False`.

    See Also
    --------
    LinearOperator

    Notes
    -----
    A preconditioner ``P`` is chosen such that ``P`` is close to `A` but easy
    to solve for. The preconditioner parameter required by this routine is
    ``M = P^-1``. The inverse should preferably not be calculated
    explicitly.  Rather, use the following template to produce ``M``:

    # Construct a linear operator that computes P^-1 @ x.
    import scipy.sparse.linalg as spla
    M_x = lambda x: spla.spsolve(P, x)
    M = spla.LinearOperator((n, n), M_x)

    Examples
    --------
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import gmresk
    >>> A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> x, exitCode = gmresk(A, b)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """

    # Check data type
    dtype = A.dtype
    if np.issubdtype(dtype, np.int64):
        dtype = float
        A = A.astype(dtype)
    if np.issubdtype(b.dtype, np.int64):
        b = b.astype(dtype)

    A, M, x, b, postprocess = make_system(A, M, x0, b)

    # Check if RHS is a zero vector
    bnorm = np.linalg.norm(b)
    if bnorm == 0.:
        x = b.copy()
        if (show):
            print("GMRES: Linear solve converged due to zero RHS iterations 0")
        return (postprocess(x), 0)

    ndofs = A.shape[0]
    if maxiter is None:
        maxiter = min(10000, ndofs * 10)
    if x0 is None:
        x0 = x.copy()
    if restart is None:
        restart = 20
    m = min(restart, ndofs)

    # Define orthonormal basis
    V = np.zeros((m+1, ndofs), dtype=dtype)
    # Define Hessenberg matrix
    Hn = np.zeros((m+1, m), dtype=dtype)
    # Define rotation matrix
    rot = np.zeros((2, m), dtype=dtype)
    g = np.zeros((m+1), dtype=dtype)
    gtmp = np.zeros((m+1), dtype=dtype)

    r = b - A.matvec(x)
    z = M.matvec(r)
    g[0] = np.linalg.norm(z)  # initial residual

    if callback_type is None:
        callback_type = 'prnorm'

    # Initial callback
    if callback is not None:
        if callback_type == 'x':
            callback(x)
        elif callback_type == 'prnorm':
            callback(g[0])

    if g[0] == 0.:
        if (show):
            print("GMRES: Linear solve converged due to zero preconditioned "
                  "residual norm iterations 0")
        return (postprocess(x), 0)

    if atol is None:
        atol = tol * bnorm
    else:
        atol = max(atol, tol * bnorm)

    # Outer iterations
    niter = 0
    for iter in range(maxiter):
        V[0] = (1./g[0]) * z  # first basis vector

        # Inner iterations
        for k in range(m):
            niter += 1

            # Modified Gram-Schmidt orthogonalization
            v = A.matvec(V[k])
            w = M.matvec(v)
            Hn[0:k+1, k] = V[0:k+1].conjugate().dot(w)
            w -= V[0:k+1].T.dot(Hn[0:k+1, k])
            Hn[k+1, k] = np.linalg.norm(w)

            # Check if Hn[k+1, k] is zero
            if Hn[k+1, k] == 0.:
                # Reconstruct the solution
                y = np.zeros((k+1), dtype=dtype)
                gtmp[0:k+1] = g[0:k+1]
                for i in range(k, -1, -1):
                    gtmp[i] -= Hn[i, i+1:k+1].dot(y[i+1:k+1])
                    if Hn[i, i] != 0:
                        y[i] = gtmp[i] / Hn[i, i]
                x += V[0:k+1].T.dot(y[0:k+1])
                if callback is not None:
                    if callback_type == 'x':
                        callback(x)
                    elif callback_type == 'prnorm':
                        r = b - A.matvec(x)
                        z = M.matvec(r)
                        g[0] = np.linalg.norm(z)
                        callback(g[0])
                if (show):
                    print("GMRES: Linear solve not converged due to BREAKDOWN "
                          f"iterations {niter}")
                return (postprocess(x), -1)
            V[k+1] = (1./Hn[k+1, k]) * w

            # QR decomposition of Hn
            for i in range(k):
                tmp0 = rot[0, i] * Hn[i, k] + rot[1, i] * Hn[i+1, k]
                tmp1 = -rot[1, i] * Hn[i, k] + rot[0, i] * Hn[i+1, k]
                Hn[i, k], Hn[i+1, k] = tmp0, tmp1
            sq = np.sqrt(Hn[k+1, k]**2 + Hn[k, k]**2)
            rot[0:2, k] = Hn[k:k+2, k] / sq
            Hn[k, k], Hn[k+1, k] = sq, 0.
            g[k+1] = -rot[1, k] * g[k]
            g[k] *= rot[0, k]

            # Convergence criterion in inner iterations
            if k < m - 1:
                if np.abs(g[k+1]) < atol:
                    break
            else:
                break

            # Inner callback
            if (callback is not None) and (callback_type == 'prnorm'):
                callback(np.abs(g[k+1]))

        # Reconstruct the solution
        y = np.zeros((k+1), dtype=dtype)
        gtmp[0:k+1] = g[0:k+1]
        for i in range(k, -1, -1):
            gtmp[i] -= Hn[i, i+1:k+1].dot(y[i+1:k+1])
            y[i] = gtmp[i] / Hn[i, i]
        x += V[0:k+1].T.dot(y[0:k+1])

        # Computing residual norm and preconditioned residual norm
        r = b - A.matvec(x)
        z = M.matvec(r)
        g[k+1] = np.linalg.norm(r)
        g[0] = np.linalg.norm(z)

        # Outer callback
        if callback is not None:
            if callback_type == 'x':
                callback(x)
            elif callback_type == 'prnorm':
                callback(g[0])

        # Convergence criterion in outer iterations
        if g[0] < atol:
            if (show):
                print("GMRES: Linear solve converged due to reach TOL "
                      f"iterations {niter}")
            return (postprocess(x), 0)

    if (show):
        print("GMRES: Linear solve not converged due to reach MAXIT "
              f"iterations {niter}")
    return (postprocess(x), maxiter)
