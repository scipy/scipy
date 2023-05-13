import numpy as np
from .utils import make_system


__all__ = ['symmlq']


def symmlq(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None,
           atol=0., verbose=False):
    """
    Use Symmetric LQ iteration to solve ``Ax = b``.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        `scipy.sparse.linalg.LinearOperator`.
    b : {ndarray}
        Right hand side of the linear system. Has shape (N,) or (N, 1).
    x0 : {ndarray}
        Starting guess for the solution.
    tol, atol : float, optional
        Tolerances for convergence, `tol` represents relative tolerance and
        `atol` represents absolute tolerance,
        ``norm(b - Ax) < max(tol*norm(b), atol)``.
        The default for `tol` is 1.0e-5.
        The default for `atol` is 0.
    maxiter : int, optional
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
        Default is ``min(10000, n * 10)``, where ``n = A.shape[0]``.
    M : {sparse matrix, ndarray, LinearOperator}
        Inverse of the preconditioner of A. M should approximate the
        inverse of A and be easy to solve for (see Notes). Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance. By default, no preconditioner is used.
    callback : function, optional
        User-supplied function to call after each iteration. It is called
        as `callback(xk)`, where `xk` is the current solution vector.
    verbose : bool, optional
        Specify ``verbose = True`` to show the convergence, ``verbose = False`` is
        to close the output of the convergence.
        Default is `False`.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:

            - 0  : successful exit
            - >0 : convergence to tolerance not achieved, number of iterations
            - <0 : illegal input or breakdown

    Notes
    -----
    The Symmetric LQ algorithm (SYMMLQ) is a variant of CG method, which is
    designed to avoid LU decomposition and breakdown, and can be applied to
    symmetric indefinite systems. Unlike MINRES, SYMMLQ keeps the residual
    orthogonal to all previous ones (Don't minimize anything).

    References
    ----------
    .. [1] C. C. Paige, M. A. Saunders, Solution of Sparse Indefinite Systems
           of Linear Equations, SIAM J. Numer. Anal., Vol. 12, No. 4, 1975.

    This file is adapted from PETSc implementation by
        https://gitlab.com/petsc/petsc/-/blob/main/src/ksp/ksp/impls/symmlq/symmlq.c

    Examples
    --------
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import symmlq
    >>> A = csc_matrix([[2, -1, 0], [-1, 0.25, -1], [0, -1, 0]], dtype=float)
    >>> b = np.array([1, 2, 3], dtype=float)
    >>> x, exitCode = symmlq(A, b)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """

    A, M, x, b, postprocess = make_system(A, M, x0, b)

    n = A.shape[0]
    if maxiter is None:
        maxiter = min(10000, n * 10)

    # Check if the R.H.S is a zero vector
    bnorm = np.linalg.norm(b)
    if bnorm == 0.:
        x = np.zeros(n)
        if verbose:
            print("SYMMLQ: Linear solve converged due to zero right-hand side "
                  "iterations 0")
        return (postprocess(x), 0)

    # Check initial guess
    if x0 is None:
        r = b
    else:
        r = b - A.matvec(x)

    # Set absolute tolerance
    atol = max(atol, tol * bnorm)

    c = cold = 1.
    s = sold = 0.
    ceta_oold = ceta_old = ceta = 0.
    w = wbar = vold = uold = np.zeros(n)
    z = M.matvec(r)

    # Check initial residual norm
    rnorm = np.linalg.norm(z)  # preconditioned residual norm
    if rnorm <= atol:
        if callback is not None:
            callback(x)
        if verbose:
            print("SYMMLQ: Linear solve converged due to reach TOL "
                  "iterations 0")
        return (postprocess(x), 0)

    # Check breakdown
    dp = np.inner(r.conjugate(), z)
    haptol = 1e-18
    if abs(dp) < haptol:
        if callback is not None:
            callback(x)
        if verbose:
            print("SYMMLQ: Linear solve converged due to HAPPY BREAKDOWN "
                  "iterations 0")
        return (postprocess(x), 0)

    if not np.iscomplex(dp) and dp < 0.:
        if callback is not None:
            callback(x)
        if verbose:
            print("SYMMLQ: Linear solve not converged due to DIVERGED "
                  "INDEFINITE PC iterations 0")
        return (postprocess(x), -1)

    beta = beta1 = np.sqrt(dp)
    s_prod = abs(beta1)
    [v, u] = [r, z] / beta
    wbar = u

    for iter in range(maxiter):
        # Update
        if iter > 0:
            vold, uold = v, u
            [v, u] = [r, z] / beta
            w = c * wbar + s * u
            wbar = - s * wbar + c * u
            x += ceta * w
            ceta_oold, ceta_old = ceta_old, ceta

        # Lanczos method
        r = A.matvec(u)
        alpha = np.inner(u.conjugate(), r)
        z = M.matvec(r)
        r = r - alpha * v - beta * vold
        z = z - alpha * u - beta * uold
        betaold = beta
        dp = np.inner(r.conjugate(), z)
        # Check breakdown
        if abs(dp) < haptol:
            dp = 0.
        if not np.iscomplex(dp) and dp < 0.:
            if c == 0.:
                ceta_bar = ceta * 1e15
            else:
                ceta_bar = ceta / c
            x += ceta_bar * wbar
            if callback is not None:
                callback(x)
            if verbose:
                print("SYMMLQ: Linear solve not converged due to DIVERGED "
                      f"INDEFINITE PC iterations {iter+1}")
            return (postprocess(x), -1)
        beta = np.sqrt(dp)

        # QR factorization
        coold, cold, soold, sold = cold, c, sold, s
        rho0 = cold * alpha - coold * sold * betaold
        rho1 = np.sqrt(rho0**2 + beta**2)
        rho2 = sold * alpha + coold * cold * betaold
        rho3 = soold * betaold

        # Givens rotation
        [c, s] = [rho0, beta] / rho1
        if iter == 0:
            ceta = beta1 / rho1
        else:
            ceta = - (rho2 * ceta_old + rho3 * ceta_oold) / rho1

        s_prod *= abs(s)
        if c == 0.:
            rnorm = s_prod * 1e15
        else:
            rnorm = s_prod / abs(c)

        # Convergence criterion
        if rnorm <= atol:
            if verbose:
                print("SYMMLQ: Linear solve converged due to reach TOL "
                      f"iterations {iter+1}")
            break

        if callback is not None and iter != maxiter-1:
            callback(x)

    if c == 0.:
        ceta_bar = ceta * 1e15
    else:
        ceta_bar = ceta / c

    x += ceta_bar * wbar

    if callback is not None:
        callback(x)

    if rnorm <= atol:
        return (postprocess(x), 0)
    else:
        if verbose:
            print("SYMMLQ: Linear solve not converged due to reach MAXIT "
                  f"iterations {maxiter}")
        return (postprocess(x), maxiter)
