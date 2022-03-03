import numpy as np
from .utils import make_system


__all__ = ['cr']


def cr(A, b, x0=None, tol=1e-5, maxiter=None, M=None,
       callback=None, atol=None, dtol=1e+5, show=False):
    """
    Use Conjugate Residual (CR) iteration to solve ``Ax = b``.
    CR is structurally similar to Conjugate Gradient method,
    but mathematically equivalent to MINRES and the system can
    be symmetric indefinite or Hermitian.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator}
        The N-by-N Hermitian matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        `scipy.sparse.linalg.LinearOperator`.
    b : {ndarray}
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : {ndarray}
        Starting guess for the solution.
    tol, atol : float, optional
        Tolerances for convergence, ``norm(residual) <= max(tol*norm(b),
        atol)``.
        The default for `tol` is 1.0e-5.
        The default for `atol` is ``tol * norm(b)``.
    dtol : float, optional
        Tolerance for divergence, ``norm(residual) > dtol * norm(b)``.
        The default for `dtol` is 1.0e+5.
    maxiter : int, optional
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
        Default is ``min(10000, n * 10)``, where ``n = A.shape[0]``.
    M : {sparse matrix, ndarray, LinearOperator}
        Inverse of the preconditioner of A.  M should approximate the
        inverse of A and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.
    callback : function, optional
        User-supplied function to call after each iteration.  It is called
        as `callback(xk)`, where `xk` is the current solution vector.
    show : bool, optional
        Specify ``show = True`` to show the convergence, ``show = False`` is
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

    References
    ----------
    .. [1] Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd edition,
           SIAM, Philadelphia, 2003.
    .. [2] G. H. Golub, C. F. van Loan, Matrix Computations, 3rd ed., The Johns
           Hopkins University Press, Baltimore and London, 1996.

    Examples
    --------
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import cr
    >>> A = csc_matrix([[2, -1, 0], [-1, 2, -1], [0, -1, 2]], dtype=float)
    >>> b = np.ones(3, dtype=float)
    >>> x, exitCode = cr(A, b)
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
        if show:
            print("CR: Linear solve converged due to zero right-hand side "
                  "iterations 0")
        return (postprocess(x), 0)

    # Check if the initial guess is a zero vector
    if np.linalg.norm(x) == 0.:
        r = b.copy()
        rnorm = bnorm
    else:
        r = b - A.matvec(x)
        rnorm = np.linalg.norm(r)

    # Check if the residual norm is zero
    if rnorm == 0.:
        if show:
            print("CR: Linear solve converged due to zero residual norm "
                  "iterations 0")
        return (postprocess(x), 0)

    if atol is None:
        atol = tol * bnorm
    else:
        atol = max(atol, tol * bnorm)

    beta = 0.
    alpha = 0.

    Mr = M.matvec(r)
    AMr = A.matvec(Mr)
    # p is a search direction
    p = Mr.copy()
    Ap = A.matvec(p)
    for iter in range(maxiter):
        MAp = M.matvec(Ap)

        # Check breakdown
        MApTAp = np.inner(MAp.conjugate(), Ap)
        if MApTAp != 0.:
            # Compute step size on the search direction `p`
            alpha = np.inner(Mr.conjugate(), AMr) / MApTAp
        else:
            if show:
                print("CR: Linear solve not converged due to BREAKDOWN "
                      f"iterations {iter+1}")
            return (postprocess(x), -1)

        # Obtain the solution of each iteration
        x += alpha * p

        if callback is not None:
            callback(x)

        # Convergence criterion
        r = r - alpha * Ap
        rnorm = np.linalg.norm(r)
        if rnorm < atol:
            if show:
                print("CR: Linear solve converged due to reach TOL "
                      f"iterations {iter+1}")
            return (postprocess(x), 0)
        if rnorm > dtol * bnorm:
            if show:
                print("CR: Linear solve not converged due to reach "
                      f"DIVERGENCE TOL iterations {iter+1}")
            return (postprocess(x), iter+1)

        # Compute Mr and AMr
        Mr_old = Mr.copy()
        AMr_old = AMr.copy()
        Mr -= alpha * MAp
        AMr = A.matvec(Mr)

        # Check breakdown
        MrTAMr_old = np.inner(Mr_old.conjugate(), AMr_old)
        if MrTAMr_old != 0.:
            beta = np.inner(Mr.conjugate(), AMr) / MrTAMr_old
        else:
            if show:
                print("CR: Linear solve not converged due to BREAKDOWN "
                      f"iterations {iter+1}")
            return (postprocess(x), -1)

        # Update the search direction
        p = Mr + beta * p
        Ap = AMr + beta * Ap

    if show:
        print("CR: Linear solve not converged due to reach MAXIT "
              f"iterations {maxiter}")
    return (postprocess(x), maxiter)
