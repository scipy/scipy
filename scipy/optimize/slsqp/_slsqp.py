import numpy as np
import scipy.linalg as la
from scipy.optimize import _nnls
from scipy.linalg.blas import dtrsm, dtrsv


def ldl_update(a, z, sigma):
    """
    Updates the LDL' factors of a matrix by a rank-one or dyadic product
    sigma*(z @ z.T). Arrays are modified in-place.

    Parameters
    ----------
    a : ndarray
        1D array, encoding the lower triangle of a positive definite array in
        LAPACK packed format with UPLO=True. The array is updated in place.
    z : ndarray
        1D update vector with compatible size that of 'a'.
    sigma : float
        Scalar factor that multiplies the rank-1 update.

    """
    if sigma == 0.0:
        return

    ij = 0
    t = 1.0/sigma
    n = len(z)

    w = z.copy()
    if sigma < 0.0:
        # Negative update
        for i in range(n):
            v = w[i]
            t += v*v/a[ij]
            for j in range(i+1, n):
                ij += 1
                w[j] -= v*a[ij]
            ij += 1
        if t >= 0.0:
            t = np.finfo(a.dtype).eps / sigma

        for i in range(n):
            j = n - i - 1
            ij -= i + 1
            u = w[j]
            w[j] = t
            t -= u*u / a[ij]

    for i in range(n):
        v = z[i]
        delta = v / a[ij]
        tp = w[i] if sigma < 0.0 else t + delta*v
        alpha = tp / t
        a[ij] *= alpha
        if i == n - 1:
            return
        beta = delta / tp
        if alpha <= 4.0:
            for j in range(i+1, n):
                ij += 1
                z[j] -= v*a[ij]
                a[ij] += beta*z[j]
        else:
            gamma = t / tp
            for j in range(i+1, n):
                ij += 1
                u = a[ij]
                a[ij] = gamma*u + beta*z[j]
                z[j] -= v*u
        ij += 1
        t = tp
    return


def ldp(G, h):
    """Solves the constrained problem, also known as "least distance problem"
               1
        min   --- |x|**2  subject to Gx >= h
         x     2

    via solving the dual problem using Nonnegative Least Squares (nnls)
    """
    E = np.empty([G.shape[1]+1, G.shape[0]], dtype=G.dtype)
    E[:-1, :] = G.T
    E[-1, :] = h
    # E = np.vstack([G.T, h])
    f = np.zeros([E.shape[0]], dtype=E.dtype)
    f[-1] = 1
    u, unorm, mode = _nnls(E, f)
    if mode == -1:
        return np.zeros_like(u), u, 3
    if unorm == 0.:
        return np.zeros_like(u), u, 4
    r = E @ u
    r -= f
    if np.abs(r[-1]) < np.finfo(np.float64).eps:
        return np.zeros_like(u), u, 4

    return -r[:-1]/r[-1], -u/r[-1], 1


def lsi(A, b, G, h):
    """Solve inequality constrained least squares problem

        min |Ax - b|  subject to Gx >= h
         x

    Assumes A to be full rank. Also b, h is expected to be 1D
    """
    eps = np.spacing(1.)
    Q, R = la.qr(A)
    m, n = A.shape
    # Only deal with full rank problems
    if (np.abs(np.diag(R)) < m*n*eps).any():
        return np.zeros(n, dtype=A.dtype), np.zeros(n, dtype=A.dtype), 5

    bb = Q.T @ b
    # Define x = inv(R) (z - ff)
    # Transform G and h:  X R = G -> R.T X.T = G.T
    GG = dtrsm(1., R[:n, :], G, side=1)
    hh = h - GG @ bb[:n]

    z, mult, mode = ldp(GG, hh)
    if mode != 1:
        return np.zeros([R.shape[1]], dtype=R.dtype), mult, mode

    x = dtrsm(1., R[:n, :], z + bb[:n])
    return x, mult, 1


def lsei(A, b, E, f, G, h):
    """Solve inequality and equality constrained least squares problem

        min |Ax - b|  subject to Gx >= h,  Ex = f
         x

    E for (e)quality, and G for "(g)reater than" inequality constraints.

    Converts the problem into only inequality constrained problem and
    uses lsi()

    Parameters
    ----------
    A, E, G: ndarray
        2D arrays for the left hand side data of the optimization problem. The
        second dimension should be all equal. E cannot have more rows than the
        number of states.
    b, f, h: ndarray
        1D arrays for the right hand side data of the optimization problem,
        sizes compatible with A, E, G.

    Returns
    -------
    x: ndarray
        1D array that contains the solution or the all zero entries depending
        on the success.
    mult : ndarray
        If applicable (E and/or G supplied), returns the multiplier vector.
    mode: int
        If 1 then problem is solved successfully otherwise some error occured,
        error codes following the former fortran code.
          2 : Size mismatches (returns dummy empty arrays to save memory)
          5 : Coefficient array A is not full rank.
          6 : Equality constraint array E is not full rank.
          7 : NNLS problem failed to find a solution or if no constraints given
              A is not full rank, up to lstsq precision.
    """
    eps = np.spacing(1.)
    mineq = G.shape[0]
    meq, n = E.shape
    if meq > n:
        return (np.array([]), np.array([]), 2)

    if meq == 0 and mineq == 0:
        xx, _, krank, _ = la.lstsq(A, b)
        if krank != n:
            return xx, np.array([]), 7
        # No constraints just least squares, no multipliers
        return (xx, np.array([]), 1)

    elif meq == 0 and mineq > 0:
        # No equality constraints, it is an LSI problem
        return lsi(A, b, G, h)

    nd = n - meq

    # We right triangulize E and, apply Q.T to A and, if given, to G;
    #
    #     [E]         [0  Et]    [  ]
    #     [A] @ Q.T = [A2 A1] := [AA]
    #     [G]         [G2 G1]    [GG]
    #
    # the indices are swapped due to 'scipy.linalg.rq' returning Et on the
    # right end, as upper triangular.

    Et, Q = la.rq(E)

    # Remove the zero part coming from RQ decomposition
    Et = Et[:, nd:]
    AA = A @ Q.T
    if mineq > 0:
        GG = G @ Q.T

    # Only deal with full rank problems
    if (np.abs(np.diag(Et)) < meq*n*eps).any():
        return (np.zeros(len(b), dtype=A.dtype),
                np.zeros(E.shape[0], dtype=E.dtype),
                6)

    # Solve Ex = f
    x1 = dtrsv(Et, f)

    if mineq == 0:  # No inequality constraints, basic least squares
        x2, _, _, _ = la.lstsq(AA[:, :nd], b - AA[:, nd:]@x1)
    else:

        x2, mult_ls, mode = lsi(AA[:, :nd], b - AA[:, nd:] @ x1,
                                GG[:, :nd], h - GG[:, nd:] @ x1)
        # If failed
        if mode != 1:
            return (x2, mult_ls, mode)

    # Solution of the original problem and the multipliers
    xx = np.hstack([x2, x1])
    mult = (AA @ xx)
    mult -= b
    mult = AA[:, nd:].T @ mult
    if mineq != 0:
        mult -= GG[:, nd:].T @ mult_ls
    dtrsv(Et, mult, trans=1, overwrite_x=True)
    if mineq != 0:
        # combine the multipliers of eq and ineq constraints
        mult = np.hstack([mult, mult_ls])
    # Note: it is [x2 x1] due to RQ swap
    xx = Q.T @ xx

    return xx, mult, 1
