# Author : Rondall E. Jones, August 2020

from math import sqrt
import numpy as np
from numpy import atleast_1d, atleast_2d
import scipy.linalg.decomp_svd
from scipy.linalg.decomp import _asarray_validated
from scipy.linalg.misc import LinAlgError, _datacopied, LinAlgWarning

__all__ = ['autosolve', 'autosolve_nonneg']


def two_norm(x):
    """ Computes two-norm of 1-D array x without overflow"""
    ln = len(x)
    if ln == 0:
        return 0.0
    y = abs(x)
    mx = max(y)
    if mx == 0.0:
        return 0.0
    sm = 0.0
    for i in range(0, ln):
        sm += (x[i] / mx) ** 2
    return mx * sqrt(sm)


def rms(x):
    """ Computes root-mean-square of 1-D array x without overflow"""
    s = two_norm(x)
    return s / sqrt(len(x))


def compute_mov_sums(g, w, m):
    """ Computes all moving sums of width w in g[0] to g[m-1]. """
    numsums = m - w + 1
    sums = np.zeros(numsums)
    for i in range(0, numsums):
        s = 0.0
        for j in range(i, i + w):
            s += g[j]
        sums[i] = s
        return sums


def splita(mg, g):
    """ Determines a usable rank based on large rise in Picard Vector g"""
    # initialize
    sensitivity = g[0]
    small = sensitivity
    local = sensitivity
    urank = 1
    for i in range(1, mg):  # start with 2nd row; (i=1; i<mg; i++)
        sensitivity = g[i]
        if sensitivity > 15.0 * small and sensitivity > local:
            break
        if sensitivity < small:
            small = small + 0.40 * (sensitivity - small)
        else:
            small = small + 0.02 * (sensitivity - small)
        local = local + 0.40 * (sensitivity - local)
        urank = i + 1
    return urank


def splitb(mg, g):
    """ Determines a usable rank based on small rise in Picard Vector g"""
    gg = np.zeros(mg)
    for i in range(0, mg):
        gg[i] = g[i] * g[i]
    w = min(int((mg + 3) / 4), 6)
    sums = compute_mov_sums(gg, w, mg)
    ilow = np.where(sums == min(sums))[0][0]
    # estimate a nominal value that should see a big rise later
    sum = 0.0
    for i in range(0, ilow):
        sum += abs(gg[i])
    gnom = sum / float(ilow + 1)
    # see if the moving average ever gets much larger
    bad = 10.0 * gnom
    ibad = 0
    for i in range(ilow + 1, mg - w + 1):
        if sums[i] > bad:
            ibad = i
            break
    if ibad <= 0:
        urank = mg  # leave urank alone
    else:
        urank = ibad + w - 1
    return urank


def rmslambdah(A, b, U, S, Vt, ur, lamb):
    """ Computes solution to Ax=b its residual, given rank & Tik.lambda."""
    mn = S.shape[0]
    ps = np.zeros(mn)
    for i in range(0, ur):
        ps[i] = 1.0 / (S[i] + lamb ** 2 / S[i]) if S[i] > 0.0 else 0.0
    for i in range(ur, mn):
        ps[i] = 0.0
    # best to do multiplies from right end....
    xa = np.transpose(Vt) @ (np.diag(ps) @ (np.transpose(U) @ b))
    res = b - A @ xa
    r = rms(res)
    return xa, r


def discrep(A, b, U, S, Vt, ur, mysigma):
    """ Computes Tikhonov's lambda using b's estimated RMS error, mysigma"""
    lo = 0.0  # for minimum achievable residual
    hi = 0.33 * float(S[0])  # for ridiculously large residual
    lamb = 0.0
    # bisect until we get the residual we want...but quit eventually
    for k in range(0, 50):
        lamb = (lo + hi) * 0.5
        xa, check = rmslambdah(A, b, U, S, Vt, ur, lamb)
        if abs(check - mysigma) < 0.0001 * mysigma:
            break  # close enough!
        if check > mysigma:
            hi = lamb
        else:
            lo = lamb
    return lamb


def autosolve(A, b):
    """Solves the linear system of equation, A*x = b, for any shape matrix.

    The system can be underdetermined, square, or over-determined.
    That is, A(m,n) can be such that m < n, m = n, or m > n.
    b is a vector of length m.
    This solver automatically detects if the system is ill-conditioned or not.
    Then
     -- If the equations are consistent then the solution will be
        exact within round-off error.
     -- If the equations are inconsistent then the the solution will be
        by least-squares. That is, it solves ``min ||b - Ax||_2``.
     -- If the equations are inconsistent and diagnosable as ill-conditioned
        using the principles of the first reference below, the system will be
        automatically regularized and the residual will be larger than minimum.

    Parameters
    ----------
    A : (M, N) array_like "Coefficient" matrix of type float.
    b : (M) 1-D array_like Ordinate or "dependent variable" values, type float.

    Returns
    -------
    x : (N) ndarray of float.
        The solution, as explained above.
        To return only this solution, call x = autosolve(A,b)[0]
    ur: int
        The "usable rank", which is usually smaller than the conventional rank.
    sigma : float
        Estimated Right hand Side root-mean-square error
    lambda : float
         The estimated Tikhonov regularization parameter, lambda.

    Raises
    ------
    If the row size of the input matrix A is not the same as the length of the
    right hand side vector, b, then autosolve raises this:

       LinAlgError: Inconsistent array sizes.

    Autosolve is protected from virtually all calculation faults.
    However, Autoreg calls scipy's package's SVD algorithm, which may raise:

       LinAlgError: Computation does not converge.

    Examples
    --------
    See documentation of scipy's lstsq for a standard least-square example.
    Autosolve will behave like lstsq when the system is well conditioned.
    Here is a tiny example of an ill-conditioned system as handled by autoreg.

    x + y = 2
    x + 1.01 y =3

    Then A = array([[ 1.,  1.],
                   [ 1.,  1.01.]])
    and b = array([2.0, 3.0])

    Then standard solvers like lstsq will return:
    x = [-98. , 100.]

    But autosolve() will see the violation of the Picard Condition and return
    x = [1.12216 , 1.12779]

    Notes:
    -----
    1. When the system is ill-conditioned, the process works best when the rows
       of A are scaled so that the elements of b have similar estimated errors.
    2. Autosolve occasionally may produce a smoother (i.e., more regularized)
       solution than desired. In this case please try scipy routine lsmr.
    3. With any linear equation solver, check that the solution is reasonable.
       In particular, you check the residual vector, A*x - b.
    4. Autosolve neither needs or accepts optional parameters such as iteration
       limits, error estimates, variable bounds, condition number limits, etc.
       It also does not return any error flags as there are no error states.
       As long as the SVD converges (and svd failure is remarkably rare)
       then autosolve and autosolve_nonneg will complete normally.

    References
    ----------
    About the Picard Condition: "The discrete picard condition for discrete
    ill-posed problems", Per Christian Hansen, 1990.
    https://link.springer.com/article/10.1007/BF01933214

    About Nonnegative solutions: "Solving Least Squares Problems",
    by Charles L. Lawson and Richard J. Hanson. Prentice-Hall 1974

    About our algorithm: "Solving Linear Algebraic Systems Arising in the
    Solution of Integral Equations of the First Kind",
    Dissertation by Rondall E. Jones, 1985, U. of N.M.
    Advisor: Cleve B. Moler, creator of MatLab and co-founder of MathWorks.

    For algorithmic details: http://www.rejones7.net/REJTRIX/indextutorials.htm

    For further examples: See the appropriate links in http://www.rejones7.net/
    """
    A = atleast_2d(_asarray_validated(A, check_finite=True))
    b = atleast_1d(_asarray_validated(b, check_finite=True))
    m = A.shape[0]
    n = A.shape[1]
    mn = min(m, n)
    if b.shape[0]!= m: raise LinAlgError('Inconsistent array sizes.')
    if np.count_nonzero(A) == 0:
        return 0.0, 0, 0.0, 0.0
    if np.count_nonzero(b) == 0:
        return 0.0, 0, 0.0, 0.0

    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    beta = np.transpose(U) @ b
    # compute contributions to norm of solution
    k = 0  # rank of A so far
    g = np.zeros(mn)
    sense = 0.0
    si = 0.0
    for i in range(0, mn):
        si = S[i]
        if si == 0.0:
            break
        sense = beta[i] / si
        if sense < 0.0:
            sense = -sense
        g[i] = sense
        k = i + 1
    if k <= 0:
        return np.zeros(n)  # zero system
    # two-stage search for break in Picard Condition Vector
    ura = splita(k, g)
    ur = splitb(ura, g)
    if ur >= mn:
        # problem is not ill-conditioned
        x, check = rmslambdah(A, b, U, S, Vt, ur, 0.0)
        sigma = 0.0
        lambdah = 0.0
    else:
        # from urb, determine sigma, then lambda, then solution
        Utb = np.transpose(U) @ b
        sigma = rms(Utb[ur:mn])
        lambdah = discrep(A, b, U, S, Vt, ur, sigma)
        x, check = rmslambdah(A, b, U, S, Vt, ur, lambdah)
    return x, ur, sigma, lambdah


def autosolve_nonneg(A, b):
    """ Computes a nonnegative solution of A*x = b, for A of any shape.

    Autosolve_nonneg uses autosolve, above, and iteratively zeroes
    variables that violates the nonnegativity constraint.
    
    Parameters
    ----------
    A : (M, N) array_like "Coefficient" matrix of float.
    b : (M) 1-D array_like Ordinate or "dependent variable" values of float.

    Returns
    -------
    x : (N) ndarray of float.
        The solution, as explained above.
        To return only this solution, call x = autosolve(A,b)[0]
    """
    A = atleast_2d(_asarray_validated(A, check_finite=True))
    b = atleast_1d(_asarray_validated(b, check_finite=True))
    m = A.shape[0]
    n = A.shape[1]
    mn = min(m, n)
    if b.shape[0]!= m: raise LinAlgError('Inconsistent array sizes.')
    x = np.zeros(n)
    if np.count_nonzero(A) == 0:
        return 0.0, 0, 0.0, 0.0
    if np.count_nonzero(b) == 0:
        return 0.0, 0, 0.0, 0.0
        
    x, ur, sigma, lambdah = autosolve(A, b)
    # see if unconstrained solution is already non-negative
    if min(x) >= 0.0:
        return x
    # the approach here is to actually delete columns, for SVD speed,
    # rather than just zero out columns and thereby complicate the SVD.
    C = A
    cols = [0] * n  # list of active column numbers
    for i in range(1, n):
        cols[i] = i
    xt = x
    nn = n
    for i in range(1, nn):
        # choose a column to zero
        p = -1
        worst = 0.0
        for j in range(0, nn):
            if xt[j] < worst:
                p = j
                worst = xt[p]
        if p < 0:
            break
        # remove column p and resolve
        C = np.delete(C, p, 1)
        cols.pop(p)
        nn -= 1
        U, S, Vt = np.linalg.svd(C, full_matrices=False)
        ms = len(S)
        ps = np.zeros(ms)
        for i in range(0, ms):
            ps[i] = 1.0 / (S[i] + lambdah ** 2 / S[i]) \
                if S[i] > 0.0 else 0.0
        xt = np.transpose(Vt) @ (np.diag(ps) @ (np.transpose(U) @ b))

    # rebuild full solution vector
    if xt[0] < 0.0: xt[0]=0.0 #degenerate case 
    xn = np.zeros(n)
    for j in range(0, nn):
        xn[int(cols[j])] = xt[j]
    return xn
