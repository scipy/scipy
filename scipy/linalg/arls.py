from math import sqrt
import numpy as np
from numpy import atleast_1d, atleast_2d
from scipy.linalg import qr, solve_triangular
from scipy.linalg.decomp import _asarray_validated
from scipy.linalg.misc import LinAlgError

# Author : Rondall E. Jones, October 2020


def checkAb(A, b, maxcols):
    if len(A.shape) != 2:
        raise LinAlgError("Input array should be 2-D.")
    m, n = A.shape
    if m == 0 or n == 0:
        raise LinAlgError("Matrix is empty.")
    if b.shape[0] != m:
        raise LinAlgError(
            "Matrix and RHS do not have the same number of rows.")
    if maxcols < 2 and len(b.shape) > 1:
        raise LinAlgError(
            "Right hand side must have only one column.")
    return


def mynorm(x):
    return sqrt(np.dot(x, x))


def myrms(x):
    return sqrt(np.dot(x, x)) / sqrt(len(x))


def decide_width(mg):
    if mg < 3:
        return 1
    elif mg <= 8:  # 4 spans
        return 2
    elif mg <= 20:  # 5 spans
        return 4
    elif mg <= 36:  # 6 spans
        return 6
    elif mg <= 64:  # 8 spans
        return 8
    elif mg <= 100:  # 10 spans
        return 10
    else:
        w = int(mg / 10)
        return 2 * int(w / 2)  # 10 spans


def splita(mg, g):
    """ Determines a usable rank based on large rise in Picard Vector"""
    # initialize
    if mg < 2:
        return mg
    w = decide_width(mg)
    sensitivity = g[0]
    small = sensitivity
    local = sensitivity
    urank = 1
    for i in range(1, mg):
        sensitivity = g[i]
        if i >= w and sensitivity > 25.0 * small and sensitivity > local:
            break
        if sensitivity < small:
            small = small + 0.40 * (sensitivity - small)
        else:
            small = small + 0.10 * (sensitivity - small)
        local = local + 0.40 * (sensitivity - local)
        urank = i + 1
    return urank


def compute_mov_sums(g, w, m):
    numsums = m - w + 1
    sums = np.zeros(numsums)
    for i in range(0, numsums):
        s = 0.0
        for j in range(i, i + w):
            s += g[j]
        sums[i] = s
    return sums


def splitb(mg, g):
    """ Determines a usable rank based on modest rise in Picard Vector
    after the low point in the PCV."""
    w = decide_width(mg)
    if w < 2:
        return mg  # splitb needs w>=2 to be reliable

    # magnify any divergence by squaring
    gg = np.zeros(mg)
    for i in range(0, mg):
        gg[i] = g[i] * g[i]

    # ignore dropouts
    for i in range(1, mg - 1):
        if gg[i] < 0.2 * gg[i - 1] and gg[i] < 0.2 * gg[i + 1]:
            gg[i] = 0.5 * min(gg[i - 1], gg[i + 1])

    # choose breakpoint as multiple of lowest moving average
    sums = compute_mov_sums(gg, w, mg)
    ilow = np.where(sums == min(sums))[0][0]
    bad = 20.0 * sums[ilow]

    # look for unexpected rise
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
    """ Computes a regularized solution to Ax=b, given the usable rank
    and the Tikhonov lambda value."""
    mn = S.shape[0]
    ps = np.zeros(mn)
    for i in range(0, ur):
        ps[i] = 1.0 / (S[i] + lamb ** 2 / S[i]) if S[i] > 0.0 else 0.0
    for i in range(ur, mn):
        ps[i] = 0.0
    # best to do multiplies from right end....
    xa = np.transpose(Vt) @ (np.diag(ps) @ (np.transpose(U) @ b))
    res = b - A @ xa
    r = myrms(res)
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
        if abs(check - mysigma) < 0.0000001 * mysigma:
            break  # close enough!
        if check > mysigma:
            hi = lamb
        else:
            lo = lamb
    return lamb


def arlsusv(A, b, U, S, Vt):
    """ core solver when SVD is already available """
    if np.count_nonzero(A) == 0 or np.count_nonzero(b) == 0:
        return np.zeros(A.shape[1])
    m, n = A.shape
    mn = min(m, n)
    # compute contributions to norm of solution
    beta = np.transpose(U) @ b
    k = 0
    g = np.zeros(mn)
    sense = 0.0
    si = 0.0
    cond = max(A.shape) * np.spacing(A.real.dtype.type(1))
    eps = S[0] * cond
    for i in range(0, mn):
        si = S[i]
        if si <= eps:
            break
        sense = beta[i] / si
        if sense < 0.0:
            sense = -sense
        g[i] = sense
        k = i + 1
    nr = k  # traditional numeric rank
    if k <= 0:
        return np.zeros(n)  # failsave check

    # two-stage search for divergence in Picard Condition Vector
    ura = splita(k, g)
    urb = splitb(ura, g)
    ur = min(ura, urb)
    if ur >= mn:
        # problem is not ill-conditioned
        x, check = rmslambdah(A, b, U, S, Vt, ur, 0.0)
        sigma = 0.0
        lambdah = 0.0
    else:
        # from ur, determine sigma
        Utb = np.transpose(U) @ b
        sigma = myrms(Utb[ur:mn])
        # from sigma, determine lambda
        lambdah = discrep(A, b, U, S, Vt, ur, sigma)
        # from lambda, determine solution
        x, check = rmslambdah(A, b, U, S, Vt, ur, lambdah)
    return x, nr, ur, sigma, lambdah


def arls(A, b):
    """Solves the linear system of equation, Ax = b, for any shape matrix.

    The system can be underdetermined, square, or over-determined.
    That is, A(m,n) can be such that m < n, m = n, or m > n.
    Argument b is a matrix of size(n,p) of p right-hand-side columns.
    This solver automatically detects if each system is ill-conditioned or not.

    Then
     -- If the equations are consistent then the solution will usually be
        exact within round-off error.
     -- If the equations are inconsistent then the the solution will be
        by least-squares. That is, it solves ``min ||b - Ax||_2``.
     -- If the equations are inconsistent and diagnosable as ill-conditioned
        using the principles of the first reference below, the system will be
        automatically regularized and the residual will be larger than minimum.
     -- If either A or b is all zeros then the solution will be all zeros.

    Parameters
    ----------
    A : (m, n) array_like
        Coefficient matrix
    b : (m, p) array_like
        Set of columns of dependent variables.

    Returns
    -------
    x : (n, p) array_like set of columns, type float.
        Each column will be the solution corresponding to a column of b.
    nr : int
        The Numerical Rank of A.
    mur : int
        The Minimum Usable Rank seen in solving all the problems.
        Note that "numerical rank" is an attribute of a matrix
        but the "usable rank" that arls computes is an attribute
        of the problem, Ax=b.
    msigma : float
        The largest estimated right-hand-side root-mean-square error seen.
    mlambda : float
        The largest estimated Tikhonov regularization parameter seen.

    Raises
    ------
    LinAlgError
        If A is not 2-D.
        If A is empty.
        If A and b do not have the same row size.
        If b has more than one column.
        If SCIPY's SVD() does not converge.

    Examples
    --------
    Arls() will behave like any good least-squares solver when the system
    is well conditioned.
    Here is a tiny example of an ill-conditioned system as handled by arls(),

       x + y = 2
       x + 1.01 y =3

    Then A = array([[ 1., 1.],
                    [ 1., 1.01.]])
    and  b = array([2.0, 3.0])

    Then standard solvers will return:
       x = [-98. , 100.]

    But arls() will see the violation of the Picard Condition and return
       x = [1.12216 , 1.12779]

    Notes:
    -----
    1. When the system is ill-conditioned, the process works best when the rows
       of A are scaled so that the elements of b have similar estimated errors.
    2. Arls() occasionally may produce a smoother (i.e., more regularized)
       solution than desired. In this case please try scipy routine lsmr.
    3. With any linear equation solver, check that the solution is reasonable.
       In particular, you should check the residual vector, Ax - b.
    4. Arls() neither needs nor accepts optional parameters such as iteration
       limits, error estimates, variable bounds, condition number limits, etc.
       It also does not return any error flags as there are no error states.
       As long as the SVD converges (and SVD failure is remarkably rare)
       then arls() and other routines in this package will complete normally.
    5. Arls()'s intent (and the intent of all routines in this module)
       is to find a reasonable solution even in the midst of excessive
       inaccuracy, ill-conditioning, singularities, duplicated data, etc.
       Its performance is often very like that of lsmr, but from a completely
       different approach.
    6. In view of note 5, arls() is not appropriate for situations
       where the requirements are more for high accuracy rather than
       robustness. So, we assume, in the coding, where needed, that no data
       needs to be considered more accurate than 8 significant figures.

    References
    ----------
    The auto-regularization algorithm in this software arose from the research
    for my dissertation, "Solving Linear Algebraic Systems Arising in the
    Solution of Integral Equations of the First Kind", University of
    New Mexico, Albuquerque, NM, 1985.

    Many thanks to Cleve B. Moler, MatLab creater and co-founder of MathWorks
    for his energy and insights in guiding my dissertation research.

    My thanks also to Richard Hanson (deceased), co-author of the classic
    "Solving Least Squares Problems", co-creater of BLAS, and co-advisor
    for the last year of my dissertation work.

    My thanks also to Per Christian Hansen for reviewing an early version
    of this software which resulted in my creating the crucial two-phase
    matrix "split" process.

    And my thanks to the central computing department at Sandia National Labs
    where I had the opportunity to contribute to the "MATHLIB" Fortran library
    in the 1960's and 1970's. MATLIB evolved into part of the
    Fortran "SLATEC" library.
    See http://www.netlib.org/slatec/src/

    For for a short presentation on the Picard Condition which is at the heart
    of this package's algorithms, please see
    See www.rejones7.net/Arls/What_Is_The_Picard_Condition.htm

    For a complete description, see "The Discrete Picard Condition for Discrete
    Ill-posed Problems", Per Christian Hansen, 1990.
    See link.springer.com/article/10.1007/BF01933214

    For discussion of incorporating equality and inequality constraints
    (including nonnegativity) in solving linear algebraic problems, see
    "Solving Least Squares Problems", by Charles L. Lawson and
    Richard J. Hanson, Prentice-Hall 1974.
    My implementation of these features has evolved somewhat
    from that fine book, but is based on those algorithms.

    Rondall E. Jones, Ph.D.
    rejones7@msn.com
    """
    A = atleast_2d(_asarray_validated(A, check_finite=True))
    b = atleast_1d(_asarray_validated(b, check_finite=True))
    if np.count_nonzero(A) == 0 or np.count_nonzero(b) == 0:
        return np.zeros(A.shape[1]), 0, 0, 0.0, 0.0
    checkAb(A, b, 2)
    m, n = A.shape
    nrhs = 1
    if len(b.shape) == 2:
        nrhs = b.shape[1]

    U, S, Vt = np.linalg.svd(A, full_matrices=False)

    # one right hand side
    if nrhs == 1:
        return arlsusv(A, b, U, S, Vt)

    # multiple right hand sides
    xx = np.zeros((n, nrhs))
    nr = min(A.shape)   # track minimum numeric rank
    mur = nr            # track minimum usable rank
    msigma = 0.0        # track maximum estimated RHS error
    mlamda = 0.0        # track maximum Tikhonov parameter
    for p in range(0, nrhs):
        xx[:, p], nr, ur, sigma, lambdah = arlsusv(A, b[:, p], U, S, Vt)
        mur = min(mur, ur)
        msigma = max(msigma, sigma)
        mlamda = max(mlamda, lambdah)
    return xx, nr, mur, msigma, mlamda


def strange(A, b):
    """determines whether the A*x=b is likely ill-conditioned."""
    m, n = A.shape
    epsilon = 10.0 * max(m, n) * np.spacing(b.real.dtype.type(1))
    rn = np.zeros(m)
    for k in range(0, m):
        rn[k] = mynorm(A[k, :])
    tol = max(rn) * epsilon
    odd = False
    for k in range(0, m):
        if rn[k] < tol:
            rn[k] = tol
            odd = True
    bb = abs(b)/rn
    maxnorm = mynorm(bb) * 1.5
    return odd, maxnorm


def arlsqr(A, b):
    """Solves the linear system of equation, Ax = b, for any shape matrix.

    Arlsqr() is called exactly like arls() above, and the returns are
    exactly the same.

    The difference is that if b is only one column then arlsqr()
    first performs a quick assessment of the system to see if it appears
    to be ill-conditioned. If not, qr() is called for a solution.
    If the resulting solution from qr() appears to be good, then it
    is returned. Otherwise, arls() continues itself to to get a solution
    based on the SVD.

    The purpose of this difference is to execute signficantly faster than
    arlsqr() does if the problems given to it are usually well behaved.
    """
    A = atleast_2d(_asarray_validated(A, check_finite=True))
    b = atleast_1d(_asarray_validated(b, check_finite=True))
    nz = np.count_nonzero(A)
    if nz == 0 or np.count_nonzero(b) == 0:
        return np.zeros(A.shape[1]), 0, 0, 0.0, 0.0
    checkAb(A, b, 2)
    m, n = A.shape
    nr = min(m, n)
    # is b a single column?
    if len(b.shape) > 1:
        return arls(A, b)
    odd, maxnorm = strange(A, b)
    if not odd or nz < int((m*n)/2):  # max norm tests fails for sparse-like
        if (m >= n):
            Q, R = qr(A, mode='economic')
            Qb = Q.T @ b
            try:
                x = solve_triangular(R, Qb)
            except LinAlgError:
                return arls(A, b)
        else:
            Q, R = qr(A.T, mode='economic')
            try:
                y = solve_triangular(R, b, 'T')
            except LinAlgError:
                return arls(A, b)
            x = Q @ y
        if mynorm(x) < maxnorm or nz < int((m*n)/2):
            return x, nr, nr, 0.0, 0.0
    return arls(A, b)


def find_max_row_norm(A):
    """ determine max row norm of A """
    m = A.shape[0]
    rnmax = 0.0
    for i in range(0, m):
        rn = mynorm(A[i, :])
        if rn > rnmax:
            rnmax = rn
    return rnmax


def find_max_sense(E, f):
    """ find the row of Ex=f which his the highest ratio of f[i]
    to the norm of the row. """
    snmax = -1.0
    ibest = 0  # default
    m = E.shape[0]
    for i in range(0, m):
        rn = mynorm(E[i, :])
        if rn > 0.0:
            s = abs(f[i]) / rn
            if s > snmax:
                snmax = s
                ibest = i
    return ibest


def prepeq(E, f, neglect):
    """ a utility routine for arlseq() below that prepares the equality
    constraints for use"""
    E = atleast_2d(_asarray_validated(E, check_finite=True))
    f = atleast_1d(_asarray_validated(f, check_finite=True))
    EE = E.copy()
    ff = f.copy()
    m, n = EE.shape
    for i in range(0, m):
        # determine new best row and put it next
        if i == 0:
            imax = find_max_sense(EE, ff)
        else:
            rnmax = -1.0
            imax = -1
            for k in range(i, m):
                rn = mynorm(EE[k, :])
                if rn > rnmax:
                    rnmax = rn
                    imax = k
        EE[[i, imax], :] = EE[[imax, i], :]
        ff[[i, imax]] = ff[[imax, i]]

        # normalize
        rin = mynorm(EE[i, :])
        if rin > 0.0:
            EE[i, :] /= rin
            ff[i] /= rin
        else:
            ff[i] = 0.0

        # subtract projections onto EE[i,:]
        for k in range(i + 1, m):
            d = np.dot(EE[k, :], EE[i, :])
            EE[k, :] -= d * EE[i, :]
            ff[k] -= d * ff[i]

    # reject ill-conditioned rows
    if m > 2:
        g = np.zeros(m)
        for k in range(0, m):
            g[k] = abs(ff[k])
        m1 = splita(m, g)
        mm = splitb(m1, g)
        if mm < m:
            EE = np.resize(EE, (mm, n))
            ff = np.resize(ff, mm)
    return EE, ff


def arlseq(A, b, E, f):
    """Solves the double linear system of equations

       Ax = b  (least squares)
       Ex = f  (exact)

    Both Ax=b and Ex=f system can be underdetermined, square,
    or over-determined. Arguments b and f must be single columns.

    Ex=f is treated as a set of equality constraints.
    These constraints are usually few in number and well behaved.
    But clearly the caller can easily provide equations in Ex=f that
    are impossible to satisfy as a group. For example, there could be
    one equation requiring x[0]=0, and another requiring x[0]=1.
    And, the solver must deal with there being redundant or other
    pathological situations within the E matrix.
    So the solution process will either solve each equation in Ex=f exactly
    (within roundoff) or if that is impossible, arlseq() will discard
    one or more equations until the remaining equations are solvable
    exactly (within roundoff).
    We will refer below to the solution of this reduced system as "xe".

    After Ex=f is processed as above, the rows of Ax=b will have their
    projections onto every row of Ex=f subtracted from them.
    We will call this reduced set of equations A'x = b'.
    (Thus, the rows of A' will all be orthogonal to the rows of E.)
    This reduced problem A'x = b', will then be solved with arls().
    We will refer to the solution of this system as "xt".

    The final solution will be x = xe + xt.
    (Since E and A' are mutually orthogonal matrices.)

    If A'x = b' is NOT ill-conditioned (or singular) then x = xe + xt
    will satisfy all the equality constraints that were not rejected.
    However, if A'x = b' is ill-conditioned, then xt will be a regularized
    solution, not an exact solution, and thus x = xe + xt will not in general
    satisfy the equality constraints exactly.

    Parameters
    ----------
    A : (m, n)  array_like "Coefficient" matrix, type float.
    b : (m)     array_like column of dependent variables, type float.
    E : (me, n) array_like "Coefficient" matrix, type float.
    f : (me)    array_like column of dependent variables, type float.

    Returns
    -------
    x : (n) array_like column, type float.
    nr : int
        The numerical rank of the matrix, A, after its projection onto the rows
        of E are subtracted.
    ur : int
        The "usable" rank of the "reduced" problem, Ax=b, after its projection
        onto the rows of Ex=f are subtracted.
        Note that "numerical rank" is an attribute of a matrix
        but the "usable rank" that arls computes is an attribute
        of the problem, Ax=b.
    sigma : float
        The estimated right-hand-side root-mean-square error.
    lambda : float
        The estimated Tikhonov regularization.

    Raises
    ------
    LinAlgError
        If A is not 2-D.
        If A is empty.
        If A and b do not have the same row size.
        If b has more than one column.
        If E is not 2-D.
        If E is empty.
        If E and f do not have the same row size.
        If f has more than one column.
        If A and E do not have the same number of columns.
        If SCIPY's SVD() does not converge.

    Examples
    --------
    Here is a tiny example of a problem which has an "unknown" amount
    of error in the right hand side, but for which the user knows that the
    correct SUM of the unknowns must be 3:

         x + 2 y = 5.3   (Least Squares)
       2 x + 3 y = 7.8
           x + y = 3     ( Exact )

    Then the arrays for arlseq are:

       A = array([[ 1.,  2.0],
                  [ 2.,  3.0]])
       b = array([5.3, 7.8])
       E = array([[ 1.0,  1.0]])
       f = array([3.0])

    Without using the equality constraint we are given here,
    standard solvers will return [x,y] = [-.3 , 2.8].
    Even arls() will return the same [x,y] = [-.3 , 2.8].
    The residual for this solution is [0.0 , 0.0] (within roundoff).
    But of course x + y = 2.5, not the 3.0 we really want.

    Arlsnn() could help here by disallowing presumably unacceptable
    negative values, producing [x,y] = [0. , 2.6].
    The residual for this solution is [-0.1 , 0.] which is of course
    an increase from zero, but this is natural since we have forced
    the solution away from being the "exact" result, for good reason.
    Note that x + y = 2.6, which is a little better.

    If we solve with arlseq(A,b,E,f) then we get [x,y] = [1.004, 1.996].
    This answer is close to the "correct" answer of [x,y] = [1.0 , 2.0]
    if the right hand side had been the correct [5.,8.] instead of [5.3,7.8].
    The residual for this solution is [-0.3 , 0.2] which is yet larger.
    Again, when adding constraints to the problem the residual
    typically increases, but the solution becomes more acceptable.
    Note that x + y = 3 exactly.

    Notes:
    -----
    See arls() above for notes and references.
    """
    A = atleast_2d(_asarray_validated(A, check_finite=True))
    b = atleast_1d(_asarray_validated(b, check_finite=True))
    if np.count_nonzero(A) == 0 or np.count_nonzero(b) == 0:
        return np.zeros(A.shape[1]), 0, 0, 0.0, 0.0
    AA = A.copy()
    bb = b.copy()
    checkAb(AA, bb, 1)
    m, n = AA.shape
    rnmax = find_max_row_norm(AA)
    neglect = rnmax * 0.00000001  # see Note 6. for arls()

    E = atleast_2d(_asarray_validated(E, check_finite=True))
    f = atleast_1d(_asarray_validated(f, check_finite=True))
    EE = E.copy()
    ff = f.copy()
    checkAb(EE, ff, 1)
    me, ne = EE.shape

    if n != ne:
        raise LinAlgError(
            "The two matrices do not have the same number of unknowns.")

    EE, ff = prepeq(EE, ff, neglect)
    mEE = EE.shape[0]

    # decouple AAx=bb from EEx=ff
    i = 0
    while i < m:
        for j in range(0, mEE):
            d = np.dot(AA[i, :], EE[j, :])
            AA[i, :] -= d * EE[j, :]
            bb[i] -= d * ff[j]
        nm = mynorm(AA[i, :])
        if nm < neglect:
            AA = np.delete(AA, i, 0)
            bb = np.delete(bb, i, 0)
            m = AA.shape[0]
        else:
            AA[i, :] = AA[i, :] / nm
            bb[i] = bb[i] / nm
            i += 1

    # final solution
    xe = np.transpose(EE) @ ff
    if AA.shape[0] > 0:
        xt, nr, ur, sigma, lambdah = arls(AA, bb)
        return xt + xe, nr, ur, sigma, lambdah
    else:
        return xe, 0, 0, 0., 0.


def arlsgt(A, b, G, h):
    """Solves the double linear system of equations
       Ax = b  (least squares)
       Gx >= h ("greater than" inequality constraints)
    Both Ax=b and Gx>=h can be underdetermined, square, or over-determined.
    Arguments b and h must be single columns.
    Arlsgt() uses arls(), above, as the core solver, and iteratively selects
    rows of Gx>=h to move to a growing list of equality constraints, choosing
    first whatever equation in Gx>=h most violates its requirement.

    Note that "less than" equations can be included by negating
    both sides of the equation, thus turning it into a "greater than".

    If either A or b is all zeros then the solution will be all zeros.

    Parameters
    ----------
    A : (m, n)  array_like "Coefficient" matrix, type float.
    b : (m)     array_like column of dependent variables, type float.
    G : (mg, n) array_like "Coefficient" matrix, type float.
    b : (mg)    array_like column of dependent variables, type float.

    Returns
    -------
    x : (n) array_like column, type float.
    nr : int
        The numerical rank of the matrix, A.
    ur : int
        The usable rank of the problem, Ax=b.
        Note that "numerical rank" is an attribute of a matrix
        but the "usable rank" that arls computes is an attribute
        of the problem, Ax=b.
    sigma : float
        The estimated right-hand-side root-mean-square error.
    lambda : float
        The estimated Tikhonov regularization.

    Raises
    ------
    LinAlgError
        If A is not 2-D.
        If A is empty.
        If A and b do not have the same row size.
        If b has more than one column.
        If G is not 2-D.
        If G is empty.
        If G and h do not have the same row size.
        If h has more than one column.
        If A and G do not have the same number of columns.
        If SCIPY's SVD() does not converge.

    Example
    -------
    Let A = [[1,1,1],
             [0,1,1],
             [1,0,1]]
    and b = [5.9, 5.0, 3.9]

    Then any least-squares solver would produce x = [0.9, 2., 3.]
    The residual for this solution is zero within roungoff.

    But if we happen to know that all the answers should be at least 1.0
    then we can add inequalites to insure that:
        x[0] >= 1
        x[1] >= 1
        x[2] >= 1

    This can be expressed in the matrix equation Gx>=h where
        G = [[1,0,0],
             [0,1,0],
             [0,0,1]]
        h =  [1,1,1]

    Then arlsgt(A,b,G,h) produces x = [1., 2.05, 2.9].
    The residual vector and its norm are then:
       res = [-0.05, 0.05, 0.0]
       norm(res) = 0.0707

    If the user had just adjusted the least-squares answer of [0.9, 2., 3.]
    to [1., 2., 3.] without re-solving then the residual vector
    and its norm would be
       res = [0.1, 0, 0.1]
       norm(res) = 0.141
    which is significantly larger.
    """
    A = atleast_2d(_asarray_validated(A, check_finite=True))
    b = atleast_1d(_asarray_validated(b, check_finite=True))
    if np.count_nonzero(A) == 0 or np.count_nonzero(b) == 0:
        return np.zeros(A.shape[1]), 0, 0, 0.0, 0.0
    checkAb(A, b, 1)
    m, n = A.shape

    G = atleast_2d(_asarray_validated(G, check_finite=True))
    h = atleast_1d(_asarray_validated(h, check_finite=True))
    GG = G.copy()
    hh = h.copy()
    checkAb(GG, hh, 1)
    mg, ng = GG.shape
    if n != ng:
        raise LinAlgError(
            "The two matrices do not have the same number of unknowns.")

    EE = []
    ff = []
    me = 0
    ne = 0

    # get initial solution... it might actually be right
    x, nr, ur, sigma, lambdah = arls(A, b)
    nx = mynorm(x)
    if nx <= 0.0:
        return np.zeros(n), 0, 0, 0., 0.

    # while constraints are not fully satisfied:
    while True:
        # assess state of inequalities
        p = -1
        mg = GG.shape[0]
        rhs = GG @ x
        worst = 0.0
        for i in range(0, mg):
            if rhs[i] < hh[i]:
                diff = hh[i] - rhs[i]
                if p < 0 or diff > worst:
                    p = i
                    worst = diff
        if p < 0:
            break

        # delete row from GGx=hh
        row = GG[p, :]
        rhs = hh[p]
        GG = np.delete(GG, p, 0)
        hh = np.delete(hh, p, 0)

        # add row to Ex>=f
        if me == 0:
            EE = np.zeros((1, ng))
            EE[0, :] = row
            ff = np.zeros(1)
            ff[0] = rhs
            me = 1
            ne = ng
        else:
            me += 1
            EE = np.resize(EE, (me, ne))
            EE[me - 1, :] = row[:]
            ff = np.resize(ff, me)
            ff[me - 1] = rhs
        # re-solve modified system
        x = arlseq(A, b, EE, ff)[0]
    return x, nr, ur, sigma, lambdah


def arlsnn(A, b):
    """Solves Ax = b in the least squares sense, with the solution
       constrained to be non-negative.

       For a nonpositive solution, use
          x = -arlsnn(A,-b)

    Parameters
    ----------
    A : (m, n) array_like "Coefficient" matrix, type float.
    b : (m) array_like column of dependent variables, type float.

    Returns
    -------
    x : (n) array_like column, type float.
    nr : int
        The numerical rank of the matrix, A.
    ur : int
        The usable rank of the problem, Ax=b.
        Note that "numerical rank" is an attribute of a matrix
        but the "usable rank" that arls computes is an attribute
        of the problem, Ax=b.
    sigma : float
        The estimated right-hand-side root-mean-square error.
    lambda : float
        The estimated Tikhonov regularization.

    Raises
    ------
    LinAlgError
        If A is not 2-D.
        If A is empty.
        If A and b do not have the same row size.
        If b has more than one column.
        If SCIPY's SVD() does not converge.

    Example
    -------
    Let A = [[2., 2., 1.],
             [2., 1., 0.],
             [1., 1., 0.]]
    and b =  [3.9, 3., 2.]
    Then any least-squares solver will produce
       x =  [1. ,1., -0.1]
    But arlsnn() produces  x = [ 1.0322, 0.9093, 0.].

    arlsnn() tries to produce a small residual for the final solution,
    while being based toward making the fewest changes feasible
    to the problem. Most older solvers try to minimize the residual
    at the expense of extra interference with the user's model.
    Arls, arlsgt, and arlsnn seek a better balance.
    """
    A = atleast_2d(_asarray_validated(A, check_finite=True))
    b = atleast_1d(_asarray_validated(b, check_finite=True))
    if np.count_nonzero(A) == 0 or np.count_nonzero(b) == 0:
        return np.zeros(A.shape[1]), 0, 0, 0.0, 0.0
    checkAb(A, b, 1)
    n = A.shape[1]
    G = np.eye(n)
    h = np.zeros(n)
    return arlsgt(A, b, G, h)
