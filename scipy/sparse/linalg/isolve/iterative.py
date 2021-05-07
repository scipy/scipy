"""Iterative methods for solving linear systems"""

__all__ = ['bicg','bicgstab','cg','cgs','gmres','qmr']

import warnings
import numpy as np

from . import _iterative

from scipy.sparse.linalg.interface import LinearOperator
from .utils import make_system
from scipy._lib._util import _aligned_zeros
from scipy._lib._threadsafety import non_reentrant

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}


# Part of the docstring common to all iterative solvers
common_doc1 = \
"""
Parameters
----------
A : {sparse matrix, dense matrix, LinearOperator}"""

common_doc2 = \
"""b : {array, matrix}
    Right hand side of the linear system. Has shape (N,) or (N,1).

Returns
-------
x : {array, matrix}
    The converged solution.
info : integer
    Provides convergence information:
        0  : successful exit
        >0 : convergence to tolerance not achieved, number of iterations
        <0 : illegal input or breakdown

Other Parameters
----------------
x0  : {array, matrix}
    Starting guess for the solution.
tol, atol : float, optional
    Tolerances for convergence, ``norm(residual) <= max(tol*norm(b), atol)``.
    The default for ``atol`` is ``'legacy'``, which emulates
    a different legacy behavior.

    .. warning::

       The default value for `atol` will be changed in a future release.
       For future compatibility, specify `atol` explicitly.
maxiter : integer
    Maximum number of iterations.  Iteration will stop after maxiter
    steps even if the specified tolerance has not been achieved.
M : {sparse matrix, dense matrix, LinearOperator}
    Preconditioner for A.  The preconditioner should approximate the
    inverse of A.  Effective preconditioning dramatically improves the
    rate of convergence, which implies that fewer iterations are needed
    to reach a given error tolerance.
callback : function
    User-supplied function to call after each iteration.  It is called
    as callback(xk), where xk is the current solution vector.

"""


def _stoptest(residual, atol):
    """
    Successful termination condition for the solvers.
    """
    resid = np.linalg.norm(residual)
    if resid <= atol:
        return resid, 1
    else:
        return resid, 0


def _get_atol(tol, atol, bnrm2, get_residual, routine_name):
    """
    Parse arguments for absolute tolerance in termination condition.

    Parameters
    ----------
    tol, atol : object
        The arguments passed into the solver routine by user.
    bnrm2 : float
        2-norm of the rhs vector.
    get_residual : callable
        Callable ``get_residual()`` that returns the initial value of
        the residual.
    routine_name : str
        Name of the routine.
    """

    if atol is None:
        warnings.warn("scipy.sparse.linalg.{name} called without specifying `atol`. "
                      "The default value will be changed in a future release. "
                      "For compatibility, specify a value for `atol` explicitly, e.g., "
                      "``{name}(..., atol=0)``, or to retain the old behavior "
                      "``{name}(..., atol='legacy')``".format(name=routine_name),
                      category=DeprecationWarning, stacklevel=4)
        atol = 'legacy'

    tol = float(tol)

    if atol == 'legacy':
        # emulate old legacy behavior
        resid = get_residual()
        if resid <= tol:
            return 'exit'
        if bnrm2 == 0:
            return tol
        else:
            return tol * float(bnrm2)
    else:
        return max(float(atol), tol * float(bnrm2))


def set_docstring(header, Ainfo, footer='', atol_default='0'):
    def combine(fn):
        fn.__doc__ = '\n'.join((header, common_doc1,
                                '    ' + Ainfo.replace('\n', '\n    '),
                                common_doc2, footer))
        return fn
    return combine


@set_docstring('Use BIConjugate Gradient iteration to solve ``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` and ``A^T x`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.',
               footer="""
               
               Examples
               --------
               >>> from scipy.sparse import csc_matrix
               >>> from scipy.sparse.linalg import bicg
               >>> A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
               >>> b = np.array([2, 4, -1], dtype=float)
               >>> x, exitCode = bicg(A, b)
               >>> print(exitCode)            # 0 indicates successful convergence
               0
               >>> np.allclose(A.dot(x), b)
               True
               
               """
               )
@non_reentrant()
def bicg(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None, atol=None):
    A,M,x,b,postprocess = make_system(A, M, x0, b)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec, rmatvec = A.matvec, A.rmatvec
    psolve, rpsolve = M.matvec, M.rmatvec
    ltr = _type_conv[x.dtype.char]
    revcom = getattr(_iterative, ltr + 'bicgrevcom')

    get_residual = lambda: np.linalg.norm(matvec(x) - b)
    atol = _get_atol(tol, atol, np.linalg.norm(b), get_residual, 'bicg')
    if atol == 'exit':
        return postprocess(x), 0

    resid = atol
    ndx1 = 1
    ndx2 = -1
    # Use _aligned_zeros to work around a f2py bug in Numpy 1.9.1
    work = _aligned_zeros(6*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            if callback is not None:
                callback(x)
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            work[slice2] *= sclr2
            work[slice2] += sclr1*rmatvec(work[slice1])
        elif (ijob == 3):
            work[slice1] = psolve(work[slice2])
        elif (ijob == 4):
            work[slice1] = rpsolve(work[slice2])
        elif (ijob == 5):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 6):
            if ftflag:
                info = -1
                ftflag = False
            resid, info = _stoptest(work[slice1], atol)
        ijob = 2

    if info > 0 and iter_ == maxiter and not (resid <= atol):
        # info isn't set appropriately otherwise
        info = iter_

    return postprocess(x), info


@set_docstring('Use BIConjugate Gradient STABilized iteration to solve '
               '``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def bicgstab(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None, atol=None):
    A, M, x, b, postprocess = make_system(A, M, x0, b)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom = getattr(_iterative, ltr + 'bicgstabrevcom')

    get_residual = lambda: np.linalg.norm(matvec(x) - b)
    atol = _get_atol(tol, atol, np.linalg.norm(b), get_residual, 'bicgstab')
    if atol == 'exit':
        return postprocess(x), 0

    resid = atol
    ndx1 = 1
    ndx2 = -1
    # Use _aligned_zeros to work around a f2py bug in Numpy 1.9.1
    work = _aligned_zeros(7*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            if callback is not None:
                callback(x)
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            work[slice1] = psolve(work[slice2])
        elif (ijob == 3):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            resid, info = _stoptest(work[slice1], atol)
        ijob = 2

    if info > 0 and iter_ == maxiter and not (resid <= atol):
        # info isn't set appropriately otherwise
        info = iter_

    return postprocess(x), info


@set_docstring('Use Conjugate Gradient iteration to solve ``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               '``A`` must represent a hermitian, positive definite matrix.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def cg(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None, atol=None):
    A, M, x, b, postprocess = make_system(A, M, x0, b)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom = getattr(_iterative, ltr + 'cgrevcom')

    get_residual = lambda: np.linalg.norm(matvec(x) - b)
    atol = _get_atol(tol, atol, np.linalg.norm(b), get_residual, 'cg')
    if atol == 'exit':
        return postprocess(x), 0

    resid = atol
    ndx1 = 1
    ndx2 = -1
    # Use _aligned_zeros to work around a f2py bug in Numpy 1.9.1
    work = _aligned_zeros(4*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            if callback is not None:
                callback(x)
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            work[slice1] = psolve(work[slice2])
        elif (ijob == 3):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            resid, info = _stoptest(work[slice1], atol)
            if info == 1 and iter_ > 1:
                # recompute residual and recheck, to avoid
                # accumulating rounding error
                work[slice1] = b - matvec(x)
                resid, info = _stoptest(work[slice1], atol)
        ijob = 2

    if info > 0 and iter_ == maxiter and not (resid <= atol):
        # info isn't set appropriately otherwise
        info = iter_

    return postprocess(x), info


@set_docstring('Use Conjugate Gradient Squared iteration to solve ``Ax = b``.',
               'The real-valued N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def cgs(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None, atol=None):
    A, M, x, b, postprocess = make_system(A, M, x0, b)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom = getattr(_iterative, ltr + 'cgsrevcom')

    get_residual = lambda: np.linalg.norm(matvec(x) - b)
    atol = _get_atol(tol, atol, np.linalg.norm(b), get_residual, 'cgs')
    if atol == 'exit':
        return postprocess(x), 0

    resid = atol
    ndx1 = 1
    ndx2 = -1
    # Use _aligned_zeros to work around a f2py bug in Numpy 1.9.1
    work = _aligned_zeros(7*n,dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            if callback is not None:
                callback(x)
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
        elif (ijob == 2):
            work[slice1] = psolve(work[slice2])
        elif (ijob == 3):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            resid, info = _stoptest(work[slice1], atol)
            if info == 1 and iter_ > 1:
                # recompute residual and recheck, to avoid
                # accumulating rounding error
                work[slice1] = b - matvec(x)
                resid, info = _stoptest(work[slice1], atol)
        ijob = 2

    if info == -10:
        # termination due to breakdown: check for convergence
        resid, ok = _stoptest(b - matvec(x), atol)
        if ok:
            info = 0

    if info > 0 and iter_ == maxiter and not (resid <= atol):
        # info isn't set appropriately otherwise
        info = iter_

    return postprocess(x), info


@non_reentrant()
def gmres(A, b, x0=None, tol=1e-5, restart=None, maxiter=None, M=None, callback=None,
          restrt=None, atol=None, callback_type=None, *, orth=None):
    """
    Use Generalized Minimal RESidual iteration to solve ``Ax = b``.

    Parameters
    ----------
    A : {sparse matrix, dense matrix, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, ``A`` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : {array, matrix}
        Right hand side of the linear system. Has shape (N,) or (N,1).

    Returns
    -------
    x : {array, matrix}
        The converged solution.
    info : int
        Provides convergence information:
          * 0  : successful exit
          * >0 : convergence to tolerance not achieved, number of iterations
          * <0 : illegal input or breakdown

    Other parameters
    ----------------
    x0 : {array, matrix}
        Starting guess for the solution (a vector of zeros by default).
    tol, atol : float, optional
        Tolerances for convergence, ``norm(residual) <= max(tol*norm(b), atol)``.
        The default for ``atol`` is ``'legacy'``, which emulates
        a different legacy behavior.

        .. warning::

           The default value for `atol` will be changed in a future release.
           For future compatibility, specify `atol` explicitly.
    restart : int, optional
        Number of iterations between restarts. Larger values increase
        iteration cost, but may be necessary for convergence.
        Default is 20.
    maxiter : int, optional
        Maximum number of iterations (restart cycles).  Iteration will stop
        after maxiter steps even if the specified tolerance has not been
        achieved.
    M : {sparse matrix, dense matrix, LinearOperator}
        Inverse of the preconditioner of A.  M should approximate the
        inverse of A and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as `callback(args)`, where `args` are selected by `callback_type`.
    callback_type : {'x', 'pr_norm', 'legacy'}, optional
        Callback function argument requested:
          - ``x``: current iterate (ndarray), called on every restart
          - ``pr_norm``: relative (preconditioned) residual norm (float),
            called on every inner iteration
          - ``legacy`` (default): same as ``pr_norm``, but also changes the
            meaning of 'maxiter' to count inner iterations instead of restart
            cycles.
    orth : 'mgs', optional
        A faster implementation with restarted GMRES based on MGS
    restrt : int, optional
        DEPRECATED - use `restart` instead.

    See Also
    --------
    LinearOperator

    Notes
    -----
    A preconditioner, P, is chosen such that P is close to A but easy to solve
    for. The preconditioner parameter required by this routine is
    ``M = P^-1``. The inverse should preferably not be calculated
    explicitly.  Rather, use the following template to produce M::

      # Construct a linear operator that computes P^-1 * x.
      import scipy.sparse.linalg as spla
      M_x = lambda x: spla.spsolve(P, x)
      M = spla.LinearOperator((n, n), M_x)

    Examples
    --------
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import gmres
    >>> A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> x, exitCode = gmres(A, b)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """

    # Faster restarted GMRES based on modified Gram-Schmidt orthogonalization
    if orth == 'mgs':
        if atol is None:
            x, info = _gmres_mgs(A, b, x0=x0, restart=restart, tol=tol, maxiter=maxiter, M=M, 
                                 callback=callback, restrt=restrt)
            return x, info
        else:
            raise ValueError("The option 'atol' is not used in the case of 'orth == \'mgs\''")

    # Change 'restrt' keyword to 'restart'
    if restrt is None:
        restrt = restart
    elif restart is not None:
        raise ValueError("Cannot specify both restart and restrt keywords. "
                         "Preferably use 'restart' only.")

    if callback is not None and callback_type is None:
        # Warn about 'callback_type' semantic changes.
        # Probably should be removed only in far future, Scipy 2.0 or so.
        warnings.warn("scipy.sparse.linalg.gmres called without specifying `callback_type`. "
                      "The default value will be changed in a future release. "
                      "For compatibility, specify a value for `callback_type` explicitly, e.g., "
                      "``{name}(..., callback_type='pr_norm')``, or to retain the old behavior "
                      "``{name}(..., callback_type='legacy')``",
                      category=DeprecationWarning, stacklevel=3)

    if callback_type is None:
        callback_type = 'legacy'

    if callback_type not in ('x', 'pr_norm', 'legacy'):
        raise ValueError("Unknown callback_type: {!r}".format(callback_type))

    if callback is None:
        callback_type = 'none'

    A, M, x, b,postprocess = make_system(A, M, x0, b)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    if restrt is None:
        restrt = 20
    restrt = min(restrt, n)

    matvec = A.matvec
    psolve = M.matvec
    ltr = _type_conv[x.dtype.char]
    revcom = getattr(_iterative, ltr + 'gmresrevcom')

    bnrm2 = np.linalg.norm(b)
    Mb_nrm2 = np.linalg.norm(psolve(b))
    get_residual = lambda: np.linalg.norm(matvec(x) - b)
    atol = _get_atol(tol, atol, bnrm2, get_residual, 'gmres')
    if atol == 'exit':
        return postprocess(x), 0

    if bnrm2 == 0:
        return postprocess(b), 0

    # Tolerance passed to GMRESREVCOM applies to the inner iteration
    # and deals with the left-preconditioned residual.
    ptol_max_factor = 1.0
    ptol = Mb_nrm2 * min(ptol_max_factor, atol / bnrm2)
    resid = np.nan
    presid = np.nan
    ndx1 = 1
    ndx2 = -1
    # Use _aligned_zeros to work around a f2py bug in Numpy 1.9.1
    work = _aligned_zeros((6+restrt)*n,dtype=x.dtype)
    work2 = _aligned_zeros((restrt+1)*(2*restrt+2),dtype=x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    iter_ = maxiter
    old_ijob = ijob
    first_pass = True
    resid_ready = False
    iter_num = 1
    while True:
        olditer = iter_
        x, iter_, presid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, restrt, work, work2, iter_, presid, info, ndx1, ndx2, ijob, ptol)
        if callback_type == 'x' and iter_ != olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):  # gmres success, update last residual
            if callback_type in ('pr_norm', 'legacy'):
                if resid_ready:
                    callback(presid / bnrm2)
            elif callback_type == 'x':
                callback(x)
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(x)
        elif (ijob == 2):
            work[slice1] = psolve(work[slice2])
            if not first_pass and old_ijob == 3:
                resid_ready = True

            first_pass = False
        elif (ijob == 3):
            work[slice2] *= sclr2
            work[slice2] += sclr1*matvec(work[slice1])
            if resid_ready:
                if callback_type in ('pr_norm', 'legacy'):
                    callback(presid / bnrm2)
                resid_ready = False
                iter_num = iter_num+1

        elif (ijob == 4):
            if ftflag:
                info = -1
                ftflag = False
            resid, info = _stoptest(work[slice1], atol)

            # Inner loop tolerance control
            if info or presid > ptol:
                ptol_max_factor = min(1.0, 1.5 * ptol_max_factor)
            else:
                # Inner loop tolerance OK, but outer loop not.
                ptol_max_factor = max(1e-16, 0.25 * ptol_max_factor)

            if resid != 0:
                ptol = presid * min(ptol_max_factor, atol / resid)
            else:
                ptol = presid * ptol_max_factor

        old_ijob = ijob
        ijob = 2

        if callback_type == 'legacy':
            # Legacy behavior
            if iter_num > maxiter:
                info = maxiter
                break

    if info >= 0 and not (resid <= atol):
        # info isn't set appropriately otherwise
        info = maxiter
        
    return postprocess(x), info


@non_reentrant()
def qmr(A, b, x0=None, tol=1e-5, maxiter=None, M1=None, M2=None, callback=None,
        atol=None):
    """Use Quasi-Minimal Residual iteration to solve ``Ax = b``.

    Parameters
    ----------
    A : {sparse matrix, dense matrix, LinearOperator}
        The real-valued N-by-N matrix of the linear system.
        Alternatively, ``A`` can be a linear operator which can
        produce ``Ax`` and ``A^T x`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : {array, matrix}
        Right hand side of the linear system. Has shape (N,) or (N,1).

    Returns
    -------
    x : {array, matrix}
        The converged solution.
    info : integer
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Other Parameters
    ----------------
    x0  : {array, matrix}
        Starting guess for the solution.
    tol, atol : float, optional
        Tolerances for convergence, ``norm(residual) <= max(tol*norm(b), atol)``.
        The default for ``atol`` is ``'legacy'``, which emulates
        a different legacy behavior.

        .. warning::

           The default value for `atol` will be changed in a future release.
           For future compatibility, specify `atol` explicitly.
    maxiter : integer
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M1 : {sparse matrix, dense matrix, LinearOperator}
        Left preconditioner for A.
    M2 : {sparse matrix, dense matrix, LinearOperator}
        Right preconditioner for A. Used together with the left
        preconditioner M1.  The matrix M1*A*M2 should have better
        conditioned than A alone.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as callback(xk), where xk is the current solution vector.

    See Also
    --------
    LinearOperator

    Examples
    --------
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import qmr
    >>> A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> x, exitCode = qmr(A, b)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A_ = A
    A, M, x, b, postprocess = make_system(A, None, x0, b)

    if M1 is None and M2 is None:
        if hasattr(A_,'psolve'):
            def left_psolve(b):
                return A_.psolve(b,'left')

            def right_psolve(b):
                return A_.psolve(b,'right')

            def left_rpsolve(b):
                return A_.rpsolve(b,'left')

            def right_rpsolve(b):
                return A_.rpsolve(b,'right')
            M1 = LinearOperator(A.shape, matvec=left_psolve, rmatvec=left_rpsolve)
            M2 = LinearOperator(A.shape, matvec=right_psolve, rmatvec=right_rpsolve)
        else:
            def id(b):
                return b
            M1 = LinearOperator(A.shape, matvec=id, rmatvec=id)
            M2 = LinearOperator(A.shape, matvec=id, rmatvec=id)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    ltr = _type_conv[x.dtype.char]
    revcom = getattr(_iterative, ltr + 'qmrrevcom')

    get_residual = lambda: np.linalg.norm(A.matvec(x) - b)
    atol = _get_atol(tol, atol, np.linalg.norm(b), get_residual, 'qmr')
    if atol == 'exit':
        return postprocess(x), 0

    resid = atol
    ndx1 = 1
    ndx2 = -1
    # Use _aligned_zeros to work around a f2py bug in Numpy 1.9.1
    work = _aligned_zeros(11*n,x.dtype)
    ijob = 1
    info = 0
    ftflag = True
    iter_ = maxiter
    while True:
        olditer = iter_
        x, iter_, resid, info, ndx1, ndx2, sclr1, sclr2, ijob = \
           revcom(b, x, work, iter_, resid, info, ndx1, ndx2, ijob)
        if callback is not None and iter_ > olditer:
            callback(x)
        slice1 = slice(ndx1-1, ndx1-1+n)
        slice2 = slice(ndx2-1, ndx2-1+n)
        if (ijob == -1):
            if callback is not None:
                callback(x)
            break
        elif (ijob == 1):
            work[slice2] *= sclr2
            work[slice2] += sclr1*A.matvec(work[slice1])
        elif (ijob == 2):
            work[slice2] *= sclr2
            work[slice2] += sclr1*A.rmatvec(work[slice1])
        elif (ijob == 3):
            work[slice1] = M1.matvec(work[slice2])
        elif (ijob == 4):
            work[slice1] = M2.matvec(work[slice2])
        elif (ijob == 5):
            work[slice1] = M1.rmatvec(work[slice2])
        elif (ijob == 6):
            work[slice1] = M2.rmatvec(work[slice2])
        elif (ijob == 7):
            work[slice2] *= sclr2
            work[slice2] += sclr1*A.matvec(x)
        elif (ijob == 8):
            if ftflag:
                info = -1
                ftflag = False
            resid, info = _stoptest(work[slice1], atol)
        ijob = 2

    if info > 0 and iter_ == maxiter and not (resid <= atol):
        # info isn't set appropriately otherwise
        info = iter_

    return postprocess(x), info


def _gmres_mgs(A, b, x0=None, restart=None, tol=1e-5, maxiter=None, M=None,
               callback=None, restrt=None):
    """
    Use Restarted Generalized Minimal RESidual based on modified Gram-Schmidt
    iteration to solve ``Ax = b``.

    Parameters
    ----------
    A : {sparse matrix, dense matrix, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        `scipy.sparse.linalg.LinearOperator`.
    b : {array, matrix}
        Right hand side of the linear system. Has shape (N,) or (N,1).

    Returns
    -------
    x : {array, matrix}
        The converged solution.
    info : int
        Provides convergence information:
          * 0  : successful exit
          * >0 : convergence to tolerance not achieved, number of iterations
          * <0 : illegal input or breakdown

    Other parameters
    ----------------
    x0 : {array, matrix}
        Starting guess for the solution (a vector of zeros by default).
    restart : int, optional
        Number of iterations between restarts. Larger values increase
        iteration cost, but may be necessary for convergence.
        Default is 30.
    tol : float, optional
        Tolerances for convergence, ``norm(residual) <= tol*norm(b)``.
    maxiter : int, optional
        Maximum number of iterations (restart cycles).  Iteration will stop
        after maxiter steps even if the specified tolerance has not been
        achieved.
    M : {sparse matrix, dense matrix, LinearOperator}
        Inverse of the preconditioner of A.  M should approximate the
        inverse of A and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.
    callback : function
        User-supplied function to call after each iteration.
    restrt : int, optional
        DEPRECATED - use `restart` instead.

    See Also
    --------
    LinearOperator

    Notes
    -----
    A preconditioner, P, is chosen such that P is close to A but easy to solve
    for. The preconditioner parameter required by this routine is
    ``M = P^-1``. The inverse should preferably not be calculated
    explicitly.  Rather, use the following template to produce M::

      # Construct a linear operator that computes P^-1 * x.
      import scipy.sparse.linalg as spla
      M_x = lambda x: spla.spsolve(P, x)
      M = spla.LinearOperator((n, n), M_x)
    """

    # Change 'restrt' keyword to 'restart'
    if restrt is None:
        restrt = restart
    elif restart is not None:
        raise ValueError("Cannot specifiy both restart and restrt keywords. "
                         "Preferably use 'restart' only.")

    # Check data type
    dtype = A.dtype
    if dtype == int:
        dtype = float
        A = A.astype(dtype)
        if b.dtype == int:
            b = b.astype(dtype)

    A, M, x, b, postprocess = make_system(A, M, x0, b)

    # Check if b is a zero vector
    if np.linalg.norm(b) == 0.:
        x = b.copy()
        return (postprocess(x), 0)

    ndofs = A.shape[0]
    if maxiter is None:
        maxiter = min(10000, ndofs*10)
    maxOuterIter = maxiter
    if x0 is None:
        x0 = x.copy()
    if restrt is None:
        restrt = 20
    m = min(restrt, ndofs)

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
    if g[0] == 0.:
        return (postprocess(x), 0)
    r0norm = np.linalg.norm(r)
    if r0norm == 0.:
        return (postprocess(x), 0)

    # Outer iterations
    for iter in range(maxOuterIter):
        V[0] = (1./g[0]) * z  # first basis vector

        # Inner iterations
        for k in range(m):
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
                rres = np.linalg.norm(b - A.matvec(x)) / r0norm
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

            # Compute relative residual
            if k < m - 1:
                rres = np.abs(g[k+1]) / r0norm
                if rres < tol:
                    break

        # Reconstruct the solution
        y = np.zeros((k+1), dtype=dtype)
        gtmp[0:k+1] = g[0:k+1]
        for i in range(k, -1, -1):
            gtmp[i] -= Hn[i, i+1:k+1].dot(y[i+1:k+1])
            y[i] = gtmp[i] / Hn[i, i]
        x += V[0:k+1].T.dot(y[0:k+1])

        # Compute relative residual
        r = b - A.matvec(x)
        z = M.matvec(r)
        g[k+1] = np.linalg.norm(r)
        g[0] = np.linalg.norm(z)
        rres = g[k+1] / r0norm
        if callback is not None:
            callback(x)
        if rres < tol:
            return (postprocess(x), 0)

    return (postprocess(x), maxiter)
