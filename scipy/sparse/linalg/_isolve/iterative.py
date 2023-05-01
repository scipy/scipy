"""Iterative methods for solving linear systems"""

__all__ = ['bicg', 'bicgstab', 'cg', 'cgs', 'gmres', 'qmr']

import warnings
from textwrap import dedent
import numpy as np

from . import _iterative

from scipy.sparse.linalg._interface import LinearOperator
from .utils import make_system
from scipy._lib._util import _aligned_zeros
from scipy._lib._threadsafety import non_reentrant

_type_conv = {'f': 's', 'd': 'd', 'F': 'c', 'D': 'z'}


# Part of the docstring common to all iterative solvers
common_doc1 = \
"""
Parameters
----------
A : {sparse matrix, ndarray, LinearOperator}"""

common_doc2 = \
"""b : ndarray
    Right hand side of the linear system. Has shape (N,) or (N,1).

Returns
-------
x : ndarray
    The converged solution.
info : integer
    Provides convergence information:
        0  : successful exit
        >0 : convergence to tolerance not achieved, number of iterations
        <0 : illegal input or breakdown

Other Parameters
----------------
x0 : ndarray
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
M : {sparse matrix, ndarray, LinearOperator}
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
                                common_doc2, dedent(footer)))
        return fn
    return combine


@set_docstring('Use BIConjugate Gradient iteration to solve ``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` and ``A^T x`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.',
               footer="""\
               Examples
               --------
               >>> import numpy as np
               >>> from scipy.sparse import csc_matrix
               >>> from scipy.sparse.linalg import bicg
               >>> A = csc_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1.]])
               >>> b = np.array([2., 4., -1.])
               >>> x, exitCode = bicg(A, b)
               >>> print(exitCode)  # 0 indicates successful convergence
               0
               >>> np.allclose(A.dot(x), b)
               True

               """
               )
def bicg(A, b, x0=None, tol=1e-5, maxiter=None, M=None,
         callback=None, atol=0.):
    A, M, x, b, postprocess = make_system(A, M, x0, b)
    n = len(b)
    if (np.iscomplexobj(A) or np.iscomplexobj(b)):
        dotprod = np.vdot
    else:
        dotprod = np.dot

    if maxiter is None:
        maxiter = n*10

    matvec, rmatvec = A.matvec, A.rmatvec
    psolve, rpsolve = M.matvec, M.rmatvec

    if A.dtype.char in 'fF':
        rhotol = np.finfo(np.float32).eps ** 2
    else:
        rhotol = np.spacing(1.) ** 2

    # Deprecate the legacy usage
    if isinstance(atol, str):
        warnings.warn("scipy.sparse.linalg.cg called with `atol` set to "
                      "string, possibly with value 'legacy'. This behavior "
                      "is deprecated and atol parameter only excepts floats."
                      " In SciPy 1.13, this will result with an error.",
                      category=DeprecationWarning, stacklevel=3)
        tol = float(tol)

        if np.linalg.norm(b) == 0:
            atol = tol
        else:
            atol = tol * float(np.linalg.norm(b))
    else:
        atol = max(float(atol), tol * float(np.linalg.norm(b)))

    # Is there any tolerance set? since b can be all 0.
    if atol == 0.:
        atol = 10*n*np.finfo(A.dtype).eps

    # Dummy values to initialize vars, silence linter warnings
    rho_prev, p, ptilde = None, None, None

    r = b - matvec(x) if x0 is not None else b.copy()
    rtilde = r.copy()

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return postprocess(x), 0

        z = psolve(r)
        ztilde = rpsolve(rtilde)
        rho_cur = dotprod(z, rtilde)

        if np.abs(rho_cur) < rhotol:  # Breakdown case
            # It brokedown but maybe converged?
            if np.linalg.norm(r) < atol:
                return postprocess(x), 0
            else:
                return postprocess, -10

        if iteration > 0:
            beta = rho_cur / rho_prev
            p *= beta
            p += z
            ptilde *= beta.conj()
            ptilde += ztilde
        else:  # First spin
            p = np.empty_like(r)
            ptilde = np.empty_like(r)
            p[:] = z[:]
            ptilde[:] = ztilde[:]

        q = matvec(p)
        qtilde = rmatvec(ptilde)
        rv = dotprod(ptilde, q)

        if rv == 0:
            return postprocess(x), -11

        alpha = rho_cur / rv
        x += alpha*p
        r -= alpha*q
        rtilde -= alpha.conj()*qtilde
        rho_prev = rho_cur

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return postprocess(x), maxiter


@set_docstring('Use BIConjugate Gradient STABilized iteration to solve '
               '``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.',
               footer="""\
               Examples
               --------
               >>> import numpy as np
               >>> from scipy.sparse import csc_matrix
               >>> from scipy.sparse.linalg import bicgstab
               >>> R = np.array([[4, 2, 0, 1],
               ...               [3, 0, 0, 2],
               ...               [0, 1, 1, 1],
               ...               [0, 2, 1, 0]])
               >>> A = csc_matrix(R)
               >>> b = np.array([-1, -0.5, -1, 2])
               >>> x, exit_code = bicgstab(A, b)
               >>> print(exit_code)  # 0 indicates successful convergence
               0
               >>> np.allclose(A.dot(x), b)
               True
               """)
def bicgstab(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None,
             atol=0.):
    A, M, x, b, postprocess = make_system(A, M, x0, b)
    n = len(b)

    if (np.iscomplexobj(A) or np.iscomplexobj(b)):
        dotprod = np.vdot
    else:
        dotprod = np.dot

    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec

    # These values make no sense but coming from original Fortran code
    # sqrt might have been meant instead.
    if A.dtype.char in 'fF':
        rhotol = np.finfo(np.float32).eps ** 2
    else:
        rhotol = np.spacing(1.) ** 2
    omegatol = rhotol

    # Deprecate the legacy usage
    if isinstance(atol, str):
        warnings.warn("scipy.sparse.linalg.cg called with `atol` set to "
                      "string, possibly with value 'legacy'. This behavior "
                      "is deprecated and atol parameter only excepts floats."
                      " In SciPy 1.13, this will result with an error.",
                      category=DeprecationWarning, stacklevel=3)
        tol = float(tol)

        if np.linalg.norm(b) == 0:
            atol = tol
        else:
            atol = tol * float(np.linalg.norm(b))
    else:
        atol = max(float(atol), tol * float(np.linalg.norm(b)))

    # Is there any tolerance set? since b can be all 0.
    if atol == 0.:
        atol = 10*n*np.finfo(A.dtype).eps

    # Dummy values to initialize vars, silence linter warnings
    rho_prev, omega, alpha, p, v = None, None, None, None, None

    r = b - matvec(x) if x0 is not None else b.copy()
    rtilde = r.copy()

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return postprocess(x), 0

        rho = dotprod(rtilde, r)
        if np.abs(rho) < rhotol:  # rho breakdown
            return postprocess, -10

        if iteration > 0:
            if np.abs(omega) < omegatol:  # omega breakdown
                return postprocess(x), -11

            beta = (rho / rho_prev) * (alpha / omega)
            p -= omega*v
            p *= beta
            p += r
        else:  # First spin
            s = np.empty_like(r)
            p = r.copy()

        phat = psolve(p)
        v = matvec(phat)
        rv = dotprod(rtilde, v)
        if rv == 0:
            return postprocess(x), -11
        alpha = rho / rv
        r -= alpha*v
        s[:] = r[:]

        if np.linalg.norm(s) < atol:
            x += alpha*phat
            return postprocess(x), 0

        shat = psolve(s)
        t = matvec(shat)
        omega = dotprod(t, s) / dotprod(t, t)
        x += alpha*phat
        x += omega*shat
        r -= omega*t
        rho_prev = rho

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return postprocess(x), maxiter


@set_docstring('Use Conjugate Gradient iteration to solve ``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               '``A`` must represent a hermitian, positive definite matrix.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.',
               footer="""\
               Examples
               --------
               >>> import numpy as np
               >>> from scipy.sparse import csc_matrix
               >>> from scipy.sparse.linalg import cg
               >>> P = np.array([[4, 0, 1, 0],
               ...               [0, 5, 0, 0],
               ...               [1, 0, 3, 2],
               ...               [0, 0, 2, 4]])
               >>> A = csc_matrix(P)
               >>> b = np.array([-1, -0.5, -1, 2])
               >>> x, exit_code = cg(A, b)
               >>> print(exit_code)    # 0 indicates successful convergence
               0
               >>> np.allclose(A.dot(x), b)
               True

               """)
def cg(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None, atol=0.):
    A, M, x, b, postprocess = make_system(A, M, x0, b)
    n = len(b)

    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec
    r = b - matvec(x) if x0 is not None else b.copy()

    # Deprecate the legacy usage
    if isinstance(atol, str):
        warnings.warn("scipy.sparse.linalg.cg called with `atol` set to "
                      "string, possibly with value 'legacy'. This behavior "
                      "is deprecated and atol parameter only excepts floats."
                      " In SciPy 1.13, this will result with an error.",
                      category=DeprecationWarning, stacklevel=3)
        tol = float(tol)

        if np.linalg.norm(b) == 0:
            atol = tol
        else:
            atol = tol * float(np.linalg.norm(b))
    else:
        atol = max(float(atol), tol * float(np.linalg.norm(b)))

    # Is there any tolerance set? since b can be all 0.
    if atol == 0.:
        atol = 10*n*np.finfo(A.dtype).eps

    # Dummy value to initialize var, silences warnings
    rho_prev, p = None, None
    if (np.iscomplexobj(A) or np.iscomplexobj(b)):
        dotprod = np.vdot
    else:
        dotprod = np.dot

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return postprocess(x), 0

        z = psolve(r)
        rho_cur = dotprod(r, z)
        if iteration > 0:
            beta = rho_cur / rho_prev
            p *= beta
            p += z
        else:  # First spin
            p = np.empty_like(r)
            p[:] = z[:]

        q = matvec(p)
        alpha = rho_cur / dotprod(p, q)
        x += alpha*p
        r -= alpha*q
        rho_prev = rho_cur

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return postprocess(x), maxiter


@set_docstring('Use Conjugate Gradient Squared iteration to solve ``Ax = b``.',
               'The real-valued N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.',
               footer="""\
               Examples
               --------
               >>> import numpy as np
               >>> from scipy.sparse import csc_matrix
               >>> from scipy.sparse.linalg import cgs
               >>> R = np.array([[4, 2, 0, 1],
               ...               [3, 0, 0, 2],
               ...               [0, 1, 1, 1],
               ...               [0, 2, 1, 0]])
               >>> A = csc_matrix(R)
               >>> b = np.array([-1, -0.5, -1, 2])
               >>> x, exit_code = cgs(A, b)
               >>> print(exit_code)  # 0 indicates successful convergence
               0
               >>> np.allclose(A.dot(x), b)
               True
               """
               )
def cgs(A, b, x0=None, tol=1e-5, maxiter=None, M=None, callback=None, atol=0.):
    A, M, x, b, postprocess = make_system(A, M, x0, b)
    n = len(b)

    if (np.iscomplexobj(A) or np.iscomplexobj(b)):
        dotprod = np.vdot
    else:
        dotprod = np.dot

    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec

    if A.dtype.char in 'fF':
        rhotol = np.finfo(np.float32).eps ** 2
    else:
        rhotol = np.spacing(1.) ** 2

    r = b - matvec(x) if x0 is not None else b.copy()

    # Deprecate the legacy usage
    if isinstance(atol, str):
        warnings.warn("scipy.sparse.linalg.cg called with `atol` set to "
                      "string, possibly with value 'legacy'. This behavior "
                      "is deprecated and atol parameter only excepts floats."
                      " In SciPy 1.13, this will result with an error.",
                      category=DeprecationWarning, stacklevel=3)
        tol = float(tol)

        if np.linalg.norm(b) == 0:
            atol = tol
        else:
            atol = tol * float(np.linalg.norm(b))
    else:
        atol = max(float(atol), tol * float(np.linalg.norm(b)))

    # Is there any tolerance set? since b can be all 0.
    if atol == 0.:
        atol = 10*n*np.finfo(A.dtype).eps

    rtilde = r.copy()
    bnorm = np.linalg.norm(b)
    if bnorm == 0:
        bnorm = 1

    # Dummy values to initialize vars, silence linter warnings
    rho_prev, p, u, q = None, None, None, None

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return postprocess(x), 0

        rho_cur = dotprod(rtilde, r)
        if np.abs(rho_cur) < rhotol:  # Breakdown case
            # It brokedown but maybe converged?
            if np.linalg.norm(r) < atol:
                return postprocess(x), 0
            else:
                return postprocess, -10

        if iteration > 0:
            beta = rho_cur / rho_prev

            # u = r + beta * q
            # p = u + beta * (q + beta * p);
            u[:] = r[:]
            u += beta*q

            p *= beta
            p += q
            p *= beta
            p += u

        else:  # First spin
            p = np.empty_like(r)
            u = np.empty_like(r)
            q = np.empty_like(r)
            p[:] = r[:]
            u[:] = r[:]

        phat = psolve(p)
        vhat = matvec(phat)
        rv = dotprod(rtilde, vhat)

        if rv == 0:  # Dot product breakdown
            return postprocess(x), -11

        alpha = rho_cur / rv
        q[:] = u[:]
        q -= alpha*vhat
        uhat = psolve(u + q)
        x += alpha*uhat
        qhat = matvec(uhat)
        r -= alpha*qhat
        rho_prev = rho_cur

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return postprocess(x), maxiter


@non_reentrant()
def gmres(A, b, x0=None, tol=1e-5, restart=None, maxiter=None, M=None, callback=None,
          restrt=None, atol=None, callback_type=None):
    """
    Use Generalized Minimal RESidual iteration to solve ``Ax = b``.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, ``A`` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:
          * 0  : successful exit
          * >0 : convergence to tolerance not achieved, number of iterations
          * <0 : illegal input or breakdown

    Other parameters
    ----------------
    x0 : ndarray
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
    M : {sparse matrix, ndarray, LinearOperator}
        Inverse of the preconditioner of A.  M should approximate the
        inverse of A and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.
        In this implementation, left preconditioning is used,
        and the preconditioned residual is minimized.
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
    restrt : int, optional, deprecated

        .. deprecated:: 0.11.0
           `gmres` keyword argument `restrt` is deprecated infavour of
           `restart` and will be removed in SciPy 1.12.0.

    See Also
    --------
    LinearOperator

    Notes
    -----
    A preconditioner, P, is chosen such that P is close to A but easy to solve
    for. The preconditioner parameter required by this routine is
    ``M = P^-1``. The inverse should preferably not be calculated
    explicitly.  Rather, use the following template to produce M::

      # Construct a linear operator that computes P^-1 @ x.
      import scipy.sparse.linalg as spla
      M_x = lambda x: spla.spsolve(P, x)
      M = spla.LinearOperator((n, n), M_x)

    Examples
    --------
    >>> import numpy as np
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

    # Change 'restrt' keyword to 'restart'
    if restrt is None:
        restrt = restart
    elif restart is not None:
        raise ValueError("Cannot specify both restart and restrt keywords. "
                         "Preferably use 'restart' only.")
    else:
        msg = ("'gmres' keyword argument 'restrt' is deprecated infavour of "
               "'restart' and will be removed in SciPy 1.12.0.")
        warnings.warn(msg, DeprecationWarning, stacklevel=2)

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
        raise ValueError(f"Unknown callback_type: {callback_type!r}")

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
    def get_residual():
        return np.linalg.norm(matvec(x) - b)
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


def qmr(A, b, x0=None, tol=1e-5, maxiter=None, M1=None, M2=None, callback=None,
        atol=0.):
    """Use Quasi-Minimal Residual iteration to solve ``Ax = b``.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator}
        The real-valued N-by-N matrix of the linear system.
        Alternatively, ``A`` can be a linear operator which can
        produce ``Ax`` and ``A^T x`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).

    Returns
    -------
    x : ndarray
        The converged solution.
    info : integer
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : illegal input or breakdown

    Other Parameters
    ----------------
    x0 : ndarray
        Starting guess for the solution.
    tol, atol : float, optional
        Tolerances for convergence,
        ``norm(residual) <= max(tol*norm(b), atol)``. The default for ``atol``
        is ``'legacy'``, which emulates a different legacy behavior.

        .. warning::

           The default value for `atol` will be changed in a future release.
           For future compatibility, specify `atol` explicitly.
    maxiter : integer
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M1 : {sparse matrix, ndarray, LinearOperator}
        Left preconditioner for A.
    M2 : {sparse matrix, ndarray, LinearOperator}
        Right preconditioner for A. Used together with the left
        preconditioner M1.  The matrix M1@A@M2 should have better
        conditioned than A alone.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as callback(xk), where xk is the current solution vector.

    See Also
    --------
    LinearOperator

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import qmr
    >>> A = csc_matrix([[3., 2., 0.], [1., -1., 0.], [0., 5., 1.]])
    >>> b = np.array([2., 4., -1.])
    >>> x, exitCode = qmr(A, b)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A_ = A
    A, M, x, b, postprocess = make_system(A, None, x0, b)

    if M1 is None and M2 is None:
        if hasattr(A_, 'psolve'):
            def left_psolve(b):
                return A_.psolve(b, 'left')

            def right_psolve(b):
                return A_.psolve(b, 'right')

            def left_rpsolve(b):
                return A_.rpsolve(b, 'left')

            def right_rpsolve(b):
                return A_.rpsolve(b, 'right')
            M1 = LinearOperator(A.shape,
                                matvec=left_psolve,
                                rmatvec=left_rpsolve)
            M2 = LinearOperator(A.shape,
                                matvec=right_psolve,
                                rmatvec=right_rpsolve)
        else:
            def id(b):
                return b
            M1 = LinearOperator(A.shape, matvec=id, rmatvec=id)
            M2 = LinearOperator(A.shape, matvec=id, rmatvec=id)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    if (np.iscomplexobj(A) or np.iscomplexobj(b)):
        dotprod = np.vdot
    else:
        dotprod = np.dot

    rhotol = np.finfo(A.dtype.char).eps
    betatol = rhotol
    gammatol = rhotol
    deltatol = rhotol
    epsilontol = rhotol
    xitol = rhotol

    # Deprecate the legacy usage
    if isinstance(atol, str):
        warnings.warn("scipy.sparse.linalg.cg called with `atol` set to "
                      "string, possibly with value 'legacy'. This behavior "
                      "is deprecated and atol parameter only excepts floats."
                      " In SciPy 1.13, this will result with an error.",
                      category=DeprecationWarning, stacklevel=3)
        tol = float(tol)

        if np.linalg.norm(b) == 0:
            atol = tol
        else:
            atol = tol * float(np.linalg.norm(b))
    else:
        atol = max(float(atol), tol * float(np.linalg.norm(b)))

    r = b - A.matvec(x) if x0 is not None else b.copy()
    vtilde = r.copy()
    y = M1.matvec(vtilde)
    rho = np.linalg.norm(y)
    wtilde = r.copy()
    z = M2.rmatvec(wtilde)
    xi = np.linalg.norm(z)
    gamma, eta, theta = 1, -1, 0
    v = np.empty_like(vtilde)
    w = np.empty_like(wtilde)

    # Dummy values to initialize vars, silence linter warnings
    epsilon, q, d, p, s = None, None, None, None, None

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return postprocess(x), 0
        if np.abs(rho) < rhotol:  # rho breakdown
            return postprocess(x), -10
        if np.abs(xi) < xitol:  # xi breakdown
            return postprocess(x), -15

        v[:] = vtilde[:]
        v *= (1 / rho)
        y *= (1 / rho)
        w[:] = wtilde[:]
        w *= (1 / xi)
        z *= (1 / xi)
        delta = dotprod(z, y)

        if np.abs(delta) < deltatol:  # delta breakdown
            return postprocess(x), -13

        ytilde = M2.matvec(y)
        ztilde = M1.rmatvec(z)

        if iteration > 0:
            ytilde -= (xi * delta / epsilon) * p
            p[:] = ytilde[:]
            ztilde -= (rho * (delta / epsilon).conj()) * q
            q[:] = ztilde[:]
        else:  # First spin
            p = ytilde.copy()
            q = ztilde.copy()

        ptilde = A.matvec(p)
        epsilon = dotprod(q, ptilde)
        if np.abs(epsilon) < epsilontol:  # epsilon breakdown
            return postprocess(x), -14

        beta = epsilon / delta
        if np.abs(beta) < betatol:  # beta breakdown
            return postprocess(x), -11

        vtilde[:] = ptilde[:]
        vtilde -= beta*v
        y = M1.matvec(vtilde)

        rho_prev = rho
        rho = np.linalg.norm(y)
        wtilde[:] = w[:]
        wtilde *= - beta.conj()
        wtilde += A.rmatvec(q)
        z = M2.rmatvec(wtilde)
        xi = np.linalg.norm(z)
        gamma_prev = gamma
        theta_prev = theta
        theta = rho / (gamma_prev * np.abs(beta))
        gamma = 1 / np.sqrt(1 + theta**2)

        if np.abs(gamma) < gammatol:  # gamma breakdown
            return postprocess(x), -12

        eta *= -(rho_prev / beta) * (gamma / gamma_prev)**2

        if iteration > 0:
            d *= (theta_prev * gamma) ** 2
            d += eta*p
            s *= (theta_prev * gamma) ** 2
            s += eta*ptilde
        else:
            d = p.copy()
            d *= eta
            s = ptilde.copy()
            s *= eta

        x += d
        r -= s

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return postprocess(x), maxiter
