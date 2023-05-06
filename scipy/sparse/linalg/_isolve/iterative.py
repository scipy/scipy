import warnings
from textwrap import dedent
import numpy as np

from scipy.sparse.linalg._interface import LinearOperator
from .utils import make_system
from scipy._lib._threadsafety import non_reentrant
from scipy.linalg import get_blas_funcs, get_lapack_funcs

__all__ = ['bicg', 'bicgstab', 'cg', 'cgs', 'gmres', 'qmr']

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
def gmres(A, b, x0=None, tol=None, restart=None, maxiter=None, M=None,
          callback=None, restrt=None, atol=None, rtol=1e-5,
          callback_type=None):
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

    Other parameters
    ----------------
    x0 : ndarray
        Starting guess for the solution (a vector of zeros by default).
    atol, rtol : float
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(tol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-5``.
    restart : int, optional
        Number of iterations between restarts. Larger values increase
        iteration cost, but may be necessary for convergence.
        If omitted, ``min(20, n)`` is used.
    maxiter : int, optional
        Maximum number of iterations (restart cycles).  Iteration will stop
        after maxiter steps even if the specified tolerance has not been
        achieved. See `callback_type`.
    M : {sparse matrix, ndarray, LinearOperator}
        Inverse of the preconditioner of A.  M should approximate the
        inverse of A and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.
        In this implementation, left preconditioning is used,
        and the preconditioned residual is minimized. However, the final
        convergence is tested with respect to the ``b - A @ x`` residual.
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
            This keyword has no effect if `callback` is not set.
    restrt : int, optional, deprecated

        .. deprecated:: 0.11.0
           `gmres` keyword argument `restrt` is deprecated in favor of
           `restart` and will be removed in SciPy 1.12.0.
    tol : float, optional, deprecated

        .. deprecated 1.11.0
           `gmres` keyword argument `tol` is deprecated in favor of `rtol` and
           will be removed in SciPy 1.13.0

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

    # Handle the deprecation frenzy
    if restrt and restart:
        raise ValueError("Cannot specify both restart and restrt"
                         " keywords. Also 'rstrt' is deprecated."
                         " and will be removed in SciPy 1.12.0. Use "
                         "'restart' instad.")
    if restrt is not None:
        msg = ("'gmres' keyword argument 'restrt' is deprecated "
               "in favor of 'restart' and will be removed in SciPy"
               " 1.12.0. Until then, if set, 'rstrt' will override 'restart'."
               )
        warnings.warn(msg, DeprecationWarning, stacklevel=3)
        restart = restrt

    if tol is not None:
        msg = ("'gmres' keyword argument 'tol' is deprecated "
               "in favor of 'rtol' and will be removed in SciPy "
               " and will be removed in SciPy v.1.13.0.")
        warnings.warn(msg, category=DeprecationWarning, stacklevel=3)
        rtol = float(tol)

    if atol == 'legacy':
        msg = ("Use of strings for 'atol' keyword is deprecated "
               "and will be removed in SciPy 1.13. To keep legacy "
               "behavior set 'atol=0., rtol=1e-5'. See "
               " ``scipy.sparse.linalg.gmres`` documentation for "
               "more details."
               )
        warnings.warn(msg, category=DeprecationWarning, stacklevel=3)

    if callback is not None and callback_type is None:
        # Warn about 'callback_type' semantic changes.
        # Probably should be removed only in far future, Scipy 2.0 or so.
        msg = ("scipy.sparse.linalg.gmres called without specifying "
               "`callback_type`. The default value will be changed in"
               " a future release. For compatibility, specify a value "
               "for `callback_type` explicitly, e.g., "
               "``gmres(..., callback_type='pr_norm')``, or to retain the "
               "old behavior ``gmres(..., callback_type='legacy')``"
               )
        warnings.warn(msg, category=DeprecationWarning, stacklevel=3)

    if callback_type is None:
        callback_type = 'legacy'

    if callback_type not in ('x', 'pr_norm', 'legacy'):
        raise ValueError(f"Unknown callback_type: {callback_type!r}")

    if callback is None:
        callback_type = None

    A, M, x, b, postprocess = make_system(A, M, x0, b)

    if np.iscomplexobj(x):
        dotprod = np.vdot
    else:
        dotprod = np.dot

    eps = np.finfo(x.dtype.char).eps

    n = len(b)

    if maxiter is None:
        maxiter = n*10

    if restart is None:
        restart = 20
    restart = min(restart, n)

    matvec = A.matvec
    psolve = M.matvec

    bnrm2 = np.linalg.norm(b)
    if bnrm2 == 0:
        return postprocess(b), 0

    atol = max(float(atol), float(rtol)*bnrm2)

    lartg = get_lapack_funcs('lartg', dtype=x.dtype)
    trsv = get_blas_funcs('trsv', dtype=x.dtype)

    # ====================================================
    # =========== Tolerance control from gh-8400 =========
    # ====================================================
    Mb_nrm2 = np.linalg.norm(psolve(b))
    ptol_max_factor = 1.0
    ptol = Mb_nrm2 * min(ptol_max_factor, atol / bnrm2)
    presid = np.nan
    # ====================================================

    # allocate internal variables
    v = np.empty([restart+1, n], dtype=x.dtype)
    h = np.zeros([restart, restart+1], dtype=x.dtype)
    givens = np.zeros([restart, 2], dtype=x.dtype)
    r = b - matvec(x) if x0 is None else b.copy()

    # legacy iteration count
    inner_iters = 0

    for iteration in range(maxiter):
        rnorm = np.linalg.norm(r)
        if rnorm < atol:  # Are we done?
            return postprocess(x), 0

        # Tightening of the inner loop tol, see gh-8400 for the reasoning.

        # The case "true residual passes tol-test but preconditioned
        # residual does not" is not considered. Tolerance problem is too
        # sophisticated for manipulating two crude knobs hence will at some
        # point require a more structured handling.
        if iteration > 0:
            if presid <= ptol:
                # Inner loop passed but we are still here
                # hence tighten the inner tolerance
                ptol_max_factor = max(1e-16, 0.25 * ptol_max_factor)
                ptol = presid * min(ptol_max_factor, atol / rnorm)

        v[0, :] = psolve(r)
        v[0, :] *= (1 / rnorm)
        # RHS of the Hessenberg problem
        S = np.zeros(restart+1, dtype=x.dtype)
        S[0] = rnorm

        # Arnoldi process
        breakdown = False
        for col in range(restart):
            av = matvec(v[col, :])
            w = psolve(av)

            # Modified Gram-Schmidt
            h0 = np.linalg.norm(w)
            for k in range(col+1):
                tmp = dotprod(v[k, :], w)
                h[col, k] = tmp
                w -= tmp*v[k, :]

            h1 = np.linalg.norm(w)
            h[col, col + 1] = h1
            v[col + 1, :] = w[:]

            # Exact solution indicator
            if h1 <= eps*h0:
                h[col, col + 1] = 0
                breakdown = True
            else:
                v[col + 1, :] *= (1 / h1)

            # apply past Givens rotations to current h column
            for k in range(col):
                c, s = givens[k, 0], givens[k, 1]
                n0, n1 = h[col, [k, k+1]]
                h[col, [k, k + 1]] = [c*n0 + s*n1, -s.conj()*n0 + c*n1]

            # get and apply current rotation to h and S
            c, s, mag = lartg(*h[col, [col, col+1]])
            givens[col, :] = [c, s]
            h[col, [col, col+1]] = mag, 0

            # S[col+1] component is always 0
            tmp = -np.conjugate(s)*S[col]
            S[[col, col + 1]] = [c*S[col], tmp]
            presid = np.abs(tmp)

            # Legacy callback behavior per outer
            if callback_type in ('pr_norm', 'legacy'):
                callback(presid / bnrm2)

            inner_iters += 1
            if callback_type == 'legacy':
                # Legacy behavior
                if inner_iters >= maxiter:
                    break

            if presid <= ptol or breakdown:
                break

        # Outer loop continues

        # Solve (col, col) upper triangular system
        # allow trsv to pseudo-solve singular cases
        if breakdown:
            S[col] = 0

        y = trsv(h[:col+1, :col+1].T, S[:col+1])
        g = y @ v[:col+1, :]
        x += g
        r = b - matvec(x)

        if callback is not None:
            if callback_type in ('pr_norm', 'legacy'):
                callback(presid / bnrm2)

                # Exit if ran out of inner iterations
                if callback_type == 'legacy':
                    if inner_iters >= maxiter:
                        # Return incomplete progress
                        return postprocess(x), maxiter
            else:
                callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return postprocess(x), maxiter


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
