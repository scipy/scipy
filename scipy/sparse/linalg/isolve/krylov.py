"""Iterative methods for solving linear systems"""

__all__ = ['cg','cgs','icgs','bicg','bicgstab','cgnr','cgne','gmres','qmr']

import warnings
import numpy as np
import time

from scipy.sparse.linalg.isolve import _iterative

from scipy.sparse.linalg.interface import LinearOperator
#from .utils import make_system
from scipy.sparse.linalg.isolve.utils import make_system
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


def set_docstring(header, Ainfo, footer='', atol_default='0'):
    def combine(fn):
        fn.__doc__ = '\n'.join((header, common_doc1,
                                '    ' + Ainfo.replace('\n', '\n    '),
                                common_doc2, footer))
        return fn
    return combine


#======================
#   CG
#======================
@set_docstring('Use Conjugate Gradient iteration to solve ``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               '``A`` must represent a hermitian, positive definite matrix.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def cg(A, b, x0=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	r = b - A.matvec(x)  # A is a LinearOperator
	z = M.matvec(r)  # z = M * r
	p = z.copy()
	r0norm = np.linalg.norm(r)
	if r0norm == 0.:
		return (postprocess(x), 0)
	for iter in range(maxit):
		rho0 = np.inner(r.conjugate(), z)
		q = A.matvec(p)
		alpha = rho0 / np.inner(p.conjugate(), q)
		x += alpha * p
		r -= alpha * q
		alpha = np.linalg.norm(r) / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format(iter+1, alpha))
		if alpha < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format(iter+1))
			print(" --- system solved with SciPy/CG (in {} s)".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format(iter+1))
			return (postprocess(x), 1)
		z = M.matvec(r)
		rho1 = np.inner(r.conjugate(), z)
		alpha = rho1 / rho0
		p = z + alpha * p


#======================
#   Traditional CGS
#======================
@set_docstring('Use Conjugate Gradient Squared iteration to solve ``Ax = b``.',
               'The real-valued N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def cgs(A, b, x0=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	r = b - A.matvec(x)
	rhat = r.copy()
	d_old = np.inner(rhat.conjugate(), r)
	p = r.copy()
	q = r.copy()
	beta = 0.
	r0norm = np.linalg.norm(r)
	if r0norm == 0.:
		return (postprocess(x), 0)
	for iter in range(maxit):
		w = r + beta * q
		p = w + beta * q + (beta**2) * p
		z = M.matvec(p)
		Az = A.matvec(z)
		alpha = d_old / np.inner(rhat.conjugate(), Az)
		q = w - alpha * Az
		w += q
		z = alpha * M.matvec(w)
		x += z
		r -= A.matvec(z)
		alpha = np.linalg.norm(r) / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format(iter+1, alpha))
		if alpha < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format(iter+1))
			print(" --- system solved with SciPy/CGS (in {} s)".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format(iter+1))
			return (postprocess(x), 1)
		d = np.inner(rhat.conjugate(), r)
		beta = d / d_old
		d_old = d


#======================
#   Improved CGS
#======================
@set_docstring('Use Improved Conjugate Gradient Squared iteration to solve ``Ax = b``.',
               'The real-valued N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def icgs(A, b, x0=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	r = b - A.matvec(x)
	# why not set rhat = M.matvec(r) directly ?
	tmp = r.copy()
	rhat = M.matvec(tmp)
	z = rhat.copy()
	d_old = np.inner(rhat.conjugate(), z)
	p = r.copy()
	q = r.copy()
	beta = 0.
	r0norm = np.linalg.norm(r)
	if r0norm == 0.:
		return (postprocess(x), 0)
	for iter in range(maxit):
		w = z + beta * q
		p = w + beta * q + (beta**2) * p
		Ap = A.matvec(p)
		z = M.matvec(Ap)
		alpha = d_old / np.inner(rhat.conjugate(), z)
		q = w - alpha * z
		w += q
		x += alpha * w
		r -= alpha * A.matvec(w)
		alpha = np.linalg.norm(r) / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format(iter+1, alpha))
		if alpha < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format(iter+1))
			print(" --- system solved with SciPy/Improved CGS (in {} s)".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format(iter+1))
			return (postprocess(x), 1)
		z = M.matvec(r)
		d = np.inner(rhat.conjugate(), z)
		beta = d / d_old
		d_old = d


#======================
#   BiCG
#======================
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
def bicg(A, b, x0=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	At = A.H
	Mt = M.H
	r = b - A.matvec(x)
	z = M.matvec(r)
	rhat = z.copy()
	z_hat = Mt.matvec(rhat)
	d_old = np.inner(rhat.conjugate(), z)
	p = r.copy()
	p_hat = z.copy()
	beta = 0.
	r0norm = np.linalg.norm(r)
	if r0norm == 0.:
		return (postprocess(x), 0)
	for iter in range(maxit):
		p = z + beta * p
		p_hat = z_hat + beta * p_hat
		Ap = A.matvec(p)
		alpha = d_old / np.inner(p_hat.conjugate(), Ap)
		x += alpha * p
		r -= alpha * Ap
		beta = np.linalg.norm(r) / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format(iter+1, beta))
		if beta < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format(iter+1))
			print(" --- system solved with SciPy/BiCG (in {} s)".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format(iter+1))
			return (postprocess(x), 1)
		Ap = At.matvec(p_hat)
		rhat -= alpha * Ap  # alpha * At * p_hat
		z = M.matvec(r)
		d = np.inner(rhat.conjugate(), z)
		beta = d / d_old
		d_old = d.copy()
		z_hat = Mt.matvec(rhat)


#======================
#   BiCGSTAB
#======================
@set_docstring('Use BIConjugate Gradient STABilized iteration to solve '
               '``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def bicgstab(A, b, x0=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	r = b - A.matvec(x)
	rhat = r.copy()
	r0norm = np.linalg.norm(r)
	if r0norm == 0.:
		return (postprocess(x), 0)
	for iter in range(maxit):
		rho1 = np.inner(rhat.conjugate(), r)
		if rho1 == 0.:
			print("Iterative method failed due to inner product of residuals is zero")
		if iter == 0:
			p = r.copy()
		else:
			beta = (rho1/rho0) * (alpha/omega)
			p = r + beta * (p - omega*v)
		p_hat = M.matvec(p)
		v = A.matvec(p_hat)
		alpha = rho1 / np.inner(rhat.conjugate(), v)
		s = r - alpha * v
		if np.linalg.norm(s) < 1e-14:
			x += alpha * p_hat
			return (postprocess(x), 0)
		s_hat = M.matvec(s)
		t = A.matvec(s_hat)
		omega = np.inner(t.conjugate(), s) / np.inner(t.conjugate(), t)
		x += alpha * p_hat + omega * s_hat
		r = s - omega * t
		beta = np.linalg.norm(r) / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format(iter+1, beta))
		if beta < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format(iter+1))
			print(" --- system solved with SciPy/BiCGSTAB (in {} s)".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format(iter+1))
			return (postprocess(x), 1)
		rho0 = rho1;


#======================
#   CGNR
#======================
@set_docstring('Use Conjugate Gradient Normal Residual iteration to solve '
               '``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def cgnr(A, b, x0=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	r = b - A.matvec(x)
	At = A.H
	rhat = At.matvec(r)
	z = M.matvec(rhat)
	p = z.copy()
	ztrhat_old = np.inner(z.conjugate(), rhat)
	r0norm = np.linalg.norm(r)
	if r0norm == 0.:
		return (postprocess(x), 0)
	for iter in range(maxit):
		w = A.matvec(p)
		alpha = ztrhat_old / np.inner(w.conjugate(), w)
		x += alpha * p
		r -= alpha * w
		alpha = np.linalg.norm(r) / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format(iter+1, alpha))
		if alpha < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format(iter+1))
			print(" --- system solved with SciPy/CGNR (in {} s)".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format(iter+1))
			return (postprocess(x), 1)
		rhat = At.matvec(r)
		z = M.matvec(rhat)
		ztrhat = np.inner(z.conjugate(), rhat)
		beta = ztrhat / ztrhat_old 
		p = z + beta * p
		ztrhat_old = ztrhat


#======================
#   CGNE
#======================
@set_docstring('Use Conjugate Gradient Normal Error iteration to solve '
               '``Ax = b``.',
               'The real or complex N-by-N matrix of the linear system.\n'
               'Alternatively, ``A`` can be a linear operator which can\n'
               'produce ``Ax`` using, e.g.,\n'
               '``scipy.sparse.linalg.LinearOperator``.')
@non_reentrant()
def cgne(A, b, x0=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	r = b - A.matvec(x)
	At = A.H
	z = M.matvec(r)
	p = At.matvec(z)
	ztr_old = np.inner(z.conjugate(), r)
	r0norm = np.linalg.norm(r)
	if r0norm == 0.:
		return (postprocess(x), 0)
	for iter in range(maxit):
		w = A.matvec(p)
		alpha = ztr_old / np.inner(p.conjugate(), p)
		x += alpha * p
		r -= alpha * w
		alpha = np.linalg.norm(r) / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format(iter+1, alpha))
		if alpha < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format(iter+1))
			print(" --- system solved with SciPy/CGNR (in {})".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format(iter+1))
			return (postprocess(x), 1)
		z = M.matvec(r)
		ztr = np.inner(z.conjugate(), r)
		beta = ztr / ztr_old
		p = At.matvec(z) + beta * p
		ztr_old = ztr


#======================
#   Restarted GMRES
#======================
@non_reentrant()
def gmres(A, b, x0=None, restart=None, tol=1e-5, maxit=10000, dtype=None, sparams=None, M=None):
	"""
    Use Restarted Generalized Minimal RESidual iteration to solve ``Ax = b``.

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
    x0: {array, matrix}
        Starting guess for the solution (a vector of zeros by default).
    restart: int, optional
        Number of iterations between restarts. Larger values increase
        iteration cost, but may be necessary for convergence.
        Default is 30.
    tol: float, optional
        Tolerances for convergence, ``norm(residual) <= tol*norm(b)``.
    maxit: int, optional
        Maximum number of iterations (restart cycles).  Iteration will stop
        after maxit steps even if the specified tolerance has not been
        achieved.
	dtype: float, optional
		Default is float64.
	sparams: string, optional
		The solver parameter option is used to output iteration information.
		For example, if check the convergence history, let ``sparams = '-scipy_ksp_monitor'``
		Default is None
    M : {sparse matrix, dense matrix, LinearOperator}
        Inverse of the preconditioner of A.  M should approximate the
        inverse of A and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.

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
	>>> from scipy.sparse import csr_matrix
	>>> from scipy.sparse.linalg import gmres
	>>> A = csr_matrix([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
	>>> b = np.array([2, 4, -1], dtype=float)
	>>> x, exitCode = gmres(A, b)
	>>> print(exitCode)            # 0 indicates successful convergence
    0
	>>> np.allclose(A.dot(x), b)
	True
	"""

	assert A.shape[0] == A.shape[1] & A.shape[0] == len(b)
	cpuKSPSolve = time.time()
	A, M, x, b, postprocess = make_system(A, M, x0, b)
	ndofs = A.shape[0]
	if x0 is None:
		x0 = x.copy()
	if restart is None:
		restart = 30
	if dtype is None:
		dtype = 'float64'
	m = restart
	# define orthonormal basis
	V = np.zeros((m+1, ndofs), dtype=dtype)
	# define Hessenberg matrix
	Hn = np.zeros((m+1, m), dtype=dtype)
	# define rotation matrix
	rot = np.zeros((2, m), dtype=dtype)

	g = np.zeros((m+1), dtype=dtype)
	gtmp = np.zeros((m+1), dtype=dtype)
	# Outer iterations
	for iter in range(maxit):
		r = b - A.matvec(x)
		z = M.matvec(r)
		g[0] = np.linalg.norm(z)  # initial residual 
		if iter == 0:
			r0norm = g[0]
		V[0] = (1./g[0]) * z  # first basis vector
		# Inner iterations
		for k in range(m):
			v = A.matvec(V[k])
			w = M.matvec(v)
			Hn[0:k+1, k] = V[0:k+1].dot(w)
			w -= V[0:k+1].T.dot(Hn[0:k+1, k])
			Hn[k+1, k] = np.linalg.norm(w)
			V[k+1] = (1./Hn[k+1, k]) * w
			# QR decomposition of Hn
			for i in range(k):
				tmp0 = rot[0, i] * Hn[i, k] + rot[1, i] * Hn[i+1, k]
				tmp1 = -rot[1, i] * Hn[i, k] + rot[0, i] * Hn[i+1, k]
				Hn[i, k] = tmp0
				Hn[i+1, k] = tmp1
			sq = np.sqrt(Hn[k+1, k]**2 + Hn[k, k]**2)
			rot[0, k] = Hn[k, k] / sq
			rot[1, k] = Hn[k+1, k] / sq
			Hn[k, k] = sq
			Hn[k+1, k] = 0.
			g[k+1] = -rot[1, k] * g[k]
			g[k] *= rot[0, k]
			# compute relative residual
			rres = np.abs(g[k+1]) / r0norm
			if k < m-1 and rres >= tol and sparams == '-scipy_ksp_monitor':
				print("{}: unpreconditioned residual = {}".format(iter*m + (k+1), rres))
			if rres < tol:
				# reconstruct the solution
				y = np.zeros((k+1), dtype=dtype)
				gtmp[0:k+1] = g[0:k+1]
				for i in range(k, -1, -1):
					gtmp[i] -= Hn[i, i+1:k+1].dot(y[i+1:k+1])
					y[i] = gtmp[i] / Hn[i, i]
				x += V[0:k+1].T.dot(y[0:k+1])
				print("{}: unpreconditioned residual = {}".format(iter*m + (k+1), rres))
				print("Linear solve converged due to reaching TOL iterations {}".format(iter*m + k+1))
				print(" --- system solved with SciPy/GMRES (in {} s)".format(time.time()-cpuKSPSolve))
				return (postprocess(x), 0)
		# reconstruct the solution
		y = np.zeros((m), dtype=dtype)
		gtmp[0:m] = g[0:m]
		for i in range(m-1, -1, -1):
			gtmp[i] -= Hn[i, i+1:m].dot(y[i+1:m])
			y[i] = gtmp[i] / Hn[i, i]
		x += V[0:m].T.dot(y[0:m])
		# compute relative residual
		r = b - A.matvec(x)
		g[k+1] = np.linalg.norm(r)
		rres = g[k+1] / r0norm
		if sparams == '-scipy_ksp_monitor':
			print("{}: unpreconditioned residual = {}".format((iter+1)*m, rres))
		if rres < tol:
			print("Linear solve converged due to reaching TOL iterations {}".format((iter+1)*m))
			print(" --- system solved with SciPy/GMRES (in {} s)".format(time.time()-cpuKSPSolve))
			return (postprocess(x), 0)
		elif iter == maxit - 1:
			print("Linear solve converged due to reaching MAXIT iterations {}".format((iter+1)*m))
			return (postprocess(x), 1)


#======================
#   QMR (Original)
#======================
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
	
    cpuKSPSolve = time.time()
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
    print("Linear solve converged due to reaching TOL iterations {}".format(iter_))
    print(" --- system solved with SciPy/QMR (in {} s)".format(time.time()-cpuKSPSolve))

    return postprocess(x), info
