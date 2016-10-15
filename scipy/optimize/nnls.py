from __future__ import division, print_function, absolute_import

from . import _nnls
from ..linalg import norm
import numpy as np

__all__ = ['nnls']


def nnls(A, b, method="lawson", x0=None, tol=1e-7, n_iter=10000,
         return_n_iter=False):
    """Non-negative linear least squares solver.

    This function finds and returns a vector :math`x` of non-negative real
    numbers that minimizes the linear least-squares objective

    ..math::
        \| Ax - b \|_2

    Parameters
    ----------
    A : ndarray
        Matrix ``A`` as shown above.
    b : ndarray
        Right-hand side vector.
    method : string
        Either "lawson" for the Lawson-Hanson algorithm, or "sca" for the
        sequential coordinate-wise algorithm. See Notes, below.
    x0 : ndarray
        Initial guess at the solution. Only applicable to ``method="sca"``.
    tol : float
        Tolerance when determining convergence when ``method="sca"``.
        Convergence is determined by checking, for each component of ``x``,
        whether it meets the KKT conditions to within ``tol``.
    n_iter : int
        Maximum number of iterations for ``method="sca"``. The number of
        iterations in the Lawson-Hanson algorithm is hardwired.
    return_n_iter : bool
        Whether to return the actual number of iterations needed.
        Only applicable to ``method="sca"``.

    Returns
    -------
    x : ndarray
        Solution vector.
    residual : float
        The residual norm, :math:`\| Ax-b \|_2`.
    n_iter : int
        Actual number of iterations needed to reach the solution, if
        ``return_n_iter=True`` was given.

    Notes
    -----
    The algorithm is that of Franc et al. [1].

    References
    ----------
    .. [1] V. Franc, M. Navara, V. Hlavac (2005). Sequential coordinate-wise
        algorithm for non-negative least squares problem. Proc. Int'l Conf. on
        Computer Analysis of Images and Patterns (CIAP).
        http://cmp.felk.cvut.cz/ftp/articles/franc/Franc-TR-2005-06.pdf

    """

    A, b = map(np.asarray_chkfinite, (A, b))

    if len(A.shape) != 2:
        raise ValueError("expected matrix")
    if len(b.shape) != 1:
        raise ValueError("expected vector")

    if method == "lawson":
        if x0 is not None or return_n_iter:
            raise TypeError("unexpected argument")
        return _lawson(A, b)
    elif method == "sca":
        return _sca(A, b, x0, tol, n_iter, return_n_iter)
    else:
        raise ValueError("unknown method %r" % method)


# Wrapper for nnls.f.
def _lawson(A, b):
    m, n = A.shape

    if m != b.shape[0]:
        raise ValueError("incompatible dimensions")

    w = np.zeros((n,), dtype=np.double)
    zz = np.zeros((m,), dtype=np.double)
    index = np.zeros((n,), dtype=int)

    x, rnorm, mode = _nnls.nnls(A, m, n, b, w, zz, index)
    if mode != 1:
        raise RuntimeError("too many iterations")

    return x, rnorm


def _sca(A, b, x0=None, tol=1e-7, n_iter=1000, return_n_iter=False):
    if x0 is not None:
        x = x0
        if np.any(x0 < 0):
            raise ValueError("x0 should be non-negative")
    else:
        x = np.zeros(A.shape[1])

    f = -np.dot(A.T, b)
    grad = f.copy()
    H = np.dot(A.T, A)

    for it in range(n_iter):
        for k in range(x.shape[0]):
            xk = x[k]
            x[k] = max(0, xk - grad[k] / H[k, k])
            grad += (x[k] - xk) * H[:, k]
        if _converged(H, x, f, tol):
            break

    residual = norm(A.dot(x) - b)
    if return_n_iter:
        return x, residual, it + 1
    else:
        return x, residual


def _converged(H, x, f, tol):
    # Evaluate the relaxed KKT conditions for the QP version of NNLS.
    Hx_f = H.dot(x) + f
    return (Hx_f >= -tol).all() and (Hx_f[x > 0] <= tol).all()
