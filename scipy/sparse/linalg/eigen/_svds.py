import numpy as np

from .arpack import _arpack  # type: ignore[attr-defined]
from . import eigsh

from scipy._lib._util import check_random_state
from scipy.sparse.linalg.interface import LinearOperator, aslinearoperator
from scipy.sparse.linalg.eigen.lobpcg import lobpcg  # type: ignore[no-redef]
from scipy.sparse.linalg._svdp import _svdp

arpack_int = _arpack.timing.nbx.dtype
__all__ = ['svds']


def _herm(x):
    return x.T.conj()


def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True,
         solver='arpack'):
    """Compute the largest or smallest k singular values/vectors for a sparse matrix. The order of the singular values is not guaranteed.

    Parameters
    ----------
    A : {sparse matrix, LinearOperator}
        Array to compute the SVD on, of shape (M, N)
    k : int, optional
        Number of singular values and vectors to compute.
        Must be 1 <= k < min(A.shape).
    ncv : int, optional
        The number of Lanczos vectors generated
        ncv must be greater than k+1 and smaller than n;
        it is recommended that ncv > 2*k
        Default: ``min(n, max(2*k + 1, 20))``
    tol : float, optional
        Tolerance for singular values. Zero (default) means machine precision.
    which : str, ['LM' | 'SM'], optional
        Which `k` singular values to find:

            - 'LM' : largest singular values
            - 'SM' : smallest singular values

        .. versionadded:: 0.12.0
    v0 : ndarray, optional
        Starting vector for iteration, of length min(A.shape). Should be an
        (approximate) left singular vector if N > M and a right singular
        vector otherwise.
        Default: random

        .. versionadded:: 0.12.0
    maxiter : int, optional
        Maximum number of iterations.

        .. versionadded:: 0.12.0
    return_singular_vectors : bool or str, optional
        - True: return singular vectors (True) in addition to singular values.

        .. versionadded:: 0.12.0

        - "u": only return the u matrix, without computing vh (if N > M).
        - "vh": only return the vh matrix, without computing u (if N <= M).

        .. versionadded:: 0.16.0
    solver : str, optional
            Eigenvalue solver to use. Should be 'arpack' or 'lobpcg'.
            Default: 'arpack'

    Returns
    -------
    u : ndarray, shape=(M, k)
        Unitary matrix having left singular vectors as columns.
        If `return_singular_vectors` is "vh", this variable is not computed,
        and None is returned instead.
    s : ndarray, shape=(k,)
        The singular values.
    vt : ndarray, shape=(k, N)
        Unitary matrix having right singular vectors as rows.
        If `return_singular_vectors` is "u", this variable is not computed,
        and None is returned instead.


    Notes
    -----
    This is a naive implementation using ARPACK or LOBPCG as an eigensolver
    on A.H * A or A * A.H, depending on which one is more efficient.

    Examples
    --------
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import svds, eigs
    >>> A = csc_matrix([[1, 0, 0], [5, 0, 2], [0, -1, 0], [0, 0, 3]], dtype=float)
    >>> u, s, vt = svds(A, k=2)
    >>> s
    array([ 2.75193379,  5.6059665 ])
    >>> np.sqrt(eigs(A.dot(A.T), k=2)[0]).real
    array([ 5.6059665 ,  2.75193379])
    """
    if which == 'LM':
        largest = True
    elif which == 'SM':
        largest = False
    else:
        raise ValueError("which must be either 'LM' or 'SM'.")

    if not (isinstance(A, LinearOperator) or isspmatrix(A) or is_pydata_spmatrix(A)):
        A = np.asarray(A)

    n, m = A.shape

    if k <= 0 or k >= min(n, m):
        raise ValueError("k must be between 1 and min(A.shape), k=%d" % k)

    if isinstance(A, LinearOperator):
        if n > m:
            X_dot = A.matvec
            X_matmat = A.matmat
            XH_dot = A.rmatvec
            XH_mat = A.rmatmat
            transpose = False
        else:
            X_dot = A.rmatvec
            X_matmat = A.rmatmat
            XH_dot = A.matvec
            XH_mat = A.matmat

            dtype = getattr(A, 'dtype', None)
            if dtype is None:
                dtype = A.dot(np.zeros([m, 1])).dtype
            transpose = True

    else:
        if n > m:
            X_dot = X_matmat = A.dot
            XH_dot = XH_mat = _herm(A).dot
            transpose = False
        else:
            XH_dot = XH_mat = A.dot
            X_dot = X_matmat = _herm(A).dot
            transpose = True

    def matvec_XH_X(x):
        return XH_dot(X_dot(x))

    def matmat_XH_X(x):
        return XH_mat(X_matmat(x))

    XH_X = LinearOperator(matvec=matvec_XH_X, dtype=A.dtype,
                          matmat=matmat_XH_X,
                          shape=(min(A.shape), min(A.shape)))

    # Get a low rank approximation of the implicitly defined gramian matrix.
    # This is not a stable way to approach the problem.
    if solver == 'lobpcg':

        if k == 1 and v0 is not None:
            X = np.reshape(v0, (-1, 1))
        else:
            X = np.random.RandomState(52).randn(min(A.shape), k)

        eigvals, eigvec = lobpcg(XH_X, X, tol=tol, maxiter=maxiter,
                                 largest=largest)

    elif solver == 'arpack' or solver is None:
        eigvals, eigvec = eigsh(XH_X, k=k, tol=tol, maxiter=maxiter,
                                ncv=ncv, which=which, v0=v0)

    else:
        raise ValueError("solver must be either 'arpack', or 'lobpcg'.")

    u = X_matmat(eigvec)
    if not return_singular_vectors:
        s = svd(u, compute_uv=False)
        return s[::-1]

    # compute the right singular vectors of X and update the left ones accordingly
    u, s, vh = svd(u, full_matrices=False)
    u = u[:, ::-1]
    s = s[::-1]
    vh = vh[::-1]
    return_u = (return_singular_vectors == 'u')
    return_vh = (return_singular_vectors == 'vh')
    if not transpose:
        if return_vh:
            u = None
        if return_u:
            vh = None
        else:
            vh = vh @ _herm(eigvec)
        return u, s, vh
    else:
        if return_u:
            u = eigvec @ _herm(vh)
            return u, s, None
        if return_vh:
            return None, s, _herm(u)
        u, vh = eigvec @ _herm(vh), _herm(u)
        return u, s, vh