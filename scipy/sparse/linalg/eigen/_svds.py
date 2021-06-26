import numpy as np

from .arpack import _arpack  # type: ignore[attr-defined]
from . import eigsh

from scipy.sparse.linalg.interface import LinearOperator, aslinearoperator
from scipy.sparse import isspmatrix
from scipy.sparse.sputils import is_pydata_spmatrix
from scipy.sparse.linalg.eigen.lobpcg import lobpcg  # type: ignore[no-redef]

arpack_int = _arpack.timing.nbx.dtype
__all__ = ['svds']


def _augmented_orthonormal_cols(x, k):
    # extract the shape of the x array
    n, m = x.shape
    # create the expanded array and copy x into it
    y = np.empty((n, m+k), dtype=x.dtype)
    y[:, :m] = x
    # do some modified gram schmidt to add k random orthonormal vectors
    for i in range(k):
        # sample a random initial vector
        v = np.random.randn(n)
        if np.iscomplexobj(x):
            v = v + 1j*np.random.randn(n)
        # subtract projections onto the existing unit length vectors
        for j in range(m+i):
            u = y[:, j]
            v -= (np.dot(v, u.conj()) / np.dot(u, u.conj())) * u
        # normalize v
        v /= np.sqrt(np.dot(v, v.conj()))
        # add v into the output array
        y[:, m+i] = v
    # return the expanded array
    return y


def _augmented_orthonormal_rows(x, k):
    return _augmented_orthonormal_cols(x.T, k).T


def _herm(x):
    return x.T.conj()


def _iv(A, k, ncv, tol, which, v0, maxiter, return_singular, solver):

    # input validation/standardization for `A`
    A = aslinearoperator(A)  # this takes care of some input validation
    if not (np.issubdtype(A.dtype, np.complexfloating)
            or np.issubdtype(A.dtype, np.floating)):
        message = "`A` must be of floating or complex floating data type."
        raise ValueError(message)
    if np.prod(A.shape) == 0:
        message = "`A` must not be empty."
        raise ValueError(message)

    # input validation/standardization for `k`
    if int(k) != k or not (0 < k < min(A.shape)):
        message = "`k` must be an integer satisfying `0 < k < min(A.shape)`."
        raise ValueError(message)
    k = int(k)

    # input validation/standardization for `solver`
    # out of order because it's needed for `ncv`
    solver = str(solver).lower()
    solvers = {"arpack", "lobpcg", "propack"}
    if solver not in solvers:
        raise ValueError(f"solver must be one of {solvers}.")

    # input validation/standardization for `ncv`
    if solver == "arpack" and ncv is not None:
        if int(ncv) != ncv or not (k < ncv < min(A.shape)):
            message = ("`ncv` must be an integer satisfying "
                       "`k < ncv < min(A.shape)`.")
            raise ValueError(message)
        ncv = int(ncv)

    # input validation/standardization for `tol`
    if tol < 0 or not np.isfinite(tol):
        message = "`tol` must be a non-negative floating point value."
        raise ValueError(message)
    tol = float(tol)

    # input validation/standardization for `which`
    which = str(which).upper()
    whichs = {'LM', 'SM'}
    if which not in {'LM', 'SM'}:
        raise ValueError(f"`which` must be in {whichs}.")

    # input validation/standardization for `v0`
    if v0 is not None:
        v0 = np.atleast_1d(v0)
        if not (np.issubdtype(v0.dtype, np.complexfloating)
                or np.issubdtype(v0.dtype, np.floating)):
            message =  ("`v0` must be of floating or complex floating "
                        "data type.")
            raise ValueError(message)
        if v0.shape != (min(A.shape),):
            message = "`v0` must have shape `(min(A.shape),)`."
            raise ValueError(message)

    return A, k, ncv, tol, which, v0, maxiter, return_singular, solver


def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True,
         solver='arpack', options=None):
    """
    Partial singular value decomposition of a sparse matrix.

    Compute the largest or smallest `k` singular values and corresponding
    singular vectors of a sparse matrix `A`. The order in which the singular
    values are returned is not guaranteed.

    In the descriptions below, let ``M, N = A.shape``.

    Parameters
    ----------
    A : sparse matrix or LinearOperator
        Matrix to decompose.
    k : int, default: 6
        Number of singular values and singular vectors to compute.
        Must satisfy ``1 <= k < min(M, N)``.
    ncv : int, optional
        When ``solver='arpack'``, this is the number of Lanczos vectors
        generated. See :ref:`'arpack' <sparse.linalg.svds-arpack>` for details.
        When ``solver='lobpcg'`` or ``solver='propack'``, this parameter is
        ignored.
    tol : float, optional
        Tolerance for singular values. Zero (default) means machine precision.
    which : {'LM', 'SM'}
        Which `k` singular values to find: either the largest magnitude ('LM')
        or smallest magnitude ('SM') singular values.
    v0 : ndarray, optional
        The starting vector for iteration; see method-specific
        documentation (:ref:`'arpack' <sparse.linalg.svds-arpack>`,
        :ref:`'lobpcg' <sparse.linalg.svds-lobpcg>`), or
        :ref:`'propack' <sparse.linalg.svds-propack>` for details.
    maxiter : int, optional
        Maximum number of iterations; see method-specific
        documentation (:ref:`'arpack' <sparse.linalg.svds-arpack>`,
        :ref:`'lobpcg' <sparse.linalg.svds-lobpcg>`), or
        :ref:`'propack' <sparse.linalg.svds-propack>` for details.
    return_singular_vectors : bool or str, optional
        Singular values are always computed and returned; this parameter
        controls the computation and return of singular vectors.

        - ``True``: return singular vectors.
        - ``False``: do not return singular vectors.
        - ``"u"``: only return the left singular values, without computing the
          right singular vectors (if ``N > M``).
        - ``"vh"``: only return the right singular values, without computing
          the left singular vectors (if ``N <= M``).

    solver : str, optional
            The solver used.
            :ref:`'arpack' <sparse.linalg.svds-arpack>`,
            :ref:`'lobpcg' <sparse.linalg.svds-lobpcg>`, and
            :ref:`'propack' <sparse.linalg.svds-propack>` are supported.
            Default: `'arpack'`.
    options : dict, optional
        A dictionary of solver-specific options. No solver-specific options
        are currently supported; this parameter is reserved for future use.

    Returns
    -------
    u : ndarray, shape=(M, k)
        Unitary matrix having left singular vectors as columns.
        If `return_singular_vectors` is ``"vh"``, this variable is not
        computed, and ``None`` is returned instead.
    s : ndarray, shape=(k,)
        The singular values.
    vh : ndarray, shape=(k, N)
        Unitary matrix having right singular vectors as rows.
        If `return_singular_vectors` is ``"u"``, this variable is not computed,
        and ``None`` is returned instead.

    Notes
    -----
    This is a naive implementation using ARPACK or LOBPCG as an eigensolver
    on ``A.conj().T @ A`` or ``A @ A.conj().T``, depending on which one is more
    efficient.

    Examples
    --------
    Construct a matrix ``A`` from singular values and vectors.

    >>> from scipy.stats import ortho_group
    >>> from scipy.sparse import csc_matrix, diags
    >>> from scipy.sparse.linalg import svds
    >>> rng = np.random.default_rng()
    >>> orthogonal = csc_matrix(ortho_group.rvs(10, random_state=rng))
    >>> s = [0.0001, 0.001, 3, 4, 5]  # singular values
    >>> u = orthogonal[:, :5]         # left singular vectors
    >>> vT = orthogonal[:, 5:].T      # right singular vectors
    >>> A = u @ diags(s) @ vT

    With only three singular values/vectors, the SVD approximates the original
    matrix.

    >>> u2, s2, vT2 = svds(A, k=3)
    >>> A2 = u2 @ np.diag(s2) @ vT2
    >>> np.allclose(A2, A.todense(), atol=1e-3)
    True

    With all five singular values/vectors, we can reproduce the original
    matrix.

    >>> u3, s3, vT3 = svds(A, k=5)
    >>> A3 = u3 @ np.diag(s3) @ vT3
    >>> np.allclose(A3, A.todense())
    True

    The singular values match the expected singular values, and the singular
    values are as expected up to a difference in sign. Consequently, the
    returned arrays of singular vectors must also be orthogonal.

    >>> (np.allclose(s3, s) and
    ...  np.allclose(np.abs(u3), np.abs(u.todense())) and
    ...  np.allclose(np.abs(vT3), np.abs(vT.todense())))
    True

    """

    args = _iv(A, k, ncv, tol, which, v0, maxiter, return_singular_vectors,
               solver)
    A, k, ncv, tol, which, v0, maxiter, return_singular_vectors, solver = args

    largest = (which == 'LM')
    n, m = A.shape

    if n > m:
        X_dot = A.matvec
        X_matmat = A.matmat
        XH_dot = A.rmatvec
        XH_mat = A.rmatmat
    else:
        X_dot = A.rmatvec
        X_matmat = A.rmatmat
        XH_dot = A.matvec
        XH_mat = A.matmat

        dtype = getattr(A, 'dtype', None)
        if dtype is None:
            dtype = A.dot(np.zeros([m, 1])).dtype

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

        eigvals, eigvec = lobpcg(XH_X, X, tol=tol ** 2, maxiter=maxiter,
                                 largest=largest)

    elif solver == 'propack':
        from scipy.sparse.linalg import svdp
        jobu = n <= m or return_singular_vectors != 'vh'
        jobv = n > m or return_singular_vectors != 'u'
        res = svdp(A, k=k, tol=tol**2, which=which, maxiter=None,
                   compute_u=jobu, compute_v=jobv, irl_mode=True,
                   kmax=maxiter, v0=v0)

        u, s, vh, _ = res  # but we'll ignore bnd, the last output

        # PROPACK order appears to be largest first. Output order is not
        # guaranteed, according to documentation, but for ARPACK and LOBPCG
        # they actually are ordered smallest to largest.
        s = s[::-1]
        u = u[:, ::-1]
        vh = vh[::-1]

        if return_singular_vectors is False:
            return s
        elif return_singular_vectors == 'u' and n <= m:
            return u, s, None
        elif return_singular_vectors == 'vh' and n >= m:
            return None, s, vh
        else:
            return u, s, vh

    elif solver == 'arpack' or solver is None:
        eigvals, eigvec = eigsh(XH_X, k=k, tol=tol ** 2, maxiter=maxiter,
                                ncv=ncv, which=which, v0=v0)

    # Gramian matrices have real non-negative eigenvalues.
    eigvals = np.maximum(eigvals.real, 0)

    # Use the sophisticated detection of small eigenvalues from pinvh.
    t = eigvec.dtype.char.lower()
    factor = {'f': 1E3, 'd': 1E6}
    cond = factor[t] * np.finfo(t).eps
    cutoff = cond * np.max(eigvals)

    # Get a mask indicating which eigenpairs are not degenerately tiny,
    # and create the re-ordered array of thresholded singular values.
    above_cutoff = (eigvals > cutoff)
    nlarge = above_cutoff.sum()
    nsmall = k - nlarge
    slarge = np.sqrt(eigvals[above_cutoff])
    s = np.zeros_like(eigvals)
    s[:nlarge] = slarge
    if not return_singular_vectors:
        return np.sort(s)

    if n > m:
        vlarge = eigvec[:, above_cutoff]
        ularge = (X_matmat(vlarge) / slarge
                  if return_singular_vectors != 'vh' else None)
        vhlarge = _herm(vlarge)
    else:
        ularge = eigvec[:, above_cutoff]
        vhlarge = (_herm(X_matmat(ularge) / slarge)
                   if return_singular_vectors != 'u' else None)

    u = (_augmented_orthonormal_cols(ularge, nsmall)
         if ularge is not None else None)
    vh = (_augmented_orthonormal_rows(vhlarge, nsmall)
          if vhlarge is not None else None)

    indexes_sorted = np.argsort(s)
    s = s[indexes_sorted]
    if u is not None:
        u = u[:, indexes_sorted]
    if vh is not None:
        vh = vh[indexes_sorted]

    return u, s, vh
