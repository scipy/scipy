import os
import numpy as np

from .arpack import _arpack  # type: ignore[attr-defined]
from . import eigsh

from scipy._lib._util import check_random_state
from scipy.sparse.linalg._interface import LinearOperator, aslinearoperator
from scipy.sparse.linalg._eigen.lobpcg import lobpcg  # type: ignore[no-redef]
if os.environ.get("USE_PROPACK"):
    from scipy.sparse.linalg._svdp import _svdp
    HAS_PROPACK = True
else:
    HAS_PROPACK = False

arpack_int = _arpack.timing.nbx.dtype
__all__ = ['svds']


def _augmented_orthonormal_cols(x, k, random_state):
    # extract the shape of the x array
    n, m = x.shape
    # create the expanded array and copy x into it
    y = np.empty((n, m+k), dtype=x.dtype)
    y[:, :m] = x
    # do some modified gram schmidt to add k random orthonormal vectors
    for i in range(k):
        # sample a random initial vector
        v = random_state.standard_normal(size=n)
        if np.iscomplexobj(x):
            v = v + 1j*random_state.standard_normal(size=n)
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


def _augmented_orthonormal_rows(x, k, random_state):
    return _augmented_orthonormal_cols(x.T, k, random_state).T


def _herm(x):
    return x.T.conj()


def _iv(A, k, ncv, tol, which, v0, maxiter,
        return_singular, solver, random_state):

    # input validation/standardization for `solver`
    # out of order because it's needed for other parameters
    solver = str(solver).lower()
    solvers = {"arpack", "lobpcg", "propack"}
    if solver not in solvers:
        raise ValueError(f"solver must be one of {solvers}.")

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
    kmax = min(A.shape) if solver == 'propack' else min(A.shape) - 1
    if int(k) != k or not (0 < k <= kmax):
        message = "`k` must be an integer satisfying `0 < k < min(A.shape)`."
        raise ValueError(message)
    k = int(k)

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
    if which not in whichs:
        raise ValueError(f"`which` must be in {whichs}.")

    # input validation/standardization for `v0`
    if v0 is not None:
        v0 = np.atleast_1d(v0)
        if not (np.issubdtype(v0.dtype, np.complexfloating)
                or np.issubdtype(v0.dtype, np.floating)):
            message = ("`v0` must be of floating or complex floating "
                       "data type.")
            raise ValueError(message)

        shape = (A.shape[0],) if solver == 'propack' else (min(A.shape),)
        if v0.shape != shape:
            message = "`v0` must have shape {shape}."
            raise ValueError(message)

    # input validation/standardization for `maxiter`
    if maxiter is not None and (int(maxiter) != maxiter or maxiter <= 0):
        message = "`maxiter` must be a positive integer."
        raise ValueError(message)
    maxiter = int(maxiter) if maxiter is not None else maxiter

    # input validation/standardization for `return_singular_vectors`
    # not going to be flexible with this; too complicated for little gain
    rs_options = {True, False, "vh", "u"}
    if return_singular not in rs_options:
        raise ValueError(f"`return_singular_vectors` must be in {rs_options}.")

    random_state = check_random_state(random_state)

    return (A, k, ncv, tol, which, v0, maxiter,
            return_singular, solver, random_state)


def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True,
         solver='arpack', random_state=None, options=None):
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
        Must satisfy ``1 <= k <= kmax``, where ``kmax=min(M, N)`` for
        ``solver='propack'`` and ``kmax=min(M, N) - 1`` otherwise.
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
    return_singular_vectors : {True, False, "u", "vh"}
        Singular values are always computed and returned; this parameter
        controls the computation and return of singular vectors.

        - ``True``: return singular vectors.
        - ``False``: do not return singular vectors.
        - ``"u"``: if ``M <= N``, compute only the left singular vectors and
          return ``None`` for the right singular vectors. Otherwise, compute
          all singular vectors.
        - ``"vh"``: if ``M > N``, compute only the right singular vectors and
          return ``None`` for the left singular vectors. Otherwise, compute
          all singular vectors.

        If ``solver='propack'``, the option is respected regardless of the
        matrix shape.

    solver :  {'arpack', 'propack', 'lobpcg'}, optional
            The solver used.
            :ref:`'arpack' <sparse.linalg.svds-arpack>`,
            :ref:`'lobpcg' <sparse.linalg.svds-lobpcg>`, and
            :ref:`'propack' <sparse.linalg.svds-propack>` are supported.
            Default: `'arpack'`.
    random_state : {None, int, `numpy.random.Generator`,
                    `numpy.random.RandomState`}, optional

        Pseudorandom number generator state used to generate resamples.

        If `random_state` is ``None`` (or `np.random`), the
        `numpy.random.RandomState` singleton is used.
        If `random_state` is an int, a new ``RandomState`` instance is used,
        seeded with `random_state`.
        If `random_state` is already a ``Generator`` or ``RandomState``
        instance then that instance is used.
    options : dict, optional
        A dictionary of solver-specific options. No solver-specific options
        are currently supported; this parameter is reserved for future use.

    Returns
    -------
    u : ndarray, shape=(M, k)
        Unitary matrix having left singular vectors as columns.
    s : ndarray, shape=(k,)
        The singular values.
    vh : ndarray, shape=(k, N)
        Unitary matrix having right singular vectors as rows.

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
    >>> np.allclose(A2, A.toarray(), atol=1e-3)
    True

    With all five singular values/vectors, we can reproduce the original
    matrix.

    >>> u3, s3, vT3 = svds(A, k=5)
    >>> A3 = u3 @ np.diag(s3) @ vT3
    >>> np.allclose(A3, A.toarray())
    True

    The singular values match the expected singular values, and the singular
    vectors are as expected up to a difference in sign.

    >>> (np.allclose(s3, s) and
    ...  np.allclose(np.abs(u3), np.abs(u.toarray())) and
    ...  np.allclose(np.abs(vT3), np.abs(vT.toarray())))
    True

    The singular vectors are also orthogonal.
    >>> (np.allclose(u3.T @ u3, np.eye(5)) and
    ...  np.allclose(vT3 @ vT3.T, np.eye(5)))
    True

    """
    rs_was_None = random_state is None  # avoid changing v0 for arpack/lobpcg

    args = _iv(A, k, ncv, tol, which, v0, maxiter, return_singular_vectors,
               solver, random_state)
    (A, k, ncv, tol, which, v0, maxiter,
     return_singular_vectors, solver, random_state) = args

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
            if rs_was_None:
                X = np.random.RandomState(52).randn(min(A.shape), k)
            else:
                X = random_state.uniform(size=(min(A.shape), k))

        eigvals, eigvec = lobpcg(XH_X, X, tol=tol ** 2, maxiter=maxiter,
                                 largest=largest, )

    elif solver == 'propack':
        if not HAS_PROPACK:
            raise ValueError("`solver='propack'` is opt-in due to potential issues on Windows, "
                             "it can be enabled by setting the `USE_PROPACK` environment "
                             "variable before importing scipy")
        jobu = return_singular_vectors in {True, 'u'}
        jobv = return_singular_vectors in {True, 'vh'}
        irl_mode = (which == 'SM')
        res = _svdp(A, k=k, tol=tol**2, which=which, maxiter=None,
                    compute_u=jobu, compute_v=jobv, irl_mode=irl_mode,
                    kmax=maxiter, v0=v0, random_state=random_state)

        u, s, vh, _ = res  # but we'll ignore bnd, the last output

        # PROPACK order appears to be largest first. `svds` output order is not
        # guaranteed, according to documentation, but for ARPACK and LOBPCG
        # they actually are ordered smallest to largest, so reverse for
        # consistency.
        s = s[::-1]
        u = u[:, ::-1]
        vh = vh[::-1]

        u = u if jobu else None
        vh = vh if jobv else None

        if return_singular_vectors:
            return u, s, vh
        else:
            return s

    elif solver == 'arpack' or solver is None:
        if v0 is None and not rs_was_None:
            v0 = random_state.uniform(size=(min(A.shape),))
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

    u = (_augmented_orthonormal_cols(ularge, nsmall, random_state)
         if ularge is not None else None)
    vh = (_augmented_orthonormal_rows(vhlarge, nsmall, random_state)
          if vhlarge is not None else None)

    indexes_sorted = np.argsort(s)
    s = s[indexes_sorted]
    if u is not None:
        u = u[:, indexes_sorted]
    if vh is not None:
        vh = vh[indexes_sorted]

    return u, s, vh
