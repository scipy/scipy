
def _svds_arpack_doc(A, k=6, ncv=None, tol=0, which='LM', v0=None,
                     maxiter=None, return_singular_vectors=True,
                     solver='arpack'):
    """
    Partial singular value decomposition of a sparse matrix using ARPACK.

    Compute the largest or smallest `k` singular values and corresponding
    singular vectors of a sparse matrix `A`. The order in which the singular
    values are returned is not guaranteed.

    In the descriptions below, let ``M, N = A.shape``.

    Parameters
    ----------
    A : {sparse matrix, LinearOperator}
        Matrix to decompose.
    k : int, default: 6
        Number of singular values and singular vectors to compute.
        Must satisfy ``1 <= k < min(M, N)``.
    ncv : int, optional
        The number of Lanczos vectors generated.
        The default is ``min(n, max(2*k + 1, 20))``.
        If specified, must satistify ``k + 1 < ncv < N``; ``ncv > 2*k`` is
        recommended.
    tol : float, optional
        Tolerance for singular values. Zero (default) means machine precision.
    which : {'LM', 'SM'}
        Which `k` singular values to find: either the largest magnitude ('LM')
        or smallest magnitude ('SM') singular values.
    v0 : ndarray, optional
        The starting vector for iteration:
        an (approximate) left singular vector if ``N > M`` and a right singular
        vector otherwise. Must be of length ``min(M, N)``.
        Default: random
    maxiter : int, optional
        Maximum number of Arnoldi update iterations allowed;
        default is ``min(M, N) * 10``.
    return_singular_vectors : bool or str, default: True
        Singular values are always computed and returned; this parameter
        controls the computation and return of singular vectors.

        - True: return singular vectors.
        - False: do not return singular vectors.
        - "u": only return the left singular values, without computing the
          right singular vectors (if ``N > M``).
        - "vh": only return the right singular values, without computing the
          left singular vectors (if ``N <= M``).

    solver : str, optional
            This is the solver-specific documentation for ``solver='arpack'``.
            :ref:`'lobpcg' <sparse.linalg.svds-lobpcg>` is also supported.

    Returns
    -------
    u : ndarray, shape=(M, k)
        Unitary matrix having left singular vectors as columns.
        If `return_singular_vectors` is "vh", this variable is not computed,
        and ``None`` is returned instead.
    s : ndarray, shape=(k,)
        The singular values.
    vt : ndarray, shape=(k, n)
        Unitary matrix having right singular vectors as rows.
        If `return_singular_vectors` is "u", this variable is not computed,
        and None is returned instead.


    Notes
    -----
    This is a naive implementation using ARPACK as an eigensolver
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
    pass

def _svds_lobpcg_doc(A, k=6, ncv=None, tol=0, which='LM', v0=None,
                     maxiter=None, return_singular_vectors=True,
                     solver='lobpcg'):
    """
    Partial singular value decomposition of a sparse matrix using LOBPCG.

    Compute the largest or smallest `k` singular values and corresponding
    singular vectors of a sparse matrix `A`. The order in which the singular
    values are returned is not guaranteed.

    In the descriptions below, let ``M, N = A.shape``.

    Parameters
    ----------
    A : {sparse matrix, LinearOperator}
        Matrix to decompose.
    k : int, default: 6
        Number of singular values and singular vectors to compute.
        Must satisfy ``1 <= k < min(M, N)``.
    ncv : int, optional
        Ignored; only used for ``solver='arpack'``.
    tol : float, optional
        Tolerance for singular values. Zero (default) means machine precision.
    which : {'LM', 'SM'}
        Which `k` singular values to find: either the largest magnitude ('LM')
        or smallest magnitude ('SM') singular values.
    v0 : ndarray, optional
        If `k` is 1, the starting vector for iteration:
        an (approximate) left singular vector if ``N > M`` and a right singular
        vector otherwise. Must be of length ``min(M, N)``.
        Ignored otherwise.
        Default: random
    maxiter : int, default: 20
        Maximum number of iterations.
    return_singular_vectors : bool or str, default: True
        Singular values are always computed and returned; this parameter
        controls the computation and return of singular vectors.

        - True: return singular vectors.
        - False: do not return singular vectors.
        - "u": only return the left singular values, without computing the
          right singular vectors (if ``N > M``).
        - "vh": only return the right singular values, without computing the
          left singular vectors (if ``N <= M``).

    solver : str, optional
            This is the solver-specific documentation for ``solver='lobpcg'``.
            :ref:`'arpack' <sparse.linalg.svds-lobpcg>` is also supported.

    Returns
    -------
    u : ndarray, shape=(M, k)
        Unitary matrix having left singular vectors as columns.
        If `return_singular_vectors` is "vh", this variable is not computed,
        and ``None`` is returned instead.
    s : ndarray, shape=(k,)
        The singular values.
    vt : ndarray, shape=(k, n)
        Unitary matrix having right singular vectors as rows.
        If `return_singular_vectors` is "u", this variable is not computed,
        and None is returned instead.


    Notes
    -----
    This is a naive implementation using LOBPCG as an eigensolver
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
    pass