"""SVD decomposition functions."""
from __future__ import division, print_function, absolute_import

import numpy
from numpy import zeros, r_, diag

# Local imports.
from .misc import LinAlgError, _datacopied
from .lapack import get_lapack_funcs, _compute_lwork
from .decomp import _asarray_validated
from scipy._lib.six import string_types

__all__ = ['svd', 'svdvals', 'diagsvd', 'orth']


def svd(a, full_matrices=True, compute_uv=True, overwrite_a=False,
        check_finite=True, lapack_driver='gesdd'):
    """
    Singular Value Decomposition.

    Factorizes the matrix a into two unitary matrices U and Vh, and
    a 1-D array s of singular values (real, non-negative) such that
    ``a == U*S*Vh``, where S is a suitably shaped matrix of zeros with
    main diagonal s.

    Parameters
    ----------
    a : (M, N) array_like
        Matrix to decompose.
    full_matrices : bool, optional
        If True, `U` and `Vh` are of shape ``(M,M)``, ``(N,N)``.
        If False, the shapes are ``(M,K)`` and ``(K,N)``, where
        ``K = min(M,N)``.
    compute_uv : bool, optional
        Whether to compute also `U` and `Vh` in addition to `s`.
        Default is True.
    overwrite_a : bool, optional
        Whether to overwrite `a`; may improve performance.
        Default is False.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
    lapack_driver : {'gesdd', 'gesvd'}, optional
        Whether to use the more efficient divide-and-conquer approach
        (``'gesdd'``) or general rectangular approach (``'gesvd'``)
        to compute the SVD. MATLAB and Octave use the ``'gesvd'`` approach.
        Default is ``'gesdd'``.

        .. versionadded:: 0.18

    Returns
    -------
    U : ndarray
        Unitary matrix having left singular vectors as columns.
        Of shape ``(M,M)`` or ``(M,K)``, depending on `full_matrices`.
    s : ndarray
        The singular values, sorted in non-increasing order.
        Of shape (K,), with ``K = min(M, N)``.
    Vh : ndarray
        Unitary matrix having right singular vectors as rows.
        Of shape ``(N,N)`` or ``(K,N)`` depending on `full_matrices`.

    For ``compute_uv=False``, only `s` is returned.

    Raises
    ------
    LinAlgError
        If SVD computation does not converge.

    See also
    --------
    svdvals : Compute singular values of a matrix.
    diagsvd : Construct the Sigma matrix, given the vector s.

    Examples
    --------
    >>> from scipy import linalg
    >>> a = np.random.randn(9, 6) + 1.j*np.random.randn(9, 6)
    >>> U, s, Vh = linalg.svd(a)
    >>> U.shape, Vh.shape, s.shape
    ((9, 9), (6, 6), (6,))

    >>> U, s, Vh = linalg.svd(a, full_matrices=False)
    >>> U.shape, Vh.shape, s.shape
    ((9, 6), (6, 6), (6,))
    >>> S = linalg.diagsvd(s, 6, 6)
    >>> np.allclose(a, np.dot(U, np.dot(S, Vh)))
    True

    >>> s2 = linalg.svd(a, compute_uv=False)
    >>> np.allclose(s, s2)
    True

    """
    a1 = _asarray_validated(a, check_finite=check_finite)
    if len(a1.shape) != 2:
        raise ValueError('expected matrix')
    m, n = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))

    if not isinstance(lapack_driver, string_types):
        raise TypeError('lapack_driver must be a string')
    if lapack_driver not in ('gesdd', 'gesvd'):
        raise ValueError('lapack_driver must be "gesdd" or "gesvd", not "%s"'
                         % (lapack_driver,))
    funcs = (lapack_driver, lapack_driver + '_lwork')
    gesXd, gesXd_lwork = get_lapack_funcs(funcs, (a1,))

    # compute optimal lwork
    lwork = _compute_lwork(gesXd_lwork, a1.shape[0], a1.shape[1],
                           compute_uv=compute_uv, full_matrices=full_matrices)

    # perform decomposition
    u, s, v, info = gesXd(a1, compute_uv=compute_uv, lwork=lwork,
                          full_matrices=full_matrices, overwrite_a=overwrite_a)

    if info > 0:
        raise LinAlgError("SVD did not converge")
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal gesdd'
                         % -info)
    if compute_uv:
        return u, s, v
    else:
        return s


def svdvals(a, overwrite_a=False, check_finite=True):
    """
    Compute singular values of a matrix.

    Parameters
    ----------
    a : (M, N) array_like
        Matrix to decompose.
    overwrite_a : bool, optional
        Whether to overwrite `a`; may improve performance.
        Default is False.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    s : (min(M, N),) ndarray
        The singular values, sorted in decreasing order.

    Raises
    ------
    LinAlgError
        If SVD computation does not converge.

    Notes
    -----
    ``svdvals(a)`` only differs from ``svd(a, compute_uv=False)`` by its
    handling of the edge case of empty ``a``, where it returns an
    empty sequence:

    >>> a = np.empty((0, 2))
    >>> from scipy.linalg import svdvals
    >>> svdvals(a)
    array([], dtype=float64)

    See Also
    --------
    svd : Compute the full singular value decomposition of a matrix.
    diagsvd : Construct the Sigma matrix, given the vector s.

    Examples
    --------
    >>> from scipy.linalg import svdvals
    >>> m = np.array([[1.0, 0.0],
    ...               [2.0, 3.0],
    ...               [1.0, 1.0],
    ...               [0.0, 2.0],
    ...               [1.0, 0.0]])
    >>> svdvals(m)
    array([ 4.28091555,  1.63516424])

    We can verify the maximum singular value of `m` by computing the maximum
    length of `m.dot(u)` over all the unit vectors `u` in the (x,y) plane.
    We approximate "all" the unit vectors with a large sample.  Because
    of linearity, we only need the unit vectors with angles in [0, pi].

    >>> t = np.linspace(0, np.pi, 2000)
    >>> u = np.array([np.cos(t), np.sin(t)])
    >>> np.linalg.norm(m.dot(u), axis=0).max()
    4.2809152422538475

    `p` is a projection matrix with rank 1.  With exact arithmetic,
    its singular values would be [1, 0, 0, 0].

    >>> v = np.array([0.1, 0.3, 0.9, 0.3])
    >>> p = np.outer(v, v)
    >>> svdvals(p)
    array([  1.00000000e+00,   2.02021698e-17,   1.56692500e-17,
             8.15115104e-34])

    The singular values of an orthogonal matrix are all 1.  Here we
    create a random orthogonal matrix by using the `rvs()` method of
    `scipy.stats.ortho_group`.

    >>> from scipy.stats import ortho_group
    >>> np.random.seed(123)
    >>> orth = ortho_group.rvs(4)
    >>> svdvals(orth)
    array([ 1.,  1.,  1.,  1.])

    """
    a = _asarray_validated(a, check_finite=check_finite)
    if a.size:
        return svd(a, compute_uv=0, overwrite_a=overwrite_a,
                   check_finite=False)
    elif len(a.shape) != 2:
        raise ValueError('expected matrix')
    else:
        return numpy.empty(0)


def diagsvd(s, M, N):
    """
    Construct the sigma matrix in SVD from singular values and size M, N.

    Parameters
    ----------
    s : (M,) or (N,) array_like
        Singular values
    M : int
        Size of the matrix whose singular values are `s`.
    N : int
        Size of the matrix whose singular values are `s`.

    Returns
    -------
    S : (M, N) ndarray
        The S-matrix in the singular value decomposition

    """
    part = diag(s)
    typ = part.dtype.char
    MorN = len(s)
    if MorN == M:
        return r_['-1', part, zeros((M, N-M), typ)]
    elif MorN == N:
        return r_[part, zeros((M-N, N), typ)]
    else:
        raise ValueError("Length of s must be M or N.")


# Orthonormal decomposition

def orth(A):
    """
    Construct an orthonormal basis for the range of A using SVD

    Parameters
    ----------
    A : (M, N) array_like
        Input array

    Returns
    -------
    Q : (M, K) ndarray
        Orthonormal basis for the range of A.
        K = effective rank of A, as determined by automatic cutoff

    See also
    --------
    svd : Singular value decomposition of a matrix

    """
    u, s, vh = svd(A, full_matrices=False)
    M, N = A.shape
    eps = numpy.finfo(float).eps
    tol = max(M, N) * numpy.amax(s) * eps
    num = numpy.sum(s > tol, dtype=int)
    Q = u[:, :num]
    return Q
