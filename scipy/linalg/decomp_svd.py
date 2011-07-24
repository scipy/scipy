"""SVD decomposition functions."""

import numpy
from numpy import asarray_chkfinite, zeros, r_, diag
from scipy.linalg import calc_lwork

# Local imports.
from misc import LinAlgError, _datacopied
from lapack import get_lapack_funcs

__all__ = ['svd', 'svdvals', 'diagsvd', 'orth']


def svd(a, full_matrices=True, compute_uv=True, overwrite_a=False):
    """Singular Value Decomposition.

    Factorizes the matrix a into two unitary matrices U and Vh and
    an 1d-array s of singular values (real, non-negative) such that
    a == U S Vh  if S is an suitably shaped matrix of zeros whose
    main diagonal is s.

    Parameters
    ----------
    a : array, shape (M, N)
        Matrix to decompose
    full_matrices : boolean
        If true,  U, Vh are shaped  (M,M), (N,N)
        If false, the shapes are    (M,K), (K,N) where K = min(M,N)
    compute_uv : boolean
        Whether to compute also U, Vh in addition to s (Default: true)
    overwrite_a : boolean
        Whether data in a is overwritten (may improve performance)

    Returns
    -------
    U:  array, shape (M,M) or (M,K) depending on full_matrices
    s:  array, shape (K,)
        The singular values, sorted so that s[i] >= s[i+1]. K = min(M, N)
    Vh: array, shape (N,N) or (K,N) depending on full_matrices

    For compute_uv = False, only s is returned.

    Raises LinAlgError if SVD computation does not converge

    Examples
    --------
    >>> from scipy import random, linalg, allclose, dot
    >>> a = random.randn(9, 6) + 1j*random.randn(9, 6)
    >>> U, s, Vh = linalg.svd(a)
    >>> U.shape, Vh.shape, s.shape
    ((9, 9), (6, 6), (6,))

    >>> U, s, Vh = linalg.svd(a, full_matrices=False)
    >>> U.shape, Vh.shape, s.shape
    ((9, 6), (6, 6), (6,))
    >>> S = linalg.diagsvd(s, 6, 6)
    >>> allclose(a, dot(U, dot(S, Vh)))
    True

    >>> s2 = linalg.svd(a, compute_uv=False)
    >>> allclose(s, s2)
    True

    See also
    --------
    svdvals : return singular values of a matrix
    diagsvd : return the Sigma matrix, given the vector s

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError('expected matrix')
    m,n = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))
    gesdd, = get_lapack_funcs(('gesdd',), (a1,))
    if gesdd.module_name[:7] == 'flapack':
        lwork = calc_lwork.gesdd(gesdd.prefix, m, n, compute_uv)[1]
        u,s,v,info = gesdd(a1,compute_uv = compute_uv, lwork = lwork,
                           full_matrices=full_matrices, overwrite_a = overwrite_a)
    else: # 'clapack'
        raise NotImplementedError('calling gesdd from %s' % gesdd.module_name)
    if info > 0:
        raise LinAlgError("SVD did not converge")
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal gesdd'
                                                                    % -info)
    if compute_uv:
        return u, s, v
    else:
        return s

def svdvals(a, overwrite_a=False):
    """Compute singular values of a matrix.

    Parameters
    ----------
    a : array, shape (M, N)
        Matrix to decompose
    overwrite_a : boolean
        Whether data in a is overwritten (may improve performance)

    Returns
    -------
    s:  array, shape (K,)
        The singular values, sorted so that s[i] >= s[i+1]. K = min(M, N)

    Raises LinAlgError if SVD computation does not converge

    See also
    --------
    svd : return the full singular value decomposition of a matrix
    diagsvd : return the Sigma matrix, given the vector s

    """
    return svd(a, compute_uv=0, overwrite_a=overwrite_a)

def diagsvd(s, M, N):
    """Construct the sigma matrix in SVD from singular values and size M,N.

    Parameters
    ----------
    s : array, shape (M,) or (N,)
        Singular values
    M : integer
    N : integer
        Size of the matrix whose singular values are s

    Returns
    -------
    S : array, shape (M, N)
        The S-matrix in the singular value decomposition

    """
    part = diag(s)
    typ = part.dtype.char
    MorN = len(s)
    if MorN == M:
        return r_['-1', part, zeros((M, N-M), typ)]
    elif MorN == N:
        return r_[part, zeros((M-N,N), typ)]
    else:
        raise ValueError("Length of s must be M or N.")


# Orthonormal decomposition

def orth(A):
    """Construct an orthonormal basis for the range of A using SVD

    Parameters
    ----------
    A : array, shape (M, N)

    Returns
    -------
    Q : array, shape (M, K)
        Orthonormal basis for the range of A.
        K = effective rank of A, as determined by automatic cutoff

    See also
    --------
    svd : Singular value decomposition of a matrix

    """
    u, s, vh = svd(A)
    M, N = A.shape
    eps = numpy.finfo(float).eps
    tol = max(M,N) * numpy.amax(s) * eps
    num = numpy.sum(s > tol, dtype=int)
    Q = u[:,:num]
    return Q
