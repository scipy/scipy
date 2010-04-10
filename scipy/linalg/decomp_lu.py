"""LU decomposition functions."""

from warnings import warn

from numpy import asarray, asarray_chkfinite

# Local imports
from misc import _datanotshared
from lapack import get_lapack_funcs
from flinalg import get_flinalg_funcs


def lu_factor(a, overwrite_a=False):
    """Compute pivoted LU decomposition of a matrix.

    The decomposition is::

        A = P L U

    where P is a permutation matrix, L lower triangular with unit
    diagonal elements, and U upper triangular.

    Parameters
    ----------
    a : array, shape (M, M)
        Matrix to decompose
    overwrite_a : boolean
        Whether to overwrite data in A (may increase performance)

    Returns
    -------
    lu : array, shape (N, N)
        Matrix containing U in its upper triangle, and L in its lower triangle.
        The unit diagonal elements of L are not stored.
    piv : array, shape (N,)
        Pivot indices representing the permutation matrix P:
        row i of matrix was interchanged with row piv[i].

    See also
    --------
    lu_solve : solve an equation system using the LU factorization of a matrix

    Notes
    -----
    This is a wrapper to the *GETRF routines from LAPACK.

    """
    a1 = asarray(a)
    if len(a1.shape) != 2 or (a1.shape[0] != a1.shape[1]):
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1, a))
    getrf, = get_lapack_funcs(('getrf',), (a1,))
    lu, piv, info = getrf(a, overwrite_a=overwrite_a)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of '
                                'internal getrf (lu_factor)' % -info)
    if info > 0:
        warn("Diagonal number %d is exactly zero. Singular matrix." % info,
                    RuntimeWarning)
    return lu, piv


def lu_solve((lu, piv), b, trans=0, overwrite_b=False):
    """Solve an equation system, a x = b, given the LU factorization of a

    Parameters
    ----------
    (lu, piv)
        Factorization of the coefficient matrix a, as given by lu_factor
    b : array
        Right-hand side
    trans : {0, 1, 2}
        Type of system to solve:

        =====  =========
        trans  system
        =====  =========
        0      a x   = b
        1      a^T x = b
        2      a^H x = b
        =====  =========

    Returns
    -------
    x : array
        Solution to the system

    See also
    --------
    lu_factor : LU factorize a matrix

    """
    b1 = asarray_chkfinite(b)
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b, '__array__'))
    if lu.shape[0] != b1.shape[0]:
        raise ValueError("incompatible dimensions.")

    getrs, = get_lapack_funcs(('getrs',), (lu, b1))
    x,info = getrs(lu, piv, b1, trans=trans, overwrite_b=overwrite_b)
    if info == 0:
        return x
    raise ValueError('illegal value in %d-th argument of internal gesv|posv'
                                                                    % -info)


def lu(a, permute_l=False, overwrite_a=False):
    """Compute pivoted LU decompostion of a matrix.

    The decomposition is::

        A = P L U

    where P is a permutation matrix, L lower triangular with unit
    diagonal elements, and U upper triangular.

    Parameters
    ----------
    a : array, shape (M, N)
        Array to decompose
    permute_l : boolean
        Perform the multiplication P*L  (Default: do not permute)
    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)

    Returns
    -------
    (If permute_l == False)
    p : array, shape (M, M)
        Permutation matrix
    l : array, shape (M, K)
        Lower triangular or trapezoidal matrix with unit diagonal.
        K = min(M, N)
    u : array, shape (K, N)
        Upper triangular or trapezoidal matrix

    (If permute_l == True)
    pl : array, shape (M, K)
        Permuted L matrix.
        K = min(M, N)
    u : array, shape (K, N)
        Upper triangular or trapezoidal matrix

    Notes
    -----
    This is a LU factorization routine written for Scipy.

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError('expected matrix')
    overwrite_a = overwrite_a or (_datanotshared(a1, a))
    flu, = get_flinalg_funcs(('lu',), (a1,))
    p, l, u, info = flu(a1, permute_l=permute_l, overwrite_a=overwrite_a)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of '
                                            'internal lu.getrf' % -info)
    if permute_l:
        return l, u
    return p, l, u
