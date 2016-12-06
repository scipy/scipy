#
# Author: Ilhan Polat, December, 2016
#

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy import atleast_2d
from .decomp import _asarray_validated
from .misc import LinAlgError
from .lapack import get_lapack_funcs

__all__ = ['ldl']


def ldl(A, lower=True, rook_pivoting=False, only_sym=False,
        overwrite_a=False, check_finite=True):
    """ Computes the ``ldl`` or Bunch-Kaufman factorization of a symmetric/
    hermitian matrix.

    This function returns a block diagonal matrix D consisting blocks of size
    at most 2x2 and also a possibly permuted unit lower triangular matrix
    ``L`` such that the factorization ``A = L D L^H`` holds.

    Depending on the value of the boolean ``lower``, only upper or lower
    triangular part of the input array is referenced. Hence a triangular
    part of matrix on entry would give the same result as if the full matrix
    is supplied.

    The permutation array can be used to triangulize the outer factors simply
    by a row shuffle, i.e., ``lu[perm, :]`` is a upper/lower triangular
    matrix. This is also equivalent to seeing the shuffle as a permutation
    matrix ``P.dot(lu)`` where ``P`` is an identity matrix permuted by
    ``perm``.

    Parameters
    ----------
    A : array_like
        Symmetric/hermitian input matrix
    lower : bool, optional
        This switches between the lower and upper triangular outer factors of
        the factorization. As the name implies lower triangular that is
        (``lower=True``) is the default.
    rook_pivoting : bool, optional
        Switches between the pivoting algorithms. See Notes section below.
    only_sym : bool, optional
        In case the input matrix is complex valued and symmetric, by setting
        this, the matrix is assumed to be only symmetric and not hermitian.
    overwrite_a : bool, optional
        Allow overwriting data in `a` (may enhance performance). The default
        is False.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    lu : ndarray
        The permuted square outer factor of the factorization.
    d : ndarray
        The block diagonal multiplier of the factorization.
    perm : ndarray
        The permutation indices array that brings lu into triangular form.

    Notes
    -----
    The main computational routines are LAPACK's ?SYTRF and, if an Hermitian
    matrix is factorized, ?HETRF routines.

    If ``only_sym`` is set instead of assuming the complex-valued array to
    be hermitian and choosing the (C/Z)HETRF solvers, symmetric solvers
    (C/Z)SYTRF are selected.

    If ``rook_pivoting`` keyword is set to True the so-called bounded
    Bunch-Kaufman algorithm given in [2]_ is used. Otherwise, the standard
    algorithm given in [1]_ is used. Internally, this selects between
    ``?(HE/SY)TRF`` and ``?(HE/SY)TRF_ROOK`` solvers.

    See also
    --------
    cholesky, lu

    .. versionadded:: 0.19.0

    References
    ----------
    .. [1] J.R. Bunch, L. Kaufman, Some stable methods for calculating
       inertia and solving symmetric linear systems, Math. Comput. Vol.31,
       1977. DOI: 10.2307/2005787

    .. [2]  S. H. Cheng and N. J. Higham, A modified cholesky algorithm
       based on a symmetric indefinite factorization, SIAM J. Matrix Anal.
       Appl., Vol.19, 1998, DOI: 10.1137/S0895479896302898
    """
    a = atleast_2d(_asarray_validated(A, check_finite=check_finite))
    if a.shape[0] != a.shape[1]:
        raise ValueError('The input matrix should be square.')
    n = a.shape[0]
    r_or_c = complex if np.iscomplexobj(a) else float

    # Get the low-level functions
    if r_or_c is complex and not only_sym:
        if np.any(np.imag(np.diag(a))):
            raise ValueError('The matrix diagonal entries should be real '
                             'for the hermitian decomposition. Otherwise '
                             'set only_sym keyword to True for symmetric '
                             'solution.')

        s, sl = ('hetrf_rook', 'hetrf_rook_lwork') if rook_pivoting else (
                                                    'hetrf', 'hetrf_lwork')
    else:
        s, sl = ('sytrf_rook', 'sytrf_rook_lwork') if rook_pivoting else (
                                                    'sytrf', 'sytrf_lwork')

    solver, solver_lwork = get_lapack_funcs((s, sl), (a,))
    lwork, linfo = solver_lwork(n, lower=lower)
    if linfo < 0:
        raise ValueError('{} exited with the internal error '
                         '"illegal value in argument number {}." while '
                         'requesting the optimal block size LWORK. See '
                         'LAPACK documentation for the error codes.'
                         ''.format(s, -linfo)
                         )

    ldu, piv, sinfo = solver(a, lwork=lwork, lower=lower)
    if sinfo < 0:
        raise ValueError('{} exited with the internal error '
                         '"illegal value in argument number {}". See '
                         'LAPACK documentation for the error codes.'
                         ''.format(s, -sinfo)
                         )

    swap_arr, pivot_arr = _ldl_sanitize_ipiv(piv, lower=lower,
                                             rook=rook_pivoting)
    d, l = _ldl_get_d_and_l(ldu, pivot_arr, lower=lower, only_sym=only_sym)
    lu, perm = _ldl_construct_tri_factor(l, swap_arr, pivot_arr, lower=lower)

    return lu, d, perm


def _ldl_sanitize_ipiv(a, lower=True, rook=False):
    """
    This helper function takes the strangely encoded permutation array
    returned by the LAPACK routines ?(HE/SY)TRF, ?(HE/SY)TRF_ROOK and
    converts it into regularized permutation and diagonal pivot size format.

    Since FORTRAN uses 1-indexing and LAPACK uses different start points for
    upper and lower formats there are certain offsets in the indices used
    below.

    Let's assume a result where the matrix is 6x6 and there are two 2x2
    and two 1x1 blocks reported by the routine. To ease the coding efforts,
    we still populate a 6-sized array and fill zeros as the following ::

        pivots = [2, 0, 2, 0, 1, 1]

    This denotes a diagonal matrix of the form ::

        [x x        ]
        [x x        ]
        [    x x    ]
        [    x x    ]
        [        x  ]
        [          x]

    In other words, we write 2 when the 2x2 block is first encountered and
    automatically write 0 to the next entry and skip the next spin of the
    loop. Thus, a separate counter or array appends to keep track of block
    sizes are avoided. Zeros can be removed much easier instead.

    Parameters
    ----------
    a : ndarray
        The permutation array ipiv returned by LAPACK
    lower : bool, optional
        The switch to select whether upper or lower triangle is chosen in
        the LAPACK call.
    rook : bool, optional
        The switch to select whether the partial or the rook pivoting method
        is chosen in the LAPACK call.
    Returns
    -------
    swap_ : ndarray
        The array that defines the row/column swap operations. For example,
        if row two is swapped with four the result is [0, 3, 2, 3].

    """
    n = a.size
    swap_ = np.arange(n)
    pivots = np.zeros_like(swap_, dtype=int)
    skip_2x2 = False

    # Some upper/lower dependent offset values
    x, y, rs, re, ri = (1, 0, 0, n, 1) if lower else (-1, -1, n-1, -1, -1)

    if rook:
        for ind in range(rs, re, ri):
            if skip_2x2:
                skip_2x2 = False
                continue

            cur_val = a[ind]
            # do have a 1x1 block or not?
            if cur_val > 0:
                if cur_val != ind + 1:
                    swap_[ind] = swap_[cur_val - 1]
                pivots[ind] = 1
            # not.
            elif cur_val < 0 and a[ind + x] < 0:
                # first neg entry of 2x2 block identifier
                if -cur_val != ind + 1:
                    swap_[ind] = swap_[-cur_val - 1]
                if -a[ind + x] != ind + x + 1:
                    swap_[ind + x] = swap_[-a[ind + x] - 1]
                pivots[ind + y] = 2
                skip_2x2 = True
            else:  # Doesn't make sense, give up
                raise LinAlgError('While parsing the permutation array '
                                  'in "scipy.linalg.ldl", invalid entries'
                                  ' found.')
    # rook = False
    else:
        for ind in range(rs, re, ri):
            if skip_2x2:
                skip_2x2 = False
                continue

            cur_val = a[ind]
            # do we have a 1x1 block or not?
            if cur_val > 0:
                if cur_val != ind+1:
                    # Index value != array value --> permutation required
                    swap_[ind] = swap_[cur_val-1]
                pivots[ind] = 1
            # Not.
            elif cur_val < 0 and cur_val == a[ind+x]:
                # first neg entry of 2x2 block identifier
                if -cur_val != ind+2:
                    # Index value != array value --> permutation required
                    swap_[ind+x] = swap_[-cur_val-1]
                pivots[ind+y] = 2
                skip_2x2 = True
            else:  # Doesn't make sense, give up
                raise LinAlgError('While parsing the permutation array '
                                  'in "scipy.linalg.ldl", invalid entries'
                                  ' found. Either "rook" key is not set'
                                  ' properly or the array is invalid.')
    return swap_, pivots


def _ldl_get_d_and_l(ldu, pivs, lower=True, only_sym=False):
    """
    Helper function to extract the diagonal and triangular matrices for
    LDL.T factorization.

    Parameters
    ----------
    ldu : ndarray
        The compact output returned by the LAPACK routing
    pivs : ndarray
        The sanitized array of {0, 1, 2} denoting the sizes of the pivots. For
        every 2 there is a succeeding 0.
    lower : bool, optional
        If set to False, upper triangular part is considered.

    Returns
    -------
    d : ndarray
        The block diagonal matrix.
    lu : ndarray
        The upper/lower triangular matrix
    """
    is_c = np.iscomplexobj(ldu)
    d = np.diag(np.diag(ldu))
    n = d.shape[0]
    blk_i = 0  # block index

    # row/column offsets for selecting sub-, super-diagonal
    x, y = (1, 0) if lower else (0, 1)

    lu = np.tril(ldu, -1) if lower else np.triu(ldu, 1)
    diag_inds = np.arange(n)
    lu[diag_inds, diag_inds] = 1

    for blk in pivs[pivs != 0]:
        # increment the block index and check for 2s
        # if 2 then copy the off diagonals depending on uplo
        inc = blk_i + blk

        if blk == 2:
            # If Hermitian matrix is factorized, the offdiagonal
            # should be conjugated.
            if is_c and not only_sym:
                d[[blk_i+x, blk_i+y], [blk_i+y, blk_i+x]] = [
                    ldu[blk_i+x, blk_i+y], ldu[blk_i+x, blk_i+y].conj()]
            else:
                d[[blk_i+x, blk_i+y], [blk_i+y, blk_i+x]] = \
                                                        ldu[blk_i+x, blk_i+y]
            lu[blk_i+x, blk_i+y] = 0
        blk_i = inc

    return d, lu


def _ldl_construct_tri_factor(lu, swap_vec, pivs, lower=True):
    """
    Helper function to construct explicit outer factors of LDL factorization.

    If lower is True the permuted factors are multiplied as L(1)*L(2)*...*L(k).
    Otherwise, the permuted factors are multiplied as L(k)*...*L(2)*L(1). See
    LAPACK documentation for more details.

    Parameters
    ----------
    lu : ndarray
        The triangular array that is extracted from LAPACK routine call with
        ones on the diagonals.
    swap_vec : ndarray
        The array that defines the row swapping indices. If k'th entry is m
        then rows k,m are swapped. Notice that m'th entry is not necessarily
        k to avoid undoing the swapping.
    lower : bool, optional
        The boolean to switch between lower and upper triangular structure.

    Returns
    -------
    lu : ndarray
        The square outer factor which satisfies the L * D * L.T = A
    perm : ndarray
        The permutation vector that brings the lu to the triangular form

    Notes
    -----
    It should be noted that we overwrite the original argument "lu".
    """
    n = lu.shape[0]
    perm = np.arange(n)
    # Setup the reading order of the permutation matrix for upper/lower
    rs, re, ri = (n-1, -1, -1) if lower else (0, n, 1)

    for ind in range(rs, re, ri):
        s_ind = swap_vec[ind]
        if s_ind != ind:
            # Column start and end positions
            col_s = ind if lower else 0
            col_e = n if lower else ind+1

            # If we stumble upon a 2x2 block include both cols in the perm.
            if pivs[ind] == (0 if lower else 2):
                col_s += -1 if lower else 0
                col_e += 0 if lower else 1
            lu[[s_ind, ind], col_s:col_e] = lu[[ind, s_ind], col_s:col_e]
            perm[[s_ind, ind]] = perm[[ind, s_ind]]

    return lu, np.argsort(perm)
