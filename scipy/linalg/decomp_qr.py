"""QR decomposition functions."""

import numpy
from numpy import asarray_chkfinite

# Local imports
import special_matrices
from blas import get_blas_funcs
from lapack import get_lapack_funcs, find_best_lapack_type
from misc import _datacopied

# XXX: what is qr_old, should it be kept?
__all__ = ['qr', 'rq', 'qr_multiply', 'qr_old']

def _raw_qr(a1, overwrite_a=False, lwork=None, mode='full', pivoting=False):
    """encapsulates common blocks from qr and qr_multiply. Returns the qr
    matrix as a product of elementary reflectors that can be used by both
    {or,un}gqr and {or,un}mqr.
    """
    jpvt = None
    M, N = a1.shape

    if pivoting:
        geqp3, = get_lapack_funcs(('geqp3',), (a1,))
        if lwork is None or lwork == -1:
            # get optimal work array
            qr, jpvt, tau, work, info = geqp3(a1, lwork=-1, overwrite_a=1)
            lwork = work[0].real.astype(numpy.int)

        qr, jpvt, tau, work, info = geqp3(a1, lwork=lwork,
            overwrite_a=overwrite_a)
        jpvt -= 1 # geqp3 returns a 1-based index array, so subtract 1
        if info < 0:
            raise ValueError("illegal value in %d-th argument of internal geqp3"
                                                                        % -info)
    else:
        geqrf, = get_lapack_funcs(('geqrf',), (a1,))
        if lwork is None or lwork == -1:
            # get optimal work array
            qr, tau, work, info = geqrf(a1, lwork=-1, overwrite_a=1)
            lwork = work[0].real.astype(numpy.int)

        qr, tau, work, info = geqrf(a1, lwork=lwork, overwrite_a=overwrite_a)
        if info < 0:
            raise ValueError("illegal value in %d-th argument of internal geqrf"
                                                                        % -info)

    if mode != 'economic' or M < N:
        R = special_matrices.triu(qr)
    else:
        R = special_matrices.triu(qr[0:N, 0:N])

    if mode == 'r':
        pass
    elif M <= N:
        qr = qr[:, :M]
    elif mode == 'full':
        tmp = numpy.empty((M, M), dtype=qr.dtype)
        tmp[:,0:N] = qr
        qr = tmp

    return qr, R, jpvt, tau, work, info


def qr(a, overwrite_a=False, lwork=None, mode='full', pivoting=False):
    """Compute QR decomposition of a matrix.

    Calculate the decomposition :lm:`A = Q R` where Q is unitary/orthogonal
    and R upper triangular.

    Parameters
    ----------
    a : array, shape (M, N)
        Matrix to be decomposed
    overwrite_a : bool, optional
        Whether data in a is overwritten (may improve performance)
    lwork : int, optional
        Work array size, lwork >= a.shape[1]. If None or -1, an optimal size
        is computed.
    mode : {'full', 'r', 'economic'}
        Determines what information is to be returned: either both Q and R
        ('full', default), only R ('r') or both Q and R but computed in
        economy-size ('economic', see Notes).
    pivoting : bool, optional
        Whether or not factorization should include pivoting for rank-revealing
        qr decomposition. If pivoting, compute the decomposition
        :lm:`A P = Q R` as above, but where P is chosen such that the diagonal
        of R is non-increasing.

    Returns
    -------
    Q : double or complex ndarray
        Of shape (M, M), or (M, K) for ``mode='economic'``.  Not returned if
        ``mode='r'``.
    R : double or complex ndarray
        Of shape (M, N), or (K, N) for ``mode='economic'``.  ``K = min(M, N)``.
    P : integer ndarray
        Of shape (N,) for ``pivoting=True``. Not returned if ``pivoting=False``.

    Raises
    ------
    LinAlgError
        Raised if decomposition fails

    Notes
    -----
    This is an interface to the LAPACK routines dgeqrf, zgeqrf,
    dorgqr, zungqr, dgeqp3, and zgeqp3.

    If ``mode=economic``, the shapes of Q and R are (M, K) and (K, N) instead
    of (M,M) and (M,N), with ``K=min(M,N)``.

    Examples
    --------
    >>> from scipy import random, linalg, dot, diag, all, allclose
    >>> a = random.randn(9, 6)

    >>> q, r = linalg.qr(a)
    >>> allclose(a, dot(q, r))
    True
    >>> q.shape, r.shape
    ((9, 9), (9, 6))

    >>> r2 = linalg.qr(a, mode='r')
    >>> allclose(r, r2)
    True

    >>> q3, r3 = linalg.qr(a, mode='economic')
    >>> q3.shape, r3.shape
    ((9, 6), (6, 6))

    >>> q4, r4, p4 = linalg.qr(a, pivoting=True)
    >>> d = abs(diag(r4))
    >>> all(d[1:] <= d[:-1])
    True
    >>> allclose(a[:, p4], dot(q4, r4))
    True
    >>> q4.shape, r4.shape, p4.shape
    ((9, 9), (9, 6), (6,))

    >>> q5, r5, p5 = linalg.qr(a, mode='economic', pivoting=True)
    >>> q5.shape, r5.shape, p5.shape
    ((9, 6), (6, 6), (6,))

    """
    if mode == 'qr':
        # 'qr' was the old default, equivalent to 'full'. Neither 'full' nor
        # 'qr' are used below, but set to 'full' anyway to be sure
        mode = 'full'
    if not mode in ['full', 'qr', 'r', 'economic']:
        raise ValueError(\
                 "Mode argument should be one of ['full', 'r', 'economic']")

    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError("expected 2D array")
    M, N = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))

    qr, R, jpvt, tau, work, info = _raw_qr(
        a1, overwrite_a=overwrite_a, lwork=lwork, mode=mode, pivoting=pivoting)

    if mode == 'r':
        if pivoting:
            return R, jpvt
        else:
            return R

    # LAPACK requires to generate Q matrix from elementary reflectors
    if find_best_lapack_type((a1,))[0] in ('s', 'd'):
        gor_un_gqr, = get_lapack_funcs(('orgqr',), (qr,))
    else:
        gor_un_gqr, = get_lapack_funcs(('ungqr',), (qr,))

    Q, work, info = gor_un_gqr(qr, tau, lwork=-1, overwrite_a=1)
    lwork = work[0].real.astype(numpy.int)
    Q, work, info = gor_un_gqr(qr, tau, lwork=lwork, overwrite_a=1)

    if info < 0:
        raise ValueError("illegal value in %d-th argument of internal gorgqr"
                                                                    % -info)
    if pivoting:
        return Q, R, jpvt
    return Q, R


def qr_multiply(a, *args, **kwargs):
    """Perform dot(q, b), dot(q, c) ... where q is given by the QR-decompostion of a.

    Parameters
    ----------
    a : array, shape (M, N)
       array to perform the QR decomposition.

    b, c, d ... : array, up to 2D
       arrays to be multiplied by Q.

    mode : {'full', 'economic'}, optional
        Wether to compute Q and R in full (mode='full') or economy-size
        (mode='economy'). If none given, 'full' will be used.

    trans : boolean, optional.
        Type of opeartion to perform. If a is a complex array, Q^H is used
        instead of Q^T. If none given, False will be used

        ========  ===========
        trans     opeartion
        ========  ===========
        False     dot(Q, A)
        True      dot(Q^T, A)
        ========  ===========

    Returns
    -------
    qb, qc, .. : sequence of arrays
        Arrays dot(q, b), dot(q, c), ...

    R : array, shape (M, N)
        Triangular matrix from the qr decompostion of a.

    P : integer array, shape (N,)
        Pivoting matrix. Onlye returned if ``pivoting=True``.

    Examples
    --------
    q_multiply can be used to compute the least squares solution ax = b:

    >>> from scipy import linalg
    >>> a = [[1, 1], [2, 1], [2, 2]]
    >>> b = [0, 1, 2]
    >>> qb, r = linalg.qr_multiply(a, b, trans=True, mode='economic')
    >>> print linalg.solve_triangular(r, qb)
    [[ 0.2],
     [ 0.6]])
    """
    trans = kwargs.pop('trans', False)
    overwrite_b = kwargs.pop('overwrite_b', False)
    a1 = asarray_chkfinite(a)
    result = []
    if trans:
        N, M = a1.shape
    else:
        M, N = a1.shape

    q, R, jpvt, tau, work, info = _raw_qr(a1, **kwargs)

    for b in args:
        b = asarray_chkfinite(b)
        if b.ndim == 1:
            b = b[:, numpy.newaxis]

        if b.dtype.char in numpy.typecodes['Complex'] or \
           q.dtype.char in numpy.typecodes['Complex']:
            ormqr, = get_lapack_funcs(('unmqr',), (q,b))
        else:
            ormqr, = get_lapack_funcs(('ormqr',), (q,b))

        if (trans == False and q.shape[1] != b.shape[0]) or \
           (trans == True and q.shape[0] != b.shape[0]):
            raise ValueError("Incompatible dimensions in qr_multiply")

        if b.shape[0] < M:
            tmp = numpy.empty((M, b.shape[1]), dtype = q.dtype)
            tmp[:b.shape[0]] = b
            b = tmp
            overwrite_b = True

        c, work, info = ormqr(q, tau, b, lwork=-1)
        lwork = work[0].real.astype(numpy.int32)
        c, work, info = ormqr(q, tau, b, lwork=lwork, trans=trans,
                              overwrite_b=overwrite_b)
        if info < 0:
            raise ValueError("illegal value in %d-th argument of internal gorgqr"
                                                                        % -info)
        result.append(c[:M])

    if jpvt is not None:
        result = result + [R, jpvt]
    else:
        result = result + [R]
    return result


def qr_old(a, overwrite_a=False, lwork=None):
    """Compute QR decomposition of a matrix.

    Calculate the decomposition :lm:`A = Q R` where Q is unitary/orthogonal
    and R upper triangular.

    Parameters
    ----------
    a : array, shape (M, N)
        Matrix to be decomposed
    overwrite_a : boolean
        Whether data in a is overwritten (may improve performance)
    lwork : integer
        Work array size, lwork >= a.shape[1]. If None or -1, an optimal size
        is computed.

    Returns
    -------
    Q : double or complex array, shape (M, M)
    R : double or complex array, shape (M, N)
        Size K = min(M, N)

    Raises LinAlgError if decomposition fails

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError('expected matrix')
    M,N = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))
    geqrf, = get_lapack_funcs(('geqrf',), (a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        qr, tau, work, info = geqrf(a1, lwork=-1, overwrite_a=1)
        lwork = work[0]
    qr, tau, work, info = geqrf(a1, lwork=lwork, overwrite_a=overwrite_a)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal geqrf'
                                                                    % -info)
    gemm, = get_blas_funcs(('gemm',), (qr,))
    t = qr.dtype.char
    R = special_matrices.triu(qr)
    Q = numpy.identity(M, dtype=t)
    ident = numpy.identity(M, dtype=t)
    zeros = numpy.zeros
    for i in range(min(M, N)):
        v = zeros((M,), t)
        v[i] = 1
        v[i+1:M] = qr[i+1:M, i]
        H = gemm(-tau[i], v, v, 1+0j, ident, trans_b=2)
        Q = gemm(1, Q, H)
    return Q, R


def rq(a, overwrite_a=False, lwork=None, mode='full'):
    """Compute RQ decomposition of a square real matrix.

    Calculate the decomposition :lm:`A = R Q` where Q is unitary/orthogonal
    and R upper triangular.

    Parameters
    ----------
    a : array, shape (M, M)
        Matrix to be decomposed
    overwrite_a : boolean
        Whether data in a is overwritten (may improve performance)
    lwork : integer
        Work array size, lwork >= a.shape[1]. If None or -1, an optimal size
        is computed.
    mode : {'full', 'r', 'economic'}
        Determines what information is to be returned: either both Q and R
        ('full', default), only R ('r') or both Q and R but computed in
        economy-size ('economic', see Notes).

    Returns
    -------
    R : double array, shape (M, N)
    Q : double or complex array, shape (M, M)

    Raises LinAlgError if decomposition fails

    Examples
    --------
    >>> from scipy import linalg
    >>> from numpy import random, dot, allclose
    >>> a = random.randn(6, 9)
    >>> r, q = linalg.rq(a)
    >>> allclose(a, dot(r, q))
    True
    >>> r.shape, q.shape
    ((6, 9), (9, 9))
    >>> r2 = linalg.rq(a, mode='r')
    >>> allclose(r, r2)
    True
    >>> r3, q3 = linalg.rq(a, mode='economic')
    >>> r3.shape, q3.shape
    ((6, 6), (6, 9))

    """
    if not mode in ['full', 'r', 'economic']:
        raise ValueError(\
                 "Mode argument should be one of ['full', 'r', 'economic']")

    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError('expected matrix')
    M, N = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))

    gerqf, = get_lapack_funcs(('gerqf',), (a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        rq, tau, work, info = gerqf(a1, lwork=-1, overwrite_a=1)
        lwork = work[0].real.astype(numpy.int)
    rq, tau, work, info = gerqf(a1, lwork=lwork, overwrite_a=overwrite_a)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal gerqf'
                                                                    % -info)
    if not mode == 'economic' or N < M:
        R = special_matrices.triu(rq, N-M)
    else:
        R = special_matrices.triu(rq[-M:, -M:])

    if mode == 'r':
        return R

    if find_best_lapack_type((a1,))[0] in ('s', 'd'):
        gor_un_grq, = get_lapack_funcs(('orgrq',), (rq,))
    else:
        gor_un_grq, = get_lapack_funcs(('ungrq',), (rq,))

    if N < M:
        # get optimal work array
        Q, work, info = gor_un_grq(rq[-N:], tau, lwork=-1, overwrite_a=1)
        lwork = work[0].real.astype(numpy.int)
        Q, work, info = gor_un_grq(rq[-N:], tau, lwork=lwork, overwrite_a=1)
    elif mode == 'economic':
        # get optimal work array
        Q, work, info = gor_un_grq(rq, tau, lwork=-1, overwrite_a=1)
        lwork = work[0].real.astype(numpy.int)
        Q, work, info = gor_un_grq(rq, tau, lwork=lwork, overwrite_a=1)
    else:
        rq1 = numpy.empty((N, N), dtype=rq.dtype)
        rq1[-M:] = rq
        # get optimal work array
        Q, work, info = gor_un_grq(rq1, tau, lwork=-1, overwrite_a=1)
        lwork = work[0].real.astype(numpy.int)
        Q, work, info = gor_un_grq(rq1, tau, lwork=lwork, overwrite_a=1)

    if info < 0:
        raise ValueError("illegal value in %d-th argument of internal orgrq"
                                                                    % -info)
    return R, Q
