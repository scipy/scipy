"""QR decomposition functions."""

import numpy

# Local imports
from blas import get_blas_funcs
from lapack import get_lapack_funcs, find_best_lapack_type
from misc import _datacopied

# XXX: what is qr_old, should it be kept?
__all__ = ['qr', 'rq', 'qr_old']

def safecall(f, name, *args, **kwargs):
    lwork = kwargs.pop("lwork", None)
    if lwork is None:
        ret = f(*args, lwork=-1, **kwargs)
        lwork = ret[-2][0].real.astype(numpy.int)
    ret = f(*args, lwork=lwork, **kwargs)
    if ret[-1] < 0:
        raise ValueError("illegal value in %d-th argument of internal %s"
                         % (-ret[-1], name))
    return ret[:-2]

def qr(a, overwrite_a=False, lwork=None, mode='full', pivoting=False, c=None,
    overwrite_c=False):
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
    mode : {'full', 'r', 'economic', 'reflectors'} if c is None,
        or {'left', 'right'} if c is given.
        Determines what information is to be returned: either both Q and R
        ('full', default), only R ('r'), both Q and R but computed in
        economy-size ('economic', see Notes), or the elementary
        reflectors of Q, and R ('reflectors', see Notes).
        In case c is given, dot(Q, c) is returned if mode is left,
        dot(c, Q) is returned if mode is right.
        
    pivoting : bool, optional
        Whether or not factorization should include pivoting for rank-revealing
        qr decomposition. If pivoting, compute the decomposition
        :lm:`A P = Q R` as above, but where P is chosen such that the diagonal
        of R is non-increasing.
    c : array, two-dimensional
        calculate the product of c and q, depending on the mode:
        'left': dot(Q, c)
        'right': dot(c, Q)
        the shape of c must be appropriate for the matrix multiplications
    overwrite_c: bool, optional
        Whether data in c is overwritten (may improve performance)
        

    Returns
    -------
    Q : double or complex ndarray
        Of shape (M, M), or (M, K) for ``mode='economic'``.  Not returned if
        ``mode='r'``.
    CQ : double or complex ndarray
        the product of Q and c, as defined in mode
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

    If ``mode=reflectors``, the matrix Q is represented as a product of
    elementary reflectors

      Q = H[1] H[2] . . . H[k], where k = min(m,n).

    Each H[i] has the form

     H[i] = I - (1 / v[i]) * v * v.T

    where v is a vector with v[:i] = 0. v is returned as Q[:, i]

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
    if (not mode in ['full', 'qr', 'r', 'economic', 'reflectors']
                 ) and c is None:
        raise ValueError(
                 "Mode argument should be one of "
                 "['full', 'r', 'economic', 'reflectors']")
    if (not mode in ['left', 'right']
                 ) and c is not None:
        raise ValueError(
                 "Mode argument should be one of ['left', 'right']")

    a1 = numpy.asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError("expected 2D array")
    M, N = a1.shape
    overwrite_a = overwrite_a or (_datacopied(a1, a))

    if pivoting:
        geqp3, = get_lapack_funcs(('geqp3',), (a1,))
        qr, jpvt, tau = safecall(geqp3, "geqp3", a1, overwrite_a=overwrite_a)
        jpvt -= 1 # geqp3 returns a 1-based index array, so subtract 1
    else:
        geqrf, = get_lapack_funcs(('geqrf',), (a1,))
        qr, tau = safecall(geqrf, "geqrf", a1, lwork=lwork,
            overwrite_a=overwrite_a)

    if not mode == 'economic' or M < N:
        R = numpy.triu(qr)
    else:
        R = numpy.triu(qr[:N, :N])

    if mode == 'r':
        if pivoting:
            return R, jpvt
        else:
            return R

    if mode == 'reflectors':
        K = min(M, N)
        Q = numpy.tril(qr[:, :K], -1)
        numpy.fill_diagonal(Q[:K, :K], 1)
        Q *= tau
        if pivoting:
            return Q, R, jpvt
        else:
            return Q, R

    if c is None:
        if find_best_lapack_type((a1,))[0] in ('s', 'd'):
            gor_un_gqr, = get_lapack_funcs(('orgqr',), (qr,))
        else:
            gor_un_gqr, = get_lapack_funcs(('ungqr',), (qr,))

        if M < N:
            Q, = safecall(gor_un_gqr, "gorgqr/gungqr", qr[:,0:M], tau,
                lwork=lwork, overwrite_a=1)
        elif mode == 'economic':
            Q, = safecall(gor_un_gqr, "gorgqr/gungqr", qr, tau, lwork=lwork,
                overwrite_a=1)
        else:
            t = qr.dtype.char
            qqr = numpy.empty((M, M), dtype=t)
            qqr[:, :N] = qr
            Q, = safecall(gor_un_gqr, "gorgqr/gungqr", qqr, tau, lwork=lwork,
                overwrite_a=1)
    else:
        c = numpy.asarray_chkfinite(c)
        if c.ndim == 1:
            if mode == "left":
                c = c[:, numpy.newaxis]
            else:
                c = c[numpy.newaxis, :]
        if find_best_lapack_type((a1,))[0] in ('s', 'd'):
            gor_un_mqr, = get_lapack_funcs(('ormqr',), (qr,))
            trans = "T"
        else:
            gor_un_mqr, = get_lapack_funcs(('unmqr',), (qr,))
            trans = "C"

        qr = qr[:, :min(M, N)]
        if M > N and mode == "left":
            if overwrite_c:
                cc = c.T.conjugate()
            else:
                cc = numpy.zeros((c.shape[1], max(M, N)),
                    dtype=c.dtype, order="F")
                cc[:, :c.shape[0]] = c.T.conjugate()
                overwrite_c = True
            lr = "R"
        elif c.flags["C_CONTIGUOUS"] and (
                mode == "left" or M <= N) and trans == "T":
            cc = c.T
            lr = "R" if mode == "left" else "L"
        else: 
            trans = "N"
            cc = c
            lr = "L" if mode == "left" else "R"
        Q, = safecall(gor_un_mqr, "gormqr/gunmqr", lr, trans, qr, tau, cc,
            lwork=lwork, overwrite_c=overwrite_c)
        if trans != "N":
            Q = Q.T.conjugate()
        if mode == "right":
            Q = Q[:, :min(M, N)]

    if pivoting:
        return Q, R, jpvt
    return Q, R



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
    a1 = numpy.asarray_chkfinite(a)
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
    R = numpy.triu(qr)
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

    a1 = numpy.asarray_chkfinite(a)
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
        R = numpy.triu(rq, N-M)
    else:
        R = numpy.triu(rq[-M:, -M:])

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
