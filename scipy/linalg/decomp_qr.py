"""QR decomposition functions."""

from warnings import warn

import numpy
from numpy import asarray_chkfinite, complexfloating

# Local imports
import special_matrices
from blas import get_blas_funcs
from lapack import get_lapack_funcs, find_best_lapack_type
from misc import _datanotshared


def qr(a, overwrite_a=False, lwork=None, econ=None, mode='qr'):
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
    econ : boolean
        Whether to compute the economy-size QR decomposition, making shapes
        of Q and R (M, K) and (K, N) instead of (M,M) and (M,N). K=min(M,N).
        Default is False.
    mode : {'qr', 'r'}
        Determines what information is to be returned: either both Q and R
        or only R.

    Returns
    -------
    (if mode == 'qr')
    Q : double or complex array, shape (M, M) or (M, K) for econ==True

    (for any mode)
    R : double or complex array, shape (M, N) or (K, N) for econ==True
        Size K = min(M, N)

    Raises LinAlgError if decomposition fails

    Notes
    -----
    This is an interface to the LAPACK routines dgeqrf, zgeqrf,
    dorgqr, and zungqr.

    Examples
    --------
    >>> from scipy import random, linalg, dot
    >>> a = random.randn(9, 6)
    >>> q, r = linalg.qr(a)
    >>> allclose(a, dot(q, r))
    True
    >>> q.shape, r.shape
    ((9, 9), (9, 6))

    >>> r2 = linalg.qr(a, mode='r')
    >>> allclose(r, r2)

    >>> q3, r3 = linalg.qr(a, econ=True)
    >>> q3.shape, r3.shape
    ((9, 6), (6, 6))

    """
    if econ is None:
        econ = False
    else:
        warn("qr econ argument will be removed after scipy 0.7. "
             "The economy transform will then be available through "
             "the mode='economic' argument.", DeprecationWarning)

    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError("expected 2D array")
    M, N = a1.shape
    overwrite_a = overwrite_a or (_datanotshared(a1, a))

    geqrf, = get_lapack_funcs(('geqrf',), (a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        qr, tau, work, info = geqrf(a1, lwork=-1, overwrite_a=1)
        lwork = work[0]

    qr, tau, work, info = geqrf(a1, lwork=lwork, overwrite_a=overwrite_a)
    if info < 0:
        raise ValueError("illegal value in %d-th argument of internal geqrf"
                                                                    % -info)
    if not econ or M < N:
        R = special_matrices.triu(qr)
    else:
        R = special_matrices.triu(qr[0:N, 0:N])

    if mode == 'r':
        return R

    if find_best_lapack_type((a1,))[0] == 's' or \
                find_best_lapack_type((a1,))[0] == 'd':
        gor_un_gqr, = get_lapack_funcs(('orgqr',), (qr,))
    else:
        gor_un_gqr, = get_lapack_funcs(('ungqr',), (qr,))

    if M < N:
        # get optimal work array
        Q, work, info = gor_un_gqr(qr[:,0:M], tau, lwork=-1, overwrite_a=1)
        lwork = work[0]
        Q, work, info = gor_un_gqr(qr[:,0:M], tau, lwork=lwork, overwrite_a=1)
    elif econ:
        # get optimal work array
        Q, work, info = gor_un_gqr(qr, tau, lwork=-1, overwrite_a=1)
        lwork = work[0]
        Q, work, info = gor_un_gqr(qr, tau, lwork=lwork, overwrite_a=1)
    else:
        t = qr.dtype.char
        qqr = numpy.empty((M, M), dtype=t)
        qqr[:,0:N] = qr
        # get optimal work array
        Q, work, info = gor_un_gqr(qqr, tau, lwork=-1, overwrite_a=1)
        lwork = work[0]
        Q, work, info = gor_un_gqr(qqr, tau, lwork=lwork, overwrite_a=1)

    if info < 0:
        raise ValueError("illegal value in %d-th argument of internal gorgqr"
                                                                    % -info)
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
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError, 'expected matrix'
    M,N = a1.shape
    overwrite_a = overwrite_a or (_datanotshared(a1, a))
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


def rq(a, overwrite_a=False, lwork=None):
    """Compute RQ decomposition of a square real matrix.

    Calculate the decomposition :lm:`A = R Q` where Q is unitary/orthogonal
    and R upper triangular.

    Parameters
    ----------
    a : array, shape (M, M)
        Square real matrix to be decomposed
    overwrite_a : boolean
        Whether data in a is overwritten (may improve performance)
    lwork : integer
        Work array size, lwork >= a.shape[1]. If None or -1, an optimal size
        is computed.
    econ : boolean

    Returns
    -------
    R : double array, shape (M, N) or (K, N) for econ==True
        Size K = min(M, N)
    Q : double or complex array, shape (M, M) or (M, K) for econ==True

    Raises LinAlgError if decomposition fails

    """
    # TODO: implement support for non-square and complex arrays
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError('expected matrix')
    M,N = a1.shape
    if M != N:
        raise ValueError('expected square matrix')
    if issubclass(a1.dtype.type, complexfloating):
        raise ValueError('expected real (non-complex) matrix')
    overwrite_a = overwrite_a or (_datanotshared(a1, a))
    gerqf, = get_lapack_funcs(('gerqf',), (a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        rq, tau, work, info = gerqf(a1, lwork=-1, overwrite_a=1)
        lwork = work[0]
    rq, tau, work, info = gerqf(a1, lwork=lwork, overwrite_a=overwrite_a)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal geqrf'
                                                                    % -info)
    gemm, = get_blas_funcs(('gemm',), (rq,))
    t = rq.dtype.char
    R = special_matrices.triu(rq)
    Q = numpy.identity(M, dtype=t)
    ident = numpy.identity(M, dtype=t)
    zeros = numpy.zeros

    k = min(M, N)
    for i in range(k):
        v = zeros((M,), t)
        v[N-k+i] = 1
        v[0:N-k+i] = rq[M-k+i, 0:N-k+i]
        H = gemm(-tau[i], v, v, 1+0j, ident, trans_b=2)
        Q = gemm(1, Q, H)
    return R, Q
