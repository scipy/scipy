import warnings

import numpy as np
from numpy import asarray_chkfinite, single

from misc import LinAlgError, _datacopied
from lapack import get_lapack_funcs

__all__ = ['qz']

_double_precision = ['i','l','d']

def _select_function(sort, typ):
    if typ in ['F','D']:
        if callable(sort):
            #assume the user knows what they're doing
            sfunction = sort
        elif sort == 'lhp':
            sfunction = lambda x,y: (np.real(x/y) < 0.0)
        elif sort == 'rhp':
            sfunction = lambda x,y: (np.real(x/y) >= 0.0)
        elif sort == 'iuc':
            sfunction = lambda x,y: (abs(x/y) <= 1.0)
        elif sort == 'ouc':
            sfunction = lambda x,y: (abs(x/y) > 1.0)
        else:
            raise ValueError("sort parameter must be None, a callable, or "
                "one of ('lhp','rhp','iuc','ouc')")
    elif typ in ['f','d']:
        if callable(sort):
            #assume the user knows what they're doing
            sfunction = sort
        elif sort == 'lhp':
            sfunction = lambda x,y,z: (np.real((x+y*1j)/z) < 0.0)
        elif sort == 'rhp':
            sfunction = lambda x,y,z: (np.real((x+y*1j)/z) >= 0.0)
        elif sort == 'iuc':
            sfunction = lambda x,y,z: (abs((x+y*1j)/z) <= 1.0)
        elif sort == 'ouc':
            sfunction = lambda x,y,z: (abs((x+y*1j)/z) > 1.0)
        else:
            raise ValueError("sort parameter must be None, a callable, or "
                "one of ('lhp','rhp','iuc','ouc')")
    else: # to avoid an error later
        raise ValueError("dtype %s not understood" % typ)
    return sfunction


def qz(A, B, output='real', lwork=None, sort=None, overwrite_a=False,
       overwrite_b=False):
    """
    QZ decompostion for generalized eigenvalues of a pair of matrices.

    The QZ, or generalized Schur, decomposition for a pair of N x N
    nonsymmetric matrices (A,B) is

        (A,B) = (Q*AA*Z', Q*BB*Z')

    where AA, BB is in generalized Schur form if BB is upper-triangular
    with non-negative diagonal and AA is upper-triangular, or for real QZ
    decomposition (output='real') block upper triangular with 1x1
    and 2x2 blocks. In this case, the 1x1 blocks correpsond to real
    generalized eigenvalues and 2x2 blocks are 'standardized' by making
    the correpsonding elements of BB have the form::

        [ a 0 ]
        [ 0 b ]

    and the pair of correpsonding 2x2 blocks in AA and BB will have a complex
    conjugate pair of generalized eigenvalues. If (output='complex') or A
    and B are complex matrices, Z' denotes the conjugate-transpose of Z.
    Q and Z are unitary matrices.

    Parameters
    ----------
    A : array_like, shape (N,N)
        2d array to decompose
    B : array_like, shape (N,N)
        2d array to decompose
    output : str {'real','complex'}
        Construct the real or complex QZ decomposition for real matrices.
    lwork : integer, optional
        Work array size. If None or -1, it is automatically computed.
    sort : {None, callable, 'lhp', 'rhp', 'iuc', 'ouc'}
        Specifies whether the upper eigenvalues should be sorted.  A callable
        may be passed that, given a eigenvalue, returns a boolean denoting
        whether the eigenvalue should be sorted to the top-left (True). For
        real matrix pairs, the sort function takes three real arguments
        (alphar, alphai, beta). The eigenvalue x = (alphar + alphai*1j)/beta.
        For complex matrix pairs or output='complex', the sort function
        takes two complex arguments (alpha, beta). The eigenvalue
        x = (alpha/beta).
        Alternatively, string parameters may be used:
            'lhp'   Left-hand plane (x.real < 0.0)
            'rhp'   Right-hand plane (x.real > 0.0)
            'iuc'   Inside the unit circle (x*x.conjugate() <= 1.0)
            'ouc'   Outside the unit circle (x*x.conjugate() > 1.0)
        Defaults to None (no sorting).

    Returns
    -------
    AA : ndarray, shape (N,N)
        Generalized Schur form of A.
    BB : ndarray, shape (N,N)
        Generalized Schur form of B.
    Q : ndarray, shape (N,N)
        The left Schur vectors.
    Z : ndarray, shape (N,N)
        The right Schur vectors.
    sdim : int
        If sorting was requested, a fifth return value will contain the
        number of eigenvalues for which the sort condition was True.

    Notes
    -----
    Q is transposed versus the equivalent function in Matlab.

    .. versionadded:: 0.11.0

    """
    if not output in ['real','complex','r','c']:
        raise ValueError("argument must be 'real', or 'complex'")

    a1 = asarray_chkfinite(A)
    b1 = asarray_chkfinite(B)

    a_m, a_n = a1.shape
    b_m, b_n = b1.shape
    try:
        assert a_m == a_n == b_m == b_n
    except AssertionError:
        raise ValueError("Array dimensions must be square and agree")

    typa = a1.dtype.char
    if output in ['complex', 'c'] and typa not in ['F','D']:
        if typa in _double_precision:
            a1 = a1.astype('D')
            typa = 'D'
        else:
            a1 = a1.astype('F')
            typa = 'F'
    typb = b1.dtype.char
    if output in ['complex', 'c'] and typb not in ['F','D']:
        if typb in _double_precision:
            b1 = b1.astype('D')
            typb = 'D'
        else:
            b1 = b1.astype('F')
            typb = 'F'

    overwrite_a = overwrite_a or (_datacopied(a1,A))
    overwrite_b = overwrite_b or (_datacopied(b1,B))

    gges, = get_lapack_funcs(('gges',), (a1,b1))

    if lwork is None or lwork == -1:
        # get optimal work array size
        result = gges(lambda x: None, a1, b1, lwork=-1)
        lwork = result[-2][0].real.astype(np.int)

    if sort is None:
        sort_t = 0
        sfunction = lambda x : None
    else:
        sort_t = 1
        sfunction = _select_function(sort, typa)

    result = gges(sfunction, a1, b1, lwork=lwork, overwrite_a=overwrite_a,
                  overwrite_b=overwrite_b, sort_t=sort_t)

    info = result[-1]
    if info < 0:
        raise ValueError("Illegal value in argument %d of gges" % -info)
    elif info > 0 and info <= a_n:
        warnings.warn("The QZ iteration failed. (a,b) are not in Schur "
                "form, but ALPHAR(j), ALPHAI(j), and BETA(j) should be correct"
                "for J=%d,...,N" % info-1, UserWarning)
    elif info == a_n+1:
        raise LinAlgError("Something other than QZ iteration failed")
    elif info == a_n+2:
        raise LinAlgError("After reordering, roundoff changed values of some"
                "complex eigenvalues so that leading eigenvalues in the"
                "Generalized Schur form no longer satisfy sort=True."
                "This could also be caused due to scaling.")
    elif info == a_n+3:
        raise LinAlgError("Reordering failed in <s,d,c,z>tgsen")

    # output for real
    #AA, BB, sdim, alphar, alphai, beta, vsl, vsr, work, info
    # output for complex
    #AA, BB, sdim, alphai, beta, vsl, vsr, work, info
    if sort_t == 0:
        return result[0], result[1], result[-4], result[-3]
    else:
        return result[0], result[1], result[-4], result[-3], result[2]
