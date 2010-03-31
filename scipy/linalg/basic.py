## Automatically adapted for scipy Oct 18, 2005 by

## Automatically adapted for scipy Oct 18, 2005 by

#
# Author: Pearu Peterson, March 2002
#
# w/ additions by Travis Oliphant, March 2002

__all__ = ['solve', 'inv', 'det', 'lstsq', 'pinv', 'pinv2',
           'cholesky_banded', 'solveh_banded', 'lu_solve', 'cho_solve',
           'solve_banded',
           #  From special_matrices:
           'tri','tril', 'triu', 'toeplitz', 'circulant', 'hankel', 'kron',
           'hadamard', 'block_diag', 'all_mat',
           #  From misc:
           'LinAlgError', 'norm',
            ]

from numpy import asarray, zeros, sum, conjugate, dot, transpose, \
        asarray_chkfinite,  single
import numpy

#from blas import get_blas_funcs
from flinalg import get_flinalg_funcs
from lapack import get_lapack_funcs
from misc import LinAlgError, norm
from special_matrices import tri, tril, triu, toeplitz, circulant, hankel, \
        hadamard, kron, block_diag, all_mat
from scipy.linalg import calc_lwork
import decomp


def lu_solve((lu, piv), b, trans=0, overwrite_b=0):
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
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))
    if lu.shape[0] != b1.shape[0]:
        raise ValueError, "incompatible dimensions."
    getrs, = get_lapack_funcs(('getrs',),(lu,b1))
    x,info = getrs(lu,piv,b1,trans=trans,overwrite_b=overwrite_b)
    if info==0:
        return x
    raise ValueError,\
          'illegal value in %-th argument of internal gesv|posv'%(-info)

def cho_solve((c, lower), b, overwrite_b=0):
    """Solve an equation system, a x = b, given the Cholesky factorization of a

    Parameters
    ----------
    (c, lower)
        Cholesky factorization of a, as given by cho_factor
    b : array
        Right-hand side

    Returns
    -------
    x : array
        The solution to the system a x = b

    See also
    --------
    cho_factor : Cholesky factorization of a matrix

    """
    b1 = asarray_chkfinite(b)
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))
    if c.shape[0] != b1.shape[0]:
        raise ValueError, "incompatible dimensions."
    potrs, = get_lapack_funcs(('potrs',),(c,b1))
    x,info = potrs(c,b1,lower=lower,overwrite_b=overwrite_b)
    if info==0:
        return x
    raise ValueError,\
          'illegal value in %-th argument of internal gesv|posv'%(-info)

# Linear equations
def solve(a, b, sym_pos=0, lower=0, overwrite_a=0, overwrite_b=0,
          debug = 0):
    """Solve the equation a x = b for x

    Parameters
    ----------
    a : array, shape (M, M)
    b : array, shape (M,) or (M, N)
    sym_pos : boolean
        Assume a is symmetric and positive definite
    lower : boolean
        Use only data contained in the lower triangle of a, if sym_pos is true.
        Default is to use upper triangle.
    overwrite_a : boolean
        Allow overwriting data in a (may enhance performance)
    overwrite_b : boolean
        Allow overwriting data in b (may enhance performance)

    Returns
    -------
    x : array, shape (M,) or (M, N) depending on b
        Solution to the system a x = b

    Raises LinAlgError if a is singular

    """
    a1, b1 = map(asarray_chkfinite,(a,b))
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError, 'expected square matrix'
    if a1.shape[0] != b1.shape[0]:
        raise ValueError, 'incompatible dimensions'
    overwrite_a = overwrite_a or (a1 is not a and not hasattr(a,'__array__'))
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))
    if debug:
        print 'solve:overwrite_a=',overwrite_a
        print 'solve:overwrite_b=',overwrite_b
    if sym_pos:
        posv, = get_lapack_funcs(('posv',),(a1,b1))
        c,x,info = posv(a1,b1,
                        lower = lower,
                        overwrite_a=overwrite_a,
                        overwrite_b=overwrite_b)
    else:
        gesv, = get_lapack_funcs(('gesv',),(a1,b1))
        lu,piv,x,info = gesv(a1,b1,
                             overwrite_a=overwrite_a,
                             overwrite_b=overwrite_b)

    if info==0:
        return x
    if info>0:
        raise LinAlgError, "singular matrix"
    raise ValueError,\
          'illegal value in %-th argument of internal gesv|posv'%(-info)

def solve_banded((l,u), ab, b, overwrite_ab=0, overwrite_b=0,
          debug = 0):
    """Solve the equation a x = b for x, assuming a is banded matrix.

    The matrix a is stored in ab using the matrix diagonal orded form::

        ab[u + i - j, j] == a[i,j]

    Example of ab (shape of a is (6,6), u=1, l=2)::

        *    a01  a12  a23  a34  a45
        a00  a11  a22  a33  a44  a55
        a10  a21  a32  a43  a54   *
        a20  a31  a42  a53   *    *

    Parameters
    ----------
    (l, u) : (integer, integer)
        Number of non-zero lower and upper diagonals
    ab : array, shape (l+u+1, M)
        Banded matrix
    b : array, shape (M,) or (M, K)
        Right-hand side
    overwrite_ab : boolean
        Discard data in ab (may enhance performance)
    overwrite_b : boolean
        Discard data in b (may enhance performance)

    Returns
    -------
    x : array, shape (M,) or (M, K)
        The solution to the system a x = b

    """
    a1, b1 = map(asarray_chkfinite,(ab,b))
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))

    gbsv, = get_lapack_funcs(('gbsv',),(a1,b1))
    a2 = zeros((2*l+u+1,a1.shape[1]), dtype=gbsv.dtype)
    a2[l:,:] = a1
    lu,piv,x,info = gbsv(l,u,a2,b1,
                         overwrite_ab=1,
                         overwrite_b=overwrite_b)
    if info==0:
        return x
    if info>0:
        raise LinAlgError, "singular matrix"
    raise ValueError,\
          'illegal value in %-th argument of internal gbsv'%(-info)

def solveh_banded(ab, b, overwrite_ab=0, overwrite_b=0,
                  lower=0):
    """Solve equation a x = b. a is Hermitian positive-definite banded matrix.

    The matrix a is stored in ab either in lower diagonal or upper
    diagonal ordered form:

        ab[u + i - j, j] == a[i,j]        (if upper form; i <= j)
        ab[    i - j, j] == a[i,j]        (if lower form; i >= j)

    Example of ab (shape of a is (6,6), u=2)::

        upper form:
        *   *   a02 a13 a24 a35
        *   a01 a12 a23 a34 a45
        a00 a11 a22 a33 a44 a55

        lower form:
        a00 a11 a22 a33 a44 a55
        a10 a21 a32 a43 a54 *
        a20 a31 a42 a53 *   *

    Cells marked with * are not used.

    Parameters
    ----------
    ab : array, shape (M, u + 1)
        Banded matrix
    b : array, shape (M,) or (M, K)
        Right-hand side
    overwrite_ab : boolean
        Discard data in ab (may enhance performance)
    overwrite_b : boolean
        Discard data in b (may enhance performance)
    lower : boolean
        Is the matrix in the lower form. (Default is upper form)

    Returns
    -------
    c : array, shape (M, u+1)
        Cholesky factorization of a, in the same banded format as ab
    x : array, shape (M,) or (M, K)
        The solution to the system a x = b

    """
    ab, b = map(asarray_chkfinite,(ab,b))

    pbsv, = get_lapack_funcs(('pbsv',),(ab,b))
    c,x,info = pbsv(ab,b,
                    lower=lower,
                    overwrite_ab=overwrite_ab,
                    overwrite_b=overwrite_b)
    if info==0:
        return c, x
    if info>0:
        raise LinAlgError, "%d-th leading minor not positive definite" % info
    raise ValueError,\
          'illegal value in %d-th argument of internal pbsv'%(-info)

def cholesky_banded(ab, overwrite_ab=0, lower=0):
    """Cholesky decompose a banded Hermitian positive-definite matrix

    The matrix a is stored in ab either in lower diagonal or upper
    diagonal ordered form:

        ab[u + i - j, j] == a[i,j]        (if upper form; i <= j)
        ab[    i - j, j] == a[i,j]        (if lower form; i >= j)

    Example of ab (shape of a is (6,6), u=2)::

        upper form:
        *   *   a02 a13 a24 a35
        *   a01 a12 a23 a34 a45
        a00 a11 a22 a33 a44 a55

        lower form:
        a00 a11 a22 a33 a44 a55
        a10 a21 a32 a43 a54 *
        a20 a31 a42 a53 *   *

    Parameters
    ----------
    ab : array, shape (M, u + 1)
        Banded matrix
    overwrite_ab : boolean
        Discard data in ab (may enhance performance)
    lower : boolean
        Is the matrix in the lower form. (Default is upper form)

    Returns
    -------
    c : array, shape (M, u+1)
        Cholesky factorization of a, in the same banded format as ab

    """
    ab = asarray_chkfinite(ab)

    pbtrf, = get_lapack_funcs(('pbtrf',),(ab,))
    c,info = pbtrf(ab,
                   lower=lower,
                   overwrite_ab=overwrite_ab)

    if info==0:
        return c
    if info>0:
        raise LinAlgError, "%d-th leading minor not positive definite" % info
    raise ValueError,\
          'illegal value in %d-th argument of internal pbtrf'%(-info)


# matrix inversion
def inv(a, overwrite_a=0):
    """Compute the inverse of a matrix.

    Parameters
    ----------
    a : array-like, shape (M, M)
        Matrix to be inverted

    Returns
    -------
    ainv : array-like, shape (M, M)
        Inverse of the matrix a

    Raises LinAlgError if a is singular

    Examples
    --------
    >>> a = array([[1., 2.], [3., 4.]])
    >>> inv(a)
    array([[-2. ,  1. ],
           [ 1.5, -0.5]])
    >>> dot(a, inv(a))
    array([[ 1.,  0.],
           [ 0.,  1.]])

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or (a1 is not a and not hasattr(a,'__array__'))
    #XXX: I found no advantage or disadvantage of using finv.
##     finv, = get_flinalg_funcs(('inv',),(a1,))
##     if finv is not None:
##         a_inv,info = finv(a1,overwrite_a=overwrite_a)
##         if info==0:
##             return a_inv
##         if info>0: raise LinAlgError, "singular matrix"
##         if info<0: raise ValueError,\
##            'illegal value in %-th argument of internal inv.getrf|getri'%(-info)
    getrf,getri = get_lapack_funcs(('getrf','getri'),(a1,))
    #XXX: C ATLAS versions of getrf/i have rowmajor=1, this could be
    #     exploited for further optimization. But it will be probably
    #     a mess. So, a good testing site is required before trying
    #     to do that.
    if getrf.module_name[:7]=='clapack'!=getri.module_name[:7]:
        # ATLAS 3.2.1 has getrf but not getri.
        lu,piv,info = getrf(transpose(a1),
                            rowmajor=0,overwrite_a=overwrite_a)
        lu = transpose(lu)
    else:
        lu,piv,info = getrf(a1,overwrite_a=overwrite_a)
    if info==0:
        if getri.module_name[:7] == 'flapack':
            lwork = calc_lwork.getri(getri.prefix,a1.shape[0])
            lwork = lwork[1]
            # XXX: the following line fixes curious SEGFAULT when
            # benchmarking 500x500 matrix inverse. This seems to
            # be a bug in LAPACK ?getri routine because if lwork is
            # minimal (when using lwork[0] instead of lwork[1]) then
            # all tests pass. Further investigation is required if
            # more such SEGFAULTs occur.
            lwork = int(1.01*lwork)
            inv_a,info = getri(lu,piv,
                               lwork=lwork,overwrite_lu=1)
        else: # clapack
            inv_a,info = getri(lu,piv,overwrite_lu=1)
    if info>0: raise LinAlgError, "singular matrix"
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal getrf|getri'%(-info)
    return inv_a


### Determinant

def det(a, overwrite_a=0):
    """Compute the determinant of a matrix

    Parameters
    ----------
    a : array, shape (M, M)

    Returns
    -------
    det : float or complex
        Determinant of a

    Notes
    -----
    The determinant is computed via LU factorization, LAPACK routine z/dgetrf.
    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or (a1 is not a and not hasattr(a,'__array__'))
    fdet, = get_flinalg_funcs(('det',),(a1,))
    a_det,info = fdet(a1,overwrite_a=overwrite_a)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal det.getrf'%(-info)
    return a_det

### Linear Least Squares

def lstsq(a, b, cond=None, overwrite_a=0, overwrite_b=0):
    """Compute least-squares solution to equation :m:`a x = b`

    Compute a vector x such that the 2-norm :m:`|b - a x|` is minimised.

    Parameters
    ----------
    a : array, shape (M, N)
    b : array, shape (M,) or (M, K)
    cond : float
        Cutoff for 'small' singular values; used to determine effective
        rank of a. Singular values smaller than rcond*largest_singular_value
        are considered zero.
    overwrite_a : boolean
        Discard data in a (may enhance performance)
    overwrite_b : boolean
        Discard data in b (may enhance performance)

    Returns
    -------
    x : array, shape (N,) or (N, K) depending on shape of b
        Least-squares solution
    residues : array, shape () or (1,) or (K,)
        Sums of residues, squared 2-norm for each column in :m:`b - a x`
        If rank of matrix a is < N or > M this is an empty array.
        If b was 1-d, this is an (1,) shape array, otherwise the shape is (K,)
    rank : integer
        Effective rank of matrix a
    s : array, shape (min(M,N),)
        Singular values of a. The condition number of a is abs(s[0]/s[-1]).

    Raises LinAlgError if computation does not converge

    """
    a1, b1 = map(asarray_chkfinite,(a,b))
    if len(a1.shape) != 2:
        raise ValueError, 'expected matrix'
    m,n = a1.shape
    if len(b1.shape)==2: nrhs = b1.shape[1]
    else: nrhs = 1
    if m != b1.shape[0]:
        raise ValueError, 'incompatible dimensions'
    gelss, = get_lapack_funcs(('gelss',),(a1,b1))
    if n>m:
        # need to extend b matrix as it will be filled with
        # a larger solution matrix
        b2 = zeros((n,nrhs), dtype=gelss.dtype)
        if len(b1.shape)==2: b2[:m,:] = b1
        else: b2[:m,0] = b1
        b1 = b2
    overwrite_a = overwrite_a or (a1 is not a and not hasattr(a,'__array__'))
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))
    if gelss.module_name[:7] == 'flapack':
        lwork = calc_lwork.gelss(gelss.prefix,m,n,nrhs)[1]
        v,x,s,rank,info = gelss(a1,b1,cond = cond,
                                lwork = lwork,
                                overwrite_a = overwrite_a,
                                overwrite_b = overwrite_b)
    else:
        raise NotImplementedError,'calling gelss from %s' % (gelss.module_name)
    if info>0: raise LinAlgError, "SVD did not converge in Linear Least Squares"
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal gelss'%(-info)
    resids = asarray([], dtype=x.dtype)
    if n<m:
        x1 = x[:n]
        if rank==n: resids = sum(x[n:]**2,axis=0)
        x = x1
    return x,resids,rank,s


def pinv(a, cond=None, rcond=None):
    """Compute the (Moore-Penrose) pseudo-inverse of a matrix.

    Calculate a generalized inverse of a matrix using a least-squares
    solver.

    Parameters
    ----------
    a : array, shape (M, N)
        Matrix to be pseudo-inverted
    cond, rcond : float
        Cutoff for 'small' singular values in the least-squares solver.
        Singular values smaller than rcond*largest_singular_value are
        considered zero.

    Returns
    -------
    B : array, shape (N, M)

    Raises LinAlgError if computation does not converge

    Examples
    --------
    >>> from numpy import *
    >>> a = random.randn(9, 6)
    >>> B = linalg.pinv(a)
    >>> allclose(a, dot(a, dot(B, a)))
    True
    >>> allclose(B, dot(B, dot(a, B)))
    True

    """
    a = asarray_chkfinite(a)
    b = numpy.identity(a.shape[0], dtype=a.dtype)
    if rcond is not None:
        cond = rcond
    return lstsq(a, b, cond=cond)[0]


eps = numpy.finfo(float).eps
feps = numpy.finfo(single).eps

_array_precision = {'f': 0, 'd': 1, 'F': 0, 'D': 1}

def pinv2(a, cond=None, rcond=None):
    """Compute the (Moore-Penrose) pseudo-inverse of a matrix.

    Calculate a generalized inverse of a matrix using its
    singular-value decomposition and including all 'large' singular
    values.

    Parameters
    ----------
    a : array, shape (M, N)
        Matrix to be pseudo-inverted
    cond, rcond : float or None
        Cutoff for 'small' singular values.
        Singular values smaller than rcond*largest_singular_value are
        considered zero.

        If None or -1, suitable machine precision is used.

    Returns
    -------
    B : array, shape (N, M)

    Raises LinAlgError if SVD computation does not converge

    Examples
    --------
    >>> from numpy import *
    >>> a = random.randn(9, 6)
    >>> B = linalg.pinv2(a)
    >>> allclose(a, dot(a, dot(B, a)))
    True
    >>> allclose(B, dot(B, dot(a, B)))
    True

    """
    a = asarray_chkfinite(a)
    u, s, vh = decomp.svd(a)
    t = u.dtype.char
    if rcond is not None:
        cond = rcond
    if cond in [None,-1]:
        cond = {0: feps*1e3, 1: eps*1e6}[_array_precision[t]]
    m,n = a.shape
    cutoff = cond*numpy.maximum.reduce(s)
    psigma = zeros((m,n),t)
    for i in range(len(s)):
        if s[i] > cutoff:
            psigma[i,i] = 1.0/conjugate(s[i])
    #XXX: use lapack/blas routines for dot
    return transpose(conjugate(dot(dot(u,psigma),vh)))
