#
# Author: Pearu Peterson, March 2002
#
# w/ additions by Travis Oliphant, March 2002

__all__ = ['solve','inv','det','lstsq','norm','pinv','pinv2',
           'tri','tril','triu','toeplitz','hankel','lu_solve',
           'cho_solve','solve_banded','LinAlgError','kron',
           'all_mat']

#from blas import get_blas_funcs
from lapack import get_lapack_funcs
from flinalg import get_flinalg_funcs
from scipy_base import asarray,zeros,sum,NewAxis,greater_equal,subtract,arange,\
     conjugate,ravel,r_,mgrid,take,ones,dot,transpose,diag,sqrt,add,real
import Matrix
import scipy_base
from scipy_base import asarray_chkfinite
import calc_lwork

class LinAlgError(Exception):
    pass

def lu_solve((lu, piv), b, trans=0, overwrite_b=0):
    """ lu_solve((lu, piv), b, trans=0, overwrite_b=0) -> x

    Solve a system of equations given a previously factored matrix

    Inputs:

      (lu,piv) -- The factored matrix, a (the output of lu_factor)
      b        -- a set of right-hand sides
      trans    -- type of system to solve:
                  0 : a   * x = b   (no transpose)
                  1 : a^T * x = b   (transpose) 
                  2   a^H * x = b   (conjugate transpose)

    Outputs:

       x -- the solution to the system
    """
    b1 = asarray_chkfinite(b)
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))
    if lu.shape[0] != b1.shape[0]:
        raise ValuError, "incompatible dimensions."
    getrs, = get_lapack_funcs(('getrs',),(lu,b1))
    x,info = getrs(lu,piv,b1,trans=trans,overwrite_b=overwrite_b)
    if info==0:
        return x
    raise ValueError,\
          'illegal value in %-th argument of internal gesv|posv'%(-info)

def cho_solve((c, lower), b, overwrite_b=0):
    """ cho_solve((c, lower), b, overwrite_b=0) -> x

    Solve a system of equations given a previously cholesky factored matrix

    Inputs:

      (c,lower) -- The factored matrix, a (the output of cho_factor)
      b        -- a set of right-hand sides

    Outputs:

       x -- the solution to the system a*x = b
    """
    b1 = asarray_chkfinite(b)
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))
    if c.shape[0] != b1.shape[0]:
        raise ValuError, "incompatible dimensions."
    potrs, = get_lapack_funcs(('potrs',),(c,b1))
    x,info = potrs(c,b1,lower=lower,overwrite_b=overwrite_b)
    if info==0:
        return x
    raise ValueError,\
          'illegal value in %-th argument of internal gesv|posv'%(-info)

# Linear equations
def solve(a, b, sym_pos=0, lower=0, overwrite_a=0, overwrite_b=0,
          debug = 0):
    """ solve(a, b, sym_pos=0, lower=0, overwrite_a=0, overwrite_b=0) -> x

    Solve a linear system of equations a * x = b for x.

    Inputs:

      a -- An N x N matrix.
      b -- An N x nrhs matrix or N vector.
      sym_pos -- Assume a is symmetric and positive definite.
      lower -- Assume a is lower triangular, otherwise upper one.
               Only used if sym_pos is true.
      overwrite_y - Discard data in y, where y is a or b.

    Outputs:

      x -- The solution to the system a * x = b
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
        posv, = get_lapack_funcs(('posv',),(a1,b1),debug=debug)
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
    """ solve_banded((l,u), ab, b, overwrite_ab=0, overwrite_b=0) -> x

    Solve a linear system of equations a * x = b for x where
    a is a banded matrix stored in diagonal orded form

     *   *     a1u
     
     *  a12 a23 ...
    a11 a22 a33 ...
    a21 a32 a43 ...
    .
    al1 ..         *

    Inputs:

      (l,u) -- number of non-zero lower and upper diagonals, respectively.
      a -- An N x (l+u+1) matrix.
      b -- An N x nrhs matrix or N vector.
      overwrite_y - Discard data in y, where y is ab or b.

    Outputs:

      x -- The solution to the system a * x = b
    """
    a1, b1 = map(asarray_chkfinite,(ab,b))
    overwrite_b = overwrite_b or (b1 is not b and not hasattr(b,'__array__'))

    gbsv, = get_lapack_funcs(('gbsv',),(a1,b1))
    a2 = zeros((2*l+u+1,a1.shape[1]),gbsv.typecode)
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


# matrix inversion
def inv(a, overwrite_a=0):
    """ inv(a, overwrite_a=0) -> a_inv

    Return inverse of square matrix a.
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


## matrix and Vector norm
import decomp
def norm(x, ord=2):
    """ norm(x, ord=2) -> n

    Matrix and vector norm.

    Inputs:

      x -- a rank-1 (vector) or rank-2 (matrix) array
      ord -- the order of norm.

     Comments:

       For vectors ord can be any real number including Inf or -Inf.
         ord = Inf, computes the maximum of the magnitudes
         ord = -Inf, computes minimum of the magnitudes
         ord is finite, computes sum(abs(x)**ord)**(1.0/ord)

       For matrices ord can only be + or - 1, 2, Inf.
         ord = 2 computes the largest singular value
         ord = -2 computes the smallest singular value
         ord = 1 computes the largest column sum of absolute values
         ord = -1 computes the smallest column sum of absolute values
         ord = Inf computes the largest row sum of absolute values
         ord = -Inf computes the smallest row sum of absolute values
         ord = 'fro' computes the frobenius norm sqrt(sum(diag(X.H * X)))
    """
    x = asarray_chkfinite(x)
    nd = len(x.shape)
    Inf = scipy_base.Inf
    if nd == 1:
        if ord == Inf:
            return scipy_base.amax(abs(x))
        elif ord == -Inf:
            return scipy_base.amin(abs(x))
        else:
            return scipy_base.sum(abs(x)**ord)**(1.0/ord)
    elif nd == 2:
        if ord == 2:
            return scipy_base.amax(decomp.svd(x,compute_uv=0))
        elif ord == -2:
            return scipy_base.amin(decomp.svd(x,compute_uv=0))
        elif ord == 1:
            return scipy_base.amax(scipy_base.sum(abs(x)))
        elif ord == Inf:
            return scipy_base.amax(scipy_base.sum(abs(x),axis=1))
        elif ord == -1:
            return scipy_base.amin(scipy_base.sum(abs(x)))
        elif ord == -Inf:
            return scipy_base.amin(scipy_base.sum(abs(x),axis=1))
        elif ord in ['fro','f']:
            val = real((conjugate(x)*x).flat)
            return sqrt(add.reduce(val))
        else:
            raise ValueError, "Invalid norm order for matrices."
    else:
        raise ValueError, "Improper number of dimensions to norm."

### Determinant

def det(a, overwrite_a=0):
    """ det(a, overwrite_a=0) -> d

    Return determinant of a square matrix.
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
    """ lstsq(a, b, cond=None, overwrite_a=0, overwrite_b=0) -> x,resids,rank,s

    Return least-squares solution of a * x = b.

    Inputs:

      a -- An M x N matrix.
      b -- An M x nrhs matrix or M vector.
      cond -- Used to determine effective rank of a.

    Outputs:

      x -- The solution (N x nrhs matrix) to the minimization problem:
                  2-norm(| b - a * x |) -> min
      resids -- The residual sum-of-squares for the solution matrix x
                (only if M>N and rank==N).
      rank -- The effective rank of a.
      s -- Singular values of a in decreasing order. The condition number
           of a is abs(s[0]/s[-1]).
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
        b2 = zeros((n,nrhs),gelss.typecode)
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
    resids = asarray([],x.typecode())
    if n<m:
        x1 = x[:n]
        if rank==n: resids = sum(x[n:]**2)
        x = x1
    return x,resids,rank,s


def pinv(a, cond=None):
    """ pinv(a, cond=None) -> a_pinv

    Compute generalized inverse of A using least-squares solver.
    """
    a = asarray_chkfinite(a)
    t = a.typecode()
    b = scipy_base.identity(a.shape[0],t)
    return lstsq(a, b, cond=cond)[0]


eps = scipy_base.limits.double_epsilon
feps = scipy_base.limits.float_epsilon
_array_precision = {'f': 0, 'd': 1, 'F': 0, 'D': 1}
def pinv2(a, cond=None):
    """ pinv2(a, cond=None) -> a_pinv

    Compute the generalized inverse of A using svd.
    """
    a = asarray_chkfinite(a)
    u, s, vh = decomp.svd(a)
    t = u.typecode()
    if cond in [None,-1]:
        cond = {0: feps*1e3, 1: eps*1e6}[_array_precision[t]]
    m,n = a.shape
    cutoff = cond*scipy_base.maximum.reduce(s)
    psigma = zeros((m,n),t)
    for i in range(len(s)):
        if s[i] > cutoff:
            psigma[i,i] = 1.0/conjugate(s[i])
    #XXX: use lapack/blas routines for dot
    return transpose(conjugate(dot(dot(u,psigma),vh)))

#-----------------------------------------------------------------------------
# matrix construction functions
#-----------------------------------------------------------------------------

def tri(N, M=None, k=0, typecode=None):
    """ returns a N-by-M matrix where all the diagonals starting from 
        lower left corner up to the k-th are all ones.
    """
    if M is None: M = N
    if type(M) == type('d'): 
        typecode = M
        M = N
    m = greater_equal(subtract.outer(arange(N), arange(M)),-k)
    if typecode is None:
        return m
    else:
        return m.astype(typecode)

def tril(m, k=0):
    """ returns the elements on and below the k-th diagonal of m.  k=0 is the
        main diagonal, k > 0 is above and k < 0 is below the main diagonal.
    """
    svsp = m.spacesaver()
    m = asarray(m,savespace=1)
    out = tri(m.shape[0], m.shape[1], k=k, typecode=m.typecode())*m
    out.savespace(svsp)
    return out

def triu(m, k=0):
    """ returns the elements on and above the k-th diagonal of m.  k=0 is the
        main diagonal, k > 0 is above and k < 0 is below the main diagonal.
    """
    svsp = m.spacesaver()
    m = asarray(m,savespace=1)
    out = (1-tri(m.shape[0], m.shape[1], k-1, m.typecode()))*m
    out.savespace(svsp)
    return out

def toeplitz(c,r=None):
    """ Construct a toeplitz matrix (i.e. a matrix with constant diagonals).

        Description:
    
           toeplitz(c,r) is a non-symmetric Toeplitz matrix with c as its first
           column and r as its first row.
    
           toeplitz(c) is a symmetric (Hermitian) Toeplitz matrix (r=c). 
    
        See also: hankel
    """
    isscalar = scipy_base.isscalar
    if isscalar(c) or isscalar(r):
        return c   
    if r is None:
        r = c
        r[0] = conjugate(r[0])
        c = conjugate(c)
    r,c = map(asarray_chkfinite,(r,c))
    r,c = map(ravel,(r,c))
    rN,cN = map(len,(r,c))
    if r[0] != c[0]:
        print "Warning: column and row values don't agree; column value used."
    vals = r_[r[rN-1:0:-1], c]
    cols = mgrid[0:cN]
    rows = mgrid[rN:0:-1]
    indx = cols[:,NewAxis]*ones((1,rN)) + \
           rows[NewAxis,:]*ones((cN,1)) - 1
    return take(vals, indx)


def hankel(c,r=None):
    """ Construct a hankel matrix (i.e. matrix with constant anti-diagonals).
    
        Description:
    
          hankel(c,r) is a Hankel matrix whose first column is c and whose
          last row is r.
    
          hankel(c) is a square Hankel matrix whose first column is C.
          Elements below the first anti-diagonal are zero.
    
        See also:  toeplitz
    """
    isscalar = scipy_base.isscalar
    if isscalar(c) or isscalar(r):
        return c   
    if r is None:
        r = zeros(len(c))
    elif r[0] != c[-1]:
        print "Warning: column and row values don't agree; column value used."
    r,c = map(asarray_chkfinite,(r,c))
    r,c = map(ravel,(r,c))
    rN,cN = map(len,(r,c))
    vals = r_[c, r[1:rN]]
    cols = mgrid[1:cN+1]
    rows = mgrid[0:rN]
    indx = cols[:,NewAxis]*ones((1,rN)) + \
           rows[NewAxis,:]*ones((cN,1)) - 1
    return take(vals, indx)

def all_mat(*args):
    return map(Matrix.Matrix,args)

def kron(a,b):
    """kronecker product of a and b

    Kronecker product of two matrices is block matrix
    [[ a[ 0 ,0]*b, a[ 0 ,1]*b, ... , a[ 0 ,n-1]*b  ],
     [ ...                                   ...   ],
     [ a[m-1,0]*b, a[m-1,1]*b, ... , a[m-1,n-1]*b  ]]
    """
    if not a.iscontiguous():
        a = reshape(a, a.shape)
    if not b.iscontiguous():
        b = reshape(b, b.shape)
    o = outerproduct(a,b)
    o.shape = a.shape + b.shape
    return concatenate(concatenate(o, axis=1), axis=1)
