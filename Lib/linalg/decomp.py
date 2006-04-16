## Automatically adapted for scipy Oct 18, 2005 by

#
# Author: Pearu Peterson, March 2002
#
# additions by Travis Oliphant, March 2002
# additions by Eric Jones,      June 2002

__all__ = ['eig','eigh','eigvals','eigvalsh','lu','svd','svdvals','diagsvd',
           'cholesky','qr',
           'schur','rsf2csf','lu_factor','cho_factor','cho_solve','orth',
           'hessenberg']

from basic import LinAlgError
import basic

from warnings import warn
from lapack import get_lapack_funcs
from blas import get_blas_funcs
from flinalg import get_flinalg_funcs
import calc_lwork
import numpy
from numpy import array, asarray_chkfinite, asarray, diag, zeros, ones, single, typecodes, isfinite

cast = numpy.cast
r_ = numpy.r_
c_ = numpy.c_

_I = cast['F'](1j)
def _make_complex_eigvecs(w,vin,cmplx_tcode):
    v = numpy.array(vin,dtype=cmplx_tcode)
    ind = numpy.nonzero(numpy.not_equal(w.imag,0.0))
    vnew = numpy.zeros((v.shape[0],len(ind)>>1),cmplx_tcode)
    vnew.real = numpy.take(vin,ind[::2],1)
    vnew.imag = numpy.take(vin,ind[1::2],1)
    count = 0
    conj = numpy.conjugate
    for i in range(len(ind)/2):
        v[:,ind[2*i]] = vnew[:,count]
        v[:,ind[2*i+1]] = conj(vnew[:,count])
        count += 1
    return v

def _datanotshared(a1,a):
    if a1 is a:
        return 0
    if hasattr(a,"__array_data__"):
        return a1.__array_data__[0] != a.__array_data__[0]
    return 1

def _geneig(a1,b,left,right,overwrite_a,overwrite_b):
    b1 = asarray(b)
    overwrite_b = overwrite_b or _datanotshared(b1,b)
    if len(b1.shape) != 2 or b1.shape[0] != b1.shape[1]:
        raise ValueError, 'expected square matrix'
    ggev, = get_lapack_funcs(('ggev',),(a1,b1))
    cvl,cvr = left,right
    if ggev.module_name[:7] == 'clapack':
        raise NotImplementedError,'calling ggev from %s' % (ggev.module_name)
    res = ggev(a1,b1,lwork=-1)
    lwork = res[-2][0]
    if ggev.prefix in 'cz':
        alpha,beta,vl,vr,work,info = ggev(a1,b1,cvl,cvr,lwork,
                                          overwrite_a,overwrite_b)
        w = alpha / beta
    else:
        alphar,alphai,beta,vl,vr,work,info = ggev(a1,b1,cvl,cvr,lwork,
                                                  overwrite_a,overwrite_b)
        w = (alphar+_I*alphai)/beta
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal ggev'%(-info)
    if info>0: raise LinAlgError,"generalized eig algorithm did not converge"

    only_real = numpy.logical_and.reduce(numpy.equal(w.imag,0.0))
    if not (ggev.prefix in 'cz' or only_real):
        t = w.dtype.char
        if left:
            vl = _make_complex_eigvecs(w, vl, t)
        if right:
            vr = _make_complex_eigvecs(w, vr, t)
    if not (left or right):
        return w
    if left:
        if right:
            return w, vl, vr
        return w, vl
    return w, vr

def eig(a,b=None,left=0,right=1,overwrite_a=0,overwrite_b=0):
    """ Solve ordinary and generalized eigenvalue problem
    of a square matrix.

    Inputs:

      a     -- An N x N matrix.
      b     -- An N x N matrix [default is identity(N)].
      left  -- Return left eigenvectors [disabled].
      right -- Return right eigenvectors [enabled].

    Outputs:

      w      -- eigenvalues [left==right==0].
      w,vr   -- w and right eigenvectors [left==0,right=1].
      w,vl   -- w and left eigenvectors [left==1,right==0].
      w,vl,vr  -- [left==right==1].

    Definitions:

      a * vr[:,i] = w[i] * b * vr[:,i]

      a^H * vl[:,i] = conjugate(w[i]) * b^H * vl[:,i]

    where a^H denotes transpose(conjugate(a)).
    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    if b is not None:
        b = asarray_chkfinite(b)
        return _geneig(a1,b,left,right,overwrite_a,overwrite_b)
    geev, = get_lapack_funcs(('geev',),(a1,))
    compute_vl,compute_vr=left,right
    if geev.module_name[:7] == 'flapack':
        lwork = calc_lwork.geev(geev.prefix,a1.shape[0],
                                compute_vl,compute_vr)[1]
        if geev.prefix in 'cz':
            w,vl,vr,info = geev(a1,lwork = lwork,
                                compute_vl=compute_vl,
                                compute_vr=compute_vr,
                                overwrite_a=overwrite_a)
        else:
            wr,wi,vl,vr,info = geev(a1,lwork = lwork,
                                    compute_vl=compute_vl,
                                    compute_vr=compute_vr,
                                    overwrite_a=overwrite_a)
            t = {'f':'F','d':'D'}[wr.dtype.char]
            w = wr+_I*wi
    else: # 'clapack'
        if geev.prefix in 'cz':
            w,vl,vr,info = geev(a1,
                                compute_vl=compute_vl,
                                compute_vr=compute_vr,
                                overwrite_a=overwrite_a)
        else:
            wr,wi,vl,vr,info = geev(a1,
                                    compute_vl=compute_vl,
                                    compute_vr=compute_vr,
                                    overwrite_a=overwrite_a)
            t = {'f':'F','d':'D'}[wr.dtype.char]
            w = wr+_I*wi
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal geev'%(-info)
    if info>0: raise LinAlgError,"eig algorithm did not converge"

    only_real = numpy.logical_and.reduce(numpy.equal(w.imag,0.0))
    if not (geev.prefix in 'cz' or only_real):
        t = w.dtype.char
        if left:
            vl = _make_complex_eigvecs(w, vl, t)
        if right:
            vr = _make_complex_eigvecs(w, vr, t)
    if not (left or right):
        return w
    if left:
        if right:
            return w, vl, vr
        return w, vl
    return w, vr

def eigh(a,lower=1,eigvals_only=0,overwrite_a=0):
    """ Solve real symmetric or complex hermitian eigenvalue problem.

    Inputs:

      a            -- A hermitian N x N matrix.
      lower        -- values in a are read from lower triangle [1: UPLO='L' (default) / 0: UPLO='U']
      eigvals_only -- don't compute eigenvectors.
      overwrite_a  -- content of a may be destroyed

    Outputs:

      w,v     -- w: eigenvalues, v: eigenvectors [for eigvals_only == False]
      w       -- eigenvalues [for eigvals_only == True (default)].

    Definitions:

      a * v[:,i] = w[i] * vr[:,i]
      v.H * v = identity

    """
    if eigvals_only or overwrite_a:
	a1 = asarray_chkfinite(a)
	overwrite_a = overwrite_a or (_datanotshared(a1,a))
    else:
	a1 = array(a)
        if (a1.dtype.char in typecodes['AllFloat']) and not isfinite(a1).all():
	    raise ValueError, "array must not contain infs or NaNs"
	overwrite_a = 1

    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError, 'expected square matrix'

    if a1.dtype.char in 'FD':
	heev, = get_lapack_funcs(('heev',),(a1,))
        if heev.module_name[:7] == 'flapack':
	    lwork = calc_lwork.heev(heev.prefix,a1.shape[0],lower)
	    w,v,info = heev(a1,lwork = lwork,
	                    compute_v = not eigvals_only,
			    lower = lower,
			    overwrite_a = overwrite_a)
	else: # 'clapack'
	    w,v,info = heev(a1,
	                    compute_v = not eigvals_only,
			    lower = lower,
			    overwrite_a = overwrite_a)
        if info<0: raise ValueError,\
	   'illegal value in %-th argument of internal heev'%(-info)
        if info>0: raise LinAlgError,"eig algorithm did not converge"
    else: # a1.dtype.char in 'fd':
	syev, = get_lapack_funcs(('syev',),(a1,))
        if syev.module_name[:7] == 'flapack':
	    lwork = calc_lwork.syev(syev.prefix,a1.shape[0],lower)
	    w,v,info = syev(a1,lwork = lwork,
	                    compute_v = not eigvals_only,
			    lower = lower,
			    overwrite_a = overwrite_a)
	else: # 'clapack'
	    w,v,info = syev(a1,
	                    compute_v = not eigvals_only,
			    lower = lower,
			    overwrite_a = overwrite_a)
        if info<0: raise ValueError,\
	   'illegal value in %-th argument of internal syev'%(-info)
        if info>0: raise LinAlgError,"eig algorithm did not converge"

    if eigvals_only:
        return w
    return w, v


def eigvals(a,b=None,overwrite_a=0):
    """Return eigenvalues of square matrix."""
    return eig(a,b=b,left=0,right=0,overwrite_a=overwrite_a)

def eigvalsh(a,lower=1,overwrite_a=0):
    """Return eigenvalues of hermitean or real symmetric matrix."""
    return eigh(a,lower=lower,eigvals_only=1,overwrite_a=overwrite_a)

def lu_factor(a, overwrite_a=0):
    """Return raw LU decomposition of a matrix and pivots, for use in solving
    a system of linear equations.

    Inputs:

      a  --- an NxN matrix

    Outputs:

      lu --- the lu factorization matrix
      piv --- an array of pivots
    """
    a1 = asarray(a)
    if len(a1.shape) != 2 or (a1.shape[0] != a1.shape[1]):
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    getrf, = get_lapack_funcs(('getrf',),(a1,))
    lu, piv, info = getrf(a,overwrite_a=overwrite_a)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal getrf (lu_factor)'%(-info)
    if info>0: warn("Diagonal number %d is exactly zero. Singular matrix." % info,
                    RuntimeWarning)
    return lu, piv

def lu_solve(a_lu_pivots,b):
    """Solve a previously factored system.  First input is a tuple (lu, pivots)
    which is the output to lu_factor.  Second input is the right hand side.
    """
    a_lu, pivots = a_lu_pivots
    a_lu = asarray_chkfinite(a_lu)
    pivots = asarray_chkfinite(pivots)
    b = asarray_chkfinite(b)
    _assert_squareness(a_lu)

    getrs, = get_lapack_funcs(('getrs',),(a_lu,))
    b, info = getrs(a_lu,pivots,b)
    if info < 0:
        msg = "Argument %d to lapack's ?getrs() has an illegal value." % info
        raise TypeError, msg
    if info > 0:
        msg = "Unknown error occured int ?getrs(): error code = %d" % info
        raise TypeError, msg
    return b


def lu(a,permute_l=0,overwrite_a=0):
    """Return LU decompostion of a matrix.

    Inputs:

      a     -- An M x N matrix.
      permute_l  -- Perform matrix multiplication p * l [disabled].

    Outputs:

      p,l,u  -- LU decomposition matrices of a [permute_l=0]
      pl,u   -- LU decomposition matrices of a [permute_l=1]

    Definitions:

      a = p * l * u    [permute_l=0]
      a = pl * u       [permute_l=1]

      p   -  An M x M permutation matrix
      l   -  An M x K lower triangular or trapezoidal matrix
             with unit-diagonal
      u   -  An K x N upper triangular or trapezoidal matrix
             K = min(M,N)
    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError, 'expected matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    flu, = get_flinalg_funcs(('lu',),(a1,))
    p,l,u,info = flu(a1,permute_l=permute_l,overwrite_a = overwrite_a)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal lu.getrf'%(-info)
    if permute_l:
        return l,u
    return p,l,u

def svd(a,compute_uv=1,overwrite_a=0):
    """Compute singular value decomposition (SVD) of matrix a.

    Description:

      Singular value decomposition of a matrix a is
        a = u * sigma * v^H,
      where v^H denotes conjugate(transpose(v)), u,v are unitary
      matrices, sigma is zero matrix with a main diagonal containing
      real non-negative singular values of the matrix a.

    Inputs:

      a -- An M x N matrix.
      compute_uv -- If zero, then only the vector of singular values
                    is returned.

    Outputs:

      u -- An M x M unitary matrix [compute_uv=1].
      s -- An min(M,N) vector of singular values in descending order,
           sigma = diagsvd(s).
      vh -- An N x N unitary matrix [compute_uv=1], vh = v^H.

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError, 'expected matrix'
    m,n = a1.shape
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    gesdd, = get_lapack_funcs(('gesdd',),(a1,))
    if gesdd.module_name[:7] == 'flapack':
        lwork = calc_lwork.gesdd(gesdd.prefix,m,n,compute_uv)[1]
        u,s,v,info = gesdd(a1,compute_uv = compute_uv, lwork = lwork,
                      overwrite_a = overwrite_a)
    else: # 'clapack'
        raise NotImplementedError,'calling gesdd from %s' % (gesdd.module_name)
    if info>0: raise LinAlgError, "SVD did not converge"
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal gesdd'%(-info)
    if compute_uv:
        return u,s,v
    else:
        return s

def svdvals(a,overwrite_a=0):
    """Return singular values of a matrix."""
    return svd(a,compute_uv=0,overwrite_a=overwrite_a)

def diagsvd(s,M,N):
    """Return sigma from singular values and original size M,N."""
    part = diag(s)
    typ = part.dtype.char
    MorN = len(s)
    if MorN == M:
        return c_[part,zeros((M,N-M),typ)]
    elif MorN == N:
        return r_[part,zeros((M-N,N),typ)]
    else:
        raise ValueError, "Length of s must be M or N."

def cholesky(a,lower=0,overwrite_a=0):
    """Compute Cholesky decomposition of matrix.

    Description:

      For a hermitian positive-definite matrix a return the
      upper-triangular (or lower-triangular if lower==1) matrix,
      u such that u^H * u = a (or l * l^H = a).

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or _datanotshared(a1,a)
    potrf, = get_lapack_funcs(('potrf',),(a1,))
    c,info = potrf(a1,lower=lower,overwrite_a=overwrite_a,clean=1)
    if info>0: raise LinAlgError, "matrix not positive definite"
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal potrf'%(-info)
    return c

def cho_factor(a, lower=0, overwrite_a=0):
    """ Compute Cholesky decomposition of matrix and return an object
    to be used for solving a linear system using cho_solve.
    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    potrf, = get_lapack_funcs(('potrf',),(a1,))
    c,info = potrf(a1,lower=lower,overwrite_a=overwrite_a,clean=0)
    if info>0: raise LinAlgError, "matrix not positive definite"
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal potrf'%(-info)
    return c, lower

def cho_solve(clow, b):
    """Solve a previously factored symmetric system of equations.
    First input is a tuple (LorU, lower) which is the output to cho_factor.
    Second input is the right-hand side.
    """
    c, lower = clow
    c = asarray_chkfinite(c)
    _assert_squareness(c)
    b = asarray_chkfinite(b)
    potrs, = get_lapack_funcs(('potrs',),(c,))
    b, info = potrs(c,b,lower)
    if info < 0:
        msg = "Argument %d to lapack's ?potrs() has an illegal value." % info
        raise TypeError, msg
    if info > 0:
        msg = "Unknown error occured int ?potrs(): error code = %d" % info
        raise TypeError, msg
    return b


def qr(a,overwrite_a=0,lwork=None):
    """QR decomposition of an M x N matrix a.

    Description:

      Find a unitary matrix, q, and an upper-trapezoidal matrix r
      such that q * r = a

    Inputs:

      a -- the matrix
      overwrite_a=0 -- if non-zero then discard the contents of a,
                     i.e. a is used as a work array if possible.

      lwork=None -- >= shape(a)[1]. If None (or -1) compute optimal
                    work array size.

    Outputs:

      q, r -- matrices such that q * r = a

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError, 'expected matrix'
    M,N = a1.shape
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    geqrf, = get_lapack_funcs(('geqrf',),(a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        qr,tau,work,info = geqrf(a1,lwork=-1,overwrite_a=1)
        lwork = work[0]
    qr,tau,work,info = geqrf(a1,lwork=lwork,overwrite_a=overwrite_a)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal geqrf'%(-info)
    gemm, = get_blas_funcs(('gemm',),(qr,))
    t = qr.dtype.char
    R = basic.triu(qr)
    Q = numpy.identity(M,dtype=t)
    ident = numpy.identity(M,dtype=t)
    zeros = numpy.zeros
    for i in range(min(M,N)):
        v = zeros((M,),t)
        v[i] = 1
        v[i+1:M] = qr[i+1:M,i]
        H = gemm(-tau[i],v,v,1+0j,ident,trans_b=2)
        Q = gemm(1,Q,H)
    return Q, R

_double_precision = ['i','l','d']

def schur(a,output='real',lwork=None,overwrite_a=0):
    """Compute Schur decomposition of matrix a.

    Description:

      Return T, Z such that a = Z * T * (Z**H) where Z is a
      unitary matrix and T is either upper-triangular or quasi-upper
      triangular for output='real'
    """
    if not output in ['real','complex','r','c']:
        raise ValueError, "argument must be 'real', or 'complex'"
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2 or (a1.shape[0] != a1.shape[1]):
        raise ValueError, 'expected square matrix'
    typ = a1.dtype.char
    if output in ['complex','c'] and typ not in ['F','D']:
        if typ in _double_precision:
            a1 = a1.astype('D')
            typ = 'D'
        else:
            a1 = a1.astype('F')
            typ = 'F'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    gees, = get_lapack_funcs(('gees',),(a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        result = gees(lambda x: None,a,lwork=-1)
        lwork = result[-2][0]
    result = gees(lambda x: None,a,lwork=result[-2][0],overwrite_a=overwrite_a)
    info = result[-1]
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal gees'%(-info)
    elif info>0: raise LinAlgError, "Schur form not found.  Possibly ill-conditioned."
    return result[0], result[-3]

eps = numpy.finfo(float).eps.item()
feps = numpy.finfo(single).eps.item()

_array_kind = {'b':0, 'h':0, 'B': 0, 'i':0, 'l': 0, 'f': 0, 'd': 0, 'F': 1, 'D': 1}
_array_precision = {'i': 1, 'l': 1, 'f': 0, 'd': 1, 'F': 0, 'D': 1}
_array_type = [['f', 'd'], ['F', 'D']]
def _commonType(*arrays):
    kind = 0
    precision = 0
    for a in arrays:
        t = a.dtype.char
        kind = max(kind, _array_kind[t])
        precision = max(precision, _array_precision[t])
    return _array_type[kind][precision]

def _castCopy(type, *arrays):
    cast_arrays = ()
    for a in arrays:
        if a.dtype.char == type:
            cast_arrays = cast_arrays + (a.copy(),)
        else:
            cast_arrays = cast_arrays + (a.astype(type),)
    if len(cast_arrays) == 1:
        return cast_arrays[0]
    else:
        return cast_arrays

def _assert_squareness(*arrays):
    for a in arrays:
        if max(a.shape) != min(a.shape):
            raise LinAlgError, 'Array must be square'

def rsf2csf(T, Z):
    """Convert real schur form to complex schur form.

    Description:

      If A is a real-valued matrix, then the real schur form is
      quasi-upper triangular.  2x2 blocks extrude from the main-diagonal
      corresponding to any complex-valued eigenvalues.

      This function converts this real schur form to a complex schur form
      which is upper triangular.
    """
    Z,T = map(asarray_chkfinite, (Z,T))
    if len(Z.shape) !=2 or Z.shape[0] != Z.shape[1]:
        raise ValueError, "matrix must be square."
    if len(T.shape) !=2 or T.shape[0] != T.shape[1]:
        raise ValueError, "matrix must be square."
    if T.shape[0] != Z.shape[0]:
        raise ValueError, "matrices must be same dimension."
    N = T.shape[0]
    arr = numpy.array
    t = _commonType(Z, T, arr([3.0],'F'))
    Z, T = _castCopy(t, Z, T)
    conj = numpy.conj
    dot = numpy.dot
    r_ = numpy.r_
    transp = numpy.transpose
    for m in range(N-1,0,-1):
        if abs(T[m,m-1]) > eps*(abs(T[m-1,m-1]) + abs(T[m,m])):
            k = slice(m-1,m+1)
            mu = eigvals(T[k,k]) - T[m,m]
            r = basic.norm([mu[0], T[m,m-1]])
            c = mu[0] / r
            s = T[m,m-1] / r
            G = r_[arr([[conj(c),s]],dtype=t),arr([[-s,c]],dtype=t)]
            Gc = conj(transp(G))
            j = slice(m-1,N)
            T[k,j] = dot(G,T[k,j])
            i = slice(0,m+1)
            T[i,k] = dot(T[i,k], Gc)
            i = slice(0,N)
            Z[i,k] = dot(Z[i,k], Gc)
        T[m,m-1] = 0.0;
    return T, Z


# Orthonormal decomposition

def orth(A):
    """Return an orthonormal basis for the range of A using svd"""
    u,s,vh = svd(A)
    M,N = A.shape
    tol = max(M,N)*numpy.amax(s)*eps
    num = numpy.sum(s > tol,dtype=int)
    Q = u[:,:num]
    return Q

def hessenberg(a,calc_q=0,overwrite_a=0):
    """ Compute Hessenberg form of a matrix.

    Inputs:

      a -- the matrix
      calc_q -- if non-zero then calculate unitary similarity
                transformation matrix q.
      overwrite_a=0 -- if non-zero then discard the contents of a,
                     i.e. a is used as a work array if possible.

    Outputs:

      h    -- Hessenberg form of a                [calc_q=0]
      h, q -- matrices such that a = q * h * q^T  [calc_q=1]

    """
    a1 = asarray(a)
    if len(a1.shape) != 2 or (a1.shape[0] != a1.shape[1]):
        raise ValueError, 'expected square matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    gehrd,gebal = get_lapack_funcs(('gehrd','gebal'),(a1,))
    ba,lo,hi,pivscale,info = gebal(a,permute=1,overwrite_a = overwrite_a)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal gebal (hessenberg)'%(-info)
    n = len(a1)
    lwork = calc_lwork.gehrd(gehrd.prefix,n,lo,hi)
    hq,tau,info = gehrd(ba,lo=lo,hi=hi,lwork=lwork,overwrite_a=1)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal gehrd (hessenberg)'%(-info)

    if not calc_q:
        for i in range(lo,hi):
            hq[i+2:hi+1,i] = 0.0
        return hq

    # XXX: Use ORGHR routines to compute q.
    ger,gemm = get_blas_funcs(('ger','gemm'),(hq,))
    typecode = hq.dtype.char
    q = None
    for i in range(lo,hi):
        if tau[i]==0.0:
            continue
        v = zeros(n,dtype=typecode)
        v[i+1] = 1.0
        v[i+2:hi+1] = hq[i+2:hi+1,i]
        hq[i+2:hi+1,i] = 0.0
        h = ger(-tau[i],v,v,a=diag(ones(n,dtype=typecode)),overwrite_a=1)
        if q is None:
            q = h
        else:
            q = gemm(1.0,q,h)
    if q is None:
        q = diag(ones(n,dtype=typecode))
    return hq,q
