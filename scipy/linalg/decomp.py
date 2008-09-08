## Automatically adapted for scipy Oct 18, 2005 by

#
# Author: Pearu Peterson, March 2002
#
# additions by Travis Oliphant, March 2002
# additions by Eric Jones,      June 2002
# additions by Johannes Loehnert, June 2006
# additions by Bart Vandereycken, June 2006
# additions by Andrew D Straw, May 2007

__all__ = ['eig','eigh','eig_banded','eigvals','eigvalsh', 'eigvals_banded',
           'lu','svd','svdvals','diagsvd','cholesky','qr','qr_old','rq',
           'schur','rsf2csf','lu_factor','cho_factor','cho_solve','orth',
           'hessenberg']

from basic import LinAlgError
import basic

from warnings import warn
from lapack import get_lapack_funcs, find_best_lapack_type
from blas import get_blas_funcs
from flinalg import get_flinalg_funcs
from scipy.linalg import calc_lwork
import numpy
from numpy import array, asarray_chkfinite, asarray, diag, zeros, ones, \
        single, isfinite, inexact, complexfloating

cast = numpy.cast
r_ = numpy.r_

_I = cast['F'](1j)
def _make_complex_eigvecs(w,vin,cmplx_tcode):
    v = numpy.array(vin,dtype=cmplx_tcode)
    #ind = numpy.flatnonzero(numpy.not_equal(w.imag,0.0))
    ind = numpy.flatnonzero(numpy.logical_and(numpy.not_equal(w.imag,0.0),
                            numpy.isfinite(w)))
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
        return False
    else:
        #try comparing data pointers
        try:
            return a1.__array_interface__['data'][0] != a.__array_interface__['data'][0]
        except:
            return True


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

def eig(a,b=None, left=False, right=True, overwrite_a=False, overwrite_b=False):
    """Solve an ordinary or generalized eigenvalue problem of a square matrix.

    Find eigenvalues w and right or left eigenvectors of a general matrix::

        a   vr[:,i] = w[i]        b   vr[:,i]
        a.H vl[:,i] = w[i].conj() b.H vl[:,i]

    where .H is the Hermitean conjugation.

    Parameters
    ----------
    a : array, shape (M, M)
        A complex or real matrix whose eigenvalues and eigenvectors
        will be computed.
    b : array, shape (M, M)
        Right-hand side matrix in a generalized eigenvalue problem.
        If omitted, identity matrix is assumed.
    left : boolean
        Whether to calculate and return left eigenvectors
    right : boolean
        Whether to calculate and return right eigenvectors

    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)
    overwrite_b : boolean
        Whether to overwrite data in b (may improve performance)

    Returns
    -------
    w : double or complex array, shape (M,)
        The eigenvalues, each repeated according to its multiplicity.

    (if left == True)
    vl : double or complex array, shape (M, M)
        The normalized left eigenvector corresponding to the eigenvalue w[i]
        is the column v[:,i].

    (if right == True)
    vr : double or complex array, shape (M, M)
        The normalized right eigenvector corresponding to the eigenvalue w[i]
        is the column vr[:,i].

    Raises LinAlgError if eigenvalue computation does not converge

    See Also
    --------
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays

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

def eigh(a, lower=True, eigvals_only=False, overwrite_a=False):
    """Solve the eigenvalue problem for a Hermitian or real symmetric matrix.

    Find eigenvalues w and optionally right eigenvectors v of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    Parameters
    ----------
    a : array, shape (M, M)
        A complex Hermitian or symmetric real matrix whose eigenvalues
        and eigenvectors will be computed.
    lower : boolean
        Whether the pertinent array data is taken from the lower or upper
        triangle of a. (Default: lower)
    eigvals_only : boolean
        Whether to calculate only eigenvalues and no eigenvectors.
        (Default: both are calculated)
    overwrite_a : boolean
        Whether data in a is overwritten (may improve performance).

    Returns
    -------
    w : double array, shape (M,)
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.

    (if eigvals_only == False)
    v : double or complex double array, shape (M, M)
        The normalized eigenvector corresponding to the eigenvalue w[i] is
        the column v[:,i].

    Raises LinAlgError if eigenvalue computation does not converge

    See Also
    --------
    eig : eigenvalues and right eigenvectors for non-symmetric arrays

    """
    if eigvals_only or overwrite_a:
        a1 = asarray_chkfinite(a)
        overwrite_a = overwrite_a or (_datanotshared(a1,a))
    else:
        a1 = array(a)
        if issubclass(a1.dtype.type, inexact) and not isfinite(a1).all():
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



def eig_banded(a_band, lower=0, eigvals_only=0, overwrite_a_band=0,
               select='a', select_range=None, max_ev = 0):
    """Solve real symmetric or complex hermetian band matrix eigenvalue problem.

    Find eigenvalues w and optionally right eigenvectors v of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

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
    a_band : array, shape (M, u+1)
        Banded matrix whose eigenvalues to calculate
    lower : boolean
        Is the matrix in the lower form. (Default is upper form)
    eigvals_only : boolean
        Compute only the eigenvalues and no eigenvectors.
        (Default: calculate also eigenvectors)
    overwrite_a_band:
        Discard data in a_band (may enhance performance)
    select: {'a', 'v', 'i'}
        Which eigenvalues to calculate

        ======  ========================================
        select  calculated
        ======  ========================================
        'a'     All eigenvalues
        'v'     Eigenvalues in the interval (min, max]
        'i'     Eigenvalues with indices min <= i <= max
        ======  ========================================
    select_range : (min, max)
        Range of selected eigenvalues
    max_ev : integer
        For select=='v', maximum number of eigenvalues expected.
        For other values of select, has no meaning.

        In doubt, leave this parameter untouched.

    Returns
    -------
    w : array, shape (M,)
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.

    v : double or complex double array, shape (M, M)
        The normalized eigenvector corresponding to the eigenvalue w[i] is
        the column v[:,i].

    Raises LinAlgError if eigenvalue computation does not converge

    """
    if eigvals_only or overwrite_a_band:
        a1 = asarray_chkfinite(a_band)
        overwrite_a_band = overwrite_a_band or (_datanotshared(a1,a_band))
    else:
        a1 = array(a_band)
        if issubclass(a1.dtype.type, inexact) and not isfinite(a1).all():
            raise ValueError, "array must not contain infs or NaNs"
        overwrite_a_band = 1

    if len(a1.shape) != 2:
        raise ValueError, 'expected two-dimensional array'
    if select.lower() not in [0, 1, 2, 'a', 'v', 'i', 'all', 'value', 'index']:
        raise ValueError, 'invalid argument for select'
    if select.lower() in [0, 'a', 'all']:
        if a1.dtype.char in 'GFD':
            bevd, = get_lapack_funcs(('hbevd',),(a1,))
            # FIXME: implement this somewhen, for now go with builtin values
            # FIXME: calc optimal lwork by calling ?hbevd(lwork=-1)
            #        or by using calc_lwork.f ???
            # lwork = calc_lwork.hbevd(bevd.prefix, a1.shape[0], lower)
            internal_name = 'hbevd'
        else: # a1.dtype.char in 'fd':
            bevd, = get_lapack_funcs(('sbevd',),(a1,))
            # FIXME: implement this somewhen, for now go with builtin values
            #         see above
            # lwork = calc_lwork.sbevd(bevd.prefix, a1.shape[0], lower)
            internal_name = 'sbevd'
        w,v,info = bevd(a1, compute_v = not eigvals_only,
                        lower = lower,
                        overwrite_ab = overwrite_a_band)
    if select.lower() in [1, 2, 'i', 'v', 'index', 'value']:
        # calculate certain range only
        if select.lower() in [2, 'i', 'index']:
            select = 2
            vl, vu, il, iu = 0.0, 0.0, min(select_range), max(select_range)
            if min(il, iu) < 0 or max(il, iu) >= a1.shape[1]:
                raise ValueError, 'select_range out of bounds'
            max_ev = iu - il + 1
        else:  # 1, 'v', 'value'
            select = 1
            vl, vu, il, iu = min(select_range), max(select_range), 0, 0
            if max_ev == 0:
                max_ev = a_band.shape[1]
        if eigvals_only:
            max_ev = 1
        # calculate optimal abstol for dsbevx (see manpage)
        if a1.dtype.char in 'fF':  # single precision
            lamch, = get_lapack_funcs(('lamch',),(array(0, dtype='f'),))
        else:
            lamch, = get_lapack_funcs(('lamch',),(array(0, dtype='d'),))
        abstol = 2 * lamch('s')
        if a1.dtype.char in 'GFD':
            bevx, = get_lapack_funcs(('hbevx',),(a1,))
            internal_name = 'hbevx'
        else: # a1.dtype.char in 'gfd'
            bevx, = get_lapack_funcs(('sbevx',),(a1,))
            internal_name = 'sbevx'
        # il+1, iu+1: translate python indexing (0 ... N-1) into Fortran
        # indexing (1 ... N)
        w, v, m, ifail, info = bevx(a1, vl, vu, il+1, iu+1,
                                    compute_v = not eigvals_only,
                                    mmax = max_ev,
                                    range = select, lower = lower,
                                    overwrite_ab = overwrite_a_band,
                                    abstol=abstol)
        # crop off w and v
        w = w[:m]
        if not eigvals_only:
            v = v[:, :m]
    if info<0: raise ValueError,\
    'illegal value in %-th argument of internal %s'%(-info, internal_name)
    if info>0: raise LinAlgError,"eig algorithm did not converge"

    if eigvals_only:
        return w
    return w, v

def eigvals(a,b=None,overwrite_a=0):
    """Compute eigenvalues from an ordinary or generalized eigenvalue problem.

    Find eigenvalues of a general matrix::

        a   vr[:,i] = w[i]        b   vr[:,i]

    Parameters
    ----------
    a : array, shape (M, M)
        A complex or real matrix whose eigenvalues and eigenvectors
        will be computed.
    b : array, shape (M, M)
        Right-hand side matrix in a generalized eigenvalue problem.
        If omitted, identity matrix is assumed.
    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)

    Returns
    -------
    w : double or complex array, shape (M,)
        The eigenvalues, each repeated according to its multiplicity,
        but not in any specific order.

    Raises LinAlgError if eigenvalue computation does not converge

    See Also
    --------
    eigvalsh : eigenvalues of symmetric or Hemitiean arrays
    eig : eigenvalues and right eigenvectors of general arrays
    eigh : eigenvalues and eigenvectors of symmetric/Hermitean arrays.

    """
    return eig(a,b=b,left=0,right=0,overwrite_a=overwrite_a)

def eigvalsh(a,lower=1,overwrite_a=0):
    """Solve the eigenvalue problem for a Hermitian or real symmetric matrix.

    Find eigenvalues w of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    Parameters
    ----------
    a : array, shape (M, M)
        A complex Hermitian or symmetric real matrix whose eigenvalues
        and eigenvectors will be computed.
    lower : boolean
        Whether the pertinent array data is taken from the lower or upper
        triangle of a. (Default: lower)
    overwrite_a : boolean
        Whether data in a is overwritten (may improve performance).

    Returns
    -------
    w : double array, shape (M,)
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.

    Raises LinAlgError if eigenvalue computation does not converge

    See Also
    --------
    eigvals : eigenvalues of general arrays
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eig : eigenvalues and right eigenvectors for non-symmetric arrays

    """
    return eigh(a,lower=lower,eigvals_only=1,overwrite_a=overwrite_a)

def eigvals_banded(a_band,lower=0,overwrite_a_band=0,
                   select='a', select_range=None):
    """Solve real symmetric or complex hermetian band matrix eigenvalue problem.

    Find eigenvalues w of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

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
    a_band : array, shape (M, u+1)
        Banded matrix whose eigenvalues to calculate
    lower : boolean
        Is the matrix in the lower form. (Default is upper form)
    overwrite_a_band:
        Discard data in a_band (may enhance performance)
    select: {'a', 'v', 'i'}
        Which eigenvalues to calculate

        ======  ========================================
        select  calculated
        ======  ========================================
        'a'     All eigenvalues
        'v'     Eigenvalues in the interval (min, max]
        'i'     Eigenvalues with indices min <= i <= max
        ======  ========================================
    select_range : (min, max)
        Range of selected eigenvalues

    Returns
    -------
    w : array, shape (M,)
        The eigenvalues, in ascending order, each repeated according to its
        multiplicity.

    Raises LinAlgError if eigenvalue computation does not converge

    See Also
    --------
    eig_banded : eigenvalues and right eigenvectors for symmetric/Hermitian band matrices
    eigvals : eigenvalues of general arrays
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eig : eigenvalues and right eigenvectors for non-symmetric arrays

    """
    return eig_banded(a_band,lower=lower,eigvals_only=1,
                      overwrite_a_band=overwrite_a_band, select=select,
                      select_range=select_range)

def lu_factor(a, overwrite_a=0):
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
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    getrf, = get_lapack_funcs(('getrf',),(a1,))
    lu, piv, info = getrf(a,overwrite_a=overwrite_a)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal getrf (lu_factor)'%(-info)
    if info>0: warn("Diagonal number %d is exactly zero. Singular matrix." % info,
                    RuntimeWarning)
    return lu, piv

def lu_solve(a_lu_pivots,b):
    """Solve an equation system, a x = b, given the LU factorization of a

    Parameters
    ----------
    (lu, piv)
        Factorization of the coefficient matrix a, as given by lu_factor
    b : array
        Right-hand side

    Returns
    -------
    x : array
        Solution to the system

    See also
    --------
    lu_factor : LU factorize a matrix

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
        raise ValueError, 'expected matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    flu, = get_flinalg_funcs(('lu',),(a1,))
    p,l,u,info = flu(a1,permute_l=permute_l,overwrite_a = overwrite_a)
    if info<0: raise ValueError,\
       'illegal value in %-th argument of internal lu.getrf'%(-info)
    if permute_l:
        return l,u
    return p,l,u

def svd(a,full_matrices=1,compute_uv=1,overwrite_a=0):
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
    # A hack until full_matrices == 0 support is fixed here.
    if full_matrices == 0:
        import numpy.linalg
        return numpy.linalg.svd(a, full_matrices=0, compute_uv=compute_uv)
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
    return svd(a,compute_uv=0,overwrite_a=overwrite_a)

def diagsvd(s,M,N):
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
        return r_['-1',part,zeros((M,N-M),typ)]
    elif MorN == N:
        return r_[part,zeros((M-N,N),typ)]
    else:
        raise ValueError, "Length of s must be M or N."

def cholesky(a,lower=0,overwrite_a=0):
    """Compute the Cholesky decomposition of a matrix.

    Returns the Cholesky decomposition, :lm:`A = L L^*` or :lm:`A = U^* U`
    of a Hermitian positive-definite matrix :lm:`A`.

    Parameters
    ----------
    a : array, shape (M, M)
        Matrix to be decomposed
    lower : boolean
        Whether to compute the upper or lower triangular Cholesky factorization
        (Default: upper-triangular)
    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)

    Returns
    -------
    B : array, shape (M, M)
        Upper- or lower-triangular Cholesky factor of A

    Raises LinAlgError if decomposition fails

    Examples
    --------
    >>> from scipy import array, linalg, dot
    >>> a = array([[1,-2j],[2j,5]])
    >>> L = linalg.cholesky(a, lower=True)
    >>> L
    array([[ 1.+0.j,  0.+0.j],
           [ 0.+2.j,  1.+0.j]])
    >>> dot(L, L.T.conj())
    array([[ 1.+0.j,  0.-2.j],
           [ 0.+2.j,  5.+0.j]])

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
    """Compute the Cholesky decomposition of a matrix, to use in cho_solve

    Returns a matrix containing the Cholesky decomposition,
    ``A = L L*`` or ``A = U* U`` of a Hermitian positive-definite matrix `a`.
    The return value can be directly used as the first parameter to cho_solve.

    .. warning::
        The returned matrix also contains random data in the entries not
        used by the Cholesky decomposition. If you need to zero these
        entries, use the function `cholesky` instead.

    Parameters
    ----------
    a : array, shape (M, M)
        Matrix to be decomposed
    lower : boolean
        Whether to compute the upper or lower triangular Cholesky factorization
        (Default: upper-triangular)
    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)

    Returns
    -------
    c : array, shape (M, M)
        Matrix whose upper or lower triangle contains the Cholesky factor
        of `a`. Other parts of the matrix contain random data.
    lower : boolean
        Flag indicating whether the factor is in the lower or upper triangle

    Raises
    ------
    LinAlgError
        Raised if decomposition fails.

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

    The equation system is

        A x = b,  A = U^H U = L L^H

    and A is real symmetric or complex Hermitian.

    Parameters
    ----------
    clow : tuple (c, lower)
        Cholesky factor and a flag indicating whether it is lower triangular.
        The return value from cho_factor can be used.
    b : array
        Right-hand side of the equation system

    First input is a tuple (LorU, lower) which is the output to cho_factor.
    Second input is the right-hand side.

    Returns
    -------
    x : array
        Solution to the equation system

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

def qr(a,overwrite_a=0,lwork=None,econ=False,mode='qr'):
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
        of Q and R (M, K) and (K, N) instead of (M,M) and (M,N). K=min(M,N)
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
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2:
        raise ValueError("expected 2D array")
    M, N = a1.shape
    overwrite_a = overwrite_a or (_datanotshared(a1,a))

    geqrf, = get_lapack_funcs(('geqrf',),(a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        qr,tau,work,info = geqrf(a1,lwork=-1,overwrite_a=1)
        lwork = work[0]

    qr,tau,work,info = geqrf(a1,lwork=lwork,overwrite_a=overwrite_a)
    if info<0:
        raise ValueError("illegal value in %-th argument of internal geqrf"
            % -info)

    if not econ or M<N:
        R = basic.triu(qr)
    else:
        R = basic.triu(qr[0:N,0:N])

    if mode=='r':
        return R

    if find_best_lapack_type((a1,))[0]=='s' or find_best_lapack_type((a1,))[0]=='d':
        gor_un_gqr, = get_lapack_funcs(('orgqr',),(qr,))
    else:
        gor_un_gqr, = get_lapack_funcs(('ungqr',),(qr,))


    if M<N:
        # get optimal work array
        Q,work,info = gor_un_gqr(qr[:,0:M],tau,lwork=-1,overwrite_a=1)
        lwork = work[0]
        Q,work,info = gor_un_gqr(qr[:,0:M],tau,lwork=lwork,overwrite_a=1)
    elif econ:
        # get optimal work array
        Q,work,info = gor_un_gqr(qr,tau,lwork=-1,overwrite_a=1)
        lwork = work[0]
        Q,work,info = gor_un_gqr(qr,tau,lwork=lwork,overwrite_a=1)
    else:
        t = qr.dtype.char
        qqr = numpy.empty((M,M),dtype=t)
        qqr[:,0:N]=qr
        # get optimal work array
        Q,work,info = gor_un_gqr(qqr,tau,lwork=-1,overwrite_a=1)
        lwork = work[0]
        Q,work,info = gor_un_gqr(qqr,tau,lwork=lwork,overwrite_a=1)

    if info < 0:
        raise ValueError("illegal value in %-th argument of internal gorgqr"
            % -info)

    return Q, R



def qr_old(a,overwrite_a=0,lwork=None):
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



def rq(a,overwrite_a=0,lwork=None):
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
        raise ValueError, 'expected matrix'
    M,N = a1.shape
    if M != N:
        raise ValueError, 'expected square matrix'
    if issubclass(a1.dtype.type,complexfloating):
        raise ValueError, 'expected real (non-complex) matrix'
    overwrite_a = overwrite_a or (_datanotshared(a1,a))
    gerqf, = get_lapack_funcs(('gerqf',),(a1,))
    if lwork is None or lwork == -1:
        # get optimal work array
        rq,tau,work,info = gerqf(a1,lwork=-1,overwrite_a=1)
        lwork = work[0]
    rq,tau,work,info = gerqf(a1,lwork=lwork,overwrite_a=overwrite_a)
    if info<0: raise ValueError, \
       'illegal value in %-th argument of internal geqrf'%(-info)
    gemm, = get_blas_funcs(('gemm',),(rq,))
    t = rq.dtype.char
    R = basic.triu(rq)
    Q = numpy.identity(M,dtype=t)
    ident = numpy.identity(M,dtype=t)
    zeros = numpy.zeros

    k = min(M,N)
    for i in range(k):
        v = zeros((M,),t)
        v[N-k+i] = 1
        v[0:N-k+i] = rq[M-k+i,0:N-k+i]
        H = gemm(-tau[i],v,v,1+0j,ident,trans_b=2)
        Q = gemm(1,Q,H)
    return R, Q

_double_precision = ['i','l','d']

def schur(a,output='real',lwork=None,overwrite_a=0):
    """Compute Schur decomposition of a matrix.

    The Schur decomposition is

        A = Z T Z^H

    where Z is unitary and T is either upper-triangular, or for real
    Schur decomposition (output='real'), quasi-upper triangular.  In
    the quasi-triangular form, 2x2 blocks describing complex-valued
    eigenvalue pairs may extrude from the diagonal.

    Parameters
    ----------
    a : array, shape (M, M)
        Matrix to decompose
    output : {'real', 'complex'}
        Construct the real or complex Schur decomposition (for real matrices).
    lwork : integer
        Work array size. If None or -1, it is automatically computed.
    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)

    Returns
    -------
    T : array, shape (M, M)
        Schur form of A. It is real-valued for the real Schur decomposition.
    Z : array, shape (M, M)
        An unitary Schur transformation matrix for A.
        It is real-valued for the real Schur decomposition.

    See also
    --------
    rsf2csf : Convert real Schur form to complex Schur form

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

eps = numpy.finfo(float).eps
feps = numpy.finfo(single).eps

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
    """Convert real Schur form to complex Schur form.

    Convert a quasi-diagonal real-valued Schur form to the upper triangular
    complex-valued Schur form.

    Parameters
    ----------
    T : array, shape (M, M)
        Real Schur form of the original matrix
    Z : array, shape (M, M)
        Schur transformation matrix

    Returns
    -------
    T : array, shape (M, M)
        Complex Schur form of the original matrix
    Z : array, shape (M, M)
        Schur transformation matrix corresponding to the complex form

    See also
    --------
    schur : Schur decompose a matrix

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
    u,s,vh = svd(A)
    M,N = A.shape
    tol = max(M,N)*numpy.amax(s)*eps
    num = numpy.sum(s > tol,dtype=int)
    Q = u[:,:num]
    return Q

def hessenberg(a,calc_q=0,overwrite_a=0):
    """Compute Hessenberg form of a matrix.

    The Hessenberg decomposition is

        A = Q H Q^H

    where Q is unitary/orthogonal and H has only zero elements below the first
    subdiagonal.

    Parameters
    ----------
    a : array, shape (M,M)
        Matrix to bring into Hessenberg form
    calc_q : boolean
        Whether to compute the transformation matrix
    overwrite_a : boolean
        Whether to ovewrite data in a (may improve performance)

    Returns
    -------
    H : array, shape (M,M)
        Hessenberg form of A

    (If calc_q == True)
    Q : array, shape (M,M)
        Unitary/orthogonal similarity transformation matrix s.t. A = Q H Q^H

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
