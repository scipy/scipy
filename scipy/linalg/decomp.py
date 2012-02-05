#
# Author: Pearu Peterson, March 2002
#
# additions by Travis Oliphant, March 2002
# additions by Eric Jones,      June 2002
# additions by Johannes Loehnert, June 2006
# additions by Bart Vandereycken, June 2006
# additions by Andrew D Straw, May 2007
# additions by Tiziano Zito, November 2008
#
# April 2010: Functions for LU, QR, SVD, Schur and Cholesky decompositions were
# moved to their own files.  Still in this file are functions for eigenstuff
# and for the Hessenberg form.

__all__ = ['eig','eigh','eig_banded','eigvals','eigvalsh', 'eigvals_banded',
           'hessenberg']

import numpy
from numpy import array, asarray_chkfinite, asarray, diag, zeros, ones, \
        isfinite, inexact, nonzero, iscomplexobj, cast, flatnonzero, conj

# Local imports
from scipy.linalg import calc_lwork
from misc import LinAlgError, _datacopied
from lapack import get_lapack_funcs
from blas import get_blas_funcs


_I = cast['F'](1j)

def _make_complex_eigvecs(w, vin, dtype):
    """
    Produce complex-valued eigenvectors from LAPACK DGGEV real-valued output
    """
    # - see LAPACK man page DGGEV at ALPHAI
    v = numpy.array(vin, dtype=dtype)
    m = (w.imag > 0)
    m[:-1] |= (w.imag[1:] < 0) # workaround for LAPACK bug, cf. ticket #709
    for i in flatnonzero(m):
        v.imag[:,i] = vin[:,i+1]
        conj(v[:,i], v[:,i+1])
    return v

def _geneig(a1, b, left, right, overwrite_a, overwrite_b):
    b1 = asarray(b)
    overwrite_b = overwrite_b or _datacopied(b1, b)
    if len(b1.shape) != 2 or b1.shape[0] != b1.shape[1]:
        raise ValueError('expected square matrix')
    ggev, = get_lapack_funcs(('ggev',), (a1, b1))
    cvl, cvr = left, right
    if ggev.module_name[:7] == 'clapack':
        raise NotImplementedError('calling ggev from %s' % ggev.module_name)
    res = ggev(a1, b1, lwork=-1)
    lwork = res[-2][0].real.astype(numpy.int)
    if ggev.prefix in 'cz':
        alpha, beta, vl, vr, work, info = ggev(a1, b1, cvl, cvr, lwork,
                                                    overwrite_a, overwrite_b)
        w = alpha / beta
    else:
        alphar, alphai, beta, vl, vr, work, info = ggev(a1, b1, cvl, cvr, lwork,
                                                        overwrite_a,overwrite_b)
        w = (alphar + _I * alphai) / beta
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal ggev'
                                                                    % -info)
    if info > 0:
        raise LinAlgError("generalized eig algorithm did not converge (info=%d)"
                                                                    % info)

    only_real = numpy.logical_and.reduce(numpy.equal(w.imag, 0.0))
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

def eig(a, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False):
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
        raise ValueError('expected square matrix')
    overwrite_a = overwrite_a or (_datacopied(a1, a))
    if b is not None:
        b = asarray_chkfinite(b)
        if b.shape != a1.shape:
            raise ValueError('a and b must have the same shape')
        return _geneig(a1, b, left, right, overwrite_a, overwrite_b)
    geev, = get_lapack_funcs(('geev',), (a1,))
    compute_vl, compute_vr = left, right
    if geev.module_name[:7] == 'flapack':
        lwork = calc_lwork.geev(geev.prefix, a1.shape[0],
                                    compute_vl, compute_vr)[1]
        if geev.prefix in 'cz':
            w, vl, vr, info = geev(a1, lwork=lwork,
                                        compute_vl=compute_vl,
                                        compute_vr=compute_vr,
                                        overwrite_a=overwrite_a)
        else:
            wr, wi, vl, vr, info = geev(a1, lwork=lwork,
                                        compute_vl=compute_vl,
                                        compute_vr=compute_vr,
                                        overwrite_a=overwrite_a)
            t = {'f':'F','d':'D'}[wr.dtype.char]
            w = wr + _I * wi
    else: # 'clapack'
        if geev.prefix in 'cz':
            w, vl, vr, info = geev(a1,
                                    compute_vl=compute_vl,
                                    compute_vr=compute_vr,
                                    overwrite_a=overwrite_a)
        else:
            wr, wi, vl, vr, info = geev(a1,
                                        compute_vl=compute_vl,
                                        compute_vr=compute_vr,
                                        overwrite_a=overwrite_a)
            t = {'f':'F','d':'D'}[wr.dtype.char]
            w = wr + _I * wi
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal geev'
                                                                    % -info)
    if info > 0:
        raise LinAlgError("eig algorithm did not converge (only eigenvalues "
                            "with order >= %d have converged)" % info)

    only_real = numpy.logical_and.reduce(numpy.equal(w.imag, 0.0))
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

def eigh(a, b=None, lower=True, eigvals_only=False, overwrite_a=False,
         overwrite_b=False, turbo=True, eigvals=None, type=1):
    """Solve an ordinary or generalized eigenvalue problem for a complex
    Hermitian or real symmetric matrix.

    Find eigenvalues w and optionally eigenvectors v of matrix a, where
    b is positive definite::

                      a v[:,i] = w[i] b v[:,i]
        v[i,:].conj() a v[:,i] = w[i]
        v[i,:].conj() b v[:,i] = 1


    Parameters
    ----------
    a : array, shape (M, M)
        A complex Hermitian or real symmetric matrix whose eigenvalues and
        eigenvectors will be computed.
    b : array, shape (M, M)
        A complex Hermitian or real symmetric definite positive matrix in.
        If omitted, identity matrix is assumed.
    lower : boolean
        Whether the pertinent array data is taken from the lower or upper
        triangle of a. (Default: lower)
    eigvals_only : boolean
        Whether to calculate only eigenvalues and no eigenvectors.
        (Default: both are calculated)
    turbo : boolean
        Use divide and conquer algorithm (faster but expensive in memory,
        only for generalized eigenvalue problem and if eigvals=None)
    eigvals : tuple (lo, hi)
        Indexes of the smallest and largest (in ascending order) eigenvalues
        and corresponding eigenvectors to be returned: 0 <= lo < hi <= M-1.
        If omitted, all eigenvalues and eigenvectors are returned.
    type: integer
        Specifies the problem type to be solved:
           type = 1: a   v[:,i] = w[i] b v[:,i]
           type = 2: a b v[:,i] = w[i]   v[:,i]
           type = 3: b a v[:,i] = w[i]   v[:,i]
    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)
    overwrite_b : boolean
        Whether to overwrite data in b (may improve performance)

    Returns
    -------
    w : real array, shape (N,)
        The N (1<=N<=M) selected eigenvalues, in ascending order, each
        repeated according to its multiplicity.

    (if eigvals_only == False)
    v : complex array, shape (M, N)
        The normalized selected eigenvector corresponding to the
        eigenvalue w[i] is the column v[:,i]. Normalization:
        type 1 and 3:       v.conj() a      v  = w
        type 2:        inv(v).conj() a  inv(v) = w
        type = 1 or 2:      v.conj() b      v  = I
        type = 3     :      v.conj() inv(b) v  = I

    Raises LinAlgError if eigenvalue computation does not converge,
    an error occurred, or b matrix is not definite positive. Note that
    if input matrices are not symmetric or hermitian, no error is reported
    but results will be wrong.

    See Also
    --------
    eig : eigenvalues and right eigenvectors for non-symmetric arrays

    """
    a1 = asarray_chkfinite(a)
    if len(a1.shape) != 2 or a1.shape[0] != a1.shape[1]:
        raise ValueError('expected square matrix')
    overwrite_a = overwrite_a or (_datacopied(a1, a))
    if iscomplexobj(a1):
        cplx = True
    else:
        cplx = False
    if b is not None:
        b1 = asarray_chkfinite(b)
        overwrite_b = overwrite_b or _datacopied(b1, b)
        if len(b1.shape) != 2 or b1.shape[0] != b1.shape[1]:
            raise ValueError('expected square matrix')

        if b1.shape != a1.shape:
            raise ValueError("wrong b dimensions %s, should "
                             "be %s" % (str(b1.shape), str(a1.shape)))
        if iscomplexobj(b1):
            cplx = True
        else:
            cplx = cplx or False
    else:
        b1 = None

    # Set job for fortran routines
    _job = (eigvals_only and 'N') or 'V'

    # port eigenvalue range from python to fortran convention
    if eigvals is not None:
        lo, hi = eigvals
        if lo < 0 or hi >= a1.shape[0]:
            raise ValueError('The eigenvalue range specified is not valid.\n'
                             'Valid range is [%s,%s]' % (0, a1.shape[0]-1))
        lo += 1
        hi += 1
        eigvals = (lo, hi)

    # set lower
    if lower:
        uplo = 'L'
    else:
        uplo = 'U'

    # fix prefix for lapack routines
    if cplx:
        pfx = 'he'
    else:
        pfx = 'sy'

    #  Standard Eigenvalue Problem
    #  Use '*evr' routines
    # FIXME: implement calculation of optimal lwork
    #        for all lapack routines
    if b1 is None:
        (evr,) = get_lapack_funcs((pfx+'evr',), (a1,))
        if eigvals is None:
            w, v, info = evr(a1, uplo=uplo, jobz=_job, range="A", il=1,
                             iu=a1.shape[0], overwrite_a=overwrite_a)
        else:
            (lo, hi)= eigvals
            w_tot, v, info = evr(a1, uplo=uplo, jobz=_job, range="I",
                                 il=lo, iu=hi, overwrite_a=overwrite_a)
            w = w_tot[0:hi-lo+1]

    # Generalized Eigenvalue Problem
    else:
        # Use '*gvx' routines if range is specified
        if eigvals is not None:
            (gvx,) = get_lapack_funcs((pfx+'gvx',), (a1,b1))
            (lo, hi) = eigvals
            w_tot, v, ifail, info = gvx(a1, b1, uplo=uplo, iu=hi,
                                        itype=type,jobz=_job, il=lo,
                                        overwrite_a=overwrite_a,
                                        overwrite_b=overwrite_b)
            w = w_tot[0:hi-lo+1]
        # Use '*gvd' routine if turbo is on and no eigvals are specified
        elif turbo:
            (gvd,) = get_lapack_funcs((pfx+'gvd',), (a1,b1))
            v, w, info = gvd(a1, b1, uplo=uplo, itype=type, jobz=_job,
                             overwrite_a=overwrite_a,
                             overwrite_b=overwrite_b)
        # Use '*gv' routine if turbo is off and no eigvals are specified
        else:
            (gv,) = get_lapack_funcs((pfx+'gv',), (a1,b1))
            v, w, info = gv(a1, b1, uplo=uplo, itype= type, jobz=_job,
                            overwrite_a=overwrite_a,
                            overwrite_b=overwrite_b)

    # Check if we had a  successful exit
    if info == 0:
        if eigvals_only:
            return w
        else:
            return w, v

    elif info < 0:
        raise LinAlgError("illegal value in %i-th argument of internal"
                          " fortran routine." % (-info))
    elif info > 0 and b1 is None:
        raise LinAlgError("unrecoverable internal error.")

    # The algorithm failed to converge.
    elif info > 0 and info <= b1.shape[0]:
        if eigvals is not None:
            raise LinAlgError("the eigenvectors %s failed to"
                              " converge." % nonzero(ifail)-1)
        else:
            raise LinAlgError("internal fortran routine failed to converge: "
                              "%i off-diagonal elements of an "
                              "intermediate tridiagonal form did not converge"
                              " to zero." % info)

    # This occurs when b is not positive definite
    else:
        raise LinAlgError("the leading minor of order %i"
                          " of 'b' is not positive definite. The"
                          " factorization of 'b' could not be completed"
                          " and no eigenvalues or eigenvectors were"
                          " computed." % (info-b1.shape[0]))

def eig_banded(a_band, lower=False, eigvals_only=False, overwrite_a_band=False,
               select='a', select_range=None, max_ev = 0):
    """Solve real symmetric or complex hermitian band matrix eigenvalue problem.

    Find eigenvalues w and optionally right eigenvectors v of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    The matrix a is stored in a_band either in lower diagonal or upper
    diagonal ordered form:

        a_band[u + i - j, j] == a[i,j]        (if upper form; i <= j)
        a_band[    i - j, j] == a[i,j]        (if lower form; i >= j)

    where u is the number of bands above the diagonal.

    Example of a_band (shape of a is (6,6), u=2)::

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
    a_band : array, shape (u+1, M)
        The bands of the M by M matrix a.
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
        overwrite_a_band = overwrite_a_band or (_datacopied(a1, a_band))
    else:
        a1 = array(a_band)
        if issubclass(a1.dtype.type, inexact) and not isfinite(a1).all():
            raise ValueError("array must not contain infs or NaNs")
        overwrite_a_band = 1

    if len(a1.shape) != 2:
        raise ValueError('expected two-dimensional array')
    if select.lower() not in [0, 1, 2, 'a', 'v', 'i', 'all', 'value', 'index']:
        raise ValueError('invalid argument for select')
    if select.lower() in [0, 'a', 'all']:
        if a1.dtype.char in 'GFD':
            bevd, = get_lapack_funcs(('hbevd',), (a1,))
            # FIXME: implement this somewhen, for now go with builtin values
            # FIXME: calc optimal lwork by calling ?hbevd(lwork=-1)
            #        or by using calc_lwork.f ???
            # lwork = calc_lwork.hbevd(bevd.prefix, a1.shape[0], lower)
            internal_name = 'hbevd'
        else: # a1.dtype.char in 'fd':
            bevd, = get_lapack_funcs(('sbevd',), (a1,))
            # FIXME: implement this somewhen, for now go with builtin values
            #         see above
            # lwork = calc_lwork.sbevd(bevd.prefix, a1.shape[0], lower)
            internal_name = 'sbevd'
        w,v,info = bevd(a1, compute_v=not eigvals_only,
                        lower=lower,
                        overwrite_ab=overwrite_a_band)
    if select.lower() in [1, 2, 'i', 'v', 'index', 'value']:
        # calculate certain range only
        if select.lower() in [2, 'i', 'index']:
            select = 2
            vl, vu, il, iu = 0.0, 0.0, min(select_range), max(select_range)
            if min(il, iu) < 0 or max(il, iu) >= a1.shape[1]:
                raise ValueError('select_range out of bounds')
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
            lamch, = get_lapack_funcs(('lamch',), (array(0, dtype='f'),))
        else:
            lamch, = get_lapack_funcs(('lamch',), (array(0, dtype='d'),))
        abstol = 2 * lamch('s')
        if a1.dtype.char in 'GFD':
            bevx, = get_lapack_funcs(('hbevx',), (a1,))
            internal_name = 'hbevx'
        else: # a1.dtype.char in 'gfd'
            bevx, = get_lapack_funcs(('sbevx',), (a1,))
            internal_name = 'sbevx'
        # il+1, iu+1: translate python indexing (0 ... N-1) into Fortran
        # indexing (1 ... N)
        w, v, m, ifail, info = bevx(a1, vl, vu, il+1, iu+1,
                                    compute_v=not eigvals_only,
                                    mmax=max_ev,
                                    range=select, lower=lower,
                                    overwrite_ab=overwrite_a_band,
                                    abstol=abstol)
        # crop off w and v
        w = w[:m]
        if not eigvals_only:
            v = v[:, :m]
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal %s'
                                                    % (-info, internal_name))
    if info > 0:
        raise LinAlgError("eig algorithm did not converge")

    if eigvals_only:
        return w
    return w, v

def eigvals(a, b=None, overwrite_a=False):
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
    return eig(a, b=b, left=0, right=0, overwrite_a=overwrite_a)

def eigvalsh(a, b=None, lower=True, overwrite_a=False,
             overwrite_b=False, turbo=True, eigvals=None, type=1):
    """Solve an ordinary or generalized eigenvalue problem for a complex
    Hermitian or real symmetric matrix.

    Find eigenvalues w of matrix a, where b is positive definite::

                      a v[:,i] = w[i] b v[:,i]
        v[i,:].conj() a v[:,i] = w[i]
        v[i,:].conj() b v[:,i] = 1


    Parameters
    ----------
    a : array, shape (M, M)
        A complex Hermitian or real symmetric matrix whose eigenvalues and
        eigenvectors will be computed.
    b : array, shape (M, M)
        A complex Hermitian or real symmetric definite positive matrix in.
        If omitted, identity matrix is assumed.
    lower : boolean
        Whether the pertinent array data is taken from the lower or upper
        triangle of a. (Default: lower)
    turbo : boolean
        Use divide and conquer algorithm (faster but expensive in memory,
        only for generalized eigenvalue problem and if eigvals=None)
    eigvals : tuple (lo, hi)
        Indexes of the smallest and largest (in ascending order) eigenvalues
        and corresponding eigenvectors to be returned: 0 <= lo < hi <= M-1.
        If omitted, all eigenvalues and eigenvectors are returned.
    type: integer
        Specifies the problem type to be solved:
           type = 1: a   v[:,i] = w[i] b v[:,i]
           type = 2: a b v[:,i] = w[i]   v[:,i]
           type = 3: b a v[:,i] = w[i]   v[:,i]
    overwrite_a : boolean
        Whether to overwrite data in a (may improve performance)
    overwrite_b : boolean
        Whether to overwrite data in b (may improve performance)

    Returns
    -------
    w : real array, shape (N,)
        The N (1<=N<=M) selected eigenvalues, in ascending order, each
        repeated according to its multiplicity.

    Raises LinAlgError if eigenvalue computation does not converge,
    an error occurred, or b matrix is not definite positive. Note that
    if input matrices are not symmetric or hermitian, no error is reported
    but results will be wrong.

    See Also
    --------
    eigvals : eigenvalues of general arrays
    eigh : eigenvalues and right eigenvectors for symmetric/Hermitian arrays
    eig : eigenvalues and right eigenvectors for non-symmetric arrays

    """
    return eigh(a, b=b, lower=lower, eigvals_only=True,
                overwrite_a=overwrite_a, overwrite_b=overwrite_b,
                turbo=turbo, eigvals=eigvals, type=type)

def eigvals_banded(a_band, lower=False, overwrite_a_band=False,
                   select='a', select_range=None):
    """Solve real symmetric or complex hermitian band matrix eigenvalue problem.

    Find eigenvalues w of a::

        a v[:,i] = w[i] v[:,i]
        v.H v    = identity

    The matrix a is stored in a_band either in lower diagonal or upper
    diagonal ordered form:

        a_band[u + i - j, j] == a[i,j]        (if upper form; i <= j)
        a_band[    i - j, j] == a[i,j]        (if lower form; i >= j)

    where u is the number of bands above the diagonal.

    Example of a_band (shape of a is (6,6), u=2)::

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
    a_band : array, shape (u+1, M)
        The bands of the M by M matrix a.
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
    return eig_banded(a_band, lower=lower, eigvals_only=1,
                      overwrite_a_band=overwrite_a_band, select=select,
                      select_range=select_range)

_double_precision = ['i','l','d']


def hessenberg(a, calc_q=False, overwrite_a=False):
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
        raise ValueError('expected square matrix')
    overwrite_a = overwrite_a or (_datacopied(a1, a))
    gehrd,gebal = get_lapack_funcs(('gehrd','gebal'), (a1,))
    ba, lo, hi, pivscale, info = gebal(a1, permute=1, overwrite_a=overwrite_a)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal gebal '
                                                    '(hessenberg)' % -info)
    n = len(a1)
    lwork = calc_lwork.gehrd(gehrd.prefix, n, lo, hi)
    hq, tau, info = gehrd(ba, lo=lo, hi=hi, lwork=lwork, overwrite_a=1)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal gehrd '
                                        '(hessenberg)' % -info)

    if not calc_q:
        for i in range(lo, hi):
            hq[i+2:hi+1, i] = 0.0
        return hq

    # XXX: Use ORGHR routines to compute q.
    typecode = hq.dtype
    ger,gemm = get_blas_funcs(('ger','gemm'), dtype=typecode)
    q = None
    for i in range(lo, hi):
        if tau[i]==0.0:
            continue
        v = zeros(n, dtype=typecode)
        v[i+1] = 1.0
        v[i+2:hi+1] = hq[i+2:hi+1, i]
        hq[i+2:hi+1, i] = 0.0
        h = ger(-tau[i], v, v,a=diag(ones(n, dtype=typecode)), overwrite_a=1)
        if q is None:
            q = h
        else:
            q = gemm(1.0, q, h)
    if q is None:
        q = diag(ones(n, dtype=typecode))
    return hq, q
