"""
Find a few eigenvectors and eigenvalues of a matrix.


Uses ARPACK: http://www.caam.rice.edu/software/ARPACK/

"""
# Wrapper implementation notes
#
# ARPACK Entry Points
# -------------------
# The entry points to ARPACK are
# - (s,d)seupd : single and double precision symmetric matrix
# - (s,d,c,z)neupd: single,double,complex,double complex general matrix
# This wrapper puts the *neupd (general matrix) interfaces in eigen()
# and the *seupd (symmetric matrix) in eigen_symmetric().
# There is no Hermetian complex/double complex interface.
# To find eigenvalues of a Hermetian matrix you
# must use eigen() and not eigen_symmetric()
# It might be desirable to handle the Hermetian case differently
# and, for example, return real eigenvalues.

# Number of eigenvalues returned and complex eigenvalues
# ------------------------------------------------------
# The ARPACK nonsymmetric real and double interface (s,d)naupd return
# eigenvalues and eigenvectors in real (float,double) arrays.
# Since the eigenvalues and eigenvectors are, in general, complex
# ARPACK puts the real and imaginary parts in consecutive entries
# in real-valued arrays.   This wrapper puts the real entries
# into complex data types and attempts to return the requested eigenvalues
# and eigenvectors.


# Solver modes
# ------------
# ARPACK and handle shifted and shift-inverse computations
# for eigenvalues by providing a shift (sigma) and a solver.
# This is currently not implemented

__docformat__ = "restructuredtext en"

__all___=['eigen','eigen_symmetric', 'svd']

import warnings

import _arpack
import numpy as np
from scipy.sparse.linalg.interface import aslinearoperator
from scipy.sparse import csc_matrix, csr_matrix

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}
_ndigits = {'f':5, 'd':12, 'F':5, 'D':12}

class _ArpackParams(object):
    def __init__(self, n, k, tp, mode="symmetric", sigma=None,
                 ncv=None, v0=None, maxiter=None, which="LM", tol=0):
        if k <= 0:
            raise ValueError("k must be positive, k=%d" % k)
        if k == n:
            raise ValueError("k must be less than rank(A), k=%d" % k)

        if maxiter is None:
            maxiter = n * 10
        if maxiter <= 0:
            raise ValueError("maxiter must be positive, maxiter=%d" % maxiter)

        if tp not in 'fdFD':
            raise ValueError("matrix type must be 'f', 'd', 'F', or 'D'")

        if v0 is not None:
            self.resid = v0
            info = 1
        else:
            self.resid = np.zeros(n, tp)
            info = 0

        if sigma is not None:
            raise NotImplementedError("shifted eigenproblem not supported yet")

        if ncv is None:
            ncv = 2 * k + 1
        ncv = min(ncv, n)

        if ncv > n or ncv < k:
            raise ValueError("ncv must be k<=ncv<=n, ncv=%s" % ncv)

        if not which in ["LM", "SM", "LR", "SR", "LI", "SI"]:
            raise ValueError("Parameter which must be one of %s" % ' '.join(whiches))

        ltr = _type_conv[tp]

        self.v = np.zeros((n, ncv), tp) # holds Ritz vectors
        self.rwork = None # Only used for unsymmetric, complex solver

        if mode == "unsymmetric":
            self.workd = np.zeros(3 * n, tp)
            self.workl = np.zeros(3 * ncv * ncv + 6 * ncv, tp)
            self.solver = _arpack.__dict__[ltr + 'naupd']
            self.extract = _arpack.__dict__[ltr + 'neupd']

            if tp in 'FD':
                self.rwork = np.zeros(ncv, tp.lower())

            self.ipntr = np.zeros(14, "int")
        elif mode == "symmetric":
            self.workd = np.zeros(3 * n, tp)
            self.workl = np.zeros(ncv * (ncv + 8), tp)
            self.solver = _arpack.__dict__[ltr + 'saupd']
            self.extract = _arpack.__dict__[ltr + 'seupd']

            self.ipntr = np.zeros(11, "int")
        else:
            raise ValueError("Unrecognized mode %s" % mode)

        self.iparam = np.zeros(11, "int")

        # set solver mode and parameters
        # only supported mode is 1: Ax=lx
        ishfts = 1
        mode1 = 1
        self.iparam[0] = ishfts
        self.iparam[2] = maxiter
        self.iparam[6] = mode1

        self.n = n
        self.mode = mode
        self.tol = tol
        self.k = k
        self.maxiter = maxiter
        self.ncv = ncv
        self.which = which
        self.tp = tp
        self.info = info
        self.bmat = 'I'

def eigen(A, k=6, M=None, sigma=None, which='LM', v0=None,
          ncv=None, maxiter=None, tol=0,
          return_eigenvectors=True):
    """Find k eigenvalues and eigenvectors of the square matrix A.

    Solves A * x[i] = w[i] * x[i], the standard eigenvalue problem for
    w[i] eigenvalues with corresponding eigenvectors x[i].


    Parameters
    ----------
    A : matrix, array, or object with matvec(x) method
        An N x N matrix, array, or an object with matvec(x) method to perform
        the matrix vector product A * x.  The sparse matrix formats
        in scipy.sparse are appropriate for A.

    k : integer
        The number of eigenvalues and eigenvectors desired

    Returns
    -------
    w : array
        Array of k eigenvalues

    v : array
        An array of k eigenvectors
        The v[i] is the eigenvector corresponding to the eigenvector w[i]

    Other Parameters
    ----------------

    M : matrix or array
        (Not implemented)
        A symmetric positive-definite matrix for the generalized
        eigenvalue problem A * x = w * M * x

    sigma : real or complex
        (Not implemented)
        Find eigenvalues near sigma.  Shift spectrum by sigma.

    v0 : array
        Starting vector for iteration.

    ncv : integer
        The number of Lanczos vectors generated
        ncv must be greater than k; it is recommended that ncv > 2*k

    which : string
        Which k eigenvectors and eigenvalues to find:
         - 'LM' : largest magnitude
         - 'SM' : smallest magnitude
         - 'LR' : largest real part
         - 'SR' : smallest real part
         - 'LI' : largest imaginary part
         - 'SI' : smallest imaginary part

    maxiter : integer
        Maximum number of Arnoldi update iterations allowed

    tol : float
        Relative accuracy for eigenvalues (stopping criterion)

    return_eigenvectors : boolean
        Return eigenvectors (True) in addition to eigenvalues

    See Also
    --------
    eigen_symmetric : eigenvalues and eigenvectors for symmetric matrix A

    Notes
    -----

    Examples
    --------

    """
    A = aslinearoperator(A)
    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix (shape=%s)' % A.shape)
    n = A.shape[0]

    params = _ArpackParams(n, k, A.dtype.char, "unsymmetric", sigma,
                           ncv, v0, maxiter, which, tol)

    if M is not None:
        raise NotImplementedError("generalized eigenproblem not supported yet")

    ido = 0

    while True:
        if params.tp in 'fd':
            ido, params.resid, params.v, params.iparam, params.ipntr, params.info = \
                params.solver(ido, params.bmat, params.which, params.k, params.tol,
                        params.resid, params.v, params.iparam, params.ipntr,
                        params.workd, params.workl, params.info)
        else:
            ido, params.resid, params.v, params.iparam, params.ipntr, params.info =\
                params.solver(ido, params.bmat, params.which, params.k, params.tol,
                        params.resid, params.v, params.iparam, params.ipntr,
                        params.workd, params.workl, params.rwork, params.info)

        xslice = slice(params.ipntr[0]-1, params.ipntr[0]-1+n)
        yslice = slice(params.ipntr[1]-1, params.ipntr[1]-1+n)
        if ido == -1:
            # initialization
            params.workd[yslice] = A.matvec(params.workd[xslice])
        elif ido == 1:
            # compute y=Ax
            params.workd[yslice] = A.matvec(params.workd[xslice])
        else:
            break

    if params.info < -1 :
        raise RuntimeError("Error info=%d in arpack" % params.info)
    elif params.info == -1:
        warnings.warn("Maximum number of iterations taken: %s" % self.iparam[2])

    # now extract eigenvalues and (optionally) eigenvectors
    rvec = return_eigenvectors
    ierr = 0
    howmny = 'A' # return all eigenvectors
    sselect = np.zeros(params.ncv, 'int') # unused
    sigmai = 0.0 # no shifts, not implemented
    sigmar = 0.0 # no shifts, not implemented
    workev = np.zeros(3 * params.ncv, params.tp)

    if params.tp in 'fd':
        dr = np.zeros(k+1, params.tp)
        di = np.zeros(k+1, params.tp)
        zr = np.zeros((n, k+1), params.tp)
        dr, di, zr, params.info=\
            params.extract(rvec, howmny, sselect, sigmar, sigmai, workev,
                   params.bmat, params.which, k, params.tol, params.resid,
                   params.v, params.iparam, params.ipntr,
                   params.workd, params.workl, params.info)

        # The ARPACK nonsymmetric real and double interface (s,d)naupd return
        # eigenvalues and eigenvectors in real (float,double) arrays.

        # Build complex eigenvalues from real and imaginary parts
        d = dr + 1.0j * di

        # Arrange the eigenvectors: complex eigenvectors are stored as
        # real,imaginary in consecutive columns
        z = zr.astype(params.tp.upper())
        eps = np.finfo(params.tp).eps
        i = 0
        while i<=k:
            # check if complex
            if abs(d[i].imag) > eps:
                # assume this is a complex conjugate pair with eigenvalues
                # in consecutive columns
                z[:,i] = zr[:,i] + 1.0j * zr[:,i+1]
                z[:,i+1] = z[:,i].conjugate()
                i +=1
            i += 1

        # Now we have k+1 possible eigenvalues and eigenvectors
        # Return the ones specified by the keyword "which"
        nreturned = params.iparam[4] # number of good eigenvalues returned
        if nreturned == k:    # we got exactly how many eigenvalues we wanted
            d = d[:k]
            z = z[:,:k]
        else:   # we got one extra eigenvalue (likely a cc pair, but which?)
            # cut at approx precision for sorting
            rd = np.round(d, decimals = _ndigits[params.tp])
            if params.which in ['LR','SR']:
                ind = np.argsort(rd.real)
            elif which in ['LI','SI']:
                # for LI,SI ARPACK returns largest,smallest abs(imaginary) why?
                ind = np.argsort(abs(rd.imag))
            else:
                ind = np.argsort(abs(rd))
            if params.which in ['LR','LM','LI']:
                d = d[ind[-k:]]
                z = z[:,ind[-k:]]
            if params.which in ['SR','SM','SI']:
                d = d[ind[:k]]
                z = z[:,ind[:k]]


    else:
        # complex is so much simpler...
        d, z, params.info =\
                params.extract(rvec, howmny, sselect, sigmar, workev,
                       params.bmat, params.which, k, params.tol, params.resid,
                       params.v, params.iparam, params.ipntr,
                       params.workd, params.workl, params.rwork, ierr)



    if ierr != 0:
        raise RuntimeError("Error info=%d in arpack"%info)
        return None
    if return_eigenvectors:
        return d,z
    return d


def eigen_symmetric(A, k=6, M=None, sigma=None, which='LM', v0=None,
                    ncv=None, maxiter=None, tol=0,
                    return_eigenvectors=True):
    """Find k eigenvalues and eigenvectors of the real symmetric
    square matrix A.

    Solves A * x[i] = w[i] * x[i], the standard eigenvalue problem for
    w[i] eigenvalues with corresponding eigenvectors x[i].


    Parameters
    ----------
    A : matrix or array with real entries or object with matvec(x) method
        An N x N real symmetric matrix or array or an object with matvec(x)
        method to perform the matrix vector product A * x.  The sparse
        matrix formats in scipy.sparse are appropriate for A.

    k : integer
        The number of eigenvalues and eigenvectors desired

    Returns
    -------
    w : array
        Array of k eigenvalues

    v : array
       An array of k eigenvectors
       The v[i] is the eigenvector corresponding to the eigenvector w[i]

    Other Parameters
    ----------------
    M : matrix or array
        (Not implemented)
        A symmetric positive-definite matrix for the generalized
        eigenvalue problem A * x = w * M * x


    sigma : real
        (Not implemented)
        Find eigenvalues near sigma.  Shift spectrum by sigma.

    v0 : array
        Starting vector for iteration.

    ncv : integer
        The number of Lanczos vectors generated
        ncv must be greater than k; it is recommended that ncv > 2*k

    which : string
        Which k eigenvectors and eigenvalues to find:
         - 'LA' : Largest (algebraic) eigenvalues
         - 'SA' : Smallest (algebraic) eigenvalues
         - 'LM' : Largest (in magnitude) eigenvalues
         - 'SM' : Smallest (in magnitude) eigenvalues
         - 'BE' : Half (k/2) from each end of the spectrum
                  When k is odd, return one more (k/2+1) from the high end

    maxiter : integer
        Maximum number of Arnoldi update iterations allowed

    tol : float
        Relative accuracy for eigenvalues (stopping criterion)

    return_eigenvectors : boolean
        Return eigenvectors (True) in addition to eigenvalues

    See Also
    --------
    eigen : eigenvalues and eigenvectors for a general (nonsymmetric) matrix A

    Notes
    -----

    Examples
    --------
    """
    A = aslinearoperator(A)
    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix (shape=%s)' % shape)
    n = A.shape[0]

    # guess type
    typ = A.dtype.char
    if typ not in 'fd':
        raise ValueError("matrix must be real valued (type must be 'f' or 'd')")

    if M is not None:
        raise NotImplementedError("generalized eigenproblem not supported yet")
    if sigma is not None:
        raise NotImplementedError("shifted eigenproblem not supported yet")

    if ncv is None:
        ncv=2*k+1
    ncv=min(ncv,n)
    if maxiter==None:
        maxiter=n*10
    # assign starting vector
    if v0 is not None:
        resid=v0
        info=1
    else:
        resid = np.zeros(n,typ)
        info=0

    # some sanity checks
    if k <= 0:
        raise ValueError("k must be positive, k=%d"%k)
    if k == n:
        raise ValueError("k must be less than rank(A), k=%d"%k)
    if maxiter <= 0:
        raise ValueError("maxiter must be positive, maxiter=%d"%maxiter)
    whiches=['LM','SM','LA','SA','BE']
    if which not in whiches:
        raise ValueError("which must be one of %s"%' '.join(whiches))
    if ncv > n or ncv < k:
        raise ValueError("ncv must be k<=ncv<=n, ncv=%s"%ncv)

    # assign solver and postprocessor
    ltr = _type_conv[typ]
    eigsolver = _arpack.__dict__[ltr+'saupd']
    eigextract = _arpack.__dict__[ltr+'seupd']

    # set output arrays, parameters, and workspace
    v = np.zeros((n,ncv),typ)
    workd = np.zeros(3*n,typ)
    workl = np.zeros(ncv*(ncv+8),typ)
    iparam = np.zeros(11,'int')
    ipntr = np.zeros(11,'int')
    ido = 0

    # set solver mode and parameters
    # only supported mode is 1: Ax=lx
    ishfts = 1
    mode1 = 1
    bmat='I'
    iparam[0] = ishfts
    iparam[2] = maxiter
    iparam[6] = mode1

    while True:
        ido,resid,v,iparam,ipntr,info =\
            eigsolver(ido,bmat,which,k,tol,resid,v,
                      iparam,ipntr,workd,workl,info)

        xslice = slice(ipntr[0]-1, ipntr[0]-1+n)
        yslice = slice(ipntr[1]-1, ipntr[1]-1+n)
        if ido == -1:
            # initialization
            workd[yslice]=A.matvec(workd[xslice])
        elif ido == 1:
            # compute y=Ax
            workd[yslice]=A.matvec(workd[xslice])
        else:
            break

    if info < -1 :
        raise RuntimeError("Error info=%d in arpack" % info)
        return None

    if info == 1:
        warnings.warn("Maximum number of iterations taken: %s" % iparam[2])

    if iparam[4] < k:
        warnings.warn("Only %d/%d eigenvectors converged" % (iparam[4], k))

    # now extract eigenvalues and (optionally) eigenvectors
    rvec = return_eigenvectors
    ierr = 0
    howmny = 'A' # return all eigenvectors
    sselect = np.zeros(ncv,'int') # unused
    sigma = 0.0 # no shifts, not implemented

    d,z,info =\
             eigextract(rvec,howmny,sselect,sigma,
                        bmat,which, k,tol,resid,v,iparam[0:7],ipntr,
                        workd[0:2*n],workl,ierr)

    if ierr != 0:
        raise RuntimeError("Error info=%d in arpack"%info)
        return None
    if return_eigenvectors:
        return d,z
    return d

def svd(A, k=6):
    """Compute a few singular values/vectors for a sparse matrix using ARPACK.

    Parameters
    ----------
    A: sparse matrix
        Array to compute the SVD on.
    k: int
        Number of singular values and vectors to compute.

    Note
    ----
    This is a naive implementation using the symmetric eigensolver on A.T * A
    or A * A.T, depending on which one is more efficient.

    Complex support is not implemented yet
    """
    # TODO: implement complex support once ARPACK-based eigen_hermitian is
    # available
    n, m = A.shape

    if np.iscomplexobj(A):
        raise NotImplementedError("Complex support for sparse SVD not " \
                                  "implemented yet")
        op = lambda x: x.T.conjugate()
    else:
        op = lambda x: x.T

    def _left(x):
        x = csc_matrix(x)
        m = op(x) * x

        eigvals, eigvec = eigen_symmetric(m, k)
        s = np.sqrt(eigvals)

        v = eigvec
        u = (x * v) / s
        return u, s, op(v)

    def _right(x):
        x = csr_matrix(x)
        m = x * op(x)

        eigvals, eigvec = eigen_symmetric(m, k)
        s = np.sqrt(eigvals)

        u = eigvec
        vh = (op(u) * x) / s[:, None]
        return u, s, vh

    if n > m:
        return _left(A)
    else:
        return _right(A)
