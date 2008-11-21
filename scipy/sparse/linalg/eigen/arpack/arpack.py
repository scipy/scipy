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

__all___=['eigen','eigen_symmetric']

import warnings

import _arpack
import numpy as np
from scipy.sparse.linalg.interface import aslinearoperator

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}
_ndigits = {'f':5, 'd':12, 'F':5, 'D':12}


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
        raise ValueError('expected square matrix (shape=%s)' % shape)
    n = A.shape[0]

    # guess type
    typ = A.dtype.char
    if typ not in 'fdFD':
        raise ValueError("matrix type must be 'f', 'd', 'F', or 'D'")

    if M is not None:
        raise NotImplementedError("generalized eigenproblem not supported yet")
    if sigma is not None:
        raise NotImplementedError("shifted eigenproblem not supported yet")


    # some defaults
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
    whiches=['LM','SM','LR','SR','LI','SI']
    if which not in whiches:
        raise ValueError("which must be one of %s"%' '.join(whiches))
    if ncv > n or ncv < k:
        raise ValueError("ncv must be k<=ncv<=n, ncv=%s"%ncv)

    # assign solver and postprocessor
    ltr = _type_conv[typ]
    eigsolver = _arpack.__dict__[ltr+'naupd']
    eigextract = _arpack.__dict__[ltr+'neupd']

    v = np.zeros((n,ncv),typ) # holds Ritz vectors
    workd = np.zeros(3*n,typ) # workspace
    workl = np.zeros(3*ncv*ncv+6*ncv,typ) # workspace
    iparam = np.zeros(11,'int') # problem parameters
    ipntr = np.zeros(14,'int') # pointers into workspaces
    ido = 0

    if typ in 'FD':
        rwork = np.zeros(ncv,typ.lower())

    # set solver mode and parameters
    # only supported mode is 1: Ax=lx
    ishfts = 1
    mode1 = 1
    bmat = 'I'
    iparam[0] = ishfts
    iparam[2] = maxiter
    iparam[6] = mode1

    while True:
        if typ in 'fd':
            ido,resid,v,iparam,ipntr,info =\
                eigsolver(ido,bmat,which,k,tol,resid,v,iparam,ipntr,
                          workd,workl,info)
        else:
            ido,resid,v,iparam,ipntr,info =\
                eigsolver(ido,bmat,which,k,tol,resid,v,iparam,ipntr,
                          workd,workl,rwork,info)

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

    if  info < -1 :
        raise RuntimeError("Error info=%d in arpack"%info)
        return None
    if info == -1:
        warnings.warn("Maximum number of iterations taken: %s"%iparam[2])
#    if iparam[3] != k:
#        warnings.warn("Only %s eigenvalues converged"%iparam[3])


    # now extract eigenvalues and (optionally) eigenvectors
    rvec = return_eigenvectors
    ierr = 0
    howmny = 'A' # return all eigenvectors
    sselect = np.zeros(ncv,'int') # unused
    sigmai = 0.0 # no shifts, not implemented
    sigmar = 0.0 # no shifts, not implemented
    workev = np.zeros(3*ncv,typ)

    if typ in 'fd':
        dr=np.zeros(k+1,typ)
        di=np.zeros(k+1,typ)
        zr=np.zeros((n,k+1),typ)
        dr,di,zr,info=\
            eigextract(rvec,howmny,sselect,sigmar,sigmai,workev,
                   bmat,which,k,tol,resid,v,iparam,ipntr,
                   workd,workl,info)

        # The ARPACK nonsymmetric real and double interface (s,d)naupd return
        # eigenvalues and eigenvectors in real (float,double) arrays.

        # Build complex eigenvalues from real and imaginary parts
        d=dr+1.0j*di

        # Arrange the eigenvectors: complex eigenvectors are stored as
        # real,imaginary in consecutive columns
        z=zr.astype(typ.upper())
        eps=np.finfo(typ).eps
        i=0
        while i<=k:
            # check if complex
            if abs(d[i].imag)>eps:
                # assume this is a complex conjugate pair with eigenvalues
                # in consecutive columns
                z[:,i]=zr[:,i]+1.0j*zr[:,i+1]
                z[:,i+1]=z[:,i].conjugate()
                i+=1
            i+=1

        # Now we have k+1 possible eigenvalues and eigenvectors
        # Return the ones specified by the keyword "which"
        nreturned=iparam[4] # number of good eigenvalues returned
        if nreturned==k:    # we got exactly how many eigenvalues we wanted
            d=d[:k]
            z=z[:,:k]
        else:   # we got one extra eigenvalue (likely a cc pair, but which?)
            # cut at approx precision for sorting
            rd=np.round(d,decimals=_ndigits[typ])
            if which in ['LR','SR']:
                ind=np.argsort(rd.real)
            elif which in ['LI','SI']:
                # for LI,SI ARPACK returns largest,smallest abs(imaginary) why?
                ind=np.argsort(abs(rd.imag))
            else:
                ind=np.argsort(abs(rd))
            if which in ['LR','LM','LI']:
                d=d[ind[-k:]]
                z=z[:,ind[-k:]]
            if which in ['SR','SM','SI']:
                d=d[ind[:k]]
                z=z[:,ind[:k]]


    else:
        # complex is so much simpler...
        d,z,info =\
              eigextract(rvec,howmny,sselect,sigmar,workev,
                         bmat,which,k,tol,resid,v,iparam,ipntr,
                         workd,workl,rwork,ierr)



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
