"""
arpack - Scipy module to find a few eigenvectors and eigenvalues of a matrix

Uses ARPACK: http://www.caam.rice.edu/software/ARPACK/

"""
__all___=['eigen','eigen_symmetric']

import _arpack 
import numpy as sb
import warnings

# inspired by iterative.py
# so inspired, in fact, that some of it was copied directly
try:
    False, True
except NameError:
    False, True = 0, 1

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}

class get_matvec:
    methname = 'matvec'
    def __init__(self, obj, *args):
        self.obj = obj
        self.args = args
        if isinstance(obj, sb.matrix):
            self.callfunc = self.type1m
            return
        if isinstance(obj, sb.ndarray):
            self.callfunc = self.type1
            return
        meth = getattr(obj,self.methname,None)
        if not callable(meth):
            raise ValueError, "Object must be an array "\
                  "or have a callable %s attribute." % (self.methname,)

        self.obj = meth
        self.callfunc = self.type2

    def __call__(self, x):
        return self.callfunc(x)

    def type1(self, x):
        return sb.dot(self.obj, x)

    def type1m(self, x):
        return sb.dot(self.obj.A, x)

    def type2(self, x):
        return self.obj(x,*self.args)


def eigen(A,k=6,M=None,ncv=None,which='LM',
          maxiter=None,tol=0, return_eigenvectors=True):
    """ Return k eigenvalues and eigenvectors of the matrix A.

    Solves A * x[i] = w[i] * x[i], the standard eigenvalue problem for
    w[i] eigenvalues with corresponding eigenvectors x[i].

    Inputs:

    A --  A matrix, array or an object with matvec(x) method to perform
          the matrix vector product A * x.  The sparse matrix formats
          in scipy.sparse are appropriate for A.

    k -- The number of eigenvalue/eigenvectors desired

    M -- (Not implemented)
          A symmetric positive-definite matrix for the generalized
          eigenvalue problem A * x = w * M * x

    Outputs:

    w -- An array of k eigenvalues

    v -- An array of k eigenvectors, k[i] is the eigenvector corresponding
         to the eigenvector w[i]

    Optional Inputs:

    ncv -- Number of Lanczos vectors generated, ncv must be greater than k
           and is recommended to be ncv > 2*k

    which -- String specifying which eigenvectors to compute.
             Compute the k eigenvalues of:
             'LM' - largest magnitude.
             'SM' - smallest magnitude.
             'LR' - largest real part.
             'SR' - smallest real part.
             'LI' - largest imaginary part.
             'SI' - smallest imaginary part.

    maxiter --  Maximum number of Arnoldi update iterations allowed

    tol -- Relative accuracy for eigenvalues (stopping criterion)

    return_eigenvectors -- True|False, return eigenvectors 

    """
    try:
        n,ny=A.shape
        n==ny
    except:
        raise AttributeError("matrix is not square")
    if M is not None:
        raise NotImplementedError("generalized eigenproblem not supported yet")
        
    # some defaults
    if ncv is None:
        ncv=2*k+1
    ncv=min(ncv,n)
    if maxiter==None:
        maxiter=n*10

    # guess type        
    resid = sb.zeros(n,'f')
    try:
        typ = A.dtype.char
    except AttributeError:
        typ = A.matvec(resid).dtype.char
    if typ not in 'fdFD':
        raise ValueError("matrix type must be 'f', 'd', 'F', or 'D'")

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
    matvec = get_matvec(A)

    v = sb.zeros((n,ncv),typ) # holds Ritz vectors
    resid = sb.zeros(n,typ) # residual
    workd = sb.zeros(3*n,typ) # workspace
    workl = sb.zeros(3*ncv*ncv+6*ncv,typ) # workspace
    iparam = sb.zeros(11,'int') # problem parameters
    ipntr = sb.zeros(14,'int') # pointers into workspaces
    info = 0
    ido = 0

    if typ in 'FD':
        rwork = sb.zeros(ncv,typ.lower())

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

        if (ido == -1 or ido == 1):
            # compute y = A * x
            xslice = slice(ipntr[0]-1, ipntr[0]-1+n)
            yslice = slice(ipntr[1]-1, ipntr[1]-1+n)
            workd[yslice]=matvec(workd[xslice])
        else: # done
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
    sselect = sb.zeros(ncv,'int') # unused 
    sigmai = 0.0 # no shifts, not implemented
    sigmar = 0.0 # no shifts, not implemented
    workev = sb.zeros(3*ncv,typ) 

    if typ in 'fd':
        dr=sb.zeros(k+1,typ)
        di=sb.zeros(k+1,typ)
        zr=sb.zeros((n,k+1),typ)
        dr,di,z,info=\
            eigextract(rvec,howmny,sselect,sigmar,sigmai,workev,
                   bmat,which,k,tol,resid,v,iparam,ipntr,
                   workd,workl,info)
        
        # make eigenvalues complex
        d=dr+1.0j*di
        # futz with the eigenvectors:
        # complex are stored as real,imaginary in consecutive columns
        z=zr.astype(typ.upper())
        for i in range(k): # fix c.c. pairs
            if di[i] > 0 :
                z[:,i]=zr[:,i]+1.0j*zr[:,i+1]
                z[:,i+1]=z[:,i].conjugate()
                     
    else:
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


def eigen_symmetric(A,k=6,M=None,ncv=None,which='LM',
                    maxiter=None,tol=0, return_eigenvectors=True):
    """ Return k eigenvalues and eigenvectors of the real symmetric matrix A.

    Solves A * x[i] = w[i] * x[i], the standard eigenvalue problem for
    w[i] eigenvalues with corresponding eigenvectors x[i].
    A must be real and symmetric.
    See eigen() for nonsymmetric or complex symmetric (Hermetian) matrices.

    Inputs:

    A --  A symmetric matrix, array or an object with matvec(x) method
          to perform the matrix vector product A * x.
          The sparse matrix formats in scipy.sparse are appropriate for A.

    k -- The number of eigenvalue/eigenvectors desired

    M -- (Not implemented)
          A symmetric positive-definite matrix for the generalized
          eigenvalue problem A * x = w * M * x

    Outputs:

    w -- An real array of k eigenvalues 

    v -- An array of k real eigenvectors, k[i] is the eigenvector corresponding
         to the eigenvector w[i]

    Optional Inputs:

    ncv -- Number of Lanczos vectors generated, ncv must be greater than k
           and is recommended to be ncv > 2*k

    which -- String specifying which eigenvectors to compute.
             Compute the k 
             'LA' - largest (algebraic) eigenvalues.
             'SA' - smallest (algebraic) eigenvalues.
             'LM' - largest (in magnitude) eigenvalues.
             'SM' - smallest (in magnitude) eigenvalues. 
             'BE' - eigenvalues, half from each end of the
                    spectrum.  When NEV is odd, compute one more from the
                    high end than from the low end.

    maxiter --  Maximum number of Arnoldi update iterations allowed

    tol -- Relative accuracy for eigenvalues (stopping criterion)

    return_eigenvectors -- True|False, return eigenvectors 

    """
    try:
        n,ny=A.shape
        n==ny
    except:
        raise AttributeError("matrix is not square")
    if M is not None:
        raise NotImplementedError("generalized eigenproblem not supported yet")
    if ncv is None:
        ncv=2*k+1
    ncv=min(ncv,n)
    if maxiter==None:
        maxiter=n*10


    # guess type        
    resid = sb.zeros(n,'f')
    try:
        typ = A.dtype.char
    except AttributeError:
        typ = A.matvec(resid).dtype.char
    if typ not in 'fd':
        raise ValueError("matrix type must be 'f' or 'd'")

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
    matvec = get_matvec(A)

    v = sb.zeros((n,ncv),typ)
    resid = sb.zeros(n,typ)
    workd = sb.zeros(3*n,typ)
    workl = sb.zeros(ncv*(ncv+8),typ)
    iparam = sb.zeros(11,'int')
    ipntr = sb.zeros(11,'int')
    info = 0
    ido = 0

    # only supported mode is 1: Ax=lx
    ishfts = 1
    mode1 = 1
    bmat='I'
    iparam[0] = ishfts
    iparam[2] = maxiter
    iparam[6] = mode1


    while True:
        ido,resid,v,iparam,ipntr,info =\
            eigsolver(ido,bmat,which,k,tol,resid,v,iparam,ipntr,
               workd,workl,info)
        if (ido == -1 or ido == 1):
            xslice = slice(ipntr[0]-1, ipntr[0]-1+n)
            yslice = slice(ipntr[1]-1, ipntr[1]-1+n)
            workd[yslice]=matvec(workd[xslice])
        else:
            break

    if  info < -1 :
        raise RuntimeError("Error info=%d in arpack"%info)
        return None
    if info == -1:
        warnings.warn("Maximum number of iterations taken: %s"%iparam[2])

    # now extract eigenvalues and (optionally) eigenvectors        
    rvec = return_eigenvectors
    ierr = 0
    howmny = 'A' # return all eigenvectors
    sselect = sb.zeros(ncv,'int') # unused 
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


