import numpy as N
import _arpack 

def eigvals(matvec, n, nev, ncv=None):
    """
    Calculate eigenvalues for system with matrix-vector product matvec, dimension n

    Arguments
    =========
    matvec -- Function that provides matrix-vector product, i.e. matvec(x) -> A*x
    n -- Matrix dimension of the problem
    nev -- Number of eigenvalues to calculate
    ncv -- Number of Arnoldi basisvectors to use. If None, default to 2*nev 

    Return Values
    =============
    Real array of nev eigenvalues, or complex array if values are complex
    """

    assert(nev <= n-4)  # ARPACK seems to cause a segfault otherwise
    if ncv is None:
        ncv = min(2*nev, n-1)
    iparam = N.zeros(11, N.int32) # Array with assorted extra paramters for F77 call
    ishfts = 1         # Some random arpack parameter
    maxitr = n*3       # Maximum number of iterations
    # Some random arpack parameter (I think it tells ARPACK to solve the
    # regular eigenproblem the "standard" way
    mode1 = 1      
    iparam[0] = ishfts
    iparam[2] = maxitr
    iparam[6] = mode1
    ipntr = N.zeros(14, N.int32) # Pointers into memory structure used by F77 calls
    d = N.zeros((ncv, 3), N.float64, order='FORTRAN') # Temp workspace
    # Temp workspace/error residuals upon iteration completion
    resid = N.zeros(n, N.float64) 
    
    workd = N.zeros(3*n, N.float64) # workspace
    workl = N.zeros(3*ncv*ncv+6*ncv, N.float64) # more workspace
    tol = 1e-16            # Desired convergence tollerance
    ido = 0                # Communication variable used by ARPACK to tell the user what to do
    info = 0               # Used for error reporting
    which = 'SM'           # Request smallest magnitude eigenvalues, see dnaupd.f for more options
    bmat  = 'I'            # Standard (not generalised) eigenproblem
    # Storage for the Arnoldi basis vectors
    v = N.zeros((n, ncv), dtype=float, order='FORTRAN') 
    
    # Causes various debug info to be printed by ARPACK
    _arpack.debug.ndigit = -3
    _arpack.debug.logfil = 6
    _arpack.debug.mnaitr = 0
    _arpack.debug.mnapps = 0
    _arpack.debug.mnaupd = 1
    _arpack.debug.mnaup2 = 0
    _arpack.debug.mneigh = 0
    _arpack.debug.mneupd = 1


    cnt = 0
    # Arnouldi iteration.
    while True:
        ido,resid,v,iparam,ipntr,info = _arpack.dnaupd(
            ido, bmat, which, nev, tol, resid, v, iparam, ipntr, workd, workl, info)
        # Exit if reverse communication flag does not request a matrix-vector
        # product
        if ido not in  (-1, 1): break 
        cnt += 1
        x = workd[ipntr[0]-1:ipntr[0]+n-1]
        workd[ipntr[1]-1:ipntr[1]+n-1] = matvec(x)    # y = A*x


    if info != 0: raise "Hell" # Indicates some error during the Arnouldi iterations

    dr = N.zeros(nev+1, N.float64) # Real part of eigenvalues
    di = N.zeros(nev+1, N.float64) # Imaginary part of eigenvalues
    # Used as workspace and to return eigenvectors if requested. Not touched if
    # eigenvectors are not requested
    z = N.zeros((n, nev+1), N.float64, order='FORTRAN')
    workev = N.zeros(3*ncv, N.float64) # More workspace
    select = N.zeros(ncv, N.int32) # Used as internal workspace since dneupd
                                   # parameter HOWMNY == 'A'

    # Postprocess the Arnouldi vectors to extract eigenvalues/vectors
    # If dneupd's first paramter is 'True' the eigenvectors are also calculated,
    # 'False' only the eigenvalues
    dr,di,z,info = _arpack.dneupd(
        True, 'A', select, 0., 0., workev, bmat, which, nev, tol,
        resid, v, iparam, ipntr, workd, workl, info)
    if N.abs(di[:-1]).max() == 0: dr = dr[:-1]
    else: dr =  dr[:-1] + 1j*di[:-1]
    return (dr, z[:,:-1])

def geneigvals(matvec, sigma_solve, n, sigma, nev, ncv=None):
    """
    Calculate eigenvalues close to sigma for generalised eigen system

    Given a system [A]x = k_i*[M]x where [A] and [M] are matrices and k_i are
    eigenvalues, nev eigenvalues close to sigma are calculated. The user needs
    to provide routines that calculate [M]*x and solve [A]-sigma*[M]*x = b for x.

    Arguments
    =========
    matvec -- Function that provides matrix-vector product, i.e. matvec(x) -> [M]*x
    sigma_solve -- sigma_solve(b) -> x, where [A]-sigma*[M]*x = b
    n -- Matrix dimension of the problem
    sigma -- Eigenvalue spectral shift real value
    nev -- Number of eigenvalues to calculate
    ncv -- Number of Arnoldi basisvectors to use. If None, default to 2*nev 

    Return Values
    =============
    Real array of nev eigenvalues, or complex array if values are complex
    """
    if ncv is None:
        ncv = min(2*nev, n-1)

    iparam = N.zeros(11, N.int32) # Array with assorted extra paramters for F77 call
    ishfts = 1         # Some random arpack parameter
    maxitr = n*3       # Maximum number of iterations
    # Some random arpack parameter (I think it tells ARPACK to solve the
    # regular eigenproblem the "standard" way
    mode = 3
    iparam[0] = ishfts
    iparam[2] = maxitr
    iparam[6] = mode
    ipntr = N.zeros(14, N.int32) # Pointers into memory structure used by F77 calls
    d = N.zeros((ncv, 3), N.float64, order='FORTRAN') # Temp workspace
    # Temp workspace/error residuals upon iteration completion
    resid = N.zeros(n, N.float64) 
    
    workd = N.zeros(3*n, N.float64) # workspace
    workl = N.zeros(3*ncv*ncv+6*ncv, N.float64) # more workspace
    tol = 1e-7             # Desired convergence tollerance
    ido = 0                # Communication variable used by ARPACK to tell the user what to do
    info = 0               # Used for error reporting
    which = 'LR'           # Request largest magnitude eigenvalues, see dnaupd.f for more options
    bmat  = 'G'            # Generalised eigenproblem
    # Storage for the Arnoldi basis vectors
    v = N.zeros((n, ncv), dtype=float, order='FORTRAN') 
    # Causes various debug info to be printed by ARPACK
    _arpack.debug.ndigit = -3
    _arpack.debug.logfil = 6
    _arpack.debug.mnaitr = 0
    _arpack.debug.mnapps = 0
    _arpack.debug.mnaupd = 1
    _arpack.debug.mnaup2 = 0
    _arpack.debug.mneigh = 0
    _arpack.debug.mneupd = 1
    
    # Arnouldi iteration.
    while True:
        ido,resid,v,iparam,ipntr,info = _arpack.dnaupd(ido, bmat, which, nev, tol, resid, v,
                            iparam, ipntr, workd, workl, info)
        if ido == -1:  # Perform y = inv[A - sigma*M]*M*x
            x = workd[ipntr[0]-1:ipntr[0]+n-1]
            Mx = matvec(x)    # Mx = [M]*x
            workd[ipntr[1]-1:ipntr[1]+n-1] = sigma_solve(Mx)
        elif ido == 1: # Perform y = inv[A - sigma*M]*M*x using saved M*x
            # Mx = [M]*x where it was saved by ARPACK
            Mx = workd[ipntr[2]-1:ipntr[2]+n-1]
            workd[ipntr[1]-1:ipntr[1]+n-1] = sigma_solve(Mx)
        elif ido == 2: # Perform y = M*x
            x = workd[ipntr[0]-1:ipntr[0]+n-1]
            workd[ipntr[1]-1:ipntr[1]+n-1] = matvec(x)
        else:          # Finished, or error
            break
        if info != 0: raise "Hell" # Indicates some error during the Arnouldi iterations
    
    # Used as workspace and to return eigenvectors if requested. Not touched if
    # eigenvectors are not requested
    z = N.zeros((n, nev+1), N.float64, order='FORTRAN')
    workev = N.zeros(3*ncv, N.float64) # More workspace
    ierr = 0
    select = N.zeros(ncv, N.int32) # Used as internal workspace since dneupd
                                   # parameter HOWMNY == 'A'
    sigmar = sigma
    sigmai = 0.
    # Postprocess the Arnouldi vectors to extract eigenvalues/vectors
    # If dneupd's first paramter is 'True' the eigenvectors are also calculated,
    # 'False' only the eigenvalues
    dr,di,z,info = _arpack.dneupd(
        True, 'A', select, sigmar, sigmai, workev, bmat, which, nev, tol, resid, v,
        iparam, ipntr, workd, workl, info)

    if N.abs(di[:-1]).max() == 0:
        return dr[:-1]
    else:
        return dr[:-1] + 1j*di[:-1]

