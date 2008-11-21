import numpy as np
import _arpack

__all___=['ArpackException','ARPACK_eigs', 'ARPACK_gen_eigs']

class ArpackException(RuntimeError):
    ARPACKErrors = { 0: """Normal exit.""",
                     3: """No shifts could be applied during a cycle of the
                     Implicitly restarted Arnoldi iteration. One possibility
                     is to increase the size of NCV relative to NEV.""",
                     -1: """N must be positive.""",
                     -2: """NEV must be positive.""",
                     -3: """NCV-NEV >= 2 and less than or equal to N.""",
                     -4: """The maximum number of Arnoldi update iteration
                     must be greater than zero.""",
                     -5: """WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'""",
                     -6: """BMAT must be one of 'I' or 'G'.""",
                     -7: """Length of private work array is not sufficient.""",
                     -8: """Error return from LAPACK eigenvalue calculation;""",
                     -9: """Starting vector is zero.""",
                     -10: """IPARAM(7) must be 1,2,3,4.""",
                     -11: """IPARAM(7) = 1 and BMAT = 'G' are incompatable.""",
                     -12: """IPARAM(1) must be equal to 0 or 1.""",
                     -9999: """Could not build an Arnoldi factorization.
                     IPARAM(5) returns the size of the current Arnoldi
                     factorization.""",
                     }
    def __init__(self, info):
        self.info = info
    def __str__(self):
        try: return self.ARPACKErrors[self.info]
        except KeyError: return "Unknown ARPACK error"

def check_init(n, nev, ncv):
    assert(nev <= n-4)  # ARPACK seems to cause a segfault otherwise
    if ncv is None:
        ncv = min(2*nev+1, n-1)
    maxitr = max(n, 1000)       # Maximum number of iterations
    return ncv, maxitr

def init_workspaces(n,nev,ncv):
    ipntr = np.zeros(14, np.int32) # Pointers into memory structure used by F77 calls
    d = np.zeros((ncv, 3), np.float64, order='FORTRAN') # Temp workspace
    # Temp workspace/error residuals upon iteration completion
    resid = np.zeros(n, np.float64)
    workd = np.zeros(3*n, np.float64) # workspace
    workl = np.zeros(3*ncv*ncv+6*ncv, np.float64) # more workspace
    # Storage for the Arnoldi basis vectors
    v = np.zeros((n, ncv), dtype=np.float64, order='FORTRAN')
    return (ipntr, d, resid, workd, workl, v)

def init_debug():
    # Causes various debug info to be printed by ARPACK
    _arpack.debug.ndigit = -3
    _arpack.debug.logfil = 6
    _arpack.debug.mnaitr = 0
    _arpack.debug.mnapps = 0
    _arpack.debug.mnaupd = 1
    _arpack.debug.mnaup2 = 0
    _arpack.debug.mneigh = 0
    _arpack.debug.mneupd = 1

def init_postproc_workspace(n, nev, ncv):
    # Used as workspace and to return eigenvectors if requested. Not touched if
    # eigenvectors are not requested
    workev = np.zeros(3*ncv, np.float64) # More workspace
    select = np.zeros(ncv, np.int32) # Used as internal workspace since dneupd
                                   # parameter HOWMNY == 'A'
    return (workev, select)

def postproc(n, nev, ncv, sigmar, sigmai, bmat, which,
             tol, resid, v, iparam, ipntr, workd, workl, info):
    workev, select = init_postproc_workspace(n, nev, ncv)
    ierr = 0
    # Postprocess the Arnouldi vectors to extract eigenvalues/vectors
    # If dneupd's first paramter is 'True' the eigenvectors are also calculated,
    # 'False' only the eigenvalues
    dr,di,z,info = _arpack.dneupd(
        True, 'A', select, sigmar, sigmai, workev, bmat, which, nev, tol, resid, v,
        iparam, ipntr, workd, workl, info)

    if np.abs(di[:-1]).max() == 0: dr = dr[:-1]
    else: dr =  dr[:-1] + 1j*di[:-1]
    return (dr, z[:,:-1])


def ARPACK_eigs(matvec, n, nev, which='SM', ncv=None, tol=1e-14):
    """
    Calculate eigenvalues for system with matrix-vector product matvec, dimension n

    Arguments
    =========
    matvec -- Function that provides matrix-vector product, i.e. matvec(x) -> A*x
    n -- Matrix dimension of the problem
    nev -- Number of eigenvalues to calculate
    which -- Spectrum selection. See details below. Defaults to 'SM'
    ncv -- Number of Arnoldi basisvectors to use. If None, default to 2*nev+1
    tol -- Numerical tollerance for Arnouldi iteration convergence. Defaults to 1e-14

    Spectrum Selection
    ==================
    which can take one of several values:

    'LM' -> Request eigenvalues with largest magnitude.
    'SM' -> Request eigenvalues with smallest magnitude.
    'LR' -> Request eigenvalues with largest real part.
    'SR' -> Request eigenvalues with smallest real part.
    'LI' -> Request eigenvalues with largest imaginary part.
    'SI' -> Request eigenvalues with smallest imaginary part.

    Return Values
    =============
    (eig_vals, eig_vecs) where eig_vals are the requested eigenvalues and
    eig_vecs the corresponding eigenvectors. If all the eigenvalues are real,
    eig_vals is a real array but if some eigenvalues are complex it is a
    complex array.

    """
    bmat = 'I'                          # Standard eigenproblem
    ncv, resid, iparam, ipntr, v, workd, workl, info = ARPACK_iteration(
        matvec, lambda x: x, n, bmat, which, nev, tol, ncv, mode=1)
    return postproc(n, nev, ncv, 0., 0., bmat, which, tol,
                    resid, v, iparam, ipntr, workd, workl, info)

def ARPACK_gen_eigs(matvec, sigma_solve, n, sigma, nev, which='LR', ncv=None, tol=1e-14):
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
    which -- Spectrum selection. See details below. Defaults to 'LR'
    ncv -- Number of Arnoldi basisvectors to use. If None, default to 2*nev+1
    tol -- Numerical tollerance for Arnouldi iteration convergence. Defaults to 1e-14

    Spectrum Shift
    ==============

    The spectrum of the orignal system is shifted by sigma. This transforms the
    original eigenvalues to be 1/(original_eig-sigma) in the shifted
    system. ARPACK then operates on the shifted system, transforming it back to
    the original system in a postprocessing step.

    The spectrum shift causes eigenvalues close to sigma to become very large
    in the transformed system. This allows quick convergence for these
    eigenvalues. This is particularly useful if a system has a number of
    trivial zero-eigenvalues that are to be ignored.

    Spectrum Selection
    ==================
    which can take one of several values:

    'LM' -> Request spectrum shifted eigenvalues with largest magnitude.
    'SM' -> Request spectrum shifted eigenvalues with smallest magnitude.
    'LR' -> Request spectrum shifted eigenvalues with largest real part.
    'SR' -> Request spectrum shifted eigenvalues with smallest real part.
    'LI' -> Request spectrum shifted eigenvalues with largest imaginary part.
    'SI' -> Request spectrum shifted eigenvalues with smallest imaginary part.

    The effect on the actual system is:
    'LM' -> Eigenvalues closest to sigma on the complex plane
    'LR' -> Eigenvalues with real part > sigma, provided they exist


    Return Values
    =============
    (eig_vals, eig_vecs) where eig_vals are the requested eigenvalues and
    eig_vecs the corresponding eigenvectors. If all the eigenvalues are real,
    eig_vals is a real array but if some eigenvalues are complex it is a
    complex array. The eigenvalues and vectors correspond to the original
    system, not the shifted system. The shifted system is only used interally.

    """
    bmat = 'G'                          # Generalised eigenproblem
    ncv, resid, iparam, ipntr, v, workd, workl, info = ARPACK_iteration(
        matvec, sigma_solve, n, bmat, which, nev, tol, ncv, mode=3)
    sigmar = sigma
    sigmai = 0.
    return postproc(n, nev, ncv, sigmar, sigmai, bmat, which, tol,
                    resid, v, iparam, ipntr, workd, workl, info)

def ARPACK_iteration(matvec, sigma_solve, n, bmat, which, nev, tol, ncv, mode):
    ncv, maxitr = check_init(n, nev, ncv)
    ipntr, d, resid, workd, workl, v = init_workspaces(n,nev,ncv)
    #init_debug()
    ishfts = 1         # Some random arpack parameter
    # Some random arpack parameter (I think it tells ARPACK to solve the
    # general eigenproblem using shift-invert
    iparam = np.zeros(11, np.int32) # Array with assorted extra paramters for F77 call
    iparam[[0,2,6]] = ishfts, maxitr, mode
    ido = 0                # Communication variable used by ARPACK to tell the user what to do
    info = 0               # Used for error reporting
    # Arnouldi iteration.
    while True:
        ido,resid,v,iparam,ipntr,info = _arpack.dnaupd(
            ido, bmat, which, nev, tol, resid, v, iparam, ipntr, workd, workl, info)
        if ido == -1 or ido == 1 and mode not in (3,4):
            # Perform y = inv[A - sigma*M]*M*x
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
        if info == 1:
            warn.warn("Maximum number of iterations taken: %s"%iparam[2])
        elif info != 0:
            raise ArpackException(info)

    return (ncv, resid, iparam, ipntr, v, workd, workl, info)
