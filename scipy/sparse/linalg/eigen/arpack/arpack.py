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
# This wrapper puts the *neupd (general matrix) interfaces in eigs()
# and the *seupd (symmetric matrix) in eigsh().
# There is no Hermetian complex/double complex interface.
# To find eigenvalues of a Hermetian matrix you
# must use eigs() and not eigsh()
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

__all___=['eigs', 'eigsh', 'svds', 'ArpackNoConvergence']

import _arpack
import numpy as np
from scipy.sparse.linalg.interface import aslinearoperator, LinearOperator
from scipy.sparse import csc_matrix, csr_matrix, isspmatrix

_type_conv = {'f':'s', 'd':'d', 'F':'c', 'D':'z'}
_ndigits = {'f':5, 'd':12, 'F':5, 'D':12}

_NAUPD_ERRORS = {
    0: "Normal exit.",
    1: "Maximum number of iterations taken. "
       "All possible eigenvalues of OP has been found.",
    2: "No longer an informational error. Deprecated starting with "
       "release 2 of ARPACK.",
    3: "No shifts could be applied during a cycle of the Implicitly "
       "restarted Arnoldi iteration. One possibility is to increase "
       "the size of NCV relative to NEV. ",
    -1: "N must be positive.",
    -2: "NEV must be positive.",
    -3: "NCV must be greater than NEV and less than or equal to N.",
    -4: "The maximum number of Arnoldi update iterations allowed "
        "must be greater than zero.",
    -5: "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.",
    -6: "BMAT must be one of 'I' or 'G'.",
    -7: "Length of private work array WORKL is not sufficient.",
    -8: "Error return from trid. eigenvalue calculation; "
        "Informational error from LAPACK routine dsteqr .",
    -9: "Starting vector is zero.",
    -10: "IPARAM(7) must be 1,2,3,4,5.",
    -11: "IPARAM(7) = 1 and BMAT = 'G' are incompatable.",
    -12: "IPARAM(1) must be equal to 0 or 1.",
    -13: "NEV and WHICH = 'BE' are incompatable. ",
    -9999: "Could not build an Arnoldi factorization. "
           "IPARAM(5) returns the size of the current Arnoldi "
           "factorization. The user is advised to check that "
           "enough workspace and array storage has been allocated.",
}

_NEUPD_ERRORS = {
    0: "Normal exit.",
    1: "The Schur form computed by LAPACK routine dlahqr "
       "could not be reordered by LAPACK routine dtrsen. "
       "Re-enter subroutine dneupd  with IPARAM(5)NCV and "
       "increase the size of the arrays DR and DI to have "
       "dimension at least dimension NCV and allocate at least NCV "
       "columns for Z. NOTE: Not necessary if Z and V share "
       "the same space. Please notify the authors if this error"
       "occurs.",
    -1: "N must be positive.",
    -2: "NEV must be positive.",
    -3: "NCV-NEV >= 2 and less than or equal to N.",
    -5: "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'",
    -6: "BMAT must be one of 'I' or 'G'.",
    -7: "Length of private work WORKL array is not sufficient.",
    -8: "Error return from calculation of a real Schur form. "
        "Informational error from LAPACK routine dlahqr .",
    -9: "Error return from calculation of eigenvectors. "
        "Informational error from LAPACK routine dtrevc.",
    -10: "IPARAM(7) must be 1,2,3,4.",
    -11: "IPARAM(7) = 1 and BMAT = 'G' are incompatible.",
    -12: "HOWMNY = 'S' not yet implemented",
    -13: "HOWMNY must be one of 'A' or 'P' if RVEC = .true.",
    -14: "DNAUPD  did not find any eigenvalues to sufficient "
         "accuracy.",
    -15: "DNEUPD got a different count of the number of converged "
         "Ritz values than DNAUPD got.  This indicates the user "
         "probably made an error in passing data from DNAUPD to "
         "DNEUPD or that the data was modified before entering "
         "DNEUPD",
}

_SEUPD_ERRORS = {
    0: "Normal exit.",
    -1: "N must be positive.",
    -2: "NEV must be positive.",
    -3: "NCV must be greater than NEV and less than or equal to N.",
    -5: "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.",
    -6: "BMAT must be one of 'I' or 'G'.",
    -7: "Length of private work WORKL array is not sufficient.",
    -8: ("Error return from trid. eigenvalue calculation; "
         "Information error from LAPACK routine dsteqr."),
    -9: "Starting vector is zero.",
    -10: "IPARAM(7) must be 1,2,3,4,5.",
    -11: "IPARAM(7) = 1 and BMAT = 'G' are incompatible.",
    -12: "NEV and WHICH = 'BE' are incompatible.",
    -14: "DSAUPD  did not find any eigenvalues to sufficient accuracy.",
    -15: "HOWMNY must be one of 'A' or 'S' if RVEC = .true.",
    -16: "HOWMNY = 'S' not yet implemented",
    -17: ("DSEUPD  got a different count of the number of converged "
          "Ritz values than DSAUPD  got.  This indicates the user "
          "probably made an error in passing data from DSAUPD  to "
          "DSEUPD  or that the data was modified before entering  "
          "DSEUPD.")
}

class ArpackError(RuntimeError):
    """
    ARPACK error
    """
    def __init__(self, info, infodict=_NAUPD_ERRORS):
        msg = infodict.get(info, "Unknown error")
        RuntimeError.__init__(self, "ARPACK error %d: %s" % (info, msg))

class ArpackNoConvergence(ArpackError):
    """
    ARPACK iteration did not converge

    Attributes
    ----------
    eigenvalues : ndarray
        Partial result. Converged eigenvalues.
    eigenvectors : ndarray
        Partial result. Converged eigenvectors.

    """
    def __init__(self, msg, eigenvalues, eigenvectors):
        ArpackError.__init__(self, -1, {-1: msg})
        self.eigenvalues = eigenvalues
        self.eigenvectors = eigenvectors

class _ArpackParams(object):
    def __init__(self, n, k, tp, matvec, sigma=None,
                 ncv=None, v0=None, maxiter=None, which="LM", tol=0):
        if k <= 0:
            raise ValueError("k must be positive, k=%d" % k)

        if maxiter is None:
            maxiter = n * 10
        if maxiter <= 0:
            raise ValueError("maxiter must be positive, maxiter=%d" % maxiter)

        if tp not in 'fdFD':
            raise ValueError("matrix type must be 'f', 'd', 'F', or 'D'")

        if v0 is not None:
            # ARPACK overwrites its initial resid,  make a copy
            self.resid = np.array(v0, copy=True)
            info = 1
        else:
            self.resid = np.zeros(n, tp)
            info = 0

        if sigma is not None:
            raise NotImplementedError("shifted eigenproblem not supported yet")

        if ncv is None:
            ncv = 2 * k + 1
        ncv = min(ncv, n)

        self.v = np.zeros((n, ncv), tp) # holds Ritz vectors
        self.iparam = np.zeros(11, "int")

        # set solver mode and parameters
        # only supported mode is 1: Ax=lx
        ishfts = 1
        mode1 = 1
        self.iparam[0] = ishfts
        self.iparam[2] = maxiter
        self.iparam[6] = mode1

        self.n = n
        self.matvec = matvec
        self.tol = tol
        self.k = k
        self.maxiter = maxiter
        self.ncv = ncv
        self.which = which
        self.tp = tp
        self.info = info
        self.bmat = 'I'

        self.converged = False
        self.ido = 0

    def _raise_no_convergence(self):
        msg = "No convergence (%d iterations, %d/%d eigenvectors converged)"
        k_ok = self.iparam[4]
        num_iter = self.iparam[2]
        try:
            ev, vec = self.extract(True)
        except ArpackError, err:
            msg = "%s [%s]" % (msg, err)
            ev = np.zeros((0,))
            vec = np.zeros((self.n, 0))
            k_ok = 0
        raise ArpackNoConvergence(msg % (num_iter, k_ok, self.k), ev, vec)

class _SymmetricArpackParams(_ArpackParams):
    def __init__(self, n, k, tp, matvec, sigma=None,
                 ncv=None, v0=None, maxiter=None, which="LM", tol=0):
        if not which in ['LM', 'SM', 'LA', 'SA', 'BE']:
            raise ValueError("which must be one of %s" % ' '.join(whiches))
        if k >= n:
            raise ValueError("k must be less than rank(A), k=%d" % k)

        _ArpackParams.__init__(self, n, k, tp, matvec, sigma,
                 ncv, v0, maxiter, which, tol)

        if self.ncv > n or self.ncv <= k:
            raise ValueError("ncv must be k<ncv<=n, ncv=%s" % self.ncv)

        self.workd = np.zeros(3 * n, self.tp)
        self.workl = np.zeros(self.ncv * (self.ncv + 8), self.tp)

        ltr = _type_conv[self.tp]
        if ltr not in ["s", "d"]:
            raise ValueError("Input matrix is not real-valued.")
        self._arpack_solver = _arpack.__dict__[ltr + 'saupd']
        self._arpack_extract = _arpack.__dict__[ltr + 'seupd']

        self.ipntr = np.zeros(11, "int")

    def iterate(self):
        self.ido, self.resid, self.v, self.iparam, self.ipntr, self.info = \
            self._arpack_solver(self.ido, self.bmat, self.which, self.k, self.tol,
                    self.resid, self.v, self.iparam, self.ipntr,
                    self.workd, self.workl, self.info)

        xslice = slice(self.ipntr[0]-1, self.ipntr[0]-1+self.n)
        yslice = slice(self.ipntr[1]-1, self.ipntr[1]-1+self.n)
        if self.ido == -1:
            # initialization
            self.workd[yslice] = self.matvec(self.workd[xslice])
        elif self.ido == 1:
            # compute y=Ax
            self.workd[yslice] = self.matvec(self.workd[xslice])
        else:
            self.converged = True

            if self.info == 0:
                pass
            elif self.info == 1:
                self._raise_no_convergence()
            else:
                raise ArpackError(self.info)

    def extract(self, return_eigenvectors):
        rvec = return_eigenvectors
        ierr = 0
        howmny = 'A' # return all eigenvectors
        sselect = np.zeros(self.ncv, 'int') # unused
        sigma = 0.0 # no shifts, not implemented

        d, z, ierr = self._arpack_extract(rvec, howmny, sselect, sigma, self.bmat,
                self.which, self.k, self.tol, self.resid, self.v,
                self.iparam[0:7], self.ipntr, self.workd[0:2*self.n],
                self.workl,ierr)

        if ierr != 0:
            raise ArpackError(ierr, infodict=_SEUPD_ERRORS)

        k_ok = self.iparam[4]
        d = d[:k_ok]
        z = z[:,:k_ok]

        if return_eigenvectors:
            return d, z
        else:
            return d

class _UnsymmetricArpackParams(_ArpackParams):
    def __init__(self, n, k, tp, matvec, sigma=None,
                 ncv=None, v0=None, maxiter=None, which="LM", tol=0):
        if not which in ["LM", "SM", "LR", "SR", "LI", "SI"]:
            raise ValueError("Parameter which must be one of %s" % ' '.join(whiches))
        if k >= n-1:
            raise ValueError("k must be less than rank(A)-1, k=%d" % k)

        _ArpackParams.__init__(self, n, k, tp, matvec, sigma,
                 ncv, v0, maxiter, which, tol)

        if self.ncv > n or self.ncv <= k+1:
            raise ValueError("ncv must be k+1<ncv<=n, ncv=%s" % self.ncv)

        self.workd = np.zeros(3 * n, self.tp)
        self.workl = np.zeros(3 * self.ncv * self.ncv + 6 * self.ncv, self.tp)

        ltr = _type_conv[self.tp]
        self._arpack_solver = _arpack.__dict__[ltr + 'naupd']
        self._arpack_extract = _arpack.__dict__[ltr + 'neupd']

        self.ipntr = np.zeros(14, "int")

        if self.tp in 'FD':
            self.rwork = np.zeros(self.ncv, self.tp.lower())
        else:
            self.rwork = None

    def iterate(self):
        if self.tp in 'fd':
            self.ido, self.resid, self.v, self.iparam, self.ipntr, self.info = \
                self._arpack_solver(self.ido, self.bmat, self.which, self.k, self.tol,
                        self.resid, self.v, self.iparam, self.ipntr,
                        self.workd, self.workl, self.info)
        else:
            self.ido, self.resid, self.v, self.iparam, self.ipntr, self.info =\
                self._arpack_solver(self.ido, self.bmat, self.which, self.k, self.tol,
                        self.resid, self.v, self.iparam, self.ipntr,
                        self.workd, self.workl, self.rwork, self.info)

        xslice = slice(self.ipntr[0]-1, self.ipntr[0]-1+self.n)
        yslice = slice(self.ipntr[1]-1, self.ipntr[1]-1+self.n)
        if self.ido == -1:
            # initialization
            self.workd[yslice] = self.matvec(self.workd[xslice])
        elif self.ido == 1:
            # compute y=Ax
            self.workd[yslice] = self.matvec(self.workd[xslice])
        else:
            self.converged = True

            if self.info == 0:
                pass
            elif self.info == 1:
                self._raise_no_convergence()
            else:
                raise ArpackError(self.info)

    def extract(self, return_eigenvectors):
        k, n = self.k, self.n

        ierr = 0
        howmny = 'A' # return all eigenvectors
        sselect = np.zeros(self.ncv, 'int') # unused
        sigmai = 0.0 # no shifts, not implemented
        sigmar = 0.0 # no shifts, not implemented
        workev = np.zeros(3 * self.ncv, self.tp)

        if self.tp in 'fd':
            dr = np.zeros(k+1, self.tp)
            di = np.zeros(k+1, self.tp)
            zr = np.zeros((n, k+1), self.tp)
            dr, di, zr, ierr=\
                self._arpack_extract(return_eigenvectors,
                       howmny, sselect, sigmar, sigmai, workev,
                       self.bmat, self.which, k, self.tol, self.resid,
                       self.v, self.iparam, self.ipntr,
                       self.workd, self.workl, self.info)

            if ierr != 0:
                raise ArpackError(ierr, infodict=_NEUPD_ERRORS)

            # The ARPACK nonsymmetric real and double interface (s,d)naupd
            # return eigenvalues and eigenvectors in real (float,double) arrays.

            # Build complex eigenvalues from real and imaginary parts
            d = dr + 1.0j * di

            # Arrange the eigenvectors: complex eigenvectors are stored as
            # real,imaginary in consecutive columns
            z = zr.astype(self.tp.upper())
            i = 0
            while i<=k:
                # check if complex
                if abs(d[i].imag) != 0:
                    # this is a complex conjugate pair with eigenvalues
                    # in consecutive columns
                    z[:,i] = zr[:,i] + 1.0j * zr[:,i+1]
                    z[:,i+1] = z[:,i].conjugate()
                    i +=1
                i += 1

            # Now we have k+1 possible eigenvalues and eigenvectors
            # Return the ones specified by the keyword "which"
            nreturned = self.iparam[4] # number of good eigenvalues returned
            if nreturned <= k:
                # we got less or equal as many eigenvalues we wanted
                d = d[:nreturned]
                z = z[:,:nreturned]
            else:
                # we got one extra eigenvalue (likely a cc pair, but which?)
                # cut at approx precision for sorting
                rd = np.round(d, decimals = _ndigits[self.tp])
                if self.which in ['LR','SR']:
                    ind = np.argsort(rd.real)
                elif self.which in ['LI','SI']:
                    # for LI,SI ARPACK returns largest,smallest
                    # abs(imaginary) why?
                    ind = np.argsort(abs(rd.imag))
                else:
                    ind = np.argsort(abs(rd))
                if self.which in ['LR','LM','LI']:
                    d = d[ind[-k:]]
                    z = z[:,ind[-k:]]
                if self.which in ['SR','SM','SI']:
                    d = d[ind[:k]]
                    z = z[:,ind[:k]]

        else:
            # complex is so much simpler...
            d, z, ierr =\
                    self._arpack_extract(return_eigenvectors,
                           howmny, sselect, sigmar, workev,
                           self.bmat, self.which, k, self.tol, self.resid,
                           self.v, self.iparam, self.ipntr,
                           self.workd, self.workl, self.rwork, ierr)

            if ierr != 0:
                raise ArpackError(ierr, infodict=_NEUPD_ERRORS)

            k_ok = self.iparam[4]
            d = d[:k_ok]
            z = z[:,:k_ok]


        if return_eigenvectors:
            return d, z
        else:
            return d

def _aslinearoperator_with_dtype(m):
    m = aslinearoperator(m)
    if not hasattr(m, 'dtype'):
        x = np.zeros(m.shape[1])
        m.dtype = (m*x).dtype
    return m

def eigs(A, k=6, M=None, sigma=None, which='LM', v0=None,
         ncv=None, maxiter=None, tol=0,
         return_eigenvectors=True):
    """
    Find k eigenvalues and eigenvectors of the square matrix A.

    Solves ``A * x[i] = w[i] * x[i]``, the standard eigenvalue problem
    for w[i] eigenvalues with corresponding eigenvectors x[i].

    Parameters
    ----------
    A : matrix, array, or object with matvec(x) method
        An N x N matrix, array, or an object with matvec(x) method to perform
        the matrix vector product A * x.  The sparse matrix formats
        in scipy.sparse are appropriate for A.
    k : integer
        The number of eigenvalues and eigenvectors desired.
        `k` must be smaller than N. It is not possible to compute all
        eigenvectors of a matrix.

    Returns
    -------
    w : array
        Array of k eigenvalues.
    v : array
        An array of `k` eigenvectors.
        ``v[:, i]`` is the eigenvector corresponding to the eigenvalue w[i].

    Other Parameters
    ----------------
    M : matrix or array
        (Not implemented)
        A symmetric positive-definite matrix for the generalized
        eigenvalue problem ``A * x = w * M * x``.
    sigma : real or complex
        (Not implemented)
        Find eigenvalues near sigma.  Shift spectrum by sigma.
    v0 : array
        Starting vector for iteration.
    ncv : integer
        The number of Lanczos vectors generated
        `ncv` must be greater than `k`; it is recommended that ``ncv > 2*k``.
    which : string
        Which `k` eigenvectors and eigenvalues to find:

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
        The default value of 0 implies machine precision.
    return_eigenvectors : boolean
        Return eigenvectors (True) in addition to eigenvalues

    Raises
    ------
    ArpackNoConvergence
        When the requested convergence is obtained.

        The currently converged eigenvalues and eigenvectors can be found
        as ``eigenvalues`` and ``eigenvectors`` attributes of the exception
        object.

    See Also
    --------
    eigsh : eigenvalues and eigenvectors for symmetric matrix A

    Notes
    -----
    This function is a wrapper to the ARPACK [1]_ SNEUPD, DNEUPD, CNEUPD,
    ZNEUPD, functions which use the Implicitly Restarted Arnoldi Method to
    find the eigenvalues and eigenvectors [2]_.

    Examples
    --------
    Find 6 eigenvectors of the identity matrix:

    >>> id = np.identity(13)
    >>> vals, vecs = sp.sparse.linalg.eigs(id, k=6)
    >>> vals
    array([ 1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j])
    >>> vecs.shape
    (13, 6)

    References
    ----------
    .. [1] ARPACK Software, http://www.caam.rice.edu/software/ARPACK/
    .. [2] R. B. Lehoucq, D. C. Sorensen, and C. Yang,  ARPACK USERS GUIDE:
       Solution of Large Scale Eigenvalue Problems by Implicitly Restarted
       Arnoldi Methods. SIAM, Philadelphia, PA, 1998.
    """
    A = _aslinearoperator_with_dtype(A)
    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix (shape=%s)' % (A.shape,))
    n = A.shape[0]

    matvec = lambda x : A.matvec(x)
    params = _UnsymmetricArpackParams(n, k, A.dtype.char, matvec, sigma,
                           ncv, v0, maxiter, which, tol)

    if M is not None:
        raise NotImplementedError("generalized eigenproblem not supported yet")

    while not params.converged:
        params.iterate()

    return params.extract(return_eigenvectors)

def eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None,
          ncv=None, maxiter=None, tol=0,
          return_eigenvectors=True):
    """
    Find k eigenvalues and eigenvectors of the real symmetric square matrix A.

    Solves A * x[i] = w[i] * x[i], the standard eigenvalue problem for
    w[i] eigenvalues with corresponding eigenvectors x[i].

    Parameters
    ----------
    A : matrix or array with real entries or object with matvec(x) method
        An N x N real symmetric matrix or array or an object with matvec(x)
        method to perform the matrix vector product A * x.  The sparse
        matrix formats in scipy.sparse are appropriate for A.
    k : integer
        The number of eigenvalues and eigenvectors desired.

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
        ncv must be greater than k and smaller than n;
        it is recommended that ncv > 2*k
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
        Relative accuracy for eigenvalues (stopping criterion).
        The default value of 0 implies machine precision.
    return_eigenvectors : boolean
        Return eigenvectors (True) in addition to eigenvalues

    Raises
    ------
    ArpackNoConvergence
        When the requested convergence is obtained.

        The currently converged eigenvalues and eigenvectors can be found
        as ``eigenvalues`` and ``eigenvectors`` attributes of the exception
        object.

    See Also
    --------
    eigs : eigenvalues and eigenvectors for a general (nonsymmetric) matrix A

    Notes
    -----
    This function is a wrapper to the ARPACK [1]_ SSEUPD and DSEUPD
    functions which use the Implicitly Restarted Lanczos Method to
    find the eigenvalues and eigenvectors [2]_.

    Examples
    --------
    >>> id = np.identity(13)
    >>> vals, vecs = sp.sparse.linalg.eigsh(id, k=6)
    >>> vals
    array([ 1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j,  1.+0.j])
    >>> vecs.shape
    (13, 6)

    References
    ----------
    .. [1] ARPACK Software, http://www.caam.rice.edu/software/ARPACK/
    .. [2] R. B. Lehoucq, D. C. Sorensen, and C. Yang,  ARPACK USERS GUIDE:
       Solution of Large Scale Eigenvalue Problems by Implicitly Restarted
       Arnoldi Methods. SIAM, Philadelphia, PA, 1998.
    """
    A = _aslinearoperator_with_dtype(A)
    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix (shape=%s)' % (A.shape,))
    n = A.shape[0]

    if M is not None:
        raise NotImplementedError("generalized eigenproblem not supported yet")

    matvec = lambda x : A.matvec(x)
    params = _SymmetricArpackParams(n, k, A.dtype.char, matvec, sigma,
                           ncv, v0, maxiter, which, tol)

    while not params.converged:
        params.iterate()

    return params.extract(return_eigenvectors)

def svds(A, k=6, ncv=None, tol=0):
    """Compute k singular values/vectors for a sparse matrix using ARPACK.

    Parameters
    ----------
    A : sparse matrix
        Array to compute the SVD on
    k : int, optional
        Number of singular values and vectors to compute.
    ncv : integer
        The number of Lanczos vectors generated
        ncv must be greater than k+1 and smaller than n;
        it is recommended that ncv > 2*k
    tol : float, optional
        Tolerance for singular values. Zero (default) means machine precision.

    Note
    ----
    This is a naive implementation using an eigensolver on A.H * A or
    A * A.H, depending on which one is more efficient.

    """
    if not (isinstance(A, np.ndarray) or isspmatrix(A)):
        A = np.asarray(A)

    n, m = A.shape

    if np.issubdtype(A.dtype, np.complexfloating):
        herm = lambda x: x.T.conjugate()
        eigensolver = eigs
    else:
        herm = lambda x: x.T
        eigensolver = eigsh

    if n > m:
        X = A
        XH = herm(A)
    else:
        XH = A
        X = herm(A)

    def matvec_XH_X(x):
        return XH.dot(X.dot(x))

    XH_X = LinearOperator(matvec=matvec_XH_X, dtype=X.dtype,
                          shape=(X.shape[1], X.shape[1]))

    eigvals, eigvec = eigensolver(XH_X, k=k, tol=tol**2)
    s = np.sqrt(eigvals)

    if n > m:
        v = eigvec
        u = X.dot(v) / s
        vh = herm(v)
    else:
        u = eigvec
        vh = herm(X.dot(u) / s)

    return u, s, vh
