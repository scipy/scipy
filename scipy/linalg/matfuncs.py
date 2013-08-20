#
# Author: Travis Oliphant, March 2002
#

from __future__ import division, print_function, absolute_import

__all__ = ['expm','expm2','expm3','cosm','sinm','tanm','coshm','sinhm',
           'tanhm','logm','funm','signm','sqrtm',
           'expm_frechet', 'fractional_matrix_power']

from numpy import asarray, Inf, dot, eye, diag, exp, \
     product, logical_not, ravel, transpose, conjugate, \
     cast, log, ogrid, imag, real, absolute, amax, sign, \
     isfinite, sqrt, single
from numpy import matrix as mat
import numpy as np


# Local imports
from .misc import norm
from .basic import solve, inv
from .lapack import ztrsyl
from .special_matrices import triu, all_mat
from .decomp import eig
from .decomp_svd import orth, svd
from .decomp_schur import schur, rsf2csf
from ._expm_frechet import expm_frechet
from ._matfuncs_sqrtm import sqrtm
import warnings

eps = np.finfo(float).eps
feps = np.finfo(single).eps


def fractional_matrix_power(A, t):
    # This fixes some issue with imports;
    # this function calls onenormest which is in scipy.sparse.
    import scipy.linalg._matfuncs_inv_ssq
    return scipy.linalg._matfuncs_inv_ssq.fractional_matrix_power(A, t)


def logm(A, disp=True):
    """
    Compute matrix logarithm.

    The matrix logarithm is the inverse of
    expm: expm(logm(`A`)) == `A`

    Parameters
    ----------
    A : (N, N) array_like
        Matrix whose logarithm to evaluate
    disp : bool, optional
        Print warning if error in the result is estimated large
        instead of returning estimated error. (Default: True)

    Returns
    -------
    logm : (N, N) ndarray
        Matrix logarithm of `A`
    errest : float
        (if disp == False)

        1-norm of the estimated error, ||err||_1 / ||A||_1

    """
    import scipy.linalg._matfuncs_inv_ssq
    A = mat(asarray(A))
    F = scipy.linalg._matfuncs_inv_ssq.logm(A)
    errtol = 1000*eps
    #TODO use a better error approximation
    errest = norm(expm(F)-A,1) / norm(A,1)
    if disp:
        if not isfinite(errest) or errest >= errtol:
            print("logm result may be inaccurate, approximate err =", errest)
        return F
    else:
        return F, errest


def expm(A, q=None):
    """
    Compute the matrix exponential using Pade approximation.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix to be exponentiated

    Returns
    -------
    expm : (N, N) ndarray
        Matrix exponential of `A`

    References
    ----------
    N. J. Higham,
    "The Scaling and Squaring Method for the Matrix Exponential Revisited",
    SIAM. J. Matrix Anal. & Appl. 26, 1179 (2005).

    """
    if q:
        warnings.warn("argument q=... in scipy.linalg.expm is deprecated.")
    import scipy.sparse.linalg
    return scipy.sparse.linalg.expm(A)


# deprecated, but probably should be left there in the long term
@np.deprecate(new_name="expm")
def expm2(A):
    """
    Compute the matrix exponential using eigenvalue decomposition.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix to be exponentiated

    Returns
    -------
    expm2 : (N, N) ndarray
        Matrix exponential of `A`

    """
    A = asarray(A)
    t = A.dtype.char
    if t not in ['f','F','d','D']:
        A = A.astype('d')
        t = 'd'
    s,vr = eig(A)
    vri = inv(vr)
    r = dot(dot(vr,diag(exp(s))),vri)
    if t in ['f', 'd']:
        return r.real.astype(t)
    else:
        return r.astype(t)


# deprecated, but probably should be left there in the long term
@np.deprecate(new_name="expm")
def expm3(A, q=20):
    """
    Compute the matrix exponential using Taylor series.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix to be exponentiated
    q : int
        Order of the Taylor series used is `q-1`

    Returns
    -------
    expm3 : (N, N) ndarray
        Matrix exponential of `A`

    """
    A = asarray(A)
    t = A.dtype.char
    if t not in ['f','F','d','D']:
        A = A.astype('d')
        t = 'd'
    A = mat(A)
    eA = eye(*A.shape,**{'dtype':t})
    trm = mat(eA, copy=True)
    castfunc = cast[t]
    for k in range(1,q):
        trm *= A / castfunc(k)
        eA += trm
    return eA

_array_precision = {'i': 1, 'l': 1, 'f': 0, 'd': 1, 'F': 0, 'D': 1}


def toreal(arr, tol=None):
    """Return as real array if imaginary part is small.

    Parameters
    ----------
    arr : array
    tol : float
        Absolute tolerance

    Returns
    -------
    arr : double or complex array
    """
    if tol is None:
        tol = {0:feps*1e3, 1:eps*1e6}[_array_precision[arr.dtype.char]]
    if (arr.dtype.char in ['F', 'D','G']) and \
       np.allclose(arr.imag, 0.0, atol=tol):
        arr = arr.real
    return arr


def cosm(A):
    """
    Compute the matrix cosine.

    This routine uses expm to compute the matrix exponentials.

    Parameters
    ----------
    A : (N, N) array_like
        Input array

    Returns
    -------
    cosm : (N, N) ndarray
        Matrix cosine of A

    """
    A = asarray(A)
    if A.dtype.char not in ['F','D','G']:
        return expm(1j*A).real
    else:
        return 0.5*(expm(1j*A) + expm(-1j*A))


def sinm(A):
    """
    Compute the matrix sine.

    This routine uses expm to compute the matrix exponentials.

    Parameters
    ----------
    A : (N, N) array_like
        Input array.

    Returns
    -------
    sinm : (N, N) ndarray
        Matrix cosine of `A`

    """
    A = asarray(A)
    if A.dtype.char not in ['F','D','G']:
        return expm(1j*A).imag
    else:
        return -0.5j*(expm(1j*A) - expm(-1j*A))


def tanm(A):
    """
    Compute the matrix tangent.

    This routine uses expm to compute the matrix exponentials.

    Parameters
    ----------
    A : (N, N) array_like
        Input array.

    Returns
    -------
    tanm : (N, N) ndarray
        Matrix tangent of `A`

    """
    A = asarray(A)
    if A.dtype.char not in ['F','D','G']:
        return toreal(solve(cosm(A), sinm(A)))
    else:
        return solve(cosm(A), sinm(A))


def coshm(A):
    """
    Compute the hyperbolic matrix cosine.

    This routine uses expm to compute the matrix exponentials.

    Parameters
    ----------
    A : (N, N) array_like
        Input array.

    Returns
    -------
    coshm : (N, N) ndarray
        Hyperbolic matrix cosine of `A`

    """
    A = asarray(A)
    if A.dtype.char not in ['F','D','G']:
        return toreal(0.5*(expm(A) + expm(-A)))
    else:
        return 0.5*(expm(A) + expm(-A))


def sinhm(A):
    """
    Compute the hyperbolic matrix sine.

    This routine uses expm to compute the matrix exponentials.

    Parameters
    ----------
    A : (N, N) array_like
        Input array.

    Returns
    -------
    sinhm : (N, N) ndarray
        Hyperbolic matrix sine of `A`

    """
    A = asarray(A)
    if A.dtype.char not in ['F','D']:
        return toreal(0.5*(expm(A) - expm(-A)))
    else:
        return 0.5*(expm(A) - expm(-A))


def tanhm(A):
    """
    Compute the hyperbolic matrix tangent.

    This routine uses expm to compute the matrix exponentials.

    Parameters
    ----------
    A : (N, N) array_like
        Input array

    Returns
    -------
    tanhm : (N, N) ndarray
        Hyperbolic matrix tangent of `A`

    """
    A = asarray(A)
    if A.dtype.char not in ['F','D']:
        return toreal(solve(coshm(A), sinhm(A)))
    else:
        return solve(coshm(A), sinhm(A))


def funm(A, func, disp=True):
    """
    Evaluate a matrix function specified by a callable.

    Returns the value of matrix-valued function ``f`` at `A`. The
    function ``f`` is an extension of the scalar-valued function `func`
    to matrices.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix at which to evaluate the function
    func : callable
        Callable object that evaluates a scalar function f.
        Must be vectorized (eg. using vectorize).
    disp : bool, optional
        Print warning if error in the result is estimated large
        instead of returning estimated error. (Default: True)

    Returns
    -------
    funm : (N, N) ndarray
        Value of the matrix function specified by func evaluated at `A`
    errest : float
        (if disp == False)

        1-norm of the estimated error, ||err||_1 / ||A||_1

    """
    # Perform Shur decomposition (lapack ?gees)
    A = asarray(A)
    if len(A.shape) != 2:
        raise ValueError("Non-matrix input to matrix function.")
    if A.dtype.char in ['F', 'D', 'G']:
        cmplx_type = 1
    else:
        cmplx_type = 0
    T, Z = schur(A)
    T, Z = rsf2csf(T,Z)
    n,n = T.shape
    F = diag(func(diag(T)))  # apply function to diagonal elements
    F = F.astype(T.dtype.char)  # e.g. when F is real but T is complex

    minden = abs(T[0,0])

    # implement Algorithm 11.1.1 from Golub and Van Loan
    #                 "matrix Computations."
    for p in range(1,n):
        for i in range(1,n-p+1):
            j = i + p
            s = T[i-1,j-1] * (F[j-1,j-1] - F[i-1,i-1])
            ksl = slice(i,j-1)
            val = dot(T[i-1,ksl],F[ksl,j-1]) - dot(F[i-1,ksl],T[ksl,j-1])
            s = s + val
            den = T[j-1,j-1] - T[i-1,i-1]
            if den != 0.0:
                s = s / den
            F[i-1,j-1] = s
            minden = min(minden,abs(den))

    F = dot(dot(Z, F),transpose(conjugate(Z)))
    if not cmplx_type:
        F = toreal(F)

    tol = {0:feps, 1:eps}[_array_precision[F.dtype.char]]
    if minden == 0.0:
        minden = tol
    err = min(1, max(tol,(tol/minden)*norm(triu(T,1),1)))
    if product(ravel(logical_not(isfinite(F))),axis=0):
        err = Inf
    if disp:
        if err > 1000*tol:
            print("funm result may be inaccurate, approximate err =", err)
        return F
    else:
        return F, err


def signm(a, disp=True):
    """
    Matrix sign function.

    Extension of the scalar sign(x) to matrices.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix at which to evaluate the sign function
    disp : bool, optional
        Print warning if error in the result is estimated large
        instead of returning estimated error. (Default: True)

    Returns
    -------
    signm : (N, N) ndarray
        Value of the sign function at `A`
    errest : float
        (if disp == False)

        1-norm of the estimated error, ||err||_1 / ||A||_1

    Examples
    --------
    >>> from scipy.linalg import signm, eigvals
    >>> a = [[1,2,3], [1,2,1], [1,1,1]]
    >>> eigvals(a)
    array([ 4.12488542+0.j, -0.76155718+0.j,  0.63667176+0.j])
    >>> eigvals(signm(a))
    array([-1.+0.j,  1.+0.j,  1.+0.j])

    """
    def rounded_sign(x):
        rx = real(x)
        if rx.dtype.char == 'f':
            c = 1e3*feps*amax(x)
        else:
            c = 1e3*eps*amax(x)
        return sign((absolute(rx) > c) * rx)
    result,errest = funm(a, rounded_sign, disp=0)
    errtol = {0:1e3*feps, 1:1e3*eps}[_array_precision[result.dtype.char]]
    if errest < errtol:
        return result

    # Handle signm of defective matrices:

    # See "E.D.Denman and J.Leyva-Ramos, Appl.Math.Comp.,
    # 8:237-250,1981" for how to improve the following (currently a
    # rather naive) iteration process:

    a = asarray(a)
    # a = result # sometimes iteration converges faster but where??

    # Shifting to avoid zero eigenvalues. How to ensure that shifting does
    # not change the spectrum too much?
    vals = svd(a,compute_uv=0)
    max_sv = np.amax(vals)
    # min_nonzero_sv = vals[(vals>max_sv*errtol).tolist().count(1)-1]
    # c = 0.5/min_nonzero_sv
    c = 0.5/max_sv
    S0 = a + c*np.identity(a.shape[0])
    prev_errest = errest
    for i in range(100):
        iS0 = inv(S0)
        S0 = 0.5*(S0 + iS0)
        Pp = 0.5*(dot(S0,S0)+S0)
        errest = norm(dot(Pp,Pp)-Pp,1)
        if errest < errtol or prev_errest == errest:
            break
        prev_errest = errest
    if disp:
        if not isfinite(errest) or errest >= errtol:
            print("signm result may be inaccurate, approximate err =", errest)
        return S0
    else:
        return S0, errest
