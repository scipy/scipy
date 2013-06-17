#
# Author: Travis Oliphant, March 2002
#

from __future__ import division, print_function, absolute_import

__all__ = ['expm','expm2','expm3','cosm','sinm','tanm','coshm','sinhm',
           'tanhm','logm','funm','signm','sqrtm',
           'expm_frechet']

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
import warnings

eps = np.finfo(float).eps
feps = np.finfo(single).eps


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


def expm3(A, q=20):
    """
    Compute the matrix exponential using Taylor series.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix to be exponentiated
    q : int
        Order of the Taylor series

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
            print("Result may be inaccurate, approximate err =", err)
        return F
    else:
        return F, err


def _almost_diagonal(M, abstol=1e-13):
    M_abs_off_diagonal = np.absolute(M - np.diag(np.diag(M)))
    return np.max(M_abs_off_diagonal) < abstol


def _inverse_squaring_helper(T0, theta):
    """
    A helper function for inverse scaling and squaring for Pade approximation.

    Parameters
    ----------
    T0 : (N, N) array_like upper triangular
        Matrix involved in inverse scaling and squaring.
    theta : indexable
        The values theta[1] .. theta[7] must be available.
        They represent bounds related to Pade approximation, and they depend
        on the matrix function which is being computed.
        For example, different values of theta are required for
        matrix logarithm than for matrix fractional power.

    Returns
    -------
    T : (N, N) array_like upper triangular
        Composition of zero or more matrix square roots of T0.
    s : non-negative integer
        Number of square roots taken.
    m : positive integer
        The degree of the Pade approximation.

    Notes
    -----
    This subroutine appears as a chunk of lines within
    a couple of published algorithms; for example it appears
    as lines 4--35 in algorithm (3.1) of [1], and
    as lines 3--34 in algorithm (4.1) of [2].

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing Lin (2013)
           "An Improved Schur-Pade Algorithm for Fractional Powers
           of a Matrix and their Frechet Derivatives."

    .. [2] Awad H. Al-Mohy and Nicholas J. Higham (2012)
           "Improved Inverse Scaling and Squaring Algorithms
           for the Matrix Logarithm."
           SIAM Journal on Scientific Computing, 34 (4). C152-C169.
           ISSN 1095-7197

    """
    T = T0

    # Find s0, the smallest s such that the spectral radius
    # of a certain diagonal matrix is at most theta[7].
    s0 = 0
    tmp_diag = np.diag(T)
    while np.max(np.absolute(tmp_diag - 1)) > theta[7]:
        tmp_diag = np.sqrt(tmp_diag)
        s0 += 1

    # Take matrix square roots of T.
    for i in range(s0):
        T = _sqrtm_triu(T)

    # Flow control in this section is a little odd.
    # This is because I am translating algorithm descriptions
    # which have GOTOs in the publication.
    s = s0
    k = 0
    d2 = onenormest_m1_power(T, 2) ** (1/2)
    d3 = onenormest_m1_power(T, 3) ** (1/3)
    a2 = max(d2, d3)
    m = None
    for i in (1, 2):
        if a2 <= theta[i]:
            m = i
            break
    while m is None:
        if s > s0:
            d3 = onenormest_m1_power(T, 3) ** (1/3)
        d4 = onenormest_m1_power(T, 4) ** (1/4)
        a3 = max(d3, d4)
        if a3 <= theta[7]:
            j1 = min(i for i in (3, 4, 5, 6, 7) if a3 <= theta[i])
            if j1 <= 6:
                m = j1
                break
            elif a3 / 2 <= theta[5] and k < 2:
                k += 1
                T = _sqrtm_triu(T)
                s += 1
                continue
        d5 = onenormest_m1_power(T, 5) ** (1/5)
        a4 = max(d4, d5)
        eta = min(a3, a4)
        for i in (6, 7):
            if eta <= theta[i]:
                m = i
                break
        if m is not None:
            break
        T = _sqrtm_triu(T)
        s += 1

    # Return the root matrix, the number of square roots, and the Pade degree.
    return T, s, m


def fractional_power(A, t):
    """
    Compute a fractional power of a matrix.

    This uses algorithm (3.1) of [1].

    Parameters
    ----------
    A : (N, N) array_like
        Matrix whose fractional power to evaluate.
    t : float
        Fractional power between -1 and 1 exclusive.

    Returns
    -------
    fractional_power : (N, N) array_like
        The fractional power of the matrix.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing Lin (2013)
           "An Improved Schur-Pade Algorithm for Fractional Powers
           of a Matrix and their Frechet Derivatives."

    """
    A = asarray(A)
    if len(A.shape) != 2:
        raise ValueError("Non-matrix input to matrix function.")
    T, Z = schur(A)
    T, Z = rsf2csf(T,Z)
    m_to_theta = {
            1 : 1.51e-5,
            2 : 2.24e-3,
            3 : 1.88e-2,
            4 : 6.04e-2,
            5 : 1.24e-1,
            6 : 2.00e-1,
            7 : 2.79e-1,
            #8 : 3.55e-1,
            #9 : 4.25e-1,
            #10 : 4.87e-1,
            #11 : 5.42e-1,
            #12 : 5.90e-1,
            #13 : 6.32e-1,
            #14 : 6.69e-1,
            #15 : 7.00e-1,
            #16 : 7.28e-1,
            #32 : 9.15e-1,
            #64 : 9.76e-1,
            }
    if _almost_diagonal(T):
        U = np.diag(np.diag(T) ** t)
        U, Z = all_mat(U, Z)
        X = (Z * U * Z.H)
        return X.A
    else:
        T0 = T
        T, s, m = _inverse_squaring_helper(T, theta)


    U, Z = all_mat(U, Z)
    X = (Z * U * Z.H)
    return X.A



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
    # Compute using general funm but then use better error estimator and
    #   make one step in improving estimate using a rotation matrix.
    A = mat(asarray(A))
    F, errest = funm(A,log,disp=0)
    errtol = 1000*eps
    # Only iterate if estimate of error is too large.
    if errest >= errtol:
        # Use better approximation of error
        errest = norm(expm(F)-A,1) / norm(A,1)
        if not isfinite(errest) or errest >= errtol:
            N,N = A.shape
            X,Y = ogrid[1:N+1,1:N+1]
            R = mat(orth(eye(N,dtype='d')+X+Y))
            F, dontcare = funm(R*A*R.H,log,disp=0)
            F = R.H*F*R
            if (norm(imag(F),1) <= 1000*errtol*norm(F,1)):
                F = mat(real(F))
            E = mat(expm(F))
            temp = mat(solve(E.T,(E-A).T))
            F = F - temp.T
            errest = norm(expm(F)-A,1) / norm(A,1)
    if disp:
        if not isfinite(errest) or errest >= errtol:
            print("Result may be inaccurate, approximate err =", errest)
        return F
    else:
        return F, errest


def logm_new(A):
    """
    Compute matrix logarithm.

    The matrix logarithm is the inverse of
    expm: expm(logm(`A`)) == `A`

    Parameters
    ----------
    A : (N, N) array_like
        Matrix whose logarithm to evaluate

    Returns
    -------
    logm : (N, N) ndarray
        Matrix logarithm of `A`

    References
    ----------
    .. [1] Awad H. Al-Mohy and Nicholas J. Higham (2012)
           "Improved Inverse Scaling and Squaring Algorithms
           for the Matrix Logarithm."
           SIAM Journal on Scientific Computing, 34 (4). C152-C169.
           ISSN 1095-7197

    .. [2] Nicholas J. Higham (2008)
           "Functions of Matrices: Theory and Computation"
           ISBN 978-0-898716-46-7

    .. [3] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    """
    A = asarray(A)
    if len(A.shape) != 2:
        raise ValueError("Non-matrix input to matrix function.")
    T, Z = schur(A)
    T, Z = rsf2csf(T,Z)
    T0 = T

    # Define bounds given in Table (2.1).
    theta = (None,
            1.59e-5, 2.31e-3, 1.94e-2, 6.21e-2,
            1.28e-1, 2.06e-1, 2.88e-1, 3.67e-1,
            4.39e-1, 5.03e-1, 5.60e-1, 6.09e-1,
            6.52e-1, 6.89e-1, 7.21e-1, 7.49e-1)

    T, s, m = _inverse_squaring_helper(T, theta)

    #TODO line 35
    #NOTE I am starting to think that fractional matrix powers 
    #NOTE should be implemented before implementing the matrix logarithm.
    # Replace diag(T - I) by diag(T0)^(1/(2^s)) - 1
    # using [1, Alg. 2] and recompute the first superdiagonal of T
    # using [18 (my [3]), eq. (5, 6)] applied to T0.

    #TODO line 36
    # Evaluate U = 2**s r_m(T - I) using the partial fraction expansion (1.1).

    # Recompute diagonal entries of U.
    U[np.diag_indices(n)] = np.log(np.diag(T0))

    # Recompute superdiagonal entries of U.
    # This indexing of this code should be renovated
    # when newer np.diagonal() becomes available.
    for i in range(n-1):
        l1 = T0[i, i]
        l2 = T0[i+1, i+1]
        t12 = T0[i, i+1]
        U[i, i+1] = _logm_superdiag_entry(l1, l2, t12)

    U, Z = all_mat(U, Z)
    X = (Z * U * Z.H)
    return X.A


def _fractional_power_superdiag_entry(l1, l2, t12, p):
    """
    Compute a superdiagonal entry of a matrix fractional power.

    This is Eq. (5.6) in [1].

    Parameters
    ----------
    l1 : complex
        A diagonal entry of the matrix.
    l2 : complex
        A diagonal entry of the matrix.
    t12 : complex
        A superdiagonal entry of the matrix.
    p : float
        A fractional power.

    Returns
    -------
    f12 : complex
        A superdiagonal entry of the matrix fractional power.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    """
    if l1 == l2:
        return t12 * p * l1**(p-1)
    elif abs(l1) < abs(l2) / 2 or abs(l2) < abs(l1) / 2:
        return t12 * ((l2**p) - (l1**p)) / (l2 - l1)
    else:
        # This is Eq. (5.5) in [1].
        z = (l2 - l1) / (l2 + l1)
        log_l1 = np.log(l1)
        log_l2 = np.log(l2)
        atanh_z = np.atanh(z)
        tmp_a = t12 * np.exp((p/2)*(log_l2 + log_l1))
        tmp_b = p * (atanh_z + np.pi * 1j * _unwindk(log_l2 - log_l1))
        tmp_c = 2 * np.sinh(tmp_b) / (l2 - l1)
        return tmp_a * tmp_c


def _logm_superdiag_entry(l1, l2, t12):
    """
    Compute a superdiagonal entry of a matrix logarithm.

    This is Eq. (11.28) in [1].

    Parameters
    ----------
    l1 : complex
        A diagonal entry of the matrix.
    l2 : complex
        A diagonal entry of the matrix.
    t12 : complex
        A superdiagonal entry of the matrix.

    Returns
    -------
    f12 : complex
        A superdiagonal entry of the matrix logarithm.

    References
    ----------
    .. [1] Nicholas J. Higham (2008)
           "Functions of Matrices: Theory and Computation"
           ISBN 978-0-898716-46-7

    """
    if l1 == l2:
        return t12 / l1
    elif abs(l1) < abs(l2) / 2 or abs(l2) < abs(l1) / 2:
        return t12 * (np.log(l2) - np.log(l1)) / (l2 - l1)
    else:
        z = (l2 - l1) / (l2 + l1)
        ua = _unwindk(np.log(l2) - np.log(l1))
        ub = _unwindk(np.log(1+z) - np.log(1-z))
        return t12 * (2*np.atanh(z) + 2*np.pi*1j*(ua + ub)) / (l2 - l1)


def _unwindk(z):
    """
    Compute the scalar unwinding number.

    Uses Eq. (5.3) in [1], and should be equal to (z - log(exp(z)) / (2 pi i).
    Note that this definition differs in sign from the original definition
    in equations (5, 6) in [2].  The sign convention is justified in [3].

    Parameters
    ----------
    z : complex
        A complex number.

    Returns
    -------
    unwinding_number : integer
        The scalar unwinding number of z.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    .. [2] Robert M. Corless and David J. Jeffrey,
           "The unwinding number." Newsletter ACM SIGSAM Bulletin
           Volume 30, Issue 2, June 1996, Pages 28-35.

    .. [3] Russell Bradford and Robert M. Corless and James H. Davenport and
           David J. Jeffrey and Stephen M. Watt,
           "Reasoning about the elementary functions of complex analysis"
           Annals of Mathematics and Artificial Intelligence,
           36: 303-318, 2002.

    """
    return int(np.ceil((z.imag - np.pi) / (2*np.pi)))


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
            print("Result may be inaccurate, approximate err =", errest)
        return S0
    else:
        return S0, errest


def _sqrtm_triu(T, blocksize=64):
    """
    Matrix square root of an upper triangular matrix.

    This is a helper function for `sqrtm` and `logm`.

    Parameters
    ----------
    T : (N, N) array_like upper triangular
        Matrix whose square root to evaluate
    blocksize : integer, optional
        If the blocksize is not degenerate with respect to the
        size of the input array, then use a blocked algorithm. (Default: 64)

    Returns
    -------
    sqrtm : (N, N) ndarray
        Value of the sqrt function at `T`

    References
    ----------
    .. [1] Edvin Deadman, Nicholas J. Higham, Rui Ralha (2013)
           "Blocked Schur Algorithms for Computing the Matrix Square Root,
           Lecture Notes in Computer Science, 7782. pp. 171-182.

    """
    n, n = T.shape
    R = np.zeros((n,n), T.dtype.char)
    
    # Take square roots of the diagonal of the triangular matrix.
    R[np.diag_indices(n)] = np.sqrt(np.diag(T))

    # Compute the number of blocks to use; use at least one block.
    nblocks = max(n // blocksize, 1)

    # Compute the smaller of the two sizes of blocks that
    # we will actually use, and compute the number of large blocks.
    bsmall, nlarge = divmod(n, nblocks)
    blarge = bsmall + 1
    nsmall = nblocks - nlarge
    if nsmall * bsmall + nlarge * blarge != n:
        raise Exception('internal inconsistency')

    # Define the index range covered by each block.
    start_stop_pairs = []
    start = 0
    for count, size in ((nsmall, bsmall), (nlarge, blarge)):
        for i in range(count):
            start_stop_pairs.append((start, start + size))
            start += size

    # Within-block interactions.
    for start, stop in start_stop_pairs:
        for j in range(start, stop):
            for i in range(j-1, start-1, -1):
                s = 0
                if j - i > 1:
                    s = R[i, i+1:j].dot(R[i+1:j, j])
                R[i,j] = (T[i,j] - s)/(R[i,i] + R[j,j])

    # Between-block interactions.
    for j in range(nblocks):
        jstart, jstop = start_stop_pairs[j]
        for i in range(j-1, -1, -1):
            istart, istop = start_stop_pairs[i]
            S = T[istart:istop, jstart:jstop]
            if j - i > 1:
                S = S - R[istart:istop, istop:jstart].dot(
                        R[istop:jstart, jstart:jstop])

            # Invoke LAPACK.
            # For more details, see the solve_sylvester implemention
            # and the fortran ztrsyl docs.
            Rii = R[istart:istop, istart:istop]
            Rjj = R[jstart:jstop, jstart:jstop]
            x, scale, info = ztrsyl(Rii, Rjj, S)
            R[istart:istop, jstart:jstop] = x * scale
    
    # Return the matrix square root.
    return R


def sqrtm(A, disp=True, blocksize=64):
    """
    Matrix square root.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix whose square root to evaluate
    disp : bool, optional
        Print warning if error in the result is estimated large
        instead of returning estimated error. (Default: True)
    blocksize : integer, optional
        If the blocksize is not degenerate with respect to the
        size of the input array, then use a blocked algorithm. (Default: 64)

    Returns
    -------
    sqrtm : (N, N) ndarray
        Value of the sqrt function at `A`

    errest : float
        (if disp == False)

        Frobenius norm of the estimated error, ||err||_F / ||A||_F

    References
    ----------
    .. [1] Edvin Deadman, Nicholas J. Higham, Rui Ralha (2013)
           "Blocked Schur Algorithms for Computing the Matrix Square Root,
           Lecture Notes in Computer Science, 7782. pp. 171-182.

    """
    A = asarray(A)
    if len(A.shape) != 2:
        raise ValueError("Non-matrix input to matrix function.")
    if blocksize < 1:
        raise ValueError("The blocksize should be at least 1.")
    T, Z = schur(A)
    T, Z = rsf2csf(T,Z)
    R = _sqrtm_triu(T, blocksize=blocksize)
    R, Z = all_mat(R,Z)
    X = (Z * R * Z.H)

    if disp:
        nzeig = np.any(diag(T) == 0)
        if nzeig:
            print("Matrix is singular and may not have a square root.")
        return X.A
    else:
        arg2 = norm(X*X - A,'fro')**2 / norm(A,'fro')
        return X.A, arg2
