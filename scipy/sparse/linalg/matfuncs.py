"""
Sparse matrix functions
"""

#
# Authors: Travis Oliphant, March 2002
#          Anthony Scopatz, August 2012 (Sparse Updates)
#          Jake Vanderplas, August 2012 (Sparse Updates)
#

from __future__ import division, print_function, absolute_import

__all__ = ['expm', 'inv']

import math

from numpy import asarray, dot, eye, ceil, log2
from numpy import matrix as mat
import numpy as np

import scipy.linalg
import scipy.special
from scipy.linalg.misc import norm
from scipy.linalg.basic import solve, inv

import scipy.sparse
from scipy.sparse.base import isspmatrix
from scipy.sparse.construct import eye as speye
from scipy.sparse.linalg import spsolve


def inv(A):
    """
    Compute the inverse of a sparse matrix

    .. versionadded:: 0.12.0

    Parameters
    ----------
    A : (M,M) ndarray or sparse matrix
        square matrix to be inverted

    Returns
    -------
    Ainv : (M,M) ndarray or sparse matrix
        inverse of `A`

    Notes
    -----
    This computes the sparse inverse of `A`.  If the inverse of `A` is expected
    to be non-sparse, it will likely be faster to convert `A` to dense and use
    scipy.linalg.inv.

    """
    I = speye(A.shape[0], A.shape[1], dtype=A.dtype, format=A.format)
    Ainv = spsolve(A, I)
    return Ainv


def expm(A):
    """
    Compute the matrix exponential using Pade approximation.

    .. versionadded:: 0.12.0

    Parameters
    ----------
    A : (M,M) array or sparse matrix
        2D Array or Matrix (sparse or dense) to be exponentiated

    Returns
    -------
    expA : (M,M) ndarray
        Matrix exponential of `A`

    References
    ----------
    N. J. Higham,
    "The Scaling and Squaring Method for the Matrix Exponential Revisited",
    SIAM. J. Matrix Anal. & Appl. 26, 1179 (2005).

    """
    n_squarings = 0
    Aissparse = isspmatrix(A)

    if Aissparse:
        A_L1 = max(abs(A).sum(axis=0).flat)
        ident = speye(A.shape[0], A.shape[1], dtype=A.dtype, format=A.format)
    else:
        A = asarray(A)
        A_L1 = norm(A,1)
        ident = eye(A.shape[0], A.shape[1], dtype=A.dtype)

    if A.dtype == 'float64' or A.dtype == 'complex128':
        if A_L1 < 1.495585217958292e-002:
            U,V = _pade3(A, ident)
        elif A_L1 < 2.539398330063230e-001:
            U,V = _pade5(A, ident)
        elif A_L1 < 9.504178996162932e-001:
            U,V = _pade7(A, ident)
        elif A_L1 < 2.097847961257068e+000:
            U,V = _pade9(A, ident)
        else:
            maxnorm = 5.371920351148152
            n_squarings = max(0, int(ceil(log2(A_L1 / maxnorm))))
            A = A / 2**n_squarings
            U,V = _pade13(A, ident)
    elif A.dtype == 'float32' or A.dtype == 'complex64':
        if A_L1 < 4.258730016922831e-001:
            U,V = _pade3(A, ident)
        elif A_L1 < 1.880152677804762e+000:
            U,V = _pade5(A, ident)
        else:
            maxnorm = 3.925724783138660
            n_squarings = max(0, int(ceil(log2(A_L1 / maxnorm))))
            A = A / 2**n_squarings
            U,V = _pade7(A, ident)
    else:
        raise ValueError("invalid type: "+str(A.dtype))

    P = U + V  # p_m(A) : numerator
    Q = -U + V  # q_m(A) : denominator

    if Aissparse:
        R = spsolve(Q, P)
    else:
        R = solve(Q,P)

    # squaring step to undo scaling
    for i in range(n_squarings):
        R = R.dot(R)

    return R


def expm_2009(A):
    """
    Compute the matrix exponential using Pade approximation.

    Parameters
    ----------
    A : (M,M) array or sparse matrix
        2D Array or Matrix (sparse or dense) to be exponentiated

    Returns
    -------
    expA : (M,M) ndarray
        Matrix exponential of `A`

    Notes
    -----
    This is algorithm (6.1) which is a simplification of algorithm (5.1).

    References
    ----------
    .. [1] Awad H. Al-Mohy and Nicholas J. Higham (2009)
           "A New Scaling and Squaring Algorithm for th Matrix Exponential."
           SIAM Journal on Matrix Analysis and Applications.
           31 (3). pp. 970-989. ISSN 1095-7162

    """

    #XXX This function intends to use a fast norm estimate,
    #XXX but because the fast norm estimation code has not yet made its
    #XXX way into scipy, we are using slow exact norm calculations.

    # Define the identity matrix depending on sparsity.
    if isspmatrix(A):
        ident = speye(A.shape[0], A.shape[1], dtype=A.dtype, format=A.format)
    else:
        A = asarray(A)
        ident = eye(A.shape[0], A.shape[1], dtype=A.dtype)

    # Try Pade order 3.
    A2 = A.dot(A)
    d6 = _onenorm_of_power(A2, 3)**(1/6.)
    eta_1 = max(_onenorm_of_power(A2, 2) **(1/4.), d6)
    if eta_1 < 1.495585217958292e-002 and _ell(A, 3) == 0:
        U, V = _pade3(A, ident, A2)
        return _solve_P_Q(U, V)

    # Try Pade order 5.
    A4 = A2.dot(A2)
    d4 = scipy.linalg.norm(A4, 1)**(1/4.)
    eta_2 = max(d4, d6)
    if eta_2 < 2.539398330063230e-001 and _ell(A, 5) == 0:
        U, V = _pade5(A, ident, A2, A4)
        return _solve_P_Q(U, V)

    # Try Pade orders 7 and 9.
    A6 = A2.dot(A4)
    d6 = scipy.linalg.norm(A6, 1)**(1/6.)
    d8 = _onenorm_of_power(A4, 2)**(1/8.)
    eta_3 = max(d6, d8)
    if eta_3 < 9.504178996162932e-001 and _ell(A, 7) == 0:
        U, V = _pade7(A, ident, A2, A4, A6)
        return _solve_P_Q(U, V)
    if eta_3 < 2.097847961257068e+000 and _ell(A, 9) == 0:
        U, V = _pade9(A, ident, A2, A4, A6)
        return _solve_P_Q(U, V)

    # Use Pade order 13.
    d10 = _onenorm_of_product(A4, A6)**(1/10.)
    eta_4 = max(d8, d10)
    eta_5 = min(eta_3, eta_4)
    theta_13 = 4.25
    s = max(int(math.ceil(math.log(eta_5 / theta_13, 2))), 0)
    s = s + _ell(2**-s * A, 13)
    A = A * 2**-s
    A2 = A2 * 2**(-2*s)
    A4 = A4 * 2**(-4*s)
    A6 = A6 * 2**(-6*s)
    U, V = _pade13(A, ident, A2, A4, A6)
    X = _solve_P_Q(U, V)
    if _is_upper_triangular(A):
        # Invoke Code Fragment 2.1.
        _fragment_2_1(X, A, s)
    else:
        # X = r_13(A)^(2^s) by repeated squaring.
        for i in range(s):
            X = X.dot(X)
    return X


def _is_upper_triangular(A):
    """
    A helper function for expm_2009.
    """
    if isspmatrix(A):
        return _sparse_count_nonzero(scipy.sparse.tril(A, -1)) == 0
    else:
        return np.count_nonzero(np.tril(A, -1)) == 0


def _sparse_count_nonzero(A):
    # XXX this is obviously not the right way to do this...
    return np.count_nonzero(A.todense())


def _solve_P_Q(U, V):
    """
    A helper function for expm_2009.
    """
    P = U + V
    Q = -U + V
    if isspmatrix(U):
        R = spsolve(Q, P)
    else:
        R = solve(Q, P)
    return R


def _sinch(x):
    """
    Stably evaluate sinch.

    Notes
    -----
    The strategy of falling back to a sixth order Taylor expansion
    was suggested by the Spallation Neutron Source docs
    which was found on the internet by google search.
    http://www.ornl.gov/~t6p/resources/xal/javadoc/gov/sns/tools/math/ElementaryFunction.html
    The details of the cutoff point and the Horner-like evaluation
    was picked without reference to anything in particular.

    Note that sinch is not currently implemented in scipy.special,
    whereas the "engineer's" definition of sinc is implemented.
    The implementation of sinc involves a scaling factor of pi
    that distinguishes it from the "mathematician's" version of sinc.

    """

    # If x is small then use sixth order Taylor expansion.
    # How small is small? I am using the point where the relative error
    # of the approximation is less than 1e-14.
    # If x is large then directly evaluate sinh(x) / x.
    x2 = x*x
    if abs(x) < 0.0135:
        return 1 + (x2/6.)*(1 + (x2/20.)*(1 + (x2/42.)))
    else:
        return scipy.special.sinh(x) / x


def _eq_10_42(lam_1, lam_2, t_12):
    """
    Equation (10.42) of Functions of Matrices: Theory and Computation.

    Notes
    -----
    This is a helper function for _fragment_2_1 of expm_2009.
    Equation (10.42) is on page 251 in the section on Schur algorithms.
    In particular, section 10.4.3 explains the Schur-Parlett algorithm.
    expm([[lam_1, t_12], [0, lam_1])
    =
    [[exp(lam_1), t_12*exp((lam_1 + lam_2)/2)*sinch((lam_1 - lam_2)/2)],
    [0, exp(lam_2)]
    """

    # The plain formula t_12 * (exp(lam_2) - exp(lam_2)) / (lam_2 - lam_1)
    # apparently suffers from cancellation, according to Higham's textbook.
    # A nice implementation of sinch, defined as sinh(x)/x,
    # will apparently work around the cancellation.
    a = 0.5 * (lam_1 + lam_2)
    b = 0.5 * (lam_1 - lam_2)
    return t_12 * math.exp(a) * _sinch(b)


def _fragment_2_1(X, T, s):
    """
    A helper function for expm_2009.
    """
    # Form X = r_m(2^-s T)
    # Replace diag(X) by exp(2^-s diag(T)).
    n = X.shape[0]
    diag_T = np.diag(T)

    # XXX Can this chunk be vectorized?
    scale = 2 ** -s
    exp_diag = np.exp(scale * diag_T)
    for k in range(n):
        X[k, k] = exp_diag[k]

    for i in range(s-1, -1, -1):
        scale = 2 ** -i
        X = X.dot(X)

        # Replace diag(X) by exp(2^-i diag(T)).
        # XXX Can this chunk be vectorized?
        exp_diag = np.exp(scale * diag_T)
        for k in range(n):
            X[k, k] = exp_diag[k]

        # Replace (first) superdiagonal of X by explicit formula
        # for superdiagonal of exp(2^-i T) from Eq (10.42) of 
        # the author's 2008 textbook
        # Functions of Matrices: Theory and Computation.
        # XXX Can this chunk be vectorized?
        for k in range(n-1):
            lam_1 = scale * diag_T[k]
            lam_2 = scale * diag_T[k+1]
            t_12 = scale * T[k, k+1]
            X[k, k+1] = _eq_10_42(lam_1, lam_2, t_12)


def _ell(A, m):
    """
    A helper function for expm_2009.

    Parameters
    ----------
    A : linear operator
        A linear operator whose norm of power we care about.
    m : int
        The power of the linear operator

    Returns
    -------
    value : int
        A value related to a bound.

    """
    power = 2*m + 1

    # The c_i are explained in (2.2) and (2.6) of the 2005 expm paper.
    # They are coefficients of terms of a generating function series expansion.
    c_numerator = pow(-1, power)*math.factorial(power)*math.factorial(power)
    c_denominator = math.factorial(power*2)*math.factorial(power*2+1)
    c_reciprocal = c_denominator / c_numerator
    c = 1 / float(c_reciprocal)

    # This is explained after Eq. (1.2) of the 2009 expm paper.
    # It is the "unit roundoff" of IEEE double precision arithmetic.
    u = 2**-53

    # This should be computed using an estimate,
    # but for now we are computing it exactly and inefficiently.
    est = _onenorm_of_power(abs(A), power)

    alpha = abs(c) * est / scipy.linalg.norm(A, 1)
    return max(int(math.ceil(math.log(alpha/u, 2) / (2 * m))), 0)


def _wrapped_onenorm(A):
    """
    Compute the 1-norm of a linear operator.

    Parameters
    ----------
    A : linear operator
        Linear operator whose 1-norm is requested.

    Returns
    -------
    value : float
        The 1-norm of the linear operator.

    Notes
    -----
    It would be nice if this function did not have to exist.

    """
    if isspmatrix(A):
        return max(abs(A).sum(axis=0).flat)
    else:
        return scipy.linalg.norm(asarray(A), 1)


def _onenorm_of_product(*args):
    """
    Compute the 1-norm of a product of linear operators.

    Parameters
    ----------
    args : sequence of linear operators
        This is the sequence of operators whose 1-norm of product
        is to be computed.

    Returns
    -------
    value : float
        The 1-norm of the product of the operators.

    Notes
    -----
    This function is inexcusably inefficient for anything but testing.
    It is a placeholder for a more efficient 1-norm estimation function.

    """
    if not args:
        return 1.0
    M_product = args[0]
    for M in args[1:]:
        M_product = M_product.dot(M)
    return _wrapped_onenorm(M_product)


def _onenorm_of_power(M, power):
    """
    Compute the 1-norm of a power of a linear operator.

    Parameters
    ----------
    M : linear operator
        The linear operator whose 1-norm of a power is to be computed.
    power : nonnegative int
        The requested power of the linear operator.

    Returns
    -------
    value : float
        The 1-norm of the power of the operator.

    Notes
    -----
    This function is extremely inefficient.
    It is a placeholder for a more efficient 1-norm estimation function.

    """
    if not power:
        return 1.0
    if int(power) != power or power < 1:
        raise ValueError('expected a nonnegative integer power')
    M_power = np.linalg.matrix_power(M, power)
    return _wrapped_onenorm(M_power)



# Implementation of Pade approximations of various degree
# using the algorithm presented in [Higham 2005].
# These should apply to both dense and sparse matricies.
# ident is the identity matrix, which matches A in being sparse or dense.
def _pade3(A, ident, A2=None):
    b = (120., 60., 12., 1.)
    if A2 is None:
        A2 = A.dot(A)
    U = A.dot(b[3]*A2 + b[1]*ident)
    V = b[2]*A2 + b[0]*ident
    return U,V

def _pade5(A, ident, A2=None, A4=None):
    b = (30240., 15120., 3360., 420., 30., 1.)
    if A2 is None:
        A2 = A.dot(A)
    if A4 is None:
        A4 = A2.dot(A2)
    U = A.dot(b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade7(A, ident, A2=None, A4=None, A6=None):
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    if A2 is None:
        A2 = A.dot(A)
    if A4 is None:
        A4 = A2.dot(A2)
    if A6 is None:
        A6 = A4.dot(A2)
    U = A.dot(b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade9(A, ident, A2=None, A4=None, A6=None):
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
                2162160., 110880., 3960., 90., 1.)
    if A2 is None:
        A2 = A.dot(A)
    if A4 is None:
        A4 = A2.dot(A2)
    if A6 is None:
        A6 = A4.dot(A2)
    A8 = A6.dot(A2)
    U = A.dot(b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade13(A, ident, A2=None, A4=None, A6=None):
    b = (64764752532480000., 32382376266240000., 7771770303897600.,
    1187353796428800., 129060195264000., 10559470521600., 670442572800.,
    33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.)
    if A2 is None:
        A2 = A.dot(A)
    if A4 is None:
        A4 = A2.dot(A2)
    if A6 is None:
        A6 = A4.dot(A2)
    U = A.dot(A6.dot(b[13]*A6 + b[11]*A4 + b[9]*A2) + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = A6.dot(b[12]*A6 + b[10]*A4 + b[8]*A2) + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V
