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

import scipy.misc
from scipy.linalg.misc import norm
from scipy.linalg.basic import solve, solve_triangular, inv

from scipy.sparse.base import isspmatrix
from scipy.sparse.construct import eye as speye
from scipy.sparse.linalg import spsolve

import scipy.sparse
import scipy.sparse.linalg
from scipy.sparse.linalg.interface import LinearOperator

UPPER_TRIANGULAR = 'upper_triangular'


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


def _exact_1_norm(A):
    # A compatibility function which should eventually disappear.
    # This is copypasted from expm_action.
    if scipy.sparse.isspmatrix(A):
        return max(abs(A).sum(axis=0).flat)
    else:
        return np.linalg.norm(A, 1)


def _ident_like(A):
    # A compatibility function which should eventually disappear.
    # This is copypasted from expm_action.
    if scipy.sparse.isspmatrix(A):
        return scipy.sparse.construct.eye(A.shape[0], A.shape[1],
                dtype=A.dtype, format=A.format)
    else:
        return np.eye(A.shape[0], A.shape[1], dtype=A.dtype)


def _count_nonzero(A):
    # A compatibility function which should eventually disappear.
    #XXX There should be a better way to do this when A is sparse
    #    in the traditional sense.
    if isspmatrix(A):
        return np.sum(A.toarray() != 0)
    else:
        return np.sum(A != 0)


def _is_upper_triangular(A):
    # This function could possibly be of wider interest.
    if isspmatrix(A):
        lower_part = scipy.sparse.tril(A, -1)
        if lower_part.nnz == 0:
            # structural upper triangularity
            return True
        else:
            # coincidental upper triangularity
            return _count_nonzero(lower_part) == 0
    else:
        return _count_nonzero(np.tril(A, -1)) == 0


class MatrixPowerOperator(LinearOperator):

    def __init__(self, A, p):
        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ValueError('expected A to be like a square matrix')
        if p < 0:
            raise ValueError('expected p to be a non-negative integer')
        self._A = A
        self._p = p
        self.ndim = A.ndim
        self.shape = A.shape

    def matvec(self, x):
        for i in range(self._p):
            x = self._A.dot(x)
        return x

    def rmatvec(self, x):
        for i in range(self._p):
            x = x.dot(self._A)
        return x

    def matmat(self, X):
        for i in range(self._p):
            X = self._A.dot(X)
        return X

    @property
    def T(self):
        return MatrixPowerOperator(self._A.T, self._p)


class ProductOperator(LinearOperator):
    """
    For now, this is limited to products of multiple square matrices.
    """

    def __init__(self, *args):
        for A in args:
            if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
                raise ValueError(
                        'For now, the ProductOperator implementation is '
                        'limited to the product of multiple square matrices.')
        if args:
            n = args[0].shape[0]
            for A in args:
                for d in A.shape:
                    if d != n:
                        raise ValueError(
                                'The square matrices of the ProductOperator '
                                'must all have the same shape.')
            self.shape = (n, n)
            self.ndim = len(self.shape)
        self._operator_sequence = args

    def matvec(self, x):
        for A in reversed(self._operator_sequence):
            x = A.dot(x)
        return x

    def rmatvec(self, x):
        for A in self._operator_sequence:
            x = x.dot(A)
        return x

    def matmat(self, X):
        for A in reversed(self._operator_sequence):
            X = A.dot(X)
        return X

    @property
    def T(self):
        T_args = [A.T for A in reversed(self._operator_sequence)]
        return ProductOperator(*T_args)


def _onenormest_matrix_power(A, p,
        t=2, itmax=5, compute_v=False, compute_w=False):
    """
    Efficiently estimate the 1-norm of A^p.

    Parameters
    ----------
    A : ndarray
        Matrix whose 1-norm of a power is to be computed.
    p : int
        Non-negative integer power.
    t : int, optional
        A positive parameter controlling the tradeoff between
        accuracy versus time and memory usage.
        Larger values take longer and use more memory
        but give more accurate output.
    itmax : int, optional
        Use at most this many iterations.
    compute_v : bool, optional
        Request a norm-maximizing linear operator input vector if True.
    compute_w : bool, optional
        Request a norm-maximizing linear operator output vector if True.

    Returns
    -------
    est : float
        An underestimate of the 1-norm of the sparse matrix.
    v : ndarray, optional
        The vector such that ||Av||_1 == est*||v||_1.
        It can be thought of as an input to the linear operator
        that gives an output with particularly large norm.
    w : ndarray, optional
        The vector Av which has relatively large 1-norm.
        It can be thought of as an output of the linear operator
        that is relatively large in norm compared to the input.

    """
    #XXX Eventually turn this into an API function in the _onenormest module,
    #XXX and remove its underscore,
    #XXX but wait until expm_action and expm_2009 go into scipy.
    return scipy.sparse.linalg.onenormest(MatrixPowerOperator(A, p))


def _onenormest_product(operator_seq,
        t=2, itmax=5, compute_v=False, compute_w=False):
    """
    Efficiently estimate the 1-norm of the matrix product of the args.

    Parameters
    ----------
    operator_seq : linear operator sequence
        Matrices whose 1-norm of product is to be computed.
    t : int, optional
        A positive parameter controlling the tradeoff between
        accuracy versus time and memory usage.
        Larger values take longer and use more memory
        but give more accurate output.
    itmax : int, optional
        Use at most this many iterations.
    compute_v : bool, optional
        Request a norm-maximizing linear operator input vector if True.
    compute_w : bool, optional
        Request a norm-maximizing linear operator output vector if True.

    Returns
    -------
    est : float
        An underestimate of the 1-norm of the sparse matrix.
    v : ndarray, optional
        The vector such that ||Av||_1 == est*||v||_1.
        It can be thought of as an input to the linear operator
        that gives an output with particularly large norm.
    w : ndarray, optional
        The vector Av which has relatively large 1-norm.
        It can be thought of as an output of the linear operator
        that is relatively large in norm compared to the input.

    """
    #XXX Eventually turn this into an API function in the _onenormest module,
    #XXX and remove its underscore,
    #XXX but wait until expm_2009 goes into scipy.
    return scipy.sparse.linalg.onenormest(ProductOperator(*operator_seq))


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

    Notes
    -----
    This is algorithm (6.1) which is a simplification of algorithm (5.1).

    References
    ----------
    .. [1] Awad H. Al-Mohy and Nicholas J. Higham (2009)
           "A New Scaling and Squaring Algorithm for the Matrix Exponential."
           SIAM Journal on Matrix Analysis and Applications.
           31 (3). pp. 970-989. ISSN 1095-7162

    """
    # Detect upper triangularity.
    structure = UPPER_TRIANGULAR if _is_upper_triangular(A) else None

    # Define the identity matrix depending on sparsity.
    ident = _ident_like(A)

    # Try Pade order 3.
    A2 = A.dot(A)
    d6 = _onenormest_matrix_power(A2, 3)**(1/6.)
    eta_1 = max(_onenormest_matrix_power(A2, 2)**(1/4.), d6)
    if eta_1 < 1.495585217958292e-002 and _ell(A, 3) == 0:
        U, V = _pade3(A, ident, A2)
        return _solve_P_Q(U, V, structure=structure)

    # Try Pade order 5.
    A4 = A2.dot(A2)
    d4 = _exact_1_norm(A4)**(1/4.)
    eta_2 = max(d4, d6)
    if eta_2 < 2.539398330063230e-001 and _ell(A, 5) == 0:
        U, V = _pade5(A, ident, A2, A4)
        return _solve_P_Q(U, V, structure=structure)

    # Try Pade orders 7 and 9.
    A6 = A2.dot(A4)
    d6 = _exact_1_norm(A6)**(1/6.)
    d8 = _onenormest_matrix_power(A4, 2)**(1/8.)
    eta_3 = max(d6, d8)
    if eta_3 < 9.504178996162932e-001 and _ell(A, 7) == 0:
        U, V = _pade7(A, ident, A2, A4, A6)
        return _solve_P_Q(U, V, structure=structure)
    if eta_3 < 2.097847961257068e+000 and _ell(A, 9) == 0:
        U, V = _pade9(A, ident, A2, A4, A6)
        return _solve_P_Q(U, V, structure=structure)

    # Use Pade order 13.
    d10 = _onenormest_product((A4, A6))**(1/10.)
    eta_4 = max(d8, d10)
    eta_5 = min(eta_3, eta_4)
    theta_13 = 4.25
    s = max(int(np.ceil(np.log2(eta_5 / theta_13))), 0)
    s = s + _ell(2**-s * A, 13)
    B = A * 2**-s
    B2 = A2 * 2**(-2*s)
    B4 = A4 * 2**(-4*s)
    B6 = A6 * 2**(-6*s)
    U, V = _pade13(B, ident, B2, B4, B6)
    X = _solve_P_Q(U, V, structure=structure)
    if structure == UPPER_TRIANGULAR:
        # Invoke Code Fragment 2.1.
        X = _fragment_2_1(X, A, s)
    else:
        # X = r_13(A)^(2^s) by repeated squaring.
        for i in range(s):
            X = X.dot(X)
    return X


def _solve_P_Q(U, V, structure=None):
    """
    A helper function for expm_2009.

    Parameters
    ----------
    U : ndarray
        Pade numerator.
    V : ndarray
        Pade denominator.
    structure : str, optional
        A string describing the structure of both matrices `U` and `V`.
        Only `upper_triangular` is currently supported.

    Notes
    -----
    The `structure` argument is inspired by similar args
    for theano and cvxopt functions.

    """
    P = U + V
    Q = -U + V
    if isspmatrix(U):
        return spsolve(Q, P)
    elif structure is None:
        return solve(Q, P)
    elif structure == UPPER_TRIANGULAR:
        return solve_triangular(Q, P)
    else:
        raise ValueError('unsupported matrix structure: ' + str(structure))


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
        return np.sinh(x) / x


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
    return t_12 * np.exp(a) * _sinch(b)


def _fragment_2_1(X, T, s):
    """
    A helper function for expm_2009.

    Notes
    -----
    The argument X is modified in-place, but this modification is not the same
    as the returned value of the function.
    This function also takes pains to do things in ways that are compatible
    with sparse matrices, for example by avoiding fancy indexing
    and by using methods of the matrices whenever possible instead of
    using functions of the numpy or scipy libraries themselves.

    """
    # Form X = r_m(2^-s T)
    # Replace diag(X) by exp(2^-s diag(T)).
    n = X.shape[0]
    diag_T = T.diagonal().copy()

    # Replace diag(X) by exp(2^-s diag(T)).
    scale = 2 ** -s
    exp_diag = np.exp(scale * diag_T)
    for k in range(n):
        X[k, k] = exp_diag[k]

    for i in range(s-1, -1, -1):
        X = X.dot(X)

        # Replace diag(X) by exp(2^-i diag(T)).
        scale = 2 ** -i
        exp_diag = np.exp(scale * diag_T)
        for k in range(n):
            X[k, k] = exp_diag[k]

        # Replace (first) superdiagonal of X by explicit formula
        # for superdiagonal of exp(2^-i T) from Eq (10.42) of
        # the author's 2008 textbook
        # Functions of Matrices: Theory and Computation.
        for k in range(n-1):
            lam_1 = scale * diag_T[k]
            lam_2 = scale * diag_T[k+1]
            t_12 = scale * T[k, k+1]
            value = _eq_10_42(lam_1, lam_2, t_12)
            X[k, k+1] = value

    # Return the updated X matrix.
    return X


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
    p = 2*m + 1

    # The c_i are explained in (2.2) and (2.6) of the 2005 expm paper.
    # They are coefficients of terms of a generating function series expansion.
    abs_c_recip = scipy.misc.comb(2*p, p, exact=True) * math.factorial(2*p + 1)

    # This is explained after Eq. (1.2) of the 2009 expm paper.
    # It is the "unit roundoff" of IEEE double precision arithmetic.
    u = 2**-53

    # Estimate the 1-norm of the matrix power.
    est = _onenormest_matrix_power(abs(A), p)

    # Treat zero norm as a special case.
    if not est:
        return 0

    alpha = est / (_exact_1_norm(A) * abs_c_recip)
    log2_alpha_div_u = np.log2(alpha/u)
    value = int(np.ceil(log2_alpha_div_u / (2 * m)))
    return max(value, 0)


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
