"""
Matrix functions that use Pade approximation with inverse scaling and squaring.

"""
from __future__ import division, print_function, absolute_import

import numpy as np

from scipy.linalg._matfuncs_sqrtm import _sqrtm_triu
from scipy.linalg.decomp_schur import schur, rsf2csf
from scipy.linalg.special_matrices import all_mat
from scipy.linalg import solve
from scipy.sparse.linalg.interface import LinearOperator
from scipy.sparse.linalg import onenormest


__all__ = ['logm_new', 'fractional_matrix_power']


#TODO renovate or move this class when scipy operators are more mature
class _MatrixM1PowerOperator(LinearOperator):
    """
    A representation of the linear operator (A - I)^p.
    """

    def __init__(self, A, p):
        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ValueError('expected A to be like a square matrix')
        if p < 0 or p != int(p):
            raise ValueError('expected p to be a non-negative integer')
        self._A = A
        self._p = p
        self.ndim = A.ndim
        self.shape = A.shape

    def matvec(self, x):
        for i in range(self._p):
            x = self._A.dot(x) - x
        return x

    def rmatvec(self, x):
        for i in range(self._p):
            x = x.dot(self._A) - x
        return x

    def matmat(self, X):
        for i in range(self._p):
            X = self._A.dot(X) - X
        return X

    @property
    def T(self):
        return _MatrixM1PowerOperator(self._A.T, self._p)


#TODO renovate or move this function when scipy operators are more mature
def _onenormest_m1_power(A, p,
        t=2, itmax=5, compute_v=False, compute_w=False):
    """
    Efficiently estimate the 1-norm of (A - I)^p.

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
    return onenormest(_MatrixM1PowerOperator(A, p))


def _almost_diagonal(M, abstol=1e-13):
    M_abs_off_diagonal = np.absolute(M - np.diag(np.diag(M)))
    return np.max(M_abs_off_diagonal) < abstol


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


def _briggs_helper_function(a, k):
    """
    Computes r = a^(1 / (2^k)) - 1.

    This is algorithm (2) of [1].
    The purpose is to avoid a danger of subtractive cancellation.
    For more computational efficiency it should probably be cythonized.

    Parameters
    ----------
    a : complex
        A complex number not belonging to the closed negative real axis.
    k : integer
        A nonnegative integer.

    Returns
    -------
    r : complex
        The value r = a^(1 / (2^k)) - 1 computed with less cancellation.

    Notes
    -----
    The algorithm as written in the publication does not handle k=0 or k=1
    correctly, so these are special-cased in this implementation.

    References
    ----------
    .. [1] Awad H. Al-Mohy (2012)
           "A more accurate Briggs method for the logarithm",
           Numerical Algorithms, 59 : 393--402.

    """
    if (a.real == a) and (a <= 0):
        raise ValueError("expected a complex number 'a' not belonging "
                "to the closed negative real axis")
    if k < 0 or int(k) != k:
        raise ValueError('expected a nonnegative integer k')
    if k == 0:
        return a - 1
    elif k == 1:
        return np.sqrt(a) - 1
    else:
        k_hat = k
        if np.angle(a) >= np.pi / 2:
            a = np.sqrt(a)
            k_hat = k - 1
        z0 = a - 1
        a = np.sqrt(a)
        r = 1 + a
        for j in range(1, k_hat):
            a = np.sqrt(a)
            r = r * (1 + a)
        r = z0 / r
        return r


def _fractional_power_superdiag_entry(l1, l2, t12, p):
    """
    Compute a superdiagonal entry of a fractional matrix power.

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
        A superdiagonal entry of the fractional matrix power.

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
        arctanh_z = np.arctanh(z)
        tmp_a = t12 * np.exp((p/2)*(log_l2 + log_l1))
        tmp_b = p * (arctanh_z + np.pi * 1j * _unwindk(log_l2 - log_l1))
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
        return t12 * (2*np.arctanh(z) + 2*np.pi*1j*(ua + ub)) / (l2 - l1)


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
        matrix logarithm than for fractional matrix power.

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
    The instances of 'goto line 38' in algorithm (3.1) of [1]
    probably mean 'goto line 36' and have been intepreted accordingly.

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
    d2 = _onenormest_m1_power(T, 2) ** (1/2)
    d3 = _onenormest_m1_power(T, 3) ** (1/3)
    a2 = max(d2, d3)
    m = None
    for i in (1, 2):
        if a2 <= theta[i]:
            m = i
            break
    while m is None:
        if s > s0:
            d3 = _onenormest_m1_power(T, 3) ** (1/3)
        d4 = _onenormest_m1_power(T, 4) ** (1/4)
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
        d5 = _onenormest_m1_power(T, 5) ** (1/5)
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


def _fractional_power_pade_constant(i, t):
    # A helper function for matrix fractional power.
    if i < 1:
        raise ValueError('expected a positive integer i')
    if not (-1 < t < 1):
        raise ValueError('expected -1 < t < 1')
    if i == 1:
        return -t
    elif i % 2 == 0:
        j = i // 2
        return (-j + t) / (2 * (2*j - 1) )
    elif i % 2 == 1:
        j = (i - 1) // 2
        return (-j - t) / (2 * (2*j + 1) )
    else:
        raise Exception('internal error')


def _fractional_power_pade(R, t, m):
    """
    Evaluate the Pade approximation of a fractional matrix power.

    Evaluate the degree-m Pade approximation of R
    to the fractional matrix power t using the continued fraction
    in bottom-up fashion using algorithm (4.1) in [1].

    Parameters
    ----------
    R : (N, N) array_like
        Matrix whose fractional power to evaluate.
    t : float
        Fractional power between -1 and 1 exclusive.
    m : positive integer
        Degree of Pade approximation.

    Returns
    -------
    U : (N, N) array_like
        The degree-m Pade approximation of R to the fractional power t.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    """
    if m < 1:
        raise ValueError('expected a positive integer m')
    if not (-1 < t < 1):
        raise ValueError('expected -1 < t < 1')
    R = np.asarray(R)
    if len(R.shape) != 2:
        raise ValueError("Non-matrix input to matrix function.")
    if R.shape[0] != R.shape[1]:
        raise ValueError("Non-square matrix input to matrix function.")
    n, n = R.shape
    ident = np.identity(n)
    Y = R * _fractional_power_pade_constant(2*m, t)
    for j in range(2*m - 1, 0, -1):
        rhs = R * _fractional_power_pade_constant(j, t)
        Y = solve(ident + Y, rhs)
    return ident + Y


def _remainder_matrix_power(A, t):
    """
    Compute the fractional power of a matrix, for fractions -1 < t < 1.

    This uses algorithm (3.1) of [1].
    The Pade approximation itself uses algorithm (4.1) of [2].

    Parameters
    ----------
    A : (N, N) array_like
        Matrix whose fractional power to evaluate.
    t : float
        Fractional power between -1 and 1 exclusive.

    Returns
    -------
    X : (N, N) array_like
        The fractional power of the matrix.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing Lin (2013)
           "An Improved Schur-Pade Algorithm for Fractional Powers
           of a Matrix and their Frechet Derivatives."

    .. [2] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    """
    A = np.asarray(A)
    if len(A.shape) != 2:
        raise ValueError("Non-matrix input to matrix function.")
    if A.shape[0] != A.shape[1]:
        raise ValueError("Non-square matrix input to matrix function.")
    n, n = A.shape
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
    T0 = T
    T0_diag = np.diag(T0)
    if _almost_diagonal(T0):
        U = np.diag(T0_diag ** t)
    else:
        T, s, m = _inverse_squaring_helper(T0, m_to_theta)
        R = np.identity(n) - T

        # Replace the diagonal and first superdiagonal of I - T0^(1/(2^s))
        # using formulas that have less subtractive cancellation.
        for j in range(n):
            a = T0[j, j]
            r = _briggs_helper_function(a, s)
            R[j, j] = -r
        p = np.exp2(-s)
        for j in range(n-1):
            l1 = T0[j, j]
            l2 = T0[j+1, j+1]
            t12 = T0[j, j+1]
            f12 = _fractional_power_superdiag_entry(l1, l2, t12, p)
            R[j, j+1] = -f12

        # Evaluate the Pade approximation.
        U = _fractional_power_pade(R, t, m)

        # Undo the inverse scaling and squaring.
        for i in range(s, -1, -1):
            if i < s:
                U = U.dot(U)
            else:
                p = t * np.exp2(-i)
                U[np.diag_indices(n)] = T0_diag ** p
                for j in range(n-1):
                    l1 = T0[j, j]
                    l2 = T0[j+1, j+1]
                    t12 = T0[j, j+1]
                    f12 = _fractional_power_superdiag_entry(l1, l2, t12, p)
                    U[j, j+1] = f12

    U, Z = all_mat(U, Z)
    X = (Z * U * Z.H)
    return X.A


def fractional_matrix_power(A, p):
    """
    Compute the fractional power of a matrix.

    Proceeds according to the discussion in section (6) of [1].

    Parameters
    ----------
    A : (N, N) array_like
        Matrix whose fractional power to evaluate.
    p : float
        Fractional power.

    Returns
    -------
    X : (N, N) array_like
        The fractional power of the matrix.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    """
    A = np.asarray(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected a square matrix')
    if p == int(p):
        return np.linalg.matrix_power(A, int(p))
    p1 = p - np.floor(p)
    p2 = p - np.ceil(p)
    k2 = np.linalg.cond(A)
    if p1 * k2 ** (1 - p1) <= -p2 * k2:
        a = int(np.floor(p))
        b = p1
    else:
        a = int(np.ceil(p))
        b = p2
    Q = np.linalg.matrix_power(A, a)
    R = _remainder_matrix_power(A, b)
    return Q.dot(R)


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
    A = np.asarray(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected a square matrix')
    n, n = A.shape
    T, Z = schur(A)
    T, Z = rsf2csf(T,Z)
    T0 = T

    # Define bounds given in Table (2.1).
    theta = (None,
            1.59e-5, 2.31e-3, 1.94e-2, 6.21e-2,
            1.28e-1, 2.06e-1, 2.88e-1, 3.67e-1,
            4.39e-1, 5.03e-1, 5.60e-1, 6.09e-1,
            6.52e-1, 6.89e-1, 7.21e-1, 7.49e-1)

    T, s, m = _inverse_squaring_helper(T0, theta)
    R = T - np.identity(n)

    # Replace the diagonal and first superdiagonal of T0^(1/(2^s)) - I
    # using formulas that have less subtractive cancellation.
    for j in range(n):
        a = T0[j, j]
        r = _briggs_helper_function(a, s)
        R[j, j] = r
    p = np.exp2(-s)
    for j in range(n-1):
        l1 = T0[j, j]
        l2 = T0[j+1, j+1]
        t12 = T0[j, j+1]
        f12 = _fractional_power_superdiag_entry(l1, l2, t12, p)
        R[j, j+1] = f12

    # Evaluate U = 2**s r_m(T - I) using the partial fraction expansion (1.1).
    # This requires the nodes and weights
    # corresponding to degree-m Gauss-Legendre quadrature.
    # These quadrature arrays need to be transformed from the [-1, 1] interval
    # to the [0, 1] interval.
    nodes, weights = np.polynomial.legendre.leggauss(m)
    if nodes.shape != (m,) or weights.shape != (m,):
        raise Exception('internal error')
    nodes = 0.5 + 0.5 * nodes
    weights = 0.5 * weights
    ident = np.identity(n)
    U = np.zeros_like(R)
    for alpha, beta in zip(weights, nodes):
        U += solve(ident + beta*R, alpha*R)
    U *= np.exp2(s)

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

