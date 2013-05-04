"""Compute the action of the matrix exponential.
"""

import math
import functools

import numpy as np
import numpy.linalg

import scipy.misc
import scipy.linalg
import scipy.optimize
import scipy.sparse.linalg

__all__ = ['expm_action']


# This table helps to compute bounds.
# They seem to have been difficult to calculate, involving symbolic
# manipulation of equations, followed by numerical root finding.
_theta = {
        # The first 30 values are from table A.3 of Computing Matrix Functions.
        1: 2.29e-16,
        2: 2.58e-8,
        3: 1.39e-8,
        4: 3.40e-4,
        5: 2.40e-3,
        6: 9.07e-3,
        7: 2.38e-2,
        8: 5.00e-2,
        9: 8.96e-2,
        10: 1.44e-1,
        # 11
        11: 2.14e-1,
        12: 3.00e-1,
        13: 4.00e-1,
        14: 5.14e-1,
        15: 6.41e-1,
        16: 7.81e-1,
        17: 9.31e-1,
        18: 1.09,
        19: 1.26,
        20: 1.44,
        # 21
        21: 1.62,
        22: 1.82,
        23: 2.01,
        24: 2.22,
        25: 2.43,
        26: 2.64,
        27: 2.86,
        28: 3.08,
        29: 3.31,
        30: 3.54,
        # The rest are from table 3.1 of
        # Computing the Action of the Matrix Exponential.
        35 : 4.7,
        40 : 6.0,
        45 : 7.2,
        50 : 8.5,
        55 : 9.9,
        }


def expm_action(A, B, start=None, stop=None, num=None, endpoint=None):
    """
    Compute the action of the matrix exponential of A on B.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose exponential is of interest.
    B : ndarray
        The matrix to be multiplied by the matrix exponential of A.
    start : scalar, optional
        The starting time point of the sequence.
    stop : scalar, optional
        The end time point of the sequence, unless `endpoint` is set to False.
        In that case, the sequence consists of all but the last of ``num + 1``
        evenly spaced time points, so that `stop` is excluded.
        Note that the step size changes when `endpoint` is False.
    num : int, optional
        Number of time points to use.
    endpoint : bool, optional
        If True, `stop` is the last time point.  Otherwise, it is not included.

    Returns
    -------
    expm_A_B : ndarray
         The result of the action :math:`e^A B`.

    Notes
    -----
    The optional arguments defining the sequence of evenly spaced time points
    are compatible with the arguments of `numpy.linspace`.

    References
    ----------
    .. [1] Awad H. Al-Mohy and Nicholas J. Higham (2011)
           "Computing the Action of the Matrix Exponential,
           with an Application to Exponential Integrators."
           SIAM Journal on Scientific Computing,
           33 (2). pp. 488-511. ISSN 1064-8275
           http://eprints.ma.man.ac.uk/1591/

    .. [2] Nicholas J. Higham and Awad H. Al-Mohy (2010)
           "Computing Matrix Functions."
           Acta Numerica,
           19. 159-208. ISSN 0962-4929
           http://eprints.ma.man.ac.uk/1451/

    """
    pass


def _expm_action_simple(A, B, t=1.0, balance=False):
    """
    Compute the action of the matrix exponential at a single time point.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose exponential is of interest.
    B : ndarray
        The matrix to be multiplied by the matrix exponential of A.
    t : float
        A time point.
    balance : bool
        Indicates whether or not to apply balancing.

    Returns
    -------
    F : ndarray
        :math:`e^{t A} B`

    Notes
    -----
    This is algorithm (3.2) in Al-Mohy and Higham (2011).

    """
    if balance:
        raise NotImplementedError
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected A to be like a square matrix')
    if A.shape[1] != B.shape[0]:
        raise ValueError('the matrices A and B have incompatible shapes')
    n = A.shape[0]
    n0 = B.shape[1]
    u_d = 2**-53
    tol = u_d
    mu = np.trace(A) / float(n)
    A = A - mu * np.identity(n)
    A_1_norm = np.linalg.norm(A, 1)
    if t*A_1_norm == 0:
        m_star, s = 0, 1
    else:
        m_star, s = _fragment_3_1(t*A, t*A_1_norm, n0, tol)
    F = B
    eta = math.exp(t*mu / float(s))
    for i in range(s):
        c1 = np.linalg.norm(B, np.inf)
        for j in range(m_star):
            B = t * np.dot(A, B) / float((s+1)*(j+1))
            c2 = np.linalg.norm(B, np.inf)
            F = F + B
            if c1 + c2 <= tol * np.linalg.norm(F, np.inf):
                break
            c1 = c2
        F = eta * F
        B = F
    return F


def _fragment_3_1(A, A_1_norm, n0, tol, m_max=55, p_max=8, ell=2):
    """
    A helper function for the _expm_action_* functions.

    Parameters
    ----------
    A : transposable linear operator
        Linear operator for which to compute the matrix exponential.
    A_1_norm : float
        The precomputed 1-norm of A.
    n0 : int
        Number of columns in the _expm_action_* B matrix.
    tol : float
        Expected to be
        :math:`2^{-24}` for single precision or
        :math:`2^{-53}` for double precision.
    m_max : int
        A value related to a bound.
    p_max : int
        A value related to a bound.
    ell : int
        The number of columns used in the 1-norm approximation.
        This is usually taken to be small, maybe between 1 and 5.

    Returns
    -------
    m_star : int
        Related to bounds for error control.
    s : int
        Amount of scaling.

    Notes
    -----
    This is code fragment (3.1) in Al-Mohy and Higham (2011).
    The discussion of default values for m_max, p_max, and ell
    is given between the definitions of equation (3.11)
    and the definition of equation (3.12).

    """
    if ell < 1:
        raise ValueError('expected ell to be a positive integer')
    if _condition_3_13(A_1_norm, n0, m_max, p_max, ell):
        argmin()
        m_star = np.argmin(what)
        s = foo
    m_star, s = foo


def _condition_3_13(A_1_norm, n0, m_max, p_max, ell):
    """
    A helper function for the _expm_action_* functions.

    Parameters
    ----------
    A_1_norm : float
        The precomputed 1-norm of A.
    n0 : int
        Number of columns in the _expm_action_* B matrix.
    m_max : int
        A value related to a bound.
    p_max : int
        A value related to a bound.
    ell : int
        The number of columns used in the 1-norm approximation.
        This is usually taken to be small, maybe between 1 and 5.

    Returns
    -------
    value : bool
        Indicates whether or not the condition has been met.
        
    Notes
    -----
    This is condition (3.13) in Al-Mohy and Higham (2011).

    """

    # This is the rhs of equation (3.12).
    a = 2 * ell * p_max * (p_max + 3)

    # Evaluate the condition (3.13).
    b = theta[m_max] / float(n0 * m_max)
    return A_1_norm <= a * b

