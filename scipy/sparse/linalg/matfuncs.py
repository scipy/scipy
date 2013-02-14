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

from numpy import asarray, dot, eye, ceil, log2
from numpy import matrix as mat
import numpy as np

from scipy.linalg.misc import norm
from scipy.linalg.basic import solve, inv

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
    Q = -U + V # q_m(A) : denominator

    if Aissparse:
        from scipy.sparse.linalg import spsolve
        R = spsolve(Q, P)
    else:
        R = solve(Q,P)

    # squaring step to undo scaling
    for i in range(n_squarings):
        R = R.dot(R)

    return R

# implementation of Pade approximations of various degree using the algorithm presented in [Higham 2005]
# These should apply to both dense and sparse matricies.
# ident is the identity matrix, which matches A in being sparse or dense.
def _pade3(A, ident):
    b = (120., 60., 12., 1.)
    A2 = A.dot(A)
    U = A.dot(b[3]*A2 + b[1]*ident)
    V = b[2]*A2 + b[0]*ident
    return U,V

def _pade5(A, ident):
    b = (30240., 15120., 3360., 420., 30., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    U = A.dot(b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade7(A, ident):
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    A6 = A4.dot(A2)
    U = A.dot(b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade9(A, ident):
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
                2162160., 110880., 3960., 90., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    A6 = A4.dot(A2)
    A8 = A6.dot(A2)
    U = A.dot(b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

def _pade13(A, ident):
    b = (64764752532480000., 32382376266240000., 7771770303897600.,
    1187353796428800., 129060195264000., 10559470521600., 670442572800.,
    33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.)
    A2 = A.dot(A)
    A4 = A2.dot(A2)
    A6 = A4.dot(A2)
    U = A.dot(A6.dot(b[13]*A6 + b[11]*A4 + b[9]*A2) + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = A6.dot(b[12]*A6 + b[10]*A4 + b[8]*A2) + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U,V

