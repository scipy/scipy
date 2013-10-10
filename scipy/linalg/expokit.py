#
# Author: Alex LB Leach, November 2012
#

import numpy as np
from scipy.linalg._expokit import dgpadm, zgpadm, dspadm, zhpadm, \
                                  dgexpv, zgexpv
from scipy.sparse.base import isspmatrix
from scipy.sparse.linalg import LinearOperator


def expm(A, t=None ):
    """Compute the matrix exponential - exp(t*A) - using Pade approximation.

    Parameters
    ----------
    A : array or sparse matrix, shape(M,M)
        2D Array or Matrix (sparse or dense) to be exponentiated
    t : array, shape(1)
        Exponent multiplier. Default=1.

    Returns
    -------
    exp(A*t) : array, shape(M,M)
        Matrix exponential of A

    References
    ----------
    Roger B. Sidje,
    "EXPOKIT: A Software Package for Computing Matrix Exponentials",
    ACM Trans. Math. Softw. 24(1), 130-156 (1998).

    http://www.maths.uq.edu.au/expokit/

    """
    # dgpadm does an internal check for square matrices,
    #  so assume ldh == m
    # constant input variables
    if hasattr( A, 'shape' ):
        m   = A.shape[0]
        #ldh = A.shape[1] 
        dtype = A.dtype
    else:
        m   = len(A)
        #ldh = len(A[0])
        dtype = type( A[0][0] )
    #if t is None:
    #    t = np.array([1.])[0]

    # output integers needed to get result:
    iexph = np.array([0])
    ns    = np.array([0])
    iflag = np.array([0])
    if isspmatrix(A):
        # Sparse matrix routines.
        matvec = lambda v : A.dot(v)
        # Please check usage of LinearOperator:-
        Av = LinearOperator( (m, m), matvec=matvec)
        itrace = np.array([0])
        tol = np.array([1e-7])
        v = np.ones(m, dtype=dtype)
        anorm = np.linalg.norm(A.todense(), np.inf) 
        if dtype in ('float64', 'float32', float,):
            # See src/expokit.f for documentation
            w, wsp = dgexpv(v, tol, anorm, Av.matvec, itrace, iflag, t=t)
        elif dtype in ('complex128', 'complex64', complex):
            w, wsp = zgexpv(v * 1 + 1j, tol, anorm, Av.matvec, itrace, iflag, t=t) 
        #print('\nA {0}: {1}'.format(A.shape, type(A)))
        #print('w {0}: {1}'.format(w.shape, type(w)))
        return np.reshape(w, (m,m), order='F')
        #return wsp[m*(m+1)+m+(m+2)**2:]
        #return np.reshape(wsp[m*(m+1)+m+(m+2)**2:], (m,m), order='F')
    else:
        iexph = np.array([0])
        ns = np.array([0])
        if dtype in ('float64', 'float32', float,):
            wsp = dgpadm(A, iexph, ns, iflag, t)
        elif dtype in ('complex128', 'complex64', complex):
            wsp = zgpadm(A, iexph, ns, iflag, t=t)
    if iflag[0] != 0:
        raise IOError("Expokit error ({0}) in routine {1}".format(iflag[0]))
    return np.reshape(wsp[iexph[0]-1 : iexph[0] + m * m - 1], (m,m), order='F')

