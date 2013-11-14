from __future__ import division, print_function, absolute_import

#
# Author: Alex LB Leach, November 2012
#

import numpy as np
from scipy.linalg._expokit_f import dgpadm, zgpadm, dspadm, zhpadm, \
                                    dgexpv, zgexpv
from scipy.sparse.base import isspmatrix
from scipy.sparse.linalg import aslinearoperator, LinearOperator


def expm(A, t=None):
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

    is_sparse_op = (isspmatrix(A) or isinstance(A, LinearOperator))

    if not is_sparse_op:
        A = np.asarray(A)

    m   = A.shape[0]
    dtype = A.dtype

    # output integers needed to get result:
    iexph = np.array([0])
    ns    = np.array([0])
    iflag = np.array([0])

    if is_sparse_op:
        # Krylov estimation
        Av = aslinearoperator(A)
        itrace = np.array([0])
        tol = np.array([1e-7])
        v = np.ones(m, dtype=dtype)
        anorm = np.linalg.norm(A.todense(), np.inf)

        if np.issubdtype(dtype, np.complexfloating):
            w, wsp = zgexpv(v * 1 + 1j, tol, anorm, Av.matvec, itrace, iflag, t=t) 
        else:
            # See src/expokit.f for documentation
            w, wsp = dgexpv(v, tol, anorm, Av.matvec, itrace, iflag, t=t)

        return np.reshape(w, (m,m), order='F')
    else:
        iexph = np.array([0])
        ns = np.array([0])
        if np.issubdtype(dtype, np.complexfloating):
            wsp = zgpadm(A, iexph, ns, iflag, t=t)
        else:
            wsp = dgpadm(A, iexph, ns, iflag, t)

    if iflag[0] != 0:
        raise IOError("Expokit error ({0})".format(iflag[0]))

    return np.reshape(wsp[iexph[0]-1 : iexph[0] + m * m - 1], (m, m), order='F')
