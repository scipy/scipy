#
# Author: Alex LB Leach, November 2012
#

from scipy.linalg.expokit import dgpadm, zgpadm, dspadm, zspadm
from scipy.sparse.base import isspmatrix


def expm(A, t=None ):
    """Compute the matrix exponential - exp(t*A) - using Padé approximation.

    Parameters
    ----------
    A : array or sparse matrix, shape(M,M)
        2D Array or Matrix (sparse or dense) to be exponentiated
    t : array, shape(1)
        Exponent multiplier. Default=1.

    Returns
    -------
    expA : array, shape(M,M)
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
    #t = np.array([1.])[0]
    #ideg = 6
    #lwsp = 4 * m*m + ideg + 1
    # output integers needed to get result.
    iexph = np.array([0])
    ns    = np.array([0])
    iflag = np.array([0])
    # Workspace / Output array
    #wsp = np.zeros( lwsp, dtype=dtype )
    if isspmatrix(A):
        if dtype in ('float64', 'float32', float,):
            wsp = dspadm( A, iexph, ns, iflag, t=t )
        elif dtype in ('complex128', 'complex64', complex):
            wsp = zhpadm( A, iexph, ns, iflag, t=t )
    else:
        if dtype in ('float64', 'float32', float,):
            wsp = dgpadm( A, iexph, ns, iflag, t=t )
        elif dtype in ('complex128', 'complex64', complex):
            wsp = zgpadm( A, iexph, ns, iflag, t=t )
    if iflag[0] != 0:
        raise IOError("Expokit routine returned error code: {0}".format(iflag[0]))
    return np.reshape(wsp[iexph[0]-1 : iexph[0] + m*m - 1], (m,m), order='F')

