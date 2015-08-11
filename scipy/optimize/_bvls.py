"""
The Fortran code for this was originally written by Charles Lawson and
Richard Hanson, who agreed to release the code under the BSD license
for inclusion in scipy.
"""

from . import _bvlslib
import numpy as np

__all__ = ['bounded_lstsq']

_error_msg_dict = {
                   1: "M <= 0 or N <= 0",
                   2: "Size or shape violation.",
                   3: "Input bounds are inconsistent",
                   4: "Exceeded maximum number of iterations."
                   }


def bounded_lstsq(A, b, bounds=()):
    """
    Bounded variable (linear) least-squares.

    Solve ``argmin_x || Ax - b ||_2`` for possibly bounded ``x``.

    Parameters
    ----------
    A : ndarray
        Matrix ``A`` as shown above.
    b : ndarray
        Right-hand side vector.
    bounds : array_like
        An (N,2) shape array or a list of tuples specifying the lower and upper
        bound for each independent variables ``[(xl0, xu0), (xl1, xu1), ...]``.
        None or -inf, inf, can be used to indicate no bounds.

    Returns
    -------
    x : ndarray
        Solution vector
    rnorm : float
        The residual, ``|| Ax-b ||_2``.
    nsetp : int
        The number of items in x that are not at their constraint value

    See also
    --------
    scipy.optimize.nnls

    Examples
    --------
    >>> from scipy import optimize

    >>> x = np.array([-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
    ...                 -2.8394297E-01,1.3416969E+00,1.3757038E+00,
    ...                 -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
    ...                 8.2065094E-01])[:,None]
    >>> x = np.column_stack((np.ones_like(x), x))
    >>> np.random.seed(12)
    >>> y = np.dot(x, [1.8, 5.4]) + np.random.random(len(x))
    >>> bounds = None

    >>> optimize.bounded_lstsq(x, y, bounds=bounds)
    (array([ 2.26327373,  5.42102807]), 1.160159078139885, 2)

    >>> bounds = [(None, None), (None, 4)]
    >>> optimize.bounded_lstsq(x, y, bounds=bounds)
    (array([ 2.39941956,  4.        ]), 5.379227084630758, 1)

    References
    ----------

    The reference for the original Fortran code ::

        Charles Lawson, Richard Hanson,
        Solving Least Squares Problems,
        SIAM, 1995,
        ISBN: 0898713560,
        LC: QA275.L38.
    """
    A = np.array(A, dtype=float, copy=True, order='F')
    b = np.array(b, dtype=float, copy=True, order='F')
    if not np.isfinite(A).all() or not np.isfinite(b).all():
        raise ValueError("A and b may not contain NaNs or infs")

    if A.ndim != 2:
        raise ValueError("A must be 2-dimensional")
    m, n = A.shape

    if b.ndim != 1:
        raise ValueError("Expected 1-dimensional b")
    if A.shape[0] != b.shape[0]:
        raise ValueError("A and b are not conformable")

    if bounds is None or len(bounds) == 0:
        bnds = np.array([[-np.inf]*n, [np.inf]*n])
    else:
        bnds = np.asarray(bounds, dtype=float).T
        infbnd = ~np.isfinite(bnds)
        bnds[0, infbnd[0]] = -np.inf
        bnds[1, infbnd[1]] = np.inf
        if bnds.shape[1] != n:
            raise ValueError("The length of bounds is not compatible with "
                             "Ax=b. Got %d. Expected %d" (len(bnds), n))

    w = np.zeros(n, dtype=float, order='F')
    index = np.zeros(n, dtype=int, order='F')
    x = np.zeros(n, dtype=float, order='F')

    rnorm, nsetp, ierr = _bvlslib.bvls(A, b, bnds, x, w, index)

    if ierr > 0:
        raise ValueError(_error_msg_dict[ierr])

    return x, rnorm, nsetp
