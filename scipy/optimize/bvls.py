from . import _bvls
from numpy import asarray_chkfinite, zeros, array, isfinite, inf

__all__ = ['bounded_lstsq']


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
    bounds : list
        A list of tuples specifying the lower and upper bound for each
        independent variables [(xl0, xu0), (xl1, xu1), ...] Infinite values
        or None will be interpreted as large floating values, negated as
        appropriate.

    Returns
    -------
    x : ndarray
        Solution vector
    rnorm : float
        The residual, ``|| Ax-b ||_2``.

    See also
    --------
    scipy.optimize.nnls

    Notes
    -----
    The Fortran code for this was originally written by Charles Lawson and
    Richard Hanson, who agreed to release the code under the BSD license
    for inclusion in scipy.
    """
    A, b = map(asarray_chkfinite, (A, b))
    A = A.astype(dtype=float, order='F')
    b = b.astype(dtype=float, order='F')

    if A.ndim != 2:
        raise ValueError("A must be 2-dimensional")
    m, n = A.shape

    if b.ndim != 1:
        raise ValueError("Expected 1-dimensional b")
    if A.shape[0] != b.shape[0]:
        raise ValueError("A and b are not conformable")

    if bounds is None or len(bounds) == 0:
        bnds = array([[-1.0E12]*n, [1.0E12]*n])
    else:
        bnds = array(bounds, float).T
        infbnd = ~isfinite(bnds)
        bnds[0, infbnd[0]] = -inf
        bnds[1, infbnd[1]] = inf
        if bnds.shape[1] != n:
            raise ValueError("The length of bounds is not compatible with "
                             "Ax=b. Got %d. Expected %d" (len(bnds), n))

    w = zeros((n,), dtype=float, order='F')
    index = zeros((n,), dtype=int, order='F')
    x = zeros((n,), dtype=float, order='F')

    rnorm, nsetp, ierr = _bvls.bvls(A, b, bnds, x, w, index)
    if ierr == 1:
        raise ValueError("M <= 0 or N <= 0")
    elif ierr == 2:
        raise ValueError("Size or shape violation.")
    elif ierr == 3:
        raise ValueError("Input bounds are inconsistent")
    elif ierr == 4:
        raise ValueError("Exceeded maximum number of iterations.")
    return x, rnorm, nsetp
