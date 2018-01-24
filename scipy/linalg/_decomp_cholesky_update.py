from __future__ import division, print_function, absolute_import
from numpy import empty_like, sqrt, spacing
from .misc import LinAlgError, _datacopied
from .decomp import _asarray_validated

__all__ = ['cholesky_update']

def cholesky_update(R, z, downdate=False, lower=False, overwrite_R=False,
               overwrite_z=False, eps=None, check_finite=True):
    """
    Performs a rank-1 update or downdate of Cholesky factors.

    Cholesky rank-1 update and downdate algorithm taken from [1]_ such that
    an update results in ``D'D = R'R + zz'`` and a downdate results in ``D'D =
    R'R - zz'``

    Parameters
    ----------
    R : (N, N) array_like
        The 2D input data from which the triangular part will be used to
        read the Cholesky factor.The remaining parts are ignored.
    z : (N,) array_like
        A 1D update/downdate vector
    downdate: bool, optional
        The type of rank-1 modification desired. False indicates an update
        while True indicates a downdate. Default is False.
    lower: bool, optional
        Whether input array is upper or lower triangular Cholesky
        factorization. Default is upper-triangular.
    overwrite_R: bool, optional
        If set to True the entries of the Cholesky factor array will be
        modified during the computations instead of creating a new array.
        Default is False.
    overwrite_z : bool, optional
        If set to True the entries of the update array will be modified during
        the computations instead of creating a new array. Default is False.
    eps : float, optional
        This determines the tolerance below which we consider ``alpha`` values
        to be effectively zero. Nominally this value should be zero, but
        numerical issues cause a tolerance about zero to be necessary.
        If the default value is unchanged then it is set to
        ``R.shape[1] * numpy.spacing(1.).`` during runtime. Default is None.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
        Default is True.

    Returns
    -------
    D : (N, N) ndarray
        The resulting modified Cholesky factor.

    References
    ----------
    .. [1] : C.T. Pan, "A modification to the LINPACK downdating algorithm",
       BIT Numerical Mathematics, Vol.30(4), 1990, DOI:10.1007/BF02165404

    """
    # Input validation
    R1 = _asarray_validated(R, check_finite=check_finite)
    z1 = _asarray_validated(z, check_finite=check_finite)
    m, n = R1.shape

    overwrite_R = overwrite_R or _datacopied(R1, R)
    overwrite_z = overwrite_z or _datacopied(z1, z)

    # Dimension check
    if len(R1.shape) != 2:
        raise ValueError("Expected 2D array to be updated.")
    if len(z1.squeeze().shape) != 1:
        raise ValueError("Expected 1D update vector")

    # Square matrix check
    if m != n:
        raise ValueError('Input needs to be a square matrix.')

    # Compatible dimensions check
    if m != len(z1):
            raise ValueError('Input z has to have same length as the number '
                             'of rows as input R')
    # Copy or not
    R = R1 if overwrite_R else R1.copy()
    z = z1 if overwrite_z else z1.copy()

    # Initializations
    if eps is None:
        eps = n * spacing(1.)  # For complex this needs modification
    alpha, beta = empty_like(z), empty_like(z)
    alpha[-1], beta[-1] = 1., 1.
    sign = -1 if downdate else 1
    R = R.T if lower else R

    for r in range(n):
        a = z[r] / R[r, r]
        alpha[r] = alpha[r - 1] + sign * a ** 2
        # Numerically zero or negative
        if alpha[r] < eps:
            # Error msg.
            raise LinAlgError('The Cholesky factor becomes nonpositive'
                              'with this downdate at the step {}'.format(r))
        beta[r] = sqrt(alpha[r])
        z[r + 1:] -= a * R[r, r + 1:]
        R[r, r:] *= beta[r] / beta[r - 1]
        R[r, r + 1:] += sign * a / (beta[r] * beta[r - 1]) * z[r + 1:]

    return R if not lower else R.T
