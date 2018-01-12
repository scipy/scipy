from __future__ import division, print_function, absolute_import
import numpy
from numpy import empty_like, sqrt, spacing, tril, triu
from scipy.linalg import LinAlgError

__all__ = ['cholesky_update']


def cholesky_update(R, z, mod_type="+", lower=False, overwrite_R=False,
               overwrite_z=False, check_finite=True):
    """
    Performs a rank-1 update or downdate of Cholesky factors.

    Cholesky rank-1 update and downdate algorithm taken from [1]_ such that
    an update results in ``D'D = R'R + zz'`` and a downdate results in ``D'D =
    R'R - zz'``

    Parameters
    ----------
    R : (N, N) ndarray
        The 2D square array from which the triangular part will be used to
        read the Cholesky factor.The remaining parts are ignored.
    z : (N,) ndarray
        A 1D modifier (update or downdate) array
    mod_type: string, optional
        The type of modification desired. A "+" indicates an update while a
        "-" indicates a downdate. Default is "+" (update).
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
    .. [1] DOI:10.1007/BF01933218
    """

    # Check for invalid mod_type argument
    if mod_type not in ['+', '-']:
        raise ValueError("'mod_type' argument should be one of ['+','-']")
    # Copy or overwrite?
    D = R if overwrite_R else R.copy()
    z = z if overwrite_z else z.copy()
    # Check finite, else force array
    if check_finite:
        R = numpy.asarray_chkfinite(R)
        z = numpy.asarray_chkfinite(z)
    else:
        R = numpy.asarray(R)
        z = numpy.asarray(z)
    # Check for appropriate array dimensions
    if len(R.shape) != 2:
        raise ValueError("Expected 2D array to be updated.")
    if len(z.squeeze().shape) != 1:
        raise ValueError("Expected 1D update vector")
    # Ensure square array
    m, n = R.shape
    if m != n:
        raise ValueError("A square array is expected.")
    # Ensure appropriately triangular array
    is_proper_tri = {True: numpy.allclose(tril(R), R),
                     False: numpy.allclose(triu(R), R)}[lower]
    if not is_proper_tri:
        raise ValueError(
            "Expected {} triangular array.".format("lower" if lower else
                                                   "upper"))
    # Cleanup
    del m, is_proper_tri
    eps = n * spacing(1.)  # For complex this needs modification
    alpha, beta = empty_like(z), empty_like(z)
    alpha[-1], beta[-1] = 1., 1.

    for r in range(n):
        a = z[r] / R[r, r]
        alpha[r] = alpha[r - 1] + a ** 2 if mod_type == "+" else \
            alpha[r - 1] - a ** 2
        # Numerically zero or negative
        if alpha[r] < eps:
            # Made up err msg.
            raise LinAlgError('The Cholesky factor becomes nonpositive'
                              'with this downdate at the step {}'.format(r))
        beta[r] = sqrt(alpha[r])
        z[r + 1:] -= a * R[r, r + 1:]
        D[r, r:] *= beta[r] / beta[r - 1]
        D[r, r + 1:] += a / (beta[r] * beta[r - 1]) * z[r + 1:] if \
            mod_type == "+" else -a / (beta[r] * beta[r - 1]) * z[r + 1:]

    return D
