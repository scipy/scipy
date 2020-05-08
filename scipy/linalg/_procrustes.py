"""
Solve the orthogonal Procrustes problem.

"""
import numpy as np
from .decomp_svd import svd


__all__ = ['orthogonal_procrustes']


def orthogonal_procrustes(A, B, check_finite=True):
    """
    Compute the matrix solution of the orthogonal (or unitary) Procrustes problem.

    Given matrices A and B of equal shape, find an orthogonal (unitary in the case of complex input matrices) matrix R
    that most closely maps A to B using the algorithm given in [1]_.

    Parameters
    ----------
    A : (M, N) array_like
        Matrix to be mapped.
    B : (M, N) array_like
        Target matrix.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    R : (N, N) ndarray
        The matrix solution of the orthogonal Procrustes problem.
        Minimizes the Frobenius norm of ``(A @ R) - B``, subject to
        ``R.T @ R = I``.
    scale : float
        Sum of the singular values of ``A.T @ B``.

    Raises
    ------
    ValueError
        If the input array shapes don't match or if check_finite is True and
        the arrays contain Inf or NaN.

    Notes
    -----
    Note that unlike higher level Procrustes analyses of spatial data, this
    function only uses orthogonal transformations like rotations and
    reflections, and it does not use scaling or translation.

    .. versionadded:: 0.15.0

    References
    ----------
    .. [1] Peter H. Schonemann, "A generalized solution of the orthogonal
           Procrustes problem", Psychometrica -- Vol. 31, No. 1, March, 1996.

    Examples
    --------
    >>> from scipy.linalg import orthogonal_procrustes
    >>> A = np.array([[ 2,  0,  1], [-2,  0,  0]])

    Flip the order of columns and check for the anti-diagonal mapping
    
    >>> R, sca = orthogonal_procrustes(A, np.fliplr(A))
    >>> R
    array([[-5.34384992e-17,  0.00000000e+00,  1.00000000e+00],
           [ 0.00000000e+00,  1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  0.00000000e+00, -7.85941422e-17]])
    >>> sca
    9.0

    # Complex example for the unitary Procrustes problem
    >>> n = 4
    >>> A = np.random.rand(n,n)+1j*np.random.rand(n,n)
    >>> Q, _ = np.linalg.qr(np.random.rand(n,)+1j*np.random.rand(n,n)) # a random unitary matrix
    >>> B = A@Q+1e-4*(np.random.randn(n,n)+1j*np.random.randn(n,n))    # B is A@Q plus complex noise
    >>> R, _ = unitary_procrustes(A,B)                                 # find best estimate of Q
    >>> print('Should be small: %g'%np.linalg.norm(R-Q))
    Should be small: 0.00111944
    
    """
    
    if check_finite:
        A = np.asarray_chkfinite(A)
        B = np.asarray_chkfinite(B)
    else:
        A = np.asanyarray(A)
        B = np.asanyarray(B)
    if A.ndim != 2:
        raise ValueError('expected ndim to be 2, but observed %s' % A.ndim)
    if A.shape != B.shape:
        raise ValueError('the shapes of A and B differ (%s vs %s)' % (
            A.shape, B.shape))
    # Be clever with transposes, with the intention to save memory.
    # The conjugate has no effect for real inputs, but gives the correct solution for complex inputs.
    u, w, vt = svd((B.T@np.conjugate(A)).T)
    R = u@vt
    scale = w.sum()
    return R, scale
