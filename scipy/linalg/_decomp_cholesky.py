"""Cholesky decomposition functions."""

import numpy as np
from numpy import asarray_chkfinite, asarray, atleast_2d, empty_like

# Local imports
from ._misc import LinAlgError, _datacopied
from .lapack import get_lapack_funcs
from ._decomp import _asarray_validated


__all__ = ['cholesky', 'cho_factor', 'cho_solve', 'cholesky_banded',
           'cho_solve_banded', 'cholesky_update']


def _cholesky(a, lower=False, overwrite_a=False, clean=True,
              check_finite=True):
    """Common code for cholesky() and cho_factor()."""

    a1 = asarray_chkfinite(a) if check_finite else asarray(a)
    a1 = atleast_2d(a1)

    # Dimension check
    if a1.ndim != 2:
        raise ValueError(f'Input array needs to be 2D but received a {a1.ndim}d-array.')
    # Squareness check
    if a1.shape[0] != a1.shape[1]:
        raise ValueError('Input array is expected to be square but has '
                         f'the shape: {a1.shape}.')

    # Quick return for square empty array
    if a1.size == 0:
        dt = cholesky(np.eye(1, dtype=a1.dtype)).dtype
        return empty_like(a1, dtype=dt), lower

    overwrite_a = overwrite_a or _datacopied(a1, a)
    potrf, = get_lapack_funcs(('potrf',), (a1,))
    c, info = potrf(a1, lower=lower, overwrite_a=overwrite_a, clean=clean)
    if info > 0:
        raise LinAlgError("%d-th leading minor of the array is not positive "
                          "definite" % info)
    if info < 0:
        raise ValueError(f'LAPACK reported an illegal value in {-info}-th argument'
                         'on entry to "POTRF".')
    return c, lower


def cholesky(a, lower=False, overwrite_a=False, check_finite=True):
    """
    Compute the Cholesky decomposition of a matrix.

    Returns the Cholesky decomposition, :math:`A = L L^*` or
    :math:`A = U^* U` of a Hermitian positive-definite matrix A.

    Parameters
    ----------
    a : (M, M) array_like
        Matrix to be decomposed
    lower : bool, optional
        Whether to compute the upper- or lower-triangular Cholesky
        factorization.  Default is upper-triangular.
    overwrite_a : bool, optional
        Whether to overwrite data in `a` (may improve performance).
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    c : (M, M) ndarray
        Upper- or lower-triangular Cholesky factor of `a`.

    Raises
    ------
    LinAlgError : if decomposition fails.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import cholesky
    >>> a = np.array([[1,-2j],[2j,5]])
    >>> L = cholesky(a, lower=True)
    >>> L
    array([[ 1.+0.j,  0.+0.j],
           [ 0.+2.j,  1.+0.j]])
    >>> L @ L.T.conj()
    array([[ 1.+0.j,  0.-2.j],
           [ 0.+2.j,  5.+0.j]])

    """
    c, lower = _cholesky(a, lower=lower, overwrite_a=overwrite_a, clean=True,
                         check_finite=check_finite)
    return c


def cho_factor(a, lower=False, overwrite_a=False, check_finite=True):
    """
    Compute the Cholesky decomposition of a matrix, to use in cho_solve

    Returns a matrix containing the Cholesky decomposition,
    ``A = L L*`` or ``A = U* U`` of a Hermitian positive-definite matrix `a`.
    The return value can be directly used as the first parameter to cho_solve.

    .. warning::
        The returned matrix also contains random data in the entries not
        used by the Cholesky decomposition. If you need to zero these
        entries, use the function `cholesky` instead.

    Parameters
    ----------
    a : (M, M) array_like
        Matrix to be decomposed
    lower : bool, optional
        Whether to compute the upper or lower triangular Cholesky factorization
        (Default: upper-triangular)
    overwrite_a : bool, optional
        Whether to overwrite data in a (may improve performance)
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    c : (M, M) ndarray
        Matrix whose upper or lower triangle contains the Cholesky factor
        of `a`. Other parts of the matrix contain random data.
    lower : bool
        Flag indicating whether the factor is in the lower or upper triangle

    Raises
    ------
    LinAlgError
        Raised if decomposition fails.

    See Also
    --------
    cho_solve : Solve a linear set equations using the Cholesky factorization
                of a matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import cho_factor
    >>> A = np.array([[9, 3, 1, 5], [3, 7, 5, 1], [1, 5, 9, 2], [5, 1, 2, 6]])
    >>> c, low = cho_factor(A)
    >>> c
    array([[3.        , 1.        , 0.33333333, 1.66666667],
           [3.        , 2.44948974, 1.90515869, -0.27216553],
           [1.        , 5.        , 2.29330749, 0.8559528 ],
           [5.        , 1.        , 2.        , 1.55418563]])
    >>> np.allclose(np.triu(c).T @ np. triu(c) - A, np.zeros((4, 4)))
    True

    """
    c, lower = _cholesky(a, lower=lower, overwrite_a=overwrite_a, clean=False,
                         check_finite=check_finite)
    return c, lower


def cho_solve(c_and_lower, b, overwrite_b=False, check_finite=True):
    """Solve the linear equations A x = b, given the Cholesky factorization of A.

    Parameters
    ----------
    (c, lower) : tuple, (array, bool)
        Cholesky factorization of a, as given by cho_factor
    b : array
        Right-hand side
    overwrite_b : bool, optional
        Whether to overwrite data in b (may improve performance)
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    x : array
        The solution to the system A x = b

    See Also
    --------
    cho_factor : Cholesky factorization of a matrix

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import cho_factor, cho_solve
    >>> A = np.array([[9, 3, 1, 5], [3, 7, 5, 1], [1, 5, 9, 2], [5, 1, 2, 6]])
    >>> c, low = cho_factor(A)
    >>> x = cho_solve((c, low), [1, 1, 1, 1])
    >>> np.allclose(A @ x - [1, 1, 1, 1], np.zeros(4))
    True

    """
    (c, lower) = c_and_lower
    if check_finite:
        b1 = asarray_chkfinite(b)
        c = asarray_chkfinite(c)
    else:
        b1 = asarray(b)
        c = asarray(c)

    if c.ndim != 2 or c.shape[0] != c.shape[1]:
        raise ValueError("The factored matrix c is not square.")
    if c.shape[1] != b1.shape[0]:
        raise ValueError(f"incompatible dimensions ({c.shape} and {b1.shape})")

    # accommodate empty arrays
    if b1.size == 0:
        dt = cho_solve((np.eye(2, dtype=b1.dtype), True),
                        np.ones(2, dtype=c.dtype)).dtype
        return empty_like(b1, dtype=dt)

    overwrite_b = overwrite_b or _datacopied(b1, b)

    potrs, = get_lapack_funcs(('potrs',), (c, b1))
    x, info = potrs(c, b1, lower=lower, overwrite_b=overwrite_b)
    if info != 0:
        raise ValueError('illegal value in %dth argument of internal potrs'
                         % -info)
    return x


def cholesky_banded(ab, overwrite_ab=False, lower=False, check_finite=True):
    """
    Cholesky decompose a banded Hermitian positive-definite matrix

    The matrix a is stored in ab either in lower-diagonal or upper-
    diagonal ordered form::

        ab[u + i - j, j] == a[i,j]        (if upper form; i <= j)
        ab[    i - j, j] == a[i,j]        (if lower form; i >= j)

    Example of ab (shape of a is (6,6), u=2)::

        upper form:
        *   *   a02 a13 a24 a35
        *   a01 a12 a23 a34 a45
        a00 a11 a22 a33 a44 a55

        lower form:
        a00 a11 a22 a33 a44 a55
        a10 a21 a32 a43 a54 *
        a20 a31 a42 a53 *   *

    Parameters
    ----------
    ab : (u + 1, M) array_like
        Banded matrix
    overwrite_ab : bool, optional
        Discard data in ab (may enhance performance)
    lower : bool, optional
        Is the matrix in the lower form. (Default is upper form)
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    c : (u + 1, M) ndarray
        Cholesky factorization of a, in the same banded format as ab

    See Also
    --------
    cho_solve_banded :
        Solve a linear set equations, given the Cholesky factorization
        of a banded Hermitian.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import cholesky_banded
    >>> from numpy import allclose, zeros, diag
    >>> Ab = np.array([[0, 0, 1j, 2, 3j], [0, -1, -2, 3, 4], [9, 8, 7, 6, 9]])
    >>> A = np.diag(Ab[0,2:], k=2) + np.diag(Ab[1,1:], k=1)
    >>> A = A + A.conj().T + np.diag(Ab[2, :])
    >>> c = cholesky_banded(Ab)
    >>> C = np.diag(c[0, 2:], k=2) + np.diag(c[1, 1:], k=1) + np.diag(c[2, :])
    >>> np.allclose(C.conj().T @ C - A, np.zeros((5, 5)))
    True

    """
    if check_finite:
        ab = asarray_chkfinite(ab)
    else:
        ab = asarray(ab)

    # accommodate square empty matrices
    if ab.size == 0:
        dt = cholesky_banded(np.array([[0, 0], [1, 1]], dtype=ab.dtype)).dtype
        return empty_like(ab, dtype=dt)

    pbtrf, = get_lapack_funcs(('pbtrf',), (ab,))
    c, info = pbtrf(ab, lower=lower, overwrite_ab=overwrite_ab)
    if info > 0:
        raise LinAlgError("%d-th leading minor not positive definite" % info)
    if info < 0:
        raise ValueError('illegal value in %d-th argument of internal pbtrf'
                         % -info)
    return c


def cho_solve_banded(cb_and_lower, b, overwrite_b=False, check_finite=True):
    """
    Solve the linear equations ``A x = b``, given the Cholesky factorization of
    the banded Hermitian ``A``.

    Parameters
    ----------
    (cb, lower) : tuple, (ndarray, bool)
        `cb` is the Cholesky factorization of A, as given by cholesky_banded.
        `lower` must be the same value that was given to cholesky_banded.
    b : array_like
        Right-hand side
    overwrite_b : bool, optional
        If True, the function will overwrite the values in `b`.
    check_finite : bool, optional
        Whether to check that the input matrices contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    x : array
        The solution to the system A x = b

    See Also
    --------
    cholesky_banded : Cholesky factorization of a banded matrix

    Notes
    -----

    .. versionadded:: 0.8.0

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import cholesky_banded, cho_solve_banded
    >>> Ab = np.array([[0, 0, 1j, 2, 3j], [0, -1, -2, 3, 4], [9, 8, 7, 6, 9]])
    >>> A = np.diag(Ab[0,2:], k=2) + np.diag(Ab[1,1:], k=1)
    >>> A = A + A.conj().T + np.diag(Ab[2, :])
    >>> c = cholesky_banded(Ab)
    >>> x = cho_solve_banded((c, False), np.ones(5))
    >>> np.allclose(A @ x - np.ones(5), np.zeros(5))
    True

    """
    (cb, lower) = cb_and_lower
    if check_finite:
        cb = asarray_chkfinite(cb)
        b = asarray_chkfinite(b)
    else:
        cb = asarray(cb)
        b = asarray(b)

    # Validate shapes.
    if cb.shape[-1] != b.shape[0]:
        raise ValueError("shapes of cb and b are not compatible.")

    # accommodate empty arrays
    if b.size == 0:
        m = cholesky_banded(np.array([[0, 0], [1, 1]], dtype=cb.dtype))
        dt = cho_solve_banded((m, True), np.ones(2, dtype=b.dtype)).dtype
        return empty_like(b, dtype=dt)

    pbtrs, = get_lapack_funcs(('pbtrs',), (cb, b))
    x, info = pbtrs(cb, b, lower=lower, overwrite_b=overwrite_b)
    if info > 0:
        raise LinAlgError("%dth leading minor not positive definite" % info)
    if info < 0:
        raise ValueError('illegal value in %dth argument of internal pbtrs'
                         % -info)
    return x


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
       BIT Numerical Mathematics, Vol.30(4), 1990, :doi:`10.1007/BF01933218`

    """
    # Input validation
    R1 = _asarray_validated(R, check_finite=check_finite)
    z1 = _asarray_validated(z, check_finite=check_finite)

    overwrite_R = overwrite_R or _datacopied(R1, R)
    overwrite_z = overwrite_z or _datacopied(z1, z)

    # Dimension check
    if len(R1.shape) != 2:
        raise ValueError("Expected 2D array to be updated.")
    if len(z1.squeeze().shape) != 1:
        raise ValueError("Expected 1D update vector.")

    # Square matrix check
    m, n = R1.shape
    if m != n:
        raise ValueError('Input needs to be a square matrix.')

    # Compatible dimensions check
    if m != len(z1):
            raise ValueError('Input z has to have same length as the number '
                             'of rows as input R.')
    # Copy or not
    R = R1 if overwrite_R else R1.copy()
    z = z1 if overwrite_z else z1.copy()

    # Initializations
    if eps is None:
        eps = n * np.spacing(1.)  # For complex this needs modification
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
        beta[r] = np.sqrt(alpha[r])
        z[r + 1:] -= a * R[r, r + 1:]
        R[r, r:] *= beta[r] / beta[r - 1]
        R[r, r + 1:] += sign * a / (beta[r] * beta[r - 1]) * z[r + 1:]

    return R if not lower else R.T
