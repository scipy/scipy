import math
import warnings

import numpy as np
from numpy.lib.stride_tricks import as_strided

__all__ = ['tri', 'tril', 'triu', 'toeplitz', 'circulant', 'hankel',
           'hadamard', 'leslie', 'kron', 'block_diag', 'companion',
           'helmert', 'hilbert', 'invhilbert', 'pascal', 'invpascal', 'dft',
           'fiedler', 'fiedler_companion', 'convolution_matrix',
           'sakurai', 'mikota_pair', 'laplacian'
          ]


# -----------------------------------------------------------------------------
#  matrix construction functions
# -----------------------------------------------------------------------------

#
# *Note*: tri{,u,l} is implemented in NumPy, but an important bug was fixed in
# 2.0.0.dev-1af2f3, the following tri{,u,l} definitions are here for backwards
# compatibility.

def tri(N, M=None, k=0, dtype=None):
    """
    .. deprecated:: 1.11.0
        `tri` is deprecated in favour of `numpy.tri` and will be removed in
        SciPy 1.13.0.    
    
    Construct (N, M) matrix filled with ones at and below the kth diagonal.

    The matrix has A[i,j] == 1 for j <= i + k

    Parameters
    ----------
    N : int
        The size of the first dimension of the matrix.
    M : int or None, optional
        The size of the second dimension of the matrix. If `M` is None,
        `M = N` is assumed.
    k : int, optional
        Number of subdiagonal below which matrix is filled with ones.
        `k` = 0 is the main diagonal, `k` < 0 subdiagonal and `k` > 0
        superdiagonal.
    dtype : dtype, optional
        Data type of the matrix.

    Returns
    -------
    tri : (N, M) ndarray
        Tri matrix.

    Examples
    --------
    >>> from scipy.linalg import tri
    >>> tri(3, 5, 2, dtype=int)
    array([[1, 1, 1, 0, 0],
           [1, 1, 1, 1, 0],
           [1, 1, 1, 1, 1]])
    >>> tri(3, 5, -1, dtype=int)
    array([[0, 0, 0, 0, 0],
           [1, 0, 0, 0, 0],
           [1, 1, 0, 0, 0]])

    """
    warnings.warn("'tri'/'tril/'triu' are deprecated as of SciPy 1.11.0 and "
                  "will be removed in v1.13.0. Please use "
                  "numpy.(tri/tril/triu) instead.",
                  DeprecationWarning, stacklevel=2)
    
    if M is None:
        M = N
    if isinstance(M, str):
        # pearu: any objections to remove this feature?
        #       As tri(N,'d') is equivalent to tri(N,dtype='d')
        dtype = M
        M = N
    m = np.greater_equal.outer(np.arange(k, N+k), np.arange(M))
    if dtype is None:
        return m
    else:
        return m.astype(dtype)


def tril(m, k=0):
    """
    .. deprecated:: 1.11.0
        `tril` is deprecated in favour of `numpy.tril` and will be removed in
        SciPy 1.13.0.

    Make a copy of a matrix with elements above the kth diagonal zeroed.

    Parameters
    ----------
    m : array_like
        Matrix whose elements to return
    k : int, optional
        Diagonal above which to zero elements.
        `k` == 0 is the main diagonal, `k` < 0 subdiagonal and
        `k` > 0 superdiagonal.

    Returns
    -------
    tril : ndarray
        Return is the same shape and type as `m`.

    Examples
    --------
    >>> from scipy.linalg import tril
    >>> tril([[1,2,3],[4,5,6],[7,8,9],[10,11,12]], -1)
    array([[ 0,  0,  0],
           [ 4,  0,  0],
           [ 7,  8,  0],
           [10, 11, 12]])

    """
    m = np.asarray(m)
    out = tri(m.shape[0], m.shape[1], k=k, dtype=m.dtype.char) * m
    return out


def triu(m, k=0):
    """
    .. deprecated:: 1.11.0
        `tril` is deprecated in favour of `numpy.triu` and will be removed in
        SciPy 1.13.0.

    Make a copy of a matrix with elements below the kth diagonal zeroed.

    Parameters
    ----------
    m : array_like
        Matrix whose elements to return
    k : int, optional
        Diagonal below which to zero elements.
        `k` == 0 is the main diagonal, `k` < 0 subdiagonal and
        `k` > 0 superdiagonal.

    Returns
    -------
    triu : ndarray
        Return matrix with zeroed elements below the kth diagonal and has
        same shape and type as `m`.

    Examples
    --------
    >>> from scipy.linalg import triu
    >>> triu([[1,2,3],[4,5,6],[7,8,9],[10,11,12]], -1)
    array([[ 1,  2,  3],
           [ 4,  5,  6],
           [ 0,  8,  9],
           [ 0,  0, 12]])

    """
    m = np.asarray(m)
    out = (1 - tri(m.shape[0], m.shape[1], k - 1, m.dtype.char)) * m
    return out


def toeplitz(c, r=None):
    """
    Construct a Toeplitz matrix.

    The Toeplitz matrix has constant diagonals, with c as its first column
    and r as its first row. If r is not given, ``r == conjugate(c)`` is
    assumed.

    Parameters
    ----------
    c : array_like
        First column of the matrix.  Whatever the actual shape of `c`, it
        will be converted to a 1-D array.
    r : array_like, optional
        First row of the matrix. If None, ``r = conjugate(c)`` is assumed;
        in this case, if c[0] is real, the result is a Hermitian matrix.
        r[0] is ignored; the first row of the returned matrix is
        ``[c[0], r[1:]]``.  Whatever the actual shape of `r`, it will be
        converted to a 1-D array.

    Returns
    -------
    A : (len(c), len(r)) ndarray
        The Toeplitz matrix. Dtype is the same as ``(c[0] + r[0]).dtype``.

    See Also
    --------
    circulant : circulant matrix
    hankel : Hankel matrix
    solve_toeplitz : Solve a Toeplitz system.

    Notes
    -----
    The behavior when `c` or `r` is a scalar, or when `c` is complex and
    `r` is None, was changed in version 0.8.0. The behavior in previous
    versions was undocumented and is no longer supported.

    Examples
    --------
    >>> from scipy.linalg import toeplitz
    >>> toeplitz([1,2,3], [1,4,5,6])
    array([[1, 4, 5, 6],
           [2, 1, 4, 5],
           [3, 2, 1, 4]])
    >>> toeplitz([1.0, 2+3j, 4-1j])
    array([[ 1.+0.j,  2.-3.j,  4.+1.j],
           [ 2.+3.j,  1.+0.j,  2.-3.j],
           [ 4.-1.j,  2.+3.j,  1.+0.j]])

    """
    c = np.asarray(c).ravel()
    if r is None:
        r = c.conjugate()
    else:
        r = np.asarray(r).ravel()
    # Form a 1-D array containing a reversed c followed by r[1:] that could be
    # strided to give us toeplitz matrix.
    vals = np.concatenate((c[::-1], r[1:]))
    out_shp = len(c), len(r)
    n = vals.strides[0]
    return as_strided(vals[len(c)-1:], shape=out_shp, strides=(-n, n)).copy()


def circulant(c):
    """
    Construct a circulant matrix.

    Parameters
    ----------
    c : (N,) array_like
        1-D array, the first column of the matrix.

    Returns
    -------
    A : (N, N) ndarray
        A circulant matrix whose first column is `c`.

    See Also
    --------
    toeplitz : Toeplitz matrix
    hankel : Hankel matrix
    solve_circulant : Solve a circulant system.

    Notes
    -----
    .. versionadded:: 0.8.0

    Examples
    --------
    >>> from scipy.linalg import circulant
    >>> circulant([1, 2, 3])
    array([[1, 3, 2],
           [2, 1, 3],
           [3, 2, 1]])

    """
    c = np.asarray(c).ravel()
    # Form an extended array that could be strided to give circulant version
    c_ext = np.concatenate((c[::-1], c[:0:-1]))
    L = len(c)
    n = c_ext.strides[0]
    return as_strided(c_ext[L-1:], shape=(L, L), strides=(-n, n)).copy()


def hankel(c, r=None):
    """
    Construct a Hankel matrix.

    The Hankel matrix has constant anti-diagonals, with `c` as its
    first column and `r` as its last row. If `r` is not given, then
    `r = zeros_like(c)` is assumed.

    Parameters
    ----------
    c : array_like
        First column of the matrix. Whatever the actual shape of `c`, it
        will be converted to a 1-D array.
    r : array_like, optional
        Last row of the matrix. If None, ``r = zeros_like(c)`` is assumed.
        r[0] is ignored; the last row of the returned matrix is
        ``[c[-1], r[1:]]``. Whatever the actual shape of `r`, it will be
        converted to a 1-D array.

    Returns
    -------
    A : (len(c), len(r)) ndarray
        The Hankel matrix. Dtype is the same as ``(c[0] + r[0]).dtype``.

    See Also
    --------
    toeplitz : Toeplitz matrix
    circulant : circulant matrix

    Examples
    --------
    >>> from scipy.linalg import hankel
    >>> hankel([1, 17, 99])
    array([[ 1, 17, 99],
           [17, 99,  0],
           [99,  0,  0]])
    >>> hankel([1,2,3,4], [4,7,7,8,9])
    array([[1, 2, 3, 4, 7],
           [2, 3, 4, 7, 7],
           [3, 4, 7, 7, 8],
           [4, 7, 7, 8, 9]])

    """
    c = np.asarray(c).ravel()
    if r is None:
        r = np.zeros_like(c)
    else:
        r = np.asarray(r).ravel()
    # Form a 1-D array of values to be used in the matrix, containing `c`
    # followed by r[1:].
    vals = np.concatenate((c, r[1:]))
    # Stride on concatenated array to get hankel matrix
    out_shp = len(c), len(r)
    n = vals.strides[0]
    return as_strided(vals, shape=out_shp, strides=(n, n)).copy()


def hadamard(n, dtype=int):
    """
    Construct an Hadamard matrix.

    Constructs an n-by-n Hadamard matrix, using Sylvester's
    construction. `n` must be a power of 2.

    Parameters
    ----------
    n : int
        The order of the matrix. `n` must be a power of 2.
    dtype : dtype, optional
        The data type of the array to be constructed.

    Returns
    -------
    H : (n, n) ndarray
        The Hadamard matrix.

    Notes
    -----
    .. versionadded:: 0.8.0

    Examples
    --------
    >>> from scipy.linalg import hadamard
    >>> hadamard(2, dtype=complex)
    array([[ 1.+0.j,  1.+0.j],
           [ 1.+0.j, -1.-0.j]])
    >>> hadamard(4)
    array([[ 1,  1,  1,  1],
           [ 1, -1,  1, -1],
           [ 1,  1, -1, -1],
           [ 1, -1, -1,  1]])

    """

    # This function is a slightly modified version of the
    # function contributed by Ivo in ticket #675.

    if n < 1:
        lg2 = 0
    else:
        lg2 = int(math.log(n, 2))
    if 2 ** lg2 != n:
        raise ValueError("n must be an positive integer, and n must be "
                         "a power of 2")

    H = np.array([[1]], dtype=dtype)

    # Sylvester's construction
    for i in range(0, lg2):
        H = np.vstack((np.hstack((H, H)), np.hstack((H, -H))))

    return H


def leslie(f, s):
    """
    Create a Leslie matrix.

    Given the length n array of fecundity coefficients `f` and the length
    n-1 array of survival coefficients `s`, return the associated Leslie
    matrix.

    Parameters
    ----------
    f : (N,) array_like
        The "fecundity" coefficients.
    s : (N-1,) array_like
        The "survival" coefficients, has to be 1-D.  The length of `s`
        must be one less than the length of `f`, and it must be at least 1.

    Returns
    -------
    L : (N, N) ndarray
        The array is zero except for the first row,
        which is `f`, and the first sub-diagonal, which is `s`.
        The data-type of the array will be the data-type of ``f[0]+s[0]``.

    Notes
    -----
    .. versionadded:: 0.8.0

    The Leslie matrix is used to model discrete-time, age-structured
    population growth [1]_ [2]_. In a population with `n` age classes, two sets
    of parameters define a Leslie matrix: the `n` "fecundity coefficients",
    which give the number of offspring per-capita produced by each age
    class, and the `n` - 1 "survival coefficients", which give the
    per-capita survival rate of each age class.

    References
    ----------
    .. [1] P. H. Leslie, On the use of matrices in certain population
           mathematics, Biometrika, Vol. 33, No. 3, 183--212 (Nov. 1945)
    .. [2] P. H. Leslie, Some further notes on the use of matrices in
           population mathematics, Biometrika, Vol. 35, No. 3/4, 213--245
           (Dec. 1948)

    Examples
    --------
    >>> from scipy.linalg import leslie
    >>> leslie([0.1, 2.0, 1.0, 0.1], [0.2, 0.8, 0.7])
    array([[ 0.1,  2. ,  1. ,  0.1],
           [ 0.2,  0. ,  0. ,  0. ],
           [ 0. ,  0.8,  0. ,  0. ],
           [ 0. ,  0. ,  0.7,  0. ]])

    """
    f = np.atleast_1d(f)
    s = np.atleast_1d(s)
    if f.ndim != 1:
        raise ValueError("Incorrect shape for f.  f must be 1D")
    if s.ndim != 1:
        raise ValueError("Incorrect shape for s.  s must be 1D")
    if f.size != s.size + 1:
        raise ValueError("Incorrect lengths for f and s.  The length"
                         " of s must be one less than the length of f.")
    if s.size == 0:
        raise ValueError("The length of s must be at least 1.")

    tmp = f[0] + s[0]
    n = f.size
    a = np.zeros((n, n), dtype=tmp.dtype)
    a[0] = f
    a[list(range(1, n)), list(range(0, n - 1))] = s
    return a


def kron(a, b):
    """
    Kronecker product.

    The result is the block matrix::

        a[0,0]*b    a[0,1]*b  ... a[0,-1]*b
        a[1,0]*b    a[1,1]*b  ... a[1,-1]*b
        ...
        a[-1,0]*b   a[-1,1]*b ... a[-1,-1]*b

    Parameters
    ----------
    a : (M, N) ndarray
        Input array
    b : (P, Q) ndarray
        Input array

    Returns
    -------
    A : (M*P, N*Q) ndarray
        Kronecker product of `a` and `b`.

    Examples
    --------
    >>> from numpy import array
    >>> from scipy.linalg import kron
    >>> kron(array([[1,2],[3,4]]), array([[1,1,1]]))
    array([[1, 1, 1, 2, 2, 2],
           [3, 3, 3, 4, 4, 4]])

    """
    if not a.flags['CONTIGUOUS']:
        a = np.reshape(a, a.shape)
    if not b.flags['CONTIGUOUS']:
        b = np.reshape(b, b.shape)
    o = np.outer(a, b)
    o = o.reshape(a.shape + b.shape)
    return np.concatenate(np.concatenate(o, axis=1), axis=1)


def block_diag(*arrs):
    """
    Create a block diagonal matrix from provided arrays.

    Given the inputs `A`, `B` and `C`, the output will have these
    arrays arranged on the diagonal::

        [[A, 0, 0],
         [0, B, 0],
         [0, 0, C]]

    Parameters
    ----------
    A, B, C, ... : array_like, up to 2-D
        Input arrays.  A 1-D array or array_like sequence of length `n` is
        treated as a 2-D array with shape ``(1,n)``.

    Returns
    -------
    D : ndarray
        Array with `A`, `B`, `C`, ... on the diagonal. `D` has the
        same dtype as `A`.

    Notes
    -----
    If all the input arrays are square, the output is known as a
    block diagonal matrix.

    Empty sequences (i.e., array-likes of zero size) will not be ignored.
    Noteworthy, both [] and [[]] are treated as matrices with shape ``(1,0)``.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import block_diag
    >>> A = [[1, 0],
    ...      [0, 1]]
    >>> B = [[3, 4, 5],
    ...      [6, 7, 8]]
    >>> C = [[7]]
    >>> P = np.zeros((2, 0), dtype='int32')
    >>> block_diag(A, B, C)
    array([[1, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 0, 0],
           [0, 0, 3, 4, 5, 0],
           [0, 0, 6, 7, 8, 0],
           [0, 0, 0, 0, 0, 7]])
    >>> block_diag(A, P, B, C)
    array([[1, 0, 0, 0, 0, 0],
           [0, 1, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0],
           [0, 0, 3, 4, 5, 0],
           [0, 0, 6, 7, 8, 0],
           [0, 0, 0, 0, 0, 7]])
    >>> block_diag(1.0, [2, 3], [[4, 5], [6, 7]])
    array([[ 1.,  0.,  0.,  0.,  0.],
           [ 0.,  2.,  3.,  0.,  0.],
           [ 0.,  0.,  0.,  4.,  5.],
           [ 0.,  0.,  0.,  6.,  7.]])

    """
    if arrs == ():
        arrs = ([],)
    arrs = [np.atleast_2d(a) for a in arrs]

    bad_args = [k for k in range(len(arrs)) if arrs[k].ndim > 2]
    if bad_args:
        raise ValueError("arguments in the following positions have dimension "
                         "greater than 2: %s" % bad_args)

    shapes = np.array([a.shape for a in arrs])
    out_dtype = np.result_type(*[arr.dtype for arr in arrs])
    out = np.zeros(np.sum(shapes, axis=0), dtype=out_dtype)

    r, c = 0, 0
    for i, (rr, cc) in enumerate(shapes):
        out[r:r + rr, c:c + cc] = arrs[i]
        r += rr
        c += cc
    return out


def companion(a):
    """
    Create a companion matrix.

    Create the companion matrix [1]_ associated with the polynomial whose
    coefficients are given in `a`.

    Parameters
    ----------
    a : (N,) array_like
        1-D array of polynomial coefficients. The length of `a` must be
        at least two, and ``a[0]`` must not be zero.

    Returns
    -------
    c : (N-1, N-1) ndarray
        The first row of `c` is ``-a[1:]/a[0]``, and the first
        sub-diagonal is all ones.  The data-type of the array is the same
        as the data-type of ``1.0*a[0]``.

    Raises
    ------
    ValueError
        If any of the following are true: a) ``a.ndim != 1``;
        b) ``a.size < 2``; c) ``a[0] == 0``.

    Notes
    -----
    .. versionadded:: 0.8.0

    References
    ----------
    .. [1] R. A. Horn & C. R. Johnson, *Matrix Analysis*.  Cambridge, UK:
        Cambridge University Press, 1999, pp. 146-7.

    Examples
    --------
    >>> from scipy.linalg import companion
    >>> companion([1, -10, 31, -30])
    array([[ 10., -31.,  30.],
           [  1.,   0.,   0.],
           [  0.,   1.,   0.]])

    """
    a = np.atleast_1d(a)

    if a.ndim != 1:
        raise ValueError("Incorrect shape for `a`.  `a` must be "
                         "one-dimensional.")

    if a.size < 2:
        raise ValueError("The length of `a` must be at least 2.")

    if a[0] == 0:
        raise ValueError("The first coefficient in `a` must not be zero.")

    first_row = -a[1:] / (1.0 * a[0])
    n = a.size
    c = np.zeros((n - 1, n - 1), dtype=first_row.dtype)
    c[0] = first_row
    c[list(range(1, n - 1)), list(range(0, n - 2))] = 1
    return c


def helmert(n, full=False):
    """
    Create an Helmert matrix of order `n`.

    This has applications in statistics, compositional or simplicial analysis,
    and in Aitchison geometry.

    Parameters
    ----------
    n : int
        The size of the array to create.
    full : bool, optional
        If True the (n, n) ndarray will be returned.
        Otherwise the submatrix that does not include the first
        row will be returned.
        Default: False.

    Returns
    -------
    M : ndarray
        The Helmert matrix.
        The shape is (n, n) or (n-1, n) depending on the `full` argument.

    Examples
    --------
    >>> from scipy.linalg import helmert
    >>> helmert(5, full=True)
    array([[ 0.4472136 ,  0.4472136 ,  0.4472136 ,  0.4472136 ,  0.4472136 ],
           [ 0.70710678, -0.70710678,  0.        ,  0.        ,  0.        ],
           [ 0.40824829,  0.40824829, -0.81649658,  0.        ,  0.        ],
           [ 0.28867513,  0.28867513,  0.28867513, -0.8660254 ,  0.        ],
           [ 0.2236068 ,  0.2236068 ,  0.2236068 ,  0.2236068 , -0.89442719]])

    """
    H = np.tril(np.ones((n, n)), -1) - np.diag(np.arange(n))
    d = np.arange(n) * np.arange(1, n+1)
    H[0] = 1
    d[0] = n
    H_full = H / np.sqrt(d)[:, np.newaxis]
    if full:
        return H_full
    else:
        return H_full[1:]


def hilbert(n):
    """
    Create a Hilbert matrix of order `n`.

    Returns the `n` by `n` array with entries `h[i,j] = 1 / (i + j + 1)`.

    Parameters
    ----------
    n : int
        The size of the array to create.

    Returns
    -------
    h : (n, n) ndarray
        The Hilbert matrix.

    See Also
    --------
    invhilbert : Compute the inverse of a Hilbert matrix.

    Notes
    -----
    .. versionadded:: 0.10.0

    Examples
    --------
    >>> from scipy.linalg import hilbert
    >>> hilbert(3)
    array([[ 1.        ,  0.5       ,  0.33333333],
           [ 0.5       ,  0.33333333,  0.25      ],
           [ 0.33333333,  0.25      ,  0.2       ]])

    """
    values = 1.0 / (1.0 + np.arange(2 * n - 1))
    h = hankel(values[:n], r=values[n - 1:])
    return h


def invhilbert(n, exact=False):
    """
    Compute the inverse of the Hilbert matrix of order `n`.

    The entries in the inverse of a Hilbert matrix are integers. When `n`
    is greater than 14, some entries in the inverse exceed the upper limit
    of 64 bit integers. The `exact` argument provides two options for
    dealing with these large integers.

    Parameters
    ----------
    n : int
        The order of the Hilbert matrix.
    exact : bool, optional
        If False, the data type of the array that is returned is np.float64,
        and the array is an approximation of the inverse.
        If True, the array is the exact integer inverse array. To represent
        the exact inverse when n > 14, the returned array is an object array
        of long integers. For n <= 14, the exact inverse is returned as an
        array with data type np.int64.

    Returns
    -------
    invh : (n, n) ndarray
        The data type of the array is np.float64 if `exact` is False.
        If `exact` is True, the data type is either np.int64 (for n <= 14)
        or object (for n > 14). In the latter case, the objects in the
        array will be long integers.

    See Also
    --------
    hilbert : Create a Hilbert matrix.

    Notes
    -----
    .. versionadded:: 0.10.0

    Examples
    --------
    >>> from scipy.linalg import invhilbert
    >>> invhilbert(4)
    array([[   16.,  -120.,   240.,  -140.],
           [ -120.,  1200., -2700.,  1680.],
           [  240., -2700.,  6480., -4200.],
           [ -140.,  1680., -4200.,  2800.]])
    >>> invhilbert(4, exact=True)
    array([[   16,  -120,   240,  -140],
           [ -120,  1200, -2700,  1680],
           [  240, -2700,  6480, -4200],
           [ -140,  1680, -4200,  2800]], dtype=int64)
    >>> invhilbert(16)[7,7]
    4.2475099528537506e+19
    >>> invhilbert(16, exact=True)[7,7]
    42475099528537378560

    """
    from scipy.special import comb
    if exact:
        if n > 14:
            dtype = object
        else:
            dtype = np.int64
    else:
        dtype = np.float64
    invh = np.empty((n, n), dtype=dtype)
    for i in range(n):
        for j in range(0, i + 1):
            s = i + j
            invh[i, j] = ((-1) ** s * (s + 1) *
                          comb(n + i, n - j - 1, exact) *
                          comb(n + j, n - i - 1, exact) *
                          comb(s, i, exact) ** 2)
            if i != j:
                invh[j, i] = invh[i, j]
    return invh


def pascal(n, kind='symmetric', exact=True):
    """
    Returns the n x n Pascal matrix.

    The Pascal matrix is a matrix containing the binomial coefficients as
    its elements.

    Parameters
    ----------
    n : int
        The size of the matrix to create; that is, the result is an n x n
        matrix.
    kind : str, optional
        Must be one of 'symmetric', 'lower', or 'upper'.
        Default is 'symmetric'.
    exact : bool, optional
        If `exact` is True, the result is either an array of type
        numpy.uint64 (if n < 35) or an object array of Python long integers.
        If `exact` is False, the coefficients in the matrix are computed using
        `scipy.special.comb` with `exact=False`. The result will be a floating
        point array, and the values in the array will not be the exact
        coefficients, but this version is much faster than `exact=True`.

    Returns
    -------
    p : (n, n) ndarray
        The Pascal matrix.

    See Also
    --------
    invpascal

    Notes
    -----
    See https://en.wikipedia.org/wiki/Pascal_matrix for more information
    about Pascal matrices.

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> from scipy.linalg import pascal
    >>> pascal(4)
    array([[ 1,  1,  1,  1],
           [ 1,  2,  3,  4],
           [ 1,  3,  6, 10],
           [ 1,  4, 10, 20]], dtype=uint64)
    >>> pascal(4, kind='lower')
    array([[1, 0, 0, 0],
           [1, 1, 0, 0],
           [1, 2, 1, 0],
           [1, 3, 3, 1]], dtype=uint64)
    >>> pascal(50)[-1, -1]
    25477612258980856902730428600
    >>> from scipy.special import comb
    >>> comb(98, 49, exact=True)
    25477612258980856902730428600

    """

    from scipy.special import comb
    if kind not in ['symmetric', 'lower', 'upper']:
        raise ValueError("kind must be 'symmetric', 'lower', or 'upper'")

    if exact:
        if n >= 35:
            L_n = np.empty((n, n), dtype=object)
            L_n.fill(0)
        else:
            L_n = np.zeros((n, n), dtype=np.uint64)
        for i in range(n):
            for j in range(i + 1):
                L_n[i, j] = comb(i, j, exact=True)
    else:
        L_n = comb(*np.ogrid[:n, :n])

    if kind == 'lower':
        p = L_n
    elif kind == 'upper':
        p = L_n.T
    else:
        p = np.dot(L_n, L_n.T)

    return p


def invpascal(n, kind='symmetric', exact=True):
    """
    Returns the inverse of the n x n Pascal matrix.

    The Pascal matrix is a matrix containing the binomial coefficients as
    its elements.

    Parameters
    ----------
    n : int
        The size of the matrix to create; that is, the result is an n x n
        matrix.
    kind : str, optional
        Must be one of 'symmetric', 'lower', or 'upper'.
        Default is 'symmetric'.
    exact : bool, optional
        If `exact` is True, the result is either an array of type
        ``numpy.int64`` (if `n` <= 35) or an object array of Python integers.
        If `exact` is False, the coefficients in the matrix are computed using
        `scipy.special.comb` with `exact=False`. The result will be a floating
        point array, and for large `n`, the values in the array will not be the
        exact coefficients.

    Returns
    -------
    invp : (n, n) ndarray
        The inverse of the Pascal matrix.

    See Also
    --------
    pascal

    Notes
    -----

    .. versionadded:: 0.16.0

    References
    ----------
    .. [1] "Pascal matrix", https://en.wikipedia.org/wiki/Pascal_matrix
    .. [2] Cohen, A. M., "The inverse of a Pascal matrix", Mathematical
           Gazette, 59(408), pp. 111-112, 1975.

    Examples
    --------
    >>> from scipy.linalg import invpascal, pascal
    >>> invp = invpascal(5)
    >>> invp
    array([[  5, -10,  10,  -5,   1],
           [-10,  30, -35,  19,  -4],
           [ 10, -35,  46, -27,   6],
           [ -5,  19, -27,  17,  -4],
           [  1,  -4,   6,  -4,   1]])

    >>> p = pascal(5)
    >>> p.dot(invp)
    array([[ 1.,  0.,  0.,  0.,  0.],
           [ 0.,  1.,  0.,  0.,  0.],
           [ 0.,  0.,  1.,  0.,  0.],
           [ 0.,  0.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  0.,  1.]])

    An example of the use of `kind` and `exact`:

    >>> invpascal(5, kind='lower', exact=False)
    array([[ 1., -0.,  0., -0.,  0.],
           [-1.,  1., -0.,  0., -0.],
           [ 1., -2.,  1., -0.,  0.],
           [-1.,  3., -3.,  1., -0.],
           [ 1., -4.,  6., -4.,  1.]])

    """
    from scipy.special import comb

    if kind not in ['symmetric', 'lower', 'upper']:
        raise ValueError("'kind' must be 'symmetric', 'lower' or 'upper'.")

    if kind == 'symmetric':
        if exact:
            if n > 34:
                dt = object
            else:
                dt = np.int64
        else:
            dt = np.float64
        invp = np.empty((n, n), dtype=dt)
        for i in range(n):
            for j in range(0, i + 1):
                v = 0
                for k in range(n - i):
                    v += comb(i + k, k, exact=exact) * comb(i + k, i + k - j,
                                                            exact=exact)
                invp[i, j] = (-1)**(i - j) * v
                if i != j:
                    invp[j, i] = invp[i, j]
    else:
        # For the 'lower' and 'upper' cases, we computer the inverse by
        # changing the sign of every other diagonal of the pascal matrix.
        invp = pascal(n, kind=kind, exact=exact)
        if invp.dtype == np.uint64:
            # This cast from np.uint64 to int64 OK, because if `kind` is not
            # "symmetric", the values in invp are all much less than 2**63.
            invp = invp.view(np.int64)

        # The toeplitz matrix has alternating bands of 1 and -1.
        invp *= toeplitz((-1)**np.arange(n)).astype(invp.dtype)

    return invp


def dft(n, scale=None):
    """
    Discrete Fourier transform matrix.

    Create the matrix that computes the discrete Fourier transform of a
    sequence [1]_. The nth primitive root of unity used to generate the
    matrix is exp(-2*pi*i/n), where i = sqrt(-1).

    Parameters
    ----------
    n : int
        Size the matrix to create.
    scale : str, optional
        Must be None, 'sqrtn', or 'n'.
        If `scale` is 'sqrtn', the matrix is divided by `sqrt(n)`.
        If `scale` is 'n', the matrix is divided by `n`.
        If `scale` is None (the default), the matrix is not normalized, and the
        return value is simply the Vandermonde matrix of the roots of unity.

    Returns
    -------
    m : (n, n) ndarray
        The DFT matrix.

    Notes
    -----
    When `scale` is None, multiplying a vector by the matrix returned by
    `dft` is mathematically equivalent to (but much less efficient than)
    the calculation performed by `scipy.fft.fft`.

    .. versionadded:: 0.14.0

    References
    ----------
    .. [1] "DFT matrix", https://en.wikipedia.org/wiki/DFT_matrix

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import dft
    >>> np.set_printoptions(precision=2, suppress=True)  # for compact output
    >>> m = dft(5)
    >>> m
    array([[ 1.  +0.j  ,  1.  +0.j  ,  1.  +0.j  ,  1.  +0.j  ,  1.  +0.j  ],
           [ 1.  +0.j  ,  0.31-0.95j, -0.81-0.59j, -0.81+0.59j,  0.31+0.95j],
           [ 1.  +0.j  , -0.81-0.59j,  0.31+0.95j,  0.31-0.95j, -0.81+0.59j],
           [ 1.  +0.j  , -0.81+0.59j,  0.31-0.95j,  0.31+0.95j, -0.81-0.59j],
           [ 1.  +0.j  ,  0.31+0.95j, -0.81+0.59j, -0.81-0.59j,  0.31-0.95j]])
    >>> x = np.array([1, 2, 3, 0, 3])
    >>> m @ x  # Compute the DFT of x
    array([ 9.  +0.j  ,  0.12-0.81j, -2.12+3.44j, -2.12-3.44j,  0.12+0.81j])

    Verify that ``m @ x`` is the same as ``fft(x)``.

    >>> from scipy.fft import fft
    >>> fft(x)     # Same result as m @ x
    array([ 9.  +0.j  ,  0.12-0.81j, -2.12+3.44j, -2.12-3.44j,  0.12+0.81j])
    """
    if scale not in [None, 'sqrtn', 'n']:
        raise ValueError("scale must be None, 'sqrtn', or 'n'; "
                         "{!r} is not valid.".format(scale))

    omegas = np.exp(-2j * np.pi * np.arange(n) / n).reshape(-1, 1)
    m = omegas ** np.arange(n)
    if scale == 'sqrtn':
        m /= math.sqrt(n)
    elif scale == 'n':
        m /= n
    return m


def fiedler(a):
    """Returns a symmetric Fiedler matrix

    Given an sequence of numbers `a`, Fiedler matrices have the structure
    ``F[i, j] = np.abs(a[i] - a[j])``, and hence zero diagonals and nonnegative
    entries. A Fiedler matrix has a dominant positive eigenvalue and other
    eigenvalues are negative. Although not valid generally, for certain inputs,
    the inverse and the determinant can be derived explicitly as given in [1]_.

    Parameters
    ----------
    a : (n,) array_like
        coefficient array

    Returns
    -------
    F : (n, n) ndarray

    See Also
    --------
    circulant, toeplitz

    Notes
    -----

    .. versionadded:: 1.3.0

    References
    ----------
    .. [1] J. Todd, "Basic Numerical Mathematics: Vol.2 : Numerical Algebra",
        1977, Birkhauser, :doi:`10.1007/978-3-0348-7286-7`

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import det, inv, fiedler
    >>> a = [1, 4, 12, 45, 77]
    >>> n = len(a)
    >>> A = fiedler(a)
    >>> A
    array([[ 0,  3, 11, 44, 76],
           [ 3,  0,  8, 41, 73],
           [11,  8,  0, 33, 65],
           [44, 41, 33,  0, 32],
           [76, 73, 65, 32,  0]])

    The explicit formulas for determinant and inverse seem to hold only for
    monotonically increasing/decreasing arrays. Note the tridiagonal structure
    and the corners.

    >>> Ai = inv(A)
    >>> Ai[np.abs(Ai) < 1e-12] = 0.  # cleanup the numerical noise for display
    >>> Ai
    array([[-0.16008772,  0.16666667,  0.        ,  0.        ,  0.00657895],
           [ 0.16666667, -0.22916667,  0.0625    ,  0.        ,  0.        ],
           [ 0.        ,  0.0625    , -0.07765152,  0.01515152,  0.        ],
           [ 0.        ,  0.        ,  0.01515152, -0.03077652,  0.015625  ],
           [ 0.00657895,  0.        ,  0.        ,  0.015625  , -0.00904605]])
    >>> det(A)
    15409151.999999998
    >>> (-1)**(n-1) * 2**(n-2) * np.diff(a).prod() * (a[-1] - a[0])
    15409152

    """
    a = np.atleast_1d(a)

    if a.ndim != 1:
        raise ValueError("Input 'a' must be a 1D array.")

    if a.size == 0:
        return np.array([], dtype=float)
    elif a.size == 1:
        return np.array([[0.]])
    else:
        return np.abs(a[:, None] - a)


def fiedler_companion(a):
    """ Returns a Fiedler companion matrix

    Given a polynomial coefficient array ``a``, this function forms a
    pentadiagonal matrix with a special structure whose eigenvalues coincides
    with the roots of ``a``.

    Parameters
    ----------
    a : (N,) array_like
        1-D array of polynomial coefficients in descending order with a nonzero
        leading coefficient. For ``N < 2``, an empty array is returned.

    Returns
    -------
    c : (N-1, N-1) ndarray
        Resulting companion matrix

    See Also
    --------
    companion

    Notes
    -----
    Similar to `companion` the leading coefficient should be nonzero. In the case
    the leading coefficient is not 1, other coefficients are rescaled before
    the array generation. To avoid numerical issues, it is best to provide a
    monic polynomial.

    .. versionadded:: 1.3.0

    References
    ----------
    .. [1] M. Fiedler, " A note on companion matrices", Linear Algebra and its
        Applications, 2003, :doi:`10.1016/S0024-3795(03)00548-2`

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import fiedler_companion, eigvals
    >>> p = np.poly(np.arange(1, 9, 2))  # [1., -16., 86., -176., 105.]
    >>> fc = fiedler_companion(p)
    >>> fc
    array([[  16.,  -86.,    1.,    0.],
           [   1.,    0.,    0.,    0.],
           [   0.,  176.,    0., -105.],
           [   0.,    1.,    0.,    0.]])
    >>> eigvals(fc)
    array([7.+0.j, 5.+0.j, 3.+0.j, 1.+0.j])

    """
    a = np.atleast_1d(a)

    if a.ndim != 1:
        raise ValueError("Input 'a' must be a 1-D array.")

    if a.size <= 2:
        if a.size == 2:
            return np.array([[-(a/a[0])[-1]]])
        return np.array([], dtype=a.dtype)

    if a[0] == 0.:
        raise ValueError('Leading coefficient is zero.')

    a = a/a[0]
    n = a.size - 1
    c = np.zeros((n, n), dtype=a.dtype)
    # subdiagonals
    c[range(3, n, 2), range(1, n-2, 2)] = 1.
    c[range(2, n, 2), range(1, n-1, 2)] = -a[3::2]
    # superdiagonals
    c[range(0, n-2, 2), range(2, n, 2)] = 1.
    c[range(0, n-1, 2), range(1, n, 2)] = -a[2::2]
    c[[0, 1], 0] = [-a[1], 1]

    return c


def convolution_matrix(a, n, mode='full'):
    """
    Construct a convolution matrix.

    Constructs the Toeplitz matrix representing one-dimensional
    convolution [1]_.  See the notes below for details.

    Parameters
    ----------
    a : (m,) array_like
        The 1-D array to convolve.
    n : int
        The number of columns in the resulting matrix.  It gives the length
        of the input to be convolved with `a`.  This is analogous to the
        length of `v` in ``numpy.convolve(a, v)``.
    mode : str
        This is analogous to `mode` in ``numpy.convolve(v, a, mode)``.
        It must be one of ('full', 'valid', 'same').
        See below for how `mode` determines the shape of the result.

    Returns
    -------
    A : (k, n) ndarray
        The convolution matrix whose row count `k` depends on `mode`::

            =======  =========================
             mode    k
            =======  =========================
            'full'   m + n -1
            'same'   max(m, n)
            'valid'  max(m, n) - min(m, n) + 1
            =======  =========================

    See Also
    --------
    toeplitz : Toeplitz matrix

    Notes
    -----
    The code::

        A = convolution_matrix(a, n, mode)

    creates a Toeplitz matrix `A` such that ``A @ v`` is equivalent to
    using ``convolve(a, v, mode)``.  The returned array always has `n`
    columns.  The number of rows depends on the specified `mode`, as
    explained above.

    In the default 'full' mode, the entries of `A` are given by::

        A[i, j] == (a[i-j] if (0 <= (i-j) < m) else 0)

    where ``m = len(a)``.  Suppose, for example, the input array is
    ``[x, y, z]``.  The convolution matrix has the form::

        [x, 0, 0, ..., 0, 0]
        [y, x, 0, ..., 0, 0]
        [z, y, x, ..., 0, 0]
        ...
        [0, 0, 0, ..., x, 0]
        [0, 0, 0, ..., y, x]
        [0, 0, 0, ..., z, y]
        [0, 0, 0, ..., 0, z]

    In 'valid' mode, the entries of `A` are given by::

        A[i, j] == (a[i-j+m-1] if (0 <= (i-j+m-1) < m) else 0)

    This corresponds to a matrix whose rows are the subset of those from
    the 'full' case where all the coefficients in `a` are contained in the
    row.  For input ``[x, y, z]``, this array looks like::

        [z, y, x, 0, 0, ..., 0, 0, 0]
        [0, z, y, x, 0, ..., 0, 0, 0]
        [0, 0, z, y, x, ..., 0, 0, 0]
        ...
        [0, 0, 0, 0, 0, ..., x, 0, 0]
        [0, 0, 0, 0, 0, ..., y, x, 0]
        [0, 0, 0, 0, 0, ..., z, y, x]

    In the 'same' mode, the entries of `A` are given by::

        d = (m - 1) // 2
        A[i, j] == (a[i-j+d] if (0 <= (i-j+d) < m) else 0)

    The typical application of the 'same' mode is when one has a signal of
    length `n` (with `n` greater than ``len(a)``), and the desired output
    is a filtered signal that is still of length `n`.

    For input ``[x, y, z]``, this array looks like::

        [y, x, 0, 0, ..., 0, 0, 0]
        [z, y, x, 0, ..., 0, 0, 0]
        [0, z, y, x, ..., 0, 0, 0]
        [0, 0, z, y, ..., 0, 0, 0]
        ...
        [0, 0, 0, 0, ..., y, x, 0]
        [0, 0, 0, 0, ..., z, y, x]
        [0, 0, 0, 0, ..., 0, z, y]

    .. versionadded:: 1.5.0

    References
    ----------
    .. [1] "Convolution", https://en.wikipedia.org/wiki/Convolution

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import convolution_matrix
    >>> A = convolution_matrix([-1, 4, -2], 5, mode='same')
    >>> A
    array([[ 4, -1,  0,  0,  0],
           [-2,  4, -1,  0,  0],
           [ 0, -2,  4, -1,  0],
           [ 0,  0, -2,  4, -1],
           [ 0,  0,  0, -2,  4]])

    Compare multiplication by `A` with the use of `numpy.convolve`.

    >>> x = np.array([1, 2, 0, -3, 0.5])
    >>> A @ x
    array([  2. ,   6. ,  -1. , -12.5,   8. ])

    Verify that ``A @ x`` produced the same result as applying the
    convolution function.

    >>> np.convolve([-1, 4, -2], x, mode='same')
    array([  2. ,   6. ,  -1. , -12.5,   8. ])

    For comparison to the case ``mode='same'`` shown above, here are the
    matrices produced by ``mode='full'`` and ``mode='valid'`` for the
    same coefficients and size.

    >>> convolution_matrix([-1, 4, -2], 5, mode='full')
    array([[-1,  0,  0,  0,  0],
           [ 4, -1,  0,  0,  0],
           [-2,  4, -1,  0,  0],
           [ 0, -2,  4, -1,  0],
           [ 0,  0, -2,  4, -1],
           [ 0,  0,  0, -2,  4],
           [ 0,  0,  0,  0, -2]])

    >>> convolution_matrix([-1, 4, -2], 5, mode='valid')
    array([[-2,  4, -1,  0,  0],
           [ 0, -2,  4, -1,  0],
           [ 0,  0, -2,  4, -1]])
    """
    if n <= 0:
        raise ValueError('n must be a positive integer.')

    a = np.asarray(a)
    if a.ndim != 1:
        raise ValueError('convolution_matrix expects a one-dimensional '
                         'array as input')
    if a.size == 0:
        raise ValueError('len(a) must be at least 1.')

    if mode not in ('full', 'valid', 'same'):
        raise ValueError(
            "'mode' argument must be one of ('full', 'valid', 'same')")

    # create zero padded versions of the array
    az = np.pad(a, (0, n-1), 'constant')
    raz = np.pad(a[::-1], (0, n-1), 'constant')

    if mode == 'same':
        trim = min(n, len(a)) - 1
        tb = trim//2
        te = trim - tb
        col0 = az[tb:len(az)-te]
        row0 = raz[-n-tb:len(raz)-tb]
    elif mode == 'valid':
        tb = min(n, len(a)) - 1
        te = tb
        col0 = az[tb:len(az)-te]
        row0 = raz[-n-tb:len(raz)-tb]
    else:  # 'full'
        col0 = az
        row0 = raz[-n:]
    return toeplitz(col0, row0)


class sakurai:
    """
    Construct a Sakurai matrix in various formats and its eigenvalues.

    Constructs the "Sakurai" matrix motivated by reference [1]_.
    The matrix is square real symmetric positive definite and banded
    with analytically known eigenvalues.
    The matrix gets ill-conditioned with its size growing.
    See the notes below for details.

    Parameters
    ----------
    n : int
        The size of the matrix.

    Returns
    -------
    sakurai_obj: custom object
        The object containing the output
    sakurai_obj.array : (n, n) ndarray, float
        The Sakurai matrix in the ndarray format
    sakurai_obj.sparse : (n, n) sparse matrix, float
        The Sakurai matrix in a DIAgonal sparse format
    sakurai_obj.banded : (3, n) ndarray, float
        The Sakurai matrix in the format for banded symmetric matrices,
        i.e., 3 upper diagonals with the main diagonal at the bottom
    sakurai_obj.callable : callable object
        The handle to a function that multiplies the Sakurai matrix
        `s` of the shape `n`-by-`n` on the right by an input matrix `x`
        of the shape `n`-by-`k` to output ``s @ x`` without constructing `s`
    sakurai_obj.eigenvalues : (n, ) ndarray, float
        Eigenvalues of the Sakurai matrix ordered ascending

    Notes
    -----
    Reference [1]_ introduces a generalized eigenproblem for the matrix pair
    `A` and `B` where `A` is the identity so we turn it into an eigenproblem
    just for the matrix `B` that this function outputs in various formats
    together with its eigenvalues.
    
    .. versionadded:: 1.11.2

    References
    ----------
    .. [1] T. Sakurai, H. Tadano, Y. Inadomi, and U. Nagashima,
       "A moment-based method for large-scale generalized
       eigenvalue problems",
       Appl. Num. Anal. Comp. Math. Vol. 1 No. 2 (2004).

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import sakurai, eig_banded
    >>> n = 6
    >>> sak = sakurai(n)
    >>> sak.array
    array([[ 5., -4.,  1.,  0.,  0.,  0.],
           [-4.,  6., -4.,  1.,  0.,  0.],
           [ 1., -4.,  6., -4.,  1.,  0.],
           [ 0.,  1., -4.,  6., -4.,  1.],
           [ 0.,  0.,  1., -4.,  6., -4.],
           [ 0.,  0.,  0.,  1., -4.,  5.]])
    >>> sak.banded
    array([[ 1.,  1.,  1.,  1.,  1.,  1.],
           [-4., -4., -4., -4., -4., -4.],
           [ 5.,  6.,  6.,  6.,  6.,  5.]])
    >>> sak.sparse
    <6x6 sparse matrix of type '<class 'numpy.float64'>'
        with 24 stored elements (5 diagonals) in DIAgonal format>
    >>> np.array_equal(sak.sparse.A, sak.array)
    True
    >>> np.array_equal(sak.callable(np.eye(n)), sak.array)
    True
    >>> sak.eigenvalues
    array([0.03922866, 0.56703972, 2.41789479, 5.97822974,
           10.54287655, 14.45473055])

    The banded form can be used in scipy functions for banded matrices, e.g.,

    >>> e = eig_banded(sak.banded, eigvals_only=True)
    >>> np.allclose(sak.eigenvalues, e, atol= n * n * n * np.finfo(float).eps)
    True

    """
    def __init__(self, n) -> None:
        from scipy.sparse import spdiags
        self.n = n
        d0 = np.r_[5, 6 * np.ones(n - 2), 5]
        d1 = -4 * np.ones(n)
        d2 = np.ones(n)
        s = spdiags([d2, d1, d0, d1, d2], [-2, -1, 0, 1, 2], n, n)
        self.sparse = s
        self.array = s.toarray()
        self.banded = np.array([d2, d1, d0])

        k = np.arange(1, n+1)
        e = np.sort(16. * np.power(np.cos(0.5 * k * np.pi / (n + 1)), 4))
        self.eigenvalues = e


    def callable(self, x):
        n = self.n
        assert n == x.shape[0]
        x = x.reshape(n, -1)
        sx = np.zeros_like(x)
        sx[0,:] = 5*x[0,:] - 4*x[1,:] + x[2,:]
        sx[-1,:] = 5*x[-1,:] - 4*x[-2,:] + x[-3,:]
        sx[1:-1,:] = (6*x[1:-1,:] - 4*(x[:-2,:] + x[2:,:])
                      + np.pad(x[:-3,:], ((1,0),(0,0)))
                      + np.pad(x[3:,:], ((0,1),(0,0))))
        return sx


class mikota_pair:
    """
    Construct the Mikota pair of matrices in various formats and
    eigenvalues of the generalized eigenproblem with them.

    The Mikota pair of matrices [1, 2]_ models a vibration problem
    of a linear mass-spring system with the ends attached where
    the stiffness of the springs and the masses increase along
    the system length such that vibration frequencies are subsequent
    integers 1, 2, ..., `n` where `n` is the number of the masses. Thus,
    eigenvalues of the generalized eigenvalue problem for
    the matrix pair `K` and `M` where `K` is he system stiffness matrix
    and `M` is the system mass matrix are the squares of the integers,
    i.e., 1, 4, 9, ..., ``n * n``.

    The stiffness matrix `K` is square real tri-diagonal symmetric
    positive definite. The mass matrix `M` is diagonal with diagonal
    entries 1, 1/2, 1/3, ...., ``1/n``. Both matrices get
    ill-conditioned with `n` growing.

    Parameters
    ----------
    n : int
        The size of the matrices of the Mikota pair.

    Returns
    -------
    mikota_obj: custom object
        The object containing the output
    mikota_obj.Karray : (n, n) ndarray, float
        The stiffness matrix in the ndarray format
    mikota_obj.Ksparse : (n, n) sparse matrix, float
        The stiffness matrix in a DIAgonal sparse format
    mikota_obj.Kbanded : (2, n) ndarray, int32
        The stiffness matrix in the format for banded symmetric matrices,
        i.e., 2 upper diagonals with the main diagonal at the bottom
    mikota_obj.Kcallable : callable object
        The handle to a function that multiplies the stiffness matrix
        `K` of the shape `n`-by-`n` on the right by an input matrix `x`
        of the shape `n`-by-`k` to output ``K @ x`` without constructing `K`
    mikota_obj.Marray : (n, n) ndarray, float
        The mass matrix in the ndarray format
    mikota_obj.Msparse : (n, n) sparse matrix, float
        The mass matrix in a DIAgonal sparse format
    mikota_obj.Mbanded : (1, n) ndarray, float
        The main diagonal of the mass matrix
    mikota_obj.Mcallable : callable object
        The handle to a function that multiplies the mass matrix
        `M` of the shape `n`-by-`n` on the right by an input matrix `x`
        of the shape `n`-by-`k` to output ``M @ x`` without constructing `M`
    mikota_obj.eigenvalues : (n, ) ndarray, float
        Eigenvalues of the Mikota matrix pair: 1, 4, 9, ..., ``n * n``
    
    .. versionadded:: 1.11.2

    References
    ----------
    .. [1] J. Mikota, "Frequency tuning of chain structure multibody oscillators
       to place the natural frequencies at omega1 and N-1 integer multiples
       omega2,..., omegaN", Z. Angew. Math. Mech. 81 (2001), S2, S201-S202.
       Appl. Num. Anal. Comp. Math. Vol. 1 No. 2 (2004).
    .. [2] Peter C. Muller and Metin Gurgoze,
       "Natural frequencies of a multi-degree-of-freedom vibration system",
       Proc. Appl. Math. Mech. 6, 319-320 (2006).
       http://dx.doi.org/10.1002/pamm.200610141.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import mikota_pair
    >>> n = 6
    >>> mik = mikota_pair(n)
    >>> mik.Karray
    array([[11., -5.,  0.,  0.,  0.,  0.],
           [-5.,  9., -4.,  0.,  0.,  0.],
           [ 0., -4.,  7., -3.,  0.,  0.],
           [ 0.,  0., -3.,  5., -2.,  0.],
           [ 0.,  0.,  0., -2.,  3., -1.],
           [ 0.,  0.,  0.,  0., -1.,  1.]])
    >>> mik.Kbanded
    array([[ 0, -5, -4, -3, -2, -1],
           [11,  9,  7,  5,  3,  1]])
    >>> mik.Mbanded
    array([1.        , 0.5       , 0.33333333, 0.25      , 0.2       ,
        0.16666667])
    >>> mik.Ksparse
    <6x6 sparse matrix of type '<class 'numpy.float64'>'
        with 16 stored elements (3 diagonals) in DIAgonal format>
    >>> mik.Msparse
    <6x6 sparse matrix of type '<class 'numpy.float64'>'
        with 6 stored elements (1 diagonals) in DIAgonal format>
    >>> np.array_equal(mik.Ksparse.A, mik.Karray)
    True
    >>> np.array_equal(mik.Msparse.A, mik.Marray)
    True
    >>> np.array_equal(mik.Kcallable(np.eye(n)), mik.Karray)
    True
    >>> np.array_equal(mik.Mcallable(np.eye(n)), mik.Marray)
    True
    >>> mik.eigenvalues
    array([ 1.,  4.,  9., 16., 25., 36.])  

    """
    def __init__(self, n) -> None:
        from scipy.sparse import diags
        self.n = n
        aranp1 = np.arange(1, n + 1)
        aranp1_inv = 1. / aranp1
        self.Mbanded = aranp1_inv
        M = diags([aranp1_inv], [0], shape=(n, n))
        self.Msparse = M
        self.Marray = M.toarray()

        y = - np.arange(n - 1, 0, -1)
        z = np.arange(2 * n - 1, 0, -2)
        K = diags([y, z, y], [-1, 0, 1], shape=(n, n))
        self.Ksparse = K
        self.Karray = K.toarray()
        self.Kbanded = np.array([np.pad(y, (1, 0), 'constant'), z])

        self.eigenvalues = aranp1 * aranp1.astype(float)


    def Mcallable(self, x):
        n = self.n
        assert n == x.shape[0]
        aranp1_inv = 1. / np.arange(1, n + 1)
        # linearoperator requires 2D array
        if len(x.shape) == 1:
            x = x.reshape(-1, 1)
        return aranp1_inv[:, np.newaxis] * x


    def Kcallable(self, x):
        n = self.n
        assert n == x.shape[0]
        x = x.reshape(n, -1)
        kx = np.zeros_like(x)
        y = - np.arange(n - 1, 0, -1)
        z = np.arange(2 * n - 1, 0, -2)
        kx[0, :] = z[0] * x[0, :] + y[0] * x[1, :]
        kx[-1, :] = y[-1] * x[-2, :] + z[-1] * x[-1, :]
        kx[1: -1, :] = (y[:-1, None] * x[: -2,:]
                        + z[1: -1, None] * x[1: -1, :]
                        + y[1:, None] * x[2:, :])
        return kx


class laplacian:
    """
    Construct a sparse matrix and a callable associated with the Laplacian
    on a uniform rectangular grid in `N`>1 dimensions given by tuple `shape`
    and output its eigenvalues.

    The Laplacian matrix `L` is square real symmetric negative definite.

    Parameters
    ----------
    shape : tuple of int
        The shape of the grid for the `N` dimensional Laplacian of length `N`
        where the shape `n` of the Laplacian matrix `L` is ``np.prod(shape)``
    bc: ``'D'`` or ``'N'`` or ``'P'``
        The type of the boundary conditions on the boundaries of the grid

    Returns
    -------
    laplacian_obj: custom object
        The object containing the output
    laplacian_obj.sparse : (n, n) sparse matrix, float
        The Laplacian matrix in a Compressed Sparse Row sparse format
    laplacian_obj.callable : callable object
        The handle to a function that multiplies the Laplacian matrix
        `L` of the shape `n`-by-`n` on the right by an input matrix `x`
        of the shape `n`-by-`k` to output ``L @ x`` without constructing `L`
    laplacian_obj.eigenvalues : (n, ) ndarray, float
        Eigenvalues of the laplacian matrix (pair) ordered ascending

    .. versionadded:: 1.11.2

    Notes
    -----
    Compared to the MATLAB/Octave implementation [1] of 1-, 2-, and 3-D
    Laplacian, this code allows the arbitrary N-D case and the matrix-free
    callable option, but is currently limited to pure Dirichlet, Neumann or
    periodic boundary conditions only and outputs just the eigenvalues,
    no eigenvectors.

    The Laplacian matrix of a graph [2] of a rectangular grid corresponds to
    the negative Laplacian with the Neumann conditions, i.e., ``bc = `N```.


    References
    ----------
    .. [1] https://github.com/lobpcg/blopex/blob/master/blopex_tools/matlab/laplacian/laplacian.m
    .. [2] https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.laplacian.html

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import laplacian
    >>> from scipy.linalg import eigh
    >>> shape = (2, 3)
    >>> n = np.prod(shape)
    >>> bc = 'D'
    >>> lap = laplacian(shape, bc)
    >>> L = lap.sparse
    >>> L
    <6x6 sparse matrix of type '<class 'numpy.float64'>'
        with 20 stored elements in Compressed Sparse Row format>
    >>> L.A
    array([[-4.,  1.,  0.,  1.,  0.,  0.],
           [ 1., -4.,  1.,  0.,  1.,  0.],
           [ 0.,  1., -4.,  0.,  0.,  1.],
           [ 1.,  0.,  0., -4.,  1.,  0.],
           [ 0.,  1.,  0.,  1., -4.,  1.],
           [ 0.,  0.,  1.,  0.,  1., -4.]])
    >>> np.array_equal(lap.callable(np.eye(n)), L.A)
    True
    >>> lap.eigenvalues
    array([-6.41421356, -5.        , -4.41421356, -3.58578644, -3.        ,
           -1.58578644])
    >>> eigvals = np.sort(eigh(L.A, eigvals_only=True))
    >>> np.allclose(lap.eigenvalues, eigvals)
    True
    >>> bc = 'P'
    >>> lap = laplacian(shape, bc)
    >>> L = lap.sparse
    >>> L
    <6x6 sparse matrix of type '<class 'numpy.float64'>'
        with 24 stored elements in Compressed Sparse Row format>
    >>> L.A
    array([[-4.,  1.,  1.,  2.,  0.,  0.],
           [ 1., -4.,  1.,  0.,  2.,  0.],
           [ 1.,  1., -4.,  0.,  0.,  2.],
           [ 2.,  0.,  0., -4.,  1.,  1.],
           [ 0.,  2.,  0.,  1., -4.,  1.],
           [ 0.,  0.,  2.,  1.,  1., -4.]])
    >>> np.array_equal(lap.callable(np.eye(n)), L.A)
    True
    >>> lap.eigenvalues
    array([-7., -7., -4., -3., -3.,  0.])
    >>> eigvals = np.sort(eigh(L.A, eigvals_only=True))
    >>> np.allclose(lap.eigenvalues, eigvals)
    True
    >>> bc = 'N'
    >>> lap = laplacian(shape, bc)
    >>> L = lap.sparse
    >>> L
    <6x6 sparse matrix of type '<class 'numpy.float64'>'
        with 20 stored elements in Compressed Sparse Row format>
    >>> L.A
    array([[-2.,  1.,  0.,  1.,  0.,  0.],
           [ 1., -3.,  1.,  0.,  1.,  0.],
           [ 0.,  1., -2.,  0.,  0.,  1.],
           [ 1.,  0.,  0., -2.,  1.,  0.],
           [ 0.,  1.,  0.,  1., -3.,  1.],
           [ 0.,  0.,  1.,  0.,  1., -2.]])
    ## >>> np.array_equal(lap.callable(np.eye(n)), L.A)
    ## True
    >>> lap.eigenvalues
    array([-5., -3., -3., -2., -1.,  0.])
    >>> eigvals = np.sort(eigh(L.A, eigvals_only=True))
    >>> np.allclose(lap.eigenvalues, eigvals)
    True

    """


    def __init__(self, shape, bc) -> None:
        def lsm(_shape, _bc):
            from scipy.sparse import spdiags, kron, eye, dia_array

            N = len(_shape)
            L = 0
            for i in range(N):
                shape_i = _shape[i]
                diagonals = np.array(
                    [
                        np.ones(shape_i),
                        -2 * np.ones(shape_i),
                        np.ones(shape_i),
                    ]
                )
                if _bc == "N":
                    diagonals[1, 0] = -1.0
                    diagonals[1, -1] = -1.0
                elif bc == "P" or "D":
                    pass
                else:
                    raise NotImplementedError(
                        f"Boundary condition {bc} is not implemented."
                    )

                offsets = np.array([-1, 0, 1])
                L_i = spdiags(diagonals, offsets, shape_i, shape_i)
                if _bc == "P":
                    t = dia_array((shape_i, shape_i))
                    t.setdiag([1.], k=-shape_i + 1)
                    t.setdiag([1.], k=shape_i - 1)
                    L_i += t
                else:
                    pass
                for j in range(i):
                    L_i = kron(eye(_shape[j]), L_i)
                for j in range(i + 1, N):
                    L_i = kron(L_i, eye(_shape[j]))
                L += L_i
            return L

        if len(shape) == 1:
            raise NotImplementedError(
                f"1D Laplacian given by shape={shape} is not implemented."
            )
        self.shape = shape
        self.bc = bc
        self.sparse = lsm(shape, bc)

        # All eigenvalues of the discrete Laplacian operator for
        # an ``N``-dimensional  regular grid of shape `shape` with ``h=1``.
        # https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors_of_the_second_derivative

        indices = np.indices(shape)
        l = np.zeros(shape)

        for j, n in zip(indices, shape):
            if bc == "D":  # pure Dirichlet boundary conditions
                l += -4 * np.sin(np.pi * (j + 1) / (2 * (n + 1))) ** 2
            elif bc == "N":  # pure Neumann boundary conditions
                l += -4 * np.sin(np.pi * j / (2 * n)) ** 2
            elif bc == "P":  # periodic boundary conditions
                l += -4 * np.sin(np.pi * np.floor((j + 1) / 2) / n) ** 2
            else:
                raise NotImplementedError(
                    f"Boundary condition {bc} is not implemented."
                )
        self.eigenvalues = np.sort(l.ravel())

    def callable(self, x):
        shape = self.shape
        bc = self.bc
        if bc == "N":
            raise NotImplementedError(
                f"Neumann boundary condition {bc} is not yet implemented for callable."
            )
        else:
            pass
        assert np.prod(shape) == x.shape[0]
        N = len(shape)
        X = x.reshape(shape + (-1,))
        Y = -2 * N * X
        for i in range(N):
            Y += np.roll(X, 1, axis=i)
            Y += np.roll(X, -1, axis=i)
            if bc == "D":
                Y[
                    (slice(None),) * i + (0,) + (slice(None),) * (N - i - 1)
                ] -= np.roll(X, 1, axis=i)[
                    (slice(None),) * i + (0,) + (slice(None),) * (N - i - 1)
                ]
                Y[
                    (slice(None),) * i + (-1,) + (slice(None),) * (N - i - 1)
                ] -= np.roll(X, -1, axis=i)[
                    (slice(None),) * i + (-1,) + (slice(None),) * (N - i - 1)
                ]
            elif bc == "P":
                pass
        return Y.reshape(-1, X.shape[-1])

