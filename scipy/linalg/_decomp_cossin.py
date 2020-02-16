from collections.abc import Iterable

import numpy as np
from numpy.core.numerictypes import genericTypeRank
from numpy.linalg import LinAlgError

from scipy._lib._util import _asarray_validated
from scipy.linalg import block_diag
from .lapack import _compute_lwork, get_lapack_funcs

__all__ = ['cossin']


def cossin(X, partitioning=None, compute_u=True, compute_vh=True):
    """
    Compute the cosine-sine (CS) decomposition of an `(m, m)` matrix.

    The matrix X is a unitary or an orthogonal matrix, partitioned to X11, X12,
    X21, X22 subblocks and the decomposition satisfies the following equality::

                                [  I11 0  0  |  0   0  0  ]
                                [  0   C  0  |  0  -S  0  ]
            [X11|X12]   [U1|  ] [  0   0  0  |  0   0 -I12] [V1|  ]**H
        X = [-------] = [-----] [-------------------------] [-----]
            [X21|X22]   [  |U2] [  0  0  0   |  I22 0  0  ] [  |V2]
                                [  0  S  0   |  0   C  0  ]
                                [  0  0  I21 |  0   0  0  ]

    `U1`, `U2`, `V1`, `V2` are orthogonal matrices.
    `C` and `S` are nonnegative diagonal matrices satisfying `C^2 + S^2 = I`.

    Dimensions:

    * `X11: (p, q)`
    * `U1: (p, p)`
    * `U2: (m-p, m-p)`
    * `V1: (q, q)`
    * `V2: (m-q, m-q)`
    * `C: (r,r)`
    * `S: (r,r)`

    where `r = min(p,m-p,q,m-q)`.

    The rank of the identity submatrices:

    * `I11: min(p, q) - r`
    * `I12: min(p, m - q) - r`
    * `I21: min(m - p, q) - r`
    * `I22: min(m - p, m - q) - r`.

    Partitioning can be given two ways:

    #. directly by listing the subblocks
    #. providing the whole matrix and providing the dimensions of the top left subblock in `partitioning`

    Parameters
    ----------
    X : (m,m) array_like or iterable
        complex unitary or real orthogonal matrix to be decomposed, or iterable
        of subblocks `X11`, `X12`, `X21`, `X22`, when ``partitioning=None``
    partitioning : (int, int), optional
        shape of the upper left block `X11`, if omitted, `X` is assumed to be given
        in iterable form
    compute_u : bool, optional
        if `False`, `u` won't be computed
    compute_vh : bool, optional
        if `False`, `vh` won't be computed

    Returns
    -------
    u : (m, m) ndarray
        When `compute_u=True`, contains the block diagonal unitary or orthogonal
        matrix consisting of the blocks U1 (p-by-p) and U2 (m-p)-by-(m-p)
    cs : (m, m) ndarray
        The cosine-sine factor with the structure described above
    vh : (m, m) ndarray
        When `compute_u=True`, contains the block diagonal unitary or orthogonal
        matrix consisting of the blocks V1 (q-by-q) and V2 (m-q)-by-(m-q)

    Examples
    --------
    >>> from scipy.linalg import cossin
    >>> from scipy.stats import unitary_group
    >>> x = unitary_group.rvs(4)
    >>> u, cs, vdh = cossin(x, (2,2))
    >>> np.allclose(x, u @ cs @ vdh)
    True

    References
    ----------
    .. [1] : Brian D. Sutton. Computing the complete CS decomposition. Numer.
           Algorithms, 50(1):33-65, 2009.

    """
    if partitioning:
        X = _asarray_validated(X, check_finite=True)
        if not np.equal(*X.shape):
            raise ValueError("Cosine Sine decomposition only supports square"
                             " matrices, got {}".format(X.shape))
        p, q = partitioning
        return _cossin(x11=X[:p, :q], x12=X[:p, q:], x21=X[p:, :q],
                       x22=X[p:, q:],
                       compute_u=compute_u, compute_vh=compute_vh)

    if not isinstance(X, Iterable):
        raise ValueError("When top_left_shape=None, X must be Iterable")

    X = list(X)
    if len(X) != 4:
        raise ValueError("When partitioning=None, "
                         "exactly four submatrices should be in X, got {}"
                         .format(len(X)))

    x11, x12, x21, x22 = X

    p, q = np.shape(np.atleast_2d(x11))
    mmp, mmq = np.shape(np.atleast_2d(x22))

    p2, mmq2 = np.shape(np.atleast_2d(x12))
    mmp2, q2 = np.shape(np.atleast_2d(x21))

    if p != p2 or q != q2 or mmp != mmp2 or mmq != mmq2 or mmq + q != mmp + p:
        shapes = [np.shape(np.atleast_2d(subx)) for subx in X]
        raise ValueError("invalid submatrix dimensions {} to build a square "
                         "matrix, maybe you forgot to use \"partitioning\" keyword?".format(shapes))
    return _cossin(x11, x12, x21, x22, compute_u, compute_vh)


def _cossin(x11, x12, x21, x22, compute_u, compute_vh):
    p, q = np.shape(x11)
    m = np.shape(x22)[0] + p

    x11, x12, x21, x22 = _cast_all_to_largest(x11, x12, x21, x22)

    iscomp = np.iscomplexobj(x11)
    func_name = "uncsd" if iscomp else "orcsd"

    csd, csd_lwork = get_lapack_funcs([func_name, func_name + "_lwork"],
                                      [x11])
    lwork = _compute_lwork(csd_lwork, m=m, p=p, q=q)

    lwork_args = (dict(lwork=lwork[0], lrwork=lwork[1]) if iscomp else
                  dict(lwork=lwork))

    csd_args = dict(x11=x11,
                    x12=x12,
                    x21=x21,
                    x22=x22,
                    compute_u1=compute_u,
                    compute_u2=compute_u,
                    compute_v1t=compute_vh,
                    compute_v2t=compute_vh,
                    trans=False,
                    signs=False)

    csd_args.update(lwork_args)

    x11o, x12o, x21o, x22o, theta, u1, u2, v1h, v2h, info = csd(
        **csd_args)

    method_name = csd.typecode + func_name
    if info < 0:
        raise ValueError('illegal value in argument {} of internal {}'
                         .format(-info, method_name))
    if info > 0:
        raise LinAlgError("{} did not converge: {}".format(method_name, info))

    U = block_diag(u1, u2)
    VDH = block_diag(v1h, v2h)
    CS = _construct_cs(theta, m, p, q, x11.dtype.type)

    return U, CS, VDH


def _cast_all_to_largest(*args):
    arg_with_largest_type = max(args, key=(
        lambda a: genericTypeRank.index(a.dtype.name)))
    return (arg_with_largest_type.dtype.type(x) for x in args)


def _construct_cs(theta, m, p, q, dtype):
    r = min(p, q, m - p, m - q)
    n11 = min(p, q) - r
    n12 = min(p, m - q) - r
    n21 = min(m - p, q) - r
    n22 = min(m - p, m - q) - r
    CS = np.zeros((m, m), dtype=dtype)
    one = dtype(1.)
    for i in range(n11):
        CS[i, i] = one
    for i in range(n22):
        CS[p + i, q + i] = one
    for i in range(n12):
        CS[i + n11 + r, i + n11 + r + n21 + n22 + r] = -one
    for i in range(n21):
        CS[p + n22 + r + i, n11 + r + i] = one

    for i in range(r):
        CS[i + n11, i + n11] = np.cos(theta[i])
        CS[p + n22 + i, i + r + n21 + n22] = np.cos(theta[i])

        CS[i + n11, i + n11 + n21 + n22 + r] = -np.sin(theta[i])
        CS[p + n22 + i, i + n11] = np.sin(theta[i])
    return CS
