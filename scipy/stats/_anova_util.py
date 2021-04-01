from itertools import zip_longest
import numpy as np
from scipy.linalg import svd


def _fmt(value):
    return np.array2string(np.array(value))


def _common_ljust(strings):
    width = max(len(t) for t in strings)
    return [t.ljust(width) for t in strings]


def _common_rjust(strings):
    width = max(len(t) for t in strings)
    return [t.rjust(width) for t in strings]


def _common_len_str(label, values):
    s = [label] + [_fmt(t) for t in values]
    return _common_rjust(s)


def _make_table(header, rowlabels, ss, df, ms, f, p):
    txt = [header]
    rowlabels = _common_ljust(['Source'] + rowlabels)
    ss = _common_len_str('SS', ss)
    df = _common_len_str('DF', df)
    ms = _common_len_str('MS', ms)
    f = _common_len_str('F', f)
    p = _common_len_str('p', p)
    fields = zip_longest(rowlabels, ss, df, ms, f, p, fillvalue='')
    for line_fields in fields:
        line = ' '.join(line_fields)
        txt.append(line)
    return '\n'.join(txt)


def _nway_groups(*factors, values, levels=None):
    """
    Parameters
    ----------
    factors : one or more 1-d sequences of values
        The factors (i.e. the independent variables) to be analyzed.
        Generally these should be integers or categorical values.
    values : 1-d sequence
        The values associated with the corresponding factors.
        A real-valued array.

    Examples
    --------
    >>> x = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2]
    >>> y = [5, 5, 5, 7, 7, 7, 7, 7, 5, 5, 7, 7, 7, 5, 5]
    >>> z = [1.4, 2.0, 1.8, 1.7, 1.6, 1.8, 2.1, 2.0, 2.1,
    ...      1.9, 2.4, 2.3, 2.3, 2.9, 2.8]
    >>> levels, groups = _nway_groups(x, y, values=z)

    `levels` is the unique values in each factor.  As is easily
    verified, we find three distinct values in `x` and two in `y`.

    >>> levels
    (array([0, 1, 2]), array([5, 7]))

    `groups` is an object array with shape (3, 2).  Each value in
    the array is a 1-d sequence of values from `z`.

    >>> for i in range(len(groups)):
    ...      for j in range(groups.shape[1]):
    ...          print(i, j, (levels[0][i], levels[1][j]), groups[i, j])
    ...
    0 0 (0, 5) [1.4 2.  1.8]
    0 1 (0, 7) [1.7 1.6 1.8]
    1 0 (1, 5) [2.1 1.9]
    1 1 (1, 7) [2.1 2. ]
    2 0 (2, 5) [2.9 2.8]
    2 1 (2, 7) [2.4 2.3 2.3]
    """
    factors = [np.asarray(a) for a in factors]
    values = np.asarray(values)
    if len(factors) == 0:
        raise TypeError("At least one input factor is required.")
    if not all(len(a) == len(factors[0]) for a in factors[1:]):
        raise ValueError("All input factors must be sequences with the same"
                         " length.")
    if len(values) != len(factors[0]):
        raise ValueError('values must have the same length as each factor.')

    if levels is None:
        # Call np.unique with return_inverse=True on each factor.
        actual_levels, inverses = zip(*[np.unique(f, return_inverse=True)
                                        for f in factors])
        shape = [len(u) for u in actual_levels]
        groups = np.empty(shape, dtype=object)
        inverses = np.array(inverses)
        u, idx = np.unique(inverses, axis=1, return_inverse=True)
        u = u.T
        for i in range(len(u)):
            groups[tuple(u[i])] = values[idx == i]
    else:
        raise NotImplementedError('specifying levels is not implemented yet.')
        # `levels` is not None...
        if len(levels) != len(factors):
            raise ValueError('len(levels) must equal the number of input '
                             'sequences')

        factors = [np.asarray(factor) for factor in factors]
        mask = np.zeros((len(factors), len(factors[0])), dtype=np.bool_)
        inv = np.zeros((len(factors), len(factors[0])), dtype=np.intp)
        actual_levels = []
        for k, (levels_list, arg) in enumerate(zip(levels, factors)):
            if levels_list is None:
                levels_list, inv[k, :] = np.unique(arg, return_inverse=True)
                mask[k, :] = True
            else:
                q = arg == np.asarray(levels_list).reshape(-1, 1)
                mask[k, :] = np.any(q, axis=0)
                qnz = q.T.nonzero()
                inv[k, qnz[0]] = qnz[1]
            actual_levels.append(levels_list)

        mask_all = mask.all(axis=0)
        shape = [len(u) for u in actual_levels]
        count = np.zeros(shape, dtype=int)
        indices = tuple(inv[:, mask_all])
        np.add.at(count, indices, 1)

    return actual_levels, groups


def _unpack_2d_data_grid(data):
    # `data` is expected to be
    # * a numpy array with shape (n, m, r), where the last
    #   dimension holds the replicates, OR
    # * a numpy array with dtype object and shape (n, m), where each
    #   value in the array (e.g. data[i, j]) is a one-dimnensional list
    #   of replicates, OR
    # * a list of list of lists, where the outer two containers make up
    #   a two-dimensional grid, and each inner list (e.g. data[i][j])
    #   is a one-dimensional list of replicates.
    #
    a = []
    b = []
    values = []
    for i, row in enumerate(data):
        for j, cell in enumerate(row):
            for value in cell:
                a.append(i)
                b.append(j)
                values.append(value)
    return a, b, values


def _encode(a):
    """
    Encode the levels found in the array `a` as "sum to zero" contrasts.

    Returns
    -------
    v : ndarray with shape (n, m - 1)
        n is len(a), and m is the number of distinct values found in `a`.
    u : ndarray
        The distinct values (or "levels") found in `a`.
    inv : ndarray
        An array with the same shape as `a` that holds the index into `u`
        of the corresponding values in `a`.  This is the result of pass
        `return_inverse=True` in `np.unique(a, return_inverse=True)``

    Examples
    --------
    Suppose `a` is [10, 20, 40, 40, 10, 10, 10, 30, 30].  There are
    four distinct values (10, 20, 30, 40).  The "sum to zero" contrast
    encoding for four levels is::

        Level    Encoding
          10    -1  -1  -1
          20     1   0   0
          30     0   1   0
          40     0   0   1

    >>> from scipy.stats._anova_util import _encode
    >>> a = [10, 20, 40, 40, 10, 10, 10, 30, 30]
    >>> v, u, inv = _encode(a)

    `u` and `inv` are the result of ``np.unique(a, return_inverse=True)`

    >>> u
    array([10, 20, 30, 40])
    >>> inv
    array([0, 1, 3, 3, 0, 0, 0, 2, 2])

    `v` is the array with shape ``(len(a), len(u) - 1))``, with each
    row holding the corresponding row from the contrast matrix that
    encodes the levels. I.e. ``a[0]`` is 10, so ``v[0]`` is [-1, -1, -1];
    ``a[2]`` is 40, so ``v[2]`` is [0, 0, 1], etc.

    >>> v
    array([[-1, -1, -1],
           [ 1,  0,  0],
           [ 0,  0,  1],
           [ 0,  0,  1],
           [-1, -1, -1],
           [-1, -1, -1],
           [-1, -1, -1],
           [ 0,  1,  0],
           [ 0,  1,  0]])
    """
    u, i = np.unique(a, return_inverse=True)
    m = len(u) - 1
    vars = np.vstack((np.full(m, fill_value=-1), np.eye(m, dtype=int)))
    return vars[i], u, i


def _encode_nested(a, b):
    """
    Create variable part of the design matrix for factor b nested in
    factor a.

    `b` must be a numpy array.
    """
    n = len(a)
    au, ai = np.unique(a, return_inverse=True)
    nu = len(au)

    masks = [ai == k for k in range(nu)]

    bsplit = [b[mask] for mask in masks]

    bui = [np.unique(b1, return_inverse=True) for b1 in bsplit]
    mb = [len(bu) - 1 for (bu, _) in bui]
    bvars = [np.vstack((np.full(len(bu) - 1, fill_value=-1),
                        np.eye(len(bu) - 1, dtype=int)))
             for (bu, _) in bui]
    result = np.zeros((n, sum(mb)), dtype=int)
    starts = np.cumsum(np.concatenate(([0], mb[:-1])))
    for k in range(nu):
        result[masks[k], starts[k]:starts[k]+mb[k]] = bvars[k][bui[k][1]]

    cell_sizes = [np.bincount(bi1) for (bu1, bi1) in bui]
    return np.array(result), cell_sizes


def _interaction(X1, X2):
    m1 = X1.shape[1]
    m2 = X2.shape[1]
    X = np.zeros((X1.shape[0], m1*m2))
    k = 0
    for i in range(m1):
        for j in range(m2):
            X[:, k] = X1[:, i]*X2[:, j]
            k += 1
    return X


def _proj(x):
    """
    Create the projection matrix (also known as the "hat matrix").

    `x` is the design matrix for a linear model.

    See, for example, https://en.wikipedia.org/wiki/Projection_matrix
    """
    # Equivalent to H = x @ inv(x.T @ x) @ x.T
    U, s, Vh = svd(x, full_matrices=False)
    H = U @ U.T
    return H
