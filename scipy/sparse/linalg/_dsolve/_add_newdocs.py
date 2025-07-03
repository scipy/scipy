from numpy.lib import add_newdoc

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU',
    """
    LU factorization of a sparse matrix.

    Factorization is represented as::

        Pr @ A @ Pc = L @ U

    To construct these `SuperLU` objects, call the `splu` and `spilu`
    functions.

    Attributes
    ----------
    shape
    nnz
    perm_c
    perm_r
    L
    U

    Methods
    -------
    solve

    Notes
    -----

    .. versionadded:: 0.14.0

    Examples
    --------
    The LU decomposition can be used to solve matrix equations. Consider:

    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import splu
    >>> A = csc_array([[1,2,0,4], [1,0,0,1], [1,0,2,1], [2,2,1,0.]])

    This can be solved for a given right-hand side:

    >>> lu = splu(A)
    >>> b = np.array([1, 2, 3, 4])
    >>> x = lu.solve(b)
    >>> A.dot(x)
    array([ 1.,  2.,  3.,  4.])

    The ``lu`` object also contains an explicit representation of the
    decomposition. The permutations are represented as mappings of
    indices:

    >>> lu.perm_r
    array([2, 1, 3, 0], dtype=int32)  # may vary
    >>> lu.perm_c
    array([0, 1, 3, 2], dtype=int32)  # may vary

    The L and U factors are sparse matrices in CSC format:

    >>> lu.L.toarray()
    array([[ 1. ,  0. ,  0. ,  0. ],  # may vary
           [ 0.5,  1. ,  0. ,  0. ],
           [ 0.5, -1. ,  1. ,  0. ],
           [ 0.5,  1. ,  0. ,  1. ]])
    >>> lu.U.toarray()
    array([[ 2. ,  2. ,  0. ,  1. ],  # may vary
           [ 0. , -1. ,  1. , -0.5],
           [ 0. ,  0. ,  5. , -1. ],
           [ 0. ,  0. ,  0. ,  2. ]])

    The permutation matrices can be constructed:

    >>> Pr = csc_array((np.ones(4), (lu.perm_r, np.arange(4))))
    >>> Pc = csc_array((np.ones(4), (np.arange(4), lu.perm_c)))

    We can reassemble the original matrix:

    >>> (Pr.T @ (lu.L @ lu.U) @ Pc.T).toarray()
    array([[ 1.,  2.,  0.,  4.],
           [ 1.,  0.,  0.,  1.],
           [ 1.,  0.,  2.,  1.],
           [ 2.,  2.,  1.,  0.]])
    """)

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('solve',
    """
    solve(rhs[, trans])

    Solves linear system of equations with one or several right-hand sides.

    Parameters
    ----------
    rhs : ndarray, shape (n,) or (n, k)
        Right hand side(s) of equation
    trans : {'N', 'T', 'H'}, optional
        Type of system to solve::

            'N':   A   @ x == rhs  (default)
            'T':   A^T @ x == rhs
            'H':   A^H @ x == rhs

        i.e., normal, transposed, and hermitian conjugate.

    Returns
    -------
    x : ndarray, shape ``rhs.shape``
        Solution vector(s)
    """))

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('L',
    """
    Lower triangular factor with unit diagonal as a
    `scipy.sparse.csc_array`.

    .. versionadded:: 0.14.0
    """))

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('U',
    """
    Upper triangular factor as a `scipy.sparse.csc_array`.

    .. versionadded:: 0.14.0
    """))

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('shape',
    """
    Shape of the original matrix as a tuple of ints.
    """))

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('nnz',
    """
    Number of nonzero elements in the matrix.
    """))

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('perm_c',
    """
    Permutation Pc represented as an array of indices.

    See the `SuperLU` docstring for details.
    """))

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('perm_r',
    """
    Permutation Pr represented as an array of indices.

    See the `SuperLU` docstring for details.
    """))

add_newdoc('scipy.sparse.linalg._dsolve._superlu', 'SuperLU', ('normest_inv',
    """Compute an estimate of the norm of the inverse of a sparse matrix.

    Parameters
    ----------
    ord : {1, inf}, optional
        Order of the norm. If ``None``, defaults to the 1-norm. inf means
        numpy's ``inf`` object.

    Returns
    -------
    est : float
        An estimate of the norm of the matrix inverse.

    See Also
    --------
    norm : Compute the norm of a sparse matrix.
    numpy.linalg.norm : Compute the norm of a dense matrix.
    cond1est : Compute an estimate for the reciprocal of
        the condition number of a sparse matrix.

    Notes
    -----
    This method runs the appropriate SuperLU ``gscon`` [0]_ procedure for the
    data type.

    .. versionadded:: 1.17.0

    References
    ----------
    .. [0] https://portal.nersc.gov/project/sparse/superlu/superlu_code_html/dgscon_8c.html

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import splu
    >>> A = csc_array([[1., 0., 0.], [5., 8., 2.], [0., -1., 0.]], dtype=float)
    >>> A.toarray()
    array([[ 1.,  0.,  0.],
           [ 5.,  8.,  2.],
           [ 0., -1.,  0.]])
    >>> splu(A).normest_inv(ord=1)
    np.float64(5.0)
    >>> np.linalg.norm(np.linalg.inv(A.toarray()), ord=1)
    np.float64(5.0)
    """))
