import numpy as np
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import kron, eye, dia_array

__all__ = ['LaplacianNd']


class LaplacianNd(LinearOperator):
    """
    The grid Laplacian in ``N`` dimensions and its eigenvalues/eigenvectors.

    Construct Laplacian on a uniform rectangular grid in `N` dimensions
    and output its eigenvalues and eigenvectors.
    The Laplacian ``L`` is square, negative definite, real symmetric array
    with signed integer entries and zeros otherwise.

    Parameters
    ----------
    grid_shape : tuple
        A tuple of integers of length ``N`` (corresponding to the dimension of
        the Lapacian), where each entry gives the size of that dimension. The
        Laplacian matrix is square of the size ``np.prod(grid_shape)``.
    boundary_conditions : {'neumann', 'dirichlet', 'periodic'}, optional
        The type of the boundary conditions on the boundaries of the grid.
        Valid values are ``'dirichlet'`` or ``'neumann'``(default) or
        ``'periodic'``.
    dtype : dtype
        Numerical type of the array. Default is ``np.int8``.

    Methods
    -------
    toarray()
        Construct a dense array from Laplacian data
    tosparse()
        Construct a sparse array from Laplacian data
    eigenvalues(m=None)
        Construct a 1D array of `m` largest (smallest in absolute value)
        eigenvalues of the Laplacian matrix in ascending order.
    eigenvectors(m=None):
        Construct the array with columns made of `m` eigenvectors (``float``)
        of the ``Nd`` Laplacian corresponding to the `m` ordered eigenvalues.

    .. versionadded:: 1.12.0

    Notes
    -----
    Compared to the MATLAB/Octave implementation [1] of 1-, 2-, and 3-D
    Laplacian, this code allows the arbitrary N-D case and the matrix-free
    callable option, but is currently limited to pure Dirichlet, Neumann or
    Periodic boundary conditions only.

    The Laplacian matrix of a graph (`scipy.sparse.csgraph.laplacian`) of a
    rectangular grid corresponds to the negative Laplacian with the Neumann
    conditions, i.e., ``boundary_conditions = 'neumann'``.

    All eigenvalues and eigenvectors of the discrete Laplacian operator for
    an ``N``-dimensional  regular grid of shape `grid_shape` with the grid
    step size ``h=1`` are analytically known [2].

    References
    ----------
    .. [1] https://github.com/lobpcg/blopex/blob/master/blopex_\
tools/matlab/laplacian/laplacian.m
    .. [2] "Eigenvalues and eigenvectors of the second derivative", Wikipedia
           https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors_\
of_the_second_derivative

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse.linalg import LaplacianNd
    >>> from scipy.sparse import diags, csgraph
    >>> from scipy.linalg import eigvalsh

    The one-dimensional Laplacian demonstrated below for pure Neumann boundary
    conditions on a regular grid with ``n=6`` grid points is exactly the
    negative graph Laplacian for the undirected linear graph with ``n``
    vertices using the sparse adjacency matrix ``G`` represented by the
    famous tri-diagonal matrix:

    >>> n = 6
    >>> G = diags(np.ones(n - 1), 1, format='csr')
    >>> Lf = csgraph.laplacian(G, symmetrized=True, form='function')
    >>> grid_shape = (n, )
    >>> lap = LaplacianNd(grid_shape, boundary_conditions='neumann')
    >>> np.array_equal(lap.matmat(np.eye(n)), -Lf(np.eye(n)))
    True

    Since all matrix entries of the Laplacian are integers, ``'int8'`` is
    the default dtype for storing matrix representations.

    >>> lap.tosparse()
    <6x6 sparse matrix of type '<class 'numpy.int8'>'
        with 16 stored elements (3 diagonals) in DIAgonal format>
    >>> lap.toarray()
    array([[-1,  1,  0,  0,  0,  0],
           [ 1, -2,  1,  0,  0,  0],
           [ 0,  1, -2,  1,  0,  0],
           [ 0,  0,  1, -2,  1,  0],
           [ 0,  0,  0,  1, -2,  1],
           [ 0,  0,  0,  0,  1, -1]], dtype=int8)
    >>> np.array_equal(lap.matmat(np.eye(n)), lap.toarray())
    True
    >>> np.array_equal(lap.tosparse().toarray(), lap.toarray())
    True

    Any number of extreme eigenvalues and/or eigenvectors can be computed.
    
    >>> lap = LaplacianNd(grid_shape, boundary_conditions='periodic')
    >>> lap.eigenvalues()
    array([-4., -3., -3., -1., -1.,  0.])
    >>> lap.eigenvalues()[-2:]
    array([-1.,  0.])
    >>> lap.eigenvalues(2)
    array([-1.,  0.])
    >>> lap.eigenvectors(1)
    array([[0.40824829],
           [0.40824829],
           [0.40824829],
           [0.40824829],
           [0.40824829],
           [0.40824829]])
    >>> lap.eigenvectors(2)
    array([[ 0.5       ,  0.40824829],
           [ 0.        ,  0.40824829],
           [-0.5       ,  0.40824829],
           [-0.5       ,  0.40824829],
           [ 0.        ,  0.40824829],
           [ 0.5       ,  0.40824829]])
    >>> lap.eigenvectors()
    array([[ 0.40824829,  0.28867513,  0.28867513,  0.5       ,  0.5       ,
             0.40824829],
           [-0.40824829, -0.57735027, -0.57735027,  0.        ,  0.        ,
             0.40824829],
           [ 0.40824829,  0.28867513,  0.28867513, -0.5       , -0.5       ,
             0.40824829],
           [-0.40824829,  0.28867513,  0.28867513, -0.5       , -0.5       ,
             0.40824829],
           [ 0.40824829, -0.57735027, -0.57735027,  0.        ,  0.        ,
             0.40824829],
           [-0.40824829,  0.28867513,  0.28867513,  0.5       ,  0.5       ,
             0.40824829]])

    The two-dimensional Laplacian is illustrated on a regular grid with
    ``grid_shape = (2, 3)`` points in each dimension.

    >>> grid_shape = (2, 3)
    >>> n = np.prod(grid_shape)

    Numeration of grid points is as follows:

    >>> np.arange(n).reshape(grid_shape + (-1,))
    array([[[0],
            [1],
            [2]],
    <BLANKLINE>
           [[3],
            [4],
            [5]]])

    Each of the boundary conditions ``'dirichlet'``, ``'periodic'``, and
    ``'neumann'`` is illustrated separately; with ``'dirichlet'``

    >>> lap = LaplacianNd(grid_shape, boundary_conditions='dirichlet')
    >>> lap.tosparse()
    <6x6 sparse array of type '<class 'numpy.int8'>'
        with 20 stored elements in Compressed Sparse Row format>
    >>> lap.toarray()
    array([[-4,  1,  0,  1,  0,  0],
           [ 1, -4,  1,  0,  1,  0],
           [ 0,  1, -4,  0,  0,  1],
           [ 1,  0,  0, -4,  1,  0],
           [ 0,  1,  0,  1, -4,  1],
           [ 0,  0,  1,  0,  1, -4]], dtype=int8)
    >>> np.array_equal(lap.matmat(np.eye(n)), lap.toarray())
    True
    >>> np.array_equal(lap.tosparse().toarray(), lap.toarray())
    True
    >>> lap.eigenvalues()
    array([-6.41421356, -5.        , -4.41421356, -3.58578644, -3.        ,
           -1.58578644])
    >>> eigvals = eigvalsh(lap.toarray().astype(np.float64))
    >>> np.allclose(lap.eigenvalues(), eigvals)
    True
    >>> np.allclose(lap.toarray() @ lap.eigenvectors(),
    ...             lap.eigenvectors() @ np.diag(lap.eigenvalues()))
    True

    with ``'periodic'``

    >>> lap = LaplacianNd(grid_shape, boundary_conditions='periodic')
    >>> lap.tosparse()
    <6x6 sparse array of type '<class 'numpy.int8'>'
        with 24 stored elements in Compressed Sparse Row format>
    >>> lap.toarray()
        array([[-4,  1,  1,  2,  0,  0],
               [ 1, -4,  1,  0,  2,  0],
               [ 1,  1, -4,  0,  0,  2],
               [ 2,  0,  0, -4,  1,  1],
               [ 0,  2,  0,  1, -4,  1],
               [ 0,  0,  2,  1,  1, -4]], dtype=int8)
    >>> np.array_equal(lap.matmat(np.eye(n)), lap.toarray())
    True
    >>> np.array_equal(lap.tosparse().toarray(), lap.toarray())
    True
    >>> lap.eigenvalues()
    array([-7., -7., -4., -3., -3.,  0.])
    >>> eigvals = eigvalsh(lap.toarray().astype(np.float64))
    >>> np.allclose(lap.eigenvalues(), eigvals)
    True
    >>> np.allclose(lap.toarray() @ lap.eigenvectors(),
    ...             lap.eigenvectors() @ np.diag(lap.eigenvalues()))
    True

    and with ``'neumann'``

    >>> lap = LaplacianNd(grid_shape, boundary_conditions='neumann')
    >>> lap.tosparse()
    <6x6 sparse array of type '<class 'numpy.int8'>'
        with 20 stored elements in Compressed Sparse Row format>
    >>> lap.toarray()
    array([[-2,  1,  0,  1,  0,  0],
           [ 1, -3,  1,  0,  1,  0],
           [ 0,  1, -2,  0,  0,  1],
           [ 1,  0,  0, -2,  1,  0],
           [ 0,  1,  0,  1, -3,  1],
           [ 0,  0,  1,  0,  1, -2]])
    >>> np.array_equal(lap.matmat(np.eye(n)), lap.toarray())
    True
    >>> np.array_equal(lap.tosparse().toarray(), lap.toarray())
    True
    >>> lap.eigenvalues()
    array([-5., -3., -3., -2., -1.,  0.])
    >>> eigvals = eigvalsh(lap.toarray().astype(np.float64))
    >>> np.allclose(lap.eigenvalues(), eigvals)
    True
    >>> np.allclose(lap.toarray() @ lap.eigenvectors(),
    ...             lap.eigenvectors() @ np.diag(lap.eigenvalues()))
    True

    """

    def __init__(self, grid_shape, *,
                 boundary_conditions='neumann',
                 dtype=np.int8):

        if boundary_conditions not in ('dirichlet', 'neumann', 'periodic'):
            raise ValueError(
                f"Unknown value {boundary_conditions!r} is given for "
                "'boundary_conditions' parameter. The valid options are "
                "'dirichlet', 'periodic', and 'neumann' (default)."
            )

        self.grid_shape = grid_shape
        self.boundary_conditions = boundary_conditions
        # LaplacianNd folds all dimensions in `grid_shape` into a single one
        N = np.prod(grid_shape)
        super().__init__(dtype=dtype, shape=(N, N))

    def _eigenvalue_ordering(self, m):
        """Compute `m` largest eigenvalues in each of the ``N`` directions,
        i.e., up to ``m * N`` total, order them and return `m` largest.
        """
        grid_shape = self.grid_shape
        if m is None:
            indices = np.indices(grid_shape)
            Leig = np.zeros(grid_shape)
        else:
            grid_shape_min = min(grid_shape,
                                 tuple(np.ones_like(grid_shape) * m))
            indices = np.indices(grid_shape_min)
            Leig = np.zeros(grid_shape_min)

        for j, n in zip(indices, grid_shape):
            if self.boundary_conditions == 'dirichlet':
                Leig += -4 * np.sin(np.pi * (j + 1) / (2 * (n + 1))) ** 2
            elif self.boundary_conditions == 'neumann':
                Leig += -4 * np.sin(np.pi * j / (2 * n)) ** 2
            else:  # boundary_conditions == 'periodic'
                Leig += -4 * np.sin(np.pi * np.floor((j + 1) / 2) / n) ** 2

        Leig_ravel = Leig.ravel()
        ind = np.argsort(Leig_ravel)
        eigenvalues = Leig_ravel[ind]
        if m is not None:
            eigenvalues = eigenvalues[-m:]
            ind = ind[-m:]

        return eigenvalues, ind

    def eigenvalues(self, m=None):
        """Return the requested number of eigenvalues.
        
        Parameters
        ----------
        m : int, optional
            The positive number of eigenvalues to return. If not provided,
            then all eigenvalues will be returned.
            
        Returns
        -------
        eigenvalues : float array
            The requested `m` or all eigenvalues, in ascending order.
        """
        eigenvalues, _ = self._eigenvalue_ordering(m)
        return eigenvalues

    def _ev1d(self, j, n):
        """Return 1 eigenvector in 1d with index `j`
        and number of grid points `n` where ``j < n``. 
        """
        if self.boundary_conditions == 'dirichlet':
            i = np.pi * (np.arange(n) + 1) / (n + 1)
            ev = np.sqrt(2. / (n + 1.)) * np.sin(i * (j + 1))
        elif self.boundary_conditions == 'neumann':
            i = np.pi * (np.arange(n) + 0.5) / n
            ev = np.sqrt((1. if j == 0 else 2.) / n) * np.cos(i * j)
        else:  # boundary_conditions == 'periodic'
            if j == 0:
                ev = np.sqrt(1. / n) * np.ones(n)
            elif j + 1 == n and n % 2 == 0:
                ev = np.sqrt(1. / n) * np.tile([1, -1], n//2)
            else:
                i = 2. * np.pi * (np.arange(n) + 0.5) / n
                ev = np.sqrt(2. / n) * np.cos(i * np.floor((j + 1) / 2))
        # make small values exact zeros correcting round-off errors
        # due to symmetry of eigenvectors the exact 0. is correct 
        ev[np.abs(ev) < np.finfo(np.float64).eps] = 0.
        return ev

    def _one_eve(self, k):
        """Return 1 eigenvector in Nd with multi-index `j`
        as a tensor product of the corresponding 1d eigenvectors. 
        """
        phi = [self._ev1d(j, n) for j, n in zip(k, self.grid_shape)]
        result = phi[0]
        for phi in phi[1:]:
            result = np.tensordot(result, phi, axes=0)
        return np.asarray(result).ravel()

    def eigenvectors(self, m=None):
        """Return the requested number of eigenvectors for ordered eigenvalues.
        
        Parameters
        ----------
        m : int, optional
            The positive number of eigenvectors to return. If not provided,
            then all eigenvectors will be returned.
            
        Returns
        -------
        eigenvectors : float array
            An array with columns made of the requested `m` or all eigenvectors.
            The columns are ordered according to the `m` ordered eigenvalues. 
        """
        _, ind = self._eigenvalue_ordering(m)
        if m is None:
            grid_shape_min = self.grid_shape
        else:
            grid_shape_min = min(self.grid_shape,
                                tuple(np.ones_like(self.grid_shape) * m))

        N_indices = np.unravel_index(ind, grid_shape_min)
        N_indices = [tuple(x) for x in zip(*N_indices)]
        eigenvectors_list = [self._one_eve(k) for k in N_indices]
        return np.column_stack(eigenvectors_list)

    def toarray(self):
        """
        Converts the Laplacian data to a dense array.

        Returns
        -------
        L : ndarray
            The shape is ``(N, N)`` where ``N = np.prod(grid_shape)``.

        """
        grid_shape = self.grid_shape
        n = np.prod(grid_shape)
        L = np.zeros([n, n], dtype=np.int8)
        # Scratch arrays
        L_i = np.empty_like(L)
        Ltemp = np.empty_like(L)

        for ind, dim in enumerate(grid_shape):
            # Start zeroing out L_i
            L_i[:] = 0
            # Allocate the top left corner with the kernel of L_i
            # Einsum returns writable view of arrays
            np.einsum("ii->i", L_i[:dim, :dim])[:] = -2
            np.einsum("ii->i", L_i[: dim - 1, 1:dim])[:] = 1
            np.einsum("ii->i", L_i[1:dim, : dim - 1])[:] = 1

            if self.boundary_conditions == 'neumann':
                L_i[0, 0] = -1
                L_i[dim - 1, dim - 1] = -1
            elif self.boundary_conditions == 'periodic':
                if dim > 1:
                    L_i[0, dim - 1] += 1
                    L_i[dim - 1, 0] += 1
                else:
                    L_i[0, 0] += 1

            # kron is too slow for large matrices hence the next two tricks
            # 1- kron(eye, mat) is block_diag(mat, mat, ...)
            # 2- kron(mat, eye) can be performed by 4d stride trick

            # 1-
            new_dim = dim
            # for block_diag we tile the top left portion on the diagonal
            if ind > 0:
                tiles = np.prod(grid_shape[:ind])
                for j in range(1, tiles):
                    L_i[j*dim:(j+1)*dim, j*dim:(j+1)*dim] = L_i[:dim, :dim]
                    new_dim += dim
            # 2-
            # we need the keep L_i, but reset the array
            Ltemp[:new_dim, :new_dim] = L_i[:new_dim, :new_dim]
            tiles = int(np.prod(grid_shape[ind+1:]))
            # Zero out the top left, the rest is already 0
            L_i[:new_dim, :new_dim] = 0
            idx = [x for x in range(tiles)]
            L_i.reshape(
                (new_dim, tiles,
                 new_dim, tiles)
                )[:, idx, :, idx] = Ltemp[:new_dim, :new_dim]

            L += L_i

        return L.astype(self.dtype)

    def tosparse(self):
        """
        Constructs a sparse array from the Laplacian data. The returned sparse
        array format is dependent on the selected boundary conditions.

        Returns
        -------
        L : scipy.sparse.sparray
            The shape is ``(N, N)`` where ``N = np.prod(grid_shape)``.

        """
        N = len(self.grid_shape)
        p = np.prod(self.grid_shape)
        L = dia_array((p, p), dtype=np.int8)

        for i in range(N):
            dim = self.grid_shape[i]
            data = np.ones([3, dim], dtype=np.int8)
            data[1, :] *= -2

            if self.boundary_conditions == 'neumann':
                data[1, 0] = -1
                data[1, -1] = -1

            L_i = dia_array((data, [-1, 0, 1]), shape=(dim, dim),
                            dtype=np.int8
                            )

            if self.boundary_conditions == 'periodic':
                t = dia_array((dim, dim), dtype=np.int8)
                t.setdiag([1], k=-dim+1)
                t.setdiag([1], k=dim-1)
                L_i += t

            for j in range(i):
                L_i = kron(eye(self.grid_shape[j], dtype=np.int8), L_i)
            for j in range(i + 1, N):
                L_i = kron(L_i, eye(self.grid_shape[j], dtype=np.int8))
            L += L_i
        return L.astype(self.dtype)

    def _matvec(self, x):
        grid_shape = self.grid_shape
        N = len(grid_shape)
        X = x.reshape(grid_shape + (-1,))
        Y = -2 * N * X
        for i in range(N):
            Y += np.roll(X, 1, axis=i)
            Y += np.roll(X, -1, axis=i)
            if self.boundary_conditions in ('neumann', 'dirichlet'):
                Y[(slice(None),)*i + (0,) + (slice(None),)*(N-i-1)
                  ] -= np.roll(X, 1, axis=i)[
                    (slice(None),) * i + (0,) + (slice(None),) * (N-i-1)
                ]
                Y[
                    (slice(None),) * i + (-1,) + (slice(None),) * (N-i-1)
                ] -= np.roll(X, -1, axis=i)[
                    (slice(None),) * i + (-1,) + (slice(None),) * (N-i-1)
                ]

                if self.boundary_conditions == 'neumann':
                    Y[
                        (slice(None),) * i + (0,) + (slice(None),) * (N-i-1)
                    ] += np.roll(X, 0, axis=i)[
                        (slice(None),) * i + (0,) + (slice(None),) * (N-i-1)
                    ]
                    Y[
                        (slice(None),) * i + (-1,) + (slice(None),) * (N-i-1)
                    ] += np.roll(X, 0, axis=i)[
                        (slice(None),) * i + (-1,) + (slice(None),) * (N-i-1)
                    ]

        return Y.reshape(-1, X.shape[-1])

    def _matmat(self, x):
        return self._matvec(x)

    def _adjoint(self):
        return self

    def _transpose(self):
        return self
