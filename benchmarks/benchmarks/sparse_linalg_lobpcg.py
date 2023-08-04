import numpy as np
from .common import Benchmark, safe_import

with safe_import():
    from scipy.linalg import (eigh,
                              cholesky_banded, cho_solve_banded, eig_banded)
    from scipy.sparse.linalg import lobpcg, eigsh, LinearOperator


class Sakurai(LinearOperator):
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
    dtype: numerial numpy type
        The `dtype` of the array, matrix, or banded output.
        The default ``dtype="float"``

    Returns
    -------
    sakurai_lo: `LinearOperator` custom object
        The object containing the output standard for `LinearOperator`:
        sakurai_lo.toarray() : (n, n) ndarray
        sakurai_lo.tosparse() : (n, n) DIAgonal sparse format
        sakurai_lo._matvec and sakurai_lo._matmat: callable objects
            The handle to a function that multiplies the Sakurai matrix
            `s` of the shape `n`-by-`n` on the right by an input matrix `x`
            of the shape `n`-by-`k` to output ``s @ x`` without constructing `s`
    sakurai_lo.tobanded() : (3, n) ndarray
        The Sakurai matrix in the format for banded symmetric matrices,
        i.e., 3 upper diagonals with the main diagonal at the bottom
    sakurai_lo.eigenvalues : (n, ) ndarray, float
        Eigenvalues of the Sakurai matrix ordered ascending

    Notes
    -----
    Reference [1]_ introduces a generalized eigenproblem for the matrix pair
    `A` and `B` where `A` is the identity so we turn it into an eigenproblem
    just for the matrix `B` that this function outputs in various formats
    together with its eigenvalues.
    
    .. versionadded:: 1.12.0

    References
    ----------
    .. [1] T. Sakurai, H. Tadano, Y. Inadomi, and U. Nagashima,
       "A moment-based method for large-scale generalized
       eigenvalue problems",
       Appl. Num. Anal. Comp. Math. Vol. 1 No. 2 (2004).

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.??? import Sakurai
    >>> from scipy.linalg import eig_banded
    >>> n = 6
    >>> sak = Sakurai(n)
    >>> sak.toarray()
    array([[ 5., -4.,  1.,  0.,  0.,  0.],
           [-4.,  6., -4.,  1.,  0.,  0.],
           [ 1., -4.,  6., -4.,  1.,  0.],
           [ 0.,  1., -4.,  6., -4.,  1.],
           [ 0.,  0.,  1., -4.,  6., -4.],
           [ 0.,  0.,  0.,  1., -4.,  5.]])
    >>> sak.tobanded()
    array([[ 1.,  1.,  1.,  1.,  1.,  1.],
           [-4., -4., -4., -4., -4., -4.],
           [ 5.,  6.,  6.,  6.,  6.,  5.]])
    >>> sak.tosparse()
    <6x6 sparse matrix of type '<class 'numpy.float64'>'
        with 24 stored elements (5 diagonals) in DIAgonal format>
    >>> np.array_equal(sak.tosparse().toarray(), sak.toarray())
    True
    >>> np.array_equal(sak(np.eye(n)), sak.tosparse().toarray())
    True
    >>> sak.eigenvalues
    array([0.03922866, 0.56703972, 2.41789479, 5.97822974,
           10.54287655, 14.45473055])

    The entries of the matrix are all integers so can use ``dtype='int'``.

    >>> sak = Sakurai(n, dtype='int')
    >>> sak.tobanded()
    array([[ 1,  1,  1,  1,  1,  1],
           [-4, -4, -4, -4, -4, -4],
           [ 5,  6,  6,  6,  6,  5]])
    >>> sak.tosparse()
    <6x6 sparse matrix of type '<class 'numpy.int32'>'
        with 24 stored elements (5 diagonals) in DIAgonal format>

    The banded form can be used in scipy functions for banded matrices, e.g.,

    >>> e = eig_banded(sak.tobanded(), eigvals_only=True)
    >>> np.allclose(sak.eigenvalues, e, atol= n * n * n * np.finfo(float).eps)
    True

    """
    def __init__(self, n, dtype=np.float64) -> None:
        self.n = n
        self.dtype = dtype
        shape = (n, n)
        super().__init__(dtype, shape)

        k = np.arange(1, n + 1)
        e = np.sort(16. * np.power(np.cos(0.5 * k * np.pi / (n + 1)), 4))
        self.eigenvalues = e


    def tosparse(self):
        from scipy.sparse import spdiags
        d0 = np.r_[5, 6 * np.ones(self.n - 2, dtype=self.dtype), 5]
        d1 = -4 * np.ones(self.n, dtype=self.dtype)
        d2 = np.ones(self.n, dtype=self.dtype)
        return spdiags([d2, d1, d0, d1, d2], [-2, -1, 0, 1, 2], n, n)


    def tobanded(self):
        d0 = np.r_[5, 6 * np.ones(self.n - 2, dtype=self.dtype), 5]
        d1 = -4 * np.ones(self.n, dtype=self.dtype)
        d2 = np.ones(self.n, dtype=self.dtype)
        return np.array([d2, d1, d0])


    def toarray(self):
        d0 = np.r_[5, 6 * np.ones(self.n - 2, dtype=self.dtype), 5]
        d1 = -4 * np.ones(self.n - 1, dtype=self.dtype)
        d2 = np.ones(self.n - 2, dtype=self.dtype)
        a = np.diag(d0)
        a += np.diag(d1, 1) + np.diag(d1, -1)
        a += np.diag(d2, 2) + np.diag(d2, -2)
        return a

    
    def _matvec(self, x):
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


    def _matmat(self, x):
        return self._matvec(x)


    def _adjoint(self):
        return self


    def _transpose(self):
        return self


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


class Bench(Benchmark):
    params = [
        [],
        ['lobpcg', 'eigsh', 'lapack']
    ]
    param_names = ['n', 'solver']

    def __init__(self):
        self.time_mikota.__func__.params = list(self.params)
        self.time_mikota.__func__.params[0] = [128, 256, 512, 1024, 2048]
        self.time_mikota.__func__.setup = self.setup_mikota

        self.time_sakurai.__func__.params = list(self.params)
        self.time_sakurai.__func__.params[0] = [50, 100]
        self.time_sakurai.__func__.setup = self.setup_sakurai

        self.time_sakuraii.__func__.params = list(self.params)
        self.time_sakuraii.__func__.params[0] = [500, 1000, 2000]
        self.time_sakuraii.__func__.setup = self.setup_sakuraii


    def setup_mikota(self, n, solver):
        self.shape = (n, n)
        mikota_pair_obj = mikota_pair(n)
        self.Ac = mikota_pair_obj.Kcallable
        self.Aa = mikota_pair_obj.Karray
        self.Bc = mikota_pair_obj.Mcallable
        self.Ba = mikota_pair_obj.Marray
        self.Ab = mikota_pair_obj.Kbanded
        self.eigenvalues = mikota_pair_obj.eigenvalues

        # if solver == 'eigh' and n >= 512:
        #     # skip: slow, and not useful to benchmark
        #     raise NotImplementedError()


    def setup_sakurai(self, n, solver):
        self.shape = (n, n)
        sakurai_obj = Sakurai(n, dtype='int')
        self.A = sakurai_obj
        self.Aa = sakurai_obj.toarray()
        self.eigenvalues = sakurai_obj.eigenvalues


    def setup_sakuraii(self, n, solver):
        self.shape = (n, n)
        sakurai_obj = Sakurai(n)
        self.A = sakurai_obj.tobanded()
        self.eigenvalues = sakurai_obj.eigenvalues

    
    def time_mikota(self, n, solver):
        def a(x):
            return cho_solve_banded((c, False), x)
        m = 10
        ee = self.eigenvalues[:m]
        tol = m * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            c = cholesky_banded(self.Ab.astype(np.float32))
            el, _ = lobpcg(self.Ac, X, self.Bc, M=a, tol=1e-4,
                                maxiter=40, largest=False)
            accuracy = max(abs(ee - el) / ee)
            assert accuracy < tol
        elif solver == 'eigsh':
            B = LinearOperator((n, n), matvec=self.Bc, matmat=self.Bc, dtype='float64')
            A = LinearOperator((n, n), matvec=self.Ac, matmat=self.Ac, dtype='float64')
            c = cholesky_banded(self.Ab)
            a_l = LinearOperator((n, n), matvec=a, matmat=a, dtype='float64')
            ea, _ = eigsh(B, k=m, M=A, Minv=a_l, which='LA', tol=1e-4, maxiter=50,
                                  v0=rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - np.sort(1./ea)) / ee)
            assert accuracy < tol
        else:
            ed, _ = eigh(self.Aa, self.Ba, subset_by_index=(0, m - 1))
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol


    def time_sakurai(self, n, solver):
        m = 3
        ee = self.eigenvalues[:m]
        tol = 100 * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            el, _ = lobpcg(self.A, X, tol=1e-9, maxiter=5000, largest=False)
            accuracy = max(abs(ee - el) / ee)
            assert accuracy < tol
        elif solver == 'eigsh':
            ea, _ = eigsh(self.A, k=m, which='SA', tol=1e-9, maxiter=15000,
                                   v0=rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - ea) / ee)
            assert accuracy < tol
        else:
            ed, _ = eigh(self.Aa, subset_by_index=(0, m - 1))
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol


    def time_sakuraii(self, n, solver):
        def a(x):
            return cho_solve_banded((c, False), x)
        m = 3
        ee = self.eigenvalues[:m]
        tol = 10 * n * n * n* np.finfo(float).eps
        rng = np.random.default_rng(0)
        X =rng.normal(size=(n, m))
        if solver == 'lobpcg':
            c = cholesky_banded(self.A)
            el, _ = lobpcg(a, X, tol=1e-9, maxiter=8)
            accuracy = max(abs(ee - 1. / el) / ee)
            assert accuracy < tol
        elif solver == 'eigsh':
            c = cholesky_banded(self.A)
            a_l = LinearOperator((n, n), matvec=a, matmat=a, dtype='float64')
            ea, _ = eigsh(a_l, k=m, which='LA', tol=1e-9, maxiter=8,
                                   v0=rng.normal(size=(n, 1)))
            accuracy = max(abs(ee - np.sort(1./ea)) / ee)
            assert accuracy < tol
        else:
            ed, _ = eig_banded(self.A, select='i', select_range=[0, m-1])
            accuracy = max(abs(ee - ed) / ee)
            assert accuracy < tol
