import math

import numpy.linalg as nl

import numpy as np
from numpy.testing import assert_
from numpy.random import rand

from .common import Benchmark, safe_import

with safe_import():
    import scipy.linalg as sl


def random(size):
    return rand(*size)


class Bench(Benchmark):
    params = [
        [20, 100, 500, 1000],
        ['contig', 'nocont', 'fcontig'],
        ['numpy', 'scipy']
    ]
    param_names = ['size', 'contiguous', 'module']

    def __init__(self):
        # likely not useful to benchmark svd for large sizes
        self.time_svd.__func__.params = [[20, 100, 500]] + self.params[1:]

    def setup(self, size, contig, module):
        if module == 'numpy' and size >= 200:
            # skip: slow, and not useful to benchmark numpy
            raise NotImplementedError()

        a = random([size, size])
        # larger diagonal ensures non-singularity:
        for i in range(size):
            a[i, i] = 10*(.1+a[i, i])
        b = random([size])

        if contig == 'nocont':
            a = a[-1::-1, -1::-1]  # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])
        elif contig == 'fcontig':
            a = np.asfortranarray(a)

        self.a = a
        self.b = b
        self.a_pos = a @ a.T + size * np.eye(size)

    def time_solve(self, size, contig, module):
        if module == 'numpy':
            nl.solve(self.a, self.b)
        else:
            sl.solve(self.a, self.b)

    def time_solve_triangular(self, size, contig, module):
        # treats self.a as a lower-triangular matrix by ignoring the strictly
        # upper-triangular part
        if module == 'numpy':
            pass
        else:
            sl.solve_triangular(self.a, self.b, lower=True)

    def time_inv(self, size, contig, module):
        if module == 'numpy':
            nl.inv(self.a)
        else:
            sl.inv(self.a)

    def time_det(self, size, contig, module):
        if module == 'numpy':
            nl.det(self.a)
        else:
            sl.det(self.a)

    def time_eigvals(self, size, contig, module):
        if module == 'numpy':
            nl.eigvals(self.a)
        else:
            sl.eigvals(self.a)

    def time_geneig(self, size, contig, module):
        if module == 'numpy':
            pass
        else:
            sl.eig(self.a, self.a, check_finite=True)

    def time_svd(self, size, contig, module):
        if module == 'numpy':
            nl.svd(self.a)
        else:
            sl.svd(self.a)

    def time_cholesky(self, size, contig, module):
        if module == 'numpy':
            nl.cholesky(self.a_pos)
        else:
            sl.cholesky(self.a_pos)



class BatchedSolveBench(Benchmark):
    params = [
        [(100, 10, 10), (100, 20, 20), (100, 100)],
        ["gen", "pos", "sym", "diagonal", "tridiagonal", "banded"],
        ["scipy/detect", "scipy/assume", "numpy"]
    ]
    param_names = ["shape", "structure" ,"module"]

    def setup(self, shape, structure, module):
        a = random(shape)
        # larger diagonal ensures non-singularity:
        for i in range(shape[-1]):
            a[..., i, i] = 10*(.1+a[..., i, i])

        if structure == "pos":
            self.a = a @ a.mT
        elif structure == "sym":
            self.a = a + a.mT
        elif structure == "diagonal":
            self.a = np.zeros_like(a)
            for i in range(shape[-1]):
                self.a[..., i, i] = a[..., i, i]
        elif structure == "tridiagonal":
            self.a = np.zeros_like(a)
            for i in range(shape[-1]):
                self.a[..., i, i] = a[..., i, i]
            for i in range(shape[-1]-1):
                self.a[..., i+1, i] = a[..., i+1, i]
            for i in range(shape[-1]-1):
                self.a[..., i, i+1] = a[..., i, i+1]
        elif structure == "banded":
            self.a = np.zeros_like(a)
            self.a += np.triu(np.tril(a, k=5), k=-5)
        else:
            self.a = a

        self.b = random([a.shape[-1]])

        self.kwd = {}
        if module.split("/")[-1] == "assume":
            self.kwd = {"assume_a": structure}

    def time_solve(self, shape, structure, module):
        if module == 'numpy':
            nl.solve(self.a, self.b)
        else:
            sl.solve(self.a, self.b, check_finite=False, **self.kwd)


class BatchedSVDBench(Benchmark):
    params = [
        [(10, 10, 10, 2), (100, 10, 10), (100, 20, 20), (100, 100, 100)],
        ["scipy", "numpy"]
    ]
    param_names = ['shape',  'module']

    def setup(self, shape, module):
        self.a = random(shape)

    def time_svd(self, shape, module):
        if module == 'numpy':
            nl.svd(self.a)
        else:
            sl.svd(self.a)


class BatchedPinvBench(Benchmark):
    params = [
        [(10, 10, 10, 2), (100, 10, 10), (100, 20, 20), (100, 100, 100)],
        ["scipy", "numpy"]
    ]
    param_names = ['shape',  'module']

    def setup(self, shape, module):
        self.a = random(shape)

    def time_pinv(self, shape, module):
        if module == 'numpy':
            nl.pinv(self.a)
        else:
            sl.pinv(self.a)


class BatchedLstsqBench(Benchmark):
    params = [
        [(10, 10, 50, 2), (100, 20, 5), (100, 10, 10), (100, 5, 20), (100, 2, 50)],
    ]
    param_names = ['shape']

    def setup(self, shape):
         self.a = random(shape)
         self.b = random((shape[-2],))

    def time_lstsq(self, shape):
        sl.lstsq(self.a, self.b, check_finite=False)


class BatchedEigBench(Benchmark):
    params = [
        [(10, 10, 3, 3), (100, 10, 10), (100, 20, 20), (100, 100, 100)],
        ["scipy", "numpy"]
    ]
    param_names = ['shape',  'module']

    def setup(self, shape, module):
        self.a = random(shape)

    def time_eig(self, shape, module):
        if module == 'numpy':
            nl.eig(self.a)
        else:
            sl.eig(self.a)


class BatchedCholeskyBench(Benchmark):
    params = [
        [(100, 3, 3), (100, 10, 10), (100, 20, 20), (100, 100, 100), (100, 100)],
        ["scipy", "numpy"]
    ]
    param_names  = ["shape", "module"]

    def setup(self, shape, module):
        x = random(shape[:-1] + (1000,))
        self.a = x @ np.swapaxes(x, axis1=-2, axis2=-1)

    def time_cholesky(self, shape, module):
        if module == "numpy":
            nl.cholesky(self.a)
        else:
            sl.cholesky(self.a)


class BatchedQRBench(Benchmark):
    params = [
        [(100, 10, 10), (100, 20, 20), (100, 100)],
        ["full/complete", "economic/reduced", "r/r", "raw/raw"],
        ["scipy", "numpy"]
    ]
    param_names = ["shape", "mode", "module"]

    def setup(self, shape, mode, module):
        self.a = random(shape)

        if module == "scipy":
            self.kwd = {"mode": mode.split("/")[0]}
        elif module == "numpy":
            self.kwd = {"mode": mode.split("/")[1]}

    def time_solve(self, shape, mode, module):
        if module == "numpy":
            nl.qr(self.a, **self.kwd)
        else:
            sl.qr(self.a, **self.kwd)


class Norm(Benchmark):
    params = [
        [(20, 20), (100, 100), (1000, 1000), (20, 1000), (1000, 20)],
        ['contig', 'nocont'],
        ['numpy', 'scipy']
    ]
    param_names = ['shape', 'contiguous', 'module']

    def setup(self, shape, contig, module):
        a = np.random.randn(*shape)
        if contig != 'contig':
            a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])
        self.a = a

    def time_1_norm(self, size, contig, module):
        if module == 'numpy':
            nl.norm(self.a, ord=1)
        else:
            sl.norm(self.a, ord=1)

    def time_inf_norm(self, size, contig, module):
        if module == 'numpy':
            nl.norm(self.a, ord=np.inf)
        else:
            sl.norm(self.a, ord=np.inf)

    def time_frobenius_norm(self, size, contig, module):
        if module == 'numpy':
            nl.norm(self.a)
        else:
            sl.norm(self.a)


class Lstsq(Benchmark):
    """
    Test the speed of four least-squares solvers on not full rank matrices.
    Also check the difference in the solutions.

    The matrix has the size ``(m, 2/3*m)``; the rank is ``1/2 * m``.
    Matrix values are random in the range (-5, 5), the same is for the right
    hand side.  The complex matrix is the sum of real and imaginary matrices.
    """

    param_names = ['dtype', 'size', 'driver']
    params = [
        [np.float64, np.complex128],
        [10, 100, 1000],
        ['gelss', 'gelsy', 'gelsd', 'numpy'],
    ]

    def setup(self, dtype, size, lapack_driver):
        if lapack_driver == 'numpy' and size >= 200:
            # skip: slow, and not useful to benchmark numpy
            raise NotImplementedError()

        rng = np.random.default_rng(1234)
        n = math.ceil(2./3. * size)
        k = math.ceil(1./2. * size)
        m = size

        if dtype is np.complex128:
            A = ((10 * rng.random((m,k)) - 5) +
                 1j*(10 * rng.random((m,k)) - 5))
            temp = ((10 * rng.random((k,n)) - 5) +
                    1j*(10 * rng.random((k,n)) - 5))
            b = ((10 * rng.random((m,1)) - 5) +
                 1j*(10 * rng.random((m,1)) - 5))
        else:
            A = (10 * rng.random((m,k)) - 5)
            temp = 10 * rng.random((k,n)) - 5
            b = 10 * rng.random((m,1)) - 5

        self.A = A.dot(temp)
        self.b = b

    def time_lstsq(self, dtype, size, lapack_driver):
        if lapack_driver == 'numpy':
            np.linalg.lstsq(self.A, self.b,
                            rcond=np.finfo(self.A.dtype).eps * 100)
        else:
            sl.lstsq(self.A, self.b, cond=None, overwrite_a=False,
                     overwrite_b=False, check_finite=False,
                     lapack_driver=lapack_driver)

    # Retain old benchmark results (remove this if changing the benchmark)
    time_lstsq.version = (
        "15ee0be14a0a597c7d1c9a3dab2c39e15c8ac623484410ffefa406bf6b596ebe"
    )


class SpecialMatrices(Benchmark):
    param_names = ['size']
    params = [[4, 128]]

    def setup(self, size):
        self.x = np.arange(1, size + 1).astype(float)
        self.small_blocks = [np.ones([2, 2])] * (size//2)
        self.big_blocks = [np.ones([size//2, size//2]),
                           np.ones([size//2, size//2])]

    def time_block_diag_small(self, size):
        sl.block_diag(*self.small_blocks)

    def time_block_diag_big(self, size):
        sl.block_diag(*self.big_blocks)

    def time_circulant(self, size):
        sl.circulant(self.x)

    def time_companion(self, size):
        sl.companion(self.x)

    def time_dft(self, size):
        sl.dft(size)

    def time_hadamard(self, size):
        sl.hadamard(size)

    def time_hankel(self, size):
        sl.hankel(self.x)

    def time_helmert(self, size):
        sl.helmert(size)

    def time_hilbert(self, size):
        sl.hilbert(size)

    def time_invhilbert(self, size):
        sl.invhilbert(size)

    def time_leslie(self, size):
        sl.leslie(self.x, self.x[1:])

    def time_pascal(self, size):
        sl.pascal(size)

    def time_invpascal(self, size):
        sl.invpascal(size)

    def time_toeplitz(self, size):
        sl.toeplitz(self.x)


class GetFuncs(Benchmark):
    def setup(self):
        self.x = np.eye(1)

    def time_get_blas_funcs(self):
        sl.blas.get_blas_funcs('gemm', dtype=float)

    def time_get_blas_funcs_2(self):
        sl.blas.get_blas_funcs(('gemm', 'axpy'), (self.x, self.x))

    def time_small_cholesky(self):
        sl.cholesky(self.x)
