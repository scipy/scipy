from __future__ import division, absolute_import, print_function

import math

import numpy.linalg as nl

import numpy as np
from numpy.testing import assert_
from numpy.random import rand

try:
    import scipy.linalg as sl
except ImportError:
    pass

from .common import Benchmark


def random(size):
    return rand(*size)


class Bench(Benchmark):
    params = [
        [20, 100, 500, 1000],
        ['contig', 'nocont'],
        ['numpy', 'scipy']
    ]
    param_names = ['size', 'contiguous', 'module']

    def setup(self, size, contig, module):
        a = random([size,size])
        # larger diagonal ensures non-singularity:
        for i in range(size):
            a[i,i] = 10*(.1+a[i,i])
        b = random([size])

        if contig != 'contig':
            a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

        self.a = a
        self.b = b

    def time_solve(self, size, contig, module):
        if module == 'numpy':
            nl.solve(self.a, self.b)
        else:
            sl.solve(self.a, self.b)

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

    def time_svd(self, size, contig, module):
        if module == 'numpy':
            nl.svd(self.a)
        else:
            sl.svd(self.a)


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
        np.random.seed(1234)
        n = math.ceil(2./3. * size)
        k = math.ceil(1./2. * size)
        m = size

        if dtype is np.complex128:
            A = ((10 * np.random.rand(m,k) - 5) +
                 1j*(10 * np.random.rand(m,k) - 5))
            temp = ((10 * np.random.rand(k,n) - 5) +
                    1j*(10 * np.random.rand(k,n) - 5))
            b = ((10 * np.random.rand(m,1) - 5) +
                 1j*(10 * np.random.rand(m,1) - 5))
        else:
            A = (10 * np.random.rand(m,k) - 5)
            temp = 10 * np.random.rand(k,n) - 5
            b = 10 * np.random.rand(m,1) - 5

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
