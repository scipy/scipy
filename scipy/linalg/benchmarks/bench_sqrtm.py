""" Benchmark linalg.sqrtm for various blocksizes.

"""

from __future__ import division, print_function, absolute_import

import time

import numpy as np
from numpy.testing import assert_allclose
import scipy.linalg


def bench_sqrtm():
    np.random.seed(1234)
    print()
    print('                       Matrix Square Root')
    print('==============================================================')
    print('      shape      |   block   |        dtype       |   time   ')
    print('                 |    size   |                    | (seconds)')
    print('--------------------------------------------------------------')
    fmt = ' %15s |   %6d  | %18s | %6.2f '
    for n in (64, 256):
        for dtype in (np.float64, np.complex128):

            # Sample a random matrix.
            # See Figs. 2 and 3 of Deadman and Higham and Ralha 2012
            # "A Recursive Blocked Schur Algorithm
            # for Computing the # Matrix Square Root"
            A = np.random.rand(n, n)
            if dtype == np.complex128:
                A = A + 1j*np.random.rand(n, n)

            # blocksize that corresponds to the block-free implementation
            blocksize = n
            tm = time.clock()
            B_1, info = scipy.linalg.sqrtm(A, disp=False, blocksize=blocksize)
            nseconds = time.clock() - tm
            print(fmt % (A.shape, blocksize, A.dtype, nseconds))

            # interesting blocksize
            blocksize = 32
            tm = time.clock()
            B_32, info = scipy.linalg.sqrtm(A, disp=False, blocksize=blocksize)
            nseconds = time.clock() - tm
            print(fmt % (A.shape, blocksize, A.dtype, nseconds))

            # bigger block size
            blocksize = 64
            tm = time.clock()
            B_64, info = scipy.linalg.sqrtm(A, disp=False, blocksize=blocksize)
            nseconds = time.clock() - tm
            print(fmt % (A.shape, blocksize, A.dtype, nseconds))

            # Check that the results are the same for all block sizes.
            assert_allclose(B_1, B_32)
            assert_allclose(B_1, B_64)
