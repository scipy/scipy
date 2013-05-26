""" Benchmark functions for linalg.decomp module

"""

from __future__ import division, print_function, absolute_import

import sys

from numpy import linalg as nl
from scipy import linalg as sl
from numpy.testing import measure, rand, assert_


def random(size):
    return rand(*size)


def bench_eigvals():
    numpy_eigvals = nl.eigvals
    scipy_eigvals = sl.eigvals
    print()
    print('           Finding matrix eigenvalues')
    print('      ==================================')
    print('      |    contiguous     |   non-contiguous ')
    print('----------------------------------------------')
    print(' size |  scipy  | numpy   |  scipy  | numpy ')

    for size,repeat in [(20,150),(100,7),(200,2)]:
        repeat *= 1
        print('%5s' % size, end=' ')
        sys.stdout.flush()

        a = random([size,size])

        print('| %6.2f ' % measure('scipy_eigvals(a)',repeat), end=' ')
        sys.stdout.flush()

        print('| %6.2f ' % measure('numpy_eigvals(a)',repeat), end=' ')
        sys.stdout.flush()

        a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
        assert_(not a.flags['CONTIGUOUS'])

        print('| %6.2f ' % measure('scipy_eigvals(a)',repeat), end=' ')
        sys.stdout.flush()

        print('| %6.2f ' % measure('numpy_eigvals(a)',repeat), end=' ')
        sys.stdout.flush()

        print('   (secs for %s calls)' % (repeat))


def bench_svd():
    numpy_svd = nl.svd
    scipy_svd = sl.svd
    print()
    print('           Finding the SVD decomposition')
    print('      ==================================')
    print('      |    contiguous     |   non-contiguous ')
    print('----------------------------------------------')
    print(' size |  scipy  | numpy   |  scipy  | numpy ')

    for size,repeat in [(20,150),(100,7),(200,2)]:
        repeat *= 1
        print('%5s' % size, end=' ')
        sys.stdout.flush()

        a = random([size,size])

        print('| %6.2f ' % measure('scipy_svd(a)',repeat), end=' ')
        sys.stdout.flush()

        print('| %6.2f ' % measure('numpy_svd(a)',repeat), end=' ')
        sys.stdout.flush()

        a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
        assert_(not a.flags['CONTIGUOUS'])

        print('| %6.2f ' % measure('scipy_svd(a)',repeat), end=' ')
        sys.stdout.flush()

        print('| %6.2f ' % measure('numpy_svd(a)',repeat), end=' ')
        sys.stdout.flush()

        print('   (secs for %s calls)' % (repeat))
