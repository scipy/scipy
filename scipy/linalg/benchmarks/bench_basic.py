from __future__ import division, print_function, absolute_import

import sys
from numpy.testing import measure, rand, assert_, TestCase
import numpy.linalg as nl
import scipy.linalg as sl


def random(size):
    return rand(*size)


class TestSolve(TestCase):

    def bench_random(self):
        numpy_solve = nl.solve
        scipy_solve = sl.solve
        print()
        print('      Solving system of linear equations')
        print('      ==================================')

        print('      |    contiguous     |   non-contiguous ')
        print('----------------------------------------------')
        print(' size |  scipy  | numpy   |  scipy  | numpy ')

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print('%5s' % size, end=' ')
            sys.stdout.flush()

            a = random([size,size])
            # larger diagonal ensures non-singularity:
            for i in range(size):
                a[i,i] = 10*(.1+a[i,i])
            b = random([size])

            print('| %6.2f ' % measure('scipy_solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            print('| %6.2f ' % measure('numpy_solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

            print('| %6.2f ' % measure('scipy_solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            print('| %6.2f ' % measure('numpy_solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            print('   (secs for %s calls)' % (repeat))


class TestInv(TestCase):

    def bench_random(self):
        numpy_inv = nl.inv
        scipy_inv = sl.inv
        print()
        print('           Finding matrix inverse')
        print('      ==================================')
        print('      |    contiguous     |   non-contiguous ')
        print('----------------------------------------------')
        print(' size |  scipy  | numpy   |  scipy  | numpy')

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print('%5s' % size, end=' ')
            sys.stdout.flush()

            a = random([size,size])
            # large diagonal ensures non-singularity:
            for i in range(size):
                a[i,i] = 10*(.1+a[i,i])

            print('| %6.2f ' % measure('scipy_inv(a)',repeat), end=' ')
            sys.stdout.flush()

            print('| %6.2f ' % measure('numpy_inv(a)',repeat), end=' ')
            sys.stdout.flush()

            a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

            print('| %6.2f ' % measure('scipy_inv(a)',repeat), end=' ')
            sys.stdout.flush()

            print('| %6.2f ' % measure('numpy_inv(a)',repeat), end=' ')
            sys.stdout.flush()

            print('   (secs for %s calls)' % (repeat))


class TestDet(TestCase):

    def bench_random(self):
        numpy_det = nl.det
        scipy_det = sl.det
        print()
        print('           Finding matrix determinant')
        print('      ==================================')
        print('      |    contiguous     |   non-contiguous ')
        print('----------------------------------------------')
        print(' size |  scipy  | numpy   |  scipy  | numpy ')

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print('%5s' % size, end=' ')
            sys.stdout.flush()

            a = random([size,size])

            print('| %6.2f ' % measure('scipy_det(a)',repeat), end=' ')
            sys.stdout.flush()

            print('| %6.2f ' % measure('numpy_det(a)',repeat), end=' ')
            sys.stdout.flush()

            a = a[-1::-1,-1::-1]  # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

            print('| %6.2f ' % measure('scipy_det(a)',repeat), end=' ')
            sys.stdout.flush()

            print('| %6.2f ' % measure('numpy_det(a)',repeat), end=' ')
            sys.stdout.flush()

            print('   (secs for %s calls)' % (repeat))


if __name__ == "__main__":
    run_module_suite()
