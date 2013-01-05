import sys
from numpy.testing import *
import numpy.linalg as linalg

from scipy.lib.six import print_

def random(size):
    return rand(*size)

class TestSolve(TestCase):

    def bench_random(self):
        basic_solve = linalg.solve
        print_()
        print_('      Solving system of linear equations')
        print_('      ==================================')

        print_('      |    contiguous     |   non-contiguous ')
        print_('----------------------------------------------')
        print_(' size |  scipy  | basic   |  scipy  | basic ')

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print_('%5s' % size, end=' ')
            sys.stdout.flush()

            a = random([size,size])
            # larger diagonal ensures non-singularity:
            for i in range(size): a[i,i] = 10*(.1+a[i,i])
            b = random([size])

            print_('| %6.2f ' % measure('solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            print_('| %6.2f ' % measure('basic_solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            a = a[-1::-1,-1::-1] # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

            print_('| %6.2f ' % measure('solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            print_('| %6.2f ' % measure('basic_solve(a,b)',repeat), end=' ')
            sys.stdout.flush()

            print_('   (secs for %s calls)' % (repeat))

class TestInv(TestCase):

    def bench_random(self):
        basic_inv = linalg.inv
        print_()
        print_('           Finding matrix inverse')
        print_('      ==================================')
        print_('      |    contiguous     |   non-contiguous ')
        print_('----------------------------------------------')
        print_(' size |  scipy  | basic   |  scipy  | basic')

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print_('%5s' % size, end=' ')
            sys.stdout.flush()

            a = random([size,size])
            # large diagonal ensures non-singularity:
            for i in range(size): a[i,i] = 10*(.1+a[i,i])

            print_('| %6.2f ' % measure('inv(a)',repeat), end=' ')
            sys.stdout.flush()

            print_('| %6.2f ' % measure('basic_inv(a)',repeat), end=' ')
            sys.stdout.flush()

            a = a[-1::-1,-1::-1] # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

            print_('| %6.2f ' % measure('inv(a)',repeat), end=' ')
            sys.stdout.flush()

            print_('| %6.2f ' % measure('basic_inv(a)',repeat), end=' ')
            sys.stdout.flush()

            print_('   (secs for %s calls)' % (repeat))


class TestDet(TestCase):

    def bench_random(self):
        basic_det = linalg.det
        print_()
        print_('           Finding matrix determinant')
        print_('      ==================================')
        print_('      |    contiguous     |   non-contiguous ')
        print_('----------------------------------------------')
        print_(' size |  scipy  | basic   |  scipy  | basic ')

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print_('%5s' % size, end=' ')
            sys.stdout.flush()

            a = random([size,size])

            print_('| %6.2f ' % measure('det(a)',repeat), end=' ')
            sys.stdout.flush()

            print_('| %6.2f ' % measure('basic_det(a)',repeat), end=' ')
            sys.stdout.flush()

            a = a[-1::-1,-1::-1] # turn into a non-contiguous array
            assert_(not a.flags['CONTIGUOUS'])

            print_('| %6.2f ' % measure('det(a)',repeat), end=' ')
            sys.stdout.flush()

            print_('| %6.2f ' % measure('basic_det(a)',repeat), end=' ')
            sys.stdout.flush()

            print_('   (secs for %s calls)' % (repeat))


if __name__ == "__main__":
    run_module_suite()
