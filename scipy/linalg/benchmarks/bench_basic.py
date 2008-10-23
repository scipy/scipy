import sys
from numpy.testing import *
import numpy.linalg as linalg

def random(size):
    return rand(*size)

class TestSolve(TestCase):

    def bench_random(self):
        basic_solve = linalg.solve
        print
        print '      Solving system of linear equations'
        print '      =================================='

        print '      |    contiguous     |   non-contiguous '
        print '----------------------------------------------'
        print ' size |  scipy  | basic   |  scipy  | basic '

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print '%5s' % size,
            sys.stdout.flush()

            a = random([size,size])
            # larger diagonal ensures non-singularity:
            for i in range(size): a[i,i] = 10*(.1+a[i,i])
            b = random([size])

            print '| %6.2f ' % measure('solve(a,b)',repeat),
            sys.stdout.flush()

            print '| %6.2f ' % measure('basic_solve(a,b)',repeat),
            sys.stdout.flush()

            a = a[-1::-1,-1::-1] # turn into a non-contiguous array
            assert not a.flags['CONTIGUOUS']

            print '| %6.2f ' % measure('solve(a,b)',repeat),
            sys.stdout.flush()

            print '| %6.2f ' % measure('basic_solve(a,b)',repeat),
            sys.stdout.flush()

            print '   (secs for %s calls)' % (repeat)

class TestInv(TestCase):

    def bench_random(self):
        basic_inv = linalg.inv
        print
        print '           Finding matrix inverse'
        print '      =================================='
        print '      |    contiguous     |   non-contiguous '
        print '----------------------------------------------'
        print ' size |  scipy  | basic   |  scipy  | basic'

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print '%5s' % size,
            sys.stdout.flush()

            a = random([size,size])
            # large diagonal ensures non-singularity:
            for i in range(size): a[i,i] = 10*(.1+a[i,i])

            print '| %6.2f ' % measure('inv(a)',repeat),
            sys.stdout.flush()

            print '| %6.2f ' % measure('basic_inv(a)',repeat),
            sys.stdout.flush()

            a = a[-1::-1,-1::-1] # turn into a non-contiguous array
            assert not a.flags['CONTIGUOUS']

            print '| %6.2f ' % measure('inv(a)',repeat),
            sys.stdout.flush()

            print '| %6.2f ' % measure('basic_inv(a)',repeat),
            sys.stdout.flush()

            print '   (secs for %s calls)' % (repeat)


class TestDet(TestCase):

    def bench_random(self):
        basic_det = linalg.det
        print
        print '           Finding matrix determinant'
        print '      =================================='
        print '      |    contiguous     |   non-contiguous '
        print '----------------------------------------------'
        print ' size |  scipy  | basic   |  scipy  | basic '

        for size,repeat in [(20,1000),(100,150),(500,2),(1000,1)][:-1]:
            repeat *= 2
            print '%5s' % size,
            sys.stdout.flush()

            a = random([size,size])

            print '| %6.2f ' % measure('det(a)',repeat),
            sys.stdout.flush()

            print '| %6.2f ' % measure('basic_det(a)',repeat),
            sys.stdout.flush()

            a = a[-1::-1,-1::-1] # turn into a non-contiguous array
            assert not a.flags['CONTIGUOUS']

            print '| %6.2f ' % measure('det(a)',repeat),
            sys.stdout.flush()

            print '| %6.2f ' % measure('basic_det(a)',repeat),
            sys.stdout.flush()

            print '   (secs for %s calls)' % (repeat)


if __name__ == "__main__":
    run_module_suite()
