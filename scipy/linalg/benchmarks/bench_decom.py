""" Benchmark functions for linalg.decomp module

"""
import sys

from numpy import linalg
from numpy.testing import *

def random(size):
    return rand(*size)

def bench_random():
    Numeric_eigvals = linalg.eigvals
    print
    print '           Finding matrix eigenvalues'
    print '      =================================='
    print '      |    contiguous     '#'|   non-contiguous '
    print '----------------------------------------------'
    print ' size |  scipy  '#'| core |  scipy  | core '

    for size,repeat in [(20,150),(100,7),(200,2)]:
        repeat *= 1
        print '%5s' % size,
        sys.stdout.flush()

        a = random([size,size])

        print '| %6.2f ' % measure('eigvals(a)',repeat),
        sys.stdout.flush()

        print '   (secs for %s calls)' % (repeat)
