""" Benchmark functions for linalg.decomp module

"""
import sys

from numpy import linalg
from numpy.testing import *

from scipy.lib.six import print_

def random(size):
    return rand(*size)

def bench_random():
    Numeric_eigvals = linalg.eigvals
    print_()
    print_('           Finding matrix eigenvalues')
    print_('      ==================================')
    print_('      |    contiguous     ')#'|   non-contiguous '
    print_('----------------------------------------------')
    print_(' size |  scipy  ')#'| core |  scipy  | core '

    for size,repeat in [(20,150),(100,7),(200,2)]:
        repeat *= 1
        print_('%5s' % size, end=' ')
        sys.stdout.flush()

        a = random([size,size])

        print_('| %6.2f ' % measure('eigvals(a)',repeat), end=' ')
        sys.stdout.flush()

        print_('   (secs for %s calls)' % (repeat))
