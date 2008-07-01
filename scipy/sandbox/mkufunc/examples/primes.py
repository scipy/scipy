#!/usr/bin/env python
from math import sqrt
import time

from numpy import arange

from mkufunc.api import mkufunc


def is_prime(n):
    if n < 2:
        return 0
    for i in xrange(2, min(n, int(sqrt(n)+2.0))):
        if n %i == 0:
            return 0
    return 1


start_time = time.time()
assert sum(is_prime(n) for n in xrange(1000000)) == 78498
print 'Python: %.6f sec' % (time.time() - start_time)


is_prime = mkufunc(int)(is_prime)


start_time = time.time()
assert is_prime(arange(1000000)).sum() == 78498
print 'Compiled: %.6f sec' % (time.time() - start_time)
