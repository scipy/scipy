#!/usr/bin/env python
from math import sqrt
import time

from mkufunc.api import mkufunc


def count_primes(N):
    res = 0
    for n in xrange(2, N):
        for i in xrange(2, min(n, int(sqrt(n)+2.0))):
            if n % i == 0:
                break
        else:
            res += 1
    return res


start_time = time.time()
assert count_primes(1000000) == 78498
print 'Python: %.6f sec' % (time.time() - start_time)


count_primes = mkufunc(int)(count_primes)


start_time = time.time()
assert count_primes(1000000) == 78498
print 'Compiled: %.6f sec' % (time.time() - start_time)
