#!/usr/bin/env python
from math import sqrt
import time, sys

src = 0
if len(sys.argv) == 2 and sys.argv[1] == 's':
    src = 1


import sys
sys.path.append('../mkufunc')
from fast_vectorize import fast_vectorize



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
print count_primes(100000)
print 'Python: %.6f sec' % (time.time() - start_time)


count_primes = fast_vectorize(int, src=src)(count_primes)


start_time = time.time()
print count_primes(100000)
print 'Compiled: %.6f sec' % (time.time() - start_time)
