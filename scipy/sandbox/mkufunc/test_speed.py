#!/usr/bin/env python
from math import sin, cos
import time

from numpy import arange, vectorize

from mkufunc import mkufunc


def f(x):
    return 4.2 * x * x + 3.7 * x + 1.5


ufunc = mkufunc([(float, float)])(f)

vfunc = vectorize(f)


x = arange(0, 1000, 0.001)    #print "x =", x, x.dtype

start_time = time.time()
y = 4.2 * x * x + 3.7 * x + 1.5
n_time = time.time() - start_time
print 'numpy: %.6f sec' % n_time

start_time = time.time()
y = vfunc(x)
v_time = time.time() - start_time
print 'vectorize: %.6f sec' % v_time

start_time = time.time()
y = ufunc(x)
u_time = time.time() - start_time
print 'mkufunc: %.6f sec' % u_time

print "speedup over numpy:", n_time/u_time
print "speedup over vectorize:", v_time/u_time
