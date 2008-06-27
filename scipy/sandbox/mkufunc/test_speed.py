#!/usr/bin/env python
from math import sin, cos
import time

from numpy import arange, vectorize, allclose
from scipy import weave

from mkufunc import mkufunc


def f(x):
    return 4.2 * x * x + 3.7 * x + 1.5


vfunc = vectorize(f)

ufunc = mkufunc([(float, float)])(f)


x = arange(0, 1000, 0.001)    #print "x =", x, x.dtype

start_time = time.time()
b_y = x.copy()
weave.blitz("b_y[:] = 4.2 * x[:] * x[:] + 3.7 * x[:] + 1.5")
b_time = time.time() - start_time
print 'blitz: %.6f sec' % b_time

start_time = time.time()
n_y = f(x)
n_time = time.time() - start_time
print 'numpy: %.6f sec' % n_time

start_time = time.time()
v_y = vfunc(x)
v_time = time.time() - start_time
print 'vectorize: %.6f sec' % v_time

start_time = time.time()
u_y = ufunc(x)
u_time = time.time() - start_time
print 'mkufunc: %.6f sec' % u_time

print "speedup over blitz:",     b_time/u_time
print "speedup over numpy:",     n_time/u_time
print "speedup over vectorize:", v_time/u_time

assert allclose(b_y, n_y)
assert allclose(v_y, n_y)
assert allclose(u_y, n_y)
