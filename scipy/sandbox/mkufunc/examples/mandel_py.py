#!/usr/bin/env python
from numpy import array
from mkufunc.api import mkufunc

D = 1000


@mkufunc([(float, float, int)])
def mandel(cr, ci):
    d = 1
    zr = cr
    zi = ci
    for d in xrange(1, D):
        zr2 = zr * zr
        zi2 = zi * zi
        if zr2 + zi2 > 16:
            return d
        zi = 2.0 * zr * zi + ci
        zr = zr2 - zi2 + cr
    else:
        return -1

    

if __name__ == '__main__':

    assert mandel(-1, .3) == 36
    assert mandel(0, 0) == -1
    assert mandel(10, 10) == 1
    assert all(mandel(array([-1, 0, 10]), array([.3, 0, 10])) ==
               array([36, -1,  1]))
