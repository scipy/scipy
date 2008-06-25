#!/usr/bin/env python
from math import sin, cos, pi

from numpy import array

from mkufunc import mkufunc

int_const = 42

float_const = 3.14

def my_sqr(x):
    return x * x / 1.23456

@mkufunc([(float, float), (int, int)])
def bar(n):
    "Bar docstring"
    if n < 10:
        for i in xrange(10):
            n += i*i
        return n
    elif n == 10:
        return int_const
    elif n == 11:
        return float_const
    elif n == 12:
        return cos(pi)
        #return 1
    elif n > 12:
        return my_sqr(n)
    else:
        return 5

    
#@mkufunc(float)
#def baz(n):
#    "Baz docstring"
#    return n * n + 1000


#print bar

x = array([0.0, 1.0, 2.5, 12.0])
print "x =", x, x.dtype
y = bar(x)
print "y =", y, y.dtype

print bar(5)
print bar(15)
print bar(10)
print bar(11)
print bar(12)
print bar(12.5)
