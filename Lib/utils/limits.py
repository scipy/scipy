""" Machine limits for Float32 and Float64.
"""

__all__ = ['float_epsilon','float_tiny','float_min',
           'float_max','float_precision','float_resolution',
           'double_epsilon','double_tiny','double_min','double_max',
           'double_precision','double_resolution']

from scipy.base import finfo, single, float_

single_epsilon = finfo(single).epsilon
single_tiny = finfo(single).tiny
single_max = finfo(single).huge
single_min = -single_max
single_precision = finfo(single).precision
single_resolution = finfo(single).resolution

float_epsilon = finfo(float_).epsilon
float_tiny = finfo(float_).tiny
float_max = finfo(float_).huge
float_min = -float_max
float_precision = finfo(float_).precision
float_resolution = finfo(float_).resolution

if __name__ == '__main__':
    print 'single epsilon:',single_epsilon
    print 'single tiny:',single_tiny
    print 'float epsilon:',float_epsilon
    print 'float tiny:',float_tiny

