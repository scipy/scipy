""" Machine limits for Float32 and Float64.
"""

import warnings
warnings.warn('limits module is deprecated, please use numpy.finfo instead',
              DeprecationWarning)


__all__ = ['float_epsilon','float_tiny','float_min',
           'float_max','float_precision','float_resolution',
           'single_epsilon','single_tiny','single_min','single_max',
           'single_precision','single_resolution',
           'double_epsilon','double_tiny','double_min','double_max',
           'double_precision','double_resolution']


from numpy import finfo, single, float_

single_epsilon = finfo(single).eps
single_tiny = finfo(single).tiny
single_max = finfo(single).max
single_min = -single_max
single_precision = finfo(single).precision
single_resolution = finfo(single).resolution

double_epsilon = float_epsilon = finfo(float_).eps
double_tiny = float_tiny = finfo(float_).tiny
double_max = float_max = finfo(float_).max
double_min = float_min = -float_max
double_precision = float_precision = finfo(float_).precision
double_resolution = float_resolution = finfo(float_).resolution

if __name__ == '__main__':
    print 'single epsilon:',single_epsilon
    print 'single tiny:',single_tiny
    print 'float epsilon:',float_epsilon
    print 'float tiny:',float_tiny
