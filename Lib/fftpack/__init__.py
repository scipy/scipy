"""
Discrete Fourier Transforms - fft.py 

The underlying code for these functions is the FFTPACK package.

fft(a, n=None, axis=-1) 
inverse_fft(a, n=None, axis=-1) 
real_fft(a, n=None, axis=-1) 
inverse_real_fft(a, n=None, axis=-1)
hermite_fft(a, n=None, axis=-1)
inverse_hermite_fft(a, n=None, axis=-1)
fftnd(a, s=None, axes=None)
inverse_fftnd(a, s=None, axes=None)
real_fftnd(a, s=None, axes=None)
inverse_real_fftnd(a, s=None, axes=None)
fft2d(a, s=None, axes=(-2,-1)) 
inverse_fft2d(a, s=None, axes=(-2, -1))
real_fft2d(a, s=None, axes=(-2,-1)) 
inverse_real_fft2d(a, s=None, axes=(-2, -1))
"""
from fft import *

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------

def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite())
    return runner

def test_suite(level=1):
    import scipy_base.testing
    import scipy_base
    this_mod = scipy_base
    # ieee_754 gets tested in the type_check module.
    # testing is the module that actually does all the testing...
    ignore = ['ieee_754','testing']
    return scipy_base.testing.harvest_test_suites(this_mod,ignore = ignore,
                                                  level=level)
