""" Benchmark functions for fftpack.pseudo_diffs module
"""

from __future__ import division, print_function, absolute_import

import sys

from numpy import arange, sin, cos, pi, exp, tanh, sign

from numpy.testing import *
from scipy.fftpack import diff, fft, ifft, tilbert, hilbert, shift, fftfreq


def random(size):
    return rand(*size)


def direct_diff(x,k=1,period=None):
    fx = fft(x)
    n = len(fx)
    if period is None:
        period = 2*pi
    w = fftfreq(n)*2j*pi/period*n
    if k < 0:
        w = 1 / w**k
        w[0] = 0.0
    else:
        w = w**k
    if n > 2000:
        w[250:n-250] = 0.0
    return ifft(w*fx).real


def direct_tilbert(x,h=1,period=None):
    fx = fft(x)
    n = len(fx)
    if period is None:
        period = 2*pi
    w = fftfreq(n)*h*2*pi/period*n
    w[0] = 1
    w = 1j/tanh(w)
    w[0] = 0j
    return ifft(w*fx)


def direct_hilbert(x):
    fx = fft(x)
    n = len(fx)
    w = fftfreq(n)*n
    w = 1j*sign(w)
    return ifft(w*fx)


def direct_shift(x,a,period=None):
    n = len(x)
    if period is None:
        k = fftfreq(n)*1j*n
    else:
        k = fftfreq(n)*2j*pi/period*n
    return ifft(fft(x)*exp(k*a)).real


class TestDiff(TestCase):

    def bench_random(self):
        print()
        print('Differentiation of periodic functions')
        print('=====================================')
        print(' size  |  convolve |    naive')
        print('-------------------------------------')
        for size,repeat in [(100,1500),(1000,300),
                            (256,1500),
                            (512,1000),
                            (1024,500),
                            (2048,200),
                            (2048*2,100),
                            (2048*4,50),
                            ]:
            print('%6s' % size, end=' ')
            sys.stdout.flush()
            x = arange(size)*2*pi/size
            if size < 2000:
                f = sin(x)*cos(4*x)+exp(sin(3*x))
            else:
                f = sin(x)*cos(4*x)
            assert_array_almost_equal(diff(f,1),direct_diff(f,1))
            assert_array_almost_equal(diff(f,2),direct_diff(f,2))
            print('| %9.2f' % measure('diff(f,3)',repeat), end=' ')
            sys.stdout.flush()
            print('| %9.2f' % measure('direct_diff(f,3)',repeat), end=' ')
            sys.stdout.flush()
            print(' (secs for %s calls)' % (repeat))


class TestTilbert(TestCase):

    def bench_random(self):
        print()
        print(' Tilbert transform of periodic functions')
        print('=========================================')
        print(' size  | optimized |    naive')
        print('-----------------------------------------')
        for size,repeat in [(100,1500),(1000,300),
                            (256,1500),
                            (512,1000),
                            (1024,500),
                            (2048,200),
                            (2048*2,100),
                            (2048*4,50),
                            ]:
            print('%6s' % size, end=' ')
            sys.stdout.flush()
            x = arange(size)*2*pi/size
            if size < 2000:
                f = sin(x)*cos(4*x)+exp(sin(3*x))
            else:
                f = sin(x)*cos(4*x)
            assert_array_almost_equal(tilbert(f,1),direct_tilbert(f,1))
            print('| %9.2f' % measure('tilbert(f,1)',repeat), end=' ')
            sys.stdout.flush()
            print('| %9.2f' % measure('direct_tilbert(f,1)',repeat), end=' ')
            sys.stdout.flush()
            print(' (secs for %s calls)' % (repeat))


class TestHilbert(TestCase):

    def bench_random(self):
        print()
        print(' Hilbert transform of periodic functions')
        print('=========================================')
        print(' size  | optimized |    naive')
        print('-----------------------------------------')
        for size,repeat in [(100,1500),(1000,300),
                            (256,1500),
                            (512,1000),
                            (1024,500),
                            (2048,200),
                            (2048*2,100),
                            (2048*4,50),
                            ]:
            print('%6s' % size, end=' ')
            sys.stdout.flush()
            x = arange(size)*2*pi/size
            if size < 2000:
                f = sin(x)*cos(4*x)+exp(sin(3*x))
            else:
                f = sin(x)*cos(4*x)
            assert_array_almost_equal(hilbert(f),direct_hilbert(f))
            print('| %9.2f' % measure('hilbert(f)',repeat), end=' ')
            sys.stdout.flush()
            print('| %9.2f' % measure('direct_hilbert(f)',repeat), end=' ')
            sys.stdout.flush()
            print(' (secs for %s calls)' % (repeat))


class TestShift(TestCase):

    def bench_random(self):
        print()
        print(' Shifting periodic functions')
        print('==============================')
        print(' size  | optimized |    naive')
        print('------------------------------')
        for size,repeat in [(100,1500),(1000,300),
                            (256,1500),
                            (512,1000),
                            (1024,500),
                            (2048,200),
                            (2048*2,100),
                            (2048*4,50),
                            ]:
            print('%6s' % size, end=' ')
            sys.stdout.flush()
            x = arange(size)*2*pi/size
            a = 1
            if size < 2000:
                f = sin(x)*cos(4*x)+exp(sin(3*x))
                sf = sin(x+a)*cos(4*(x+a))+exp(sin(3*(x+a)))
            else:
                f = sin(x)*cos(4*x)
                sf = sin(x+a)*cos(4*(x+a))
            assert_array_almost_equal(direct_shift(f,1),sf)
            assert_array_almost_equal(shift(f,1),sf)
            print('| %9.2f' % measure('shift(f,a)',repeat), end=' ')
            sys.stdout.flush()
            print('| %9.2f' % measure('direct_shift(f,a)',repeat), end=' ')
            sys.stdout.flush()
            print(' (secs for %s calls)' % (repeat))

if __name__ == "__main__":
    run_module_suite()
