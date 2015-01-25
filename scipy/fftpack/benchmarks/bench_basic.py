""" Test functions for fftpack.basic module
"""

from __future__ import division, print_function, absolute_import

import sys
from numpy.testing import *
from scipy.fftpack import ifft, fft, fftn, irfft, rfft

from numpy import arange, asarray, zeros, dot, exp, pi, double, cdouble
import numpy.fft

from numpy.random import rand


def random(size):
    return rand(*size)


def direct_dft(x):
    x = asarray(x)
    n = len(x)
    y = zeros(n,dtype=cdouble)
    w = -arange(n)*(2j*pi/n)
    for i in range(n):
        y[i] = dot(exp(i*w),x)
    return y


def direct_idft(x):
    x = asarray(x)
    n = len(x)
    y = zeros(n,dtype=cdouble)
    w = arange(n)*(2j*pi/n)
    for i in range(n):
        y[i] = dot(exp(i*w),x)/n
    return y


class TestFft(TestCase):

    def bench_random(self):
        from numpy.fft import fft as numpy_fft
        print()
        print('                 Fast Fourier Transform')
        print('=================================================')
        print('      |    real input     |   complex input    ')
        print('-------------------------------------------------')
        print(' size |  scipy  |  numpy  |  scipy  |  numpy ')
        print('-------------------------------------------------')
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print('%5s' % size, end=' ')
            sys.stdout.flush()

            for x in [random([size]).astype(double),
                      random([size]).astype(cdouble)+random([size]).astype(cdouble)*1j
                      ]:
                if size > 500:
                    y = fft(x)
                else:
                    y = direct_dft(x)
                assert_array_almost_equal(fft(x),y)
                print('|%8.2f' % measure('fft(x)',repeat), end=' ')
                sys.stdout.flush()

                assert_array_almost_equal(numpy_fft(x),y)
                print('|%8.2f' % measure('numpy_fft(x)',repeat), end=' ')
                sys.stdout.flush()

            print(' (secs for %s calls)' % (repeat))
        sys.stdout.flush()


class TestIfft(TestCase):

    def bench_random(self):
        from numpy.fft import ifft as numpy_ifft
        print()
        print('       Inverse Fast Fourier Transform')
        print('===============================================')
        print('      |     real input    |    complex input   ')
        print('-----------------------------------------------')
        print(' size |  scipy  |  numpy  |  scipy  |  numpy  ')
        print('-----------------------------------------------')
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print('%5s' % size, end=' ')
            sys.stdout.flush()

            for x in [random([size]).astype(double),
                      random([size]).astype(cdouble)+random([size]).astype(cdouble)*1j
                      ]:
                if size > 500:
                    y = ifft(x)
                else:
                    y = direct_idft(x)
                assert_array_almost_equal(ifft(x),y)
                print('|%8.2f' % measure('ifft(x)',repeat), end=' ')
                sys.stdout.flush()

                assert_array_almost_equal(numpy_ifft(x),y)
                print('|%8.2f' % measure('numpy_ifft(x)',repeat), end=' ')
                sys.stdout.flush()

            print(' (secs for %s calls)' % (repeat))
        sys.stdout.flush()


class TestRfft(TestCase):

    def bench_random(self):
        from numpy.fft import rfft as numpy_rfft
        print()
        print('Fast Fourier Transform (real data)')
        print('==================================')
        print(' size |  scipy  |  numpy  ')
        print('----------------------------------')
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print('%5s' % size, end=' ')
            sys.stdout.flush()

            x = random([size]).astype(double)
            print('|%8.2f' % measure('rfft(x)',repeat), end=' ')
            sys.stdout.flush()

            print('|%8.2f' % measure('numpy_rfft(x)',repeat), end=' ')
            sys.stdout.flush()

            print(' (secs for %s calls)' % (repeat))
        sys.stdout.flush()


class TestIrfft(TestCase):

    def bench_random(self):
        from numpy.fft import irfft as numpy_irfft

        print()
        print('Inverse Fast Fourier Transform (real data)')
        print('==================================')
        print(' size |  scipy  |  numpy  ')
        print('----------------------------------')
        for size,repeat in [(100,7000),(1000,2000),
                            (256,10000),
                            (512,10000),
                            (1024,1000),
                            (2048,1000),
                            (2048*2,500),
                            (2048*4,500),
                            ]:
            print('%5s' % size, end=' ')
            sys.stdout.flush()

            x = random([size]).astype(double)
            x1 = zeros(size/2+1,dtype=cdouble)
            x1[0] = x[0]
            for i in range(1,size//2):
                x1[i] = x[2*i-1] + 1j * x[2*i]
            if not size % 2:
                x1[-1] = x[-1]
            y = irfft(x)

            print('|%8.2f' % measure('irfft(x)',repeat), end=' ')
            sys.stdout.flush()

            assert_array_almost_equal(numpy_irfft(x1,size),y)
            print('|%8.2f' % measure('numpy_irfft(x1,size)',repeat), end=' ')
            sys.stdout.flush()

            print(' (secs for %s calls)' % (repeat))

        sys.stdout.flush()


class TestFftn(TestCase):

    def bench_random(self):
        from numpy.fft import fftn as numpy_fftn
        print()
        print('    Multi-dimensional Fast Fourier Transform')
        print('===================================================')
        print('          |    real input     |   complex input    ')
        print('---------------------------------------------------')
        print('   size   |  scipy  |  numpy  |  scipy  |  numpy ')
        print('---------------------------------------------------')
        for size,repeat in [((100,100),100),((1000,100),7),
                            ((256,256),10),
                            ((512,512),3),
                            ]:
            print('%9s' % ('%sx%s' % size), end=' ')
            sys.stdout.flush()

            for x in [random(size).astype(double),
                      random(size).astype(cdouble)+random(size).astype(cdouble)*1j
                      ]:
                y = fftn(x)
                #if size > 500: y = fftn(x)
                #else: y = direct_dft(x)
                assert_array_almost_equal(fftn(x),y)
                print('|%8.2f' % measure('fftn(x)',repeat), end=' ')
                sys.stdout.flush()

                assert_array_almost_equal(numpy_fftn(x),y)
                print('|%8.2f' % measure('numpy_fftn(x)',repeat), end=' ')
                sys.stdout.flush()

            print(' (secs for %s calls)' % (repeat))

        sys.stdout.flush()


if __name__ == "__main__":
    Tester().bench()
