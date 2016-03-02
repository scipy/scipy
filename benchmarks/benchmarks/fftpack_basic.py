""" Test functions for fftpack.basic module
"""
from __future__ import division, absolute_import, print_function

from numpy import arange, asarray, zeros, dot, exp, pi, double, cdouble
import numpy.fft

from numpy.random import rand

try:
    from scipy.fftpack import ifft, fft, fftn, irfft, rfft
except ImportError:
    pass

from .common import Benchmark


def random(size):
    return rand(*size)


def direct_dft(x):
    x = asarray(x)
    n = len(x)
    y = zeros(n, dtype=cdouble)
    w = -arange(n)*(2j*pi/n)
    for i in range(n):
        y[i] = dot(exp(i*w), x)
    return y


def direct_idft(x):
    x = asarray(x)
    n = len(x)
    y = zeros(n, dtype=cdouble)
    w = arange(n)*(2j*pi/n)
    for i in range(n):
        y[i] = dot(exp(i*w), x)/n
    return y


class Fft(Benchmark):
    params = [
        [100, 256, 512, 1000, 1024, 2048, 2048*2, 2048*4],
        ['real', 'cmplx'],
        ['scipy', 'numpy']
    ]
    param_names = ['size', 'type', 'module']

    def setup(self, size, cmplx, module):
        if cmplx == 'cmplx':
            self.x = random([size]).astype(cdouble)+random([size]).astype(cdouble)*1j
        else:
            self.x = random([size]).astype(double)

    def time_fft(self, size, cmplx, module):
        if module == 'numpy':
            numpy.fft.fft(self.x)
        else:
            fft(self.x)

    def time_ifft(self, size, cmplx, module):
        if module == 'numpy':
            numpy.fft.ifft(self.x)
        else:
            ifft(self.x)


class RFft(Benchmark):
    params = [
        [100, 256, 512, 1000, 1024, 2048, 2048*2, 2048*4],
        ['scipy', 'numpy']
    ]
    param_names = ['size', 'module']

    def setup(self, size, module):
        self.x = random([size]).astype(double)

    def time_rfft(self, size, module):
        if module == 'numpy':
            numpy.fft.rfft(self.x)
        else:
            rfft(self.x)

    def time_irfft(self, size, module):
        if module == 'numpy':
            numpy.fft.irfft(self.x)
        else:
            irfft(self.x)


class Fftn(Benchmark):
    params = [
        ["100x100", "1000x100", "256x256", "512x512"],
        ['real', 'cmplx'],
        ['scipy', 'numpy']
    ]
    param_names = ['size', 'type', 'module']

    def setup(self, size, cmplx, module):
        size = map(int, size.split("x"))

        if cmplx != 'cmplx':
            self.x = random(size).astype(double)
        else:
            self.x = random(size).astype(cdouble)+random(size).astype(cdouble)*1j

    def time_fftn(self, size, cmplx, module):
        if module == 'numpy':
            numpy.fft.fftn(self.x)
        else:
            fftn(self.x)
