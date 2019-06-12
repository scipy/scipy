""" Test functions for fftpack.basic module
"""
from __future__ import division, absolute_import, print_function

from numpy import arange, asarray, zeros, dot, exp, pi, double, cdouble

from numpy.random import rand

import scipy.fftpack
import numpy.fft
try:
    import scipy.fft
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
        [100, 256, 313, 512, 1000, 1024, 2048, 2048*2, 2048*4],
        ['real', 'cmplx'],
        ['scipy.fftpack', 'scipy.fft', 'numpy.fft']
    ]
    param_names = ['size', 'type', 'module']

    def setup(self, size, cmplx, module):
        if cmplx == 'cmplx':
            self.x = random([size]).astype(cdouble)+random([size]).astype(cdouble)*1j
        else:
            self.x = random([size]).astype(double)

        self.fft = eval(module + '.fft')
        self.ifft = eval(module + '.ifft')

    def time_fft(self, size, cmplx, module):
        self.fft(self.x)

    def time_ifft(self, size, cmplx, module):
        self.ifft(self.x)


class RFft(Benchmark):
    params = [
        [100, 256, 313, 512, 1000, 1024, 2048, 2048*2, 2048*4],
        ['scipy.fftpack', 'scipy.fft', 'numpy.fft']
    ]
    param_names = ['size', 'module']

    def setup(self, size, module):
        self.x = random([size]).astype(double)

        self.rfft = eval(module + '.rfft')
        self.irfft = eval(module + '.irfft')

        self.y = self.rfft(self.x)

    def time_rfft(self, size, module):
        self.rfft(self.x)

    def time_irfft(self, size, module):
        self.irfft(self.y)


class Fftn(Benchmark):
    params = [
        ["100x100", "313x100", "1000x100", "256x256", "512x512"],
        ['real', 'cmplx'],
        ['scipy.fftpack', 'scipy.fft', 'numpy.fft']
    ]
    param_names = ['size', 'type', 'module']

    def setup(self, size, cmplx, module):
        size = list(map(int, size.split("x")))

        if cmplx != 'cmplx':
            self.x = random(size).astype(double)
        else:
            self.x = random(size).astype(cdouble)+random(size).astype(cdouble)*1j

        self.fftn = eval(module + '.fftn')

    def time_fftn(self, size, cmplx, module):
        self.fftn(self.x)
