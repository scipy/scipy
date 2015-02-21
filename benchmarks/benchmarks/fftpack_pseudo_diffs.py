""" Benchmark functions for fftpack.pseudo_diffs module
"""
from __future__ import division, absolute_import, print_function

import sys

from numpy import arange, sin, cos, pi, exp, tanh, sign

try:
    from scipy.fftpack import diff, fft, ifft, tilbert, hilbert, shift, fftfreq
except ImportError:
    pass

from .common import Benchmark


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


class Bench(Benchmark):
    params = [
        [100, 256, 512, 1000, 1024, 2048, 2048*2, 2048*4],
        ['fft', 'direct'],
    ]
    param_names = ['size', 'type']

    def setup(self, size, type):
        size = int(size)

        x = arange(size)*2*pi/size
        a = 1
        self.a = a
        if size < 2000:
            self.f = sin(x)*cos(4*x)+exp(sin(3*x))
            self.sf = sin(x+a)*cos(4*(x+a))+exp(sin(3*(x+a)))
        else:
            self.f = sin(x)*cos(4*x)
            self.sf = sin(x+a)*cos(4*(x+a))

    def time_diff(self, size, soltype):
        if soltype == 'fft':
            diff(self.f, 3)
        else:
            direct_diff(self.f,3)

    def time_tilbert(self, size, soltype):
        if soltype == 'fft':
            tilbert(self.f, 1)
        else:
            direct_tilbert(self.f, 1)

    def time_hilbert(self, size, soltype):
        if soltype == 'fft':
            hilbert(self.f)
        else:
            direct_hilbert(self.f)

    def time_shift(self, size, soltype):
        if soltype == 'fft':
            shift(self.f, self.a)
        else:
            direct_shift(self.f, self.a)
