""" Test functions for fftpack.basic module
"""
from __future__ import division, absolute_import, print_function

from numpy import arange, asarray, zeros, dot, exp, pi, double, cdouble

from numpy.random import rand

import scipy.fftpack
import numpy.fft
try:
    import scipy.fft as scipy_fft
    has_scipy_fft = True
except ImportError:
    scipy_fft = {}
    has_scipy_fft = False

from .common import Benchmark

try:
    import pyfftw.interfaces.numpy_fft as pyfftw_fft
    import pyfftw
    pyfftw.interfaces.cache.enable()
    has_pyfftw = True
except ImportError:
    pyfftw_fft = {}
    has_pyfftw = False

class PyfftwBackend:
    """Backend for pyfftw"""
    __ua_domain__ = 'numpy.scipy.fft'

    @staticmethod
    def __ua_function__(method, args, kwargs):
        kwargs.pop('overwrite_x', None)

        fn = getattr(pyfftw_fft, method.__name__, None)
        return (NotImplemented if fn is None
                else fn(*args, **kwargs))


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


def get_module(mod_name):
    module_map = {
        'scipy.fftpack': scipy.fftpack,
        'scipy.fft': scipy_fft,
        'numpy.fft': numpy.fft
    }

    if not has_scipy_fft and mod_name == 'scipy.fft':
        raise NotImplementedError

    return module_map[mod_name]


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

        module = get_module(module)
        self.fft = getattr(module, 'fft')
        self.ifft = getattr(module, 'ifft')

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

        module = get_module(module)
        self.rfft = getattr(module, 'rfft')
        self.irfft = getattr(module, 'irfft')

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

        self.fftn = getattr(get_module(module), 'fftn')

    def time_fftn(self, size, cmplx, module):
        self.fftn(self.x)


class FftBackends(Benchmark):
    params = [
        [100, 256, 313, 512, 1000, 1024, 2048, 2048*2, 2048*4],
        ['real', 'cmplx'],
        ['pocketfft', 'pyfftw', 'numpy', 'direct']
    ]
    param_names = ['size', 'type', 'backend']

    def setup(self, size, cmplx, backend):
        import scipy.fft
        if cmplx == 'cmplx':
            self.x = random([size]).astype(cdouble)+random([size]).astype(cdouble)*1j
        else:
            self.x = random([size]).astype(double)

        self.fft = scipy.fft.fft
        self.ifft = scipy.fft.ifft

        if backend == 'pocketfft':
            scipy.fft.set_global_backend('scipy')
        elif backend == 'pyfftw':
            if not has_pyfftw:
                raise NotImplementedError
            scipy.fft.set_global_backend(PyfftwBackend)
        elif backend == 'numpy':
            from scipy.fft._debug_backends import NumPyBackend
            scipy.fft.set_global_backend(NumPyBackend)
        elif backend == 'direct':
            import scipy.fft._pocketfft
            self.fft = scipy.fft._pocketfft.fft
            self.ifft = scipy.fft._pocketfft.ifft

    def time_fft(self, size, cmplx, module):
        self.fft(self.x)

    def time_ifft(self, size, cmplx, module):
        self.ifft(self.x)


class FftnBackends(Benchmark):
    params = [
        ["100x100", "313x100", "1000x100", "256x256", "512x512"],
        ['real', 'cmplx'],
        ['pocketfft', 'pyfftw', 'numpy', 'direct']
    ]
    param_names = ['size', 'type', 'backend']

    def setup(self, size, cmplx, backend):
        import scipy.fft
        size = list(map(int, size.split("x")))

        if cmplx == 'cmplx':
            self.x = random(size).astype(double)+random(size).astype(double)*1j
        else:
            self.x = random(size).astype(double)

        self.fftn = scipy.fft.fftn
        self.ifftn = scipy.fft.ifftn

        if backend == 'pocketfft':
            scipy.fft.set_global_backend('scipy')
        elif backend == 'pyfftw':
            if not has_pyfftw:
                raise NotImplementedError
            scipy.fft.set_global_backend(PyfftwBackend)
        elif backend == 'numpy':
            from scipy.fft._debug_backends import NumPyBackend
            scipy.fft.set_global_backend(NumPyBackend)
        elif backend == 'direct':
            import scipy.fft._pocketfft
            self.fftn = scipy.fft._pocketfft.fftn
            self.ifftn = scipy.fft._pocketfft.ifftn

    def time_fft(self, size, cmplx, module):
        self.fftn(self.x)

    def time_ifft(self, size, cmplx, module):
        self.ifftn(self.x)
