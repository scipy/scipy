import numpy as np
import scipy.fft

class _MockFunction:
    def __init__(self, return_value = None):
        self.number_calls = 0
        self.return_value = return_value
        self.last_args = ([], {})


    def __call__(self, *args, **kwargs):
        self.number_calls += 1
        self.last_args = (args, kwargs)
        return self.return_value


fft = _MockFunction(np.random.random(10))
fft2 = _MockFunction(np.random.random(10))
fftn = _MockFunction(np.random.random(10))

ifft = _MockFunction(np.random.random(10))
ifft2 = _MockFunction(np.random.random(10))
ifftn = _MockFunction(np.random.random(10))

rfft = _MockFunction(np.random.random(10))
rfft2 = _MockFunction(np.random.random(10))
rfftn = _MockFunction(np.random.random(10))

irfft = _MockFunction(np.random.random(10))
irfft2 = _MockFunction(np.random.random(10))
irfftn = _MockFunction(np.random.random(10))

dct = _MockFunction(np.random.random(10))
idct = _MockFunction(np.random.random(10))
dctn = _MockFunction(np.random.random(10))
idctn = _MockFunction(np.random.random(10))

dst = _MockFunction(np.random.random(10))
idst = _MockFunction(np.random.random(10))
dstn = _MockFunction(np.random.random(10))
idstn = _MockFunction(np.random.random(10))

implemented = {}

implemented[scipy.fft.fft] = fft
implemented[scipy.fft.fft2] = fft2
implemented[scipy.fft.fftn] = fftn

implemented[scipy.fft.ifft] = ifft
implemented[scipy.fft.ifft2] = ifft2
implemented[scipy.fft.ifftn] = ifftn

implemented[scipy.fft.rfft] = rfft
implemented[scipy.fft.rfft2] = rfft2
implemented[scipy.fft.rfftn] = rfftn

implemented[scipy.fft.irfft] = irfft
implemented[scipy.fft.irfft2] = irfft2
implemented[scipy.fft.irfftn] = irfftn

implemented[scipy.fft.dct] = dct
implemented[scipy.fft.idct] = idct
implemented[scipy.fft.dctn] = dctn
implemented[scipy.fft.idctn] = idctn

implemented[scipy.fft.dst] = dst
implemented[scipy.fft.idst] = idst
implemented[scipy.fft.dstn] = dstn
implemented[scipy.fft.idstn] = idstn

class MockBackend:
    __ua_domain__ = "scipy.fft"

    @staticmethod
    def __ua_function__(method, args, kwargs):
        fn = implemented.get(method)
        if fn is None:
            return NotImplemented
        return fn(*args, **kwargs)
