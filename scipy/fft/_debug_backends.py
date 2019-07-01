import numpy as np
try:
    import pyfftw.interfaces.numpy_fft as pyfftw_fft
except ImportError:
    pyfftw_fft = {}

class NumPyBackend:
    """Backend that uses numpy.fft"""
    __ua_domain__ = "scipy.fft"

    @staticmethod
    def __ua_function__(method, args, kwargs):
        kwargs.pop("overwrite_x", None)

        fn = getattr(np.fft, method.__name__, None)
        return (NotImplemented if fn is None
                else fn(*args, **kwargs))


class PyfftwBackend:
    __ua_domain__ = 'scipy.fft'

    @staticmethod
    def __ua_function__(method, args, kwargs):
        kwargs.pop('overwrite_x', None)

        fn = getattr(pyfftw_fft, method.__name__, None)
        return (NotImplemented if fn is None
                else fn(*args, **kwargs))


class EchoBackend:
    """Backend that just prints the __ua_function__ arguements"""
    __ua_domain__ = "scipy.fft"

    @staticmethod
    def __ua_function__(method, args, kwargs):
        print(method, args, kwargs, sep='\n')
