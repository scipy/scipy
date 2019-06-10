import numpy as np

class NumPyBackend:
    """Backend that uses numpy.fft"""
    __ua_domain__ = "scipy.fft"

    @staticmethod
    def __ua_convert__(value, type, coerce):
        if type is np.array:
            if not coerce:
                return value
            return np.asarray(value)
        return NotImplemented

    @staticmethod
    def __ua_function__(method, args, kwargs):
        kwargs.pop("overwrite_x", None)

        fn = getattr(np.fft, method.__name__)
        return (NotImplemented if fn is None
                else fn(*args, **kwargs))


class EchoBackend:
    """Backend that just prints the __ua_function__ arguements"""
    __ua_domain__ = "scipy.fft"

    @staticmethod
    def __ua_convert__(value, type, coerce):
        return value

    @staticmethod
    def __ua_function__(method, args, kwargs):
        print(method, args, kwargs, sep='\n')
