import scipy._lib.uarray as ua
from scipy.linalg import _api


__all__ = ['register_backend', 'set_backend',
           'set_global_backend', 'skip_backend']


class _ScipyLinalgBackend:
    __ua_domain__ = "numpy.scipy.linalg"

    @staticmethod
    def __ua_function__(method, args, kwargs):
        fn = getattr(_api, method.__name__, None)

        if fn is None:
            return NotImplemented
        return fn(*args, **kwargs)


_named_backends = {
    'scipy': _ScipyLinalgBackend,
}


def _backend_from_arg(backend):

    if isinstance(backend, str):
        try:
            backend = _named_backends[backend]
        except KeyError as e:
            raise ValueError('Unknown backend {}'.format(backend)) from e

    if backend.__ua_domain__ != 'numpy.scipy.linalg':
        raise ValueError('Backend does not implement "numpy.scipy.linalg"')

    return backend


def set_global_backend(backend, coerce=False, only=False, try_last=False):
    backend = _backend_from_arg(backend)
    ua.set_global_backend(backend, coerce=coerce, only=only, try_last=try_last)


def register_backend(backend):
    backend = _backend_from_arg(backend)
    ua.register_backend(backend)


def set_backend(backend, coerce=False, only=False):
    backend = _backend_from_arg(backend)
    return ua.set_backend(backend, coerce=coerce, only=only)


def skip_backend(backend):
    backend = _backend_from_arg(backend)
    return ua.skip_backend(backend)


set_global_backend('scipy', try_last=True)
