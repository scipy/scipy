import scipy._lib.uarray as ua

from ._basic import (
    fft, ifft, fft2,ifft2, fftn, ifftn,
    rfft, irfft, rfft2, irfft2, rfftn, irfftn)
from ._realtransforms import dct, idct, dst, idst, dctn, idctn, dstn, idstn
import scipy.fftpack as _fftpack
from . import _pocketfft


class _ScipyBackend:
    """
    The default backend for fft calculations
    """
    __ua_domain__ = "scipy.fft"

    @staticmethod
    def __ua_function__(method, args, kwargs):
        fn = (getattr(_pocketfft, method.__name__, None)
              or getattr(_fftpack, method.__name__, None))

        if fn is None:
            return NotImplemented
        return fn(*args, **kwargs)


_named_backends = {
    'scipy': _ScipyBackend,
}

def _backend_from_arg(backend):
    """Maps strings to known backends and validates the backend"""

    if isinstance(backend, str):
        try:
            backend = _named_backends[backend]
        except KeyError:
            raise ValueError('Unknown backend {}'.format(backend))

    if backend.__ua_domain__ != 'scipy.fft':
        raise ValueError('Backend does not implement "scipy.fft"')

    return backend


def set_global_backend(backend):
    """Sets the global fft backend

    The global backend has higher priority than registered backends, but lower
    priority than context-specific backends set with `set_backend`.

    Parameters
    ----------

    backend: {object, 'scipy'}

        The backend to use.
        Can either be a ``str`` containing the name of a known backend
        {'scipy'}, or an object that implements the uarray protocol.

    Raises
    ------

    ValueError: If the backend does not implement ``scipy.fft``

    Notes
    -----

    This will overwrite the previously set global backend, which by default is
    the SciPy implementation.

    """

    backend = _backend_from_arg(backend)
    ua.set_global_backend(backend)


def register_backend(backend):
    """
    Register a backend for permanent use.

    Registered backends have the lowest priority and will be tried after the
    global backend.

    Parameters
    ----------

    backend: {object, 'scipy'}

        The backend to use.
        Can either be a ``str`` containing the name of a known backend
        {'scipy'}, or an object that implements the uarray protocol.

    Raises
    ------
    ValueError: If the backend does not implement ``scipy.fft``

    """
    backend = _backend_from_arg(backend)
    ua.register_backend(backend)


def set_backend(backend, coerce=False, only=False):
    """Context manager to set the backend within a fixed scope.

    Upon entering the ``with`` statement, given backend will be added to the
    list of available backends with the highest priority. Upon exit, the
    backend is reset to the state before entering the scope.

    Parameters
    ----------
    backend: {object, 'scipy'}

        The backend to use.
        Can either be a ``str`` containing the name of a known backend
        {'scipy'}, or an object that implements the uarray protocol.

    coerce: bool, optional

        Whether to allow expensive conversions for the ``x`` parameter. e.g.
        copying a numpy array to the GPU for a CuPy backend. Implies ``only``.

    only: bool, optional

       If only is ``True`` and this backend returns ``NotImplemented`` then a
       BackendNotImplemented error will be raised immediately. Ignoring any
       lower priority backends.

    Examples
    --------
    >>> import scipy.fft as fft
    >>> with fft.set_backend('scipy', only=True):
    >>>     fft.fft([1])  # Always calls the scipy implementation
    array([1.+0.j])

    """
    backend = _backend_from_arg(backend)
    return ua.set_backend(backend, coerce=coerce, only=only)


set_global_backend(_ScipyBackend())
