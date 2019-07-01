import scipy._uarray as ua

from ._basic import (
    fft, ifft, fft2,ifft2, fftn, ifftn,
    rfft, irfft, rfft2, irfft2, rfftn, irfftn)
from ._realtransforms import dct, idct, dst, idst, dctn, idctn, dstn, idstn
import scipy.fftpack as _fftpack
from . import _pocketfft

class BackendError(RuntimeError):
    pass

class BackendWarning(UserWarning):
    pass

'''
def _wrap_fallback(fname, on_missing):
    """Wrap fallback function with error reporting according to on_missing"""
    message = "The current backend does not implement '{}'".format(fname)

    if on_missing == 'fallback':
        return cfg._fallbacks[fname]
    elif on_missing == 'warn':
        def wrap_warn(*args, **kwargs):
            from warnings import warn
            warn(message, BackendWarning)
            return cfg._fallbacks[fname](*args, **kwargs)
        return wrap_warn
    elif on_missing == 'raise':
        def wrap_raise(*args, **kwargs):
            raise BackendError(message)
        return wrap_raise

    raise ValueError("Unrecognized on_missing type '{}'".format(on_missing))
'''



def set_backend(backend, on_missing='fallback'):
    """Sets the current fft backend

    Parameters
    ----------

    backend: string
        Can either be one of the known backends {'scipy'}, or a
        module import specification of the form 'module://example.fft'
    on_missing: {'fallback', 'warn', 'raise'}, optional
        Behavior when the backend does not provide a given function:
        - 'fallback': silently use the built-in SciPy function
        - 'warn': emit a warning, then use SciPy's default
        - 'raise': raise an error

    Raises
    ------
    ImportError: If the specified backend could not be imported
    ValueError: If an invalid parameter is given

    """

    assert(backend.__ua_domain__ == 'scipy.fft')
    ua.set_global_backend(backend)


def backend(backend, on_missing='fallback'):
    """Context manager to change the backend within a fixed scope

    Upon entering a ``with`` statement, the current backend is changed. Upon
    exit, the backend is reset to the state before entering the scope.

    Parameters
    ---------
    backend: string
        Can either be one of the known backends {'scipy'}, or a
        module import specification of the form 'module://example.fft'
    on_missing: {'fallback', 'warn', 'raise'}, optional
        Behavior when the backend does not provide a given function:

    Examples
    --------
    >>> with scipy.fft.backend('scipy'):
    >>>     pass

    """
    return ua.set_backend(backend)


class _ScipyBackend:
    __ua_domain__ = "scipy.fft"

    implemented = {
        fft: _pocketfft.fft,
        fft2: _pocketfft.fft2,
        fftn: _pocketfft.fftn,

        ifft: _pocketfft.ifft,
        ifft2: _pocketfft.ifft2,
        ifftn: _pocketfft.ifftn,

        rfft: _pocketfft.rfft,
        rfft2: _pocketfft.rfft2,
        rfftn: _pocketfft.rfftn,

        irfft: _pocketfft.irfft,
        irfft2: _pocketfft.irfft2,
        irfftn: _pocketfft.irfftn,

        dct: _fftpack.dct,
        idct: _fftpack.idct,
        dctn: _fftpack.dctn,
        idctn: _fftpack.idctn,

        dst: _fftpack.dst,
        idst: _fftpack.idst,
        dstn: _fftpack.dstn,
        idstn: _fftpack.idstn,
    }


    @staticmethod
    def __ua_function__(method, args, kwargs):
        fn = _ScipyBackend.implemented.get(method)
        if fn is None:
            return NotImplemented
        return fn(*args, **kwargs)


set_backend(_ScipyBackend())
