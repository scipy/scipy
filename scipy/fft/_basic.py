from scipy._lib._array_api import array_namespace
from . import _basic_np as npbasic
from . import _pocketfft as pfft
import numpy as np


def fft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
        plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of fft.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.fft(x, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                           workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.fft(x, n, axis, norm)
    x = np.asarray(x)
    y = pfft.fft(x, n, axis, norm)
    return xp.asarray(y)


def ifft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of ifft.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.ifft(x, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                            workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.ifft(x, n, axis, norm)
    x = np.asarray(x)
    y = pfft.ifft(x, n, axis, norm)
    return xp.asarray(y)


def fft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.fft2(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                     workers=workers, plan=plan)
    return xp.asarray(y)


def ifft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.ifft2(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                      workers=workers, plan=plan)
    return xp.asarray(y)


def fftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of fftn.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.fftn(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                            workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.fftn(x, s, axes, norm)
    x = np.asarray(x)
    y = pfft.fftn(x, s, axes, norm)
    return xp.asarray(y)


def ifftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of ifftn.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.ifftn(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                             workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.ifftn(x, s, axes, norm)
    x = np.asarray(x)
    y = pfft.ifftn(x, s, axes, norm)
    return xp.asarray(y)


def rfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of rfft.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.rfft(x, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                            workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.rfft(x, n, axis, norm)
    x = np.asarray(x)
    y = pfft.rfft(x, n, axis, norm)
    return xp.asarray(y)


def irfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of irfft.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.irfft(x, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                             workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.irfft(x, n, axis, norm)
    x = np.asarray(x)
    y = pfft.irfft(x, n, axis, norm)
    return xp.asarray(y)


def rfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.rfft2(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                      workers=workers, plan=plan)
    return xp.asarray(y)


def irfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.irfft2(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                       workers=workers, plan=plan)
    return xp.asarray(y)


def rfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of rfftn.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.rfftn(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                             workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.rfftn(x, s, axes, norm)
    x = np.asarray(x)
    y = pfft.rfftn(x, s, axes, norm)
    return xp.asarray(y)


def irfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of irfftn.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.irfftn(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                              workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.irfftn(x, s, axes, norm)
    x = np.asarray(x)
    y = pfft.irfftn(x, s, axes, norm)
    return xp.asarray(y)


def hfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of hfft.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.hfft(x, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                            workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.hfft(x, n, axis, norm)
    x = np.asarray(x)
    y = pfft.hfft(x, n, axis, norm)
    return xp.asarray(y)


def ihfft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    For non-numpy arrays, this implements the Array API specification of ihfft.
    For numpy arrays, see the documentation in _basic_np.py.
    Note that if arguments outside of those in the Array API specification
    are provided with a non-numpy array, an exception is raised.
    """
    if isinstance(x, np.ndarray):
        return npbasic.ihfft(x, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                             workers=workers, plan=plan)
    if overwrite_x is not False:
        Exception
    if workers is not None:
        Exception
    if plan is not None:
        Exception
    xp = array_namespace(x)
    if hasattr(xp, 'fft'):
        return xp.fft.ihfft(x, n, axis, norm)
    x = np.asarray(x)
    y = pfft.ihfft(x, n, axis, norm)
    return xp.asarray(y)


def hfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.hfft2(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                      workers=workers, plan=plan)
    return xp.asarray(y)


def ihfft2(x, s=None, axes=(-2, -1), norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.ihfft2(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                       workers=workers, plan=plan)
    return xp.asarray(y)


def hfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
          plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.hfftn(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                      workers=workers, plan=plan)
    return xp.asarray(y)


def ihfftn(x, s=None, axes=None, norm=None, overwrite_x=False, workers=None, *,
           plan=None):
    """
    See the documentation in _basic_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = npbasic.ihfftn(x, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                       workers=workers, plan=plan)
    return xp.asarray(y)
