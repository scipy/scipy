from scipy._lib._array_api import array_namespace
from . import _realtransforms_np as nprt
import numpy as np

__all__ = ['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn', 'idstn']


def dctn(x, type=2, s=None, axes=None, norm=None, overwrite_x=False,
         workers=None, *, orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.dctn(x, type=type, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                  workers=workers, orthogonalize=orthogonalize)
    return xp.asarray(y)


def idctn(x, type=2, s=None, axes=None, norm=None, overwrite_x=False,
          workers=None, orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.idctn(x, type=type, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                   workers=workers, orthogonalize=orthogonalize)
    return xp.asarray(y)


def dstn(x, type=2, s=None, axes=None, norm=None, overwrite_x=False,
         workers=None, orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.dstn(x, type=type, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                  workers=workers, orthogonalize=orthogonalize)
    return xp.asarray(y)


def idstn(x, type=2, s=None, axes=None, norm=None, overwrite_x=False,
          workers=None, orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.idstn(x, type=type, s=s, axes=axes, norm=norm, overwrite_x=overwrite_x,
                   workers=workers, orthogonalize=orthogonalize)
    return xp.asarray(y)


def dct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False, workers=None,
        orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.dct(x, type=type, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                 workers=workers, orthogonalize=orthogonalize)
    return xp.asarray(y)


def idct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False,
         workers=None, orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.idct(x, type=type, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                  workers=workers, orthogonalize=orthogonalize)
    return xp.asarray(y)


def dst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False, workers=None,
        orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.dst(x, type=type, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                 workers=workers, orthogonalize=orthogonalize)
    return xp.asarray(y)


def idst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False,
         workers=None, orthogonalize=None):
    """
    See the documentation in _realtransforms_np.py.
    Returns an array of the same library as the input array (where possible).
    """
    xp = array_namespace(x)
    x = np.asarray(x)
    y = nprt.idst(x, type=type, n=n, axis=axis, norm=norm, overwrite_x=overwrite_x,
                  workers=workers, orthogonalize=orthogonalize)

    return xp.asarray(y)
