from scipy._lib._array_api import _has_own_device, array_namespace, xp_device
import numpy as np
from . import _duccfft

__all__ = ['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn', 'idstn']


def _execute(duccfft_func, x, type, s, axes, norm, 
             overwrite_x, workers, orthogonalize):
    xp = array_namespace(x)
    # the NumPy round-trip must return the result on the input's device
    device = xp_device(x) if _has_own_device(x) else None
    x = np.asarray(x)
    y = duccfft_func(x, type, s, axes, norm,
                       overwrite_x=overwrite_x, workers=workers,
                       orthogonalize=orthogonalize)
    return xp.asarray(y, device=device)


def dctn(x, type=2, s=None, axes=None, norm=None,
         overwrite_x=False, workers=None, *, orthogonalize=None):
    return _execute(_duccfft.dctn, x, type, s, axes, norm, 
                    overwrite_x, workers, orthogonalize)


def idctn(x, type=2, s=None, axes=None, norm=None,
          overwrite_x=False, workers=None, *, orthogonalize=None):
    return _execute(_duccfft.idctn, x, type, s, axes, norm, 
                    overwrite_x, workers, orthogonalize)


def dstn(x, type=2, s=None, axes=None, norm=None,
         overwrite_x=False, workers=None, orthogonalize=None):
    return _execute(_duccfft.dstn, x, type, s, axes, norm, 
                    overwrite_x, workers, orthogonalize)


def idstn(x, type=2, s=None, axes=None, norm=None,
          overwrite_x=False, workers=None, *, orthogonalize=None):
    return _execute(_duccfft.idstn, x, type, s, axes, norm, 
                    overwrite_x, workers, orthogonalize)


def dct(x, type=2, n=None, axis=-1, norm=None,
        overwrite_x=False, workers=None, orthogonalize=None):
    return _execute(_duccfft.dct, x, type, n, axis, norm, 
                    overwrite_x, workers, orthogonalize)


def idct(x, type=2, n=None, axis=-1, norm=None,
         overwrite_x=False, workers=None, orthogonalize=None):
    return _execute(_duccfft.idct, x, type, n, axis, norm, 
                    overwrite_x, workers, orthogonalize)


def dst(x, type=2, n=None, axis=-1, norm=None,
        overwrite_x=False, workers=None, orthogonalize=None):
    return _execute(_duccfft.dst, x, type, n, axis, norm, 
                    overwrite_x, workers, orthogonalize)


def idst(x, type=2, n=None, axis=-1, norm=None,
         overwrite_x=False, workers=None, orthogonalize=None):
    return _execute(_duccfft.idst, x, type, n, axis, norm, 
                    overwrite_x, workers, orthogonalize)
