from scipy._lib._array_api import array_namespace, is_numpy
from . import _pocketfft
import numpy as np


def arg_err_msg(param):
    return f'Providing {param!r} is only supported for numpy arrays.'


def _validate_fft_args(workers, plan, norm):
    if workers is not None:
        raise ValueError(arg_err_msg("workers"))
    if plan is not None:
        raise ValueError(arg_err_msg("plan"))
    if norm is None:
        norm = 'backward'
    return norm


def _execute(func_str, pocketfft_func, x, **kwargs):
    xp = array_namespace(x)
    # pocketfft is used whenever SCIPY_ARRAY_API is not set,
    # or x is a NumPy array or array-like.
    # When SCIPY_ARRAY_API is set, we try to use xp.fft for CuPy arrays,
    # PyTorch arrays and other array API standard supporting objects.
    # If xp.fft does not exist, we attempt to convert to np and back to use pocketfft.
    if is_numpy(xp):
        return pocketfft_func(x, **kwargs)

    try:
        s = kwargs["n"]
    except KeyError:
        s = kwargs["s"]
    try:
        axes = kwargs["axis"]
    except KeyError:
        axes = kwargs["axes"]
    norm = kwargs["norm"]
    workers = kwargs["workers"]
    plan = kwargs["plan"]

    norm = _validate_fft_args(workers, plan, norm)
    if hasattr(xp, 'fft'):
        xp_func = getattr(xp.fft, func_str)
        return xp_func(x, s, axes, norm=norm)

    x = np.asarray(x)
    y = pocketfft_func(x, s, axes, norm=norm)
    return xp.asarray(y)


def fft(x, n=None, axis=-1, norm=None,
        overwrite_x=False, workers=None, *, plan=None):
    return _execute('fft', _pocketfft.fft, x, n=n, axis=axis, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def ifft(x, n=None, axis=-1, norm=None, overwrite_x=False, workers=None, *,
         plan=None):
    return _execute('ifft', _pocketfft.ifft, x, n=n, axis=axis, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def rfft(x, n=None, axis=-1, norm=None,
         overwrite_x=False, workers=None, *, plan=None):
    return _execute('rfft', _pocketfft.rfft, x, n=n, axis=axis, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def irfft(x, n=None, axis=-1, norm=None,
          overwrite_x=False, workers=None, *, plan=None):
    return _execute('irfft', _pocketfft.irfft, x, n=n, axis=axis, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def hfft(x, n=None, axis=-1, norm=None,
         overwrite_x=False, workers=None, *, plan=None):
    return _execute('hfft', _pocketfft.hfft, x, n=n, axis=axis, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def ihfft(x, n=None, axis=-1, norm=None,
          overwrite_x=False, workers=None, *, plan=None):
    return _execute('ihfft', _pocketfft.ihfft, x, n=n, axis=axis, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def fftn(x, s=None, axes=None, norm=None,
         overwrite_x=False, workers=None, *, plan=None):
    return _execute('fftn', _pocketfft.fftn, x, s=s, axes=axes, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)



def ifftn(x, s=None, axes=None, norm=None,
          overwrite_x=False, workers=None, *, plan=None):
    return _execute('ifftn', _pocketfft.ifftn, x, s=s, axes=axes, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def fft2(x, s=None, axes=(-2, -1), norm=None,
         overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.fft2(x, s=s, axes=axes, norm=norm,
                        overwrite_x=overwrite_x,
                        workers=workers, plan=plan)
    return xp.asarray(y)


def ifft2(x, s=None, axes=(-2, -1), norm=None,
          overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.ifft2(x, s=s, axes=axes, norm=norm,
                         overwrite_x=overwrite_x,
                         workers=workers, plan=plan)
    return xp.asarray(y)


def rfftn(x, s=None, axes=None, norm=None,
          overwrite_x=False, workers=None, *, plan=None):
    return _execute('rfftn', _pocketfft.rfftn, x, s=s, axes=axes, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def rfft2(x, s=None, axes=(-2, -1), norm=None,
         overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.rfft2(x, s=s, axes=axes, norm=norm,
                         overwrite_x=overwrite_x,
                         workers=workers, plan=plan)
    return xp.asarray(y)


def irfftn(x, s=None, axes=None, norm=None,
           overwrite_x=False, workers=None, *, plan=None):
    return _execute('irfftn', _pocketfft.irfftn, x, s=s, axes=axes, norm=norm,
                    overwrite_x=overwrite_x, workers=workers, plan=plan)


def irfft2(x, s=None, axes=(-2, -1), norm=None,
           overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.irfft2(x, s=s, axes=axes, norm=norm,
                          overwrite_x=overwrite_x,
                          workers=workers, plan=plan)
    return xp.asarray(y)


def hfftn(x, s=None, axes=None, norm=None,
          overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.hfftn(x, s=s, axes=axes, norm=norm,
                         overwrite_x=overwrite_x,
                         workers=workers, plan=plan)
    return xp.asarray(y)


def hfft2(x, s=None, axes=(-2, -1), norm=None,
          overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.hfft2(x, s=s, axes=axes, norm=norm,
                         overwrite_x=overwrite_x,
                         workers=workers, plan=plan)
    return xp.asarray(y)


def ihfftn(x, s=None, axes=None, norm=None,
           overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.ihfftn(x, s=s, axes=axes, norm=norm,
                          overwrite_x=overwrite_x,
                          workers=workers, plan=plan)
    return xp.asarray(y)


def ihfft2(x, s=None, axes=(-2, -1), norm=None,
           overwrite_x=False, workers=None, *, plan=None):
    xp = array_namespace(x)
    x = np.asarray(x)
    y = _pocketfft.ihfft2(x, s=s, axes=axes, norm=norm,
                          overwrite_x=overwrite_x,
                          workers=workers, plan=plan)
    return xp.asarray(y)
