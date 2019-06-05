from . import _fftpack

__all__ = ['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn', 'idstn']

def dctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return _fftpack.dctn(x, type, shape, axes, norm, overwrite_x)
dctn.__doc__ = _fftpack.dctn.__doc__.replace('fftpack', 'fft')

def idctn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return _fftpack.idctn(x, type, shape, axes, norm, overwrite_x)
idctn.__doc__ = _fftpack.idctn.__doc__.replace('fftpack', 'fft')

def dstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return _fftpack.dstn(x, type, shape, axes, norm, overwrite_x)
dstn.__doc__ = _fftpack.dstn.__doc__.replace('fftpack', 'fft')

def idstn(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
    return _fftpack.idstn(x, type, shape, axes, norm, overwrite_x)
idstn.__doc__ = _fftpack.idstn.__doc__.replace('fftpack', 'fft')

def dct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return _fftpack.dct(x, type, n, axis, norm, overwrite_x)
dct.__doc__ = _fftpack.dct.__doc__.replace('fftpack', 'fft')

def idct(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return _fftpack.idct(x, type, n, axis, norm, overwrite_x)
idct.__doc__ = _fftpack.idct.__doc__.replace('fftpack', 'fft')

def dst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return _fftpack.dst(x, type, n, axis, norm, overwrite_x)
dst.__doc__ = _fftpack.dst.__doc__.replace('fftpack', 'fft')

def idst(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
    return _fftpack.idst(x, type, n, axis, norm, overwrite_x)
idst.__doc__ = _fftpack.idst.__doc__.replace('fftpack', 'fft')
