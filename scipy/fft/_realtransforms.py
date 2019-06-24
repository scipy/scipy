import scipy.fftpack as _fftpack

__all__ = ['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn', 'idstn']

def _doc_wrap(transform_func, new_func):
    doc = transform_func.__doc__ or ''
    new_func.__doc__ = doc.replace('fftpack', 'fft')
    new_func.__name__ = transform_func.__name__
    return new_func


def _doc_wrap_1d(transform_func):
    def inner(x, type=2, n=None, axis=-1, norm=None, overwrite_x=False):
        return transform_func(x, type, n, axis, norm, overwrite_x)
    return _doc_wrap(transform_func, inner)


def _doc_wrap_nd(transform_func):
    def inner(x, type=2, shape=None, axes=None, norm=None, overwrite_x=False):
        return transform_func(x, type, shape, axes, norm, overwrite_x)
    return _doc_wrap(transform_func, inner)


dctn = _doc_wrap_nd(_fftpack.dctn)
idctn = _doc_wrap_nd(_fftpack.idctn)
dstn = _doc_wrap_nd(_fftpack.dstn)
idstn = _doc_wrap_nd(_fftpack.idstn)

dct = _doc_wrap_1d(_fftpack.dct)
idct = _doc_wrap_1d(_fftpack.idct)
dst = _doc_wrap_1d(_fftpack.dst)
idst = _doc_wrap_1d(_fftpack.idst)
