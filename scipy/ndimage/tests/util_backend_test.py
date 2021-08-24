import numpy as np
import scipy._lib.uarray as ua
from scipy.ndimage._backend import _ScipyImageBackend
from scipy.ndimage._multimethods import _ndimage_output, _mark_index


class _BackendWithConverter(_ScipyImageBackend):

    @ua.wrap_single_convertor
    def __ua_convert__(value, dispatch_type, coerce):
        if value is None:
            return None

        if dispatch_type is np.ndarray:
            if not coerce and not isinstance(value, np.ndarray):
                return NotImplemented

            return np.asarray(value)

        if dispatch_type is np.dtype:
            try:
                return np.dtype(str(value))
            except TypeError:
                return np.dtype(value)

        if dispatch_type is _ndimage_output:
            if isinstance(value, np.ndarray):
                return np.asarray(value)
            else:
                try:
                    return np.dtype(str(value))
                except TypeError:
                    return np.dtype(value)

        if dispatch_type is _mark_index:
            if isinstance(value, np.ndarray):
                return np.asarray(value)
            else:
                return value

        return value