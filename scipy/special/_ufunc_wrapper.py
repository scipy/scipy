
import numpy as np

class ufunc_wrapper(object):
    def __init__(self, ufuncs):
        self.ufuncs = ufuncs

        self._as_ufunc_out = None
        self._resolve_out_shapes = None
        self._resolve_ufunc = None

    def resolve_out_shapes(self, func):
        self._resolve_out_shapes = func

    def resolve_ufunc(self, func):
        self._resolve_ufunc = func

    def as_ufunc_out(self, func):
        self._as_ufunc_out = func

    def __call__(self, *args, **kwargs):
        ufunc = self._resolve_ufunc(self.ufuncs, **kwargs)

        if (self._resolve_out_shapes is None and self._as_ufunc_out is None):
            return ufunc(*args)

        resolve_out_shapes_args = args[:-ufunc.nin]
        args = args[-ufunc.nin:]

        arg_shapes = tuple(np.shape(arg) for arg in args)
        out_shapes = self._resolve_out_shapes(*resolve_out_shapes_args, arg_shapes)
        if (ufunc.nout > 1):
            out_shapes = (ufunc.nout) * (out_shapes,)

        arg_dtypes = tuple(arg.dtype if hasattr(arg, 'dtype') else type(arg) for arg in args) + ufunc.nout * (None,)
        dtypes = ufunc.resolve_dtypes(arg_dtypes)
        out_dtypes = dtypes[-ufunc.nout:]

        if (self._as_ufunc_out is None):
            f = lambda out: out
        else:
            f = self._as_ufunc_out

        if (ufunc.nout >= 1):
            f = lambda out: tuple(self._as_ufunc_out(out[i]) for i in range(len(out)))

        out = tuple(np.empty(out_shape, dtype = out_dtype) for out_shape, out_dtype in zip(out_shapes, out_dtypes))
        ufunc(*args, out = f(out)) 

        if (len(out) == 1):
            out, = out

        return out
