
import numpy as np

class ufunc_wrapper(object):
    def __init__(self, ufunc, diffs = (), diff_alls = ()):
        self.ufunc = ufunc
        self.until_diffs = diff_alls
        self._as_ufunc_out = lambda out: out

    def resolve_out_shapes(self, func):
        self._resolve_out_shapes = func

    def as_ufunc_out(self, func):
        self._as_ufunc_out = func

    def __call__(self, *args, diff = 0, diff_all = False):
        if (diff == 0):
            uf = self.ufunc
        elif diff_all:
            uf = self.until_diffs[diff - 1]
        else:
            uf = self.diffs[diff - 1]

        resolve_out_shapes_args = args[:-uf.nin]
        args = args[-uf.nin:]

        arg_shapes = tuple(np.shape(arg) for arg in args)
        out_shapes = self._resolve_out_shapes(*resolve_out_shapes_args, arg_shapes)
        if (diff > 0):
            out_shapes = (diff + 1) * (out_shapes,)

        arg_dtypes = tuple(arg.dtype if hasattr(arg, 'dtype') else type(arg) for arg in args) + uf.nout * (None,)
        dtypes = uf.resolve_dtypes(arg_dtypes)
        out_dtypes = dtypes[-uf.nout:]

        f = self._as_ufunc_out
        if (diff > 0):
            f = lambda out: tuple(self._as_ufunc_out(out[i]) for i in range(len(out)))

        out = tuple(np.empty(out_shape, dtype = out_dtype) for out_shape, out_dtype in zip(out_shapes, out_dtypes))
        uf(*args, out = f(out)) 

        if (len(out) == 1):
            out, = out

        return out
