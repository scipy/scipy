
import numpy as np

class ufunc_wrapper(object):
    def __init__(self, ufunc, diffs = (), until_diffs = ()):
        if until_diffs:
            until_diffs = tuple(ufunc_wrapper(until_diff) for until_diff in until_diffs)

        self.ufunc = ufunc
        self.until_diffs = until_diffs
        self._as_ufunc_out = lambda out: out

    @property
    def until_jac(self):
        return self.until_diffs[0]

    def resolve_out_shapes(self, func):
        for k, until_diff in enumerate(self.until_diffs, 1):
            until_diff.resolve_out_shapes(lambda *args: (k + 1) * (func(*args),))

        self._resolve_out_shapes = func

    def as_ufunc_out(self, func):
        for k, until_diff in enumerate(self.until_diffs, 1):
            until_diff.as_ufunc_out(lambda out: tuple(func(out[i]) for i in range(len(out))))

        self._as_ufunc_out = func

    def __call__(self, *args):
        resolve_out_shapes_args = args[:-self.ufunc.nin]
        args = args[-self.ufunc.nin:]

        arg_shapes = tuple(np.shape(arg) for arg in args)
        out_shapes = self._resolve_out_shapes(*resolve_out_shapes_args, arg_shapes)

        arg_dtypes = tuple(arg.dtype if hasattr(arg, 'dtype') else type(arg) for arg in args) + self.ufunc.nout * (None,)
        dtypes = self.ufunc.resolve_dtypes(arg_dtypes)
        out_dtypes = dtypes[-self.ufunc.nout:]

        out = tuple(np.empty(out_shape, dtype = out_dtype) for out_shape, out_dtype in zip(out_shapes, out_dtypes))
        self.ufunc(*args, out = self._as_ufunc_out(out)) 

        if (len(out) == 1):
            out, = out

        return out
