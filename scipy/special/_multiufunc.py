import numpy as np

class multiufunc:
    def __init__(self, ufuncs, resolve_ufunc = None, resolve_out_shapes = None, force_out_complex = False):
        self.ufuncs = ufuncs
        self._resolve_out_shapes = resolve_out_shapes
        self._resolve_ufunc = resolve_ufunc
        self.force_out_complex = force_out_complex

    def resolve_ufunc(self, func):
        self._resolve_ufunc = func

    def resolve_out_shapes(self, func):
        self._resolve_out_shapes = func

    def __call__(self, *args, **kwargs):
        ufunc = self._resolve_ufunc(self.ufuncs, **kwargs)
        if ((ufunc.nout == 0) or (self._resolve_out_shapes is None)):
            return ufunc(*args)

        resolve_out_shapes_args = args[:-ufunc.nin]
        args = args[-ufunc.nin:]

        arg_shapes = tuple(np.shape(arg) for arg in args)
        out_shapes = self._resolve_out_shapes(*resolve_out_shapes_args, *arg_shapes, ufunc.nout)

        arg_dtypes = tuple(arg.dtype if hasattr(arg, 'dtype') else np.dtype(type(arg)) for arg in args)
        if hasattr(ufunc, 'resolve_dtypes'):
            dtypes = arg_dtypes + ufunc.nout * (None,)
            dtypes = ufunc.resolve_dtypes(dtypes) 
            out_dtypes = dtypes[-ufunc.nout:]
        else:
            out_dtype = np.result_type(*arg_dtypes)
            if (not np.issubdtype(out_dtype, np.inexact)):
                out_dtype = np.float64

            out_dtypes = ufunc.nout * (out_dtype,)

        if self.force_out_complex:
            out_dtypes = tuple(np.result_type(1j, out_dtype) for out_dtype in out_dtypes)

        b_shape = np.broadcast_shapes(*arg_shapes)

        out = tuple(np.empty(out_shape, dtype = out_dtype) for out_shape, out_dtype in zip(out_shapes, out_dtypes))

        new_dims = tuple(len(out_shape) - len(b_shape) for out_shape in out_shapes)

        out_src_axes = tuple(range(new_dim))
        out_dst_axes = tuple(range(-new_dim, 0))

        out_ufunc = tuple(np.moveaxis(out, out_src_axes, out_dst_axes) for out, new_dim in zip(out, new_dims))
        ufunc(*args, out = out_ufunc)

        if (len(out) == 1):
            out, = out

        return out
