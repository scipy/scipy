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

        ufunc_args = args[-ufunc.nin:] # array arguments to be passed to the ufunc

        ufunc_arg_shapes = tuple(np.shape(ufunc_arg) for ufunc_arg in ufunc_args)
        ufunc_out_shapes = self._resolve_out_shapes(*args[:-ufunc.nin], *ufunc_arg_shapes,
            ufunc.nout)

        ufunc_arg_dtypes = tuple(arg.dtype if hasattr(arg, 'dtype') else np.dtype(type(arg)) for arg in ufunc_args)
        if hasattr(ufunc, 'resolve_dtypes'):
            ufunc_dtypes = ufunc_arg_dtypes + ufunc.nout * (None,)
            ufunc_dtypes = ufunc.resolve_dtypes(ufunc_dtypes) 
            ufunc_out_dtypes = ufunc_dtypes[-ufunc.nout:]
        else:
            ufunc_out_dtype = np.result_type(*ufunc_arg_dtypes)
            if (not np.issubdtype(ufunc_out_dtype, np.inexact)):
                ufunc_out_dtype = np.float64

            ufunc_out_dtypes = ufunc.nout * (ufunc_out_dtype,)

        if self.force_out_complex:
            ufunc_out_dtypes = tuple(np.result_type(1j, ufunc_out_dtype) for ufunc_out_dtype in ufunc_out_dtypes)

        b = np.broadcast(*ufunc_args)
        ufunc_out_new_dims = tuple(len(out_shape) - b.ndim for out_shape in ufunc_out_shapes)

        out = tuple(np.empty(ufunc_out_shape, dtype = ufunc_out_dtype)
            for ufunc_out_shape, ufunc_out_dtype in zip(ufunc_out_shapes, ufunc_out_dtypes))

        ufunc_out = tuple(np.moveaxis(out, tuple(range(new_dim)), tuple(range(-new_dim, 0)))
            for out, new_dim in zip(out, ufunc_out_new_dims))
        ufunc(*ufunc_args, out = ufunc_out)

        if (len(out) == 1):
            out, = out

        return out
