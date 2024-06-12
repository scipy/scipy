import collections
import inspect
import numpy as np

class MultiUFunc:
    __slots__ = ["_ufuncs_map", "resolve_out_shapes", "resolve_ufunc",
                 "__forces_complex_output", "_ufuncs"]
    def __init__(self, ufuncs_map, *, resolve_ufunc=None,
                 resolve_out_shapes=None, force_complex_output=False):

        # Gather leaf level ufuncs from ufuncs_map.
        ufuncs = []
        def traverse(obj):
            if isinstance(obj, collections.abc.Mapping):
                for value in obj.values():
                    traverse(value)
            elif (isinstance(obj, collections.abc.Iterable)
                and not isinstance(obj, (str, bytes))):
                for item in obj:
                    traverse(item)
            else:
                ufuncs.append(obj)
        traverse(ufuncs_map)

        # Perform input validation to ensure all ufuncs in ufuncs_map are
        # actually ufuncs, have distinct names, and all take the same input
        # types.
        seen_names = set()
        seen_input_types = set()
        for ufunc in ufuncs:
            if not isinstance(ufunc, np.ufunc):
                raise ValueError("All leaf elements of ufuncs_map must have"
                                 f" type `numpy.ufunc`. Received {ufuncs_map}")

            if ufunc.__name__ in seen_names:
                raise ValueError("ufuncs within ufuncs_map must all have"
                                 f" distinct names. Received {ufuncs_map}")
            seen_names.add(ufunc.__name__)
            seen_input_types.add(frozenset(x.split("->")[0] for x in ufunc.types))
        if len(seen_input_types) > 1:
            raise ValueError("All ufuncs in ufuncs_map must take the same"
                             " input types.")

        self._ufuncs = tuple(ufuncs)
        self._ufuncs_map = ufuncs_map
        self.resolve_out_shapes = resolve_out_shapes
        self.resolve_ufunc = resolve_ufunc
        self.__forces_complex_output = force_complex_output

    @property
    def forces_complex_output(self):
        return self.__forces_complex_output

    def as_resolve_ufunc(self, func):
        """Set `resolve_ufunc` method by decorating a function.

        The decorated function's first argument should be a JSON-like
        `ufuncs_map`, and additional arguments should be keywords which 
        are used to dispatch to the correct ufunc at the leaf level of
        `ufuncs_map`.
        """
        # Given func, we construct a wrapper which no longer takes the
        # `ufuncs_map` as input, but instead gets it from the
        # class. Use inspect to add an informative signature. Add a
        # a docstring if none exists.
        sig = inspect.signature(func)
        params = list(inspect.signature(func).parameters.values())[1:]
        new_sig = sig.replace(parameters=params)

        def resolve_ufunc(**kwargs):
            return func(self._ufuncs_map, **kwargs)

        resolve_ufunc.__signature__ = new_sig
        docstring = func.__doc__
        if docstring is not None:
            resolve_ufunc.__doc__ = docstring
        resolve_ufunc.__doc__ = \
            """Resolve to a ufunc based on keyword arguments."""
        self.resolve_ufunc = resolve_ufunc

    def as_resolve_out_shapes(self, func):
        """Set `resolve_out_shapes` method by decorating a function."""
        if func.__doc__ is None:
            func.__doc__ = \
                """Resolve to output shapes based on relevant inputs."""
        func.__name__ = "resolve_out_shapes"
        self.resolve_out_shapes = func

    def __call__(self, *args, **kwargs):
        ufunc = self.resolve_ufunc(**kwargs)
        if ((ufunc.nout == 0) or (self.resolve_out_shapes is None)):
            return ufunc(*args)

        ufunc_args = args[-ufunc.nin:] # array arguments to be passed to the ufunc

        ufunc_arg_shapes = tuple(np.shape(ufunc_arg) for ufunc_arg in ufunc_args)
        ufunc_out_shapes = self.resolve_out_shapes(*args[:-ufunc.nin],
            *ufunc_arg_shapes, ufunc.nout)

        ufunc_arg_dtypes = tuple((ufunc_arg.dtype if hasattr(ufunc_arg, 'dtype')
            else np.dtype(type(ufunc_arg))) for ufunc_arg in ufunc_args)
        if hasattr(ufunc, 'resolve_dtypes'):
            ufunc_dtypes = ufunc_arg_dtypes + ufunc.nout * (None,)
            ufunc_dtypes = ufunc.resolve_dtypes(ufunc_dtypes) 
            ufunc_out_dtypes = ufunc_dtypes[-ufunc.nout:]
        else:
            ufunc_out_dtype = np.result_type(*ufunc_arg_dtypes)
            if (not np.issubdtype(ufunc_out_dtype, np.inexact)):
                ufunc_out_dtype = np.float64

            ufunc_out_dtypes = ufunc.nout * (ufunc_out_dtype,)

        if self.forces_complex_output:
            ufunc_out_dtypes = tuple(np.result_type(1j, ufunc_out_dtype)
                for ufunc_out_dtype in ufunc_out_dtypes)

        b = np.broadcast(*ufunc_args)
        ufunc_out_new_dims = tuple(len(ufunc_out_shape) - b.ndim
            for ufunc_out_shape in ufunc_out_shapes)

        out = tuple(np.empty(ufunc_out_shape, dtype = ufunc_out_dtype)
            for ufunc_out_shape, ufunc_out_dtype
            in zip(ufunc_out_shapes, ufunc_out_dtypes))

        ufunc_out = tuple(np.moveaxis(out[i],
            tuple(range(axis)), tuple(range(-axis, 0)))
            for i, axis in enumerate(ufunc_out_new_dims))
        ufunc(*ufunc_args, out = ufunc_out)

        if (len(out) == 1):
            out, = out

        return out
