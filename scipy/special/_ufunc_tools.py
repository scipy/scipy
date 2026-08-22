"""Helpers for producing efficient wrappers of ufuncs.
"""

import re
import numpy as np
import warnings

from scipy._external.packaging_version import version


def _parse_core_ndims(signature):
    """Return tuple of num core dims per input from gufunc signature."""
    input_sig = signature.split('->')[0]
    groups = re.findall(r"\((.*?)\)", input_sig)
    return tuple(0 if not g.strip() else g.count(',') + 1 for g in groups)


def _resolve_alloc_order(args, order):
    """Determine contiguity of output when using in-ufunc caching.

    We hijack the ufunc iteration order, so NumPy's own layout logic can't be
    reused; this reproduces its rules for 'A' and 'K' from the inputs directly.

    """
    order = order.upper() if order is not None else "K"
    if order in ('K', 'A'):
        if order == 'A':
            # NumPy: F if every input is F-contiguous, C otherwise.
            if all(arg.flags.f_contiguous for arg in args):
                return 'F'
            return 'C'
        if order == 'K':
            # NumPy: F only if every input is genuinely Fortran ordered. Arrays
            # which are both C and F contiguous (1-D, or shapes with a single
            # non-unit dimension) don't by themselves make the output F ordered.
            if (
                    all(arg.flags.f_contiguous for arg in args)
                    and any(not arg.flags.c_contiguous for arg in args)
            ):
                return 'F'
            return 'C'
    return order  # Returns 'C' or 'F'


def _allocate_out(args, shape, dtype, order, subok):
    """Allocate output array, respecting subclasses and priority."""
    prototype = None
    if subok and args:
        prototype = max(args, key=lambda x: getattr(x, "__array_priority__", 0.0))

    if prototype is not None:
        return np.empty_like(
            prototype, shape=shape, dtype=dtype, order=order, subok=subok
        )

    return np.empty(shape, dtype=dtype, order=order)


def _normalize_out(out, nout):
    """Validate ``out`` and return it as a tuple of length ``nout``.

    Entries may be None, meaning that output should be allocated for the user,
    matching the behavior of a plain ufunc.

    """
    if out is None:
        return (None,)*nout
    if not isinstance(out, tuple):
        if nout > 1:
            raise TypeError("'out' must be a tuple of arrays")
        out_tuple = (out,)
    else:
        out_tuple = out
    if len(out_tuple) != nout:
        raise ValueError(
            f"The 'out' tuple must have exactly {nout} entries: one per "
            f"ufunc output"
        )
    for entry in out_tuple:
        if entry is not None and not isinstance(entry, np.ndarray):
            raise TypeError("return arrays must be of ArrayType")
    return out_tuple


def _with_cache_optimization(
        *,
        name,
        arg_names,
        docstring,
        ufunc,
        cache_arg_indices,
        module="scipy.special",
):
    """Helper to ensure optimal iteration order for ufuncs that use caching.

    This concerns internal caches which are only live over the course of one
    call to a ufunc to avoid repeated computations during the course of the
    loops. See the notes below for more information.

    Parameters
    ----------
    name : str
    arg_names : list[str]
        The function name and arg names are passed in so that the wrapper
        can be generated with the create name and argument names, improving
        documentation and autocomplete.
    docstring : str
    ufunc : numpy.ufunc
    cache_arg_indices : list[int]
       Arguments to ufunc which are used in the kernel to compute an output
       which is being cached for reuse when iterating over other arguments.
    module : str, optional
       Value to use for the ``__module__`` attribute of the wrapper.

    Returns
    -------
    callable
        A wrapper for ufunc which transposes the axes of the inputs to ensure
        iteration precedes in such a way to allow the cache within the ufunc
        kernel to eliminate redundant computation.
        The underlying ufunc can be accessed through the ``ufunc`` attribute
        of the returned callable.

    Notes
    -----
    There is a common pattern in ufunc kernels exemplified by the situation
    where some of the arguments are used to compute coefficients of an expansion
    that is taken over one or more of the other arguments. A classic example is
    Mathieu functions, which compute coefficients corresponding to the shape
    parameter q and order m which in principle could be reused for varying
    values of the parameter x.

    It had long been the case that the expensive computation of coefficients is
    repeated unnecessarily across values of x. It is possible to add a cache to
    the ufunc kernel which stores the expansion coefficients and only updates
    if the pointers into the q and m arrays advance during the course of the
    ufunc loops. Such a cache is instantiated each time a ufunc is called and
    only lives during the course of the loops that are carried out for that
    particular call.

    Whether such a cache actually helps depends on the order in which iteration
    occurs. Ideally, one would want q and m to advance most slowly and for the
    iterations over x for fixed q and m to be pushed to the inner most loops.
    This helper replaces each input array (and a pre-allocated output array)
    with a view where the axes which should vary most slowly are transposed to
    the ends and forces computation in C order. This ensures iteration proceeds
    in the optimal order.

    """

    # Need to keep track of the number of core dimensions for each input
    # since core dimensions don't participate in broadcasting.
    core_ndims = (
        _parse_core_ndims(ufunc.signature)
        if ufunc.signature is not None
        else (0,)*ufunc.nin
    )

    # The kwarg ``where`` is only supported for elementwise ufuncs, while
    # the kwarg ``axes` is only supported for non-elementwise gufuncs, so
    # keep track of this here for later.
    is_elementwise = sum(core_ndims) == 0

    def wrapper(
            *args,
            out=None,
            where=True,
            casting="same_kind",
            order="K",
            dtype=None,
            subok=True,
            signature=None,
    ):
        asarray = np.asanyarray if subok else np.asarray
        args = [asarray(arg) for arg in args]
        kwargs = dict(casting=casting, subok=subok)

        # ``where`` is normalized once, up front, so that the fast path and the
        # transposed path below agree on its meaning. ``True`` is the "no mask"
        # sentinel and is not forwarded, matching a plain ufunc call.
        if is_elementwise and where is not True:
            where = asarray(where)
            kwargs["where"] = where
        if out is not None:
            kwargs["out"] = out
        if signature is not None:
            kwargs["signature"] = signature
        if dtype is not None:
            kwargs["dtype"] = dtype

        # Fast path for when the arguments which are used in the cached
        # computation don't have batches.
        if all(args[i].ndim == core_ndims[i] for i in cache_arg_indices):
            kwargs["order"] = order
            return ufunc(*args, **kwargs)

        out_tuple = _normalize_out(out, ufunc.nout)

        if (
                out is None
                and is_elementwise
                and where is not True
                and version.parse(np.__version__) >= version.parse("2.4.0")
        ):
            # The transposed path always passes ``out`` to the ufunc, so NumPy
            # can't issue this warning itself. Match its wording.
            warnings.warn(
                "'where' used without 'out', expect unitialized memory in "
                "output. If this is intentional, use out=None.",
                UserWarning,
                stacklevel=3,
            )

        # To get batch_shapes, need to exclude core dimensions. Again, the core
        # dimensions won't participate in broadcasting.
        batch_shapes = [
            arg.shape[:-core_ndims[i]] if core_ndims[i] > 0 else arg.shape
            for i, arg in enumerate(args)
        ]
        if is_elementwise and where is not True:
            batch_shapes.append(where.shape)

        batch_shapes.extend(
            entry.shape for entry in out_tuple if entry is not None
        )

        batch_shape = np.broadcast_shapes(*batch_shapes)
        batch_ndim = len(batch_shape)

        # Broadcast each arg so that the batch shapes all agree, but the
        # core dimensions may still differ.
        args_b = [
            np.broadcast_to(arg, batch_shape + arg.shape[-core_ndims[i]:], subok=subok)
            if core_ndims[i] > 0 else np.broadcast_to(arg, batch_shape, subok=subok)
            for i, arg in enumerate(args)
        ]

        # After broadcasting, determine which axes have stride-length
        # zero for each of the args participating in the cache. The cached
        # value will remain the same for iterations across these axes.
        varying_axes = []
        constant_axes = []

        for ax in range(batch_ndim):
            is_constant = all(
                args_b[i].strides[ax] == 0 for i in cache_arg_indices
            )
            if is_constant:
                constant_axes.append(ax)
            else:
                varying_axes.append(ax)
        sorted_batch_axes = varying_axes + constant_axes

        # Push the non-varying axes to the end so that they will vary most
        # slowly when iterating in C order.
        args_t = []
        for i, arg_b in enumerate(args_b):
            axes = sorted_batch_axes + list(
                range(batch_ndim, batch_ndim + core_ndims[i])
            )
            args_t.append(np.transpose(arg_b, axes=axes))

        # Handle ``where`` kwarg.
        if is_elementwise and where is not True:
            kwargs["where"] = np.transpose(
                np.broadcast_to(where, batch_shape, subok=subok),
                axes=sorted_batch_axes,
            )

        # Allocate any outputs the user didn't supply. A plain ufunc allocates
        # for each None entry of ``out`` too.
        if any(entry is None for entry in out_tuple):
            # Respect the user's choice of order. _resolve_alloc_order
            # reproduces NumPy's layout rules from the untransposed inputs,
            # since we hijack the iteration order for the benefit of in-ufunc
            # caches and so can't let NumPy pick the layout itself.
            alloc_order = _resolve_alloc_order(args, order)
            if dtype is not None:
                out_dtypes = (np.dtype(dtype),)*ufunc.nout
            else:
                resolve_kwargs = {"casting": casting}
                if signature is not None:
                    resolve_kwargs["signature"] = signature
                given_dtypes = tuple(
                    None if entry is None else entry.dtype for entry in out_tuple
                )
                input_dtypes = tuple(arg.dtype for arg in args_t) + given_dtypes
                out_dtypes = ufunc.resolve_dtypes(
                    input_dtypes, **resolve_kwargs
                )[ufunc.nin:]
            outputs = tuple(
                _allocate_out(args, batch_shape, out_dtype, alloc_order, subok)
                if entry is None else entry
                for entry, out_dtype in zip(out_tuple, out_dtypes)
            )
        else:
            outputs = out_tuple

        # Views of the output arrays with axes sorted as needed.
        kwargs["out"] = tuple(
            np.transpose(x, axes=sorted_batch_axes) for x in outputs
        )
        kwargs["order"] = "C"

        # Set out to the above views, but return the untransposed outputs. This
        # avoids having non-contiguous output.
        ufunc(*args_t, **kwargs)

        return outputs[0] if ufunc.nout == 1 else outputs

    wrapper._ufunc = ufunc
    wrapper.__qualname__ = name
    return _make_ufunc_like_wrapper(wrapper, name, arg_names, docstring, module=module)


def _reconstruct_wrapper(name):
    """Helper to allow pickling of dynamically generated ufunc wrappers."""
    import scipy.special
    return getattr(scipy.special, name)


class _UFuncLikeWrapper:
    """Base class for ufunc wrappers to provide IDE-friendly properties."""

    def __init__(self, attributes):
        self._attributes = attributes

    @property
    def nin(self):
        """The number of inputs."""
        return self._attributes["nin"]

    @property
    def nout(self):
        """The number of outputs."""
        return self._attributes["nout"]

    @property
    def nargs(self):
        """The number of arguments."""
        return self._attributes["nargs"]

    @property
    def ntypes(self):
        """The number of types."""
        return self._attributes["ntypes"]

    @property
    def types(self):
        """A list with types grouped input->output."""
        return self._attributes["types"]

    @property
    def signature(self):
        """Definition of the core elements a generalized ufunc operates on."""
        return self._attributes["signature"]

    def resolve_dtypes(self, dtypes, *, signature=None, casting=None, reduction=False):
        # The underlying wrapper/ufunc logic handles the actual type resolution
        kwargs = {"reduction": reduction}
        # although the defaults for signature and casting are ``None``,
        # one cannot actually pass these kwargs with ``None`` values.
        if casting is not None:
            kwargs["casting"] = casting
        if signature is not None:
            kwargs["signature"] = signature
        return self._attributes["resolve_dtypes"](dtypes, **kwargs)

    def __reduce__(self):
        # Tells pickle exactly how to reconstruct this specific instance
        return (_reconstruct_wrapper, (self.__name__,))


_UFuncLikeWrapper.resolve_dtypes.__doc__ = np.ufunc.resolve_dtypes.__doc__


def _make_ufunc_like_wrapper(
        func,
        name,
        arg_names,
        docstring,
        ufunc_attributes=None,
        module="scipy.special",
):
    """Wrapper for a ufunc that preserves ufunc-like behavior.

    Parameters
    ----------
    func : callable
        The internal function to wrap. Must either have a ``_ufunc``
        attribute to infer metadata from, or all metadata arguments must be
        explicitly supplied.
    name : str
        Name of the function.
    arg_names : list[str]
        Names of the arguments. May include positional default values
        (e.g., ``["z", "k=0"]``).
    docstring : str
        The docstring for the wrapper.
    ufunc_attributes : dict, optional
        nin : int
            The number of inputs.
        nout : int
            The number of outputs.
        nargs : int
            The number of arguments.
        ntypes : int
            The number of types.
        types: list[str]
            A list with types grouped input->output.
        signature: str or None
            Definition of the core elements a generalized ufunc operates on.
            If ``None``, the function is treated as an elementwise ufunc.
        resolve_dtypes: callable
            A function determining the output dtypes based on input dtypes.
    module : str, optional
       Value to use for the ``__module__`` attribute of the wrapper.

    Returns
    -------
    callable
        A wrapper for func that maintains ufunc kwargs and attributes.
    """
    metadata_source = getattr(func, "_ufunc", None)
    attributes = None
    if metadata_source is not None:
        attributes = {
            "nin": metadata_source.nin,
            "nout": metadata_source.nout,
            "nargs": metadata_source.nargs,
            "ntypes": metadata_source.ntypes,
            "types": metadata_source.types,
            "signature": metadata_source.signature,
            "resolve_dtypes": metadata_source.resolve_dtypes
        }
    else:
        attributes = ufunc_attributes
    if attributes is None:
        raise ValueError("func has no metadata and metadata was not supplied.")

    arg_str = ", ".join(arg_names)
    clean_args = [arg.split("=")[0].strip() for arg in arg_names]
    call_args = ", ".join(clean_args)

    code = (
        f"""class {name}(_UFuncLikeWrapper):
            def __call__(self, {arg_str}, /, out=None, **kwargs):
                return func({call_args}, out=out, **kwargs)

        """
    )

    namespace = {"_UFuncLikeWrapper": _UFuncLikeWrapper, "func": func}

    exec(code, namespace)
    WrapperClass = namespace[name]
    WrapperClass.__module__ = module
    WrapperClass.__qualname__ = name

    wrapper = WrapperClass(attributes)

    wrapper.__doc__ = docstring
    wrapper.__module__ = module
    wrapper.__name__ = name
    wrapper.__qualname__ = name

    return wrapper
