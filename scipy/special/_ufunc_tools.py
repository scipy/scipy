"""Helpers for producing efficient wrappers of ufuncs.
"""

import re
import numpy as np


def _parse_core_ndims(signature):
    """Return tuple of num core dims per input from gufunc signature."""
    input_sig = signature.split('->')[0]
    groups = re.findall(r"\((.*?)\)", input_sig)
    return tuple(0 if not g.strip() else g.count(',') + 1 for g in groups)


def _resolve_alloc_order(args, order):
    """Determine contiguity of output when using in-ufunc caching.

    Since we hijack the iteration order, 'K' isn't really possible
    so it's just treated like 'A'.

    """
    order = order.upper() if order is None else "K"
    if order in ('K', 'A'):
        # If all inputs are F-contiguous, return F-contiguous.
        # Otherwise, default to C.
        if all(arg.flags.f_contiguous for arg in args):
            return 'F'
        return 'C'
    return order  # Returns 'C' or 'F'


def _with_cache_optimization(
        *,
        name,
        arg_names,
        docstring,
        ufunc,
        cache_arg_indices,
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

    Returns
    -------
    callable
        A wrapper for ufunc which transposes the axes of the inputs to ensure
        iteration precedes in such a way to allow the cache within the ufunc
        kernel to eliminate redundant computation.

    Notes
    -----
    There is a common pattern in ufunc kernels exemplified by the situation
    where some of the arguments are used to compute coeffients of an expansion
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

    Note that because the pre-allocated output array used internally is C
    contiguous, the output will be C contiguous regardless of contiguity of
    the inputs.

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

    def _wrapper(
            *args, out=None, where=True, casting="same_kind", order="K", dtype=None
    ):
        args = [np.asarray(arg) for arg in args]
        kwargs = dict(casting=casting, dtype=dtype)
        if out is not None:
            kwargs["out"] = out
            if is_elementwise:
                kwargs["where"] = where

        # Fast path for when the arguments which are used in the cached
        # computation don't have batches.
        if all(args[i].ndim == core_ndims[i] for i in cache_arg_indices):
            kwargs["order"] = order
            return ufunc(*args, **kwargs)

        # To get batch_shapes, need to exclude core dimensions. Again, the core
        # dimensions won't participate in broadcasting.
        batch_shapes = [
            arg.shape[:-core_ndims[i]] if core_ndims[i] > 0 else arg.shape
            for i, arg in enumerate(args)
        ]
        batch_shape = np.broadcast_shapes(*batch_shapes)
        batch_ndim = len(batch_shape)

        # Broadcast each arg so that the batch shapes all agree, but the
        # core dimensions may still differ.
        args_b = [
            np.broadcast_to(arg, batch_shape + arg.shape[-core_ndims[i]:])
            if core_ndims[i] > 0 else np.broadcast_to(arg, batch_shape)
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
                np.broadcast_to(where, batch_shape), axes=sorted_batch_axes
            )

        # Handle output array (use provided 'out' or pre-allocate)
        if out is not None:
            out_final = out
            if ufunc.nout == 1:
                # a view of the provided output array with axes sorted as needed.
                out_t = np.transpose(out_final, axes=sorted_batch_axes)
            else:
                out_t = tuple(
                    np.transpose(x, axes=sorted_batch_axes) for x in out_final
                )
        else:
            # Respect users choice of order, and preserve F-contiguity is all inputs
            # are F-contiguous and order="K" or "A". This departs from the typical
            # meaning of order = "K", which enables NumPy's clever logic to determine
            # a cache friendly iteration order. This is because we hijack the iteration
            # order for the benefit of in-ufunc caches.
            alloc_order = _resolve_alloc_order(args, order)
            if dtype is not None:
                out_dtype = np.dtype(dtype)
            else:
                input_dtypes = tuple(arg.dtype for arg in args_t) + (None,)*ufunc.nout
                out_dtype = ufunc.resolve_dtypes(input_dtypes, casting=casting)[-1]
            if ufunc.nout == 1:
                out_final = np.empty(batch_shape, dtype=out_dtype, order=alloc_order)
                # a view of the output array with axes sorted as needed.
                out_t = np.transpose(out_final, axes=sorted_batch_axes)
            else:
                out_final = tuple(
                    np.empty(batch_shape, dtype=out_dtype, order=alloc_order)
                    for _ in range(ufunc.nout)
                )
                out_t = tuple(
                    np.transpose(x, axes=sorted_batch_axes) for x in out_final
                )
        kwargs["out"] = out_t
        kwargs["order"] = "C"

        # Set out to the above view, but return the C contiguous output.
        # This avoids having non-contiguous output.
        ufunc(*args_t, **kwargs)

        return out_final

    # Do some metaprogramming with exec so that the arg names and func
    # name are as expected.
    arg_str = ", ".join(arg_names)

    # Handle kwargs for function definition and call to wrapper.
    kwarg_defs = ["out=None", "casting='same_kind'", "order='K'", "dtype=None"]
    kwarg_calls = ["out=out", "casting=casting", "order=order", "dtype=dtype"]
    if is_elementwise:
        # where is only available for elementwise ufuncs, not gufuncs.
        kwarg_defs.insert(1, "where=True")
        kwarg_calls.insert(1, "where=where")
    signature_kwargs = ", ".join(kwarg_defs)
    call_kwargs = ", ".join(kwarg_calls)
    code = (
        f"""def {name}({arg_str}, {signature_kwargs}):
            return _wrapper({arg_str}, {call_kwargs})
        """
        )
    namespace = {"_wrapper": _wrapper}
    exec(code, namespace)
    wrapper = namespace[name]
    wrapper.__doc__ = docstring
    return wrapper
