"""Helpers for producing efficient wrappers of ufuncs.
"""

import re
import numpy as np


def _parse_core_ndims(signature: str) -> tuple[int]:
    """Return tuple of num core dims per input from gufunc signature."""
    input_sig = signature.split('->')[0]
    groups = re.findall(r"\((.*?)\)", input_sig)
    return tuple(0 if not g.strip() else g.count(',') + 1 for g in groups)


def _with_cache_optimization(ufunc, cache_arg_indices):
    """Helper to ensure optimal iteration order for ufuncs that use caching

    Parameters
    ----------
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
    over other the arguments. A classic example is mathieu functions, which
    compute coefficients corresponding to the shape parameter q and order m
    which in principle could be reused for varying values of the parameter
    x.

    It has long been the case that the expensive computation of coefficients is
    repeated unnecessarily across values of x. It is possible to add a cache to
    the ufunc kernel in this case which stores the expansion coefficients and
    only updates it if the pointers for the q and m arrays advance, however
    whether this cache actually helps depends on the order in which iteration
    occurs. Ideally, one would want q and m to advance most slowly and for the
    iterations over x for fixed q and m to be pushed to the inner most loops.
    To work around difficulties in controlling the iteration order in ufuncs,
    this helper will transpose the axes which we want to vary most slowly to
    the ends of the arrays and force computation in C order.

    Take note that the output will be C contiguous regardless of contiguity
    of the inputs.

    """

    # Need to keep track of the number of core dimensions for each input
    # since core dimensions don't participate in broadcasting.
    core_ndims = (
        _parse_core_ndims(ufunc.signature)
        if ufunc.signature is not None
        else (0,)*ufunc.nin
    )

    def wrapper(*args):
        args = [np.asarray(arg) for arg in args]

        # Fast path for when the arguments which are used in the cached
        # computation don't have batches.
        if all(args[i].ndim == core_ndims[i] for i in cache_arg_indices):
            return ufunc(*args)

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

        # Premake the output array.
        input_dtypes = tuple(arg.dtype for arg in args_t) + (None,)*ufunc.nout
        out_dtype = ufunc.resolve_dtypes(input_dtypes)[-1]
        out_final = np.empty(batch_shape, dtype=out_dtype, order='C')
        # a view of the output array with axes sorted as needed.
        out_t = np.transpose(out_final, axes=sorted_batch_axes)

        # Set out to the above view, but return the C contiguous output.
        # This avoids having non-contiguous output.
        ufunc(*args_t, out=out_t, order='C')
        return out_final
    return wrapper
