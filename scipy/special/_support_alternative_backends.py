import functools
import operator
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import ModuleType

import numpy as np
from scipy._lib._array_api import (
    array_namespace, scipy_namespace_for, is_numpy, is_dask, is_marray, is_jax_array,
    is_jax, xp_promote, xp_capabilities, get_native_namespace_name,
    is_array_api_obj, xp_result_device,
)
import scipy._external.array_api_extra as xpx
from . import _basic
from . import _ellip_harm
from . import _mathieu
from . import _spfun_stats
from . import _ufuncs
from ._ufunc_tools import _make_ufunc_wrapper


def _special_namespace_for(xp):
    spx = scipy_namespace_for(xp)
    return getattr(spx, "special", None)


def _ufunc_kwargs_extra_note(
        name=None, out_unsupported_backends=(), additional_note=None
):
    if (name is None) != (not out_unsupported_backends):
        raise ValueError(
            "`name` and `out_unsupported_backends` must either both be supplied "
            "or both be omitted."
        )

    out_note = ""
    additional_note = "" if additional_note is None else additional_note
    if name is not None:
        backend_names = {
            "cupy": "CuPy",
            "torch": "PyTorch",
        }
        backends = [
            backend_names[backend] for backend in out_unsupported_backends
        ]

        if len(backends) == 1:
            backend_text = f"the {backends[0]} backend"
        else:
            backend_text = f"the {' and '.join(backends)} backends"

        out_note = (
            f"``{name}`` does not currently support ``out`` for {backend_text}."
        )

    return (
        "For the NumPy backend, this function supports all\n"
        "    `NumPy ufunc keyword arguments "
        "<https://numpy.org/doc/stable/reference/ufuncs.html#optional-keyword-arguments>`_.\n"
        "    Other backends may support ``out``, but typically none of the other\n"
        "    ufunc kwargs. ``out`` is typically supported for CuPy and PyTorch, but\n"
        "    not currently in cases where SciPy relies on a generic Array API\n"
        "    implementation or, for PyTorch on CPU, falls back to the NumPy backend.\n"
        "    ``out`` is never supported for JAX because JAX arrays are immutable.\n"
        f"    {out_note}\n"
        f"    {additional_note}\n\n"
    )


@dataclass
class _FuncInfo:
    # NumPy-only function. IT MUST BE ELEMENTWISE.
    func: Callable
    # List of argument names to generate useful call signature for ufuncs.
    arg_names: list[str]
    # @xp_capabilities decorator, for the purpose of
    # documentation and unit testing. Omit to indicate
    # full support for all backends.
    xp_capabilities: Callable[[Callable], Callable] | None = None
    # Generic implementation to fall back on if there is no native dispatch
    # available. This is a function that accepts (main namespace, scipy namespace)
    # and returns the final callable, or None if not available.
    generic_impl: Callable[
        [ModuleType, ModuleType | None], Callable | None
    ] | None = None
    # Handle case where a backend uses an alternative name for a function.
    # Should map backend names to alternative function names.
    alt_names_map: dict[str, str] | None = None
    # Handle case where a backend has an implementation that takes its
    # arguments in a different order.
    alt_arg_orders_map: dict[str, tuple[int, ...]] | None = None
    # For testing purposes. Whether tests should only use integer arrays
    # for some arguments. This may be because noninteger args or not supported
    # or simply because the results for a particular backend do not match those
    # for NumPy. To make backend specific choices for positive_only,
    # pass in a dict mapping backend names to bool or tuple[bool, ...].
    int_only: tuple[bool, ...] | dict[str, tuple[bool, ...]] | None = None
    # For testing purposes, whether tests should only use positive values
    # for some arguments. If bool and equal to True, restrict to positive
    # values for all arguments. To restrict only some arguments to positive
    # values, pass a tuple of bool of the same length as the number of
    # arguments, the ith entry in the tuple controls positive_only for
    # the ith argument. To make backend specific choices for positive_only,
    # pass in a dict mapping backend names to bool or tuple[bool, ...].
    positive_only: (
        bool | tuple[bool, ...] | Mapping[str, bool | tuple[bool, ...]]
    ) = False
    # Some special functions are not ufuncs and ufunc-specific tests
    # should not be applied to these.
    is_ufunc: bool = True
    # Whether or not this is an elementwise function, or behaves like a
    # gufunc with nonzero core dimensionality in one or more inputs or outputs.
    is_elementwise: bool = True
    # For gufuncs, a function determining how output shapes depend on input shapes
    # to allow gufunc fallback to work with JAX JIT through lazy_apply.
    shape_mapper: Callable[..., ...] | None = None
    # Some non-ufunc special functions take only Python ints for some arguments.
    # If so, python_int_only should be a tuple of the same length as the number
    # of arguments,with value True if the corresponding argument needs to be a
    # Python int.
    # Can also take a dict mapping backends to such tuples if an argument being
    # Python int only is backend specific.
    python_int_only: dict[str, tuple[bool, ...]] | tuple[bool, ...] | None = None
    # Some functions which seem to be scalar also accept 0d arrays.
    scalar_or_0d_only: dict[str, tuple[bool, ...]] | tuple[bool, ...] | None = None
    # Some functions may not work well with very large integer valued arguments.
    test_large_ints: bool = True
    # Some non-ufunc special functions don't decay 0d arrays to scalar.
    produces_0d: bool = False
    # Whether or not uses native PyTorch or falls back to NumPy/SciPy. This
    # is needed because in PyTorch, the default dtype affects promotion
    # rules when mixing integer and floating dtypes, so relying on a
    # NumPy/SciPy fallback when the default dtype is other than float64 can lead
    # to float64 output when native PyTorch would have e.g. float32 output. This
    # must be accounted for in tests. Not putting this in xp_capabilities for now,
    # but in the future I think it's likely we may want to add a warning to
    # xp_capabilities when not using native PyTorch on CPU.
    torch_native: bool = True
    # Place a backend in this tuple if `func` is available as `xp.func` but not
    # available in the `scipy.special` namespace for this backend.
    # One example is `jax.numpy.sinc` being available but not `jax.scipy.special.sinc`.
    backends_with_func_in_xp: tuple[str, ...] = ()

    @property
    def name(self):
        return self.func.__name__

    # These are needed by @lru_cache below
    def __hash__(self):
        return hash(self.func)

    def __eq__(self, other):
        return isinstance(other, _FuncInfo) and self.func == other.func

    @property
    def n_args(self):
        return len(self.arg_names)

    @property
    def wrapper(self):
        if self.name in globals():
            # Already initialised. We are likely in a unit test.
            # Return function potentially overridden by xpx.testing.lazy_xp_function.
            import scipy.special
            return getattr(scipy.special, self.name)

        if self.is_ufunc:
            def wrapped(*args, **kwargs):
                xp = array_namespace(*args)
                return self._wrapper_for(xp)(*args, **kwargs)

            func = _make_ufunc_wrapper(
                wrapped,
                self.func,
                self.name,
                self.arg_names,
                self.func.__doc__,
            )
        else:
            @functools.wraps(self.func)
            def func(*args, **kwargs):
                xp = array_namespace(*args, *kwargs.values())
                return self._wrapper_for(xp)(*args, **kwargs)

            # needed to allow pickling
            func.__module__ = "scipy.special"
            func.__qualname__ = self.name

        if self.name in {"mathieu_cem", "mathieu_sem"}:
            # Warn outside backend dispatch and register the public wrapper.
            func = _mathieu._deprecated_mathieu(func, self.name[:-1])
        capabilities = self.xp_capabilities or xp_capabilities()
        cap_func = capabilities(func)
        assert cap_func is func
        return func

    @functools.lru_cache(1000)
    def _wrapper_for(self, xp):
        if is_numpy(xp):
            return self.func

        # If a native implementation is available, use that
        in_xp = get_native_namespace_name(xp) in self.backends_with_func_in_xp
        namespace = xp if in_xp else _special_namespace_for(xp)
        f = _get_native_func(
            xp, namespace, self.name, alt_names_map=self.alt_names_map
        )
        if f is not None:
            order = (self.alt_arg_orders_map or {}).get(get_native_namespace_name(xp))
            if order is not None:
                native_f = f

                def f(*args, **kwargs):
                    return native_f(*(args[i] for i in order), **kwargs)

            return f

        if in_xp:
            # when namespace is passed to self.generic_impl below, we want to
            # make sure that it is the special namespace for xp and not xp
            # itself, so raise if xp was incorrectly placed in
            # `backends_with_func_in_xp`.
            raise RuntimeError(
                f"func {self.func} is not available as {xp.__name__}.{self.func}"
                f" but {xp.__name__} was passed in ``in_xp``."
            )

        # If generic Array API implementation is available, use that
        if self.generic_impl is not None:
            f_generic = self.generic_impl(xp, namespace)
            if f_generic is not None:
                return f_generic

        if is_marray(xp):
            if not self.is_elementwise:
                raise NotImplementedError(
                    f"{self.name} is not supported with marray because it is "
                    "not elementwise."
                )
            # Unwrap the array, apply the function on the wrapped namespace,
            # and then re-wrap it.
            # IMPORTANT: this only works because we exclude functions that are
            # not elementwise. Otherwise, we would not be able to define a
            # general rule for mask propagation.

            _f = globals()[self.name]  # Allow nested wrapping
            def f(*args, out=None, _f=_f, xp=xp, **kwargs):
                if out is not None:
                    raise NotImplementedError(
                        "`out` is not supported with marray."
                    )
                data_args = [getattr(arg, 'data', arg) for arg in args]
                res = _f(*data_args, **kwargs)
                mask = functools.reduce(operator.or_,
                                        (getattr(arg, 'mask', False) for arg in args))
                if isinstance(res, tuple):
                    return tuple(xp.asarray(x, mask=mask) for x in res)
                return xp.asarray(res, mask=mask)

            return f

        if is_dask(xp):
            if not self.is_elementwise:
                raise NotImplementedError(
                    f"{self.name} is not supported with dask.array because it is "
                    "not elementwise."
                )
            # Apply the function to each block of the Dask array.
            # IMPORTANT: map_blocks works only because we exclude functions which
            # are not elementwise. It would be a grave mistake to apply this to gufuncs
            # or any other function with reductions, as they would change their
            # output depending on chunking!

            _f = globals()[self.name]  # Allow nested wrapping
            def f(*args, out=None, _f=_f, xp=xp, **kwargs):
                if out is not None:
                    raise NotImplementedError(
                        "The `out` parameter is not supported for dask.array."
                    )
                # Hide dtype kwarg from map_blocks
                return xp.map_blocks(functools.partial(_f, **kwargs), *args)

            return f

        # As a final resort, use the NumPy/SciPy implementation
        _f = self.func

        if is_jax(xp) and self.is_ufunc:
            # if this is a ufunc, we can use the resolve_dtypes method to figure
            # out what the output dtype should be and use lazy_apply to make this
            # work with the JAX JIT. Since pure_callback used in lazy_apply will
            # convert all of the inputs to eager JAX arrays, we will also need to
            # get the input dtypes inferred from resolve_dtypes so that Python
            # scalar inputs can be cast to the correct dtype under NEP50 weak
            # promotion rules rather than getting promoted to the default currently
            # set in JAX. One cannot just use xp_promote for the input dtypes because
            # some ufuncs have integer only args.
            nin, nout = self.func.nin, self.func.nout
            if self.is_elementwise:
                def f(*args, _f=_f, xp=xp, **kwargs):
                    if kwargs:
                        raise NotImplementedError(
                            "ufunc keyword arguments are not supported "
                            f"for {self.name} with backend {xp.__name__}."
                        )
                    dtypes = (
                        arg.dtype if is_jax_array(arg) else type(arg) for arg in args
                    )
                    # result_dtypes needs an arg for the dtype of the optional out
                    # params. Uses None to request output dtype inference.
                    dtypes = (*dtypes, *(None,) * nout)
                    # JAX uses NumPy dtypes so we can just pass these directly to
                    # resolve_dtypes. TODO: generalize to other lazy backends.
                    dtypes = self.func.resolve_dtypes(dtypes)
                    out_dtypes = dtypes[nin:]
                    out_dtype = tuple(out_dtypes) if nout > 1 else out_dtypes[0]
                    args = [
                        xp.asarray(arg, dtype=dtype)
                        for arg, dtype in zip(args, dtypes[:nin])
                    ]

                    if nout > 1:
                        bshape = xp.broadcast_shapes(*(arg.shape for arg in args))
                        shape = (bshape,) * nout
                    else:
                        shape = None
                    return xpx.lazy_apply(
                        _f, *args, shape=shape, xp=xp, as_numpy=True, dtype=out_dtype
                    )
            else:
                # gufunc path.
                if self.shape_mapper is None:
                    raise RuntimeError(
                        f"No shape mapper is defined for non-elementwise ufunc "
                        f"{self.name!r}; set `_FuncInfo.shape_mapper`."
                    )
                shape_mapper = self.shape_mapper
                # The below function is the same as above but makes use of shape_mapper
                # to get the output shape which can no longer be inferred directly by
                # broadcasting the input shapes.
                def f(*args, _f=_f, xp=xp, shape_mapper=shape_mapper, **kwargs):
                    dtypes = (
                        arg.dtype if is_jax_array(arg) else type(arg)
                        for arg in args
                    )
                    dtypes = (*dtypes, *(None,) * nout)
                    dtypes = self.func.resolve_dtypes(dtypes)

                    out_dtypes = dtypes[nin:]
                    out_dtype = tuple(out_dtypes) if nout > 1 else out_dtypes[0]

                    args = [
                        xp.asarray(arg, dtype=dtype)
                        for arg, dtype in zip(args, dtypes[:nin])
                    ]

                    # The shape mapper will decide if the kwargs are actually supported
                    # for this backend with JAX and raise if not.
                    shape = shape_mapper(
                        *(arg.shape for arg in args), **kwargs
                    )

                    return xpx.lazy_apply(
                        _f, *args, shape=shape, xp=xp,
                        as_numpy=True, dtype=out_dtype, **kwargs
                    )
        else:
            def f(*args, out=None, _f=_f, xp=xp, **kwargs):
                if out is not None:
                    raise NotImplementedError(
                        f"`out` parameter is not supported for {self.name}"
                        f" with backend {xp.__name__}."
                    )
                if self.is_ufunc and not self.is_elementwise:
                    gufunc_kwargs = {}
                    for kwarg in ("axis", "axes", "keepdims"):
                        if kwarg in kwargs:
                            gufunc_kwargs[kwarg] = kwargs.pop(kwarg)
                    if kwargs:
                        kwarg_names = ", ".join(f"`{name}`" for name in kwargs)
                        raise NotImplementedError(
                            f"received kwargs: {kwarg_names} that are not supported"
                            f" for {self.name} with backend {xp.__name__}."
                        )
                    kwargs = gufunc_kwargs
                elif self.is_ufunc and kwargs:
                    raise NotImplementedError(
                        "ufunc keyword arguments are not supported "
                        f"for {self.name} with backend {xp.__name__}."
                    )
                # The NumPy round-trip must return results on the device of the
                # input arrays, not on the backend's default device (see gh-22680)
                device = xp_result_device(*args)

                # Check with `is_array_api_obj` to keep Python scalars untouched so that
                # NEP50 can be followed.
                args = [
                    np.asarray(arg) if is_array_api_obj(arg) else arg for arg in args
                ]
                out = _f(*args, **kwargs)

                if isinstance(out, tuple):
                    return tuple(xp.asarray(out_i, device=device) for out_i in out)
                return xp.asarray(out, device=device)

        return f


def _get_native_func(xp, namespace, f_name, *, alt_names_map=None):
    if alt_names_map is None:
        alt_names_map = {}
    f_name = alt_names_map.get(get_native_namespace_name(xp), f_name)
    f = getattr(namespace, f_name, None)
    if f is None and hasattr(xp, 'special'):
        # Currently dead branch, in anticipation of 'special' Array API extension
        # https://github.com/data-apis/array-api/issues/725
        f = getattr(xp.special, f_name, None)
    return f


def _rel_entr(xp, spsx):
    def __rel_entr(x, y, *, xp=xp):
        # https://github.com/data-apis/array-api-extra/issues/160
        mxp = array_namespace(x._meta, y._meta) if is_dask(xp) else xp
        x, y = xp_promote(x, y, broadcast=True, force_floating=True, xp=xp)
        xy_pos = (x > 0) & (y > 0)
        xy_inf = xp.isinf(x) & xp.isinf(y)
        res = xpx.apply_where(
            xy_pos & ~xy_inf,
            (x, y),
            # Note: for very large x, this can overflow.
            lambda x, y: x * (mxp.log(x) - mxp.log(y)),
            fill_value=xp.inf
        )
        res = xpx.at(res)[(x == 0) & (y >= 0)].set(0)
        res = xpx.at(res)[xp.isnan(x) | xp.isnan(y) | (xy_pos & xy_inf)].set(xp.nan)
        return res

    return __rel_entr


def _xlogy(xp, spsx):
    def __xlogy(x, y, *, xp=xp):
        x, y = xp_promote(x, y, force_floating=True, xp=xp)
        with np.errstate(divide='ignore', invalid='ignore'):
            temp = x * xp.log(y)
        return xp.where(x == 0., 0., temp)
    return __xlogy



def _chdtr(xp, spsx):
    # The difference between this and just using `gammainc`
    # defined by `get_array_special_func` is that if `gammainc`
    # isn't found, we don't want to use the SciPy version; we'll
    # return None here and use the SciPy version of `chdtr`.
    gammainc = _get_native_func(xp, spsx, 'gammainc')
    if gammainc is None:
        return None

    def __chdtr(v, x):
        res = gammainc(v / 2, x / 2)  # this is almost all we need
        # The rest can be removed when google/jax#20507 is resolved
        mask = (v == 0) & (x > 0)  # JAX returns NaN
        res = xp.where(mask, 1., res)
        mask = xp.isinf(v) & xp.isinf(x)  # JAX returns 1.0
        return xp.where(mask, xp.nan, res)
    return __chdtr


def _chdtrc(xp, spsx):
    # The difference between this and just using `gammaincc`
    # defined by `get_array_special_func` is that if `gammaincc`
    # isn't found, we don't want to use the SciPy version; we'll
    # return None here and use the SciPy version of `chdtrc`.
    gammaincc = _get_native_func(xp, spsx, 'gammaincc')
    if gammaincc is None:
        return None

    def __chdtrc(v, x):
        res = xp.where(x >= 0, gammaincc(v/2, x/2), 1)
        i_nan = ((x == 0) & (v == 0)) | xp.isnan(x) | xp.isnan(v) | (v < 0)
        res = xp.where(i_nan, xp.nan, res)
        return res
    return __chdtrc


def _betaincc(xp, spsx):
    betainc = _get_native_func(xp, spsx, 'betainc')
    if betainc is None:
        return None

    def __betaincc(a, b, x):
        # not perfect; might want to just rely on SciPy
        return betainc(b, a, 1-x)
    return __betaincc


def _stdtr(xp, spsx):
    betainc = _get_native_func(xp, spsx, 'betainc')
    if betainc is None:
        return None

    def __stdtr(df, t):
        x = df / (t ** 2 + df)
        tail = betainc(df / 2, 0.5, x) / 2
        return xp.where(t < 0, tail, 1 - tail)

    return __stdtr


def _stdtrit(xp, spsx):
    # Need either native stdtr or native betainc
    stdtr = _get_native_func(xp, spsx, 'stdtr') or _stdtr(xp, spsx)
    # If betainc is not defined, the root-finding would be done with `xp`
    # despite `stdtr` being evaluated with SciPy/NumPy `stdtr`. Save the
    # conversions: in this case, just evaluate `stdtrit` with SciPy/NumPy.
    if stdtr is None:
        return None

    from scipy.optimize.elementwise import bracket_root, find_root

    def __stdtrit(df, p):
        def fun(t, df, p):  return stdtr(df, t) - p
        res_bracket = bracket_root(fun, xp.zeros_like(p), args=(df, p))
        res_root = find_root(fun, res_bracket.bracket, args=(df, p))
        return res_root.x

    return __stdtrit


def _ndtri_exp(xp, spsx):
    # The difference between this and just using `ndtri`
    # defined by `get_array_special_func` is that if `ndtri`
    # isn't found, we don't want to use the SciPy version; we'll
    # return None here and use the SciPy version of `ndtri_exp`.
    ndtri = _get_native_func(xp, spsx, 'ndtri')
    if ndtri is None:
        return None

    def __ndtri_exp(x):
        return ndtri(xp.exp(x))  # this is almost all we need

    return __ndtri_exp


def _poisson_binom_cdf_shape_mapper(k_shape, p_shape, axis=-1, **kwargs):
    # Used to infer out shape for lazy_apply under JAX JIT.
    if kwargs:
        raise NotImplementedError(
            "`poisson_binom_cdf` does not support gufunc kwargs other than "
            "`axis` with the JAX backend."
        )
    axis = operator.index(axis)
    if axis < 0:
        axis += len(p_shape)
    if not 0 <= axis < len(p_shape):
        raise np.exceptions.AxisError(axis, ndim=len(p_shape))
    p_batch_shape = p_shape[:axis] + p_shape[axis + 1:]
    return np.broadcast_shapes(k_shape, p_batch_shape)


# Inventory of automatically dispatched functions

# PyTorch doesn't implement `betainc`.
# On torch CPU we can fall back to NumPy, but on GPU it won't work.
def _needs_betainc(name):
    return xp_capabilities(
        cpu_only=True, exceptions=["jax.numpy", "cupy"],
        extra_note=_ufunc_kwargs_extra_note(name, ["torch"]),
    )


orthogonal_poly_eval_additional_note = (
    "With PyTorch, noninteger values of the degree `n` are cast to"
    " integers before evaluation."
)

_special_funcs = (
    # Airy functions
    _FuncInfo(
        _ufuncs.airy, ["z"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("airy", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.airye, ["z"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("airye", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.itairy, ["x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("itairy", ["torch"]),
        ),
        torch_native=False,
    ),
    # Elliptic functions and integrals
    _FuncInfo(
        _ufuncs.ellipj, ["u", "m"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("ellipj", ["torch"]),
        ),
        positive_only={"cupy": (False, True)},
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ellipk, ["m"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("ellipk", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ellipkm1, ["p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("ellipkm1", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ellipkinc, ["phi", "m"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("ellipkinc", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ellipe, ["m"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("ellipe", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ellipeinc, ["phi", "m"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("ellipeinc", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.elliprc, ["x", "y"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("elliprc", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.elliprd, ["x", "y", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("elliprd", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.elliprf, ["x", "y", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("elliprf", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.elliprg, ["x", "y", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("elliprg", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.elliprj, ["x", "y", "z", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("elliprj", ["torch"]),
        ),
        torch_native=False,
    ),
    # Bessel functions
    _FuncInfo(
        _ufuncs.jv, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("jv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.jve, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("jve", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.yn, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("yn", ["torch"]),
        ),
        int_only=(True, False),
        positive_only={"cupy": (True, False)},
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.yv, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("yv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.yve, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("yve", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.iv, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("iv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ive, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("ive", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.kn, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("kn", ["torch"]),
        ),
        int_only=(True, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.kv, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("kv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.kve, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("kve", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.hankel1, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("hankel1", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.hankel1e, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("hankel1e", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.hankel2, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("hankel2", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.hankel2e, ["v", "z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("hankel2e", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.wright_bessel, ["a", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("wright_bessel", ["torch"]),
        ),
        positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.log_wright_bessel, ["a", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("log_wright_bessel", ["torch"]),
        ),
        positive_only=True,
        torch_native=False,
    ),
    # Faster versions of common Bessel functions
    _FuncInfo(
        _ufuncs.j0, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "bessel_j0"},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.j1, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "bessel_j1"},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.y0, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "bessel_y0"},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.y1, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "bessel_y1"},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.i0, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.i0e, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.i1, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.i1e, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.k0, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "modified_bessel_k0"},
    ),
    _FuncInfo(
        _ufuncs.k0e, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "scaled_modified_bessel_k0"},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.k1, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "modified_bessel_k1"},
    ),
    _FuncInfo(
        _ufuncs.k1e, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note(),
        ),
        alt_names_map={"torch": "scaled_modified_bessel_k1"},
        test_large_ints=False,
    ),
    # Integrals of Bessel functions
    _FuncInfo(
        _ufuncs.itj0y0, ["x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("itj0y0", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.it2j0y0, ["x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("it2j0y0", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.iti0k0, ["x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("iti0k0", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.it2i0k0, ["x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("it2i0k0", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.besselpoly, ["a", "lmb", "nu"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("besselpoly", ["torch"]),
        ),
        torch_native=False,
    ),
    # Spherical Bessel functions
    # --
    # Struve functions
    _FuncInfo(
        _ufuncs.struve, ["v", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("struve", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.modstruve, ["v", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("modstruve", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.itstruve0, ["x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("itstruve0", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.it2struve0, ["x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("it2struve0", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.itmodstruve0, ["x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("itmodstruve0", ["torch"]),
        ),
        torch_native=False,
    ),
    # Raw Statistical Functions
    ## Binomial Distribution
    _FuncInfo(
        _ufuncs.bdtr, ["k", "n", "p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("bdtr", ["torch"]),
        ),
        int_only=(False, True, False),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.bdtrc, ["k", "n", "p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("bdtrc", ["torch"]),
        ),
        int_only=(False, True, False),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.bdtri, ["k", "n", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("bdtri", ["torch"]),
        ),
        int_only=(False, True, False),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.bdtrik, ["y", "n", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("bdtrik", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.bdtrin, ["k", "y", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("bdtrin", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Beta Distribution
    _FuncInfo(
        _ufuncs.btdtria, ["p", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("btdtria", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.btdtrib, ["a", "p", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("btdtrib", ["torch"]),
        ),
        torch_native=False,
    ),
    ## F Distribution
    _FuncInfo(
        _ufuncs.fdtr, ["dfn", "dfd", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("fdtr", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.fdtrc, ["dfn", "dfd", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("fdtrc", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.fdtri, ["dfn", "dfd", "p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("fdtri", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.fdtridfd, ["dfn", "p", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("fdtridfd", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Gamma distribution
    _FuncInfo(
        _ufuncs.gdtr, ["a", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("gdtr", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gdtrc, ["a", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("gdtrc", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gdtria, ["p", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("gdtria", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gdtrib, ["a", "p", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("gdtrib", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gdtrix, ["a", "b", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("gdtrix", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Negative Binomial Distribution
    _FuncInfo(
        _ufuncs.nbdtr, ["k", "n", "p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("nbdtr", ["torch"]),
        ),
        int_only=(True, True, False),
        positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nbdtrc, ["k", "n", "p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("nbdtrc", ["torch"]),
        ),
        int_only=(True, True, False),
        positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nbdtri, ["k", "n", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("nbdtri", ["torch"]),
        ),
        int_only=(True, True, False),
        positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nbdtrik, ["y", "n", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nbdtrik", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nbdtrin, ["k", "y", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nbdtrin", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Noncentral F distribution
    _FuncInfo(
        _ufuncs.ncfdtr, ["dfn", "dfd", "nc", "f"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("ncfdtr", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ncfdtridfd, ["dfn", "p", "nc", "f"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("ncfdtridfd", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ncfdtridfn, ["p", "dfd", "nc", "f"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("ncfdtridfn", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ncfdtri, ["dfn", "dfd", "nc", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("ncfdtri", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.ncfdtrinc, ["dfn", "dfd", "p", "f"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("ncfdtrinc", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Noncentral t distribution
    _FuncInfo(
        _ufuncs.nctdtr, ["df", "nc", "t"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nctdtr", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nctdtridf, ["p", "nc", "t"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nctdtridf", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nctdtrit, ["df", "nc", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nctdtrit", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nctdtrinc, ["df", "p", "t"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nctdtrinc", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Normal distribution
    _FuncInfo(
        _ufuncs.nrdtrimn, ["p", "std", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nrdtrimn", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.nrdtrisd, ["mn", "p", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("nrdtrisd", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.log_ndtr, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.ndtr, ["x"],
        xp_capabilities(
            mparray=True,
            extra_note=_ufunc_kwargs_extra_note(),
        ),
    ),
    _FuncInfo(
        _ufuncs.ndtri, ["p"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.ndtri_exp, ["y"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note("ndtri_exp", ["torch"])),
        generic_impl=_ndtri_exp,
        test_large_ints=False,
    ),
    ## Poisson distribution
    _FuncInfo(
        _ufuncs.pdtr, ["k", "m"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("pdtr", ["torch"]),
        ),
        positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pdtrc, ["k", "m"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("pdtrc", ["torch"]),
        ),
        positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pdtri, ["k", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("pdtri", ["torch"]),
        ),
        int_only=(True, False),
        positive_only=True,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pdtrik, ["p", "m"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("pdtrik", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Student t distribution
    _FuncInfo(
        _ufuncs.stdtr, ["df", "t"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["jax.numpy", "cupy"],  # needs betainc
            extra_note=_ufunc_kwargs_extra_note("stdtr", ["cupy", "torch"]),
        ),
        generic_impl=_stdtr,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.stdtridf, ["p", "t"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("stdtridf", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.stdtrit, ["df", "p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],  # needs betainc
            skip_backends=[("jax.numpy", "no scipy.optimize support")],
            extra_note=_ufunc_kwargs_extra_note("stdtrit", ["cupy", "torch"]),
        ),
        generic_impl=_stdtrit,
        torch_native=False,
    ),
    # Chi square distribution
    _FuncInfo(
        _ufuncs.chdtr, ["v", "x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note("chdtr", ["torch"])),
        generic_impl=_chdtr,
    ),
    _FuncInfo(
        _ufuncs.chdtrc, ["v", "x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note("chdtrc", ["torch"])),
        generic_impl=_chdtrc,
        # scipy/scipy#20972
        positive_only={"cupy": True, "jax.numpy": True, "torch": True},
    ),
    _FuncInfo(
        _ufuncs.chdtri, ["v", "p"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("chdtri", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.chdtriv, ["p", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("chdtriv", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Non-central chi square distribution
    _FuncInfo(
        _ufuncs.chndtr, ["x", "df", "nc"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("chndtr", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.chndtridf, ["x", "p", "nc"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("chndtridf", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.chndtrinc, ["x", "df", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("chndtrinc", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.chndtrix, ["p", "df", "nc"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("chndtrix", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Kolmogorov distribution
    _FuncInfo(
        _ufuncs.smirnov, ["n", "d"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("smirnov", ["torch"]),
        ),
        int_only=(True, False),
        positive_only=(True, False),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.smirnovi, ["n", "p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("smirnovi", ["torch"]),
        ),
        int_only=(True, False),
        positive_only=(True, False),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.kolmogorov, ["x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("kolmogorov", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.kolmogi, ["p"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("kolmogi", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Box-Cox transformation
    _FuncInfo(
        _ufuncs.boxcox, ["x", "lmbda"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("boxcox", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.boxcox1p, ["x", "lmbda"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("boxcox1p", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.inv_boxcox, ["y", "lmbda"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("inv_boxcox", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.inv_boxcox1p, ["y", "lmbda"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("inv_boxcox1p", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Sigmoidal functions
    _FuncInfo(
        _ufuncs.logit, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.expit, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.log_expit, ["x"],
        xp_capabilities(
            cpu_only=True, exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("log_expit", ["torch"]),
        ),
        torch_native=False,
    ),
    ## Miscellaneous
    _FuncInfo(
        _ufuncs.tklmbda, ["x", "lmbda"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("tklmbda", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.owens_t, ["h", "a"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("owens_t", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _spfun_stats.poisson_binom_cdf, ["k", "p"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "not elementwise")],
            extra_note=_ufunc_kwargs_extra_note(
                "poisson_binom_cdf",
                ["torch"],
                "It supports the ``axis`` kwarg for all supported backends.",
            ),
        ),
        is_elementwise=False,
        shape_mapper=_poisson_binom_cdf_shape_mapper,
        int_only=(True, False),
        test_large_ints=False,
        torch_native=False,
    ),
    # Information theory functions
    _FuncInfo(
        _ufuncs.entr, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.rel_entr, ["x", "y"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note("rel_entr", ["torch"])),
        generic_impl=_rel_entr,
    ),
    _FuncInfo(
        _ufuncs.kl_div, ["x", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("kl_div", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.huber, ["delta", "r"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("huber", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pseudo_huber, ["delta", "r"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("pseudo_huber", ["torch"]),
        ),
        torch_native=False,
    ),
    # Gamma and related functions
    _FuncInfo(
        _ufuncs.gamma, ["z"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("gamma", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gammaln, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.loggamma, ["z"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("loggamma", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gammasgn, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("gammasgn", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gammainc, ["a", "x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.gammaincinv, ["a", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("gammaincinv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.gammaincc, ["a", "x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
        # google/jax#20699
        positive_only={"jax.numpy": True},
    ),
    _FuncInfo(
        _ufuncs.gammainccinv, ["a", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("gammainccinv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.log_gammainc, ["a", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("log_gammainc", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.log_gammaincc, ["a", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("log_gammaincc", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.beta, ["a", "b"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("beta", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.betaln, ["a", "b"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("betaln", ["torch"]),
        ),
        # For betaln, nan mismatches can occur at negative integer a or b of
        # sufficiently large magnitude.
        positive_only={"jax.numpy": True},
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.betainc, ["a", "b", "x"],
        _needs_betainc("betainc"),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.betaincc, ["a", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["jax.numpy", "cupy"],  # needs betainc
            extra_note=_ufunc_kwargs_extra_note("betaincc", ["cupy", "torch"]),
        ),
        generic_impl=_betaincc,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.betaincinv, ["a", "b", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("betaincinv", ["torch"]),
        ),
        positive_only=True,
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.betainccinv, ["a", "b", "y"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("betainccinv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.psi, ["z"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
        alt_names_map={"jax.numpy": "digamma"},
    ),
    _FuncInfo(
        _ufuncs.digammainv, ["y"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("digammainv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.rgamma, ["z"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("rgamma", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _basic.polygamma, ["n", "x"],
        is_ufunc=False,
        scalar_or_0d_only={"torch": (True, False)},
        produces_0d=True,
        int_only=(True, False),
        positive_only={"torch": (True, False), "jax.numpy": True},
        test_large_ints=False,
    ),
    _FuncInfo(
        _spfun_stats.multigammaln, ["a", "d"],
        is_ufunc=False,
        python_int_only={
            "cupy": (False, True),
            "jax.numpy": (False, True),
            "torch": (False, True),
        },
        scalar_or_0d_only={
            "array_api_strict": (False, True),
            "numpy": (False, True),
            "dask.array": (False, True),
            "marray": (False, True),
        },
        int_only=(False, True),
        positive_only=True,
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.poch, ["z", "m"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("poch", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    # Error function and Fresnel integrals
    _FuncInfo(
        _ufuncs.erf, ["z"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.erfc, ["x"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.erfcx, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
            extra_note=_ufunc_kwargs_extra_note("erfcx", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.erfi, ["z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("erfi", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.erfinv, ["y"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.erfcinv, ["y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("erfcinv", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.wofz, ["z"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("wofz", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.dawsn, ["x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("dawsn", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.fresnel, ["z"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("fresnel", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.modfresnelp, ["z"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("modfresnelp", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.modfresnelm, ["z"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("modfresnelm", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.voigt_profile, ["x", "sigma", "gamma"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("voigt_profile", ["torch"]),
        ),
        torch_native=False,
    ),
    # Legendre functions
    _FuncInfo(
        _ufuncs.lpmv, ["m", "v", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("lpmv", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    # Ellipsoidal harmonics
    # --
    _FuncInfo(
        _ellip_harm.ellip_harm, ["h2", "k2", "n", "p", "s", "signm=1", "signn=1"],
        xp_capabilities(cpu_only=True, jax_jit=False),
        is_ufunc=False,
        int_only=(False, False, True, True, False, True, True),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ellip_harm.ellip_harm_2, ["h2", "k2", "n", "p", "s"],
        xp_capabilities(cpu_only=True, jax_jit=False),
        is_ufunc=False,
        int_only=(False, False, True, True, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ellip_harm.ellip_normal, ["h2", "k2", "n", "p"],
        xp_capabilities(cpu_only=True, jax_jit=False),
        is_ufunc=False,
        int_only=(False, False, True, True),
        test_large_ints=False,
        torch_native=False,
    ),
    # Orthogonal polynomials
    _FuncInfo(
        _ufuncs.eval_legendre, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "legendre_polynomial_p"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only={"torch": (True, False)},
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_chebyt, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "chebyshev_polynomial_t"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only={"torch": (True, False)},
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_chebyu, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "chebyshev_polynomial_u"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only={"torch": (True, False)},
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_chebyc, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("eval_chebyc", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.eval_chebys, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("eval_chebys", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.eval_jacobi, ["n", "alpha", "beta", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("eval_jacobi", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.eval_laguerre, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "laguerre_polynomial_l"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only={"torch": (True, False)},
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_genlaguerre, ["n", "alpha", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("eval_genlaguerre", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.eval_hermite, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "hermite_polynomial_h"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only=(True, False),
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_hermitenorm, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "hermite_polynomial_he"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only=(True, False),
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_gegenbauer, ["n", "alpha", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("eval_gegenbauer", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.eval_sh_legendre, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("eval_sh_legendre", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.eval_sh_chebyt, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "shifted_chebyshev_polynomial_t"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only={"torch": (True, False)},
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_sh_chebyu, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["torch"],
            extra_note=_ufunc_kwargs_extra_note(
                additional_note=orthogonal_poly_eval_additional_note,
            ),
        ),
        alt_names_map={"torch": "shifted_chebyshev_polynomial_u"},
        alt_arg_orders_map={"torch": (1, 0)},
        int_only={"torch": (True, False)},
        positive_only={"torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.eval_sh_jacobi, ["n", "p", "q", "x"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("eval_sh_jacobi", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    # Hypergeometric functions
    _FuncInfo(
        _ufuncs.hyp2f1, ["a", "b", "c", "z"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("hyp2f1", ["torch"]),
        ),
        positive_only={"jax.numpy": True},
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.hyp1f1, ["a", "b", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("hyp1f1", ["torch"]),
        ),
        positive_only={"jax.numpy": True},
        test_large_ints=False,
        torch_native=False,
    ),
    # Parabolic cylinder functions
    # --
    _FuncInfo(
        _ufuncs.pbdv, ["v", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pbdv", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pbvv, ["v", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pbvv", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pbwa, ["a", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pbwa", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    # Mathieu and related functions
    _FuncInfo(
        _ufuncs.mathieu_a, ["m", "q"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("mathieu_a", ["torch"]),
        ),
        int_only=(True, False),
        positive_only=(True, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.mathieu_b, ["m", "q"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("mathieu_b", ["torch"]),
        ),
        int_only=(True, False),
        positive_only=(True, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _mathieu.mathieu_sem, ["m", "q", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("mathieu_sem", ["torch"]),
        ),
        int_only=(True, False, False),
        positive_only=(True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _mathieu.mathieu_cem, ["m", "q", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("mathieu_cem", ["torch"]),
        ),
        int_only=(True, False, False),
        positive_only=(True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.mathieu_modcem1, ["m", "q", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("mathieu_modcem1", ["torch"]),
        ),
        int_only=(True, False, False),
        positive_only=(True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.mathieu_modcem2, ["m", "q", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("mathieu_modcem2", ["torch"]),
        ),
        int_only=(True, False, False),
        positive_only=(True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.mathieu_modsem1, ["m", "q", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("mathieu_modsem1", ["torch"]),
        ),
        int_only=(True, False, False),
        positive_only=(True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.mathieu_modsem2, ["m", "q", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("mathieu_modsem2", ["torch"]),
        ),
        int_only=(True, False, False),
        positive_only=(True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    # Spheroidal wave functions
    # --
    _FuncInfo(
        _ufuncs.pro_ang1, ["m", "n", "c", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pro_ang1", ["torch"]),
        ),
        int_only=(True, True, False, False),
        positive_only=(True, True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pro_rad1, ["m", "n", "c", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pro_rad1", ["torch"]),
        ),
        int_only=(True, True, False, False),
        positive_only=(True, True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pro_rad2, ["m", "n", "c", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pro_rad2", ["torch"]),
        ),
        int_only=(True, True, False, False),
        positive_only=(True, True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.obl_ang1, ["m", "n", "c", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("obl_ang1", ["torch"]),
        ),
        int_only=(True, True, False, False),
        positive_only=(True, True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.obl_rad1, ["m", "n", "c", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("obl_rad1", ["torch"]),
        ),
        int_only=(True, True, False, False),
        positive_only=(True, True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.obl_rad2, ["m", "n", "c", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("obl_rad2", ["torch"]),
        ),
        int_only=(True, True, False, False),
        positive_only=(True, True, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pro_cv, ["m", "n", "c"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("pro_cv", ["torch"]),
        ),
        int_only=(True, True, False),
        positive_only=(True, True, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.obl_cv, ["m", "n", "c"],
        xp_capabilities(
            cpu_only=True,
            extra_note=_ufunc_kwargs_extra_note("obl_cv", ["torch"]),
        ),
        int_only=(True, True, False),
        positive_only=(True, True, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pro_ang1_cv, ["m", "n", "c", "cv", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pro_ang1_cv", ["torch"]),
        ),
        int_only=(True, True, False, False, False),
        positive_only=(True, True, False, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pro_rad1_cv, ["m", "n", "c", "cv", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pro_rad1_cv", ["torch"]),
        ),
        int_only=(True, True, False, False, False),
        positive_only=(True, True, False, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.pro_rad2_cv, ["m", "n", "c", "cv", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("pro_rad2_cv", ["torch"]),
        ),
        int_only=(True, True, False, False, False),
        positive_only=(True, True, False, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.obl_ang1_cv, ["m", "n", "c", "cv", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("obl_ang1_cv", ["torch"]),
        ),
        int_only=(True, True, False, False, False),
        positive_only=(True, True, False, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.obl_rad1_cv, ["m", "n", "c", "cv", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("obl_rad1_cv", ["torch"]),
        ),
        int_only=(True, True, False, False, False),
        positive_only=(True, True, False, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.obl_rad2_cv, ["m", "n", "c", "cv", "x"],
        xp_capabilities(
            cpu_only=True,
            skip_backends=[("dask.array", "multiple outputs")],
            extra_note=_ufunc_kwargs_extra_note("obl_rad2_cv", ["torch"]),
        ),
        int_only=(True, True, False, False, False),
        positive_only=(True, True, False, False, False),
        test_large_ints=False,
        torch_native=False,
    ),
    # Kelvin functions
    # --
    # Lambert w and related functions
    # --
    # Other special functions
    _FuncInfo(
        _ufuncs.binom, ["x", "y"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("binom", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.expn, ["n", "x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("expn", ["torch"]),
        ),
        int_only=(True, False),
        # Inconsistent behavior for negative n. expn is not defined here without
        # taking analytic continuation.
        positive_only=True,
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.exp1, ["z"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("exp1", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.expi, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("expi", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.spence, ["z"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["jax.numpy"],
            extra_note=_ufunc_kwargs_extra_note("spence", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _basic.zeta, ["x", "q=None"],
        is_ufunc=False,
        positive_only={"jax.numpy": True, "torch": (True, False)},
        test_large_ints=False,
    ),
    _FuncInfo(
        _ufuncs.zetac, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("zetac", ["torch"]),
        ),
        torch_native=False,
    ),
    # Convenience functions
    _FuncInfo(
        _ufuncs.cbrt, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("cbrt", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.exp10, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("exp10", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.exp2, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("exp2", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.radian, ["d", "m", "s"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("radian", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.cosdg, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("cosdg", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.sindg, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("sindg", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.tandg, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("tandg", ["torch"]),
        ),
        test_large_ints=False,
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.cotdg, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("cotdg", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.cosm1, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("cosm1", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _ufuncs.xlogy, ["x", "y"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note("xlogy", ["torch"])),
        generic_impl=_xlogy,
    ),
    _FuncInfo(
        _ufuncs.xlog1py, ["x", "y"],
        xp_capabilities(extra_note=_ufunc_kwargs_extra_note()),
    ),
    _FuncInfo(
        _ufuncs.exprel, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy"],
            extra_note=_ufunc_kwargs_extra_note("exprel", ["torch"]),
        ),
        torch_native=False,
    ),
    _FuncInfo(
        _basic.sinc, ["x"],
        xp_capabilities(
            cpu_only=True,
            exceptions=["cupy", "torch"],
        ),
        backends_with_func_in_xp=("jax.numpy",),
        is_ufunc=False,
    ),
)

# Override ufuncs.
globals().update({nfo.func.__name__: nfo.wrapper for nfo in _special_funcs})
# digamma is an alias for psi. Define here so it also has alternative backend
# support. Add noqa because the linter gets confused by the sneaky way psi
# is inserted into globals above.
digamma = psi  # type:ignore[name-defined]  # noqa: F821
__all__ = [nfo.func.__name__ for nfo in _special_funcs] + ["digamma"]
