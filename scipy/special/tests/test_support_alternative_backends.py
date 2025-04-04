from functools import partial
from types import ModuleType
import pytest
from hypothesis import given, strategies
import hypothesis.extra.numpy as npst

from scipy import special
from scipy.special._support_alternative_backends import (get_array_special_func,
                                                         array_special_func_map)
from scipy._lib._array_api_no_0d import xp_assert_close
from scipy._lib._array_api import (is_cupy, is_dask, is_jax, is_torch,
                                   xp_default_dtype, SCIPY_ARRAY_API, SCIPY_DEVICE)
from scipy._lib.array_api_compat import numpy as np
from scipy._lib.array_api_extra.testing import lazy_xp_function


special_wrapped = ModuleType("special_wrapped")
lazy_xp_modules = [special_wrapped]
for f_name in array_special_func_map:
    f = getattr(special, f_name)
    setattr(special_wrapped, f_name, f)
    lazy_xp_function(f)


@pytest.mark.skipif(not SCIPY_ARRAY_API, reason="Alternative backends must be enabled.")
def test_dispatch_to_unrecognized_library():
    xp = pytest.importorskip("array_api_strict")
    f = get_array_special_func('ndtr', xp=xp)
    x = [1, 2, 3]
    res = f(xp.asarray(x))
    ref = xp.asarray(special.ndtr(np.asarray(x)))
    xp_assert_close(res, ref)


def _skip_or_tweak_alternative_backends(xp, f_name, dtypes):
    """Skip tests for specific intersections of scipy.special functions 
    vs. backends vs. dtypes vs. devices.
    Also suggest bespoke tweaks.

    Returns
    -------
    positive_only : bool
        Whether you should exclusively test positive inputs.
    dtypes_np_ref : list[str]
        The dtypes to use for the reference NumPy arrays.
    """
    if (SCIPY_DEVICE != 'cpu'
        and is_torch(xp)
        and f_name in {'stdtr', 'stdtrit', 'betaincc', 'betainc'}
    ):
        pytest.skip(f"`{f_name}` does not have an array-agnostic implementation "
                    "and cannot delegate to PyTorch.")
    if is_jax(xp) and f_name == "stdtrit":
        pytest.skip(f"`{f_name}` requires scipy.optimize support for immutable arrays")

    if ((is_jax(xp) and f_name == 'gammaincc')  # google/jax#20699
        # gh-20972
        or ((is_cupy(xp) or is_jax(xp) or is_torch(xp)) and f_name == 'chdtrc')):
        positive_only = True
    else:
        positive_only = False

    if not any('int' in dtype for dtype in dtypes):
        return positive_only, dtypes

    # Integer-specific issues from this point onwards

    if ((is_torch(xp) and f_name in {'gammainc', 'gammaincc'})
        or (is_cupy(xp) and f_name in {'stdtr', 'i0e', 'i1e'})
        or (is_jax(xp) and f_name in {'stdtr', 'ndtr', 'ndtri', 'log_ndtr'})
    ):
        pytest.skip(f"`{f_name}` does not support integer types")

    # int/float mismatched args support is sketchy
    if (any('float' in dtype for dtype in dtypes)
        and ((is_torch(xp) and f_name in ('rel_entr', 'xlogy'))
             or (is_jax(xp) and f_name in ('gammainc', 'gammaincc',
                                           'rel_entr', 'xlogy')))
    ):
        pytest.xfail("dtypes do not match")

    dtypes_np_ref = dtypes
    if (is_torch(xp) and xp_default_dtype(xp) == xp.float32
        and f_name not in {'betainc', 'betaincc', 'stdtr', 'stdtrit'}
    ):
        # On PyTorch with float32 default dtype, sometimes ints are promoted
        # to float32, and sometimes to float64.
        # When they are promoted to float32, explicitly convert the reference
        # numpy arrays to float32 to prevent them from being automatically promoted
        # to float64 instead.
        dtypes_np_ref = ['float32' if 'int' in dtype else dtype for dtype in dtypes]

    return positive_only, dtypes_np_ref


@pytest.mark.parametrize('f_name,n_args', array_special_func_map.items())
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning:dask")
@pytest.mark.parametrize('dtype', ['float32', 'float64', 'int64'])
@pytest.mark.parametrize('shapes', [[(0,)]*4, [tuple()]*4, [(10,)]*4,
                                    [(10,), (11, 1), (12, 1, 1), (13, 1, 1, 1)]])
def test_support_alternative_backends(xp, f_name, n_args, dtype, shapes):
    positive_only, [dtype_np_ref] = _skip_or_tweak_alternative_backends(
        xp, f_name, [dtype])
    f = getattr(special, f_name)  # Unwrapped
    fw = getattr(special_wrapped, f_name)  # Wrapped by lazy_xp_function
    dtype_np = getattr(np, dtype)
    dtype_xp = getattr(xp, dtype)

    shapes = shapes[:n_args]
    rng = np.random.default_rng(984254252920492019)
    if 'int' in dtype:
        iinfo = np.iinfo(dtype_np)
        rand = partial(rng.integers, iinfo.min, iinfo.max + 1)
    else:
        rand = rng.standard_normal
    args_np = [rand(size=shape, dtype=dtype_np) for shape in shapes]
    if positive_only:
        args_np = [np.abs(arg) for arg in args_np]

    args_xp = [xp.asarray(arg, dtype=dtype_xp) for arg in args_np]
    args_np = [np.asarray(arg, dtype=dtype_np_ref) for arg in args_np]

    if is_dask(xp):
        # We're using map_blocks to dispatch the function to Dask.
        # This is the correct thing to do IF all tested functions are elementwise;
        # otherwise the output would change depending on chunking.
        # Try to trigger bugs related to having multiple chunks.
        args_xp = [arg.rechunk(5) for arg in args_xp]

    res = fw(*args_xp)
    ref = f(*args_np)

    # When dtype_np is integer, the output dtype can be float
    atol = 0 if ref.dtype.kind in 'iu' else 10 * np.finfo(ref.dtype).eps
    xp_assert_close(res, xp.asarray(ref), atol=atol)


@pytest.mark.parametrize('f_name,n_args',
                         [(f_name, n_args)
                          for f_name, n_args in array_special_func_map.items()
                          if n_args >= 2])
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning:dask")
def test_support_alternative_backends_mismatched_dtypes(xp, f_name, n_args):
    """Test mix-n-match of int and float arguments"""
    assert n_args <= 3
    dtypes = ['int64', 'float32', 'float64'][:n_args]
    dtypes_xp = [xp.int64, xp.float32, xp.float64][:n_args]
    positive_only, dtypes_np_ref = _skip_or_tweak_alternative_backends(
        xp, f_name, dtypes)
    f = getattr(special, f_name)

    rng = np.random.default_rng(984254252920492019)
    iinfo = np.iinfo(np.int64)
    randint = partial(rng.integers, iinfo.min, iinfo.max + 1)
    args_np = [
        randint(size=1, dtype=np.int64),
        rng.standard_normal(size=1, dtype=np.float32),
        rng.standard_normal(size=1, dtype=np.float64),
    ][:n_args]
    if positive_only:
        args_np = [np.abs(arg) for arg in args_np]

    args_xp = [xp.asarray(arg, dtype=dtype_xp)
               for arg, dtype_xp in zip(args_np, dtypes_xp)]
    args_np = [np.asarray(arg, dtype=dtype_np_ref)
               for arg, dtype_np_ref in zip(args_np, dtypes_np_ref)]

    res = f(*args_xp)
    ref = f(*args_np)

    atol = 10 * np.finfo(ref.dtype).eps
    xp_assert_close(res, xp.asarray(ref), atol=atol)


@pytest.mark.xslow
@given(data=strategies.data())
@pytest.mark.fail_slow(5)
# `reversed` is for developer convenience: test new function first = less waiting
@pytest.mark.parametrize('f_name,n_args', reversed(array_special_func_map.items()))
@pytest.mark.filterwarnings("ignore:invalid value encountered:RuntimeWarning:dask")
@pytest.mark.filterwarnings("ignore:divide by zero encountered:RuntimeWarning:dask")
@pytest.mark.filterwarnings(
    "ignore:overflow encountered:RuntimeWarning:array_api_strict"
)
def test_support_alternative_backends_hypothesis(xp, f_name, n_args, data):
    dtype = data.draw(strategies.sampled_from(['float32', 'float64', 'int64']))
    positive_only, [dtype_np_ref] = _skip_or_tweak_alternative_backends(
        xp, f_name, [dtype])
    f = getattr(special, f_name)
    dtype_np = getattr(np, dtype)
    dtype_xp = getattr(xp, dtype)

    elements = {'allow_subnormal': False}
    # Most failures are due to NaN or infinity; uncomment to suppress them
    # elements['allow_infinity'] = False
    # elements['allow_nan'] = False
    if positive_only:
        elements['min_value'] = 0

    shapes, _ = data.draw(npst.mutually_broadcastable_shapes(num_shapes=n_args))
    args_np = [data.draw(npst.arrays(dtype_np, shape, elements=elements))
               for shape in shapes]

    args_xp = [xp.asarray(arg, dtype=dtype_xp) for arg in args_np]
    args_np = [np.asarray(arg, dtype=dtype_np_ref) for arg in args_np]

    res = f(*args_xp)
    ref = f(*args_np)

    # When dtype_np is integer, the output dtype can be float
    atol = 0 if ref.dtype.kind in 'iu' else 10 * np.finfo(ref.dtype).eps
    xp_assert_close(res, xp.asarray(ref), atol=atol)


def test_chdtr_gh21311(xp):
    # the edge case behavior of generic chdtr was not right; see gh-21311
    # be sure to test at least these cases
    # should add `np.nan` into the mix when gh-21317 is resolved
    x = np.asarray([-np.inf, -1., 0., 1., np.inf])
    v = x.reshape(-1, 1)
    ref = special.chdtr(v, x)
    res = special.chdtr(xp.asarray(v), xp.asarray(x))
    xp_assert_close(res, xp.asarray(ref))
