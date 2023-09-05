import pytest
from hypothesis import given, strategies, reproduce_failure  # noqa
import hypothesis.extra.numpy as npst

from scipy.special._support_cupy_torch_jax import (get_array_special_func,
                                                   array_special_func_map)
from scipy.conftest import array_api_compatible
from scipy import special
from scipy._lib._array_api import assert_close
from scipy._lib.array_api_compat.array_api_compat import numpy as np
import numpy.array_api as np_array_api


def test_dispatch_to_unrecognize_library():
    xp = np_array_api
    f = get_array_special_func('ndtr', xp=xp, n_array_args=1)
    x = [1, 2, 3]
    res = f(xp.asarray(x))
    ref = xp.asarray(special.ndtr(np.asarray(x)))
    assert_close(res, ref)


@array_api_compatible
@given(data=strategies.data())
@pytest.mark.parametrize('f_name_n_args', array_special_func_map.items())
def test_cupy_torch_jax_support(xp, data, f_name_n_args):
    f_name, n_args = f_name_n_args
    f = getattr(special, f_name)

    mbs = npst.mutually_broadcastable_shapes(num_shapes=n_args)
    shapes, final_shape = data.draw(mbs)

    dtype = data.draw(strategies.sampled_from(['float32', 'float64']))
    dtype = getattr(np, dtype)

    elements = dict(min_value=dtype(-10), max_value=dtype(10),
                    allow_subnormal=False)
    args_np = [np.asarray(data.draw(npst.arrays(dtype, shape, elements=elements)))
               for shape in shapes]

    # `torch.asarray(np.asarray(1.))` produces
    # TypeError: can't convert np.ndarray of type numpy.object_.
    # So we extract the scalar from 0d arrays.
    args_xp = [xp.asarray(arg[()]) for arg in args_np]

    ref = np.asarray(f(*args_np))
    res = f(*args_xp)
    # When `xp` is NumPy, `res` is a NumPy scalar, not an array.
    if res.shape != ():
        # TODO: use `_assert_matching_namespace` when gh-19005 merges
        assert isinstance(res, type(xp.asarray([])))

    eps = np.finfo(dtype).eps
    assert_close(res, xp.asarray(ref), rtol=eps**0.5, atol=eps*10)
    assert res.shape == ref.shape
    assert res.dtype == dtype
