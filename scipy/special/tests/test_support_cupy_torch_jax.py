import pytest
from hypothesis import given, strategies, reproduce_failure  # noqa
import hypothesis.extra.numpy as npst
from numpy.testing import assert_allclose

from scipy.special._support_cupy_torch_jax import (get_array_special_func,
                                                   array_special_func_map)
from scipy.conftest import array_api_compatible
from scipy import special
from scipy._lib.array_api_compat.array_api_compat import numpy as np
import numpy.array_api as np_array_api


def test_dispatch_to_unrecognize_library():
    f = get_array_special_func('ndtr', xp=np_array_api, n_array_args=1)
    x = [1, 2, 3]
    res = f(x)
    ref = special.ndtr(np.asarray(x))
    assert_allclose(np.asarray(res), ref)


@array_api_compatible
@given(data=strategies.data())
@pytest.mark.parametrize('f_name_n_args', array_special_func_map.items())
def test_cupy_torch_jax_support(xp, data, f_name_n_args):
    f_name, n_args = f_name_n_args
    f = getattr(special, f_name)
    dtypes = ['float32', 'float64']  # good enough for me
    mbs = npst.mutually_broadcastable_shapes(num_shapes=n_args)

    shapes, final_shape = data.draw(mbs)
    dtype = data.draw(strategies.sampled_from(dtypes))
    dtype = getattr(np, dtype)
    rtol = np.finfo(dtype).eps * 10
    args_np = [np.asarray(data.draw(npst.arrays(dtype, shape)))
               for shape in shapes]
    # `torch.asarray(np.asarray(1.))` produces
    # TypeError: can't convert np.ndarray of type numpy.object_.
    # So we extract the scalar from 0d arrays.
    args_xp = [xp.asarray(arg[()]) for arg in args_np]

    ref = np.asarray(f(*args_np))
    res = f(*args_xp)
    # When `xp` is NumPy, `res` NumPy scalar, not an array.
    if res.shape != ():
        # TODO: use `_assert_matching_namespace` when gh-19005 merges
        assert type(res) == type(xp.asarray([]))
    assert res.shape == ref.shape
    # TODO: use `set_assert_allclose` when gh-19005 merges
    assert_allclose(np.asarray(res), ref, rtol=rtol)
