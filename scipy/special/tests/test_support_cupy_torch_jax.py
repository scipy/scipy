import pytest
from hypothesis import given, strategies
import hypothesis.extra.numpy as npst
from numpy.testing import assert_allclose

from scipy.special._support_cupy_torch_jax import (get_array_special_func,
                                                   array_special_func_map)
from scipy.conftest import array_api_compatible
from scipy import special
from scipy._lib.array_api_compat.array_api_compat import numpy as np


def test_dispatch_failure():
    message = ("SciPy cannot dispatch to special function 'ndtr' "
               "in array library 'object'.")
    with pytest.raises(ValueError, match=message):
        get_array_special_func('ndtr', xp=object)


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
    args_np = [np.asarray(data.draw(npst.arrays(dtype, shape)))
               for shape in shapes]
    args_xp = [xp.asarray(arg) for arg in args_np]

    ref = f(*args_np)
    res = f(*args_xp)
    # having trouble with this locally
    assert type(res) == type(xp.asarray([]))
    assert res.shape == ref.shape
    assert_allclose(np.asarray(res), ref)
