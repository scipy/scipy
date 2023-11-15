import pytest
from hypothesis import given, strategies, reproduce_failure  # noqa: F401
import hypothesis.extra.numpy as npst

from scipy.special._support_alternative_backends import (get_array_special_func,
                                                         array_special_func_map)
from scipy.conftest import array_api_compatible
from scipy import special
from scipy._lib._array_api import xp_assert_close
from scipy._lib.array_api_compat.array_api_compat import numpy as np
import numpy.array_api as np_array_api


def test_dispatch_to_unrecognize_library():
    xp = np_array_api
    f = get_array_special_func('ndtr', xp=xp, n_array_args=1)
    x = [1, 2, 3]
    res = f(xp.asarray(x))
    ref = xp.asarray(special.ndtr(np.asarray(x)))
    xp_assert_close(res, ref)


@array_api_compatible
@given(data=strategies.data())
@pytest.mark.parametrize('f_name_n_args', array_special_func_map.items())
def test_support_alternative_backends(xp, data, f_name_n_args):
    f_name, n_args = f_name_n_args
    f = getattr(special, f_name)

    mbs = npst.mutually_broadcastable_shapes(num_shapes=n_args)
    shapes, final_shape = data.draw(mbs)

    dtype = data.draw(strategies.sampled_from(['float32', 'float64']))
    dtype_np = getattr(np, dtype)
    dtype_xp = getattr(xp, dtype)

    elements = dict(min_value=dtype_np(-10), max_value=dtype_np(10),
                    allow_subnormal=False)
    args_np = [np.asarray(data.draw(npst.arrays(dtype_np, shape, elements=elements)))
               for shape in shapes]

    # `torch.asarray(np.asarray(1.))` produces
    # TypeError: can't convert np.ndarray of type numpy.object_.
    # So we extract the scalar from 0d arrays.
    args_xp = [xp.asarray(arg[()], dtype=dtype_xp) for arg in args_np]

    ref = np.asarray(f(*args_np))
    res = f(*args_xp)

    eps = np.finfo(dtype).eps
    xp_assert_close(res, xp.asarray(ref, dtype=dtype_xp), rtol=eps**0.5, atol=eps*10,
                    check_namespace=True, check_shape=True, check_dtype=True)
